#!/usr/bin/env python3
"""Map epitope peptides to parent proteins, verify positions on the
canonical UniProt sequence, and report the best available 3D structure
(experimental PDB first, AlphaFold fallback).

Stdlib-only. See ../SKILL.md for usage.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from html import escape
from pathlib import Path
from typing import Any

IEDB_URL = "https://query-api.iedb.org/epitope_search"
IEDB_MHC_URL = "https://query-api.iedb.org/mhc_search"
IEDB_TCELL_URL = "https://query-api.iedb.org/tcell_search"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
USER_AGENT = "search-epitope-host/0.1 (+https://github.com/skills)"

# Thread-safe stderr logging and per-accession caches. The script runs the
# outer peptide loop on a ThreadPoolExecutor; without locks, log lines can
# interleave mid-message and duplicate UniProt/AlphaFold round-trips fire
# from concurrent workers resolving the same parent.
_PRINT_LOCK = threading.Lock()
_UNIPROT_CACHE: dict[str, dict | None] = {}
_UNIPROT_LOCK = threading.Lock()
_ALPHAFOLD_CACHE: dict[str, list[dict]] = {}
_ALPHAFOLD_LOCK = threading.Lock()


def _log(msg: str, log=sys.stderr) -> None:
    with _PRINT_LOCK:
        print(msg, file=log)


# ---------- HTTP helpers ----------

def _get_json(url: str, timeout: float = 30.0, retries: int = 1) -> Any:
    """GET a URL and parse JSON. One retry on transient failure. Returns
    None on 404 (resource genuinely missing); raises on other errors."""
    last_exc: Exception | None = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                       "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last_exc = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_exc = e
        if attempt < retries:
            time.sleep(1.5)
    raise RuntimeError(f"GET failed: {url}: {last_exc}")


# ---------- Input parsing ----------

def parse_inputs(positional: list[str], input_file: str | None) -> list[str]:
    peptides: list[str] = []
    if input_file:
        text = Path(input_file).read_text()
        # FASTA, delimited (CSV/TSV), or plain (one per line)
        if ">" in text:
            current: list[str] = []
            for line in text.splitlines():
                if line.startswith(">"):
                    if current:
                        peptides.append("".join(current))
                        current = []
                else:
                    current.append(line.strip())
            if current:
                peptides.append("".join(current))
        elif "," in text or "\t" in text:
            # CSV/TSV. Pick the column that looks most like peptide
            # sequences. If a header is present, prefer a column called
            # "sequence"/"peptide"/"epitope" (case-insensitive); else
            # take whichever column has the highest fraction of AA-only
            # tokens.
            delim = "," if text.count(",") >= text.count("\t") else "\t"
            rows = [r for r in csv.reader(text.splitlines(), delimiter=delim) if r]
            if not rows:
                rows = []
            header = rows[0]
            header_is_real = any(
                (cell or "").strip().lower() in {"sequence", "peptide", "epitope"}
                for cell in header) or any(
                bool(re.sub(r"[\s*]", "", cell or "").strip()) and
                set(re.sub(r"[\s*]", "", (cell or "")).upper()) - AA_SET
                for cell in header)
            body = rows[1:] if header_is_real else rows
            col_idx = 0
            if header_is_real:
                lowered = [(cell or "").strip().lower() for cell in header]
                for want in ("sequence", "peptide", "epitope"):
                    if want in lowered:
                        col_idx = lowered.index(want)
                        break
                else:
                    # Pick column with highest fraction of AA-only cells.
                    best, best_score = 0, -1.0
                    for j in range(len(header)):
                        ok = 0
                        for r in body:
                            if j < len(r):
                                s = re.sub(r"[\s*]", "", r[j]).upper()
                                if s and not (set(s) - AA_SET):
                                    ok += 1
                        score = ok / max(1, len(body))
                        if score > best_score:
                            best, best_score = j, score
                    col_idx = best
            for r in body:
                if col_idx < len(r):
                    peptides.append(r[col_idx])
        else:
            peptides.extend(line.strip() for line in text.splitlines() if line.strip())
    peptides.extend(positional)
    # Normalize & validate
    out: list[str] = []
    seen: set[str] = set()
    for p in peptides:
        seq = re.sub(r"[\s*]", "", p).upper()
        if not seq or seq in seen:
            continue
        seen.add(seq)
        out.append(seq)
    return out


def validate_peptide(seq: str) -> str | None:
    """Return None if valid, else error message."""
    if len(seq) < 5:
        return "sequence shorter than 5 residues"
    bad = set(seq) - AA_SET
    if bad:
        return f"non-standard residues: {''.join(sorted(bad))}"
    return None


# ---------- IEDB ----------

_UNIPROT_NAME_DECORATION = re.compile(r"\s*\(UniProt:[A-Z0-9-]+\)\s*$")


def _parse_uniprot_iri(iri: str) -> str | None:
    """`UNIPROT:P01012.2` → `P01012`. Returns None for non-UNIPROT IRIs."""
    if not iri or not iri.startswith("UNIPROT:"):
        return None
    return iri.split(":", 1)[1].split(".", 1)[0]


def query_iedb(seq: str) -> list[dict]:
    """Return per-antigen rows from IEDB epitope_search, with epitope-level
    metadata (assay outcomes, MHC) repeated on each row for downstream
    grouping.

    Each item: {
        'uniprot_acc', 'antigen_name', 'organism',
        'starting_position', 'ending_position', 'iedb_structure_id',
        'assay_outcomes': list[str],       # qualitative_measures
        'mhc_classes': list[str],
        'mhc_alleles': list[str],
    }.

    Parent resolution: first take every `UNIPROT:` IRI under
    `curated_source_antigens` (keeps the IEDB-claimed start/end positions);
    then top up with any extra accessions from `parent_source_antigen_iris`
    (IEDB's UniProt-normalized parent), with no positions — `verify_position`
    handles that case via substring scan against the canonical sequence.
    """
    select = ",".join([
        "structure_id", "linear_sequence", "curated_source_antigens",
        "parent_source_antigen_iris", "parent_source_antigen_names",
        "parent_source_antigen_source_org_names",
        "qualitative_measures", "mhc_classes", "mhc_allele_names",
    ])
    params = {"linear_sequence": f"eq.{seq}", "select": select}
    url = f"{IEDB_URL}?{urllib.parse.urlencode(params)}"
    rows = _get_json(url) or []
    out: list[dict] = []
    for row in rows:
        sid = row.get("structure_id")
        meta = {
            "assay_outcomes": list(row.get("qualitative_measures") or []),
            "mhc_classes": list(row.get("mhc_classes") or []),
            "mhc_alleles": list(row.get("mhc_allele_names") or []),
        }

        per_acc: dict[str, dict] = {}
        for ant in row.get("curated_source_antigens") or []:
            acc = _parse_uniprot_iri(ant.get("iri") or "")
            if not acc:
                continue
            per_acc.setdefault(acc, {
                "uniprot_acc": acc,
                "antigen_name": ant.get("name"),
                "organism": ant.get("source_organism_name"),
                "starting_position": ant.get("starting_position"),
                "ending_position": ant.get("ending_position"),
                "iedb_structure_id": sid,
                **meta,
            })

        parent_iris = row.get("parent_source_antigen_iris") or []
        parent_names = row.get("parent_source_antigen_names") or []
        parent_orgs = row.get("parent_source_antigen_source_org_names") or []
        for i, p_iri in enumerate(parent_iris):
            acc = _parse_uniprot_iri(p_iri or "")
            if not acc or acc in per_acc:
                continue
            name = parent_names[i] if i < len(parent_names) else None
            if name:
                name = _UNIPROT_NAME_DECORATION.sub("", name)
            org = parent_orgs[i] if i < len(parent_orgs) else None
            per_acc[acc] = {
                "uniprot_acc": acc,
                "antigen_name": name,
                "organism": org,
                "starting_position": None,
                "ending_position": None,
                "iedb_structure_id": sid,
                **meta,
            }

        out.extend(per_acc.values())
    return out


def query_iedb_assays(seq: str) -> dict[str, dict[str, set[str]]]:
    """Pull per-assay (allele, outcome) pairs from `mhc_search` and
    `tcell_search`. Returns:

        {
            'mhc': {allele_name: {outcomes}},
            'tcell': {allele_name: {outcomes}},
        }

    where outcomes are raw IEDB `qualitative_measure` strings (`Positive`,
    `Positive-High`, …, `Negative`). Allele names use IEDB's
    `mhc_allele_name` field — typically `HLA-A*02:01`, `H2-Kb`, etc.
    Unknown / unrestricted alleles are dropped. Each endpoint is queried
    with a tight `select` so payloads stay small even for popular epitopes.
    """
    select = "mhc_allele_name,qualitative_measure,mhc_class"
    endpoints = (("mhc", IEDB_MHC_URL), ("tcell", IEDB_TCELL_URL))

    def _fetch(base: str) -> list[dict]:
        url = f"{base}?{urllib.parse.urlencode({'linear_sequence': f'eq.{seq}', 'select': select})}"
        try:
            return _get_json(url) or []
        except Exception:  # noqa: BLE001 — assay endpoint failures shouldn't kill the epitope
            return []

    # Fire mhc + tcell concurrently — both are network-bound and independent.
    # Tiny inline pool: per-call overhead is in the µs range and pays back
    # ~200 ms on every matched peptide.
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = {key: pool.submit(_fetch, base) for key, base in endpoints}
        responses = {key: fut.result() for key, fut in results.items()}

    out: dict[str, dict[str, set[str]]] = {"mhc": {}, "tcell": {}}
    for key, rows in responses.items():
        for r in rows:
            allele = r.get("mhc_allele_name")
            measure = r.get("qualitative_measure")
            if not allele or not measure:
                continue
            out[key].setdefault(allele, set()).add(measure)
    return out


# ---------- UniProt ----------

def fetch_uniprot(acc: str) -> dict | None:
    with _UNIPROT_LOCK:
        if acc in _UNIPROT_CACHE:
            return _UNIPROT_CACHE[acc]
    data = _get_json(UNIPROT_URL.format(acc=acc))
    if not data:
        with _UNIPROT_LOCK:
            _UNIPROT_CACHE[acc] = None
        return None
    seq = (data.get("sequence") or {}).get("value", "")
    name = ((data.get("proteinDescription") or {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value")) or data.get("uniProtkbId")
    organism = (data.get("organism") or {}).get("scientificName")
    pdb_xrefs: list[dict] = []
    for xref in data.get("uniProtKBCrossReferences") or []:
        if xref.get("database") != "PDB":
            continue
        props = {p["key"]: p["value"] for p in xref.get("properties") or []}
        pdb_xrefs.append({
            "pdb_id": xref.get("id"),
            "method": props.get("Method"),
            "resolution": props.get("Resolution"),
            "chains": props.get("Chains"),
        })
    result = {
        "accession": acc,
        "name": name,
        "organism": organism,
        "sequence": seq,
        "length": len(seq),
        "pdb_xrefs": pdb_xrefs,
    }
    with _UNIPROT_LOCK:
        _UNIPROT_CACHE[acc] = result
    return result


def verify_position(seq: str, epitope: str, start: int | None, end: int | None) -> dict:
    """Check whether the epitope sits at the claimed 1-indexed [start, end]
    on the full sequence, and find every occurrence in the sequence."""
    claimed_match = False
    if seq and start and end and 1 <= start <= len(seq) and end <= len(seq):
        claimed_match = seq[start - 1:end] == epitope
    occurrences: list[tuple[int, int]] = []
    i = 0
    while True:
        j = seq.find(epitope, i)
        if j < 0:
            break
        occurrences.append((j + 1, j + len(epitope)))
        i = j + 1
    return {
        "claimed_start": start,
        "claimed_end": end,
        "claimed_match": claimed_match,
        "occurrences": occurrences,
    }


# ---------- AlphaFold ----------

def fetch_alphafold(acc: str) -> list[dict]:
    with _ALPHAFOLD_LOCK:
        if acc in _ALPHAFOLD_CACHE:
            return _ALPHAFOLD_CACHE[acc]
    data = _get_json(ALPHAFOLD_URL.format(acc=acc)) or []
    out: list[dict] = []
    for entry in data:
        out.append({
            "entry_id": entry.get("entryId"),
            "pdb_url": entry.get("pdbUrl"),
            "cif_url": entry.get("cifUrl"),
            "pae_url": entry.get("paeImageUrl"),
            "version": entry.get("latestVersion"),
            "uniprot_start": entry.get("uniprotStart"),
            "uniprot_end": entry.get("uniprotEnd"),
        })
    with _ALPHAFOLD_LOCK:
        _ALPHAFOLD_CACHE[acc] = out
    return out


# ---------- Aggregation ----------

def process_epitope(seq: str, use_alphafold: bool = True,
                    log=sys.stderr) -> dict:
    err = validate_peptide(seq)
    if err:
        return {"epitope": seq, "status": "invalid", "error": err}
    _log(f"[IEDB] {seq}", log)
    antigens = query_iedb(seq)
    if not antigens:
        return {"epitope": seq, "status": "unmatched", "parents": [],
                "assay_outcomes": [], "mhc_classes": [], "mhc_alleles": [],
                "hla_outcomes": {}}
    # Per-allele (HLA, outcome) pairs from the assay endpoints. Two extra
    # round-trips per matched epitope (mhc + tcell).
    _log(f"[assays] {seq}", log)
    try:
        assay_idx = query_iedb_assays(seq)
    except Exception as e:  # noqa: BLE001
        _log(f"[error] assay query {seq}: {e}", log)
        assay_idx = {"mhc": {}, "tcell": {}}
    # Merge mhc + tcell into a single {allele: set[outcomes]} map.
    hla_outcomes: dict[str, set[str]] = {}
    for d in (assay_idx["mhc"], assay_idx["tcell"]):
        for allele, outs in d.items():
            hla_outcomes.setdefault(allele, set()).update(outs)
    # Epitope-level metadata is repeated across antigen rows; dedup across all.
    ep_outcomes: list[str] = []
    ep_mhc_classes: list[str] = []
    ep_mhc_alleles: list[str] = []
    _seen_o, _seen_c, _seen_a = set(), set(), set()
    for a in antigens:
        for v in a.get("assay_outcomes") or []:
            if v not in _seen_o:
                _seen_o.add(v); ep_outcomes.append(v)
        for v in a.get("mhc_classes") or []:
            if v not in _seen_c:
                _seen_c.add(v); ep_mhc_classes.append(v)
        for v in a.get("mhc_alleles") or []:
            if v not in _seen_a:
                _seen_a.add(v); ep_mhc_alleles.append(v)
    # Group by accession; merge claimed positions
    grouped: dict[str, dict] = {}
    for a in antigens:
        g = grouped.setdefault(a["uniprot_acc"], {
            "uniprot_acc": a["uniprot_acc"],
            "antigen_names": set(),
            "organisms": set(),
            "claimed_positions": [],
            "iedb_structure_ids": set(),
        })
        if a["antigen_name"]:
            g["antigen_names"].add(a["antigen_name"])
        if a["organism"]:
            g["organisms"].add(a["organism"])
        if a["starting_position"] and a["ending_position"]:
            g["claimed_positions"].append((a["starting_position"], a["ending_position"]))
        if a["iedb_structure_id"]:
            g["iedb_structure_ids"].add(a["iedb_structure_id"])

    parents: list[dict] = []
    for acc, g in grouped.items():
        _log(f"[UniProt] {acc}", log)
        try:
            up = fetch_uniprot(acc)
        except Exception as e:  # noqa: BLE001
            _log(f"[error] {acc}: {e}", log)
            parents.append({
                "uniprot_acc": acc,
                "status": "uniprot_fetch_failed",
                "error": str(e),
                "antigen_names": sorted(g["antigen_names"]),
                "organisms": sorted(g["organisms"]),
                "claimed_positions": g["claimed_positions"],
            })
            continue
        if not up:
            parents.append({
                "uniprot_acc": acc,
                "status": "uniprot_not_found",
                "antigen_names": sorted(g["antigen_names"]),
                "organisms": sorted(g["organisms"]),
                "claimed_positions": g["claimed_positions"],
            })
            continue
        # Verify each claimed position
        position_checks = []
        for s, e in g["claimed_positions"] or [(None, None)]:
            position_checks.append(verify_position(up["sequence"], seq, s, e))
        af: list[dict] = []
        if not up["pdb_xrefs"] and use_alphafold:
            _log(f"[AlphaFold] {acc}", log)
            af = fetch_alphafold(acc)
        parents.append({
            "uniprot_acc": acc,
            "status": "ok",
            "name": up["name"],
            "organism": up["organism"],
            "length": up["length"],
            "sequence": up["sequence"],
            "antigen_names": sorted(g["antigen_names"]),
            "organisms_iedb": sorted(g["organisms"]),
            "iedb_structure_ids": sorted(g["iedb_structure_ids"]),
            "position_checks": position_checks,
            "pdb_xrefs": up["pdb_xrefs"],
            "alphafold": af,
        })
    return {"epitope": seq, "status": "ok", "parents": parents,
            "assay_outcomes": ep_outcomes,
            "mhc_classes": ep_mhc_classes,
            "mhc_alleles": ep_mhc_alleles,
            "hla_outcomes": hla_outcomes}


# ---------- Rendering ----------

def _pos_summary(checks: list[dict]) -> str:
    parts = []
    for c in checks:
        if c["claimed_start"]:
            tag = "✓" if c["claimed_match"] else "✗"
            parts.append(f"{c['claimed_start']}-{c['claimed_end']} {tag}")
        else:
            parts.append("(no IEDB position)")
    occ = checks[0]["occurrences"] if checks else []
    extra = ""
    if occ:
        extra = " | found at " + ", ".join(f"{s}-{e}" for s, e in occ)
    return "; ".join(parts) + extra


def _summary_rows(results: list[dict]) -> list[dict]:
    """Flat one-row-per-(epitope, parent) summary used for the top table
    and the master CSV. Unmatched/invalid epitopes get a single row."""
    rows: list[dict] = []
    for r in results:
        ep = r["epitope"]
        status = r.get("status")
        if status == "ok":
            for p in r["parents"]:
                if p.get("status") != "ok":
                    rows.append({
                        "epitope": ep,
                        "length": len(ep),
                        "status": "uniprot_fetch_failed",
                        "uniprot_acc": p["uniprot_acc"],
                        "protein_name": "",
                        "organism": "",
                        "assay_outcomes": "; ".join(r.get("assay_outcomes") or []),
                        "mhc_classes": "; ".join(r.get("mhc_classes") or []),
                        "mhc_alleles": _fmt_alleles(r.get("mhc_alleles") or []),
                        "hla_outcomes": _fmt_hla_outcomes(r.get("hla_outcomes") or {}),
                        "uniprot_length": "",
                        "claimed_positions": "",
                        "position_verified": "",
                        "occurrences": "",
                        "n_pdb": 0,
                        "has_alphafold": False,
                    })
                    continue
                claimed = "; ".join(
                    f"{c['claimed_start']}-{c['claimed_end']}"
                    for c in p["position_checks"] if c["claimed_start"])
                verified = "; ".join(
                    ("yes" if c["claimed_match"] else "no")
                    for c in p["position_checks"] if c["claimed_start"])
                occ = p["position_checks"][0]["occurrences"] if p["position_checks"] else []
                rows.append({
                    "epitope": ep,
                    "length": len(ep),
                    "status": "matched",
                    "uniprot_acc": p["uniprot_acc"],
                    "protein_name": p["name"] or "",
                    "organism": p["organism"] or "",
                    "assay_outcomes": "; ".join(r.get("assay_outcomes") or []),
                    "mhc_classes": "; ".join(r.get("mhc_classes") or []),
                    "mhc_alleles": _fmt_alleles(r.get("mhc_alleles") or []),
                    "hla_outcomes": _fmt_hla_outcomes(r.get("hla_outcomes") or {}),
                    "uniprot_length": p["length"],
                    "claimed_positions": claimed,
                    "position_verified": verified,
                    "occurrences": ", ".join(f"{s}-{e}" for s, e in occ),
                    "n_pdb": len(p["pdb_xrefs"]),
                    "has_alphafold": bool(p["alphafold"]),
                })
        elif status == "unmatched":
            rows.append({
                "epitope": ep, "length": len(ep), "status": "no_iedb_match",
                "uniprot_acc": "", "protein_name": "", "organism": "",
                "assay_outcomes": "", "mhc_classes": "", "mhc_alleles": "",
                "hla_outcomes": "",
                "uniprot_length": "", "claimed_positions": "",
                "position_verified": "", "occurrences": "",
                "n_pdb": 0, "has_alphafold": False,
            })
        elif status == "invalid":
            rows.append({
                "epitope": ep, "length": len(ep),
                "status": f"invalid: {r.get('error', '')}",
                "uniprot_acc": "", "protein_name": "", "organism": "",
                "assay_outcomes": "", "mhc_classes": "", "mhc_alleles": "",
                "hla_outcomes": "",
                "uniprot_length": "", "claimed_positions": "",
                "position_verified": "", "occurrences": "",
                "n_pdb": 0, "has_alphafold": False,
            })
        elif status == "error":
            rows.append({
                "epitope": ep, "length": len(ep),
                "status": f"error: {r.get('error', '')}",
                "uniprot_acc": "", "protein_name": "", "organism": "",
                "assay_outcomes": "", "mhc_classes": "", "mhc_alleles": "",
                "hla_outcomes": "",
                "uniprot_length": "", "claimed_positions": "",
                "position_verified": "", "occurrences": "",
                "n_pdb": 0, "has_alphafold": False,
            })
    return rows


SUMMARY_COLUMNS = [
    "epitope", "length", "status", "uniprot_acc", "protein_name",
    "organism", "assay_outcomes", "mhc_classes", "mhc_alleles",
    "hla_outcomes",
    "uniprot_length", "claimed_positions", "position_verified",
    "occurrences", "n_pdb", "has_alphafold",
]


def _fmt_alleles(alleles: list[str], cap: int = 5) -> str:
    if not alleles:
        return ""
    if len(alleles) <= cap:
        return "; ".join(alleles)
    return "; ".join(alleles[:cap]) + f"; …(+{len(alleles) - cap})"


def _outcome_sort_key(o: str) -> int:
    # Order outcomes by binding strength for stable display:
    # Negative < Positive-Low < Positive-Intermediate < Positive-High < Positive
    order = {
        "Negative": 0,
        "Positive-Low": 1,
        "Positive-Intermediate": 2,
        "Positive-High": 3,
        "Positive": 4,
    }
    return order.get(o, 5)


def _fmt_hla_outcomes(hla: dict[str, set[str]] | None, cap: int = 8) -> str:
    """`{HLA-A*02:01: {Negative, Positive}}` → `HLA-A*02:01[Negative,Positive]`.
    Cap at `cap` alleles, with `…(+N)` suffix when truncated. Alleles sorted
    alphabetically; outcomes within each allele sorted by binding strength."""
    if not hla:
        return ""
    items = []
    for allele in sorted(hla.keys()):
        outs = sorted(hla[allele], key=_outcome_sort_key)
        items.append(f"{allele}[{','.join(outs)}]")
    if len(items) <= cap:
        return "; ".join(items)
    return "; ".join(items[:cap]) + f"; …(+{len(items) - cap})"


def write_master_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow(r)


STRUCTURE_COLUMNS = [
    "uniprot_acc", "source", "entry_id", "method", "resolution",
    "chains", "range", "version", "url",
]


def _structure_rows_for_parent(p: dict) -> list[dict]:
    out: list[dict] = []
    acc = p["uniprot_acc"]
    for x in p["pdb_xrefs"]:
        out.append({
            "uniprot_acc": acc,
            "source": "PDB",
            "entry_id": x["pdb_id"],
            "method": x["method"] or "",
            "resolution": x["resolution"] or "",
            "chains": x["chains"] or "",
            "range": "",
            "version": "",
            "url": f"https://www.rcsb.org/structure/{x['pdb_id']}",
        })
    if not p["pdb_xrefs"]:
        for a in p["alphafold"]:
            out.append({
                "uniprot_acc": acc,
                "source": "AlphaFold",
                "entry_id": a["entry_id"] or "",
                "method": "predicted",
                "resolution": "",
                "chains": "",
                "range": f"{a['uniprot_start']}-{a['uniprot_end']}",
                "version": f"v{a['version']}" if a["version"] else "",
                "url": a["pdb_url"] or a["cif_url"] or "",
            })
    return out


def write_parent_csvs(out_dir: Path, results: list[dict]) -> dict[str, Path]:
    """Write one structures CSV per matched parent UniProt accession.
    Returns {acc: relative_path} for use in the report."""
    sub = out_dir / "structures"
    sub.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    seen: set[str] = set()
    for r in results:
        if r.get("status") != "ok":
            continue
        for p in r["parents"]:
            if p.get("status") != "ok":
                continue
            acc = p["uniprot_acc"]
            if acc in seen:
                continue
            seen.add(acc)
            rows = _structure_rows_for_parent(p)
            path = sub / f"{acc}.csv"
            with path.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=STRUCTURE_COLUMNS)
                w.writeheader()
                for row in rows:
                    w.writerow(row)
            written[acc] = path.relative_to(out_dir)
    return written


def render_markdown(results: list[dict], summary_rows: list[dict],
                    parent_csv_map: dict[str, Path],
                    fasta_map: dict[str, Path] | None = None,
                    map_map: dict[str, Path] | None = None) -> str:
    fasta_map = fasta_map or {}
    map_map = map_map or {}
    lines = ["# Epitope → host protein report", ""]
    matched = [r for r in results if r.get("status") == "ok"]
    unmatched = [r for r in results if r.get("status") == "unmatched"]
    invalid = [r for r in results if r.get("status") == "invalid"]
    lines.append(f"- Input epitopes: **{len(results)}**")
    lines.append(f"- Matched in IEDB: **{len(matched)}**")
    lines.append(f"- Unmatched: **{len(unmatched)}**")
    if invalid:
        lines.append(f"- Invalid input: **{len(invalid)}**")
    lines.append("")

    # --- Top summary table: one row per (epitope, parent) ---
    lines.append("## Summary")
    lines.append("")
    lines.append("| Epitope | Len | Status | UniProt | Protein | Organism | "
                 "Assay outcome | MHC class | HLA outcomes (per-allele) | "
                 "Position (claimed → verified) | PDB | AF |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for row in summary_rows:
        acc = row["uniprot_acc"]
        acc_cell = (f"[{acc}](https://www.uniprot.org/uniprotkb/{acc})"
                    if acc else "–")
        pos_cell = ""
        if row["claimed_positions"]:
            verified = row["position_verified"].split("; ")
            claimed = row["claimed_positions"].split("; ")
            pos_cell = "; ".join(
                f"{c} {'✓' if v == 'yes' else '✗'}"
                for c, v in zip(claimed, verified))
            if row["occurrences"]:
                pos_cell += f" | found at {row['occurrences']}"
        lines.append(
            f"| `{row['epitope']}` | {row['length']} | {row['status']} "
            f"| {acc_cell} | {row['protein_name']} | {row['organism']} "
            f"| {row.get('assay_outcomes', '')} | {row.get('mhc_classes', '')} "
            f"| {row.get('hla_outcomes', '')} "
            f"| {pos_cell} | {row['n_pdb']} "
            f"| {'yes' if row['has_alphafold'] else ''} |")
    lines.append("")
    lines.append("_Master CSV: [`epitopes.csv`](epitopes.csv). "
                 "Per-parent structure CSVs in [`structures/`](structures/)._")
    lines.append("")

    lines.append("## Per-epitope details")
    lines.append("")
    for r in matched:
        ep = r["epitope"]
        lines.append(f"## Epitope `{ep}` ({len(ep)} aa)")
        lines.append("")
        lines.append("### Parent proteins")
        lines.append("")
        lines.append("| UniProt | Name | Organism | Length | Position (claimed → verified) |")
        lines.append("|---|---|---|---|---|")
        for p in r["parents"]:
            if p["status"] != "ok":
                lines.append(f"| {p['uniprot_acc']} | _UniProt fetch failed_ | – | – | – |")
                continue
            lines.append(f"| [{p['uniprot_acc']}](https://www.uniprot.org/uniprotkb/{p['uniprot_acc']}) "
                         f"| {p['name'] or ''} | {p['organism'] or ''} | {p['length']} "
                         f"| {_pos_summary(p['position_checks'])} |")
        lines.append("")
        for p in r["parents"]:
            if p["status"] != "ok":
                continue
            acc = p["uniprot_acc"]
            link_bits = []
            if acc in parent_csv_map:
                link_bits.append(f"[structures CSV]({parent_csv_map[acc]})")
            if acc in fasta_map:
                link_bits.append(f"[FASTA]({fasta_map[acc]})")
            if acc in map_map:
                link_bits.append(f"[epitope map]({map_map[acc]})")
            link_suffix = " — " + " · ".join(link_bits) if link_bits else ""
            lines.append(f"### Structures for {acc} – "
                         f"{p['name'] or ''}{link_suffix}")
            lines.append("")
            if p["pdb_xrefs"]:
                lines.append(f"_Experimental PDB entries ({len(p['pdb_xrefs'])}):_")
                lines.append("")
                lines.append("| PDB | Method | Resolution (Å) | Chains |")
                lines.append("|---|---|---|---|")
                for x in p["pdb_xrefs"]:
                    lines.append(f"| [{x['pdb_id']}](https://www.rcsb.org/structure/{x['pdb_id']}) "
                                 f"| {x['method'] or ''} | {x['resolution'] or ''} | {x['chains'] or ''} |")
                lines.append("")
            elif p["alphafold"]:
                lines.append("_No experimental PDB entry — using AlphaFold:_")
                lines.append("")
                lines.append("| AlphaFold entry | Version | Range | PDB | CIF |")
                lines.append("|---|---|---|---|---|")
                for a in p["alphafold"]:
                    eid = a["entry_id"]
                    page = f"https://alphafold.ebi.ac.uk/entry/{p['uniprot_acc']}"
                    lines.append(f"| [{eid}]({page}) | v{a['version']} "
                                 f"| {a['uniprot_start']}-{a['uniprot_end']} "
                                 f"| [pdb]({a['pdb_url']}) | [cif]({a['cif_url']}) |")
                lines.append("")
            else:
                lines.append("_No experimental PDB and no AlphaFold prediction available._")
                lines.append("")

    if unmatched:
        lines.append("## Unmatched epitopes (no IEDB hit)")
        lines.append("")
        for r in unmatched:
            lines.append(f"- `{r['epitope']}`")
        lines.append("")
    if invalid:
        lines.append("## Invalid input")
        lines.append("")
        for r in invalid:
            lines.append(f"- `{r['epitope']}` — {r['error']}")
        lines.append("")
    return "\n".join(lines)


def render_html(results: list[dict], summary_rows: list[dict],
                parent_csv_map: dict[str, Path],
                fasta_map: dict[str, Path] | None = None,
                map_map: dict[str, Path] | None = None) -> str:
    fasta_map = fasta_map or {}
    map_map = map_map or {}
    css = """
      body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
           max-width:1100px;margin:2em auto;padding:0 1em;color:#222;line-height:1.5}
      h1{border-bottom:2px solid #444;padding-bottom:.3em}
      h2{margin-top:2em;border-bottom:1px solid #ccc;padding-bottom:.2em}
      h3{margin-top:1.4em;color:#555}
      table{border-collapse:collapse;width:100%;margin:.6em 0;font-size:14px}
      th,td{border:1px solid #ccc;padding:.4em .6em;text-align:left;vertical-align:top}
      th{background:#f4f4f4}
      code{background:#f0f0f0;padding:.1em .35em;border-radius:3px}
      .ok{color:#197a19;font-weight:bold}
      .bad{color:#b00020;font-weight:bold}
      .muted{color:#888}
      .pill{display:inline-block;padding:.1em .55em;border-radius:10px;
            background:#eef;color:#225;font-size:12px;margin-left:.4em}
    """
    parts: list[str] = []
    parts.append(f"<!doctype html><html><head><meta charset='utf-8'>"
                 f"<title>Epitope → host protein report</title>"
                 f"<style>{css}</style></head><body>")
    parts.append("<h1>Epitope → host protein report</h1>")
    matched = [r for r in results if r.get("status") == "ok"]
    unmatched = [r for r in results if r.get("status") == "unmatched"]
    invalid = [r for r in results if r.get("status") == "invalid"]
    parts.append("<ul>"
                 f"<li>Input epitopes: <b>{len(results)}</b></li>"
                 f"<li>Matched in IEDB: <b>{len(matched)}</b></li>"
                 f"<li>Unmatched: <b>{len(unmatched)}</b></li>"
                 + (f"<li>Invalid input: <b>{len(invalid)}</b></li>" if invalid else "")
                 + "</ul>")

    # Top summary table
    parts.append("<h2>Summary</h2>")
    parts.append("<table><tr><th>Epitope</th><th>Len</th><th>Status</th>"
                 "<th>UniProt</th><th>Protein</th><th>Organism</th>"
                 "<th>Assay outcome</th><th>MHC class</th>"
                 "<th>HLA outcomes (per-allele)</th>"
                 "<th>Position (claimed → verified)</th><th>PDB</th>"
                 "<th>AF</th></tr>")
    for row in summary_rows:
        acc = row["uniprot_acc"]
        acc_cell = (f"<a href='https://www.uniprot.org/uniprotkb/{escape(acc)}'>"
                    f"{escape(acc)}</a>" if acc else "<span class='muted'>–</span>")
        pos_cell = ""
        if row["claimed_positions"]:
            verified = row["position_verified"].split("; ")
            claimed = row["claimed_positions"].split("; ")
            pieces = []
            for c, v in zip(claimed, verified):
                cls = "ok" if v == "yes" else "bad"
                sym = "✓" if v == "yes" else "✗"
                pieces.append(f"<span class='{cls}'>{escape(c)} {sym}</span>")
            pos_cell = "; ".join(pieces)
            if row["occurrences"]:
                pos_cell += (" | <span class='muted'>found at "
                             f"{escape(row['occurrences'])}</span>")
        outcome_cell = escape(row.get("assay_outcomes", ""))
        if "Positive" in row.get("assay_outcomes", "") and "Negative" in row.get("assay_outcomes", ""):
            outcome_cell = f"<span class='bad'>{outcome_cell}</span>"
        elif row.get("assay_outcomes", "").startswith("Negative"):
            outcome_cell = f"<span class='muted'>{outcome_cell}</span>"
        parts.append(
            f"<tr><td><code>{escape(row['epitope'])}</code></td>"
            f"<td>{row['length']}</td>"
            f"<td>{escape(row['status'])}</td>"
            f"<td>{acc_cell}</td>"
            f"<td>{escape(row['protein_name'])}</td>"
            f"<td>{escape(row['organism'])}</td>"
            f"<td>{outcome_cell}</td>"
            f"<td>{escape(row.get('mhc_classes', ''))}</td>"
            f"<td><small>{escape(row.get('hla_outcomes', ''))}</small></td>"
            f"<td>{pos_cell}</td>"
            f"<td>{row['n_pdb']}</td>"
            f"<td>{'yes' if row['has_alphafold'] else ''}</td></tr>")
    parts.append("</table>")
    parts.append("<p class='muted'>Master CSV: "
                 "<a href='epitopes.csv'>epitopes.csv</a>. "
                 "Per-parent structure CSVs in "
                 "<a href='structures/'>structures/</a>.</p>")

    parts.append(
        "<p class='muted'>Per-epitope details (parent tables + PDB/AlphaFold "
        "structures + per-parent links) are omitted from this HTML to keep "
        "page size manageable for large inputs. See <a href='report.md'>"
        "<code>report.md</code></a> for the full per-epitope writeup, or "
        "follow the UniProt / map links per row above.</p>")

    if unmatched:
        parts.append("<h2>Unmatched epitopes (no IEDB hit)</h2><ul>")
        for r in unmatched:
            parts.append(f"<li><code>{escape(r['epitope'])}</code></li>")
        parts.append("</ul>")
    if invalid:
        parts.append("<h2>Invalid input</h2><ul>")
        for r in invalid:
            parts.append(f"<li><code>{escape(r['epitope'])}</code> — {escape(r['error'])}</li>")
        parts.append("</ul>")
    parts.append("</body></html>")
    return "\n".join(parts)


# ---------- Per-parent FASTA + epitope map ----------

EPITOPE_COLORS = [
    "#ffd54f", "#80cbc4", "#f48fb1", "#a5d6a7", "#90caf9",
    "#ce93d8", "#ffab91", "#fff59d", "#b39ddb", "#bcaaa4",
]


def build_parent_index(results: list[dict]) -> dict[str, dict]:
    """Group successful parents across all epitopes.

    Returns {acc: {name, organism, length, sequence,
                   epitopes: [{epitope, occurrences:[(s,e),...]}, ...]}}.
    Only parents whose UniProt fetch succeeded and that have at least one
    occurrence of the epitope in the canonical sequence are kept.
    """
    index: dict[str, dict] = {}
    for r in results:
        if r.get("status") != "ok":
            continue
        ep = r["epitope"]
        for p in r["parents"]:
            if p.get("status") != "ok":
                continue
            acc = p["uniprot_acc"]
            occ = p["position_checks"][0]["occurrences"] if p["position_checks"] else []
            if not occ:
                continue
            slot = index.setdefault(acc, {
                "uniprot_acc": acc,
                "name": p["name"],
                "organism": p["organism"],
                "length": p["length"],
                "sequence": p["sequence"],
                "epitopes": [],
            })
            # If the same epitope is reported twice for the same parent
            # (e.g., from multiple IEDB rows), merge.
            existing = next((e for e in slot["epitopes"] if e["epitope"] == ep), None)
            if existing is None:
                slot["epitopes"].append({"epitope": ep, "occurrences": list(occ)})
            else:
                for pos in occ:
                    if pos not in existing["occurrences"]:
                        existing["occurrences"].append(pos)
    return index


def _wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_parent_fastas(out_dir: Path, parent_index: dict[str, dict],
                        log=sys.stderr) -> dict[str, Path]:
    sub = out_dir / "parents"
    sub.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for acc, p in parent_index.items():
        header = f">{acc} {p['name'] or ''} OS={p['organism'] or ''} LEN={p['length']}"
        path = sub / f"{acc}.fasta"
        path.write_text(header + "\n" + _wrap_fasta(p["sequence"]) + "\n")
        written[acc] = path.relative_to(out_dir)
        print(f"wrote {path}", file=log)
    return written


def _render_sequence_block(seq: str, ranges: list[tuple[int, int, int]],
                           per_line: int = 60, block: int = 10) -> str:
    """Render `seq` as an HTML <pre> block. `ranges` is a list of
    (start, end, color_idx) 1-indexed inclusive ranges. Residues in any
    range get a colored <span> background — the FIRST matching range wins
    on overlap (lowest color_idx)."""
    # marks[i-1] = color_idx for 1-indexed position i, or None
    marks: list[int | None] = [None] * len(seq)
    for s, e, idx in ranges:
        for i in range(s - 1, min(e, len(seq))):
            if marks[i] is None:
                marks[i] = idx
    lines = []
    n = len(seq)
    width = len(str(n))
    for start in range(0, n, per_line):
        end = min(start + per_line, n)
        prefix = f"{start + 1:>{width}}  "
        cells = []
        cur_color: int | None = "__init__"  # sentinel
        run: list[str] = []

        def flush():
            if not run:
                return
            text = "".join(run)
            if cur_color is None:
                cells.append(text)
            else:
                color = EPITOPE_COLORS[cur_color % len(EPITOPE_COLORS)]
                cells.append(f"<span style='background:{color}'>{text}</span>")
        for i in range(start, end):
            c = marks[i]
            if c != cur_color:
                flush()
                run = []
                cur_color = c
            run.append(seq[i])
            if (i + 1) % block == 0 and i + 1 < end:
                # insert a space between 10-aa blocks, outside any span
                flush()
                run = []
                cells.append(" ")
                # keep cur_color so the next residue starts a fresh span
                cur_color = "__init__"
        flush()
        lines.append(prefix + "".join(cells))
    return "<pre style='font-family:monospace;font-size:13px;line-height:1.45;" \
           "background:#fafafa;padding:.8em 1em;border:1px solid #ddd;" \
           "border-radius:4px;overflow-x:auto'>" + "\n".join(lines) + "</pre>"


def _render_parent_map_html(p: dict, fasta_rel: Path | None,
                            structures_rel: Path | None) -> str:
    acc = p["uniprot_acc"]
    name = p["name"] or ""
    organism = p["organism"] or ""
    seq = p["sequence"]
    # Sort epitopes by first occurrence position for deterministic colour
    # assignment.
    eps = sorted(
        p["epitopes"],
        key=lambda e: (e["occurrences"][0][0] if e["occurrences"] else 0, e["epitope"]))
    ranges: list[tuple[int, int, int]] = []
    for idx, e in enumerate(eps):
        for s, end in e["occurrences"]:
            ranges.append((s, end, idx))

    legend_items: list[str] = []
    for idx, e in enumerate(eps):
        color = EPITOPE_COLORS[idx % len(EPITOPE_COLORS)]
        positions = ", ".join(f"{s}-{end}" for s, end in e["occurrences"])
        legend_items.append(
            f"<li><span style='background:{color};padding:.1em .5em;"
            f"border-radius:3px;font-family:monospace'>"
            f"{escape(e['epitope'])}</span> &nbsp; "
            f"<span class='muted'>{len(e['epitope'])} aa</span> &nbsp; at "
            f"<b>{positions}</b></li>")
    block_html = _render_sequence_block(seq, ranges)

    links: list[str] = [
        f"<a href='https://www.uniprot.org/uniprotkb/{escape(acc)}'>UniProt page</a>"
    ]
    if fasta_rel is not None:
        # FASTA lives next to this HTML inside parents/ — strip the dir.
        links.append(f"<a href='{escape(Path(fasta_rel).name)}'>FASTA</a>")
    if structures_rel is not None:
        links.append(f"<a href='../{escape(str(structures_rel))}'>structures CSV</a>")
    links.append("<a href='../report.html'>back to report</a>")

    css = """
      body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
           max-width:980px;margin:2em auto;padding:0 1em;color:#222;line-height:1.5}
      h1{border-bottom:2px solid #444;padding-bottom:.3em}
      h2{margin-top:1.6em;border-bottom:1px solid #ddd;padding-bottom:.2em}
      .muted{color:#888;font-size:13px}
      .meta li{margin:.2em 0}
      ul.legend{padding-left:1.2em} ul.legend li{margin:.25em 0}
    """
    parts = [
        f"<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>{escape(acc)} – epitope map</title>"
        f"<style>{css}</style></head><body>",
        f"<h1>{escape(acc)} – {escape(name)}</h1>",
        "<ul class='meta'>"
        f"<li><b>Organism:</b> {escape(organism)}</li>"
        f"<li><b>Length:</b> {p['length']} aa</li>"
        f"<li><b>Links:</b> {' &middot; '.join(links)}</li>"
        "</ul>",
        f"<h2>Mapped epitopes ({len(eps)})</h2>",
        "<ul class='legend'>" + "".join(legend_items) + "</ul>",
        "<h2>Sequence</h2>",
        block_html,
        "</body></html>",
    ]
    return "\n".join(parts)


def write_parent_maps(out_dir: Path, parent_index: dict[str, dict],
                      fasta_map: dict[str, Path],
                      structures_map: dict[str, Path],
                      log=sys.stderr) -> dict[str, Path]:
    sub = out_dir / "parents"
    sub.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for acc, p in parent_index.items():
        html = _render_parent_map_html(
            p, fasta_map.get(acc), structures_map.get(acc))
        path = sub / f"{acc}.html"
        path.write_text(html)
        written[acc] = path.relative_to(out_dir)
        print(f"wrote {path}", file=log)
    return written


# ---------- CLI ----------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Map epitope peptides to parent proteins and structures.")
    ap.add_argument("peptides", nargs="*", help="Epitope peptide(s).")
    ap.add_argument("--input", "-i", help="File with one peptide per line, or FASTA.")
    ap.add_argument("--out-dir", "-o", default="epitope_report",
                    help="Output directory (default: ./epitope_report). "
                         "Will contain report.md, report.html, epitopes.csv, "
                         "structures/<UNIPROT>.csv, and report.json.")
    ap.add_argument("--md-only", action="store_true", help="Skip HTML output.")
    ap.add_argument("--html-only", action="store_true", help="Skip Markdown output.")
    ap.add_argument("--no-alphafold", action="store_true",
                    help="Don't query AlphaFold even if no experimental PDB.")
    ap.add_argument("--no-json", action="store_true",
                    help="Skip the raw report.json dump.")
    ap.add_argument("--workers", "-j", type=int, default=8,
                    help="Number of concurrent worker threads for the outer "
                         "peptide loop (default: 8). Lower to 1–4 if IEDB or "
                         "UniProt start returning HTTP 429.")
    args = ap.parse_args()

    peptides = parse_inputs(args.peptides, args.input)
    if not peptides:
        ap.error("No epitopes given. Pass one as an argument or via --input.")

    def _safe(seq: str) -> dict:
        try:
            return process_epitope(seq, use_alphafold=not args.no_alphafold)
        except Exception as e:  # noqa: BLE001
            _log(f"[error] {seq}: {e}")
            return {"epitope": seq, "status": "error", "error": str(e)}

    n_workers = max(1, args.workers)
    if n_workers == 1 or len(peptides) <= 1:
        results: list[dict] = [_safe(seq) for seq in peptides]
    else:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            # submit-then-collect preserves input order
            futures = [pool.submit(_safe, seq) for seq in peptides]
            results = [f.result() for f in futures]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = _summary_rows(results)
    parent_csv_map = write_parent_csvs(out_dir, results)

    master_csv = out_dir / "epitopes.csv"
    write_master_csv(master_csv, summary_rows)
    print(f"wrote {master_csv}", file=sys.stderr)

    for acc, rel in parent_csv_map.items():
        print(f"wrote {out_dir / rel}", file=sys.stderr)

    parent_index = build_parent_index(results)
    fasta_map = write_parent_fastas(out_dir, parent_index)
    map_map = write_parent_maps(out_dir, parent_index, fasta_map, parent_csv_map)

    if not args.html_only:
        md_path = out_dir / "report.md"
        md_path.write_text(render_markdown(
            results, summary_rows, parent_csv_map, fasta_map, map_map))
        print(f"wrote {md_path}", file=sys.stderr)
    if not args.md_only:
        html_path = out_dir / "report.html"
        html_path.write_text(render_html(
            results, summary_rows, parent_csv_map, fasta_map, map_map))
        print(f"wrote {html_path}", file=sys.stderr)
    if not args.no_json:
        json_path = out_dir / "report.json"
        json_path.write_text(json.dumps(results, indent=2, default=str))
        print(f"wrote {json_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
