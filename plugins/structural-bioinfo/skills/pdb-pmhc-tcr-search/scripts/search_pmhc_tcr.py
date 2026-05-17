#!/usr/bin/env python3
"""Search RCSB for pMHC and TCR-pMHC structures by epitope peptide sequence.

For each input epitope:
  1. Search the RCSB Search API for entries containing both (a) a polymer chain whose
     canonical sequence equals or contains the epitope and (b) an MHC chain (matched by
     UniProt accession).
  2. Pull metadata via the RCSB Data API GraphQL endpoint.
  3. Classify each entry's polymer chains as peptide / MHC heavy / β2m or class-II β /
     TCR α / TCR β by UniProt accession or description regex.
  4. Tag each entry as `pmhc` or `tcr_pmhc`; group redundant copies (same epitope, same
     MHC alleles, same TCR if present); mark the lowest-resolution copy as representative.
  5. Download every entry's PDB (CIF fallback handled by the cached helper).
  6. Split TCR-pMHC entries into derivative `<pdb>_pmhc.pdb` and `<pdb>_tcr.pdb` files.
  7. Emit `report.md`, `report.html`, and `report_data.json` in the output directory.

Stdlib-only. Reuses the HTTP helpers, download helper, and CIF-to-PDB conversion patterns
from the sibling `pdb-holostructure-search` and `pdb-extractor` skills.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import html
import json
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

# Local import (same skill, scripts/ dir)
_SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPT_DIR))
from split_tcr_pmhc import split_entry, EXCLUDE_RESIDUES  # noqa: E402


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/graphql"
RCSB_PDB_URL_FMT = "https://files.rcsb.org/download/{pdb}.pdb"
RCSB_CIF_URL_FMT = "https://files.rcsb.org/download/{pdb}.cif"

USER_AGENT = "pdb-pmhc-tcr-search/1.0"

AA = "ACDEFGHIKLMNPQRSTVWY"
EPITOPE_RE = re.compile(f"^[{AA}]+$")

# MHC Class I heavy chains (human + mouse)
MHC_CLASS_I_UNIPROTS = {
    "P04439": "HLA-A",  "P01889": "HLA-B",  "P10321": "HLA-C",
    "P13747": "HLA-E",  "P30511": "HLA-F",  "P17693": "HLA-G",
    "P01901": "H-2K",   "P01899": "H-2D",   "P14430": "H-2L",
}
# β2-microglobulin
B2M_UNIPROTS = {"P61769": "β2m (human)", "P01887": "β2m (mouse)"}
# MHC Class II α
MHC_CLASS_II_ALPHA_UNIPROTS = {
    "P01903": "HLA-DRA",  "P20036": "HLA-DPA1", "P01909": "HLA-DQA1",
    "P14434": "H2-Aa",    "P14437": "H2-Ea",
}
# MHC Class II β
MHC_CLASS_II_BETA_UNIPROTS = {
    "P04229": "HLA-DRB1", "P04440": "HLA-DPB1", "P01920": "HLA-DQB1",
    "P14483": "H2-Ab1",   "P04230": "H2-Eb1",
}
# Union: any of these on a polymer entity means "this entry contains MHC machinery"
MHC_UNIPROTS = {
    **MHC_CLASS_I_UNIPROTS,
    **B2M_UNIPROTS,
    **MHC_CLASS_II_ALPHA_UNIPROTS,
    **MHC_CLASS_II_BETA_UNIPROTS,
}

# TCR constant-region UniProts (human + mouse)
TCR_ALPHA_UNIPROTS = {"P01848", "P01849"}   # TRAC human, mouse
TCR_BETA_UNIPROTS = {"P01850", "P01852", "A0A075B6X1"}  # TRBC1 human, mouse; TRBC2 human

# Description-based fallback regexes (case-insensitive). See references/tcr_keywords.md.
_TCR_HINT = re.compile(r"t[\- ]?cell\s*receptor|\btcr\b|\btra[vjc]\b|\btrb[vdjc]\b", re.I)
_TCR_ALPHA_RE = re.compile(
    r"t[\- ]?cell\s*receptor[^.]*alpha|\btcr[^.]*alpha|\btcr[\-\s]?α|\btra[vjc]\b",
    re.I,
)
_TCR_BETA_RE = re.compile(
    r"t[\- ]?cell\s*receptor[^.]*beta|\btcr[^.]*beta|\btcr[\-\s]?β|\btrb[vdjc]\b",
    re.I,
)
# Class I heavy chain hints — order matters: we test these AFTER the Class-II hints in
# `_structural_role` so a "histocompatibility" match on a DRA/DRB chain doesn't capture
# the entity into the Class-I bucket. The patterns are deliberately Class-I-specific
# (HLA-A/B/C/E/F/G, H-2K/D/L/M; explicit "class I" prefix).
_MHC_HEAVY_HINT = re.compile(
    r"hla[\- ][abcefg]\b|h[\- ]2[kdlm]\b|class\s*i\s*histocompatibility|"
    r"mhc\s*class\s*i\b|h2[\- ]?[kdlm]\b",
    re.I,
)
_MHC_LIGHT_HINT = re.compile(r"beta[\-\s]?2[\-\s]?microglobulin|\bb2m\b", re.I)
_MHC_II_BETA_HINT = re.compile(
    r"hla[\- ]?d[pqr]b|class\s*ii\s*histocompatibility.*\bbeta\b|"
    r"class\s*ii.*beta|\bdrb\d|\bdpb\d|\bdqb\d|h2[\- ]?[ae]b",
    re.I,
)
_MHC_II_ALPHA_HINT = re.compile(
    r"hla[\- ]?d[pqr]a|class\s*ii\s*histocompatibility.*\balpha\b|"
    r"class\s*ii.*alpha|\bdra\b|\bdpa\d|\bdqa\d|h2[\- ]?[ae]a",
    re.I,
)
# Generic "histocompatibility" fallback (no class signal). Used only AFTER class-II
# hints have already been ruled out — defaults to Class I (the historical PDB default
# for unannotated HLA chains).
_HISTOCOMP_GENERIC = re.compile(r"histocompatibility", re.I)


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

def http_post_json(url: str, payload: dict, timeout: int = 60) -> dict:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8")
        if not body.strip():
            return {"result_set": [], "total_count": 0}
        return json.loads(body)


# ---------------------------------------------------------------------------
# Epitope input
# ---------------------------------------------------------------------------

def normalize_epitope(raw: str) -> str:
    seq = raw.strip().upper().replace(" ", "")
    if not seq:
        raise ValueError("empty epitope")
    if not EPITOPE_RE.match(seq):
        bad = sorted(set(seq) - set(AA))
        raise ValueError(f"epitope '{raw}' contains non-canonical residues: {bad}")
    if not (7 <= len(seq) <= 30):
        raise ValueError(f"epitope '{raw}' length {len(seq)} outside expected 7–30 range "
                         f"(Class I 8–11, Class II 13–25)")
    return seq


def collect_epitopes(args: argparse.Namespace) -> list[tuple[str, str]]:
    """Return list of (sequence, label) tuples."""
    out: list[tuple[str, str]] = []
    seen: set[str] = set()

    def _add(seq: str, label: str | None = None) -> None:
        s = normalize_epitope(seq)
        if s in seen:
            return
        seen.add(s)
        out.append((s, label or s))

    if args.epitope:
        _add(args.epitope)
    if args.epitopes:
        for chunk in args.epitopes.split(","):
            if chunk.strip():
                _add(chunk)
    if args.epitope_file:
        with open(args.epitope_file, encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                parts = ln.split(None, 1)
                seq = parts[0]
                label = parts[1] if len(parts) > 1 else None
                _add(seq, label)

    if not out:
        raise SystemExit("[error] no epitopes supplied. Use --epitope, --epitopes, or --epitope-file.")
    return out


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def build_query(epitope: str) -> dict:
    """Build a Search API query for an epitope.

    The dedicated `seqmotif` service is the only RCSB Search service that supports short
    peptide substring lookups (the `text` service does not index the canonical sequence
    field, and the `sequence` service requires a 25-aa minimum). We issue seqmotif alone
    and rely on the metadata classifier to reject hits that lack an MHC chain — this is
    a small amount of extra metadata fetching but avoids RCSB's quirky restrictions on
    which attributes are searchable in `text` terminals.

    `match_type` (free vs. fused) is determined post-hoc from the metadata.
    """
    return {
        "return_type": "entry",
        "query": {
            "type": "terminal", "service": "seqmotif",
            "parameters": {
                "value": epitope,
                "pattern_type": "simple",
                "sequence_type": "protein",
            },
        },
        "request_options": {
            "results_content_type": ["experimental"],
            "return_all_hits": True,
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }


def search_epitope(epitope: str) -> list[str]:
    """Run the seqmotif query and return a list of unique PDB IDs (unfiltered for MHC)."""
    sys.stderr.write(f"[info] {epitope}: searching RCSB (seqmotif)…\n")
    try:
        resp = http_post_json(SEARCH_URL, build_query(epitope))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")[:500]
        sys.stderr.write(f"[error] {epitope}: HTTP {exc.code}: {body}\n")
        return []
    except Exception as exc:
        sys.stderr.write(f"[error] {epitope}: {exc}\n")
        return []
    hits: list[str] = []
    seen: set[str] = set()
    for h in resp.get("result_set") or []:
        pid = (h.get("identifier") or "").upper()
        if pid and pid not in seen:
            seen.add(pid)
            hits.append(pid)
    sys.stderr.write(f"[info] {epitope}: {len(hits)} unique entries\n")
    return hits


# ---------------------------------------------------------------------------
# Metadata + classification
# ---------------------------------------------------------------------------

DATA_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    struct_keywords { text pdbx_keywords }
    rcsb_entry_info {
      resolution_combined
      experimental_method
    }
    rcsb_accession_info { deposit_date }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        entity_id
        uniprot_ids
        auth_asym_ids
      }
      entity_poly {
        pdbx_seq_one_letter_code_can
        rcsb_sample_sequence_length
        type
      }
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_name_com { name }
    }
  }
}
"""


def fetch_metadata(pdb_ids: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not pdb_ids:
        return out
    for i in range(0, len(pdb_ids), 50):
        batch = [p.upper() for p in pdb_ids[i:i + 50]]
        try:
            resp = http_post_json(DATA_URL, {"query": DATA_QUERY, "variables": {"ids": batch}})
        except Exception as exc:
            sys.stderr.write(f"[warn] Data API batch {batch[0]}..{batch[-1]} failed: {exc}\n")
            continue
        for entry in (resp.get("data") or {}).get("entries") or []:
            pid = (entry.get("rcsb_id") or "").upper()
            if not pid:
                continue
            out[pid] = entry
    return out


def _structural_role(uniprots: set[str], desc: str) -> str | None:
    """Return the protein-identity role of an entity (mhc_heavy / mhc_light / mhc_ii_alpha /
    mhc_ii_beta / tcr_alpha / tcr_beta), or None if not recognised. Peptide-hosting is
    determined separately by the caller."""
    # TCR chains — UniProt first, then regex.
    if uniprots & TCR_ALPHA_UNIPROTS:
        return "tcr_alpha"
    if uniprots & TCR_BETA_UNIPROTS:
        return "tcr_beta"
    if _TCR_HINT.search(desc):
        if _TCR_ALPHA_RE.search(desc):
            return "tcr_alpha"
        if _TCR_BETA_RE.search(desc):
            return "tcr_beta"

    # MHC chains — UniProt first.
    if uniprots & set(MHC_CLASS_I_UNIPROTS):
        return "mhc_heavy"
    if uniprots & set(B2M_UNIPROTS):
        return "mhc_light"
    if uniprots & set(MHC_CLASS_II_ALPHA_UNIPROTS):
        return "mhc_ii_alpha"
    if uniprots & set(MHC_CLASS_II_BETA_UNIPROTS):
        return "mhc_ii_beta"

    # MHC chains — description regex fallback. Check Class II first (more specific
    # patterns include "d[pqr][ab]"); fall back to Class I; only then accept a generic
    # "histocompatibility" hit and assume Class I (the historical PDB default).
    if _MHC_II_BETA_HINT.search(desc):
        return "mhc_ii_beta"
    if _MHC_II_ALPHA_HINT.search(desc):
        return "mhc_ii_alpha"
    if _MHC_LIGHT_HINT.search(desc):
        return "mhc_light"
    if _MHC_HEAVY_HINT.search(desc):
        return "mhc_heavy"
    if _HISTOCOMP_GENERIC.search(desc):
        return "mhc_heavy"

    return None


def _hosts_peptide(epitope: str, seq: str, struct_role: str | None) -> bool:
    """Does this entity carry the epitope? True if (a) entire chain == epitope, or
    (b) chain contains the epitope as substring (single-chain-trimer / fusion construct).
    To avoid mis-tagging a long unrelated protein that happens to embed the epitope by
    chance (e.g. full ovalbumin contains SIINFEKL), substring matches are only accepted
    when the entity also has a recognised MHC or TCR role — i.e. it's part of the
    pMHC/TCR machinery and the embedding is intentional."""
    if not seq:
        return False
    if seq == epitope:
        return True
    if epitope in seq and struct_role is not None:
        return True
    return False


def classify_entry(pdb_id: str, epitope: str, entry: dict) -> dict | None:
    """Convert one Data-API entry into the structured record this skill works with.

    Returns None and logs a [warn] if no peptide chain matching the epitope is found, or
    if no MHC chain is found.
    """
    chains_by_role: dict[str, list[str]] = {
        "peptide": [], "mhc_heavy": [], "mhc_light": [],
        "mhc_ii_alpha": [], "mhc_ii_beta": [],
        "tcr_alpha": [], "tcr_beta": [],
    }
    mhc_uniprots_seen: set[str] = set()
    tcr_uniprots_seen: set[str] = set()
    descriptions: list[str] = []
    mhc_descriptions: list[str] = []
    peptide_seq_seen: str | None = None

    for pe in entry.get("polymer_entities") or []:
        cid = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        uniprots = set(cid.get("uniprot_ids") or [])
        # auth_asym_ids can be multi-character in mmCIF (e.g. "DDD"); the PDB file format
        # only has one chain-ID column, and the CIF→PDB converter takes the first char.
        # Apply the same transformation so the splitter matches what's in the .pdb file.
        auth_chains = [c[0] for c in (cid.get("auth_asym_ids") or []) if c]
        ep = pe.get("entity_poly") or {}
        seq = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").upper()
        desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "")
        name_coms = [n.get("name", "") for n in (pe.get("rcsb_polymer_entity_name_com") or [])]
        if desc:
            descriptions.append(desc)

        struct_role = _structural_role(uniprots, desc)
        if struct_role is not None:
            chains_by_role[struct_role].extend(auth_chains)
            if struct_role in ("mhc_heavy", "mhc_light", "mhc_ii_alpha", "mhc_ii_beta"):
                mhc_uniprots_seen |= uniprots
                if struct_role != "mhc_light":
                    for t in (desc, *name_coms):
                        if t:
                            mhc_descriptions.append(t)
            if struct_role in ("tcr_alpha", "tcr_beta"):
                tcr_uniprots_seen |= uniprots & (TCR_ALPHA_UNIPROTS | TCR_BETA_UNIPROTS)

        # Peptide hosting is independent of structural role: a single-chain trimer chain
        # can be BOTH the MHC heavy and the peptide host.
        if _hosts_peptide(epitope, seq, struct_role):
            chains_by_role["peptide"].extend(auth_chains)
            if not peptide_seq_seen or len(seq) < len(peptide_seq_seen):
                peptide_seq_seen = seq

    has_peptide = bool(chains_by_role["peptide"])
    has_mhc = bool(chains_by_role["mhc_heavy"] or chains_by_role["mhc_ii_alpha"]
                   or chains_by_role["mhc_ii_beta"])
    has_tcr = bool(chains_by_role["tcr_alpha"] or chains_by_role["tcr_beta"])

    if not (has_peptide and has_mhc):
        reasons = []
        if not has_peptide:
            reasons.append("no chain carrying the epitope")
        if not has_mhc:
            reasons.append("no MHC chain found")
        return {
            "pdb_id": pdb_id,
            "dropped": True,
            "reason": "; ".join(reasons),
            "descriptions": descriptions,
        }

    # MHC class: Class I if heavy chain present; else Class II
    if chains_by_role["mhc_heavy"]:
        mhc_class = "I"
    elif chains_by_role["mhc_ii_alpha"] or chains_by_role["mhc_ii_beta"]:
        mhc_class = "II"
    else:
        # length-based fallback (8–11 → I, 13–25 → II)
        mhc_class = "I" if len(epitope) <= 11 else "II"

    light_chains = chains_by_role["mhc_light"] + chains_by_role["mhc_ii_beta"]

    res_list = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    resolution = res_list[0] if res_list else None
    method = (entry.get("rcsb_entry_info") or {}).get("experimental_method") or ""
    title = (entry.get("struct") or {}).get("title") or ""
    kw = entry.get("struct_keywords") or {}
    keywords_text = " ".join(filter(None, [kw.get("pdbx_keywords"), kw.get("text")]))
    deposit_date = (entry.get("rcsb_accession_info") or {}).get("deposit_date") or ""
    if title:
        mhc_descriptions.append(title)
    if keywords_text:
        mhc_descriptions.append(keywords_text)

    return {
        "pdb_id": pdb_id,
        "dropped": False,
        "complex_type": "tcr_pmhc" if has_tcr else "pmhc",
        "mhc_class": mhc_class,
        "resolution": resolution,
        "method": method,
        "title": title,
        "deposit_date": deposit_date,
        "peptide_chains": sorted(set(chains_by_role["peptide"])),
        "mhc_heavy_chains": sorted(set(chains_by_role["mhc_heavy"])),
        "mhc_ii_alpha_chains": sorted(set(chains_by_role["mhc_ii_alpha"])),
        "mhc_ii_beta_chains": sorted(set(chains_by_role["mhc_ii_beta"])),
        "light_chains": sorted(set(light_chains)),
        "tcr_alpha_chains": sorted(set(chains_by_role["tcr_alpha"])),
        "tcr_beta_chains": sorted(set(chains_by_role["tcr_beta"])),
        "mhc_uniprots": sorted(mhc_uniprots_seen),
        "tcr_uniprots": sorted(tcr_uniprots_seen),
        "mhc_descriptions": mhc_descriptions,
        "peptide_seq": peptide_seq_seen,
    }


# ---------------------------------------------------------------------------
# Group redundants
# ---------------------------------------------------------------------------

def assign_groups(records: list[dict]) -> None:
    """Group by (mhc_uniprots, tcr_uniprots); mark lowest-resolution as representative."""
    keys: dict[tuple, list[dict]] = {}
    for r in records:
        if r.get("dropped"):
            continue
        key = (tuple(r["mhc_uniprots"]), tuple(r["tcr_uniprots"]))
        keys.setdefault(key, []).append(r)

    group_id = 0
    for key in sorted(keys.keys(), key=lambda k: (-len(keys[k]),)):
        group_id += 1
        gid = f"g{group_id:02d}"
        members = keys[key]
        # Best: lowest resolution; tie-break by earliest deposit date.
        members.sort(key=lambda r: (
            r["resolution"] if isinstance(r["resolution"], (int, float)) else 999.0,
            r["deposit_date"],
        ))
        for i, r in enumerate(members):
            r["group_id"] = gid
            r["representative"] = (i == 0)
            r["group_size"] = len(members)


# ---------------------------------------------------------------------------
# Download (replicates pdb-extractor.download_coords)
# ---------------------------------------------------------------------------

def _cif_to_pdb_lines(cif_path: str):
    """Minimal mmCIF atom_site → PDB ATOM/HETATM converter. Mirrors pdb-extractor."""
    with open(cif_path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    col_map: dict[str, int] = {}
    atom_lines: list[str] = []
    in_loop = False
    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")
        if line.startswith("_atom_site.") and not in_loop:
            col_names: list[str] = []
            j = i
            while j < len(lines) and lines[j].startswith("_atom_site."):
                col_names.append(lines[j].split(".")[-1].strip())
                j += 1
            col_map = {n: idx for idx, n in enumerate(col_names)}
            in_loop = True
            i = j
            continue
        if in_loop:
            stripped = line.strip()
            if not stripped or stripped.startswith("_") or stripped.startswith("loop_"):
                break
            if not stripped.startswith("#"):
                atom_lines.append(line)
        i += 1

    def col(name: str) -> int | None:
        return col_map.get(name)

    idxs = {k: col(k) for k in (
        "id", "group_PDB", "type_symbol", "label_alt_id", "pdbx_PDB_ins_code",
        "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
        "auth_atom_id", "label_atom_id", "auth_comp_id", "auth_asym_id", "auth_seq_id",
    )}
    n_cols = max([v for v in col_map.values()]) + 1 if col_map else 0

    for raw in atom_lines:
        fields = raw.split()
        while len(fields) < n_cols:
            fields.append("?")

        def get(k: str, default: str = "") -> str:
            ix = idxs.get(k)
            if ix is None or ix >= len(fields):
                return default
            v = fields[ix]
            return v if v not in (".", "?") else ""

        x, y, z = get("Cartn_x"), get("Cartn_y"), get("Cartn_z")
        if not (x and y and z):
            continue
        group = get("group_PDB", "ATOM") or "ATOM"
        atm_id = get("id", "")
        atm_name = get("auth_atom_id") or get("label_atom_id") or ""
        alt_loc = get("label_alt_id", "")
        res_name = (get("auth_comp_id", "") or "")[:3]
        chain = (get("auth_asym_id", "") or " ")[0]
        res_seq = get("auth_seq_id", "")
        ins_code = (get("pdbx_PDB_ins_code", "") or " ")[0]
        occ = get("occupancy", "1.00") or "1.00"
        b = get("B_iso_or_equiv", "0.00") or "0.00"
        elem = get("type_symbol", "")
        if len(atm_name) < 4:
            atm_name = f"{atm_name:<4}" if len(elem) == 2 else f" {atm_name:<3}"
        atm_name = atm_name[:4]
        try:
            serial = int(float(atm_id)) % 100000 if atm_id else 0
        except (ValueError, TypeError):
            serial = 0
        try:
            res_seq_i = int(float(res_seq))
        except (ValueError, TypeError):
            res_seq_i = 0
        try:
            xf, yf, zf = float(x), float(y), float(z)
            occf = float(occ); bf = float(b)
        except ValueError:
            continue
        pdb_line = (
            f"{group[:6]:<6}{serial:>5} {atm_name:<4}{(alt_loc or ' ')[0]}{res_name:>3}"
            f" {chain}{res_seq_i:>4}{ins_code}   {xf:>8.3f}{yf:>8.3f}{zf:>8.3f}"
            f"{occf:>6.2f}{bf:>6.2f}          {elem:>2}"
        )
        yield pdb_line[:79].ljust(79) + "\n"


def download_coords(pdb_id: str, output_dir: str) -> str | None:
    pdb_id = pdb_id.upper()
    entry_dir = os.path.join(output_dir, pdb_id)
    os.makedirs(entry_dir, exist_ok=True)
    pdb_path = os.path.join(entry_dir, f"{pdb_id}.pdb")

    if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 0:
        return pdb_path

    # Try legacy PDB first.
    try:
        urllib.request.urlretrieve(RCSB_PDB_URL_FMT.format(pdb=pdb_id), pdb_path)
        if os.path.getsize(pdb_path) > 0:
            return pdb_path
    except Exception:
        pass

    # Fallback: CIF + on-the-fly conversion.
    cif_path = os.path.join(entry_dir, f"{pdb_id}.cif")
    try:
        urllib.request.urlretrieve(RCSB_CIF_URL_FMT.format(pdb=pdb_id), cif_path)
        if os.path.getsize(cif_path) == 0:
            raise RuntimeError("empty CIF")
    except Exception as exc:
        sys.stderr.write(f"[warn] {pdb_id}: download failed: {exc}\n")
        return None

    try:
        with open(pdb_path, "w", encoding="utf-8") as fh:
            for line in _cif_to_pdb_lines(cif_path):
                fh.write(line)
            fh.write("END\n")
    except Exception as exc:
        sys.stderr.write(f"[warn] {pdb_id}: CIF→PDB conversion failed: {exc}\n")
        return None

    return pdb_path


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

# HLA gene letter pattern: classical class I (A/B/C/E/F/G) or class II (DRA, DRB1, DPA1,
# DPB1, DQA1, DQB1). The `w?` after the letter handles historical notation like 'Cw6'.
_GENE_PAT = r"(?:[ABCEFG]|D[PQR][AB]\d?)"

# Pattern 1 — HLA-anchored. Handles "HLA-A*02:01" / "HLA-A 0201" / "HLA-A0201" /
# "HLA-A2.1" / "HLA-A2" / "HLA-A201" (legacy 3-digit shorthand) — the last is parsed as
# family-1-digit + subtype-2-digit, NOT family-2-digit + subtype-1-digit, since the
# legacy shorthand means A*02:01 rather than A*20:01.
_P_HLA = re.compile(
    rf"\bHLA[-_ ]?(?P<gene>{_GENE_PAT})w?"
    r"(?:\*|[ ])?"
    r"(?:"
        r"(?P<f1>\d{1,2})[.:\-](?P<s1>\d{1,3})"     # explicit separator
        r"|(?P<f2>\d{2})(?P<s2>\d{2,3})"             # 4-5 contiguous = 2+2/3
        r"|(?P<f4>\d{1})(?P<s4>\d{2})"               # 3 contiguous = 1+2 (legacy)
        r"|(?P<f3>\d{1,2})"                          # 1-2 digits, no subtype
    r")"
    r"\b",
    re.I,
)

# Pattern 2 — chain text: "A-2 alpha chain" / "B-7 ALPHA CHAIN" / "Cw-6 alpha chain"
_P_CHAIN = re.compile(
    rf"\b(?P<gene>[ABCEFG])w?-(?P<f>\d{{1,2}})\s+(?:alpha|heavy)\s+chain",
    re.I,
)

# Pattern 3 — name_com "antigen A*2" / "MHC class I antigen, B-7"
_P_ANTIGEN = re.compile(
    rf"\bantigen[,\s]+(?P<gene>[ABCEFG])w?[\*\- ](?P<f>\d{{1,2}})"
    r"(?:[.:\- ]?(?P<s>\d{1,3}))?",
    re.I,
)

# Pattern 4 — "MHC (class) I" / "MHC (class) II" anchor without HLA prefix:
# "MHC I A02", "MHC class I, A2", "MHC II DRB1*04:01" — common in titles/keywords.
_P_MHC_ANCHOR = re.compile(
    r"\bMHC(?:\s+class)?\s+I{1,2}\b[\s,:\-]*"
    rf"(?P<gene>{_GENE_PAT})w?"
    r"(?:\*|[ ])?"
    r"(?:"
        r"(?P<f1>\d{1,2})[.:\-](?P<s1>\d{1,3})"
        r"|(?P<f2>\d{2})(?P<s2>\d{2,3})"
        r"|(?P<f4>\d{1})(?P<s4>\d{2})"
        r"|(?P<f3>\d{1,2})"
    r")\b",
    re.I,
)


def _format_allele(gene: str, f: str | None, s: str | None) -> str:
    out = f"HLA-{gene.upper()}"
    if f is None:
        return out
    out = f"{out}*{f.zfill(2)}"
    if s is None:
        return out
    return f"{out}:{s.zfill(2)}"


def _pick_fs(m) -> tuple[str | None, str | None]:
    """Pull out (family, subtype) from whichever alternation arm matched in _P_HLA."""
    for tag in ("1", "2", "4", "3"):
        f = m.group("f" + tag)
        if f is not None:
            s = m.groupdict().get("s" + tag)
            return f, s
    return None, None


def _alleles_in_text(text: str) -> list[str]:
    """Find all HLA alleles in text; return an ordered, deduplicated list. When a 4-digit
    allele subsumes a 2-digit one (e.g. HLA-A*02:01 vs HLA-A*02 from different fields of
    the same entry), the less-specific form is dropped."""
    if not text:
        return []
    found: dict[str, None] = {}
    for m in _P_HLA.finditer(text):
        f, s = _pick_fs(m)
        found[_format_allele(m.group("gene"), f, s)] = None
    for m in _P_CHAIN.finditer(text):
        found[_format_allele(m.group("gene"), m.group("f"), None)] = None
    for m in _P_ANTIGEN.finditer(text):
        found[_format_allele(m.group("gene"), m.group("f"), m.group("s"))] = None
    for m in _P_MHC_ANCHOR.finditer(text):
        f, s = _pick_fs(m)
        found[_format_allele(m.group("gene"), f, s)] = None
    keys = list(found)
    # Drop strict prefixes — keep the most specific form for each gene+family.
    return [a for a in keys if not any(b != a and b.startswith(a + ":") for b in keys)]


def _allele_str(record: dict) -> str:
    """Format the MHC column. Shows the parsed full allele(s); if not parseable, falls
    back to the gene symbol with a '(allele not specified)' note. β2m is omitted unless
    it's a non-default species (anything other than human P61769)."""
    alleles = _alleles_in_text(" | ".join(record.get("mhc_descriptions", [])))
    if alleles:
        main = "; ".join(alleles)
    else:
        genes = sorted({
            MHC_UNIPROTS[u] for u in record.get("mhc_uniprots", [])
            if u in MHC_UNIPROTS and u not in B2M_UNIPROTS
        })
        main = f"{', '.join(genes)} (allele not specified)" if genes \
               else "(allele not specified)"
    # Surface non-default β2m as a suffix note.
    b2m = sorted({
        B2M_UNIPROTS[u] for u in record.get("mhc_uniprots", [])
        if u in B2M_UNIPROTS and u != "P61769"
    })
    if b2m:
        main += f" · {', '.join(b2m)}"
    return main


def _tcr_str(record: dict) -> str:
    if record["complex_type"] != "tcr_pmhc":
        return "—"
    a = record["tcr_alpha_chains"]
    b = record["tcr_beta_chains"]
    pieces = []
    if a: pieces.append(f"α:{','.join(a)}")
    if b: pieces.append(f"β:{','.join(b)}")
    return " ".join(pieces) if pieces else "TCR"


def _files_str(record: dict, root: str) -> str:
    parts: list[str] = []
    for label, key in (("full", "pdb_path"), ("pmhc", "pmhc_path"), ("tcr", "tcr_path")):
        p = record.get(key)
        if not p:
            continue
        rel = os.path.relpath(p, root)
        parts.append(f"[{label}]({rel})")
    return " ".join(parts) if parts else "—"


def write_markdown_report(records_by_epitope: dict, output_dir: str,
                          epitopes: list[tuple[str, str]]) -> str:
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    out_path = os.path.join(output_dir, "report.md")

    lines: list[str] = []
    lines.append("# pMHC / TCR-pMHC PDB search report")
    lines.append("")
    lines.append(f"_Run: {now} · RCSB Search API v2 · {len(epitopes)} epitope(s)_")
    lines.append("")
    lines.append("## Input epitopes")
    lines.append("")
    lines.append("| label | sequence | length | class (heuristic) |")
    lines.append("|---|---|---|---|")
    for seq, label in epitopes:
        cls = "I" if len(seq) <= 11 else "II"
        lines.append(f"| {label} | `{seq}` | {len(seq)} | {cls} |")
    lines.append("")

    for seq, label in epitopes:
        records = records_by_epitope.get(seq, [])
        kept = [r for r in records if not r.get("dropped")]
        n_pmhc = sum(1 for r in kept if r["complex_type"] == "pmhc")
        n_tcr = sum(1 for r in kept if r["complex_type"] == "tcr_pmhc")
        groups = {r["group_id"] for r in kept}
        lines.append(f"## {label}  (`{seq}`)")
        lines.append("")
        lines.append(f"**Summary:** {len(kept)} entries · {n_pmhc} pMHC · {n_tcr} TCR-pMHC · "
                     f"{len(groups)} redundancy group(s)")
        lines.append("")
        if not kept:
            lines.append("_No matching entries returned._")
            lines.append("")
            continue

        lines.append("| group | rep | PDB | type | class | resolution (Å) | method | "
                     "match | MHC | TCR | files |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        kept.sort(key=lambda r: (r["group_id"], not r["representative"]))
        for r in kept:
            res = f"{r['resolution']:.2f}" if isinstance(r["resolution"], (int, float)) else "—"
            pdb_link = f"[{r['pdb_id']}](https://www.rcsb.org/structure/{r['pdb_id']})"
            rep_marker = "★" if r["representative"] else ""
            lines.append(
                f"| {r['group_id']} | {rep_marker} | {pdb_link} | "
                f"`{r['complex_type']}` | {r['mhc_class']} | {res} | "
                f"{r['method']} | {r.get('match_type', '?')} | "
                f"{_allele_str(r)} | {_tcr_str(r)} | "
                f"{_files_str(r, output_dir)} |"
            )
        lines.append("")

    Path(out_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


def _html_row(r: dict, root: str) -> str:
    res = f"{r['resolution']:.2f}" if isinstance(r["resolution"], (int, float)) else "—"
    pdb_link = (f'<a href="https://www.rcsb.org/structure/{r["pdb_id"]}" target="_blank">'
                f'{r["pdb_id"]}</a> '
                f'(<a href="https://www.rcsb.org/3d-view/{r["pdb_id"]}" target="_blank">3D</a>)')
    type_tag = (f'<span class="tag tag-tcr">TCR-pMHC</span>'
                if r["complex_type"] == "tcr_pmhc" else
                f'<span class="tag tag-pmhc">pMHC</span>')
    class_tag = (f'<span class="tag tag-classI">Class I</span>' if r["mhc_class"] == "I"
                 else f'<span class="tag tag-classII">Class II</span>')
    rep_tag = ('<span class="tag tag-rep">representative</span>' if r["representative"]
               else '<span class="tag tag-dup">duplicate</span>')

    files: list[str] = []
    for label, key in (("full", "pdb_path"), ("pmhc", "pmhc_path"), ("tcr", "tcr_path")):
        p = r.get(key)
        if not p:
            continue
        rel = os.path.relpath(p, root)
        files.append(f'<a href="{html.escape(rel)}">{label}</a>')
    files_html = " ".join(files) if files else "—"

    return (
        "<tr>"
        f"<td>{html.escape(r['group_id'])}</td>"
        f"<td>{rep_tag}</td>"
        f"<td>{pdb_link}</td>"
        f"<td>{type_tag}</td>"
        f"<td>{class_tag}</td>"
        f"<td>{res}</td>"
        f"<td>{html.escape(r['method'])}</td>"
        f"<td>{html.escape(r.get('match_type', '?'))}</td>"
        f"<td>{html.escape(_allele_str(r))}</td>"
        f"<td>{html.escape(_tcr_str(r))}</td>"
        f"<td class='files'>{files_html}</td>"
        "</tr>"
    )


def write_html_report(records_by_epitope: dict, output_dir: str,
                      epitopes: list[tuple[str, str]]) -> str:
    template_path = _SCRIPT_DIR.parent / "references" / "report_template.html"
    template = template_path.read_text(encoding="utf-8")

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    subtitle = (f"Generated {now} · RCSB Search API v2 · "
                f"{len(epitopes)} epitope(s) · "
                f"{sum(len([r for r in records_by_epitope.get(s, []) if not r.get('dropped')]) for s, _ in epitopes)} entries total")
    title = ", ".join(label for _, label in epitopes)

    body_chunks: list[str] = []
    for seq, label in epitopes:
        records = records_by_epitope.get(seq, [])
        kept = [r for r in records if not r.get("dropped")]
        n_pmhc = sum(1 for r in kept if r["complex_type"] == "pmhc")
        n_tcr = sum(1 for r in kept if r["complex_type"] == "tcr_pmhc")
        groups = {r["group_id"] for r in kept}
        body_chunks.append(f"<h2>{html.escape(label)} "
                           f"&middot; <code>{html.escape(seq)}</code> "
                           f"&middot; {len(seq)} aa</h2>")
        body_chunks.append(
            f"<div class='summary'>"
            f"<div><strong>{len(kept)}</strong> entries</div>"
            f"<div><strong>{n_pmhc}</strong> pMHC</div>"
            f"<div><strong>{n_tcr}</strong> TCR-pMHC</div>"
            f"<div><strong>{len(groups)}</strong> redundancy group(s)</div>"
            f"</div>"
        )
        if not kept:
            body_chunks.append("<p><em>No matching entries returned.</em></p>")
            continue
        kept.sort(key=lambda r: (r["group_id"], not r["representative"]))
        body_chunks.append(
            "<table><thead><tr>"
            "<th>group</th><th>rep</th><th>PDB</th><th>type</th><th>class</th>"
            "<th>resolution Å</th><th>method</th><th>match</th>"
            "<th>MHC</th><th>TCR</th><th>files</th>"
            "</tr></thead><tbody>"
            + "".join(_html_row(r, output_dir) for r in kept)
            + "</tbody></table>"
        )

    caveats_html = ""

    out_html = (template
                .replace("__TITLE__", html.escape(title or "epitope search"))
                .replace("__SUBTITLE__", html.escape(subtitle))
                .replace("__BODY__", "".join(body_chunks))
                .replace("__CAVEATS__", caveats_html))

    out_path = os.path.join(output_dir, "report.html")
    Path(out_path).write_text(out_html, encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group("epitope input (combine freely)")
    g.add_argument("--epitope", help="A single peptide sequence (e.g. SIINFEKL).")
    g.add_argument("--epitopes",
                   help="Comma-separated peptide list (e.g. SIINFEKL,GILGFVFTL,GLCTLVAML).")
    g.add_argument("--epitope-file",
                   help="Path to a text file with one peptide per line. "
                        "Optional 2nd whitespace column is a human-readable label.")
    p.add_argument("--output-dir", default="collect",
                   help="Root output directory (default: collect/). One subdirectory per "
                        "epitope, with per-PDB subdirectories underneath.")
    p.add_argument("--no-download", action="store_true",
                   help="Skip downloading PDB files (still emits the reports).")
    p.add_argument("--no-split", action="store_true",
                   help="Skip writing _pmhc.pdb / _tcr.pdb derivatives.")
    p.add_argument("--query-only", action="store_true",
                   help="Print the JSON queries that would be issued and exit.")
    args = p.parse_args()

    try:
        epitopes = collect_epitopes(args)
    except SystemExit:
        raise
    except Exception as exc:
        sys.stderr.write(f"[error] {exc}\n")
        return 2

    sys.stderr.write(f"[info] {len(epitopes)} epitope(s) to process: "
                     f"{', '.join(s for s, _ in epitopes)}\n")

    if args.query_only:
        for seq, _ in epitopes:
            print(f"# Epitope: {seq}  (seqmotif + MHC UniProt filter)")
            print(json.dumps(build_query(seq), indent=2))
        return 0

    os.makedirs(args.output_dir, exist_ok=True)
    records_by_epitope: dict[str, list[dict]] = {}
    dropped_all: list[dict] = []

    for seq, label in epitopes:
        hits = search_epitope(seq)
        if not hits:
            records_by_epitope[seq] = []
            continue

        metadata = fetch_metadata(hits)
        records: list[dict] = []
        for pdb_id in hits:
            entry = metadata.get(pdb_id)
            if not entry:
                sys.stderr.write(f"[warn] {pdb_id}: no metadata returned\n")
                continue
            rec = classify_entry(pdb_id, seq, entry)
            if rec is None:
                continue
            # match_type derived post-hoc: free if the matched peptide chain's canonical
            # sequence equals the epitope; fused if it strictly contains it.
            pseq = rec.get("peptide_seq") or ""
            if pseq == seq:
                rec["match_type"] = "free"
            elif seq in pseq:
                rec["match_type"] = "fused"
            else:
                rec["match_type"] = "?"
            rec["epitope"] = seq
            rec["epitope_label"] = label
            if rec.get("dropped"):
                dropped_all.append(rec)
            else:
                records.append(rec)

        assign_groups(records)
        records_by_epitope[seq] = records

        # Download + split (per epitope so paths are nested under that epitope's dir)
        epitope_dir = os.path.join(args.output_dir, seq)
        for r in records:
            entry_dir = os.path.join(epitope_dir, r["pdb_id"])
            os.makedirs(entry_dir, exist_ok=True)
            if not args.no_download:
                pdb_path = download_coords(r["pdb_id"], epitope_dir)
                r["pdb_path"] = pdb_path
            else:
                r["pdb_path"] = None

            if r["pdb_path"] and not args.no_split:
                split_result = split_entry(
                    pdb_path=r["pdb_path"],
                    output_dir=entry_dir,
                    pdb_id=r["pdb_id"],
                    peptide_chains=r["peptide_chains"],
                    mhc_chains=r["mhc_heavy_chains"] + r["mhc_ii_alpha_chains"],
                    light_chains=r["light_chains"],
                    tcr_alpha_chains=r["tcr_alpha_chains"],
                    tcr_beta_chains=r["tcr_beta_chains"],
                )
                r["pmhc_path"] = split_result["pmhc"]
                r["tcr_path"] = split_result["tcr"]
            else:
                r["pmhc_path"] = None
                r["tcr_path"] = None

    # Reports
    md_path = write_markdown_report(records_by_epitope, args.output_dir, epitopes)
    html_path = write_html_report(records_by_epitope, args.output_dir, epitopes)

    json_path = os.path.join(args.output_dir, "report_data.json")
    Path(json_path).write_text(
        json.dumps({
            "generated": datetime.datetime.now().isoformat(timespec="seconds"),
            "epitopes": [{"sequence": s, "label": l} for s, l in epitopes],
            "records_by_epitope": records_by_epitope,
            "dropped": dropped_all,
        }, indent=2, default=str),
        encoding="utf-8",
    )

    dropped_csv_path = os.path.join(args.output_dir, "report_dropped.csv")
    with open(dropped_csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["epitope", "label", "pdb_id", "reason", "descriptions"])
        for d in dropped_all:
            w.writerow([
                d.get("epitope", ""),
                d.get("epitope_label", ""),
                d.get("pdb_id", ""),
                d.get("reason", ""),
                "; ".join(d.get("descriptions", [])),
            ])

    sys.stderr.write(f"[info] Report (md):    {md_path}\n")
    sys.stderr.write(f"[info] Report (html):  {html_path}\n")
    sys.stderr.write(f"[info] Data (json):    {json_path}\n")
    sys.stderr.write(f"[info] Dropped (csv):  {dropped_csv_path} ({len(dropped_all)} entries)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
