#!/usr/bin/env python3
"""epitope_ss.py — per-epitope DSSP secondary structure from host-mapping results.

Consumes the output directory of `search-epitope-host` (epitopes.csv + structures/<ACC>.csv
+ parents/<ACC>.fasta) and, for every (epitope, parent) row that has a verified position on the
canonical UniProt sequence, computes the secondary structure at the residues the epitope maps to.

For each position it produces a probability distribution over the 8 DSSP states
(H,G,I,E,B,T,S,C) aggregated across the parent's structures, plus a simplified 3-state (H/E/C)
distribution. Only experimentally resolved structures are used (AlphaFold-only parents are skipped
and flagged unless --use-alphafold). By default, isolated peptide chains in pMHC / TCR / antibody
complexes are excluded (--min-align-length) so the result reflects the secondary structure the
epitope adopts in its PARENT fold, not its MHC-bound conformation.

8-state DSSP is computed locally with `mkdssp` (the RCSB Data API only serves coarse 3-state SS).
The RCSB Data API is used for the UniProt->structure residue mapping (SIFTS alignment).

Stdlib only, except the external `mkdssp` binary (install: conda install -c conda-forge dssp).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from html import escape
from pathlib import Path

# --------------------------------------------------------------------------- constants
RCSB_DATA = "https://data.rcsb.org/rest/v1/core"
RCSB_FILES = "https://files.rcsb.org/download"
UA = {"User-Agent": "epitope-secondary-structure/1.0 (comp-immunology)"}

# The 8 DSSP categories. "C" is the coil/loop/none bucket (DSSP blank / '-' / 'P' all fold here).
SS8 = ["H", "G", "I", "E", "B", "T", "S", "C"]
SS8_INDEX = {s: i for i, s in enumerate(SS8)}
# 8 -> 3 collapse (standard): {H,G,I}->H, {E,B}->E, rest->C.
SS8_TO_3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E",
            "T": "C", "S": "C", "C": "C"}
SS3 = ["H", "E", "C"]

# Colours (light theme) for the HTML report.
SS8_COLORS = {"H": "#d32f2f", "G": "#f06292", "I": "#ad1457",
              "E": "#1976d2", "B": "#64b5f6",
              "T": "#fbc02d", "S": "#aed581", "C": "#bdbdbd"}
SS3_COLORS = {"H": "#d32f2f", "E": "#1976d2", "C": "#757575"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}

_log_lock = threading.Lock()


def _log(msg: str) -> None:
    with _log_lock:
        sys.stderr.write(msg.rstrip() + "\n")
        sys.stderr.flush()


# --------------------------------------------------------------------------- HTTP / IO
def _get_json(url: str, timeout: float = 30.0, retries: int = 2):
    """GET a URL and parse JSON. 404 -> None. Retries transient failures with backoff."""
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            last = e
        if attempt < retries:
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed: {url}: {last}")


def _download(url: str, dest: Path, timeout: float = 90.0, retries: int = 2) -> None:
    """Download url to dest atomically (tmp + rename). Raises on persistent failure."""
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError,
                ConnectionError, OSError) as e:
            last = e
        if attempt < retries:
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"download failed: {url}: {last}")


# --------------------------------------------------------------------------- input parsing
def parse_occurrences(field: str):
    """'300-308; 426-434' -> [(300,308),(426,434)]; '' -> []."""
    out = []
    for chunk in (field or "").split(";"):
        chunk = chunk.strip()
        m = re.match(r"^(\d+)\s*-\s*(\d+)$", chunk)
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
    return out


def parse_chains_field(s: str):
    """'A/B/C=21-221; D=5-99' -> [(['A','B','C'],21,221),(['D'],5,99)].
    Ranges are UniProt-numbered (from the structures CSV, SIFTS-derived). Used for prescreen."""
    groups = []
    for grp in (s or "").split(";"):
        grp = grp.strip()
        if not grp:
            continue
        if "=" in grp:
            chains_part, rng = grp.split("=", 1)
            m = re.match(r"^(\d+)\s*-\s*(\d+)$", rng.strip())
            rng_t = (int(m.group(1)), int(m.group(2))) if m else None
        else:
            chains_part, rng_t = grp, None
        chains = [c for c in chains_part.split("/") if c]
        groups.append((chains, rng_t[0] if rng_t else None, rng_t[1] if rng_t else None))
    return groups


def load_epitopes(csv_path: Path):
    rows = []
    with csv_path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            r["_occ"] = parse_occurrences(r.get("occurrences", ""))
            rows.append(r)
    return rows


def load_structures(structures_dir: Path, acc: str):
    path = structures_dir / f"{acc}.csv"
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def load_parent_seq(parents_dir: Path, acc: str):
    path = parents_dir / f"{acc}.fasta"
    if not path.exists():
        return None
    seq = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            continue
        seq.append(line.strip())
    return "".join(seq) or None


def select_structures(struct_rows, use_af: bool):
    """experimental PDB present -> use all PDB rows. Otherwise the parent is *skipped*:
    AlphaFold-only -> 'no_experimental_structure' (flagged, not computed) unless --use-alphafold;
    no structure at all -> 'no_structure'."""
    pdb = [r for r in struct_rows if r.get("source") == "PDB"]
    af = [r for r in struct_rows if r.get("source") == "AlphaFold"]
    if pdb:
        return "experimental", pdb
    if af:
        return ("alphafold", af) if use_af else ("no_experimental_structure", [])
    return "no_structure", []


def window_overlaps_structure(struct_row, u_start, u_end) -> bool:
    """Cheap prescreen: does the epitope UniProt window overlap any of the structure's ranges?"""
    if struct_row.get("source") == "AlphaFold":
        rng = struct_row.get("range", "")
        m = re.match(r"^(\d+)\s*-\s*(\d+)$", rng.strip()) if rng else None
        if not m:
            return True  # unknown coverage -> don't skip
        s, e = int(m.group(1)), int(m.group(2))
        return not (u_end < s or u_start > e)
    groups = parse_chains_field(struct_row.get("chains", ""))
    if not groups or all(g[1] is None for g in groups):
        return True  # no range info -> don't skip
    for _chains, s, e in groups:
        if s is None:
            return True
        if not (u_end < s or u_start > e):
            return True
    return False


# --------------------------------------------------------------------------- mmCIF parse
def parse_mmcif(cif_path: Path):
    """Header-driven parse of the _atom_site loop (CIF column order varies: RCSB vs AlphaFold).

    Returns dict with, for the first CA atom of model 1 per residue (first altloc):
      by_label[(label_asym_id, label_seq_id:int)] = (auth_asym, auth_seq:int, icode, aa1, bfac)
      by_auth[(auth_asym_id, auth_seq_id:int, icode)] = aa1            # for AF identity mapping
    """
    cols = []
    in_header = False
    by_label = {}
    by_auth = {}
    seen = set()
    with cif_path.open() as fh:
        reading = False
        ncols = 0
        ix = {}
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".", 1)[1])
                in_header = True
                continue
            if in_header and not reading:
                # first non-header line after the column list -> data starts
                ix = {c: i for i, c in enumerate(cols)}
                ncols = len(cols)
                reading = True
            if reading:
                if line.startswith(("#", "_", "loop_")) or not line.strip():
                    break
                t = line.split()
                if len(t) < ncols:
                    continue
                if t[ix["group_PDB"]] != "ATOM":
                    continue
                if t[ix["label_atom_id"]] != "CA":
                    continue
                if "pdbx_PDB_model_num" in ix and t[ix["pdbx_PDB_model_num"]] != "1":
                    continue
                alt = t[ix["label_alt_id"]] if "label_alt_id" in ix else "."
                if alt not in (".", "A"):
                    continue
                try:
                    label_seq = int(t[ix["label_seq_id"]])
                    auth_seq = int(t[ix["auth_seq_id"]])
                except (ValueError, KeyError):
                    continue
                label_asym = t[ix["label_asym_id"]]
                auth_asym = t[ix["auth_asym_id"]]
                comp = t[ix["label_comp_id"]]
                icode_raw = t[ix["pdbx_PDB_ins_code"]] if "pdbx_PDB_ins_code" in ix else "?"
                icode = "" if icode_raw in (".", "?") else icode_raw
                try:
                    bfac = float(t[ix["B_iso_or_equiv"]])
                except (ValueError, KeyError):
                    bfac = None
                aa1 = THREE_TO_ONE.get(comp.upper(), "X")
                key = (label_asym, label_seq)
                if key in seen:
                    continue
                seen.add(key)
                by_label[key] = (auth_asym, auth_seq, icode, aa1, bfac)
                by_auth[(auth_asym, auth_seq, icode)] = aa1
    return {"by_label": by_label, "by_auth": by_auth}


# --------------------------------------------------------------------------- DSSP
def compute_dssp(cif_path: Path, mkdssp: str, dssp_out: Path):
    """Run mkdssp on a structure file and parse the classic residue table.

    Returns {(auth_chain, auth_seq:int, icode): ss8} where ss8 in SS8.
    Caches the raw .dssp file at dssp_out. Raises RuntimeError on failure.
    """
    if not dssp_out.exists():
        try:
            proc = subprocess.run(
                [mkdssp, "-i", str(cif_path), "-o", str(dssp_out)],
                capture_output=True, text=True, timeout=300,
            )
        except FileNotFoundError as e:
            raise RuntimeError(f"mkdssp not found ({mkdssp}): {e}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("mkdssp timed out")
        if proc.returncode != 0 or not dssp_out.exists():
            raise RuntimeError(f"mkdssp failed (rc={proc.returncode}): {proc.stderr[:300]}")
    return _parse_dssp(dssp_out)


def _parse_dssp(path: Path):
    out = {}
    started = False
    for line in path.read_text(errors="replace").splitlines():
        if not started:
            if line.startswith("  #  RESIDUE"):
                started = True
            continue
        if len(line) < 17:
            continue
        aa = line[13]
        if aa == "!":  # chain break marker
            continue
        resnum_field = line[5:10].strip()
        if not resnum_field:
            continue
        try:
            resnum = int(resnum_field)
        except ValueError:
            continue
        icode = line[10].strip()
        chain = line[11]
        ss_raw = line[16]
        ss = ss_raw if ss_raw in SS8_INDEX else "C"
        out[(chain, resnum, icode)] = ss
    return out


# --------------------------------------------------------------------------- RCSB align
def get_entity_align(entry: str, timeout: float):
    """Return list of {asym_ids:[label...], aligns:{acc: [(ref_beg, ent_beg, length)]}} per entity."""
    ent = _get_json(f"{RCSB_DATA}/entry/{entry}", timeout=timeout)
    if not ent:
        return []
    eids = (ent.get("rcsb_entry_container_identifiers", {}) or {}).get("polymer_entity_ids", []) or []
    entities = []
    for eid in eids:
        pe = _get_json(f"{RCSB_DATA}/polymer_entity/{entry}/{eid}", timeout=timeout)
        if not pe:
            continue
        ci = pe.get("rcsb_polymer_entity_container_identifiers", {}) or {}
        asym_ids = ci.get("asym_ids", []) or []
        aligns = {}
        for al in pe.get("rcsb_polymer_entity_align", []) or []:
            if al.get("reference_database_name") != "UniProt":
                continue
            acc = al.get("reference_database_accession")
            regs = [(r["ref_beg_seq_id"], r["entity_beg_seq_id"], r["length"])
                    for r in al.get("aligned_regions", []) or []]
            if acc and regs:
                aligns.setdefault(acc, []).extend(regs)
        entities.append({"asym_ids": asym_ids, "aligns": aligns})
    return entities


def map_uniprot_to_label(u: int, regions):
    """UniProt residue -> entity label_seq_id using SIFTS aligned_regions. None if not covered."""
    for ref_beg, ent_beg, length in regions:
        if ref_beg <= u < ref_beg + length:
            return u - ref_beg + ent_beg
    return None


# --------------------------------------------------------------------------- structure prep (cached)
class Cache:
    def __init__(self):
        self.lock = threading.Lock()
        self.entries = {}          # entry_id -> prep dict or {"error": ...}
        self.entry_locks = {}

    def get_lock(self, entry_id):
        with self.lock:
            if entry_id not in self.entry_locks:
                self.entry_locks[entry_id] = threading.Lock()
            return self.entry_locks[entry_id]


def prepare_structure(entry_id, source_type, struct_row, cache: Cache,
                      cache_dir: Path, mkdssp: str, offline: bool, timeout: float):
    """Download structure, parse mmCIF, run DSSP, fetch SIFTS align (PDB). Cached by entry_id."""
    lock = cache.get_lock(entry_id)
    with lock:
        if entry_id in cache.entries:
            return cache.entries[entry_id]
        prep = {"entry_id": entry_id, "source_type": source_type, "struct_row": struct_row}
        try:
            cif_dir = cache_dir / "cif"
            cif_dir.mkdir(parents=True, exist_ok=True)
            cif_path = cif_dir / f"{entry_id}.cif"
            if source_type == "experimental":
                if not cif_path.exists():
                    if offline:
                        raise RuntimeError("offline: cif not cached")
                    _download(f"{RCSB_FILES}/{entry_id}.cif", cif_path, timeout=timeout)
            else:  # alphafold
                if not cif_path.exists():
                    if offline:
                        raise RuntimeError("offline: cif not cached")
                    url = (struct_row.get("url") or "")
                    cif_url = re.sub(r"\.pdb$", ".cif", url)
                    if not cif_url.endswith(".cif"):
                        raise RuntimeError(f"cannot derive cif url from {url!r}")
                    _download(cif_url, cif_path, timeout=timeout)
            prep["mmcif"] = parse_mmcif(cif_path)
            dssp_dir = cache_dir / "dssp"
            dssp_dir.mkdir(parents=True, exist_ok=True)
            prep["dssp"] = compute_dssp(cif_path, mkdssp, dssp_dir / f"{entry_id}.dssp")
            if source_type == "experimental":
                if offline:
                    align_path = cache_dir / "align" / f"{entry_id}.json"
                    prep["entities"] = json.loads(align_path.read_text()) if align_path.exists() else []
                else:
                    ents = get_entity_align(entry_id, timeout=timeout)
                    prep["entities"] = ents
                    apath = cache_dir / "align"
                    apath.mkdir(parents=True, exist_ok=True)
                    (apath / f"{entry_id}.json").write_text(json.dumps(ents))
            prep["error"] = None
        except Exception as e:  # noqa: BLE001 - record and continue
            prep["error"] = str(e)
        cache.entries[entry_id] = prep
        return prep


# --------------------------------------------------------------------------- per-structure observation
def observe_window(prep, source_type, parent_acc, u_start, u_end, epitope, min_align_length=0):
    """For one structure, return per-position lists of (ss8, mismatch) over modeled chains,
    plus pLDDT (AF) per position. positions index 0..L-1 over the epitope window.

    If min_align_length > 0, an experimental entity whose SIFTS alignment to the parent spans
    fewer than that many residues is treated as an isolated peptide (e.g. the epitope presented
    in a pMHC complex) and skipped with note 'peptide_only'.

    Returns: per_pos (list of {"states":[ss8...], "mismatch":bool, "plddt":[...]}) or None on failure.
    """
    L = u_end - u_start + 1
    if prep.get("error"):
        return None, prep["error"]
    by_label = prep["mmcif"]["by_label"]
    by_auth = prep["mmcif"]["by_auth"]
    dssp = prep["dssp"]
    per_pos = [{"states": [], "mismatch": False, "plddt": []} for _ in range(L)]
    covered = False

    if source_type == "experimental":
        # find entity (and asym_ids + regions) matching the parent accession
        ent = None
        for e in prep.get("entities", []):
            if parent_acc in e.get("aligns", {}):
                ent = e
                break
        if ent is None:
            return per_pos, "not_aligned"
        regions = ent["aligns"][parent_acc]
        if min_align_length and sum(length for (_r, _e, length) in regions) < min_align_length:
            return per_pos, "peptide_only"
        asym_ids = ent["asym_ids"]
        for i in range(L):
            u = u_start + i
            label_seq = map_uniprot_to_label(u, regions)
            if label_seq is None:
                continue
            for asym in asym_ids:
                rec = by_label.get((asym, label_seq))
                if rec is None:
                    continue  # residue unmodeled in this chain
                auth_asym, auth_seq, icode, aa1, _bfac = rec
                ss = dssp.get((auth_asym, auth_seq, icode))
                if ss is None:
                    continue
                covered = True
                per_pos[i]["states"].append(ss)
                if aa1 not in ("X",) and aa1 != epitope[i]:
                    per_pos[i]["mismatch"] = True
    else:  # alphafold: auth_seq_id == UniProt position, single chain
        for i in range(L):
            u = u_start + i
            # AF chain id is whatever auth_asym the model uses (usually 'A'); scan keys for this u
            hit = None
            for (c, rn, ic), aa1 in by_auth.items():
                if rn == u and ic == "":
                    hit = (c, rn, ic, aa1)
                    break
            if hit is None:
                continue
            c, rn, ic, aa1 = hit
            ss = dssp.get((c, rn, ic))
            if ss is None:
                continue
            covered = True
            per_pos[i]["states"].append(ss)
            if aa1 not in ("X",) and aa1 != epitope[i]:
                per_pos[i]["mismatch"] = True
            rec = by_label.get(None)  # plddt via by_label scan
            # pull pLDDT from by_label for this auth residue
            for (la, ls), (aasym, aseq, aic, baa, bfac) in by_label.items():
                if aasym == c and aseq == rn and aic == ic:
                    if bfac is not None:
                        per_pos[i]["plddt"].append(bfac)
                    break
    return per_pos, (None if covered else "uncovered")


# --------------------------------------------------------------------------- aggregation
def reduce_structure_vector(states, mode):
    """states: list of ss8 over a structure's chains at one position -> dict state->prob."""
    if not states:
        return None
    v = {s: 0.0 for s in SS8}
    for s in states:
        v[s] += 1.0
    n = len(states)
    return {s: v[s] / n for s in SS8}


def aggregate(struct_results, L, mode):
    """struct_results: list of per_pos (one per contributing structure).
    Returns list of position dicts with prob8, prob3, n_structures, n_chains, plddt, argmax."""
    positions = []
    for i in range(L):
        if mode == "per-chain":
            pooled = []
            n_struct = 0
            n_chains = 0
            plddt = []
            for per_pos in struct_results:
                states = per_pos[i]["states"]
                if states:
                    n_struct += 1
                    n_chains += len(states)
                    pooled.extend(states)
                    plddt.extend(per_pos[i]["plddt"])
            if pooled:
                v = {s: 0.0 for s in SS8}
                for s in pooled:
                    v[s] += 1.0
                prob8 = {s: v[s] / len(pooled) for s in SS8}
            else:
                prob8 = None
        else:  # per-structure (default)
            vecs = []
            n_chains = 0
            plddt = []
            for per_pos in struct_results:
                states = per_pos[i]["states"]
                vs = reduce_structure_vector(states, mode)
                if vs is not None:
                    vecs.append(vs)
                    n_chains += len(states)
                    plddt.extend(per_pos[i]["plddt"])
            n_struct = len(vecs)
            if vecs:
                prob8 = {s: sum(v[s] for v in vecs) / len(vecs) for s in SS8}
            else:
                prob8 = None

        if prob8 is None:
            positions.append({"prob8": None, "prob3": None, "n_structures": 0,
                              "n_chains": 0, "plddt": None, "argmax8": None, "argmax3": None})
            continue
        prob3 = {c: 0.0 for c in SS3}
        for s in SS8:
            prob3[SS8_TO_3[s]] += prob8[s]
        argmax8 = max(SS8, key=lambda s: (prob8[s], -SS8_INDEX[s]))
        argmax3 = max(SS3, key=lambda c: (prob3[c], -SS3.index(c)))
        positions.append({
            "prob8": prob8, "prob3": prob3, "n_structures": n_struct,
            "n_chains": n_chains,
            "plddt": (sum(plddt) / len(plddt)) if plddt else None,
            "argmax8": argmax8, "argmax3": argmax3,
        })
    return positions


def consensus_strings(positions):
    ss8 = "".join(p["argmax8"] if p["argmax8"] else "." for p in positions)
    ss3 = "".join(p["argmax3"] if p["argmax3"] else "-" for p in positions)
    return ss8, ss3


# --------------------------------------------------------------------------- core per-row
def process_row(row, structures_dir, parents_dir, cache, cache_dir, mkdssp,
                use_af, max_per_parent, offline, timeout, min_align_length=0):
    epitope = row["epitope"]
    acc = row["uniprot_acc"]
    L = len(epitope)
    result = {
        "epitope": epitope, "length": L, "uniprot_acc": acc,
        "protein_name": row.get("protein_name", ""), "organism": row.get("organism", ""),
        "occurrence": "", "status": "", "source_type": "", "n_structures_used": 0,
        "positions_covered": 0, "mean_plddt": None, "ss3_consensus": "", "ss8_consensus": "",
        "positions": [], "structures": [], "flags": [],
    }
    occ = row["_occ"]
    if not occ:
        result["status"] = "no_position"
        return result
    u_start, u_end = occ[0]  # primary occurrence
    result["occurrence"] = f"{u_start}-{u_end}"
    if u_end - u_start + 1 != L:
        # occurrence length disagreement (shouldn't happen) -> trust occurrence span
        L = u_end - u_start + 1

    struct_rows = load_structures(structures_dir, acc)
    source_type, chosen = select_structures(struct_rows, use_af)
    result["source_type"] = source_type
    if source_type in ("no_structure", "no_experimental_structure"):
        result["status"] = source_type
        return result

    # prescreen by UniProt-range overlap, then (experimental) sort/cap by method+resolution
    chosen = [r for r in chosen if window_overlaps_structure(r, u_start, u_end)]
    if source_type == "experimental":
        def sortkey(r):
            res = r.get("resolution", "") or ""
            m = re.match(r"([\d.]+)", res)
            return (0 if r.get("method") == "X-ray" else 1,
                    float(m.group(1)) if m else float("inf"))
        chosen.sort(key=sortkey)
        if max_per_parent > 0:
            chosen = chosen[:max_per_parent]

    if not chosen:
        result["status"] = "not_in_any_structure"
        return result

    struct_results = []
    for r in chosen:
        entry_id = r["entry_id"]
        prep = prepare_structure(entry_id, source_type, r, cache, cache_dir,
                                 mkdssp, offline, timeout)
        per_pos, note = observe_window(prep, source_type, acc, u_start, u_end, epitope,
                                       min_align_length)
        srec = {
            "entry_id": entry_id, "method": r.get("method", ""),
            "resolution": r.get("resolution", ""), "version": r.get("version", ""),
            "url": r.get("url", ""), "dssp_status": "ok", "n_modeled": 0,
            "seq_mismatch": False,
        }
        if prep.get("error"):
            srec["dssp_status"] = "prep_failed"
            result["structures"].append(srec)
            result["flags"].append(f"{entry_id}:prep_failed")
            continue
        if per_pos is None or note == "not_aligned":
            srec["dssp_status"] = "not_aligned"
            result["structures"].append(srec)
            continue
        if note == "peptide_only":
            srec["dssp_status"] = "peptide_only"
            result["structures"].append(srec)
            continue
        n_modeled = sum(1 for p in per_pos if p["states"])
        srec["n_modeled"] = n_modeled
        srec["seq_mismatch"] = any(p["mismatch"] for p in per_pos)
        if note == "uncovered" or n_modeled == 0:
            srec["dssp_status"] = "unmodeled"
            result["structures"].append(srec)
            continue
        if srec["seq_mismatch"]:
            result["flags"].append(f"{entry_id}:seq_mismatch")
        result["structures"].append(srec)
        struct_results.append(per_pos)

    if not struct_results:
        # structures existed but none yielded a modeled, aligned, non-peptide window
        ds = [s["dssp_status"] for s in result["structures"]]
        if "peptide_only" in ds and all(d in ("peptide_only", "not_aligned", "prep_failed")
                                        for d in ds):
            result["status"] = "only_pmhc"
        elif "unmodeled" in ds:
            result["status"] = "unmodeled"
        else:
            result["status"] = "not_in_any_structure"
        return result

    positions = aggregate(struct_results, L, MODE)
    ss8, ss3 = consensus_strings(positions)
    result["ss8_consensus"] = ss8
    result["ss3_consensus"] = ss3
    result["positions"] = positions
    result["status"] = "ok"
    result["n_structures_used"] = max((p["n_structures"] for p in positions), default=0)
    result["positions_covered"] = sum(1 for p in positions if p["n_structures"] > 0)
    plddts = [p["plddt"] for p in positions if p["plddt"] is not None]
    if source_type == "alphafold" and plddts:
        result["mean_plddt"] = sum(plddts) / len(plddts)
    result["flags"] = sorted(set(result["flags"]))
    return result


# --------------------------------------------------------------------------- outputs
def _frac(positions, c):
    vals = [p["prob3"][c] for p in positions if p["prob3"] is not None]
    return sum(vals) / len(vals) if vals else None


def write_epitope_ss_csv(path, results):
    cols = ["epitope", "length", "uniprot_acc", "protein_name", "organism", "occurrence",
            "status", "source_type", "n_structures_used", "positions_covered", "mean_plddt",
            "ss3_consensus", "ss8_consensus", "frac_H", "frac_E", "frac_C",
            "n_positions", "n_positions_unmodeled", "n_mismatch_structures", "flags"]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in results:
            pos = r["positions"]
            n_unmod = sum(1 for p in pos if p["n_structures"] == 0) if pos else ""
            n_mm = sum(1 for s in r["structures"] if s.get("seq_mismatch"))
            w.writerow([
                r["epitope"], r["length"], r["uniprot_acc"], r["protein_name"], r["organism"],
                r["occurrence"], r["status"], r["source_type"], r["n_structures_used"],
                r["positions_covered"],
                f"{r['mean_plddt']:.1f}" if r["mean_plddt"] is not None else "",
                r["ss3_consensus"], r["ss8_consensus"],
                _fmt(_frac(pos, "H")), _fmt(_frac(pos, "E")), _fmt(_frac(pos, "C")),
                r["length"] if pos else "", n_unmod, n_mm, ";".join(r["flags"]),
            ])


def _fmt(x):
    return f"{x:.3f}" if isinstance(x, float) else ""


def write_positions_csv(path, results):
    cols = (["epitope", "uniprot_acc", "uni_pos", "epitope_index", "aa", "n_structures",
             "n_chains"] + [f"p_{s}" for s in SS8] + [f"p3_{c}" for c in SS3]
            + ["argmax8", "argmax3", "mean_plddt"])
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in results:
            if not r["positions"]:
                continue
            occ = r["occurrence"].split("-")
            u0 = int(occ[0]) if occ and occ[0] else None
            for i, p in enumerate(r["positions"]):
                aa = r["epitope"][i] if i < len(r["epitope"]) else ""
                row = [r["epitope"], r["uniprot_acc"],
                       (u0 + i) if u0 is not None else "", i + 1, aa,
                       p["n_structures"], p["n_chains"]]
                if p["prob8"]:
                    row += [f"{p['prob8'][s]:.3f}" for s in SS8]
                    row += [f"{p['prob3'][c]:.3f}" for c in SS3]
                    row += [p["argmax8"], p["argmax3"]]
                else:
                    row += [""] * (len(SS8) + len(SS3) + 2)
                row += [f"{p['plddt']:.1f}" if p["plddt"] is not None else ""]
                w.writerow(row)


def write_structures_csv(path, results):
    cols = ["epitope", "uniprot_acc", "source_type", "entry_id", "method", "resolution",
            "version", "dssp_status", "n_modeled", "seq_mismatch"]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in results:
            for s in r["structures"]:
                w.writerow([r["epitope"], r["uniprot_acc"], r["source_type"], s["entry_id"],
                            s["method"], s["resolution"], s["version"], s["dssp_status"],
                            s["n_modeled"], s["seq_mismatch"]])


def write_json(path, results, meta):
    path.write_text(json.dumps({"metadata": meta, "epitopes": results},
                               indent=2, default=str))


# --------------------------------------------------------------------------- HTML
def _legend_html():
    items8 = " ".join(
        f'<span class="chip" style="background:{SS8_COLORS[s]}">{s}</span>'
        f'<span class="muted">{lbl}</span>'
        for s, lbl in [("H", "α-helix"), ("G", "3₁₀"), ("I", "π"), ("E", "strand"),
                       ("B", "bridge"), ("T", "turn"), ("S", "bend"), ("C", "coil")])
    items3 = " ".join(
        f'<span class="chip" style="background:{SS3_COLORS[c]}">{c}</span>'
        for c in SS3)
    return (f'<div class="legend"><b>8-state:</b> {items8} &nbsp;|&nbsp; '
            f'<b>3-state:</b> {items3}</div>')


def _bar_html(p):
    """Stacked vertical 8-state probability bar for one position."""
    if not p["prob8"]:
        return '<div class="bar empty" title="no structure coverage">·</div>'
    segs = []
    for s in SS8:
        v = p["prob8"][s]
        if v <= 0:
            continue
        segs.append(f'<div style="height:{v*100:.1f}%;background:{SS8_COLORS[s]}" '
                    f'title="{s} {v:.2f}"></div>')
    return f'<div class="bar">{"".join(segs)}</div>'


def render_html(results, meta):
    css = """
  body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
       max-width:1100px;margin:2em auto;padding:0 1em;color:#222;line-height:1.5}
  h1{border-bottom:2px solid #444;padding-bottom:.3em}
  h2{margin-top:2em;border-bottom:1px solid #ccc;padding-bottom:.2em;font-size:1.15em}
  h3{margin-top:1.2em;color:#555;font-size:1em}
  table{border-collapse:collapse;width:100%;margin:.6em 0;font-size:13px}
  th,td{border:1px solid #ccc;padding:.35em .55em;text-align:left;vertical-align:top}
  th{background:#f4f4f4}
  code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
  .muted{color:#888;font-size:12px}
  .pill{display:inline-block;padding:.05em .5em;border-radius:10px;background:#eef;color:#225;
        font-size:12px}
  .pill.af{background:#e7f0ff;color:#1a4} .pill.none{background:#f3f3f3;color:#888}
  .bad{color:#b00020;font-weight:bold}
  .chip{display:inline-block;width:1.2em;text-align:center;color:#fff;border-radius:3px;
        font-size:11px;margin:0 .1em;font-weight:bold}
  .legend{margin:.8em 0;font-size:12px}
  .ssrow{display:flex;gap:3px;align-items:flex-end;margin:.4em 0}
  .ssrow .col{display:flex;flex-direction:column;align-items:center;width:34px}
  .bar{height:70px;width:26px;border:1px solid #ddd;display:flex;flex-direction:column-reverse;
       background:#fafafa}
  .bar.empty{align-items:center;justify-content:center;color:#bbb}
  .aa{font-family:ui-monospace,monospace;font-weight:bold;font-size:14px}
  .ss3{font-family:ui-monospace,monospace;font-weight:bold;color:#fff;border-radius:3px;
       width:100%;text-align:center}
  .ncount{font-size:10px;color:#999}
  .conskey{font-family:ui-monospace,monospace;font-weight:bold;letter-spacing:1px}
"""
    h = ["<!doctype html><html><head><meta charset='utf-8'>",
         "<title>Epitope secondary structure</title><style>", css, "</style></head><body>"]
    h.append("<h1>Epitope secondary structure (DSSP)</h1>")
    m = meta
    h.append(f"<p class='muted'>Input: {escape(m['input_dir'])} &middot; aggregation: "
             f"<b>{escape(m['aggregation'])}</b> &middot; DSSP: {escape(m['mkdssp'])} "
             f"&middot; {escape(m['generated'])}</p>")

    # overall summary
    by_status = {}
    by_src = {}
    for r in results:
        by_status[r["status"]] = by_status.get(r["status"], 0) + 1
        by_src[r["source_type"]] = by_src.get(r["source_type"], 0) + 1
    n_ok = by_status.get("ok", 0)
    n_mm = sum(1 for r in results if any(s.get("seq_mismatch") for s in r["structures"]))
    h.append("<ul>")
    h.append(f"<li><b>{len(results)}</b> epitope×parent rows; <b>{n_ok}</b> with secondary "
             f"structure computed.</li>")
    h.append("<li>Source: " + ", ".join(f"{escape(k)}={v}" for k, v in sorted(by_src.items()))
             + "</li>")
    h.append("<li>Status: " + ", ".join(f"{escape(k)}={v}" for k, v in sorted(by_status.items()))
             + "</li>")
    if n_mm:
        h.append(f"<li class='bad'>{n_mm} row(s) had a sequence mismatch between a structure "
                 f"and the epitope — see flags.</li>")
    h.append("</ul>")
    h.append(_legend_html())

    # summary table
    h.append("<h2>Summary</h2><table><tr><th>Epitope</th><th>UniProt</th><th>Protein</th>"
             "<th>Position</th><th>Source</th><th>3-state</th><th>Status</th><th>Flags</th></tr>")
    for r in results:
        cons = _color_ss3(r["ss3_consensus"]) if r["ss3_consensus"] else \
            f"<span class='muted'>—</span>"
        src = _src_pill(r)
        flags = ("<span class='bad'>" + escape("; ".join(r["flags"])) + "</span>") if r["flags"] else ""
        h.append(f"<tr><td class='mono'>{escape(r['epitope'])}</td>"
                 f"<td><a href='https://www.uniprot.org/uniprotkb/{escape(r['uniprot_acc'])}'>"
                 f"{escape(r['uniprot_acc'])}</a></td>"
                 f"<td>{escape(r['protein_name'])}</td>"
                 f"<td class='mono'>{escape(r['occurrence']) or '—'}</td>"
                 f"<td>{src}</td><td>{cons}</td>"
                 f"<td>{escape(r['status'])}</td><td>{flags}</td></tr>")
    h.append("</table>")

    # per-epitope detail
    h.append("<h2>Per-epitope detail</h2>")
    notes = {
        "no_experimental_structure": "no experimentally resolved structure "
                                     "(AlphaFold model available — skipped for now)",
        "no_structure": "no structure available (no PDB, no AlphaFold)",
        "no_position": "epitope not located on the parent's canonical UniProt sequence "
                       "(likely an EBV strain variant)",
        "not_in_any_structure": "epitope position not covered by any deposited structure",
        "unmodeled": "epitope residues are unresolved (unmodeled) in all covering structures",
        "only_pmhc": "only structures are pMHC complexes (MHC-bound peptide conformation, "
                     "not the parent fold) — excluded by --min-align-length",
    }
    for r in results:
        if r["status"] != "ok":
            note = notes.get(r["status"], r["status"])
            h.append(f"<h3 class='mono'>{escape(r['epitope'])} &middot; {escape(r['uniprot_acc'])}"
                     f"</h3><p class='muted'>{escape(note)}</p>")
            continue
        h.append(f"<h3 class='mono'>{escape(r['epitope'])} &middot; "
                 f"{escape(r['uniprot_acc'])} {escape(r['protein_name'])} "
                 f"&middot; {escape(r['occurrence'])} &middot; {_src_pill(r)}</h3>")
        # per-position bars
        h.append("<div class='ssrow'>")
        for i, p in enumerate(r["positions"]):
            aa = r["epitope"][i] if i < len(r["epitope"]) else "?"
            c3 = p["argmax3"] or "-"
            col = SS3_COLORS.get(c3, "#999")
            ncell = (str(p["n_structures"]) if p["n_structures"] else "0")
            h.append("<div class='col'>"
                     f"<div class='aa'>{escape(aa)}</div>"
                     f"{_bar_html(p)}"
                     f"<div class='ss3' style='background:{col}'>{escape(c3)}</div>"
                     f"<div class='ncount'>{ncell}</div></div>")
        h.append("</div>")
        h.append(f"<p class='muted'>8-state consensus: <span class='conskey'>"
                 f"{escape(r['ss8_consensus'])}</span> &middot; n = #structures contributing "
                 f"per position</p>")
        # contributing structures
        h.append("<table><tr><th>Structure</th><th>Method</th><th>Resolution</th>"
                 "<th>Modeled</th><th>DSSP</th></tr>")
        for s in r["structures"]:
            if r["source_type"] == "alphafold":
                ident = (f"<a href='https://alphafold.ebi.ac.uk'>{escape(s['entry_id'])}</a> "
                         f"{escape(s['version'])}")
                method = "AlphaFold" + (f" (pLDDT {r['mean_plddt']:.0f})"
                                        if r["mean_plddt"] is not None else "")
            else:
                ident = (f"<a href='https://www.rcsb.org/structure/{escape(s['entry_id'])}'>"
                         f"{escape(s['entry_id'])}</a>")
                method = escape(s["method"])
            mm = " <span class='bad'>mismatch</span>" if s.get("seq_mismatch") else ""
            h.append(f"<tr><td class='mono'>{ident}</td><td>{method}</td>"
                     f"<td>{escape(s['resolution'])}</td><td>{s['n_modeled']}</td>"
                     f"<td>{escape(s['dssp_status'])}{mm}</td></tr>")
        h.append("</table>")
    h.append(f"<p class='muted'>Generated by epitope-secondary-structure &middot; "
             f"DSSP 8-state computed with mkdssp; mapping via RCSB SIFTS &middot; "
             f"{escape(m['generated'])}</p>")
    h.append("</body></html>")
    return "".join(h)


def _color_ss3(s):
    return "".join(f"<span class='ss3' style='display:inline-block;width:1.1em;"
                   f"background:{SS3_COLORS.get(c, '#999')}'>{escape(c)}</span>" for c in s)


def _src_pill(r):
    st = r["source_type"]
    if st == "experimental":
        return (f"<span class='pill'>{r['n_structures_used']} PDB · "
                f"{r['positions_covered']}/{r['length']} pos</span>")
    if st == "alphafold":
        pl = f" pLDDT {r['mean_plddt']:.0f}" if r["mean_plddt"] is not None else ""
        return f"<span class='pill af'>AlphaFold{escape(pl)}</span>"
    if st == "no_experimental_structure":
        return "<span class='pill none'>AlphaFold only — skipped</span>"
    if st == "no_structure":
        return "<span class='pill none'>no structure</span>"
    return f"<span class='pill none'>{escape(st)}</span>"


# --------------------------------------------------------------------------- main
MODE = "per-structure"


def main(argv=None):
    global MODE
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir", help="search-epitope-host output dir (host_mapping/)")
    ap.add_argument("--epitopes-csv", default=None)
    ap.add_argument("--structures-dir", default=None)
    ap.add_argument("--parents-dir", default=None)
    ap.add_argument("-o", "--out-dir", default=None)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("-j", "--workers", type=int, default=8)
    ap.add_argument("--mkdssp-path", default="mkdssp")
    ap.add_argument("--aggregation", choices=["per-structure", "per-chain"],
                    default="per-structure")
    ap.add_argument("--use-alphafold", action="store_true",
                    help="Also annotate parents that have only AlphaFold models (default: skip "
                         "and flag them). NOTE: requires DSSP 4 — DSSP 3.1.2 cannot parse "
                         "AlphaFold mmCIF.")
    ap.add_argument("--max-structures-per-parent", type=int, default=0,
                    help="0 = use all (default). >0 caps to top-N by method+resolution.")
    ap.add_argument("--min-align-length", type=int, default=25,
                    help="Default 25: skip experimental entities whose SIFTS alignment to the "
                         "parent spans fewer residues — i.e. isolated peptides in pMHC / TCR / "
                         "antibody complexes — so only genuine parent-fold structures contribute "
                         "(epitopes left with only such structures become 'only_pmhc'). "
                         "Set 0 to include them (the MHC-bound peptide conformation).")
    ap.add_argument("--offline", action="store_true", help="use only cached cif/dssp/align")
    ap.add_argument("--timeout", type=float, default=60.0)
    args = ap.parse_args(argv)

    MODE = args.aggregation
    in_dir = Path(args.input_dir)
    epi_csv = Path(args.epitopes_csv) if args.epitopes_csv else in_dir / "epitopes.csv"
    structures_dir = Path(args.structures_dir) if args.structures_dir else in_dir / "structures"
    parents_dir = Path(args.parents_dir) if args.parents_dir else in_dir / "parents"
    out_dir = Path(args.out_dir) if args.out_dir else in_dir / "secondary_structure"
    cache_dir = Path(args.cache_dir) if args.cache_dir else out_dir / "cache"
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    if not epi_csv.exists():
        sys.exit(f"epitopes.csv not found: {epi_csv}")

    # verify mkdssp unless fully offline (still want it; offline only skips network)
    try:
        subprocess.run([args.mkdssp_path, "--version"], capture_output=True, timeout=30)
    except FileNotFoundError:
        sys.exit(f"mkdssp not found ({args.mkdssp_path}).\n"
                 f"Install: conda install -c conda-forge dssp\n"
                 f"then pass --mkdssp-path \"$(conda run -n <env> which mkdssp)\".")

    rows = load_epitopes(epi_csv)
    _log(f"loaded {len(rows)} epitope×parent rows")
    cache = Cache()

    def work(row):
        try:
            return process_row(row, structures_dir, parents_dir, cache, cache_dir,
                               args.mkdssp_path, args.use_alphafold,
                               args.max_structures_per_parent, args.offline, args.timeout,
                               args.min_align_length)
        except Exception as e:  # noqa: BLE001
            _log(f"ERROR {row.get('epitope')} {row.get('uniprot_acc')}: {e}")
            return {"epitope": row.get("epitope", ""), "length": len(row.get("epitope", "")),
                    "uniprot_acc": row.get("uniprot_acc", ""),
                    "protein_name": row.get("protein_name", ""),
                    "organism": row.get("organism", ""), "occurrence": "",
                    "status": f"error: {e}", "source_type": "", "n_structures_used": 0,
                    "positions_covered": 0, "mean_plddt": None, "ss3_consensus": "",
                    "ss8_consensus": "", "positions": [], "structures": [], "flags": []}

    if args.workers <= 1:
        results = [work(r) for r in rows]
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(work, rows))

    meta = {
        "input_dir": str(in_dir), "aggregation": MODE, "mkdssp": args.mkdssp_path,
        "use_alphafold": args.use_alphafold, "min_align_length": args.min_align_length,
        "n_rows": len(results),
        "generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
    }
    write_epitope_ss_csv(out_dir / "epitope_ss.csv", results)
    write_positions_csv(out_dir / "epitope_ss_positions.csv", results)
    write_structures_csv(out_dir / "structures_used.csv", results)
    write_json(out_dir / "epitope_ss.json", results, meta)
    (out_dir / "report.html").write_text(render_html(results, meta))

    n_ok = sum(1 for r in results if r["status"] == "ok")
    _log(f"done: {n_ok}/{len(results)} rows with SS; outputs in {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
