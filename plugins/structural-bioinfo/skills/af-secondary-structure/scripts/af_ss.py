#!/usr/bin/env python3
"""af_ss.py — secondary structure of an AlphaFold (or any) protein structure, with
pLDDT-based confidence filtering.

For each sequence position it assigns the 8-state DSSP code (H,G,I,E,B,T,S,C) and the
reduced 3-state code (H/E/C), then layers on AlphaFold quality control: the per-residue
pLDDT (read from the B-factor column) is classified into the four standard confidence
bands, and a *masked* secondary-structure string is produced in which low-confidence
positions (pLDDT below a cutoff, default 70) are replaced with `X`. This is non-destructive:
the raw assignment, the masked assignment, and the per-residue pLDDT/band are all reported,
so unreliable disordered regions can be flagged or removed downstream without losing the
underlying call.

DSSP is computed with mdtraj (a self-contained reimplementation of the Kabsch–Sander
algorithm) — no external `mkdssp` binary, and it reads AlphaFold mmCIF directly (which
mkdssp 3.x cannot). Works on PDB and mmCIF, single files or batches.

Engine lives in a conda env with mdtraj (typically `mdanalysis`); the script auto-discovers
it and re-execs if the active interpreter lacks mdtraj.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from html import escape
from pathlib import Path


# --------------------------------------------------------------------------- mdtraj bootstrap
def _bootstrap_mdtraj() -> None:
    """Ensure mdtraj is importable; else find a conda env that has it and re-exec there."""
    try:
        import mdtraj  # noqa: F401
        return
    except Exception:
        pass
    if os.environ.get("AF_SS_REEXECED") == "1":
        sys.stderr.write(
            "ERROR: mdtraj is required but could not be imported.\n"
            "Install it in a conda env, e.g.:\n"
            "  conda create -n mdanalysis -c conda-forge mdtraj\n"
            "  # or: conda install -n <env> -c conda-forge mdtraj\n")
        sys.exit(2)
    homes = [Path.home() / d for d in ("miniconda3", "anaconda3", "miniforge3", "mambaforge")]
    homes.append(Path("/opt/conda"))
    prefer = ["mdanalysis", "mdtraj"]
    envs = []
    for h in homes:
        envdir = h / "envs"
        if envdir.is_dir():
            for e in sorted(envdir.iterdir()):
                py = e / "bin" / "python"
                if py.exists():
                    envs.append((e.name, py))
    envs.sort(key=lambda t: (prefer.index(t[0]) if t[0] in prefer else len(prefer), t[0]))
    for _name, py in envs:
        try:
            r = subprocess.run([str(py), "-c", "import mdtraj"],
                               capture_output=True, timeout=120)
        except Exception:
            continue
        if r.returncode == 0:
            os.environ["AF_SS_REEXECED"] = "1"
            os.execv(str(py), [str(py), os.path.abspath(__file__)] + sys.argv[1:])
    sys.stderr.write(
        "ERROR: no conda env with mdtraj found (searched miniconda3/anaconda3/miniforge3/"
        "mambaforge/opt-conda envs).\nInstall it, e.g.:\n"
        "  conda create -n mdanalysis -c conda-forge mdtraj\n")
    sys.exit(2)


_bootstrap_mdtraj()
import mdtraj as md  # noqa: E402


# --------------------------------------------------------------------------- constants
# 8 DSSP categories; "C" is the coil/loop bucket (mdtraj emits a blank for it).
SS8 = ["H", "G", "I", "E", "B", "T", "S", "C"]
SS8_SET = set(SS8)
# 8 -> 3 collapse (standard): {H,G,I}->H, {E,B}->E, rest->C.
SS8_TO_3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E",
            "T": "C", "S": "C", "C": "C"}
SS3 = ["H", "E", "C"]
MASK = "X"

# Light-theme colours for the HTML report (matched to the sibling epitope-secondary-structure
# skill for plugin consistency).
SS8_COLORS = {"H": "#d32f2f", "G": "#f06292", "I": "#ad1457",
              "E": "#1976d2", "B": "#64b5f6",
              "T": "#fbc02d", "S": "#aed581", "C": "#bdbdbd"}
SS3_COLORS = {"H": "#d32f2f", "E": "#1976d2", "C": "#757575"}
MASK_COLOR = "#9e9e9e"

# Official AlphaFold pLDDT palette + band thresholds.
BANDS = ["very_low", "low", "confident", "very_high"]
BAND_COLORS = {"very_low": "#FF7D45", "low": "#FFDB13",
               "confident": "#65CBF3", "very_high": "#0053D6"}
BAND_LABELS = {"very_low": "very low (<50)", "low": "low (50–70)",
               "confident": "confident (70–90)", "very_high": "very high (≥90)"}

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


def band_of(plddt):
    if plddt is None:
        return "unknown"
    if plddt < 50:
        return "very_low"
    if plddt < 70:
        return "low"
    if plddt < 90:
        return "confident"
    return "very_high"


# --------------------------------------------------------------------------- pLDDT (B-factor) readers
def read_ca_plddt(path: Path):
    """Return {(chain_id, res_seq): bfactor} over CA atoms. The B-factor column of an
    AlphaFold model holds the per-residue pLDDT. Handles PDB and mmCIF."""
    name = path.name.lower()
    if name.endswith((".cif", ".mmcif", ".cif.gz")):
        return _plddt_from_cif(path)
    return _plddt_from_pdb(path)


def _plddt_from_pdb(path: Path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):          # first altloc only
                continue
            chain = line[21].strip() or " "
            try:
                res_seq = int(line[22:26])
            except ValueError:
                continue
            try:
                b = float(line[60:66])
            except ValueError:
                b = None
            out.setdefault((chain, res_seq), b)
    return out


def _plddt_from_cif(path: Path):
    """Header-driven parse of the _atom_site loop (column order varies between sources)."""
    cols, out, seen = [], {}, set()
    reading, ix, ncols = False, {}, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".", 1)[1])
                continue
            if cols and not reading:
                ix = {c: i for i, c in enumerate(cols)}
                ncols = len(cols)
                reading = True
            if reading:
                if line.startswith(("#", "loop_", "_")) or not line.strip():
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
                    auth_seq = int(t[ix["auth_seq_id"]])
                except (ValueError, KeyError):
                    continue
                auth_asym = t[ix["auth_asym_id"]]
                try:
                    b = float(t[ix["B_iso_or_equiv"]])
                except (ValueError, KeyError):
                    b = None
                key = (auth_asym, auth_seq)
                if key in seen:
                    continue
                seen.add(key)
                out[key] = b
    return out


# --------------------------------------------------------------------------- core analysis
def _c(x) -> str:
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def _norm_ss8(raw: str) -> str:
    if raw in ("H", "G", "I", "E", "B", "T", "S"):
        return raw
    return "C"          # blank / '-' / 'P' / anything else -> coil


def analyze_structure(path: Path, cutoff: float, mask: bool):
    """Compute per-residue SS (8/3), pLDDT band, masking, disordered segments, summary."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        traj = md.load(str(path))
    ss8_arr = md.compute_dssp(traj, simplified=False)[0]
    ss3_arr = md.compute_dssp(traj, simplified=True)[0]
    plddt_map = read_ca_plddt(path)

    residues = []
    for res, r8, r3 in zip(traj.topology.residues, ss8_arr, ss3_arr):
        raw8 = _c(r8)
        if raw8 == "NA" or not res.is_protein:
            continue                                  # ligands, ions, water, NA
        ss8 = _norm_ss8(raw8)
        ss3 = _c(r3)
        if ss3 not in SS3:
            ss3 = SS8_TO_3.get(ss8, "C")
        chain = getattr(res.chain, "chain_id", None) or str(res.chain.index)
        resid = int(res.resSeq)
        plddt = plddt_map.get((chain, resid))
        b = band_of(plddt)
        is_masked = bool(mask and plddt is not None and plddt < cutoff)
        residues.append({
            "chain": chain, "resid": resid, "resname": res.name,
            "aa": THREE_TO_ONE.get(res.name.upper(), "X"),
            "ss8": ss8, "ss3": ss3,
            "plddt": round(plddt, 2) if plddt is not None else None,
            "band": b, "masked": is_masked,
            "ss8_masked": MASK if is_masked else ss8,
            "ss3_masked": MASK if is_masked else ss3,
        })

    n = len(residues)
    sequence = "".join(r["aa"] for r in residues)
    ss8 = "".join(r["ss8"] for r in residues)
    ss3 = "".join(r["ss3"] for r in residues)
    ss8_masked = "".join(r["ss8_masked"] for r in residues)
    ss3_masked = "".join(r["ss3_masked"] for r in residues)

    # disordered segments: contiguous runs (within a chain) of pLDDT < cutoff
    segments = []
    cur = None
    for r in residues:
        below = r["plddt"] is not None and r["plddt"] < cutoff
        if below and cur and cur["chain"] == r["chain"]:
            cur["resids"].append(r["resid"])
            cur["plddts"].append(r["plddt"])
        elif below:
            if cur:
                segments.append(cur)
            cur = {"chain": r["chain"], "resids": [r["resid"]], "plddts": [r["plddt"]]}
        else:
            if cur:
                segments.append(cur)
            cur = None
    if cur:
        segments.append(cur)
    disordered_segments = [{
        "chain": s["chain"], "start": s["resids"][0], "end": s["resids"][-1],
        "length": len(s["resids"]),
        "mean_plddt": round(sum(s["plddts"]) / len(s["plddts"]), 2),
    } for s in segments]

    # summary
    def counts(seq, alphabet):
        return {k: seq.count(k) for k in alphabet}

    ss3_counts = counts(ss3, SS3)
    ss8_counts = counts(ss8, SS8)
    band_counts = {bd: sum(1 for r in residues if r["band"] == bd) for bd in BANDS}
    plddts = [r["plddt"] for r in residues if r["plddt"] is not None]
    n_conf = sum(1 for r in residues if r["plddt"] is not None and r["plddt"] >= cutoff)
    conf_res = [r for r in residues if r["plddt"] is not None and r["plddt"] >= cutoff]
    ss3_conf = counts("".join(r["ss3"] for r in conf_res), SS3)
    n_masked = sum(1 for r in residues if r["masked"])

    def frac(d, denom):
        return {k: (v / denom if denom else 0.0) for k, v in d.items()}

    summary = {
        "n_residues": n,
        "chains": sorted({r["chain"] for r in residues}),
        "mean_plddt": round(sum(plddts) / len(plddts), 2) if plddts else None,
        "ss3_counts": ss3_counts, "ss3_frac": frac(ss3_counts, n),
        "ss8_counts": ss8_counts,
        "ss3_frac_confident": frac(ss3_conf, n_conf),
        "band_counts": band_counts, "band_frac": frac(band_counts, n),
        "n_confident": n_conf, "frac_confident": (n_conf / n if n else 0.0),
        "n_masked": n_masked, "frac_masked": (n_masked / n if n else 0.0),
        "n_disordered_segments": len(disordered_segments),
        "disordered_segments": disordered_segments,
    }

    return {
        "n_residues": n, "chains": summary["chains"], "plddt_cutoff": cutoff,
        "masking": mask,
        "sequence": sequence, "ss8": ss8, "ss3": ss3,
        "ss8_masked": ss8_masked, "ss3_masked": ss3_masked,
        "residues": residues, "summary": summary,
    }


# --------------------------------------------------------------------------- writers
def _wrap(s, width=60):
    return "\n".join(s[i:i + width] for i in range(0, len(s), width)) or s


def _pl(v, nd=1):
    """Format a pLDDT value, or an em-dash if missing."""
    return "—" if v is None else f"{v:.{nd}f}"


def write_residues_csv(path: Path, struct):
    cols = ["chain", "resid", "resname", "aa", "ss8", "ss3", "plddt", "band",
            "masked", "ss8_masked", "ss3_masked"]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in struct["residues"]:
            w.writerow([r["chain"], r["resid"], r["resname"], r["aa"], r["ss8"], r["ss3"],
                        "" if r["plddt"] is None else f'{r["plddt"]:.2f}',
                        r["band"], r["masked"], r["ss8_masked"], r["ss3_masked"]])


def write_ss_txt(path: Path, stem, struct):
    cut = struct["plddt_cutoff"]
    recs = [
        (f">{stem} sequence", struct["sequence"]),
        (f">{stem} ss8", struct["ss8"]),
        (f">{stem} ss3", struct["ss3"]),
        (f">{stem} ss8_masked plddt<{cut:g}=X", struct["ss8_masked"]),
        (f">{stem} ss3_masked plddt<{cut:g}=X", struct["ss3_masked"]),
    ]
    path.write_text("\n".join(f"{h}\n{_wrap(s)}" for h, s in recs) + "\n")


def write_struct_json(path: Path, stem, struct):
    out = {"stem": stem, **struct}
    path.write_text(json.dumps(out, indent=2))


def write_summary_csv(path: Path, rows):
    cols = ["file", "stem", "status", "n_residues", "n_chains", "pct_H", "pct_E", "pct_C",
            "mean_plddt", "pct_very_low", "pct_low", "pct_confident", "pct_very_high",
            "pct_masked", "n_disordered_segments", "plddt_cutoff"]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            s = r.get("summary")
            if not s:
                w.writerow([r["file"], r["stem"], r["status"]] + [""] * (len(cols) - 3))
                continue
            f3, bf = s["ss3_frac"], s["band_frac"]
            w.writerow([
                r["file"], r["stem"], r["status"], s["n_residues"], len(s["chains"]),
                f'{f3["H"]*100:.1f}', f'{f3["E"]*100:.1f}', f'{f3["C"]*100:.1f}',
                "" if s["mean_plddt"] is None else f'{s["mean_plddt"]:.1f}',
                f'{bf["very_low"]*100:.1f}', f'{bf["low"]*100:.1f}',
                f'{bf["confident"]*100:.1f}', f'{bf["very_high"]*100:.1f}',
                f'{s["frac_masked"]*100:.1f}', s["n_disordered_segments"],
                f'{r["plddt_cutoff"]:g}',
            ])


def write_summary_json(path: Path, rows, meta):
    path.write_text(json.dumps({"metadata": meta, "structures": [
        {k: r[k] for k in ("file", "stem", "status", "plddt_cutoff")} | {"summary": r.get("summary")}
        for r in rows]}, indent=2))


# --------------------------------------------------------------------------- reports
def _md_report(rows, meta):
    out = ["# Secondary structure of AlphaFold structures",
           "",
           f"Generated {meta['generated']} · pLDDT mask cutoff **{meta['cutoff']:g}** · "
           f"engine: mdtraj DSSP (8-state + 3-state).", ""]
    out += ["## Summary", "",
            "| File | Residues | %H | %E | %C | mean pLDDT | % confident (≥cutoff) | % masked | disordered segments |",
            "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for r in rows:
        s = r.get("summary")
        if not s:
            out.append(f"| {r['stem']} | — | | | | | | | _{r['status']}_ |")
            continue
        f3 = s["ss3_frac"]
        out.append(
            f"| {r['stem']} | {s['n_residues']} | {f3['H']*100:.0f} | {f3['E']*100:.0f} | "
            f"{f3['C']*100:.0f} | {_pl(s['mean_plddt'], 0)} | "
            f"{s['frac_confident']*100:.0f} | {s['frac_masked']*100:.0f} | "
            f"{s['n_disordered_segments']} |")
    out.append("")
    for r in rows:
        s = r.get("summary")
        out += [f"## {r['stem']}", ""]
        if not s:
            out += [f"_{r['status']}_", ""]
            continue
        st = r["struct"]
        out += [
            f"- Residues: **{s['n_residues']}** · chains: {', '.join(s['chains'])} · "
            f"mean pLDDT: {_pl(s['mean_plddt'], 1)}",
            f"- 3-state: H {s['ss3_frac']['H']*100:.1f}% · E {s['ss3_frac']['E']*100:.1f}% · "
            f"C {s['ss3_frac']['C']*100:.1f}%",
            f"- Confidence bands: very-low {s['band_frac']['very_low']*100:.1f}% · "
            f"low {s['band_frac']['low']*100:.1f}% · confident {s['band_frac']['confident']*100:.1f}% · "
            f"very-high {s['band_frac']['very_high']*100:.1f}%",
            f"- Masked (pLDDT < {r['plddt_cutoff']:g}): **{s['frac_masked']*100:.1f}%** "
            f"({s['n_masked']} residues) in {s['n_disordered_segments']} disordered segment(s)",
            "",
            "```",
            "sequence   " + st["sequence"],
            "ss8        " + st["ss8"],
            "ss3        " + st["ss3"],
            "ss8_masked " + st["ss8_masked"],
            "```", ""]
        if s["disordered_segments"]:
            out += ["| Disordered segment (chain) | Range | Length | mean pLDDT |", "|---|---|--:|--:|"]
            for seg in s["disordered_segments"]:
                out.append(f"| {seg['chain']} | {seg['start']}–{seg['end']} | "
                           f"{seg['length']} | {seg['mean_plddt']:.1f} |")
            out.append("")
    return "\n".join(out)


_CSS = """
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
     max-width:1100px;margin:2em auto;padding:0 1em;color:#222;line-height:1.5}
h1{border-bottom:2px solid #444;padding-bottom:.3em}
h2{margin-top:2em;border-bottom:1px solid #ccc;padding-bottom:.2em;font-size:1.15em}
table{border-collapse:collapse;margin:.6em 0;font-size:13px}
th,td{border:1px solid #ccc;padding:.35em .55em;text-align:left}
th{background:#f4f4f4}
td.num,th.num{text-align:right}
.muted{color:#888;font-size:12px}
.chip{display:inline-block;min-width:1.1em;text-align:center;color:#fff;border-radius:3px;
      font-size:11px;margin:0 .1em;padding:0 .2em;font-weight:bold}
.legend{margin:.8em 0;font-size:12px}
.block{margin:.5em 0;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
.block .lbl{color:#999;font-size:10px;margin-right:.4em}
.trk{white-space:nowrap;line-height:1.25}
.cell{display:inline-block;width:1.15ch;text-align:center}
.ss{color:#fff;font-weight:bold}
.pos{color:#aaa;font-size:10px}
.bad{color:#b00020;font-weight:bold}
"""


def _legend():
    s8 = " ".join(f'<span class="chip" style="background:{SS8_COLORS[s]}">{s}</span>'
                  f'<span class="muted">{lbl}</span>'
                  for s, lbl in [("H", "α"), ("G", "3₁₀"), ("I", "π"), ("E", "strand"),
                                 ("B", "bridge"), ("T", "turn"), ("S", "bend"), ("C", "coil")])
    sm = f'<span class="chip" style="background:{MASK_COLOR}">X</span><span class="muted">masked</span>'
    bands = " ".join(f'<span class="chip" style="background:{BAND_COLORS[b]}">&nbsp;</span>'
                     f'<span class="muted">{escape(BAND_LABELS[b])}</span>' for b in BANDS)
    return (f'<div class="legend"><b>SS:</b> {s8} {sm}<br>'
            f'<b>pLDDT:</b> {bands}</div>')


def _struct_blocks(struct, width=60):
    res = struct["residues"]
    parts = []
    for off in range(0, len(res), width):
        chunk = res[off:off + width]
        start_resid = chunk[0]["resid"]
        ss = "".join(
            f'<span class="cell ss" style="background:{MASK_COLOR}" title="{r["resid"]} '
            f'pLDDT {r["plddt"]} (masked)">X</span>' if r["masked"] else
            f'<span class="cell ss" style="background:{SS8_COLORS[r["ss8"]]}" '
            f'title="{r["resid"]} {r["ss8"]}/{r["ss3"]} pLDDT {r["plddt"]}">{r["ss8"]}</span>'
            for r in chunk)
        aa = "".join(f'<span class="cell">{r["aa"]}</span>' for r in chunk)
        pl = "".join(f'<span class="cell" style="background:{BAND_COLORS.get(r["band"], "#fff")}" '
                     f'title="pLDDT {r["plddt"]}">&nbsp;</span>' for r in chunk)
        parts.append(
            f'<div class="block"><div><span class="pos">{start_resid}</span></div>'
            f'<div class="trk"><span class="lbl">aa </span>{aa}</div>'
            f'<div class="trk"><span class="lbl">ss </span>{ss}</div>'
            f'<div class="trk"><span class="lbl">pLDDT</span>{pl}</div></div>')
    return "".join(parts)


def _html_report(rows, meta):
    h = ["<!doctype html><html><head><meta charset='utf-8'>",
         "<title>AlphaFold secondary structure</title><style>", _CSS, "</style></head><body>",
         "<h1>Secondary structure of AlphaFold structures</h1>",
         f"<p class='muted'>Generated {escape(meta['generated'])} &middot; pLDDT mask cutoff "
         f"<b>{meta['cutoff']:g}</b> &middot; engine: mdtraj DSSP (8-state + 3-state)</p>",
         _legend(),
         "<h2>Summary</h2>",
         "<table><tr><th>File</th><th class='num'>Residues</th><th class='num'>%H</th>"
         "<th class='num'>%E</th><th class='num'>%C</th><th class='num'>mean pLDDT</th>"
         "<th class='num'>% confident</th><th class='num'>% masked</th>"
         "<th class='num'>disordered segs</th></tr>"]
    for r in rows:
        s = r.get("summary")
        if not s:
            h.append(f"<tr><td>{escape(r['stem'])}</td><td colspan='8' class='bad'>"
                     f"{escape(r['status'])}</td></tr>")
            continue
        f3 = s["ss3_frac"]
        h.append(
            f"<tr><td><a href='#{escape(r['stem'])}'>{escape(r['stem'])}</a></td>"
            f"<td class='num'>{s['n_residues']}</td>"
            f"<td class='num'>{f3['H']*100:.0f}</td><td class='num'>{f3['E']*100:.0f}</td>"
            f"<td class='num'>{f3['C']*100:.0f}</td>"
            f"<td class='num'>{_pl(s['mean_plddt'], 0)}</td>"
            f"<td class='num'>{s['frac_confident']*100:.0f}</td>"
            f"<td class='num'>{s['frac_masked']*100:.0f}</td>"
            f"<td class='num'>{s['n_disordered_segments']}</td></tr>")
    h.append("</table>")
    for r in rows:
        s = r.get("summary")
        h.append(f"<h2 id='{escape(r['stem'])}'>{escape(r['stem'])}</h2>")
        if not s:
            h.append(f"<p class='bad'>{escape(r['status'])}</p>")
            continue
        h.append(
            f"<p class='muted'>{s['n_residues']} residues &middot; chains "
            f"{escape(', '.join(s['chains']))} &middot; mean pLDDT "
            f"{_pl(s['mean_plddt'], 1)} &middot; "
            f"masked (pLDDT&lt;{r['plddt_cutoff']:g}) {s['frac_masked']*100:.1f}% in "
            f"{s['n_disordered_segments']} segment(s)</p>")
        h.append(_struct_blocks(r["struct"]))
        if s["disordered_segments"]:
            h.append("<table><tr><th>Disordered segment</th><th>Chain</th><th>Range</th>"
                     "<th class='num'>Length</th><th class='num'>mean pLDDT</th></tr>")
            for i, seg in enumerate(s["disordered_segments"], 1):
                h.append(f"<tr><td>#{i}</td><td>{escape(seg['chain'])}</td>"
                         f"<td>{seg['start']}–{seg['end']}</td>"
                         f"<td class='num'>{seg['length']}</td>"
                         f"<td class='num'>{seg['mean_plddt']:.1f}</td></tr>")
            h.append("</table>")
    h.append("</body></html>")
    return "".join(h)


# --------------------------------------------------------------------------- input gathering
_EXTS = [".pdb.gz", ".cif.gz", ".pdb", ".pdbx", ".cif", ".mmcif", ".ent", ".bcif"]


def _stem_of(p: Path) -> str:
    nl = p.name.lower()
    for ext in _EXTS:
        if nl.endswith(ext):
            return p.name[: -len(ext)]
    return p.stem


def gather_inputs(inputs, glob_pat):
    pats = [glob_pat] if glob_pat else ["*.pdb", "*.cif", "*.mmcif", "*.ent", "*.pdbx"]
    files = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            for pat in pats:
                files += sorted(p.glob(pat))
        elif p.exists():
            files.append(p)
        else:
            _log(f"WARN: input not found, skipping: {inp}")
    seen, uniq = set(), []
    for f in files:
        rp = f.resolve()
        if rp not in seen:
            seen.add(rp)
            uniq.append(f)
    return uniq


# --------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="One or more .pdb/.cif files, or directories (use --glob).")
    ap.add_argument("-o", "--out-dir", required=True, help="Output directory.")
    ap.add_argument("--plddt-cutoff", type=float, default=70.0,
                    help="Mask SS where pLDDT < cutoff (default 70). The 4 AF bands are "
                         "always reported regardless.")
    ap.add_argument("--no-mask", action="store_true",
                    help="Do not mask low-pLDDT positions (still report pLDDT + bands).")
    ap.add_argument("--glob", default=None,
                    help="Glob for directory inputs (default: *.pdb,*.cif,*.mmcif,*.ent,*.pdbx).")
    ap.add_argument("-j", "--workers", type=int, default=8, help="Parallel workers (default 8).")
    args = ap.parse_args(argv)

    files = gather_inputs(args.inputs, args.glob)
    if not files:
        sys.exit("No input structures found.")
    out_dir = Path(args.out_dir)
    per_dir = out_dir / "per_structure"
    per_dir.mkdir(parents=True, exist_ok=True)

    # assign unique stems (disambiguate collisions, e.g. same name .pdb and .cif)
    used = {}
    jobs = []
    for f in files:
        stem = _stem_of(f)
        if stem in used:
            used[stem] += 1
            stem = f"{stem}__{f.suffix.lstrip('.') or used[stem]}"
        else:
            used[stem] = 1
        jobs.append((f, stem))

    cutoff = args.plddt_cutoff
    mask = not args.no_mask
    _log(f"analyzing {len(jobs)} structure(s); pLDDT cutoff {cutoff:g}, masking {'on' if mask else 'off'}")

    def work(job):
        f, stem = job
        rec = {"file": f.name, "stem": stem, "plddt_cutoff": cutoff}
        try:
            struct = analyze_structure(f, cutoff, mask)
            if struct["n_residues"] == 0:
                rec["status"] = "error: no protein residues"
                rec["summary"] = None
                return rec
            rec["struct"] = struct
            rec["summary"] = struct["summary"]
            rec["status"] = "ok"
            write_residues_csv(per_dir / f"{stem}.residues.csv", struct)
            write_ss_txt(per_dir / f"{stem}.ss.txt", stem, struct)
            write_struct_json(per_dir / f"{stem}.json", stem, struct)
        except Exception as e:                                   # noqa: BLE001
            _log(f"ERROR {f.name}: {e}")
            rec["status"] = f"error: {e}"
            rec["summary"] = None
        return rec

    if args.workers <= 1 or len(jobs) == 1:
        rows = [work(j) for j in jobs]
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            rows = list(pool.map(work, jobs))

    meta = {"generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
            "cutoff": cutoff, "masking": mask, "n_structures": len(rows)}
    write_summary_csv(out_dir / "summary.csv", rows)
    write_summary_json(out_dir / "summary.json", rows, meta)
    (out_dir / "report.md").write_text(_md_report(rows, meta))
    (out_dir / "report.html").write_text(_html_report(rows, meta))

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    _log(f"done: {n_ok}/{len(rows)} ok; outputs in {out_dir}")
    return 0 if n_ok else 1


if __name__ == "__main__":
    sys.exit(main())
