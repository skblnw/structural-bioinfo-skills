#!/usr/bin/env python3
"""
Featurize protein sequences with ESM Cambrian (ESM C 300M).

One sequence in -> one .npz out. Contains per-residue embedding (L, 960),
mean-pooled embedding (960,), the input sequence, the id, and the model name.

Usage:
    python featurize.py --sequence MPRTKEINDAGLIVHSPQWFYK --id myprot
    python featurize.py --fasta seqs.fasta --out features/
    python featurize.py --list seqs.txt --out features/
"""

import argparse
import datetime
import glob
import html
import os
import re
import sys
import time

# ── esm bootstrap ────────────────────────────────────────────────────
# The `esm` package is rarely on the system Python. If we can't import it
# here, look for a conda env named `esmc` (preferred) or any sibling env that
# has `esm`, and re-exec the script under that env's Python. A loop guard
# prevents infinite re-exec.

def _candidate_pythons():
    home = os.path.expanduser("~")
    roots = [
        f"{home}/miniconda3", f"{home}/anaconda3",
        f"{home}/miniforge3", f"{home}/mambaforge", f"{home}/conda",
        "/opt/conda", "/opt/miniconda3", "/opt/anaconda3",
        "/opt/homebrew/Caskroom/miniconda/base",
    ]
    cp = os.environ.get("CONDA_PREFIX", "")
    if cp:
        if os.path.basename(os.path.dirname(cp)) == "envs":
            roots.insert(0, os.path.dirname(os.path.dirname(cp)))
        else:
            roots.insert(0, cp)
    seen = set()
    # Prefer an env literally named "esmc"
    for root in roots:
        py = os.path.join(root, "envs", "esmc", "bin", "python")
        if py not in seen and os.path.isfile(py):
            seen.add(py); yield py
    # Then any sibling env
    for root in roots:
        for env in sorted(glob.glob(os.path.join(root, "envs", "*"))):
            py = os.path.join(env, "bin", "python")
            if py not in seen and os.path.isfile(py):
                seen.add(py); yield py


def _bootstrap_esm():
    try:
        import esm  # noqa: F401
        return
    except ImportError:
        pass
    if os.environ.get("ESM_FEATURIZE_REEXECED"):
        sys.stderr.write(
            "[error] esm is not importable even after auto-discovery.\n"
            "        Run: conda create -n esmc python=3.11 -y \\\n"
            "             && conda activate esmc \\\n"
            "             && pip install esm numpy\n"
        )
        sys.exit(2)
    import subprocess
    for py in _candidate_pythons():
        try:
            r = subprocess.run([py, "-c", "import esm"],
                               capture_output=True, timeout=15)
        except (subprocess.TimeoutExpired, OSError):
            continue
        if r.returncode == 0:
            sys.stderr.write(f"[info] Re-executing under {py} (auto-discovered esm env)\n")
            os.environ["ESM_FEATURIZE_REEXECED"] = "1"
            os.execv(py, [py] + sys.argv)
    sys.stderr.write(
        "[error] esm is not importable and no conda env with esm was found.\n"
        "        Run: conda create -n esmc python=3.11 -y \\\n"
        "             && conda activate esmc \\\n"
        "             && pip install esm numpy\n"
    )
    sys.exit(2)


_bootstrap_esm()

import numpy as np
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein


VALID_AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


# ── Input parsing ────────────────────────────────────────────────────

def parse_fasta(path):
    records = []
    cur_id, cur_seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(cur_seq)))
                cur_id = line[1:].split()[0] if line[1:].strip() else f"seq_{len(records)+1:04d}"
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id is not None:
            records.append((cur_id, "".join(cur_seq)))
    if not records:
        sys.stderr.write(f"[error] no records found in FASTA: {path}\n")
        sys.exit(1)
    return records


def parse_list(path):
    records = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            records.append((f"seq_{len(records)+1:04d}", s))
    if not records:
        sys.stderr.write(f"[error] no sequences found in list: {path}\n")
        sys.exit(1)
    return records


def validate(records, max_len):
    cleaned = []
    for sid, seq in records:
        seq = seq.upper().replace(" ", "").replace("\t", "")
        if not seq:
            sys.stderr.write(f"[error] {sid}: empty sequence\n")
            sys.exit(1)
        if not VALID_AA.match(seq):
            bad = sorted({c for c in seq if c not in "ACDEFGHIKLMNPQRSTVWY"})
            sys.stderr.write(
                f"[error] {sid}: sequence contains non-canonical residues {bad}. "
                f"Only the 20 canonical amino acids (ACDEFGHIKLMNPQRSTVWY) are accepted.\n"
            )
            sys.exit(1)
        if len(seq) > max_len:
            sys.stderr.write(
                f"[warn] {sid}: length {len(seq)} > --max-len {max_len}, proceeding anyway\n"
            )
        cleaned.append((sid, seq))
    return cleaned


# ── Device ───────────────────────────────────────────────────────────

def pick_device(flag):
    if flag in ("cpu", "mps"):
        return flag
    return "mps" if torch.backends.mps.is_available() else "cpu"


# ── Embedding ────────────────────────────────────────────────────────

def extract_embeddings(out_obj):
    """Pull the per-token embedding tensor out of an ESMC forward output,
    regardless of whether the SDK returned a dataclass, a namedtuple, or a
    raw tensor."""
    for attr in ("embeddings", "hidden_states", "last_hidden_state"):
        if hasattr(out_obj, attr):
            t = getattr(out_obj, attr)
            if torch.is_tensor(t):
                return t
    if torch.is_tensor(out_obj):
        return out_obj
    if isinstance(out_obj, (tuple, list)) and out_obj and torch.is_tensor(out_obj[0]):
        return out_obj[0]
    raise RuntimeError(
        f"Could not locate embedding tensor in ESMC forward output of type "
        f"{type(out_obj).__name__}; attributes: {dir(out_obj)}"
    )


def embed_one(model, sequence, device):
    protein = ESMProtein(sequence=sequence)
    protein_tensor = model.encode(protein)        # ESMProteinTensor with .sequence (1D)
    tokens = protein_tensor.sequence.unsqueeze(0) # (1, L+special)
    with torch.inference_mode():
        out = model.forward(sequence_tokens=tokens)
    emb = extract_embeddings(out)                 # (1, L+special, H)
    emb = emb.squeeze(0).to(torch.float32).cpu().numpy()

    L = len(sequence)
    if emb.shape[0] == L + 2:
        emb = emb[1:-1]                         # strip BOS/EOS
    elif emb.shape[0] == L + 1:
        emb = emb[1:]                           # strip one leading special
    elif emb.shape[0] != L:
        raise RuntimeError(
            f"Unexpected embedding length {emb.shape[0]} for sequence length {L}; "
            f"cannot reliably strip special tokens."
        )
    return emb                                  # (L, H) float32


# ── Main ─────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Featurize protein sequences with ESM Cambrian 300M."
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--sequence", type=str, help="inline amino-acid sequence")
    src.add_argument("--fasta", type=str, help="path to a FASTA file")
    src.add_argument("--list", dest="list_path", type=str,
                     help="path to a text file, one sequence per line")
    p.add_argument("--id", type=str, default="query",
                   help="id for --sequence mode (default: query)")
    p.add_argument("--out", type=str, default="./esmc_features",
                   help="output directory (default: ./esmc_features)")
    p.add_argument("--device", choices=("auto", "cpu", "mps"), default="auto",
                   help="device (default: auto)")
    p.add_argument("--model", type=str, default="esmc_300m",
                   help="ESM C model id (default: esmc_300m)")
    p.add_argument("--max-len", type=int, default=2048,
                   help="warn if a sequence exceeds this length (default: 2048)")
    p.add_argument("--overwrite", action="store_true",
                   help="re-run even if the output .npz already exists")
    args = p.parse_args()

    if args.sequence is not None:
        records = [(args.id, args.sequence)]
    elif args.fasta is not None:
        records = parse_fasta(args.fasta)
    else:
        records = parse_list(args.list_path)

    records = validate(records, args.max_len)
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)

    device = pick_device(args.device)
    sys.stderr.write(f"[info] device={device}, model={args.model}\n")
    sys.stderr.write(f"[info] loading {args.model} (first run downloads ~600 MB)...\n")
    t0 = time.time()
    model = ESMC.from_pretrained(args.model).to(device).eval()
    sys.stderr.write(f"[info] model loaded in {time.time()-t0:.1f}s\n")

    t_start = time.time()
    n_done, n_skip = 0, 0
    entries = []      # one dict per sequence (written or skipped), for the report
    pooled_vecs = []  # mean_pooled vectors of *processed* sequences, in order
    pooled_ids = []
    for sid, seq in records:
        out_path = os.path.join(out_dir, f"{sid}.npz")
        if os.path.exists(out_path) and not args.overwrite:
            sys.stderr.write(f"[skip] {sid} (exists; use --overwrite to redo)\n")
            n_skip += 1
            with np.load(out_path) as d:
                mp = d["mean_pooled"]
                pr_shape = d["per_residue"].shape
            entries.append({
                "id": sid, "length": len(seq), "npz": f"{sid}.npz",
                "status": "skipped",
                "mean_norm": float(np.linalg.norm(mp)),
                "per_residue_shape": pr_shape,
            })
            pooled_vecs.append(mp); pooled_ids.append(sid)
            continue
        sys.stderr.write(f"[run ] {sid} (L={len(seq)})...\n")
        per_residue = embed_one(model, seq, device).astype(np.float32)
        mean_pooled = per_residue.mean(axis=0).astype(np.float32)
        np.savez_compressed(
            out_path,
            per_residue=per_residue,
            mean_pooled=mean_pooled,
            sequence=np.array(seq),
            id=np.array(sid),
            model=np.array(args.model),
        )
        n_done += 1
        entries.append({
            "id": sid, "length": len(seq), "npz": f"{sid}.npz",
            "status": "written",
            "mean_norm": float(np.linalg.norm(mean_pooled)),
            "per_residue_shape": per_residue.shape,
        })
        pooled_vecs.append(mean_pooled); pooled_ids.append(sid)

    elapsed = time.time() - t_start
    hidden = pooled_vecs[0].shape[0] if pooled_vecs else 0
    sys.stderr.write(
        f"[done] {n_done} written, {n_skip} skipped -> {out_dir} "
        f"(hidden_dim={hidden}, {elapsed:.1f}s on {device})\n"
    )

    # ── Report ────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sim = None
    if len(pooled_vecs) >= 2:
        V = np.stack(pooled_vecs).astype(np.float32)
        n = np.linalg.norm(V, axis=1, keepdims=True)
        n[n == 0] = 1.0
        sim = (V / n) @ (V / n).T

    md_path = os.path.join(out_dir, "report.md")
    html_path = os.path.join(out_dir, "report.html")
    write_markdown_report(md_path, args, device, hidden, elapsed,
                          n_done, n_skip, entries, sim, pooled_ids, timestamp)
    write_html_report(html_path, args, device, hidden, elapsed,
                      n_done, n_skip, entries, sim, pooled_ids, timestamp)
    sys.stderr.write(f"[done] report: {md_path}\n")
    sys.stderr.write(f"[done] report: {html_path}\n")


# ── Reports ──────────────────────────────────────────────────────────

def _sim_table_rows(sim, ids, fmt):
    """Yield (id, [(other_id, value)]) for at most 20 rows/cols."""
    k = min(20, len(ids))
    truncated = len(ids) > k
    for i in range(k):
        yield ids[i], [(ids[j], fmt(sim[i, j])) for j in range(k)]
    return truncated


def write_markdown_report(path, args, device, hidden, elapsed,
                          n_done, n_skip, entries, sim, ids, timestamp):
    lines = []
    lines.append(f"# ESM-Featurize Report")
    lines.append("")
    lines.append(f"- **Timestamp:** {timestamp}")
    lines.append(f"- **Model:** `{args.model}`  (hidden dim = {hidden})")
    lines.append(f"- **Device:** `{device}`")
    lines.append(f"- **Output directory:** `{os.path.abspath(args.out)}`")
    lines.append(f"- **Written:** {n_done}    **Skipped:** {n_skip}    "
                 f"**Elapsed:** {elapsed:.2f}s")
    lines.append("")
    lines.append("## Sequences")
    lines.append("")
    lines.append("| # | id | length | status | npz | ‖mean_pooled‖₂ |")
    lines.append("|---:|----|---:|----|----|---:|")
    for i, e in enumerate(entries, 1):
        lines.append(f"| {i} | `{e['id']}` | {e['length']} | {e['status']} | "
                     f"`{e['npz']}` | {e['mean_norm']:.3f} |")
    lines.append("")

    if sim is not None:
        k = min(20, len(ids))
        lines.append("## Cosine similarity (mean-pooled embeddings)")
        lines.append("")
        if len(ids) > k:
            lines.append(f"_Showing first {k} of {len(ids)} sequences._")
            lines.append("")
        header = "| | " + " | ".join(f"`{ids[j]}`" for j in range(k)) + " |"
        sep = "|---|" + "|".join(["---:"] * k) + "|"
        lines.append(header)
        lines.append(sep)
        for i in range(k):
            row = f"| `{ids[i]}` | " + " | ".join(f"{sim[i,j]:.3f}" for j in range(k)) + " |"
            lines.append(row)
        lines.append("")

    lines.append("## How to load")
    lines.append("")
    lines.append("```python")
    lines.append("import numpy as np")
    lines.append(f"d = np.load(\"{entries[0]['npz'] if entries else 'out.npz'}\")")
    lines.append("mean_vec = d['mean_pooled']   # (H,) for similarity / classifiers")
    lines.append("per_res  = d['per_residue']   # (L, H) for residue-level work")
    lines.append("```")
    lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))


def write_html_report(path, args, device, hidden, elapsed,
                      n_done, n_skip, entries, sim, ids, timestamp):
    esc = html.escape
    css = """
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
           max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #222; }
    h1, h2 { border-bottom: 1px solid #eee; padding-bottom: 0.3em; }
    table { border-collapse: collapse; margin: 1em 0; font-size: 0.92em; }
    th, td { border: 1px solid #ddd; padding: 4px 8px; }
    th { background: #f6f8fa; text-align: left; }
    td.num { text-align: right; font-variant-numeric: tabular-nums; }
    code, pre { font-family: 'SF Mono', Menlo, Consolas, monospace; font-size: 0.9em; }
    pre { background: #f6f8fa; padding: 0.8em 1em; border-radius: 4px; overflow-x: auto; }
    .meta li { margin: 0.15em 0; }
    .sim td { background: #fff; }
    """
    parts = ["<!doctype html><html><head><meta charset='utf-8'>",
             "<title>ESM-Featurize Report</title>",
             f"<style>{css}</style></head><body>"]
    parts.append("<h1>ESM-Featurize Report</h1>")
    parts.append("<ul class='meta'>")
    parts.append(f"<li><b>Timestamp:</b> {esc(timestamp)}</li>")
    parts.append(f"<li><b>Model:</b> <code>{esc(args.model)}</code> "
                 f"(hidden dim = {hidden})</li>")
    parts.append(f"<li><b>Device:</b> <code>{esc(device)}</code></li>")
    parts.append(f"<li><b>Output directory:</b> <code>{esc(os.path.abspath(args.out))}</code></li>")
    parts.append(f"<li><b>Written:</b> {n_done} &nbsp; <b>Skipped:</b> {n_skip} "
                 f"&nbsp; <b>Elapsed:</b> {elapsed:.2f}s</li>")
    parts.append("</ul>")

    parts.append("<h2>Sequences</h2>")
    parts.append("<table><thead><tr>"
                 "<th>#</th><th>id</th><th>length</th><th>status</th>"
                 "<th>npz</th><th>&#x2016;mean_pooled&#x2016;<sub>2</sub></th>"
                 "</tr></thead><tbody>")
    for i, e in enumerate(entries, 1):
        parts.append(
            f"<tr><td class='num'>{i}</td>"
            f"<td><code>{esc(e['id'])}</code></td>"
            f"<td class='num'>{e['length']}</td>"
            f"<td>{esc(e['status'])}</td>"
            f"<td><code>{esc(e['npz'])}</code></td>"
            f"<td class='num'>{e['mean_norm']:.3f}</td></tr>"
        )
    parts.append("</tbody></table>")

    if sim is not None:
        k = min(20, len(ids))
        parts.append("<h2>Cosine similarity (mean-pooled embeddings)</h2>")
        if len(ids) > k:
            parts.append(f"<p><em>Showing first {k} of {len(ids)} sequences.</em></p>")
        parts.append("<table class='sim'><thead><tr><th></th>")
        for j in range(k):
            parts.append(f"<th><code>{esc(ids[j])}</code></th>")
        parts.append("</tr></thead><tbody>")
        for i in range(k):
            parts.append(f"<tr><th><code>{esc(ids[i])}</code></th>")
            for j in range(k):
                v = sim[i, j]
                # color cell by similarity: 1.0 -> green, 0.0 -> white, <0 -> red
                if v >= 0:
                    g = int(255 - 120 * v); r = 255; b = int(255 - 120 * v)
                else:
                    r = 255; g = int(255 + 120 * v); b = int(255 + 120 * v)
                parts.append(f"<td class='num' style='background:rgb({r},{g},{b})'>{v:.3f}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table>")

    sample = entries[0]['npz'] if entries else "out.npz"
    parts.append("<h2>How to load</h2>")
    parts.append("<pre>import numpy as np\n"
                 f"d = np.load(\"{esc(sample)}\")\n"
                 "mean_vec = d['mean_pooled']   # (H,) for similarity / classifiers\n"
                 "per_res  = d['per_residue']   # (L, H) for residue-level work</pre>")
    parts.append("</body></html>")

    with open(path, "w") as f:
        f.write("\n".join(parts))


if __name__ == "__main__":
    main()
