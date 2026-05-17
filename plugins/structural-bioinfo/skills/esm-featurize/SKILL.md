---
name: esm-featurize
description: Use this skill when the user asks to "embed a protein sequence with ESM C", "compute ESM Cambrian embeddings", "featurize sequences with esmc_300m", "get protein representations for FASTA", "esmc embeddings from a query sequence", "vectorize protein sequences", "generate per-residue ESM embeddings", "mean-pooled protein embedding", or supplies one or more amino-acid sequences (single string, FASTA file, or one-per-line list) and wants a numeric representation written to disk. Triggers on phrases like "embed this sequence with ESM", "ESM C 300M representation", "protein language model features for these sequences", "save ESM embeddings as npz", "featurize my FASTA". This skill is REPRESENTATION-ONLY — it does NOT do generation, masking, structure prediction, inverse folding, fine-tuning, attention extraction, or call any cloud/Forge API. Use the broader scientific-agent-skills:esm skill for any of those. Outputs one compressed `.npz` per sequence containing per-residue embeddings (L×960 float32), mean-pooled sequence embedding (960 float32), the input sequence, the id, and the model name.
---

# ESM-Featurize

## Purpose

Turn one or more protein sequences into ESM Cambrian (ESM C 300M) embeddings,
written to disk as compressed NumPy archives. One sequence in → one `.npz`
out. Nothing more.

## When to Use

- "Embed this sequence with ESM C: MPRTKEINDAGLIVHSP…"
- "Featurize my FASTA with esmc_300m"
- "Get protein representations for these 50 sequences"
- "I need per-residue embeddings for downstream clustering"
- "Mean-pooled ESM embedding for similarity search"

Do **not** use this skill for: generation, masking/completion, structure
prediction, inverse folding, logits, attention weights, fine-tuning, or
Forge API calls. The broader `scientific-agent-skills:esm` skill covers
those.

## Prerequisites

A conda env named `esmc` with the `esm` package installed:

```bash
conda create -n esmc python=3.11 -y
conda activate esmc
pip install esm numpy httpx
```

`httpx` isn't a declared dependency of `esm` but the package imports it
unconditionally at top level — installing it explicitly avoids an
`ImportError` on first use.

First run will download the esmc_300m weights (~600 MB) into
`~/.cache/huggingface/`. After that the env runs offline.

The script auto-discovers the `esmc` env: you can invoke it with any Python
(including the system Python) and it will re-exec itself under the right
interpreter.

## Quick Start

> **Script location:** `scripts/featurize.py` lives in the skill directory.
> Use the full path or `cd` to the skill base before running.

### Single sequence

```bash
python scripts/featurize.py --sequence MPRTKEINDAGLIVHSPQWFYK --id myprot
# writes ./esmc_features/myprot.npz
```

### FASTA file (one .npz per record)

```bash
python scripts/featurize.py --fasta seqs.fasta --out features/
# writes features/<header1>.npz, features/<header2>.npz, ...
```

### Plain list (one sequence per line)

```bash
python scripts/featurize.py --list seqs.txt --out features/
# writes features/seq_0001.npz, features/seq_0002.npz, ...
```

## Flags

| Flag | Default | Meaning |
|------|---------|---------|
| `--sequence STR` | — | Inline single sequence (mutually exclusive with `--fasta`, `--list`) |
| `--fasta PATH`   | — | FASTA input; ids come from the first whitespace token of each header |
| `--list PATH`    | — | One sequence per non-blank, non-`#` line; ids auto-numbered `seq_0001`… |
| `--id NAME`      | `query` | Override id for `--sequence` mode |
| `--out DIR`      | `./esmc_features/` | Output directory (created if missing) |
| `--device {auto,cpu,mps}` | `auto` | `auto` picks `mps` if available, else `cpu` |
| `--model NAME`   | `esmc_300m` | Exposed for future swap to `esmc-600m`; 300M is the default |
| `--max-len N`    | `2048` | Warn (don't error) if a sequence exceeds N residues |
| `--overwrite`    | off | Re-run even if the output `.npz` already exists |

## Output Schema

Each `.npz` contains exactly these keys:

| key           | dtype     | shape       | meaning |
|---------------|-----------|-------------|---------|
| `per_residue` | float32   | `(L, 960)`  | per-position embedding, special tokens stripped so length matches the input sequence |
| `mean_pooled` | float32   | `(960,)`    | mean over the L residues |
| `sequence`    | `<U`      | scalar str  | the input AA string |
| `id`          | `<U`      | scalar str  | FASTA header / list index / `--id` |
| `model`       | `<U`      | scalar str  | `esmc_300m` |

Load with:

```python
import numpy as np
d = np.load("esmc_features/myprot.npz")
mean_vec   = d["mean_pooled"]    # (960,) for similarity / classifiers
per_res    = d["per_residue"]    # (L, 960) for residue-level work
```

## Reports

At the end of every run, the script also writes:

- `<out>/report.md`   — GitHub-flavored Markdown summary
- `<out>/report.html` — self-contained HTML (no JS, no CDN, opens offline)

Each report contains:

- Run metadata (model, device, timestamp, elapsed, written / skipped counts)
- One-row-per-sequence table: id, length, npz path, ‖mean_pooled‖₂
- When N ≥ 2 sequences are processed: a cosine-similarity matrix between the
  mean-pooled embeddings (HTML version is heat-mapped; truncated to 20×20 if
  more sequences were run)
- A copy-paste snippet showing how to `np.load` an `.npz`

Reports overwrite on every run.

## Sequence Validation

Sequences are uppercased and must contain only the canonical 20 amino acids:
`ACDEFGHIKLMNPQRSTVWY`. Ambiguity codes (`B/Z/X/U/O/*`) and gaps (`-`) are
rejected with an error naming the offending id. Empty sequences are rejected.
Sequences longer than `--max-len` produce a warning but proceed.
