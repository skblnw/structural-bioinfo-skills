---
name: af-secondary-structure
description: Assign per-residue SECONDARY STRUCTURE to an AlphaFold structure (or any PDB/mmCIF) — the 8-state DSSP code (H,G,I,E,B,T,S,C) and the reduced 3-state code (H/E/C) for every sequence position — with built-in pLDDT confidence filtering so unreliable disordered regions can be flagged or masked. Use this skill whenever the user has an AlphaFold model (a .pdb or .cif, e.g. AF-XXXXXX-F1) and wants its secondary structure per residue, asks "what's the DSSP / fold / helix-sheet-coil of this predicted structure", "assign secondary structure to my AlphaFold model", "which regions are reliable helices/strands vs disordered loops", "give me an 8-state or 3-state SS string", "mask the low-confidence parts", or wants to find/trim low-pLDDT disordered segments. Reads pLDDT from the B-factor column, classifies the 4 AlphaFold confidence bands, and emits a masked SS string (low-pLDDT → X) alongside the raw call. Handles a single file or a whole batch/directory. Outputs a per-residue CSV, FASTA-like SS strings, a JSON dump, and self-contained Markdown + HTML reports. For mapping epitope peptides onto secondary structure aggregated across EXPERIMENTAL PDB depositions, use the sibling `epitope-secondary-structure` skill instead — this skill is for assigning SS directly on a (predicted) structure file.
---

# AlphaFold secondary structure (DSSP + pLDDT QC)

Given **a structure file** — typically an AlphaFold model — assign, **for each sequence
position**, the 8-state DSSP secondary structure and the reduced 3-state code, then layer on
AlphaFold quality control. AlphaFold models routinely contain long, low-confidence regions that
are *predicted* to be coil/extended but are really just disordered and unreliable; this skill keeps
the raw assignment **and** produces a confidence-filtered ("masked") version so those regions can be
flagged or removed downstream without losing information.

## When to use this skill

- "Assign secondary structure to this AlphaFold model / predicted structure."
- "What's the per-residue DSSP (8-state) of `AF-P04637-F1.pdb`? Give me a 3-state string too."
- "Which parts of this prediction are reliable helices/strands and which are disordered?"
- "Mask out the low-pLDDT regions of the secondary structure" / "find the disordered segments."
- "Run secondary-structure + confidence QC over this folder of AlphaFold models."

The unit of analysis is **one structure → one per-residue table**. Batch mode just repeats it over
many files and adds a combined summary.

**Not this skill:** to ask "what secondary structure does this *epitope peptide* adopt in its parent
protein, pooled over the parent's experimentally solved PDB structures," use
`epitope-secondary-structure` (epitope-centric, experimental-first, mkdssp). This skill assigns SS
directly on a given structure file and is the right tool for AlphaFold models specifically.

## Prerequisites

The DSSP engine is **mdtraj** (a self-contained reimplementation — no external `mkdssp` binary, and
unlike mkdssp 3.x it parses AlphaFold mmCIF fine). It must be importable from some conda env. The
script **auto-discovers** an env that has mdtraj (preferring one named `mdanalysis`, then `mdtraj`)
and re-execs itself there, so you can launch it with any Python. If none is found it prints an
install hint and exits:

```bash
conda create -n mdanalysis -c conda-forge mdtraj      # or: conda install -n <env> -c conda-forge mdtraj
```

## How to invoke

```bash
S=~/.claude/plugins/structural-bioinfo/skills/af-secondary-structure/scripts/af_ss.py

# single AlphaFold model (PDB or mmCIF)
python "$S" AF-P04637-F1.pdb -o p53_ss/

# a whole directory of models (batch)
python "$S" ./af_models/ --glob "*.pdb" -o af_models_ss/

# several explicit files; keep SS for the 50–70 band too (mask only < 50)
python "$S" modelA.cif modelB.cif --plddt-cutoff 50 -o out/

# annotate confidence but never alter the SS letters
python "$S" model.pdb --no-mask -o out/
```

### Options

| Flag | Default | Purpose |
|---|---|---|
| `inputs` (positional) | (required) | One or more `.pdb`/`.cif`/`.mmcif` files, or directories. |
| `--out-dir`, `-o` | (required) | Output directory. |
| `--plddt-cutoff` | `70` | Mask SS where pLDDT < cutoff. The 4 AF bands are reported regardless. |
| `--no-mask` | off | Report pLDDT + bands but never replace SS with `X`. |
| `--glob` | `*.pdb,*.cif,*.mmcif,*.ent,*.pdbx` | Glob applied to directory inputs. |
| `--workers`, `-j` | `8` | Parallel workers across input files. |

## What it does (per structure)

1. **Load** the structure with mdtraj (PDB or mmCIF; first model of multi-model files).
2. **DSSP** twice: `simplified=False` → 8-state, `simplified=True` → 3-state. These align 1:1 with
   the structure's protein residues; non-protein residues (ligands, ions, water) are dropped.
3. **pLDDT** is read from the **B-factor column** (the AlphaFold convention) per residue and matched
   by `(chain, residue number)`. Each residue is binned into the standard AlphaFold bands:
   very-low `<50`, low `50–70`, confident `70–90`, very-high `≥90`.
4. **Mask**: a copy of the SS string where every position with pLDDT `< cutoff` becomes `X`
   (`--no-mask` disables this; `--plddt-cutoff` changes the threshold).
5. **Disordered segments**: contiguous runs (within a chain) of pLDDT `< cutoff` are reported as
   ranges with their length and mean pLDDT — the "unreliable disordered regions."
6. **Summary**: SS composition (overall and among confident residues), the band distribution, mean
   pLDDT, percent masked, and the segment list.

The 8→3 collapse is the standard `{H,G,I}→H, {E,B}→E, {T,S,C}→C`. See
`references/method_notes.md` for the full state/band tables, engine rationale, and edge cases.

## Output structure

```
<out>/
├── report.html              # self-contained: summary table + per-residue SS/pLDDT strip + legend
├── report.md                # same content as GitHub-flavored markdown (SS strings in code fences)
├── summary.csv              # one row per structure (SS %, band %, % masked, # disordered segments)
├── summary.json             # metadata + per-structure summaries
└── per_structure/
    ├── <stem>.residues.csv  # chain,resid,resname,aa,ss8,ss3,plddt,band,masked,ss8_masked,ss3_masked
    ├── <stem>.ss.txt        # FASTA-like: sequence / ss8 / ss3 / ss8_masked / ss3_masked
    └── <stem>.json          # full per-residue dump + disordered_segments + summary
```

`<stem>` is the input filename without extension. If a batch contains two inputs with the same stem
(e.g. a `.pdb` and a `.cif` of the same model), the second is disambiguated with a suffix.

### Secondary-structure alphabets

- **8-state** (`ss8`): `H` α-helix · `G` 3₁₀-helix · `I` π-helix · `E` β-strand · `B` β-bridge ·
  `T` turn · `S` bend · `C` coil/loop.
- **3-state** (`ss3`): `H` helix · `E` strand · `C` coil.
- **Masked** (`ss8_masked` / `ss3_masked`): same as raw, but `X` wherever pLDDT < cutoff.

## Reading the results

- Lead with the headline per structure: dominant SS composition **among confident residues**
  (`ss3_frac_confident`), mean pLDDT, and how much was masked. A model that is "70% coil" overall is
  usually just disordered — the confident-only composition is the honest fold description.
- Point the user at the disordered-segment ranges: these are exactly the regions to trim before
  docking/MD/analysis. Don't paste the HTML into chat; point to `report.html`.
- The masked SS string is the safe input for downstream consumers that shouldn't trust AF coil in
  low-pLDDT stretches.

## Gotchas

- **pLDDT lives in B-factors only for predictions.** On an experimental PDB the B-factor column is a
  crystallographic B-factor, *not* pLDDT — the band/mask columns are then meaningless. This skill is
  meant for AlphaFold (or other predictors that write pLDDT to B-factors). It still computes correct
  SS on experimental structures; just ignore the confidence columns there (or use `--no-mask`).
- **Engine differences.** mdtraj's DSSP is a faithful Kabsch–Sander reimplementation but can differ
  slightly from `mkdssp` 4 in rare states (it has no PPII `P` state and can be liberal with π-helix
  `I`). For per-residue helix/strand/coil this is immaterial.
- **Multi-chain / multimer** models work; chain breaks are respected and disordered segments never
  bridge chains. Coordinates are taken from the first model of multi-model files.
- A residue whose pLDDT can't be matched (no CA B-factor) gets band `unknown` and is left unmasked.

## Dependencies

Python 3.8+ and **mdtraj** (in a conda env; `conda install -c conda-forge mdtraj`). No other
third-party packages; everything else is stdlib. No network access.
