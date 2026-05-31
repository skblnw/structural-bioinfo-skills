# Method notes for `af-secondary-structure`

How a structure file is turned into per-residue secondary structure + pLDDT confidence.

## 1. Why mdtraj (and not mkdssp / biotite / MDAnalysis)

The skill needs **8-state DSSP** on **AlphaFold** files, with **no manual setup**.

- **mdtraj `compute_dssp`** — a self-contained reimplementation of the Kabsch–Sander DSSP
  algorithm. `simplified=False` returns the 8 states, `simplified=True` returns 3. It reads PDB and
  mmCIF (`pdb`/`pdbx`) directly, needs no external binary, and parses AlphaFold mmCIF without issue.
  **Chosen.**
- **`mkdssp`** (the classic DSSP binary) — gold standard, but (a) often not installed, and (b)
  the common `mkdssp` 3.1.2 **cannot parse AlphaFold mmCIF** ("empty protein…"); AF needs DSSP 4.
  The sibling `epitope-secondary-structure` skill uses mkdssp on *experimental* structures for
  exactly that reason; this skill avoids the binary so it works out-of-the-box on AF.
- **`MDAnalysis.analysis.dssp`** and **`biotite.structure.annotate_sse`** — both give **3-state
  only** (MDAnalysis ports pydssp; biotite uses the CA-only P-SEA geometry). Insufficient when the
  user wants 8-state.

mdtraj's DSSP can differ slightly from `mkdssp` 4 in rare cases: it has no PPII `P` state, and it can
assign π-helix (`I`) a little more liberally. For helix/strand/coil this is immaterial.

## 2. Secondary-structure states

8-state (DSSP), with mdtraj's blank/loop mapped to `C`:

| Code | Meaning |
|---|---|
| `H` | α-helix |
| `G` | 3₁₀-helix |
| `I` | π-helix |
| `E` | extended β-strand (in a ladder) |
| `B` | isolated β-bridge |
| `T` | hydrogen-bonded turn |
| `S` | bend |
| `C` | coil / loop / irregular (mdtraj emits a blank here) |

### The 8 → 3 collapse (standard)

```
H, G, I  -> H   (helix)
E, B     -> E   (strand)
T, S, C  -> C   (coil; DSSP blank / '-' / DSSP4 'P' all fold here)
```

This matches the collapse used by the sibling `epitope-secondary-structure` skill, so 3-state calls
are consistent across the plugin. The 3-state string is taken from mdtraj's own `simplified=True`
output (which applies the same reduction); if a position is anomalous it falls back to collapsing the
8-state code.

## 3. pLDDT — read from the B-factor column

AlphaFold writes the per-residue confidence (**pLDDT**, 0–100) into the **B-factor** column of every
atom. mdtraj does not carry B-factors through its topology, so pLDDT is read separately from the
file's CA atoms and matched to mdtraj residues by `(chain id, author residue number)`:

- **PDB**: columns 61–66 of the `CA` `ATOM` line (first altloc only).
- **mmCIF**: `_atom_site.B_iso_or_equiv` of the `CA` atom, model 1, first altloc — parsed
  header-driven because column order varies between sources.

For an AlphaFold model (single chain `A`, author numbering `1..N`, label == author), this match is
exact. On an **experimental** structure the same column is a crystallographic B-factor, not pLDDT —
so the band/mask columns are meaningless there (see SKILL.md gotchas).

### Confidence bands (official AlphaFold scale)

| Band | pLDDT | Report colour |
|---|---|---|
| `very_low`  | `< 50`     | `#FF7D45` (orange) |
| `low`       | `50 – <70` | `#FFDB13` (yellow) |
| `confident` | `70 – <90` | `#65CBF3` (cyan)   |
| `very_high` | `≥ 90`     | `#0053D6` (blue)   |

These boundaries reproduce the AlphaFold DB's own `fractionPlddt{VeryLow,Low,Confident,VeryHigh}`
values (verified to within rounding on real models), so the reported band distribution matches what
the AlphaFold website shows.

## 4. Masking and disordered segments

- **Mask cutoff** (`--plddt-cutoff`, default **70**): every residue with `pLDDT < cutoff` is replaced
  with `X` in `ss8_masked` / `ss3_masked`. The raw `ss8` / `ss3` are always preserved — masking is
  non-destructive. `--no-mask` disables replacement entirely.
- **Disordered segments**: maximal contiguous runs (within one chain) of residues with
  `pLDDT < cutoff`, reported as `{chain, start, end, length, mean_plddt}` (start/end are author
  residue numbers). These are the "unreliable disordered regions" — the natural thing to trim before
  docking, MD, or downstream structural analysis.

Choosing the cutoff: `70` keeps only the confident + very-high bands (the standard threshold for
trusting AlphaFold local structure). Use `50` to keep the `low` band too and only discard the
genuinely very-low/disordered stretches.

## 5. Edge cases

- **Multi-model files** (NMR-style): coordinates and DSSP use model 1.
- **Multi-chain / multimer**: each chain is handled independently; chain breaks are respected by DSSP
  and disordered segments never span a chain boundary.
- **Non-protein residues** (ligands, ions, water): mdtraj marks them `NA`; they are excluded from the
  sequence and all SS/QC outputs.
- **Unmatched pLDDT**: if a residue has no CA B-factor, its band is `unknown` and it is left unmasked.
