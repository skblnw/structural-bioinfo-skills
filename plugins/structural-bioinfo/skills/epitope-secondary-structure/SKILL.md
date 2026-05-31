---
name: epitope-secondary-structure
description: Given the output of `search-epitope-host` (epitopes mapped to parent UniProt proteins + their PDB/AlphaFold structures), compute the SECONDARY STRUCTURE each epitope adopts in its parent — a per-position probability over the 8 DSSP states (H,G,I,E,B,T,S,C) aggregated across the parent's experimentally solved structures, plus a simplified 3-state (H/E/C) consensus. Use whenever a user asks "what secondary structure do these epitopes sit in", "is this epitope helical / a strand / a loop in the parent", "map epitopes to DSSP", "per-position helix/sheet/coil probabilities for my peptides", or wants to add structural context (fold/SS) to an epitope→host mapping. Uses experimentally resolved structures (parents with only an AlphaFold model, or none, are skipped and flagged); by default it also excludes isolated peptide chains in pMHC / TCR / antibody complexes so the result reflects the secondary structure the epitope adopts in its PARENT fold, not its MHC-bound conformation. 8-state DSSP is computed locally with mkdssp (the RCSB Data API serves only coarse 3-state SS); the RCSB Data API is used for the UniProt→structure residue mapping (SIFTS). Builds a self-contained HTML report plus CSV/JSON. Stdlib-only except the external `mkdssp` binary.
---

# Epitope secondary structure (DSSP)

Take a set of epitopes already mapped to their parent proteins and structures (the
`search-epitope-host` output) and answer: **for each epitope, what secondary structure does it
adopt in the parent?** The result is, per epitope position, a **probability over the 8 DSSP
states** (`H` α-helix, `G` 3₁₀-helix, `I` π-helix, `E` β-strand, `B` β-bridge, `T` turn,
`S` bend, `C` coil) pooled across the parent's experimentally solved structures, and a simplified
**3-state (H/E/C)** consensus derived from it.

## When to use this skill

Invoke after `search-epitope-host` (or any pipeline that produced an `epitopes.csv` +
`structures/<ACC>.csv` + `parents/<ACC>.fasta` bundle), whenever the user wants structural context:

- "What secondary structure are these epitopes in?"
- "Is `GILGFVFTL` helical or a strand in its parent?"
- "Give me per-position helix/sheet/coil probabilities for my HLA-A\*02:01 9-mers."
- "Map each epitope onto DSSP across all the parent's PDB structures."
- "Which epitopes sit in loops vs. ordered secondary structure?"

The unit of analysis is an **(epitope, parent) row** — the same grain as `search-epitope-host`'s
`epitopes.csv`. An epitope mapped to two parents gets two rows.

## Why DSSP is computed locally (not pulled from RCSB)

The RCSB Data API exposes only **coarse** secondary structure (`HELIX_P` / `SHEET` /
`UNASSIGNED_SEC_STRUCT`, provenance PROMOTIF) — not the 8-state DSSP this skill reports. So 8-state
DSSP is computed locally with `mkdssp`. The RCSB Data API **is** used, but only for the
**UniProt→structure residue mapping**: `rcsb_polymer_entity_align` (provenance SIFTS) maps each
UniProt residue to the structure's `label_seq_id`, which is then reconciled to the author numbering
that DSSP reports via the mmCIF `_atom_site` table.

## Structure selection (experimental-first; parent-fold; AlphaFold skipped)

Per parent:

1. **Experimental PDB present** → use **all** of the parent's PDB structures (after a cheap
   UniProt-range prescreen drops structures that don't cover the epitope window).
2. **AlphaFold-only** → the parent is **skipped** and flagged `no_experimental_structure`
   (opt back in with `--use-alphafold`; note that requires DSSP 4 — the common conda `mkdssp`
   3.1.2 cannot parse AlphaFold mmCIF).
3. **No structure at all** → flagged `no_structure`.

**By default, isolated peptide chains in pMHC / TCR / antibody complexes are excluded**
(`--min-align-length 25`): a structure where the only thing aligned to the parent is the epitope
peptide itself reflects the **MHC-bound** conformation, not the parent fold, so it is dropped. An
epitope whose only covering structures are such complexes is flagged `only_pmhc`. Set
`--min-align-length 0` to include them (e.g. to study the bound conformation).

## How to invoke the bundled script

Stdlib-only Python; needs the external `mkdssp` binary.

```bash
# install the DSSP engine once (provides mkdssp):
conda install -c conda-forge dssp

S=~/.claude/plugins/structural-bioinfo/skills/epitope-secondary-structure/scripts/epitope_ss.py
python "$S" path/to/host_mapping/ \
    --mkdssp-path "$(conda run -n <env> which mkdssp)" \
    -o path/to/host_mapping/secondary_structure/
```

The positional argument is a `search-epitope-host` output directory; the script auto-discovers
`epitopes.csv`, `structures/`, and `parents/` inside it.

### Options

| Flag | Default | Purpose |
|---|---|---|
| `input_dir` | (required) | A `search-epitope-host` output dir (`host_mapping/`). |
| `--epitopes-csv` / `--structures-dir` / `--parents-dir` | (auto under `input_dir`) | Override individual inputs. |
| `--out-dir`, `-o` | `{input}/secondary_structure` | Output directory. |
| `--cache-dir` | `{out}/cache` | Persistent cache of downloaded mmCIF, DSSP, and SIFTS-align (keyed by entry id). |
| `--mkdssp-path` | `mkdssp` | Path to the `mkdssp` binary (e.g. a conda env's `bin/mkdssp`). |
| `--aggregation` | `per-structure` | `per-structure` (each deposition = one vote) or `per-chain` (each modeled chain = one vote). |
| `--use-alphafold` | off | Also annotate AlphaFold-only parents (requires DSSP 4). |
| `--max-structures-per-parent` | `0` (all) | Cap experimental structures per parent (sorted X-ray-then-resolution). |
| `--min-align-length` | `25` | Skip experimental entities whose SIFTS alignment to the parent spans fewer residues — i.e. isolated peptides in pMHC / TCR / antibody complexes — so only genuine parent-fold structures contribute (epitopes left with only such structures become `only_pmhc`). This is **on by default** to report the SS an epitope adopts *in its parent*; set `0` to include the MHC-bound conformation. |
| `--workers`, `-j` | `8` | Worker threads over epitope rows (structure download+DSSP is cached per entry). |
| `--offline` | off | Use only cached cif/dssp/align; never hit the network. |

### Pipeline (per epitope×parent row)

1. **Position** — take the verified `occurrences` window on the parent's canonical UniProt
   sequence (rows with no occurrence → `no_position`).
2. **Select structures** — experimental PDB if any (else skip/flag); prescreen by UniProt-range overlap.
3. **Per structure (cached by entry id):** download mmCIF (`files.rcsb.org`), run `mkdssp`, fetch
   SIFTS alignment (RCSB `rcsb_polymer_entity_align`), parse `_atom_site` (label↔auth + residue
   identity + B-factor).
4. **Map** each UniProt residue → `label_seq_id` (SIFTS) → author residue → DSSP 8-state, for every
   modeled chain of the matching entity. The residue identity is cross-checked against the epitope;
   mismatches are flagged (guards against repeat regions / engineered constructs).
5. **Aggregate** to a per-position 8-state probability (see below), collapse to 3-state, and emit
   consensus strings + reports.

## Aggregation

For each epitope position, observations are pooled across the parent's structures.

- **`per-structure`** (default): each structure is first reduced to its own per-position
  distribution (fraction of that structure's modeled chains in each state), then structures are
  averaged with equal weight. This prevents a single large homo-oligomer (e.g. a 12-chain crystal)
  from dominating — it answers "how many independent depositions agree."
- **`per-chain`**: every modeled chain across all structures is one equal observation.

The **3-state probability** is the 8→3 collapse `{H,G,I}→H, {E,B}→E, {T,S,C}→C`. The 3-state
consensus letter at each position is the argmax of the 3-state probability (not the collapse of the
8-state argmax — these can differ; the 3-state argmax is the more faithful "simplified" call).
Positions modeled in no structure get `.` (8-state) / `-` (3-state) and `n_structures = 0`.

## Output schema

```
<out>/
├── report.html                  # self-contained: summary + per-epitope stacked-probability bars + legend
├── epitope_ss.csv               # one row per (epitope, parent): consensus strings, fractions, status
├── epitope_ss_positions.csv     # one row per (epitope, parent, position): 8 + 3 state probs, n, argmax
├── structures_used.csv          # one row per (epitope, parent, structure): method/resolution/dssp_status
├── epitope_ss.json              # full nested dump (metadata + per-epitope + per-position + per-structure)
└── cache/                       # cif/, dssp/, align/ — reused on re-run (and by --offline)
```

**`epitope_ss.csv` columns:** `epitope, length, uniprot_acc, protein_name, organism, occurrence,
status, source_type, n_structures_used, positions_covered, mean_plddt, ss3_consensus, ss8_consensus,
frac_H, frac_E, frac_C, n_positions, n_positions_unmodeled, n_mismatch_structures, flags`.

**`status` values:** `ok` (SS computed) · `no_position` (epitope not located on the canonical
sequence) · `no_experimental_structure` (AlphaFold-only parent, skipped) · `no_structure` ·
`not_in_any_structure` (no deposited structure covers the window) · `unmodeled` (window unresolved
in all covering structures) · `only_pmhc` (only covering structures are pMHC complexes; excluded by
`--min-align-length`) · `error: …`.

**`source_type`:** `experimental` / `no_experimental_structure` / `no_structure` (or `alphafold`
with `--use-alphafold`).

## HTML report

Light-theme, self-contained, CSS-only (no JS/libraries). A summary section (counts by status and
source) and a per-epitope table with the colored 3-state consensus, then per-epitope detail: a
**stacked vertical probability bar per residue** (8 states colour-coded), the 3-state consensus
strip, the number of structures contributing per position, and a table of contributing structures
(PDB id + method + resolution, DSSP status, mismatch flags). A legend maps every colour.

## Workflow guidance for the assistant

1. Run `search-epitope-host` first; point this skill at that output directory.
2. Always pass `--mkdssp-path` unless `mkdssp` is on PATH.
3. After the run, scan `epitope_ss.csv` for non-`ok` statuses and report them: how many epitopes
   got SS, how many parents were skipped for lacking an experimental structure (`no_experimental_structure`),
   and any `no_position` / `unmodeled` / sequence-mismatch flags.
4. Lead with the headline (e.g. "N/M epitopes resolved; dominant fold call per epitope") and a
   couple of standout cases; point the user at `report.html`. Don't paste the HTML into chat.

## Gotchas

- **DSSP engine:** `conda install -c conda-forge dssp` typically gives **mkdssp 3.1.2**, which
  produces standard 8-state output and reads RCSB mmCIF fine, but **cannot parse AlphaFold mmCIF**
  ("empty protein…"). AlphaFold support (`--use-alphafold`) therefore needs **DSSP 4**.
- A parent with experimental PDB is used **exclusively** — AlphaFold is never mixed in.
- Probabilities are empirical frequencies over deposited structures, **not** a predictor; a
  position seen in one structure has a one-hot distribution.
- The same epitope can map to two parents (e.g. a B95-8 and a GD1 strain entry); each is a separate
  row, and one may be `ok` while the other is `no_structure`.
- A `no_position` row means the epitope sequence wasn't found verbatim on the parent's canonical
  UniProt sequence (often an EBV strain variant) — there is nothing to map onto a structure.

## Extended API patterns

For the RCSB SIFTS-alignment recipe, the mmCIF `_atom_site` label↔auth reconciliation, the `mkdssp`
invocation/column layout, and the 8→3 collapse, see `references/api_examples.md`.

## Dependencies

Python 3.8+ (stdlib only) and the external `mkdssp` binary (`conda install -c conda-forge dssp`).
No other Python packages. Network access to `data.rcsb.org` and `files.rcsb.org` (cached; `--offline`
reuses the cache).
