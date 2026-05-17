---
name: pdb-pmhc-tcr-search
description: Search the RCSB Protein Data Bank for pMHC and TCR-pMHC complex structures by epitope peptide sequence. Use whenever the user asks to "find PDB structures of an epitope", "search PDB for a pMHC structure", "get TCR-pMHC complexes for a peptide", "collect HLA-bound peptide crystals", "find pMHC structures of SIINFEKL / GILGFVFTL / KRAS G12D / any 8–25mer peptide", "build a structural dataset for an MHC-I or MHC-II epitope", or supplies one or several short peptide sequences and wants the corresponding PDB entries. Accepts a single peptide, a comma-separated list, or a file of peptides (one per line, optional label column). Covers both MHC Class I and Class II — class is auto-classified per entry. Distinguishes pMHC-only entries from TCR-pMHC complexes; for TCR-pMHC entries it splits the structure into derivative `<pdb>_pmhc.pdb` (peptide + MHC heavy + β2m / class-II β) and `<pdb>_tcr.pdb` (TCR α + TCR β) files for downstream docking / MD / pharmacophore work. Redundant crystal copies are KEPT (not deduplicated) and grouped in the report — the lowest-resolution member of each group is flagged as `representative`. Produces three artifacts in the output directory — `report.md` (GitHub-flavored markdown), `report.html` (self-contained, sortable), and `report_data.json` (raw records for re-rendering). Stdlib Python only — no external dependencies.
---

# pdb-pmhc-tcr-search

Search the RCSB Search API for **pMHC** and **TCR-pMHC** complex structures by epitope
sequence, download them, split TCR-pMHC complexes into derivative pMHC-only and TCR-only
PDB files, and emit a markdown + HTML report.

## When to use this skill

- The user provides one or more **short peptide sequences** (epitopes, 7–30 aa) and wants
  the corresponding PDB structures.
- The user is doing **TCR engineering**, **immunopeptidomics structural analysis**,
  **neoantigen modelling**, **HLA-bound peptide visualisation**, or assembling a
  structural training set for an immunology ML model.
- The user explicitly asks for "pMHC structures", "TCR-pMHC complexes", "HLA structures
  with peptide X", "MHC-bound peptide PDBs", or similar.

This skill is **complementary** to two siblings in the `structural-bioinfo` plugin:

- `pdb-holostructure-search` — protein-ligand co-crystals by *protein* identity.
- `pdb-extractor` — single-monomer chain extraction with bound small-molecule ligands.

Use this skill instead of either when the input is an **epitope peptide** and the output
must distinguish pMHC vs. TCR-pMHC and produce derivative split structures.

## How to invoke

Main entrypoint:

```bash
SKILL=$HOME/.claude/plugins/structural-bioinfo/skills/pdb-pmhc-tcr-search
python $SKILL/scripts/search_pmhc_tcr.py --epitope SIINFEKL --output-dir collect/
```

Three equivalent input modes (combine them freely):

```bash
# Single peptide
python $SKILL/scripts/search_pmhc_tcr.py --epitope SIINFEKL

# Comma-separated list
python $SKILL/scripts/search_pmhc_tcr.py --epitopes SIINFEKL,GILGFVFTL,GLCTLVAML

# File: one peptide per line; optional whitespace-separated label column
python $SKILL/scripts/search_pmhc_tcr.py --epitope-file peptides.txt
```

Example `peptides.txt`:

```
SIINFEKL    OVA257-264  (mouse H-2Kb)
GILGFVFTL   flu_M1
GLCTLVAML   EBV_BMLF1
```

### Options

| Flag | Purpose |
|---|---|
| `--output-dir DIR` | Root output directory (default `collect/`). |
| `--no-download` | Skip downloading PDBs (still emits reports). |
| `--no-split` | Skip writing the `_pmhc.pdb` / `_tcr.pdb` derivatives. |
| `--query-only` | Print the RCSB Search API JSON queries and exit. |

### Per-PDB output layout

```
<output-dir>/
├── report.md
├── report.html
├── report_data.json
└── <epitope>/
    └── <PDB_ID>/
        ├── <PDB_ID>.pdb         # cached full structure (CIF→PDB fallback if needed)
        ├── <PDB_ID>_pmhc.pdb    # peptide + MHC heavy + β2m / class-II β chains
        └── <PDB_ID>_tcr.pdb     # TCR α + TCR β chains (TCR-pMHC entries only)
```

For pMHC-only entries (no TCR present), `<PDB_ID>_pmhc.pdb` is still written but no
`_tcr.pdb` is produced.

## How the search works

For each epitope, a single **RCSB Search API `seqmotif` query** is issued. `seqmotif` is
the only RCSB search service that handles short peptides — `sequence` requires ≥25 aa
input (BLAST minimum) and the `text` service does not index the canonical sequence
attribute. See `references/rcsb_query_examples.md` for the verbatim query and the reasons
for this choice.

No server-side MHC filter is applied (the relevant UniProt attribute also isn't text-
searchable in RCSB). Instead, every seqmotif hit is fetched via the Data API GraphQL
endpoint, and the **classifier** drops any entry that doesn't have at least one MHC
chain. For 8+ aa epitopes the over-fetch is small — most pMHC-relevant epitopes hit
fewer than 100 entries before classification.

Each polymer entity is tagged with a structural role (`mhc_heavy`, `mhc_light`,
`mhc_ii_alpha`, `mhc_ii_beta`, `tcr_alpha`, `tcr_beta`) by UniProt match first,
description regex second (see `references/tcr_keywords.md`). **Peptide-hosting** is
checked *independently*: an entity is flagged as carrying the epitope if its canonical
sequence equals the epitope (free peptide) OR contains it AND the entity has a
recognised MHC/TCR role (single-chain trimer or fusion construct — substring matches on
unrelated proteins like full ovalbumin are rejected). The same chain can therefore be
both the MHC heavy chain AND the peptide host.

An entry is dropped (and listed under "Caveats" in the report) if it has no MHC chain or
no peptide-hosting chain after this pass.

Entries are then grouped by `(MHC UniProts, TCR UniProts)`; the lowest-resolution member
of each group (tie-broken by deposit date) is marked as `representative`. **All entries
are kept** — no entries are silently dropped on redundancy grounds, per design choice.

## Report contents

`report.md` and `report.html` both contain, per epitope:

- Summary line: `N entries · M pMHC · K TCR-pMHC · R redundancy groups`
- A table: `group | rep | PDB | type | class | resolution Å | method | match | MHC | TCR | files`
  - `match` is `free` / `fused` / `both` (which RCSB query strategy hit the entry)
  - `MHC` shows allele names + UniProt(s)
  - `files` links to the full PDB, the pMHC derivative, and the TCR derivative

The HTML report adds an inline column-sorter (no external JS), and links each PDB ID to
both the RCSB structure page and the 3D viewer.

## Examples

```bash
# Mouse Class I — OVA SIINFEKL
python $SKILL/scripts/search_pmhc_tcr.py --epitope SIINFEKL --output-dir ova_collect

# Human Class I — flu M1
python $SKILL/scripts/search_pmhc_tcr.py --epitope GILGFVFTL --output-dir m1_collect

# Multi-epitope cancer neoantigen sweep
python $SKILL/scripts/search_pmhc_tcr.py \
  --epitopes KLVALGINAV,ASNENMETM,FLYALALLL \
  --output-dir neoag_collect

# Class II — HIV gag p24 epitope
python $SKILL/scripts/search_pmhc_tcr.py --epitope PEVIPMFSALSEGATP --output-dir gag_p24

# Inspect the query without hitting RCSB
python $SKILL/scripts/search_pmhc_tcr.py --epitope SIINFEKL --query-only
```

## References

- `references/mhc_uniprot.md` — curated MHC + β2m UniProt accessions used in the filter.
- `references/tcr_keywords.md` — UniProt + description regex heuristics for TCR α/β
  classification.
- `references/rcsb_query_examples.md` — verbatim RCSB Search and Data API payloads.
- `references/report_template.html` — the HTML template (placeholders + embedded sort JS).

## Dependencies

Standard library only (`urllib.request`, `argparse`, `json`, `re`, `pathlib`, `html`,
`datetime`). No `requests`, no `biopython`, no `pandas`.

## Caveats

- **γδ TCRs** are not classified as TCR-pMHC. They will land in the "dropped" bucket if
  no αβ TCR signature is present.
- **Engineered chains** (chimeric MHC, soluble TCR with non-native constants) may evade
  both the UniProt and description regex. The report's "Caveats" section lists every
  dropped entry with its chain descriptions, making such cases easy to spot.
- The search is **case-sensitive** at the API level for the epitope sequence (RCSB
  one-letter codes are uppercase). Input is uppercased automatically.
- Splitting a **single-chain trimer** produces a `_pmhc.pdb` that still contains the
  fusion linker as part of the MHC chain. The peptide is not surgically excised from the
  fused chain — only entire chains are partitioned. If a clean peptide-only chain is
  needed, prefer free-peptide hits (`match = free`).
