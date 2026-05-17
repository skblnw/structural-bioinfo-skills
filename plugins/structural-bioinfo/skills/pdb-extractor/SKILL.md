---
name: pdb-extractor
description: This skill should be used when the user asks to "download PDB", "download and extract PDB", "extract clean PDB", "get PDB structure", "extract monomer with ligands", "extract chain from PDB", "get protein structure from RCSB", "process PDB holostructures", "download protein structure", provides a PDB ID or a CSV file of PDB entries, mentions extracting a specific chain and bound ligands from a PDB/mmCIF file, or wants to batch-download and clean PDB structures for a target protein.
---

# PDB Extractor

## Purpose

Download protein structures from the RCSB Protein Data Bank, resolve target
protein chains via UniProt accession (RCSB GraphQL API), and extract clean
PDB files containing only the target monomer chain and bound small-molecule
ligands. Supports single-PDB and CSV batch workflows, with automatic
mmCIF-to-PDB conversion when legacy PDB format is unavailable.

## When to Use

This skill triggers when the user wants to:

- "Download and extract PDB 3ODU for CXCR4"
- "Get the PDB structure for 8U4N with bound ligands"
- "Extract chain A and ligands from 1AKE"
- "Process this CSV of PDB holostructures into clean monomer+ligand PDBs"
- "Download PDB structures for my target protein"
- "Extract clean monomer from a PDB entry"
- "Batch download PDBs from a holostructures table"

The script handles mmCIF-to-PDB conversion automatically. This is essential
for recent cryo-EM structures that are only available in mmCIF format.

## Quick Start

> **Script:** `$HOME/.claude/plugins/structural-bioinfo/skills/pdb-extractor/scripts/download_and_extract.py`

### Single PDB ID with auto-ligand detection (simplest)

```bash
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-extractor/scripts/download_and_extract.py --pdb-id 8U4N --uniprot <UNIPROT_ACCESSION>
```

### Single PDB with explicit ligands

```bash
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-extractor/scripts/download_and_extract.py --pdb-id 3ODU --uniprot <UNIPROT_ACCESSION> --ligands STI,CLR
```

### Single PDB with manual chain (no UniProt)

```bash
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-extractor/scripts/download_and_extract.py --pdb-id 3ODU --chain A --ligands STI
```

### CSV batch processing

```bash
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-extractor/scripts/download_and_extract.py --csv holostructures.csv --uniprot <UNIPROT_ACCESSION>
```

## Usage Modes

### 1. Single PDB ID Mode (`--pdb-id`)

Downloads one structure, resolves chains via RCSB GraphQL, and writes
`<output-dir>/<pdb_id>/<pdb_id>_clean.pdb`. Requires either `--uniprot`
(for automatic chain resolution) or `--chain` (manual override).

### 2. CSV Batch Mode (`--csv`)

Processes all rows in a CSV file. The CSV must contain at minimum a `pdb_id`
column. Optional columns:

| Column | Description |
|--------|-------------|
| `pdb_id` | 4-character PDB identifier (required) |
| `ligand_comp_ids` | Semicolon-separated residue codes (3–5 chars accepted; matched by first 3 chars) |
| `uniprot_ids` | Semicolon-separated UniProt accessions for multi-protein entries |

When `ligand_comp_ids` is missing or empty for a row, auto-ligand detection
is used for that entry.

### 3. Chain Override (`--chain`)

Skip the RCSB GraphQL lookup entirely. Use when:
- RCSB GraphQL is unavailable
- The target chain is known in advance
- The protein maps to multiple chains and the first match is wrong

### 4. Auto-Ligand Detection

When `--ligands` is omitted (single mode) or when the CSV lacks a
`ligand_comp_ids` column, the script enters **auto-ligand mode**: all
HETATM residues are kept **except** common non-ligand artifacts (water,
ions, buffers, cryoprotectants, membrane-mimetic lipids, and unknown
residues).

For the current exclusion list, see the `EXCLUDE_RESIDUES` set at the top
of `scripts/download_and_extract.py`.

## Full CLI Reference

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--pdb-id` | * | — | Single PDB ID (4-character code) |
| `--csv` | * | — | CSV file path for batch processing |
| `--uniprot` | no | — | UniProt accession (e.g. Q13535) |
| `--chain` | no | — | Manual chain ID, skip GraphQL |
| `--ligands` | no | — | Comma/semicolon-separated codes (3–5 chars; matched by 3-char prefix) |
| `--output-dir` | no | `raw/` | Output root directory |
| `--csv-uniprot-col` | no | `uniprot_ids` | CSV column for UniProt IDs |

\* Either `--pdb-id` or `--csv` is required.

## Output Structure

```
<output-dir>/
├── <pdb_id>/
│   ├── <pdb_id>.pdb          # Full downloaded structure (cached)
│   ├── <pdb_id>.cif          # CIF if PDB unavailable (cached)
│   └── <pdb_id>_clean.pdb    # Extracted monomer + ligands
├── <pdb_id_2>/
│   └── ...
```

Cached PDB/CIF files are reused on subsequent runs. Delete the cache file
to force a fresh download.

## How It Works (Internal Flow)

1. **Chain resolution**: Queries `data.rcsb.org/graphql` for polymer entities
   matching the target UniProt accession. Extracts `pdbx_strand_id` values.
   Uses the first matching chain unless multiple are found.

2. **Download**: Fetches from `files.rcsb.org/download/<pdb_id>.pdb`. On
   failure (common for recent EM structures), downloads `<pdb_id>.cif` and
   converts to PDB format using author-assigned (`auth_*`) chain, residue,
   and atom identifiers for maximum compatibility.

3. **Extraction**: Parses the PDB file line-by-line:
   - Keeps all `ATOM` records matching the target chain (column 22)
   - Keeps `HETATM` records for requested/auto-detected ligands **only on
     the target chain** (multimer copies on other chains are discarded)
   - Ligand matching uses the first 3 characters of the residue code
     (e.g., `A1EIE` matches `A1E` in the PDB columns 18–20)
   - Writes a trailing `END` record

4. **CIF-to-PDB conversion**: Parses the `_atom_site` loop from mmCIF,
   maps each column by name, and formats lines to the 80-column PDB
   specification with correct column alignments. Non-standard residue
   names from `auth_comp_id` (up to 5 chars) are truncated to the
   PDB-standard 3 characters to prevent column overrun.

## Troubleshooting

**"No chains matching <uniprot> found"**
The RCSB entry may not list the target UniProt, or the protein is not
present. Check the entry on rcsb.org. Use `--chain` to manually specify
the chain if known.

**GraphQL errors or timeouts**
Use `--chain` to bypass GraphQL. Look up chain IDs on the RCSB website.

**Wrong ligand kept or missing**
Use `--ligands` to explicitly specify which 3-letter codes to retain.
For auto-ligand mode, check whether the residue is in `EXCLUDE_RESIDUES`.

**Multi-chain structures (e.g., GPCR/G-protein complexes)**
The script uses the first matching chain. If the complex has multiple
copies (dimer/tetramer), use `--chain` to pick the desired copy. Ligands
are automatically restricted to the selected target chain — copies bound
to other chains are excluded.

**Ligand appears to have wrong chain ID in PDB viewer**
Non-standard residue codes (e.g., 5-char `A1EIE`) are truncated to 3 chars
in the PDB output to maintain column alignment. The chain ID at column 22
is always the target chain. If a viewer shows a different chain, check
that it parses columns 18–22 correctly.

## Additional Resources

### Scripts

- **`scripts/download_and_extract.py`** — Main extraction script. Execute
  directly with Python. Contains all logic: RCSB download, CIF-to-PDB
  conversion, UniProt chain resolution via GraphQL, monomer extraction,
  and ligand filtering. Run with `--help` for full options.
