# structural-bioinfo

Structural biology toolkit for Claude Code — search the PDB for protein-ligand holostructures, extract clean monomer+ligand PDB files, and map immune epitopes to their parent proteins and structures.

## Skills

### `pdb-holostructure-search`

Search the RCSB Protein Data Bank for structures of a target protein (or its homologs) bound to small-molecule ligands. Outputs a CSV table ready for batch extraction.

**Triggers:** "Find PDB structures of EGFR with bound inhibitors", "Get holo structures of P00533", "Co-crystal structures of kinase with drugs"

**Output CSV columns:** `pdb_id`, `resolution_A`, `method`, `title`, `uniprot_ids`, `ligand_comp_ids`, `ligand_names`, `ligand_mw_Da`

### `pdb-extractor`

Download PDB/CIF structures from RCSB, resolve protein chains via UniProt (GraphQL), and extract clean PDB files containing only the target monomer and its bound ligands. Handles mmCIF-to-PDB conversion automatically. Supports single-PDB and CSV batch modes.

**Triggers:** "Download and extract PDB 3ODU", "Extract chain A and ligands from 1AKE", "Process this CSV into clean PDBs"

**Output:** `<pdb_id>/<pdb_id>_clean.pdb` per entry — protein chain + ligand HETATM, PDB-format correct columns.

### `search-epitope-host`

Map an epitope peptide (or a list of them) to its parent host protein via the IEDB IQ-API, verify the position on the canonical UniProt sequence, and report the best 3D structure: experimental PDB entries first, AlphaFold prediction as fallback. Outputs Markdown + HTML reports.

**Triggers:** "What protein does SIINFEKL come from?", "Find the parent antigen of these peptides", "Where is GILGFVFTL on flu M1?", "Map these epitopes to PDB or AlphaFold structures"

**Output:** `<prefix>.md` and `<prefix>.html` — parent protein table (UniProt, name, organism, position + verification) and structure table (PDB or AlphaFold) per epitope.

### `pdb-pmhc-tcr-search`

Search the RCSB Protein Data Bank for **pMHC** and **TCR-pMHC** complex structures by epitope peptide sequence. Auto-classifies MHC class I vs II; splits TCR-pMHC complexes into derivative `<pdb>_pmhc.pdb` and `<pdb>_tcr.pdb` files for downstream docking / MD / pharmacophore work. Redundant copies are kept and grouped; the lowest-resolution member is flagged `representative`.

**Triggers:** "Find pMHC structures of SIINFEKL", "Get TCR-pMHC complexes for GILGFVFTL", "Collect HLA-bound peptide crystals", "Build a structural dataset for an MHC-I/II epitope"

**Output:** `report.md`, `report.html`, and `report_data.json` in the output dir, plus downloaded PDBs and TCR/pMHC split files.

## Workflow

The holostructure skills form a pipeline:

1. **Search** — `pdb-holostructure-search` → `holostructures.csv`
2. **Extract** — `pdb-extractor --csv holostructures.csv --uniprot <ACCESSION>` → clean PDBs

```bash
# 1. Find all holostructures for your target
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py Q13535 -o atr_holo.csv

# 2. Download and extract clean monomer+ligand PDBs
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-extractor/scripts/download_and_extract.py --csv atr_holo.csv --uniprot Q13535 --output-dir raw/
```

`search-epitope-host` is a standalone entry point for the epitope side — given peptides, it produces the parent-protein/position/structure report directly.

```bash
python $HOME/.claude/plugins/structural-bioinfo/skills/search-epitope-host/scripts/search_epitope_host.py SIINFEKL GILGFVFTL -o epitopes
```

## Installation

### Claude Code

```bash
cc --plugin-dir ~/.claude/plugins/structural-bioinfo
```

### OpenCode (via OpenPackage)

```bash
# Install opkg once (skip if already installed: `which opkg`)
command -v opkg >/dev/null || npm install -g opkg

# Install skills to your workspace
opkg install ~/.claude/plugins/structural-bioinfo --platforms opencode
```

## Requirements

- Python 3.6+ (stdlib only — no pip install needed)
- Internet access (RCSB PDB APIs, UniProt REST)

## License

MIT
