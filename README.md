# Structural Bioinfo Skills

Four Claude Code skills for everyday structural-biology grunt work: fishing the right PDB entries out of RCSB, splitting them into clean monomer+ligand files, mapping immune epitopes back to their parent proteins, and collecting pMHC / TCR-pMHC complexes. Stdlib Python only — no pip install, no API keys.

## Quickstart

In Claude Code:

```
/plugin marketplace add skblnw/structural-bioinfo-skills
/plugin install structural-bioinfo
```

Then just ask Claude in natural language — the right skill triggers on intent. No CLI to memorise.

## Why these exist

Every structural-bio project starts the same way: paste the same RCSB GraphQL query, the same UniProt chain-resolution lookup, the same IEDB peptide-search hack, into yet another notebook. These skills are the versions actually used in my pipeline, cleaned up and given triggers Claude can recognise. Hand it a peptide, a UniProt accession, or a PDB ID — get back a CSV, a clean PDB, or a Markdown+HTML report.

The skills compose: holostructure search → extractor for drug-target work; epitope-host → pMHC/TCR search for immunology work. Each is also fine on its own.

## Skills

- **[pdb-holostructure-search](plugins/structural-bioinfo/skills/pdb-holostructure-search/SKILL.md)** — find PDB entries of a target protein (UniProt / PDB ID / sequence / gene name) with bound small-molecule ligands; water, ions, buffers, and cryoprotectants auto-filtered. Outputs a CSV ready for batch download.
- **[pdb-extractor](plugins/structural-bioinfo/skills/pdb-extractor/SKILL.md)** — download PDB / mmCIF, resolve the right chain via UniProt GraphQL, write a clean monomer + bound-ligand PDB with correct columns. Single ID or CSV batch.
- **[search-epitope-host](plugins/structural-bioinfo/skills/search-epitope-host/SKILL.md)** — map an epitope peptide to its parent UniProt protein via IEDB, verify the position on the canonical sequence, and return the best 3D structure: experimental PDB if any exist, AlphaFold otherwise.
- **[pdb-pmhc-tcr-search](plugins/structural-bioinfo/skills/pdb-pmhc-tcr-search/SKILL.md)** — collect pMHC and TCR-pMHC crystal structures for one or many peptides; auto-classifies MHC class I vs II and splits TCR-pMHC complexes into derivative `<pdb>_pmhc.pdb` and `<pdb>_tcr.pdb` for docking / MD / pharmacophore work.

## Typical pipeline

```bash
# 1. Find every holostructure of ATR
python plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py Q13535 -o atr_holo.csv

# 2. Download and extract clean monomer+ligand PDBs
python plugins/structural-bioinfo/skills/pdb-extractor/scripts/download_and_extract.py --csv atr_holo.csv --uniprot Q13535 --output-dir raw/
```

For an epitope-side workflow, `search-epitope-host` and `pdb-pmhc-tcr-search` are the entry points.

## Requirements

- Python 3.6+ (stdlib only — no pip install)
- Internet access for RCSB PDB, UniProt, IEDB, and the AlphaFold DB

## Updating

This repo is a downstream mirror of my local plugin at `~/.claude/plugins/structural-bioinfo/`. Sync is automated via the [`sync-structural-bioinfo`](https://github.com/skblnw/structural-bioinfo-skills) skill living at `~/.claude/skills/sync-structural-bioinfo/` — ask Claude to "sync structural-bioinfo" and it mirrors, commits, and pushes in one shot.

## License

MIT — see [LICENSE](LICENSE).
