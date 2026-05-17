---
name: pdb-holostructure-search
description: Search the RCSB Protein Data Bank for holostructures — PDB entries containing a specified protein or its close homologs bound to a meaningful small-molecule ligand. Use this whenever the user asks for PDB entries with bound ligands, drug-design starting structures, ligand-bound co-crystals, holo structures, or co-complex crystal structures of a target protein, even if they don't explicitly say "holostructure". The skill automatically excludes water, ions, buffers, cryoprotectants, and crystallization additives. Accepts UniProt accession, PDB ID, raw protein sequence (FASTA or plain), or protein/gene name as input. Outputs a CSV table with PDB ID, resolution, experimental method, UniProt IDs, ligand chemical-component IDs, names, and molecular weights.
---

# PDB Holostructure Search

Find PDB entries containing a target protein (and its homologs at a chosen sequence-identity cutoff) bound to a biologically meaningful small-molecule ligand. Built on the RCSB Search API and Data API.

## When to use this skill

Invoke this skill whenever the user wants ligand-bound PDB structures of a protein — even if they don't use the word "holostructure". Common phrasings:

- "Find PDB structures of EGFR with bound inhibitors"
- "Co-crystal structures of <kinase> with drugs"
- "Get holo structures of P00533"
- "PDB entries of thymidine kinase homologs that have a ligand"
- "Drug-design starting structures for <target>"
- "Ligand-bound structures of <protein>"

If the user explicitly wants apo (ligand-free) structures, this skill is not the right fit.

## How to invoke the bundled script

The skill ships one self-contained Python script (stdlib only). Run it from this skill directory:

```bash
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py <TARGET> [options]
```

`<TARGET>` is auto-detected as one of:

| Input form | Example | Detection rule |
|---|---|---|
| UniProt accession | `P00533` | Matches the canonical UniProt accession regex |
| PDB ID | `1ATP` | 4 chars, first is a digit. The script fetches the **longest** protein chain in the entry — important because many entries (e.g., kinase + bound peptide inhibitor) have a short peptide as entity 1. The chain selection is logged to stderr so the user can verify. |
| Sequence | `MKTAYIAKQRQI...` or FASTA with `>` header | All amino-acid letters, length ≥ 20 |
| Name / gene symbol | `"EGFR"`, `"ATR"`, `"thymidine kinase"` | Anything else; resolved via UniProt REST. For **short, no-space queries** (≤7 chars — gene symbols like `EGFR`, `ATR`, `TP53`) the lookup tries (1) gene-name match, (2) exact protein-name match, (3) general text search. For **multi-word queries** (`"thymidine kinase"`) it tries protein-name match first. Results are logged to stderr with the gene symbol shown. If the top match's protein name doesn't contain the query string and the gene name doesn't match either, a warning is printed — common names like "thymidine kinase" can match adjacent enzymes (thymidylate kinase) on a loose query. |

You can override detection with `--target-type uniprot|pdb|sequence|name`.

### Defaults (do not change unless the user asks)

- **Homolog scope**: 90% sequence identity (`--identity 0.90`). At 90%, results are essentially the same protein with point mutations / close orthologs. Lower the cutoff (e.g., `0.50` or `0.30`) only if the user asks for distant homologs or fold-family members.
- **Ligand filter**: MW > 150 Da plus a manual exclusion list applied per-ligand at the metadata stage. The exclusion list covers waters, ions, common buffers (HEPES/Tris/MES/PIPES), cryoprotectants (glycerol/PEG/MPD/EDO), crystallization additives, and frequently-non-functional sugars. This filter works on every PDB entry regardless of deposition date.
- **Quality filters**: none. The user adds `--method` and/or `--max-resolution` if they want them.

### Why these defaults

There's a tradeoff between two ligand-filtering strategies:

| Strategy | Coverage | Precision | Default? |
|---|---|---|---|
| MW > 150 Da + manual exclusion list | Works on every entry | Very good — the exclusion list catches the common offenders, MW catches the rest | ✅ default |
| `SUBJECT_OF_INVESTIGATION` (RCSB curated flag) | Sparse on pre-2020 entries (the annotation programme rolled out ~2020 and isn't fully backfilled) | Excellent — manually curated by RCSB | opt-in via `--loi-only` |

In practice the MW + exclusion approach yields **~10× more hits** than LOI-only on real targets while keeping the output clean. The LOI flag is genuinely better-curated, but coverage gaps make it the wrong default for "give me everything I can dock against". Use `--loi-only` only when the user explicitly wants the curated subset.

### Options

| Flag | Default | Purpose |
|---|---|---|
| `--identity FLOAT` | `0.90` | Sequence-identity cutoff for homolog inclusion (0.0–1.0). |
| `--loi-only` | off | Use RCSB's curated `SUBJECT_OF_INVESTIGATION` flag for ligand filtering. More precise but undercounts pre-2020 entries that lack the annotation. Default is MW > 150 Da + exclusion list, which works on every entry. |
| `--method TEXT` | none | Filter by experimental method, e.g., `"X-RAY DIFFRACTION"`, `"ELECTRON MICROSCOPY"`, `"SOLUTION NMR"`. |
| `--max-resolution FLOAT` | none | Maximum resolution in Å. |
| `--output PATH`, `-o` | stdout | CSV output path. |
| `--max-results INT` | all | Cap number of returned entries. |
| `--query-only` | off | Print the JSON search query without executing. Useful for letting the user inspect or tweak the query. |

### Output schema

CSV with columns:

| Column | Meaning |
|---|---|
| `pdb_id` | PDB entry ID |
| `resolution_A` | Resolution in Å (blank for non-diffraction methods) |
| `method` | Experimental method |
| `title` | Entry title |
| `uniprot_ids` | UniProt accessions of polymer entities, semicolon-separated |
| `ligand_comp_ids` | Chemical-component IDs of the meaningful ligands, semicolon-separated |
| `ligand_names` | Human-readable ligand names, semicolon-separated |
| `ligand_mw_Da` | Molecular weights in Da, semicolon-separated |

Entries whose ligands all get filtered out at the metadata stage are dropped from the output, so every row is guaranteed to have at least one meaningful ligand.

## Examples

```bash
# All EGFR holostructures and close mutants (default 90% identity, all methods)
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py P00533 -o egfr_holo.csv

# Distant kinase-family homologs of EGFR at 30%, X-ray ≤ 2.5 Å only
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py P00533 --identity 0.30 \
  --method "X-RAY DIFFRACTION" --max-resolution 2.5 -o kinase_family.csv

# Use a PDB entry as the sequence seed
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py 1ATP --identity 0.50 -o pka_homologs.csv

# Use a raw protein sequence
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py MKTAYIAKQRQISFVKSHFSRQLEERL... -o by_seq.csv

# Resolve a name to UniProt automatically
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py "thymidine kinase" -o tk.csv

# Use RCSB's curated SUBJECT_OF_INVESTIGATION subset (smaller but more precise)
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py P00533 --loi-only -o egfr_loi.csv

# Just inspect the JSON query the skill would send
python $HOME/.claude/plugins/structural-bioinfo/skills/pdb-holostructure-search/scripts/search_holostructures.py P00533 --query-only
```

## Workflow guidance for the assistant

1. Decide what `<TARGET>` to pass. If the user supplied a UniProt ID, PDB ID, or sequence directly, use it verbatim. If they named a protein:
   - **Watch the stderr log for the "UniProt name lookup" block** — it shows the top 5 candidates, their gene symbols, and which one was picked. For short gene-symbol queries (`EGFR`, `ATR`, `TP53`) the gene-name search runs first to avoid substring-match pitfalls (e.g., `"ATR"` matching `"ATR-interacting protein"` instead of ATR kinase). If the picked accession's protein name doesn't match what the user asked for and the gene doesn't match either (the script prints a warning), stop and ask the user which accession they meant.
   - The first chain selected in a PDB entry is also logged; verify it matches the protein the user actually wants (kinase vs. bound peptide, etc.).
2. Pick non-default flags only when the user requests them. Do not silently add resolution or method filters.
3. After running, briefly summarize the result: number of entries, top 1–3 by resolution, and the ligand spread (e.g., "37 entries, 12 unique ligand chemical components, top hit 5UG9 at 1.20 Å with TKI molecule X"). Don't dump the full CSV into chat — point to the file.
4. If the user wants ligand structural data (SDF, SMILES, binding affinity), tell them to follow up with the Data API; this skill stops at the CSV.

### Distant-homolog gotcha

Sequence search at low identity (e.g., 0.30) seeded from a *single* sequence can miss family branches that have diverged in different directions. Example: human cytosolic thymidine kinase (TK1, P04183) and HSV-1 thymidine kinase share only ~15% identity and don't show up in each other's MMseqs2 hit lists, even at 30% cutoff. If the user asks for "distant homologs" or "the whole family", warn them that one seed may miss divergent branches, and offer to run a second pass with a different seed.

## Extended API patterns

For niche needs — exact ligand chemical-component lookup, binding-affinity filters, covalent-vs-noncovalent ligand splits, ligand-quality (RSCC/RSR) thresholds — see `references/api_examples.md`. Read it only when those patterns are needed; the bundled script covers the holostructure base case.
