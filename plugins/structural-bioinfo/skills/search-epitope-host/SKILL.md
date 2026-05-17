---
name: search-epitope-host
description: Map an epitope peptide (or a list of them) to its parent host protein, verify the epitope sits at the reported position on the canonical UniProt sequence, and report the best available 3D structure — experimental PDB entries if any exist, otherwise the AlphaFold prediction. Use whenever a user provides peptide sequence(s) and asks where the epitope comes from, what protein it belongs to, where it is located in the parent, or what structure(s) cover it. Outputs everything into a single output directory: Markdown + HTML report (with a top-of-report summary table), a master `epitopes.csv`, and one `structures/<UNIPROT>.csv` per matched parent. Uses the IEDB IQ-API (`epitope_search`), UniProt REST, and the AlphaFold Database API; no external dependencies.
---

# Search epitope host

Given one or more epitope peptide sequences, return — for each — the parent protein(s) curated in IEDB, the epitope's location on the canonical UniProt sequence (with a verification check), and the best available 3D structure: experimental PDB first, AlphaFold prediction as fallback.

## When to use this skill

Invoke this skill whenever the user supplies one or more peptide sequences and asks any of:

- "What protein does this epitope come from?"
- "Find the parent antigen of `SIINFEKL`."
- "Where is `GILGFVFTL` located in flu M1?"
- "Do we have a structure containing this epitope?"
- "Map these MHC ligands back to their host proteins."
- "Give me UniProt IDs for these peptides and their PDB structures (or AlphaFold if no PDB)."

The skill is the right fit for **short linear peptides (8–25 aa)** that look like T- or B-cell epitopes. For longer protein fragments, prefer a BLAST workflow instead.

## How to invoke the bundled script

The script is stdlib-only. Run it from the skill directory:

```bash
python scripts/search_epitope_host.py <PEPTIDE> [PEPTIDE ...] [options]
python scripts/search_epitope_host.py --input epitopes.txt
```

Inputs accepted:

| Form | Example |
|---|---|
| One peptide as a positional arg | `SIINFEKL` |
| Many peptides as positional args | `SIINFEKL GILGFVFTL NLVPMVATV` |
| `--input file.txt` — one peptide per line | a plain list |
| `--input file.fasta` — FASTA | sequences across multiple lines per `>` block are joined |
| `--input file.csv` / `.tsv` — delimited | auto-detects header; picks a `sequence`/`peptide`/`epitope` column, otherwise the column with the most AA-only tokens |

Duplicate peptides are silently collapsed before any network call, so a 35 k-row CSV with many repeats only costs you the unique count in IEDB requests. Empty cells and `*` (stop codon) are stripped during normalization.

### Options

| Flag | Default | Purpose |
|---|---|---|
| `--out-dir`, `-o` | `epitope_report` | Output directory. Created if missing. All artifacts (report, CSVs, JSON) are written inside it. |
| `--md-only` | off | Skip the HTML report. |
| `--html-only` | off | Skip the Markdown report. |
| `--no-alphafold` | off | Do not query AlphaFold even when the parent has no PDB entry. |
| `--no-json` | off | Skip the raw `report.json` dump. |
| `--workers`, `-j N` | `8` | Number of concurrent worker threads on the outer peptide loop. Stdlib `ThreadPoolExecutor`. Drop to `1`–`4` if IEDB/UniProt start returning HTTP 429. For 100-peptide inputs the speedup over `-j 1` is ~7×; scales linearly to a few thousand peptides. |

### Pipeline (per epitope)

1. **IEDB `/epitope_search`** — exact match on `linear_sequence`. Pulls every curated source antigen for the peptide plus the epitope-level fields `qualitative_measures`, `mhc_classes`, `mhc_allele_names`, and `parent_source_antigen_iris` / `_names` / `_source_org_names`.
2. **Resolve UniProt parents.** First pass: every `UNIPROT:` IRI in `curated_source_antigens` (keeps the IEDB-claimed start/end positions). Second pass: any additional accessions from `parent_source_antigen_iris` (IEDB's UniProt-normalized parent — this is what surfaces peptides whose only curated source antigen is in GenPept/NCBI but which still map to a UniProt entry). Versioned accessions (`P01012.2`) are normalized to bare accessions (`P01012`).
2b. **Per-allele assay outcomes.** For every matched epitope, also query `/mhc_search` and `/tcell_search` to get one row per assay with `mhc_allele_name` + `qualitative_measure`. These are merged into `hla_outcomes`: a `{HLA-A*02:01: {Negative, Positive}, …}` dict that reveals when the same peptide is reported Positive for one allele and Negative for another. Two extra GETs per matched epitope.
3. **UniProt JSON fetch** for each unique parent accession. Pulls the canonical sequence, name, organism, and all PDB cross-references with method/resolution/chains.
4. **Position verification.** For each (start, end) IEDB reports, the script slices the canonical sequence and confirms it equals the input peptide (`✓` / `✗`). Independent of the claim, it also scans the full sequence with `find` and reports every occurrence — this catches the occasional case where the IEDB position is wrong but the peptide still exists at the right spot.
5. **Structure decision.** If the UniProt entry has any PDB cross-references, those are reported and **AlphaFold is not queried**. Only if PDB cross-refs are empty does the script call `https://alphafold.ebi.ac.uk/api/prediction/{acc}` and surface `pdbUrl`/`cifUrl`/version.

## Output schema

Everything is written inside the directory passed to `--out-dir` (default `./epitope_report/`):

```
<out-dir>/
├── report.md             # Markdown report — summary table first, then per-epitope detail
├── report.html           # HTML report with clickable links and inline CSS
├── report.json           # Raw result list (skipped with --no-json)
├── epitopes.csv          # Master mapping: one row per (epitope, parent UniProt)
├── structures/
│   └── <UNIPROT_ACC>.csv     # One file per matched parent — PDB or AlphaFold rows
└── parents/
    ├── <UNIPROT_ACC>.fasta   # Full-length canonical UniProt sequence
    └── <UNIPROT_ACC>.html    # Per-parent epitope map — sequence with epitopes highlighted
```

**`report.md`** opens with a **summary section** (one row per (epitope, parent) — same data as `epitopes.csv`) followed by **Per-epitope details**: parent tables and full PDB / AlphaFold tables per parent, with links to each parent's structures CSV, FASTA, and epitope map.

**`report.html`** has the same Summary table only — the Per-epitope details section is omitted because at >100 input peptides the HTML balloons past tens of MB and browsers struggle to render. For per-epitope detail, point users to `report.md` or to the per-parent files under `parents/`. The summary table itself carries all the high-signal columns including `HLA outcomes`.

**`epitopes.csv` columns:**
`epitope, length, status, uniprot_acc, protein_name, organism, assay_outcomes, mhc_classes, mhc_alleles, hla_outcomes, uniprot_length, claimed_positions, position_verified, occurrences, n_pdb, has_alphafold`

Possible `status` values: `matched`, `no_iedb_match`, `uniprot_fetch_failed`, `invalid: <reason>`, `error: <message>`. Unmatched and invalid epitopes still get one row each so every input is accounted for.

**Assay outcomes.** `assay_outcomes` is the `;`-joined deduplicated list of IEDB `qualitative_measures` for that epitope, copied verbatim. Possible raw values: `Positive`, `Positive-High`, `Positive-Intermediate`, `Positive-Low`, `Negative`. A single peptide can carry both `Positive*` **and** `Negative` when different assays disagree — common, not a bug. `mhc_classes` is `I`, `II`, or `I; II`; `mhc_alleles` is capped at 5 (with a `…(+N)` suffix when truncated) to keep cells readable.

**Per-allele outcomes (`hla_outcomes`).** The same outcomes broken down by HLA allele, format `HLA-A*02:01[Negative,Positive]; HLA-B*07:02[Negative]`. Built from `/mhc_search` + `/tcell_search` rows. Use this when investigating *why* `assay_outcomes` is mixed — sometimes a peptide is a strong binder to one allele and a non-binder to another, which is not real label noise but real biology. Capped at 8 alleles per row.

For datasets labeled "negative", the leakage query is **`status == matched` AND `assay_outcomes` contains any `Positive*`** — those are peptides the user labeled negative that IEDB has positive evidence for. `hla_outcomes` tells you whether the positive was on the *same* allele as the user's target restriction (real leakage) or a different allele (could be excusable depending on the experiment).

**`structures/<UNIPROT_ACC>.csv` columns:**
`uniprot_acc, source, entry_id, method, resolution, chains, range, version, url`
`source` is `PDB` or `AlphaFold`. PDB entries fill method/resolution/chains; AlphaFold entries fill range/version. A CSV is written only for parents whose UniProt fetch succeeded.

**`parents/<UNIPROT_ACC>.fasta`** — the canonical UniProt sequence, wrapped at 60 columns, with a header line `>ACC name OS=organism LEN=n`. Written for every parent whose UniProt fetch succeeded and that has at least one verified epitope occurrence in the canonical sequence.

**`parents/<UNIPROT_ACC>.html`** — a single-page per-parent epitope map. Lists every epitope that lands in this parent (with the IEDB position, the verified occurrence(s), and a colored chip) and shows the full sequence in 60-aa lines, 10-aa blocks, with each epitope highlighted in its own background color. When two or more input epitopes map to the same parent (e.g., two HLA-A*02:01 epitopes from SARS-CoV-2 Spike), they appear together on the same page so the user can see their relative locations at a glance. The FASTA sits next to the HTML in the same directory so no dir-hopping is needed.

A `✓` next to a position means the canonical sequence really contains the epitope at that position. `✗` means the IEDB-claimed position no longer matches the canonical sequence; the report still surfaces every occurrence found by direct substring search so the user can spot drift.

## Examples

```bash
# Classic positive control — ovalbumin SIINFEKL
python scripts/search_epitope_host.py SIINFEKL -o sieg/
#   -> sieg/report.md, sieg/report.html, sieg/epitopes.csv,
#      sieg/structures/P01012.csv (12 experimental PDBs),
#      sieg/parents/P01012.fasta, sieg/parents/P01012.html

# Influenza M1 epitope, ranks the canonical strain P03485 plus near-identical variants
python scripts/search_epitope_host.py GILGFVFTL -o m1/

# Batch from file, write only HTML (still writes CSVs and JSON)
python scripts/search_epitope_host.py --input my_peptides.txt --html-only -o batch/

# Force AlphaFold-free output (only experimental PDB)
python scripts/search_epitope_host.py SIINFEKL --no-alphafold -o pdb_only/
```

## Workflow guidance for the assistant

1. Always run the script — don't hand-craft IEDB / UniProt / AlphaFold URLs yourself.
2. Pass `-o <dir>/` so the output bundle lives in its own folder (the directory will be created if missing). Point the user there when reporting back — the top-of-report summary table and `epitopes.csv` are the main entrypoints; per-parent CSVs in `structures/` are for downstream scripting.
3. After the run, **scan `epitopes.csv` (or the report's Summary table) for non-`matched` statuses**: `no_iedb_match`, `uniprot_fetch_failed`, `invalid: …`, `error: …`. Mention any to the user. `uniprot_fetch_failed` typically means a stale IRI; the other parents for the same epitope are still resolved correctly.
4. **Scan for `✗` verification marks** in the position column: a stale IEDB position. The "found at …" suffix gives the corrected location.
5. If a single epitope returns **multiple parents from different organisms** (common for conserved viral motifs), summarize the list and ask the user which they care about before any downstream work.
6. If a parent has **no PDB and no AlphaFold** (small viral peptides, fragments, partial sequences in TrEMBL), flag it — the user may want to switch to the full-length reference proteome accession instead.
7. Don't dump the full Markdown into chat; point the user to the directory and quote one or two highlights.
8. **For large inputs (>500 peptides)** raise `--workers` (default 8 is safe; 12–16 also fine on a quiet network). Most of the wall-clock is the per-peptide IEDB round-trip, which parallelizes well — concurrency is far more effective than truncating the input. If you see HTTP 429 / "Connection reset" in stderr, drop to `--workers 4`. UniProt/AlphaFold results are cached per process, so a follow-up rerun is cheaper.
9. **Auditing a "negative" or "non-binder" dataset for IEDB leakage:** filter `epitopes.csv` to `status == matched` then look at `assay_outcomes`. Pure-`Negative` rows mean IEDB and the user agree the peptide is a non-binder. Rows containing any `Positive*` value mean IEDB has positive immunogenicity evidence — surface those to the user as labeling disagreements (likely candidates for removal from a negative training set).

## Gotchas

- IEDB's `epitope_search` is exact-match on `linear_sequence`. Off-by-one peptides won't be found — confirm the user's exact sequence before assuming "not in IEDB".
- IRIs can carry a version suffix (`P01012.2`). The script strips it; the UniProt accession resolved is always the bare form.
- Non-UniProt **curated** source antigens (`GENPEPT:`, `NCBI:`, etc.) still resolve to a UniProt parent via IEDB's `parent_source_antigen_iris` field. Only epitope rows with no UniProt parent at all are skipped. As a result the hit-rate of this script is materially higher than older versions that filtered strictly on `curated_source_antigens[*].iri`.
- TrEMBL entries (Q*-prefixed) often lack PDB and AlphaFold predictions; that's not a bug. SwissProt entries are usually richer.
- AlphaFold's API returns an empty list (not 404) for accessions it doesn't model — the script treats that as "no AlphaFold model".

## Extended API patterns

For one-off needs beyond the bundled pipeline — substring/`ilike` matching in IEDB, peptide-to-MHC linkage queries, fetching CIF files in bulk, BFVD or other predicted-structure databases — see `references/api_examples.md`.
