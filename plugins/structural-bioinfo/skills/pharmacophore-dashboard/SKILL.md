---
name: pharmacophore-dashboard
description: Generate a self-contained interactive HTML dashboard from structure-based pharmacophore analysis outputs (RDKit *_filtered / *_all pharmacophore JSON + cleaned PDB receptors, optional consensus query). Use when a user wants to visualize, report on, summarize, or share pharmacophore analysis results for a drug target — producing an interactive 3D pharmacophore viewer (3Dmol.js: receptor pocket + ligand + colored feature spheres + H-bond/salt-bridge lines), feature-retention charts, sortable per-structure and residue-engagement tables, and a consensus virtual-screening query. Triggers include: pharmacophore report or dashboard, pharmacophore visualization, virtual screening pharmacophore, interactive pharmacophore HTML, hit-discovery pharmacophore summary, "make a report for the pharmacophore analysis". This is the interactive-HTML counterpart to `pharmacophore-report-generator` (which emits a static markdown writeup from the same per-PDB JSONs): pick this skill when the user wants an interactive, visual, or shareable HTML deliverable with a 3D viewer, and `pharmacophore-report-generator` when they want a plain-text/markdown report. It consumes the per-structure outputs of the `pharmacophore-analyzer` skill.
---

# Pharmacophore Dashboard

Turn structure-based pharmacophore analysis outputs into a single self-contained, interactive
HTML report — one per target. The page embeds all data (no server, no build step to open it) and
renders a 3D pharmacophore explorer, feature-retention charts, sortable tables, a receptor
residue-engagement view, and (when present) a consensus virtual-screening query.

## When to Use This Skill

Use when the user wants to **visualize / report on / share pharmacophore analysis results**, e.g.:
- "Make an interactive report/dashboard for the pharmacophore analysis of <target>"
- "Visualize the pharmacophore features and which residues they engage"
- "Turn these `*_filtered.pharmacophore.json` files into an HTML report"
- "Show the consensus pharmacophore query in 3D for virtual screening"

This skill **consumes existing pharmacophore-detection outputs** (it does not detect features
itself). The upstream step is a structure-based pharmacophore detector (e.g. RDKit
`BaseFeatures.fdef` with an interaction filter) that emits, per analyzed ligand copy:

```
<prefix>_filtered.pharmacophore.json   # retained features (have a complementary protein partner)
<prefix>_all.pharmacophore.json        # all detected features (optional, used for retention %)
```

plus the cleaned receptor PDBs they reference (`input_pdb`), and optionally a
`*consensus*.json` query. See `references/INPUT_FORMAT.md` for the exact schema.

## Approach 1: One-shot (Recommended)

Point `generate_report.py` at the directory containing the `*_filtered.pharmacophore.json`
files. It auto-discovers structures/chains/ligands, locates the cleaned PDBs, extracts each
binding pocket, computes statistics, auto-detects a consensus query, and writes the HTML.

```bash
python3 scripts/generate_report.py /path/to/<target>/pharmacophore \
    --raw-dir /path/to/<target> \
    --target-name "C-X-C chemokine receptor type 4" --gene CXCR4 --uniprot P61073 \
    --accent "#0d9488" \
    --out CXCR4_pharmacophore_report.html
```

- `--raw-dir` is the base used to resolve each JSON's `input_pdb` (typically the target folder
  that contains `raw/<PDB>/<PDB>_clean.pdb`). If omitted, the parent of the pharmacophore dir is
  tried, plus a recursive fallback search. Broken symlinks are skipped automatically.
- A consensus query (`*consensus*.json`) in the directory is detected and overlaid on its
  alignment-reference structure. Pass `--consensus FILE` to specify one explicitly.
- With no metadata flags and no config, the report still works: the executive summary,
  limitations, and manifest are **auto-generated** from the data.

**Requirements:** Python 3 with `numpy`. The HTML loads `3Dmol.js` from a CDN at view time —
internet is needed **only for the 3D viewers**; charts, tables, and all data are fully offline.

## Approach 2: Modular (data bundle, then HTML)

Run the two stages separately when you want to inspect or post-process the intermediate data,
batch many targets, or regenerate the HTML without re-extracting:

```bash
# 1) Extract a portable JSON data bundle (pockets, features, stats, consensus)
python3 scripts/extract_pharmacophore_data.py <target>/pharmacophore \
    --raw-dir <target> --out <target>_data.json

# 2) Build the dashboard from the bundle (+ optional narrative config)
python3 scripts/build_dashboard.py <target>_data.json \
    --config my_report_config.json --out <target>_report.html
```

## Customization: report config (optional)

Pass `--config report_config.json` to add narrative/branding and to override how each ligand is
classified. Anything you omit is auto-generated. Supported keys:

- **Branding:** `accent`, `accent_dark`, `accent_soft` (hex), `kicker`, `meta`.
- **Target:** `target: {name, gene, uniprot, organism, n_structures}`.
- **`entry_overrides[]`** — per-structure or per-entry labels and classification, keyed by
  `pdb_id` (whole structure) or `key` (one ligand copy = filtered filename without
  `.pharmacophore.json`). Fields: `label`, `chemotype_label`, `role`, `note`, and `category`.
  Set `category` to `excluded` / `control` / `negative` to treat a ligand as a **negative
  control** (kept in the table, dropped from retention & residue aggregation) — this is how you
  separate orthosteric binders from incidental binders such as cholesterol. Categories
  `orthosteric` / `binding` / `productive` count toward the aggregates.
- **Narrative:** `summary_html`, `modes[]` (query-recommendation cards), `strategies[]`
  (screening strategy cards shown with a consensus query), `interaction_map` +
  `interaction_map_note`, `caveat`, `limitations[]` (`{issue, impact, mitigation}`), `manifest`,
  `footer`. HTML is allowed in the text fields.

Worked examples that reproduce polished, publication-style reports are in `examples/`:

```bash
python3 scripts/generate_report.py atr/pharmacophore --raw-dir atr \
    --config examples/atr_report_config.json   --out ATR_report.html
python3 scripts/generate_report.py cxcr4/pharmacophore --raw-dir cxcr4 \
    --config examples/cxcr4_report_config.json --out CXCR4_report.html
```

## What the report contains

1. **Header** — target name, UniProt link, and stat tiles (structures, entries, productive poses, retained features, consensus size).
2. **Executive summary** — auto-generated from the data, or supplied via config.
3. **Feature retention** — grouped bar chart (all vs retained per family) with INCLUDE/EXCLUDE verdicts. Aromatic/hydrophobic families that never engage the receptor surface here.
4. **Structures & 3D explorer** — sortable per-structure table (detected/retained/%, family composition) wired to a 3Dmol.js viewer: ligand sticks, labeled engaged residues, feature spheres colored by family, dashed H-bond/salt-bridge lines, with family/H-bond/label/tolerance-sphere/spin toggles. A synced feature-detail table lists every retained feature and its partner residue, atom, and distance.
5. **Residue engagement** — stacked bar chart + table of which receptor residues the features engage.
6. **Consensus query** (if a consensus JSON is present) — mandatory/optional feature tables, a 3D overlay on the reference pocket, and screening-strategy cards. Or **query recommendations** cards if `modes[]` is supplied.
7. **Global interaction map**, **data-quality/limitations**, and a **data manifest** (collapsible).

## Outputs

- A single `*.html` file (typically 130–280 KB) — self-contained and shareable.
- Optionally the intermediate `*_data.json` bundle (`--data-out` on `generate_report.py`, or the
  `extract_*` step) for reuse or external tooling.

When sharing with a user, copy the HTML to their workspace folder and present that single file.

## Best Practices

1. **Resolve the PDBs.** If the run reports `no clean PDB for <id>`, set `--raw-dir` to the folder
   holding the cleaned receptors (the 3D viewer for those entries is disabled otherwise).
2. **Classify non-orthosteric ligands.** Cholesterol, lipids, buffer components and other
   incidental binders may retain a stray feature. Mark them `category: excluded` via
   `entry_overrides` so retention and residue stats reflect the true binding chemotypes; they
   remain in the table as built-in negative controls.
3. **Trust the JSON, not prose.** All counts, coordinates, partners and distances are read
   directly from the per-structure filtered JSON files, so the report stays consistent with the
   underlying analysis. Cross-check low-confidence ligands (`smiles_source` not a clean template
   match) in PyMOL.
4. **One report per target.** Run the skill once per target directory.

## Reference Materials

`references/INPUT_FORMAT.md` documents the expected input file layout, the per-feature JSON
schema, the consensus-query schema, the full report-config schema, and troubleshooting
(missing PDBs, multi-chain/multi-copy ligands, numbering offsets, custom feature families).
