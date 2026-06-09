# Input format & configuration reference

This skill renders the outputs of a **structure-based pharmacophore detector** (e.g. RDKit
`BaseFeatures.fdef` with an interaction filter). It does not detect features itself.

## 1. Directory layout

Point the skill at the directory holding the per-structure pharmacophore JSON files. Two common
layouts are both supported (auto-discovered recursively):

```
# Flat
<target>/pharmacophore/
    9L45_filtered.pharmacophore.json
    9L45_all.pharmacophore.json
    9L46_filtered.pharmacophore.json
    ...
<target>/raw/9L45/9L45_clean.pdb          # receptor referenced by input_pdb

# Per-structure subfolders (also fine)
<target>/pharmacophore/
    8U4P/8U4P_chainR_filtered.pharmacophore.json
    8U4P/8U4P_chainR_all.pharmacophore.json
    consensus_pharmacophore.json
<target>/raw/8U4P/8U4P_clean.pdb
```

- The leading token of each filename is treated as the **structure id** (`9L45`, `8U4P`, `22XC`).
- `<prefix>_filtered.pharmacophore.json` is required; the matching `_all` file is optional (used
  for the all-vs-retained retention chart).
- Duplicate representations of the same ligand copy (same structure id + ligand resname + chain +
  resnum) are de-duplicated; the shallower path wins.
- The cleaned receptor PDB is located from each JSON's `input_pdb` field, resolved against
  `--raw-dir`, the parent of the pharmacophore dir, and a recursive basename search. **Broken
  symlinks are skipped** (pharmacophore dirs often symlink `*_clean.pdb` to `../raw/...`).

## 2. Per-structure pharmacophore JSON schema

```json
{
  "input_pdb": "raw/9L45/9L45_clean.pdb",
  "ligand": {
    "resname": "A1E", "chain": "A", "resnum": 2701,
    "smiles_used": "O=c1ccc2...", "smiles_source": "template_failed_best_effort"
  },
  "features": [
    {
      "id": 0, "family": "Donor", "type": "SingleAtomDonor",
      "position": [137.396, 211.369, 187.255], "tolerance": 1.0,
      "ligand_atom_ids": [27],
      "interaction": {
        "partner_residue": "ASP2494", "partner_chain": "A",
        "partner_atom": "OD1", "distance_A": 3.297, "type": "hbond"
      }
    }
  ],
  "metadata": { "n_features_before_filter": 17, "n_features_after_filter": 4 }
}
```

- `family` is one of `Donor, Acceptor, PosIonizable, NegIonizable, Aromatic, LumpedHydrophobe`.
- `position` is the ligand-side feature coordinate (Å) in the **same frame** as `input_pdb`. The
  skill renders a sphere here and draws a dashed line to `interaction.partner_atom` (looked up in
  the receptor by residue + atom name), so partners must exist in the PDB.
- `interaction.type` (`hbond`, `salt_bridge`, …) colors the dashed line.
- `metadata.n_features_*` drive the retention counts; if absent they are inferred from the
  feature lists.

Each entry's binding pocket = the ligand heavy atoms + every protein residue with an atom within
`--pocket-cutoff` Å (default 5.0) of the ligand. Hydrogens are dropped to keep the embedded PDB
small.

## 3. Consensus query schema (optional)

A file matching `*consensus*.json` in the directory (or `--consensus FILE`) is embedded and
overlaid in 3D on its alignment-reference structure.

```json
{
  "target": {"name": "...", "gene": "CXCR4", "uniprot": "P61073"},
  "alignment": {"reference": "8ZPN chain R C-alpha backbone"},
  "source_entries": [{"pdb_id": "8U4P", "ligand": "VH6", "n_filtered_features": 12}],
  "query_summary": {"total_features": 22, "mandatory": 11, "optional": 11},
  "features": [
    {"id": 0, "family": "Donor", "position": [112.95,133.73,122.42],
     "tolerance": 1.0, "classification": "mandatory", "n_sources": 2}
  ]
}
```

- Each feature's `classification` is `mandatory` or `optional`. If absent, a boolean
  `required`/`mandatory` field is honored; otherwise it defaults to `optional`.
- The overlay frame is chosen from `alignment.reference` (matched to a structure id present in the
  dataset); the consensus coordinates must be in that structure's frame.

## 4. Report config schema (optional, `--config`)

All keys optional; omitted values are auto-generated.

| Key | Type | Purpose |
|-----|------|---------|
| `accent`, `accent_dark`, `accent_soft` | hex string | Theme colors |
| `kicker`, `meta` | string | Header eyebrow text / provenance line |
| `target` | object | `{name, gene, uniprot, organism, n_structures}` |
| `entry_overrides` | array | Per-structure/entry labels + classification (below) |
| `summary_html` | HTML string | Executive summary (else auto-generated) |
| `modes` | array | Query-recommendation cards (below) |
| `strategies` | array | `{name, criteria, behavior}` cards (shown with a consensus query) |
| `interaction_map` | string | Monospace global interaction map (preformatted) |
| `interaction_map_note` | HTML string | Note under the map |
| `caveat` | HTML string | Highlighted callout (e.g. residue-numbering offset) |
| `limitations` | array | `{issue, impact, mitigation}` rows |
| `manifest` | string | File manifest (preformatted) |
| `footer` | HTML string | Footer text |

### `entry_overrides[]`

Keyed by `pdb_id` (applies to all that structure's entries) or `key` (one ligand copy = the
filtered filename without `.pharmacophore.json`, e.g. `9L40_filtered_chainA_res2702`).

```json
{"pdb_id": "8U4N", "chemotype_label": "CLR (Cholesterol)",
 "role": "Allosteric/membrane lipid", "category": "excluded"}
```

- `category`: `orthosteric` / `binding` / `productive` → counts toward retention & residue
  aggregation; `excluded` / `control` / `negative` → kept in the table but dropped from the
  aggregates (built-in negative control).
- `label`, `chemotype_label`, `role`, `note` are display-only.

### `modes[]` (query recommendation cards)

```json
{"name": "Primary query — ASP2494 / GLY2385", "confidence": "Moderate confidence",
 "conf_cls": "mod", "support": "9L40 + 9L45", "desc": "…",
 "features": [{"fam": "Donor", "tol": "1.0", "partner": "ASP2494 OD1"}],
 "strategy": "Require 3 of 4 features; the Donor/PosIonizable pair is mandatory."}
```

`conf_cls` ∈ `mod | low | high` (badge color). `features` and `strategy` are optional.

## 5. Troubleshooting

- **`no clean PDB for <id>`** — set `--raw-dir` to the folder that contains the cleaned receptors;
  the 3D viewer is disabled for entries whose PDB can't be found (tables/charts still render).
- **Incidental binders inflate the stats** — mark lipids/cholesterol/buffer ligands
  `category: excluded` in `entry_overrides`.
- **Multi-chain / multi-copy ligands** — handled automatically; each chain/copy is a separate
  entry, labeled `<PDB> · chain <X>`.
- **Residue-numbering offsets between structures** — the skill reports residues verbatim from the
  JSON. Explain the mapping with a `caveat` in the config (e.g. a construct-truncation offset).
- **Custom feature families** — the six standard RDKit families are styled. Others render in
  tables but have no dedicated sphere color.
- **Offline 3D** — if `3Dmol.js` can't load, the viewer shows a clear message and the rest of the
  report is unaffected; reconnect and reload to view structures.
