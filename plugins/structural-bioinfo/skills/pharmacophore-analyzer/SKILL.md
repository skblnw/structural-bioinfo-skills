---
name: pharmacophore-analyzer
description: This skill should be used when the user asks to "extract pharmacophore from PDB", "compute pharmacophore features", "structure-based pharmacophore", "interaction-filtered pharmacophore", "ligand pharmacophore from co-crystal", "HBD/HBA/aromatic features from a bound ligand", "build a pharmacophore query for virtual screening", "PML pharmacophore visualization", "RDKit pharmacophore analysis", or provides a protein-ligand PDB file (with or without a SMILES) and wants 3D feature coordinates with tolerance radii in JSON plus a PyMOL .pml visual object. Triggers on phrases like "pharmacophore from this co-crystal", "what HBD/HBA does the ligand have", "pharmacophore for VS query", "pharmacophore PML", "BaseFeatures.fdef", "ChemicalFeatures.GetFeaturesForMol", or any task that takes a protein-ligand 3D structure and returns a pharmacophore description. Auto-fetches the ligand's idealized SMILES from the RCSB Chemical Component Dictionary so bond orders are correct, applies an interaction filter against pocket residues by default (LigandScout-style), and writes both a JSON pharmacophore (feature, coordinates, tolerance, partner residue) and a self-contained PML scene (ligand sticks + colored feature spheres).
---

# Pharmacophore Analyzer

## Purpose

Compute a 3D pharmacophore from a protein-ligand complex PDB file using RDKit.
The skill identifies the bound ligand, restores correct bond orders and
hydrogens via a SMILES template (auto-fetched from the RCSB Chemical Component
Dictionary by default), runs the standard `BaseFeatures.fdef` feature factory
on the ligand, and then keeps only those features that engage the protein
(structure-based / interaction-filtered pharmacophore). The output is a
JSON file listing each feature's type, 3D coordinates, tolerance radius and
binding partner, plus a PyMOL `.pml` scene that loads the input PDB and draws
colored feature spheres alongside the ligand.

## When to Use

This skill triggers when the user wants to:

- "Compute a pharmacophore from co-crystal 1IEP"
- "Extract HBD, HBA, aromatic and hydrophobic features from this bound ligand"
- "Build a structure-based pharmacophore query for virtual screening"
- "Generate a PML visualization of the pharmacophore for PDB X"
- "What pharmacophore does imatinib display in the Abl pocket?"
- "Run RDKit ChemicalFeatures on the ligand in this PDB and write a PML"
- "Interaction-filtered pharmacophore for figure 2"

The skill is the natural follow-on to `pdb-extractor` (which produces the
clean monomer + ligand PDB used as input here) and a natural input to
`visualfactory` (which can render the resulting PML to a publication figure).

## Quick Start

> **Script location:** `scripts/pharmacophore_analysis.py` lives in the skill
> directory. Use the full path or `cd` to the skill base before running.
>
> **RDKit auto-discovery:** the script's first action is to make sure RDKit is
> importable. If the active Python doesn't have it, the script searches common
> conda locations (`~/miniconda3`, `~/anaconda3`, `~/miniforge3`, `~/mambaforge`,
> `/opt/conda`, the parent of `$CONDA_PREFIX`) for an env that does — preferring
> one literally named `rdkit` — and re-execs itself under that env's Python.
> So `python scripts/pharmacophore_analysis.py ...` works from any Python on
> the machine as long as *some* conda env has rdkit installed. A loop guard
> prevents infinite re-exec; if no env has rdkit, the script prints a clear
> install hint and exits.

### Auto-detect ligand, auto-fetch SMILES from RCSB (simplest)

```bash
python scripts/pharmacophore_analysis.py 1IEP_clean.pdb
# writes pharmacophore/1IEP_clean.pharmacophore.{json,pml} and symlinks 1IEP_clean.pdb in there
```

### Specify the ligand explicitly + supply SMILES (most accurate, fully offline)

```bash
python scripts/pharmacophore_analysis.py complex.pdb \
    --ligand-resname STI \
    --smiles "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1" \
    --output-prefix imatinib
```

### Skip the interaction filter (return every ligand feature)

```bash
python scripts/pharmacophore_analysis.py 1IEP_clean.pdb --include-all-features
```

### Offline best-effort (no network, no SMILES — accuracy degrades)

```bash
python scripts/pharmacophore_analysis.py complex.pdb --no-fetch
```

## Usage Modes

### 1. Auto-fetch mode (default)

Auto-detects the ligand (largest non-water/non-ion HETATM residue with
≥ 6 heavy atoms), then resolves its SMILES from the RCSB Chemical Component
Dictionary at `https://files.rcsb.org/ligands/download/<RESN>_ideal.sdf`.
Network is required. Use this whenever the bound ligand has a real
PDB-CCD entry (almost all crystallographic ligands do).

### 2. User-supplied template (`--smiles` or `--template-sdf`)

Provide the canonical SMILES (or an SDF file) for the ligand. This is the
fully-offline path and is also the most accurate when the PDB CCD entry
disagrees with the actual covalent state in the structure (e.g. observed
tautomers, charged forms, or atypical protonation).

### 3. Best-effort mode (`--no-fetch`)

Skip the RCSB fetch and run pharmacophore detection on whatever RDKit can
sanitize from the PDB-derived ligand alone. **Bond orders and aromaticity
are usually wrong** for drug-like ligands without a template, which corrupts
Donor/Acceptor/Aromatic detection. The script warns clearly on stderr and
records `"smiles_source": "best_effort"` in the JSON metadata.

### 4. Multi-copy ligands

If the input PDB contains multiple copies of the same ligand (symmetry mates,
multi-chain biological assemblies), each copy is processed independently and
output files are suffixed with the chain letter:
`<prefix>_chainA.pharmacophore.json`, `<prefix>_chainB.pharmacophore.json`,
etc.

## Full CLI Reference

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `pdb_path` | yes | — | Positional: path to the input protein-ligand PDB |
| `--ligand-resname` | no | auto | 3-letter ligand residue code (e.g. `STI`) |
| `--ligand-chain` | no | auto | Restrict to a specific chain |
| `--smiles` | no | — | Canonical SMILES template for bond-order fix |
| `--template-sdf` | no | — | SDF file with the ligand template |
| `--no-fetch` | no | off | Skip RCSB CCD fetch, do best-effort sanitize |
| `--include-all-features` | no | off | Skip interaction filter, emit all ligand features |
| `--output-dir` / `-d` | no | `pharmacophore` | Output directory (auto-created). Input PDB is symlinked here. |
| `--output-prefix` / `-o` | no | PDB stem | File stem within `--output-dir`. Pass an absolute path to bypass `--output-dir`. |
| `--cutoff-hbond` | no | 3.5 | Donor↔Acceptor cutoff in Å |
| `--cutoff-hydrophobic` | no | 4.5 | Hydrophobic↔Hydrophobic cutoff in Å |
| `--cutoff-aromatic` | no | 5.5 | Aromatic↔Aromatic / π-cation cutoff in Å |
| `--cutoff-ionic` | no | 5.0 | PosIonizable↔NegIonizable cutoff in Å |
| `--neighbor-radius` | no | 5.0 | Pocket-residue inclusion radius in Å |

## Output Structure

By default everything lands inside a single `pharmacophore/` subdir of the
current working directory (override with `-d/--output-dir`). The input PDB is
symlinked into that subdir so the PML's relative `load <basename>` resolves
correctly when PyMOL is launched from inside the dir:

```
pharmacophore/
├── <basename>.pdb                         # symlink to the input PDB
├── <prefix>.pharmacophore.json            # feature list (type, position, tolerance, partner)
└── <prefix>.pharmacophore.pml             # PyMOL scene: ligand sticks + feature spheres
```

For multi-copy ligands (symmetry mates / multimers), the JSON/PML files are
suffixed `_chainA`, `_chainB`, etc. — but the single PDB symlink is shared
across all of them since `load <basename>` is the same:

```
pharmacophore/
├── 1IEP.pdb                               # symlink (one)
├── 1IEP_chainA.pharmacophore.json
├── 1IEP_chainA.pharmacophore.pml
├── 1IEP_chainB.pharmacophore.json
└── 1IEP_chainB.pharmacophore.pml
```

To open a result in PyMOL:

```bash
cd pharmacophore
pymol 1IEP_chainA.pharmacophore.pml
```

JSON shape:

```json
{
  "input_pdb": "1IEP_clean.pdb",
  "ligand": {
    "resname": "STI", "chain": "A", "resnum": 1001,
    "smiles_used": "...",
    "smiles_source": "rcsb_ccd"
  },
  "tolerance_defaults_A": {"Donor":1.0,"Acceptor":1.0,"Aromatic":1.0,
                           "LumpedHydrophobe":1.5,"PosIonizable":1.5,
                           "NegIonizable":1.5,"ZnBinder":1.0},
  "interaction_cutoffs_A": {"Donor-Acceptor":3.5,"Aromatic-Aromatic":5.5,
                            "PosIonizable-NegIonizable":5.0,
                            "LumpedHydrophobe-Hydrophobe":4.5,
                            "ZnBinder-Zn":3.0},
  "features": [
    {
      "id": 0,
      "family": "Donor", "type": "SingleAtomDonor",
      "position": [12.345, 8.210, -3.117],
      "tolerance": 1.0,
      "ligand_atom_ids": [3],
      "interaction": {"partner_residue":"ASP381","partner_chain":"A",
                      "partner_atom":"OD2","distance_A":2.78,"type":"hbond"}
    }
  ],
  "metadata": {
    "rdkit_version":"...",
    "fdef":"BaseFeatures.fdef",
    "interaction_filtered": true,
    "n_features_before_filter": 14,
    "n_features_after_filter": 7
  }
}
```

PML shape:

```pml
load 1IEP_clean.pdb, complex
hide everything
show cartoon, polymer
color grey80, polymer
show sticks, resn STI
color cyan, (resn STI and elem C)
util.cnc resn STI

pseudoatom phar_Donor_0, pos=[12.345, 8.210, -3.117], color=blue, label="HBD"
set sphere_scale, 1.0, phar_Donor_0
# ... one block per feature ...

show spheres, phar_*
set sphere_transparency, 0.5, phar_*
set label_size, 14
bg_color white
zoom (resn STI), 8
```

## How It Works (Internal Flow)

1. **Parse PDB with RDKit**: `Chem.MolFromPDBFile(removeHs=False, sanitize=False)`
   then `GetMolFrags(asMols=True)`.
2. **Identify ligand fragment**: HETATM, not in `EXCLUDE_RESIDUES`, ≥ 6 heavy
   atoms; honor `--ligand-resname` / `--ligand-chain` overrides; report all
   matching copies and process each.
3. **Resolve SMILES template**: explicit `--smiles` or `--template-sdf`, else
   auto-fetch from `https://files.rcsb.org/ligands/download/<RESN>_ideal.sdf`,
   else best-effort sanitize (with stderr warning).
4. **Fix bond orders + add hydrogens**: `AllChem.AssignBondOrdersFromTemplate`
   then `Chem.AddHs(addCoords=True)` and `Chem.SanitizeMol`.
5. **Build feature factory**: `ChemicalFeatures.BuildFeatureFactory(
   $RDDataDir/BaseFeatures.fdef)`. Drop raw `Hydrophobe`; keep
   `LumpedHydrophobe` (one centroid per hydrophobic region).
6. **Compute ligand features**: `factory.GetFeaturesForMol(ligand)` →
   family / type / position / contributing atom indices.
7. **Build protein-side feature points**: textually re-parse the PDB
   `ATOM` lines within `--neighbor-radius` of any ligand heavy atom and
   emit feature points using a residue-atom-name lookup (Ser/Thr/Tyr
   side-chain hydroxyls, Asn/Gln amides, Asp/Glu carboxylates, Lys NZ,
   Arg guanidinium, His imidazole, backbone N/O for HBD/HBA, Phe/Tyr/Trp/His
   ring centroids for Aromatic, Leu/Ile/Val/Met/Ala/Pro/Cys side-chain
   atoms for Hydrophobe). The full table lives in
   `references/interaction_criteria.md`.
8. **Interaction filter**: each ligand feature is kept only if a
   complementary protein feature is within the cutoff for that pair-type.
   `--include-all-features` skips this step.
9. **Tolerance assignment**: from `references/feature_definitions.md`
   defaults (1.0 Å for directional features, 1.5 Å for hydrophobic /
   ionizable). Per-family overrides are easy to add inline if needed.
10. **Write outputs**: `<prefix>.pharmacophore.json` and
    `<prefix>.pharmacophore.pml`.

## Defaults & Conventions

Tolerance radii (Å) — LigandScout / MOE convention:

| Family | Tolerance |
|---|---|
| Donor (HBD) | 1.0 |
| Acceptor (HBA) | 1.0 |
| Aromatic | 1.0 |
| LumpedHydrophobe | 1.5 |
| PosIonizable | 1.5 |
| NegIonizable | 1.5 |
| ZnBinder | 1.0 |

PML colors:

| Family | PyMOL color |
|---|---|
| Donor | `blue` |
| Acceptor | `red` |
| Aromatic | `orange` |
| LumpedHydrophobe | `forest` |
| PosIonizable | `marine` |
| NegIonizable | `firebrick` |
| ZnBinder | `purple` |

Interaction cutoffs (Å) used by the structure-based filter:

| Ligand family | Protein partner | Cutoff |
|---|---|---|
| Donor | Acceptor | 3.5 |
| Acceptor | Donor | 3.5 |
| PosIonizable | NegIonizable | 5.0 |
| NegIonizable | PosIonizable | 5.0 |
| Aromatic | Aromatic / PosIonizable | 5.5 |
| LumpedHydrophobe | Hydrophobe | 4.5 |
| ZnBinder | Zn (HETATM) | 3.0 |

## Troubleshooting

**`ImportError: rdkit`**
The script auto-discovers a conda env with rdkit and re-execs there, so this
error means *no* conda env on the machine has rdkit. Install with
`conda install -c conda-forge rdkit -n rdkit` (or any env name) or
`pip install rdkit` into the active Python. To force a specific Python and
skip discovery, just invoke that Python directly:
`/path/to/python scripts/pharmacophore_analysis.py ...`

**"No ligand fragment detected"**
Auto-detection requires a HETATM residue with ≥ 6 heavy atoms that is not in
the exclusion list. Pass `--ligand-resname <CODE>` to override. If your input
is a ligand-only PDB (no protein), the skill still works — pass
`--include-all-features` to skip the interaction filter that would otherwise
remove everything.

**"RCSB CCD fetch failed"**
Network may be down or the ligand has no PDB-CCD entry. Pass `--smiles "..."`
or `--template-sdf path.sdf` instead, or `--no-fetch` to proceed best-effort.

**"AssignBondOrdersFromTemplate failed"**
The PDB ligand and the SMILES template have different heavy-atom counts,
typically because of alternate locations, partial occupancy, or a mismatch
between the deposited ligand and the canonical CCD entry. Inspect the raw
ligand fragment, supply a custom SMILES that matches the deposited atoms,
or use `--no-fetch`.

**Covalent inhibitor merged into protein**
RDKit groups covalently bonded ligands with the protein chain; auto-detection
won't see them. Either supply an unbound SMILES + `--ligand-resname` and the
script will reach the residue by name, or first run a separate cleanup pass
that breaks the covalent bond.

**Multiple ligand copies (symmetry mates / multimer)**
The script processes each copy and writes per-chain outputs. Use
`--ligand-chain A` to restrict to a single copy.

**Feature count looks low**
The default behaviour is the **interaction-filtered** structure-based mode —
ligand features that don't engage the protein within the cutoffs are dropped.
Pass `--include-all-features` to see the full ligand-only pharmacophore.

## Workflow Guidance for the Assistant

The natural pipeline is:

```
pdb-holostructure-search → pdb-extractor → pharmacophore-analyzer → visualfactory
```

When invoked alongside `pdb-extractor`, point this skill at the cleaned PDB
that `pdb-extractor` writes (`raw/<pdb_id>/<pdb_id>_clean.pdb`). The cleaned
file already has only the target chain plus its bound ligands, which is the
ideal input — protein-side feature detection is restricted to the same chain,
so a clean PDB removes ambiguity from multi-copy assemblies.

After running, the resulting `.pml` is suitable for direct use with the
`visualfactory` skill (e.g., `MasterOfMorphing` or any scene that takes a
PML scene file) to produce a 4K figure. The JSON file can be passed to
downstream pharmacophore screening tools (e.g., `Pharmer`, `align-it`) or
loaded into Python for custom shape comparison.

## Additional Resources

### Scripts

- **`scripts/pharmacophore_analysis.py`** — Single entry-point script.
  Runs the full pipeline: PDB parse → ligand identification → SMILES
  template fetch → bond-order fix → RDKit ChemicalFeatures → protein-side
  feature inference → interaction filter → JSON + PML output. Run with
  `--help` for the full flag list.

### References

- **`references/feature_definitions.md`** — RDKit feature families, tolerance
  radii defaults, color conventions, and notes on `BaseFeatures.fdef` vs
  `MinimalFeatures.fdef`.
- **`references/interaction_criteria.md`** — Distance cutoffs per pair-type
  and the residue-atom-name → feature lookup table used to compute
  protein-side feature points.
