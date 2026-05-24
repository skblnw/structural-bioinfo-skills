# Feature Definitions

This skill uses the standard RDKit pharmacophore feature factory built from
`$RDDataDir/BaseFeatures.fdef`. The fdef file contains SMARTS-based
definitions for each feature family — RDKit applies these on a sanitized
3D molecule and returns one feature per match.

## Default families (BaseFeatures.fdef)

| Family | Meaning | Tolerance (Å) | PML color | Label |
|---|---|---|---|---|
| `Donor` | Hydrogen-bond donor (heavy atom + H) | 1.0 | `blue` | `HBD` |
| `Acceptor` | Hydrogen-bond acceptor (lone-pair heavy atom) | 1.0 | `red` | `HBA` |
| `Aromatic` | Aromatic ring centroid | 1.0 | `orange` | `ARO` |
| `LumpedHydrophobe` | One centroid per hydrophobic region | 1.5 | `forest` | `HYD` |
| `Hydrophobe` | (suppressed — granular per-atom) | — | — | — |
| `PosIonizable` | Positively ionizable (amine, guanidinium) | 1.5 | `marine` | `POS` |
| `NegIonizable` | Negatively ionizable (carboxylate, phosphate) | 1.5 | `firebrick` | `NEG` |
| `ZnBinder` | Zinc-binding group (thiols, hydroxamates) | 1.0 | `purple` | `ZN` |

The skill drops the granular `Hydrophobe` family (one feature per
hydrophobic atom) in favour of `LumpedHydrophobe`, which collapses a
contiguous hydrophobic region into a single centroid feature — this is
what pharmacophore tools (LigandScout, MOE, Pharmer) typically display.

## Tolerance radii

RDKit does **not** assign per-feature tolerance radii — `Pharm3D.Pharmacophore`
only stores pairwise distance bounds between features. The defaults used
here follow the LigandScout / MOE convention:

| Tolerance | Used for |
|---|---|
| 1.0 Å | Directional features: HBD, HBA, Aromatic, ZnBinder |
| 1.5 Å | Larger / less-directional: hydrophobic regions, ionizables |

These values are reasonable starting points for shape-matching virtual
screening. Tightening to ~0.7 Å gives a more selective query (fewer false
hits, more false negatives); loosening to 2.0 Å is rarely useful for
HBD/HBA but acceptable for hydrophobic features.

To change the defaults, edit `TOL_DEFAULTS` in
`scripts/pharmacophore_analysis.py`.

## Why BaseFeatures.fdef and not MinimalFeatures.fdef

- `BaseFeatures.fdef` — the standard pharmacophore set; matches what most
  external tools expect when consuming pharmacophore queries.
- `MinimalFeatures.fdef` — a denser variant intended for similarity work,
  with finer-grained subtypes. Returns more features but breaks
  interoperability with downstream consumers that expect the standard
  family names.

The skill uses `BaseFeatures.fdef`. If you want to swap, change the
`fdef_path` line in `main()`.

## Color choices

The PyMOL color names used in the PML are picked for high contrast on a
white background and to match the LigandScout / Pharmer convention as
closely as PyMOL's named-color palette allows. The `ZnBinder` family uses
`purple` because true metal-binding features are rare and need to stand
out from the more common HBD/HBA spheres.
