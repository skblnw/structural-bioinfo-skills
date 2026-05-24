# Interaction Criteria

The skill computes a **structure-based** pharmacophore by default: ligand
features are kept only when a complementary functional group on the protein
side lies within a type-specific distance cutoff. This document records the
cutoffs used and the residue-atom-name lookup that turns a PDB into a list
of protein-side feature points.

## Distance cutoffs

These are heavy-atom-to-heavy-atom (or centroid-to-centroid for aromatic)
distances. Angles are not enforced in v1 of the skill — the implementation
uses simple distance criteria, which matches what tools like PLIP use as a
coarse filter.

| Ligand family | Protein partner | Cutoff (Å) | Interaction type |
|---|---|---|---|
| `Donor` | `Acceptor` | 3.5 | `hbond` |
| `Acceptor` | `Donor` | 3.5 | `hbond` |
| `PosIonizable` | `NegIonizable` | 5.0 | `salt_bridge` |
| `NegIonizable` | `PosIonizable` | 5.0 | `salt_bridge` |
| `Aromatic` | `Aromatic` | 5.5 | `pi_stack` |
| `Aromatic` | `PosIonizable` | 5.5 | `pi_cation` |
| `LumpedHydrophobe` | `Hydrophobe` | 4.5 | `hydrophobic` |
| `ZnBinder` | `Zn` (HETATM) | 3.0 | `metal` |

The cutoffs are configurable per run via `--cutoff-hbond`,
`--cutoff-hydrophobic`, `--cutoff-aromatic`, `--cutoff-ionic`. The pocket
radius (which atoms count as "near the ligand") is `--neighbor-radius`,
default 5.0 Å.

## Protein-side feature lookup

For each ATOM record within the pocket, the skill emits one or more feature
points using the table below. Backbone N/O contribute HBD/HBA features for
every residue regardless of identity. Aromatic rings emit a single feature
at the ring centroid (also reused as a hydrophobic point).

### Side-chain atoms

| Residue | Atom | Features |
|---|---|---|
| SER | OG | Donor, Acceptor |
| THR | OG1 | Donor, Acceptor |
| TYR | OH | Donor, Acceptor |
| ASN | OD1 | Acceptor |
| ASN | ND2 | Donor |
| GLN | OE1 | Acceptor |
| GLN | NE2 | Donor |
| ASP | OD1 / OD2 | Acceptor, NegIonizable |
| GLU | OE1 / OE2 | Acceptor, NegIonizable |
| LYS | NZ | Donor, PosIonizable |
| ARG | NH1 / NH2 / NE | Donor, PosIonizable |
| HIS | ND1 / NE2 | Donor, Acceptor, PosIonizable |
| TRP | NE1 | Donor |
| MET | SD | Hydrophobe |
| ALA | CB | Hydrophobe |
| VAL | CB | Hydrophobe |
| LEU | CG | Hydrophobe |
| ILE | CG1 | Hydrophobe |
| PRO | CG | Hydrophobe |
| CYS | SG | Hydrophobe |

### Backbone atoms (every residue)

| Atom | Features |
|---|---|
| N | Donor |
| O | Acceptor |
| OXT | Acceptor (C-terminal carboxylate; also acts as NegIonizable) |

### Aromatic ring centroids

For PHE, TYR, TRP, HIS within the pocket, the skill computes the centroid
of the indicated ring atoms and emits one `Aromatic` feature plus one
`Hydrophobe` feature at the centroid:

| Residue | Ring atoms used |
|---|---|
| PHE | CG, CD1, CD2, CE1, CE2, CZ |
| TYR | CG, CD1, CD2, CE1, CE2, CZ |
| HIS | CG, ND1, CD2, CE1, NE2 |
| TRP | CG, CD1, NE1, CE2, CD2 (5-membered ring) |

Note: TRP has both a 5-membered and a 6-membered aromatic ring. The skill
currently uses only the 5-membered ring. For most pharmacophore purposes
this is sufficient since the 6-membered ring centroid is within ~1.5 Å of
the 5-membered one.

## Why distance-only (no angles)

Angle filtering (e.g. donor-H...acceptor angle ≥ 120° for H-bonds, ring-normal
deviations for π-π) is more rigorous but adds a meaningful amount of
parsing logic for marginal accuracy improvement on the structure-based-filter
step. The skill is conservative on the distance side (slightly tight cutoffs)
to compensate. Add angle filtering in a v2 if the user reports false-positive
features.

## Calibration notes

These cutoffs match the defaults used by:
- **PLIP** (Adasme et al., 2021): H-bond 3.9 Å, salt bridge 5.5 Å, π-stacking 7.5 Å (centroid). PLIP is more permissive; we use slightly tighter values to keep the pharmacophore concise.
- **LigandScout**: H-bond 3.5 Å, salt bridge 5.0 Å, π-stacking 5.5 Å. We match LigandScout.
- **PoseView / 2D ligand interaction diagrams**: same as LigandScout.

If you find the filter is too aggressive (drops features that should clearly
be there), bump the cutoffs by 0.5 Å. If too permissive (keeps features that
look unphysical), tighten by 0.5 Å.
