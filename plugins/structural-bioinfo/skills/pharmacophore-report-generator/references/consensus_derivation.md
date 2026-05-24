# Consensus Pharmacophore Derivation — Method Reference

This file documents the algorithm `extract_report_data.py --derive-consensus`
uses to build a multi-structure consensus query from per-PDB filtered
pharmacophore JSONs. Read this when you need to explain the method in §4.2 /
§5 of the report, or when you're tuning the clustering threshold for a
receptor where the default doesn't fit.

## Overview

The consensus pharmacophore answers the question: *"Across N orthosteric
holostructures of the same receptor, which feature positions are conserved
enough to demand of any future hit, and which are nice-to-haves?"*

The derivation has four stages:

1. **Pick a reference structure** — the orthosteric PDB with the most
   filtered features. A feature-rich reference reduces the chance that a
   conserved feature in one structure has no nearby counterpart in the
   reference frame.
2. **Kabsch-align** every other orthosteric structure's Cα backbone to the
   reference. Use only Cα atoms common to both structures (matched by chain
   + residue number).
3. **Transform features** of every aligned structure into the reference
   frame using the rotation/translation from step 2.
4. **Cluster** features by family using greedy single-linkage with cutoff =
   `cluster_factor × tolerance` (default factor 2.0). Classify each cluster
   as **mandatory** if its source PDBs ≥ ⌈N/2⌉, else **optional**.

## Why Kabsch and not RMSD-fit-the-ligand?

Aligning ligands directly would mask the very thing we want to measure: how
much the ligand pose drifts because the receptor flexes. Receptor Cα
alignment fixes the protein's reference frame so feature drift reflects
real ligand-pose differences, not coordinate-system drift.

For an N-residue Cα set common to two structures (centered at their
centroids), the Kabsch rotation R minimizes ‖p·R − q‖² and is given by:

```
H = (p_centered).T @ q_centered          # 3×3 covariance
U, S, V.T = SVD(H)
d = sign(det(V.T @ U.T))                  # handedness fix to avoid reflection
R = V.T @ diag(1, 1, d) @ U.T
t = q_centroid − p_centroid @ R
```

Apply `aligned_p = p @ R + t` to transform any point from the source frame
to the target frame. The script reuses `R, t` from each pairwise alignment
to transform that source's feature coordinates as well.

## Why greedy single-linkage clustering?

Two reasons:

1. **It matches user intuition.** Two pharmacophore features are "the same
   point" iff one falls inside the other's tolerance sphere. Greedy
   single-linkage propagates this: feature C joins {A, B} if C is within
   tolerance of *any* member of {A, B}.
2. **It's order-independent for typical cases.** For sparse,
   well-separated features (the common case in orthosteric pockets),
   greedy and exact hierarchical clustering produce identical clusters. For
   dense regions (e.g., when 4+ features overlap on the same ligand atom),
   borderline ambiguities may emerge — see "Limitations" below.

The default threshold `2 × tolerance` is calibrated empirically: it groups
features that are obvious matches (1.0–2.5 Å apart for Donor, 1.5–3.0 Å for
PosIonizable) without merging chemically distinct features.

## Source-count classification

Let N = number of orthosteric source PDBs. A cluster's classification is:

```
mandatory   if  n_sources >= ceil(N / 2)
optional    otherwise
```

Examples:

| N (sources) | mandatory threshold | rationale |
|---|---|---|
| 2 | 2 | Both must agree — a single source isn't enough to anchor a query |
| 3 | 2 | Majority rule (CXCR4 case) |
| 4 | 2 | Half-and-half is enough |
| 5 | 3 | Majority rule |
| 6 | 3 | |
| 7 | 4 | |

For N=2, mandatory features are exactly the intersection: both structures
must contribute a feature to the cluster. This is strict but appropriate —
two-structure consensus is otherwise too brittle.

## Tuning the clustering threshold

The default `--cluster-factor 2.0` works for most receptor families. Tune it
when:

- **Receptor has many close-spaced amine groups** (polyamines like
  AMD3100 / Plerixafor): may merge features that should stay separate. Try
  `--cluster-factor 1.5`.
- **Receptor has very flexible loops** (e.g., kinase activation loop) that
  swing 5+ Å between structures: may fail to merge same-residue features
  across activation states. Try `--cluster-factor 3.0`. But beware:
  larger factors drift toward "everything merges into one feature".
- **Receptor has only N=2 sources**: stick with 2.0. The mandatory threshold
  is already 2 (= both sources), so loosening clustering only adds noise.

After re-clustering, eyeball the consensus PML in PyMOL: clusters should
look like spheres, not ellipsoids of weird shape. If a cluster spans more
than ~3 Å center-to-edge for a Donor (1.0 Å tolerance), the cluster is
suspect.

## Output schema (canonical)

The consensus JSON written by this skill has this shape (top-level
`classification` field, not nested under `conservation` — the latter is a
legacy schema some older projects use, and the script reads both):

```json
{
  "description": "<receptor> orthosteric pharmacophore consensus for virtual screening",
  "target": "<receptor> (UniProt <accession>)",
  "source_entries": [
    {"pdb_id": "...", "ligand": "...", "n_filtered_features": <int>}
  ],
  "excluded_entries": [
    {"pdb_id": "...", "ligand": "...", "reason": "..."}
  ],
  "alignment": [
    {"pdb_id": "<reference>", "role": "reference", "common_ca": ..., "rmsd_before_A": 0, "rmsd_after_A": 0},
    {"pdb_id": "<other>", "role": "aligned",   "common_ca": ..., "rmsd_before_A": ..., "rmsd_after_A": ...}
  ],
  "method": "Kabsch alignment of N orthosteric Cα + greedy single-linkage clustering (2.0× tolerance)",
  "tolerances_A": {"Donor": 1.0, "Acceptor": 1.0, "PosIonizable": 1.5, ...},
  "family_summary": {
    "<family>": {"mandatory": <int>, "optional": <int>}
  },
  "features": [
    {
      "id": <int>,
      "family": "Donor|Acceptor|PosIonizable|...",
      "type": "<RDKit feature type>",
      "position": [x, y, z],
      "tolerance": <float>,
      "radius_from_centroid": <float>,
      "classification": "mandatory|optional",
      "n_sources": <int>,
      "n_instances": <int>,
      "members": [
        {"pdb": "...", "id": <int>, "position": [x, y, z]}
      ]
    }
  ],
  "query_summary": {
    "total_features": <int>,
    "mandatory": <int>,
    "optional": <int>,
    "screening_strategy": {
      "strict": "Match all <N> mandatory features",
      "relaxed": "Match >=ceil(N/2) mandatory + >=3 optional",
      "scoring": "Score = mandatory_matches + 0.5 * optional_matches; threshold >=0.75 * total_mandatory"
    }
  }
}
```

## Limitations

- **Order-dependence in dense clusters.** When 4+ features are within
  cutoff of each other, the greedy assignment may produce slightly
  different clusters depending on iteration order. The script processes
  features in source-PDB order then feature-id order for stability.
- **No outlier handling.** A feature 4 Å from any other in its family
  becomes a singleton cluster. This is correct behavior for a heterogeneous
  pocket but can inflate the optional-feature count when one structure has a
  genuinely anomalous pose.
- **No per-feature tolerance scaling.** All clusters inherit the source
  tolerance; if a cluster spans many sources, you might want to widen the
  tolerance to capture the spread. The script reports the per-cluster
  `radius_from_centroid` so a downstream tool can do this.
- **Cross-receptor consensus is out of scope.** The Kabsch alignment
  assumes both structures are the same protein with shared Cα residue
  numbering.
