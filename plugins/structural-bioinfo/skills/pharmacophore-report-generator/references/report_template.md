# Pharmacophore Report — Canonical Template

This file is the source of truth for report structure, table formats, and tone.
Read it before composing the report. Placeholders are wrapped in `{{...}}`.

The data values come from `<project>/.report_data.json` (output of
`scripts/extract_report_data.py`). Don't compute counts in your head — pull
them from the JSON.

---

## Top of file

```markdown
# {{Receptor Name}} Orthosteric Pharmacophore for Virtual Screening

**Report generated:** {{YYYY-MM-DD}} | **RDKit:** {{version, e.g. 2026.03.1}} | **Feature definitions:** `BaseFeatures.fdef`
**Target:** {{Receptor Name}} (UniProt {{accession}})
**Goal:** Derive a validated structure-based pharmacophore query for prospective virtual screening

---
```

The header block is fixed. If RDKit version isn't known, look in any per-PDB
JSON's `metadata.rdkit_version`. The "Goal" line stays as written — this is
always the deliverable.

---

## Executive Summary

Always present. Goes immediately under the header, before §1.

Structure: one prose paragraph framing the dataset → chemotype table →
"Key finding" prose paragraph → one-line pointer to the deliverable JSON.

```markdown
## Executive Summary

{{N}} {{Receptor}} holostructures from the RCSB PDB ({{method, e.g. all cryo-EM, 2.7–3.4 Å}}) were analyzed to derive a structure-based pharmacophore for virtual screening. {{N_orthosteric_chemotypes}} orthosteric chemotypes anchor the query; {{any excluded ligand summary, e.g. cholesterol (CLR, X copies across Y structures) and a co-bound lipid (D21)}} are excluded by the interaction filter and serve as built-in negative controls.

| Chemotype | Structures | Role | Filtered features | Bond-order confidence |
|-----------|-----------|------|-------------------|------------------------|
| **{{ligand_a}}** ({{drug name if known}}) | {{pdb_ids}} | Orthosteric inhibitor | {{counts}} | {{High — CCD template matched / Lower — smiles_used null}} |
| **{{ligand_b}}** | {{pdb_ids}} | Orthosteric inhibitor | {{counts}} | {{...}} |
| {{excluded_ligand}} | {{N structures}}, {{M copies}} | {{Allosteric / membrane / Crystallographic lipid / etc.}} | {{0–N}} | — (excluded) |

**Key finding for screening:** {{One sentence stating the dominant feature families (e.g. "Donor + PosIonizable–dominated pharmacophore")}} anchored on {{the top 1–3 conserved residues with their secondary-structure context: e.g. "ASP270 in TM6 and GLU296 in TM7"}}. {{Numbering caveat sentence if present in numbering_offsets, e.g. "(= ASP262 in 8U4P, an 8-residue construct numbering offset)"}} {{One sentence on what makes the pharmacophore robust — conformational range, chemotype diversity, etc.}}

The deliverable is `{{path to consensus_pharmacophore.json}}` — a {{N}}-feature query ({{n_mandatory}} mandatory, {{n_optional}} optional) directly importable into Pharmer / align-it / RDKit shape-screening tools.

---
```

Keep the executive summary under ~250 words. The "Key finding" sentence is
the most important line in the entire report — make it commit to a specific
claim; don't hedge.

---

## §1. Dataset and Preprocessing

```markdown
## 1. Dataset and Preprocessing

{{N}} {{Receptor}} holostructures from the RCSB PDB were processed with `pdb-extractor` to isolate one {{Receptor}} monomer (chain {{primary chain, e.g. R}}) with all bound ligands. Multi-copy ligands from symmetry mates or multimeric assemblies are preserved per-chain to avoid discarding valid binding poses.

| PDB ID | Ligand(s) | Copies | Source |
|--------|-----------|--------|--------|
| {{pdb}}   | {{ligand_resname}}       | {{n_copies (chain list if multi-chain)}} | {{Cryo-EM / X-ray / —}} |
| {{...}}   | {{...}}       | {{...}}                    | {{...}} |

---
```

Sources go in display order matching the metadata CSV. Use "Cryo-EM" /
"X-ray" labels only for the orthosteric source PDBs; leave other rows as `—`
unless the user supplies the data.

---

## §2. Pharmacophore Feature Detection

```markdown
## 2. Pharmacophore Feature Detection

Performed with the `pharmacophore-analyzer` skill: RDKit `BaseFeatures.fdef` factory on the hydrogen-enriched and bond-order-corrected ligand (SMILES template auto-fetched from the RCSB Chemical Component Dictionary).

### Interaction Filter

Each ligand feature is retained only if a **complementary protein-side feature** (derived from residue-atom-name lookup within 5.0 Å of any ligand heavy atom) falls within the family-pair distance cutoff:

| Ligand Feature      | Protein Partner     | Cutoff (Å) |
|---------------------|---------------------|------------|
| Donor (HBD)         | Acceptor (HBA)      | 3.5        |
| Acceptor (HBA)      | Donor (HBD)         | 3.5        |
| PosIonizable        | NegIonizable        | 5.0        |
| NegIonizable        | PosIonizable        | 5.0        |
| Aromatic            | Aromatic, PosIonizable | 5.5     |
| LumpedHydrophobe    | Hydrophobe          | 4.5        |

### Feature Tolerances (LigandScout/MOE convention)

| Family             | Radius (Å) |
|--------------------|------------|
| Donor, Acceptor, Aromatic | 1.0 |
| LumpedHydrophobe, PosIonizable, NegIonizable | 1.5 |

---
```

Both tables are fixed (they describe the analysis convention, not
receptor-specific data). Don't modify unless the per-PDB JSON's
`interaction_cutoffs_A` or `tolerance_defaults_A` differ from these — in
which case, replace with whatever the JSON says.

---

## §3. Ligand Classification

```markdown
## 3. Ligand Classification: Orthosteric vs Non-Orthosteric

The interaction filter quantitatively separates orthosteric binders from non-orthosteric molecules based on the fraction of retained features.

### 3.1 Results by Ligand

| Ligand | Entries (copies) | All-Ligand Features | Filtered Features | Mean Retention |
|--------|-----------------|---------------------|-------------------|----------------|
| **{{ortho_lig_a}}** | {{N entries}} × {{copies}} | {{all_total}}                    | {{filtered_total}}             | {{XX%}}            |
| {{excluded_lig}} | {{...}} | {{...}} | {{...}} | {{X%}} |

(Sort orthosteric ligands first; bold their resnames.)

### 3.2 Exclusion Rationale

For each excluded ligand, supply 2–3 numbered rationales drawing on:
1. **Quantitative** — retention statistics (e.g., "11 of 16 copies retain zero filtered features").
2. **Feature-type mismatch** — what features the ligand has vs. what the orthosteric pocket demands (cite §4 of this report).
3. **Biological** — known binding site (CCM for cholesterol, phospholipid groove for D21, etc.).

### 3.3 Retained for Query Derivation

One-paragraph statement of which ligands form the orthosteric set, with PDB IDs.

---
```

The exclusion rationale is the most prose-heavy section. Lean on the data:
specific feature counts ("retains 0 filtered features in CXCR4 chain R"),
specific known biology ("cholesterol binds the canonical CCM groove on the
lateral surface"), and side-by-side comparisons when one structure has both
an orthosteric and a non-orthosteric ligand co-bound (this is a powerful
argument).

---

## §4. Orthosteric Pharmacophore Derivation

This is the longest section. Subsections §4.1–§4.5 use derived counts and
alignment math; §4.6–§4.7 surface the per-residue partner data.

### §4.1 Per-Entry Feature Profiles

Per-family count table for each orthosteric source PDB. Bold the dominant families.

### §4.2 Spatial Alignment

Two- or three-row table: "Alignment Pair | Common Cα Residues | RMSD Before (Å) | RMSD After (Å)". One row per non-reference orthosteric PDB. Add 1–2 prose sentences interpreting the RMSD (does it reflect activation-state divergence? cross-chemotype packing differences?).

### §4.3 Aligned Pharmacophore Similarity

Per-pair table of aligned-feature distances: "Pair | Matched Features | Mean Distance | Min | Max". Bold the mean. Add one prose sentence highlighting the "key finding" (e.g., "despite a 7.2 Å protein conformational shift, features align within 1.99 Å mean").

### §4.4 Per-Feature Match (one pair)

Pick the most informative pair (typically the same-ligand cross-state pair, or two structurally diverse chemotypes that converge). Show every feature in a table: "Feature | Family | Closest Match | Distance". Conclude with whether all features have a counterpart and what that means.

### §4.5 Consensus Feature Types

Per-family check-mark table across all orthosteric entries: "Family | {{pdb_a}} | {{pdb_b}} | ... | Consensus". The Consensus column says "{{n_min–n_max}}" or "absent". One prose paragraph naming the dominant families and tying it to the receptor's known ligand preferences.

### §4.6 Receptor Residue Engagement

**This is the most distinctive section** — the per-PDB JSONs already contain
partner_residue/partner_atom/distance/interaction_type, but this data is
almost never reported. Surface it aggressively.

For each orthosteric source PDB, write a sub-subsection with:

```markdown
#### {{pdb}} ({{ligand}}, chain {{X}}, {{n}} features)

| # | Family | Type | Partner | Atom | Dist (Å) | Anchor |
|---|--------|------|---------|------|:--------:|--------|
| 0 | Donor | hbond | **{{RES}}{{NUM}}** | {{atom}} | {{dist}} | {{TM6 / TM7 / ECL2 / hinge / etc.}} |
| ... | ... | ... | ... | ... | ... | ... |

{{One-paragraph interpretation: which residues absorb the most contacts, anything notable about feature distribution.}}
```

The "Anchor" column needs receptor-class context — TM/ECL/loop labels for
GPCRs, hinge/DFG/αC for kinases, H3/H12/AF-2 for nuclear receptors. If the
class is unknown, use "primary acidic anchor" / "secondary anchor" /
"peripheral acidic" instead.

After the per-PDB tables, add an aggregate residue engagement table:

```markdown
#### Aggregate residue engagement

| Residue (per-PDB numbering) | Role | Engaged in | Contacts (Donor + PosIon) |
|---|---|---|:---:|
| **{{primary_anchor}}** ({{aliases if any}}) | {{primary anchor descriptor}} | {{all 3 entries}} | **{{N}}** ({{ND D + NP P}}) |
| {{secondary_anchor}} ({{aliases}}) | {{descriptor}} | {{entries}} | {{N}} |
| ... | | | |

> **Numbering caveat.** {{If numbering_offsets is non-empty: explain that PDB X uses {{ALIAS_A}} while PDBs Y, Z use {{ALIAS_B}}, with offset = {{N}} residues — same conserved residue, different construct numbering. Recommend mapping to the canonical UniProt accession in any downstream cross-reference.}}
```

### §4.7 Global Interaction Map

Always present. ASCII map summarizing residue × feature family across all
orthosteric structures.

```markdown
### 4.7 Global Interaction Map

\`\`\`
{{REGION}}  {{RES_NAME}} ({{aliases if any}})  ← [Donor × {{N}}]   [PosIonizable × {{M}}]   ←  {{which entries}}
{{REGION}}  {{RES_NAME}}                         ← [Donor × {{N}}]   [PosIonizable × {{M}}]   ←  {{which entries}}
{{REGION}}  {{RES_NAME}}                         ← [Acceptor]                                  ←  {{which entry, if any}} — sole Acceptor partner
\`\`\`

{{One prose paragraph interpreting the map: which residues are universal vs which are chemotype-specific, what this implies for screening prioritization.}}

---
```

Sort rows by total contact count (most-engaged residue first). Use the same
region labels as in §4.6 for consistency.

---

## §5. Consensus Pharmacophore Query

```markdown
## 5. Consensus Pharmacophore Query

Features were clustered across the {{N}} aligned entries using greedy single-linkage clustering (cutoff = 2 × feature tolerance radius). Features present in **≥{{ceil(N/2)}} of {{N}} source PDBs** are classified as **mandatory**; single-source features are **optional**.

### 5.1 Query Summary

|                       | Count |
|-----------------------|-------|
| **Total features**    | {{n_total}}    |
| **Mandatory**         | {{n_mandatory}}    |
| **Optional**          | {{n_optional}}    |

| Family        | Mandatory | Optional |
|---------------|-----------|----------|
| {{Donor}}        | {{...}}        | {{...}}        |
| {{...}}  | {{...}}        | {{...}}        |

### 5.2 Mandatory Features ({{n_mandatory}})

| ID | Family          | Position (x, y, z) Å            | Tol (Å) | Sources |
|----|-----------------|----------------------------------|---------|---------|
| {{id}}  | {{family}}           | [{{x}}, {{y}}, {{z}}]        | {{tol}}     | {{n_sources}} |
| {{...}} | {{...}} | {{...}} | {{...}} | {{...}} |

### 5.3 Optional Features ({{n_optional}})

(Same table as §5.2 for optional features.)

### 5.4 Screening Strategies

| Strategy  | Criteria                                                  | Expected Behavior              |
|-----------|-----------------------------------------------------------|--------------------------------|
| **Strict**  | Match all {{n_mandatory}} mandatory features              | High precision, low recall     |
| **Relaxed** | Match ≥{{ceil(n_mandatory/2)}} mandatory **+** ≥3 optional | Balanced precision/recall      |
| **Scoring** | Match = mandatory_matches + 0.5 × optional_matches; threshold ≥{{0.75 × n_mandatory}} | Tunable via threshold          |

A standalone `consensus_pharmacophore.json` file (in this directory) is provided for direct import into screening tools (Pharmer, align-it, RDKit shape screening).

### 5.5 Feature Retention Analysis

Per-family retention across the {{N}} orthosteric source PDBs ({{pdb_a}}, {{pdb_b}}, {{pdb_c}}), summing the unfiltered (`_all`) and filtered counts:

| Family             | All-ligand total | Filtered total | Retention | Inclusion in query |
|--------------------|:---:|:---:|:---:|---|
| **Donor**          | {{all}} | {{filtered}} | {{XX%}}   | **{{Mandatory / Optional / Exclude}}** — {{1-line rationale tied to inclusion decision}} |
| **PosIonizable**   | {{...}} | {{...}} | {{...}}   | **{{...}}** — {{...}} |
| ... | | | | |

> **Build screening queries from {{primary families}} features only.** {{secondary family}} is conditional. Drop {{excluded families}} — they {{specific reason: appear in some ligands but never engage / have no instances at all}}.

---
```

Coordinate columns in §5.2/§5.3 come from the consensus JSON. Keep three
decimal places; pad with two-digit subscript IDs for sortability if needed.

---

## §6. Validation Summary

```markdown
## 6. Validation Summary

| Criterion                                | Evidence                                                                 |
|------------------------------------------|---------------------------------------------------------------------------|
| **Internal consistency**                 | {{Same-ligand pair: mean feature match Å, # unmatched features}} |
| **Chemotype independence**               | {{Cross-chemotype mean distance}} |
| **Filter specificity**                   | {{Excluded ligands' retention vs orthosteric retention}} |
| **Biological plausibility**              | {{Headline feature dominance matches known receptor biology}} |
| **Reproducibility across conformational states** | {{Pharmacophore preserved across X Å Cα RMSD shift}} |

### Limitations

- {{Number of source structures and chemotypes}}
- {{Any low-confidence feature, e.g. unique Acceptor with single source}}
- {{Reference frame caveat — features are in the aligned reference's frame}}
- {{Clustering threshold caveat — borderline cases}}

---
```

The Validation table is fixed-format. The Limitations are receptor-specific
but follow the same 4-bullet pattern: source diversity, low-confidence
features, frame caveat, clustering threshold.

---

## §7. Data Quality Assessment

```markdown
## 7. Data Quality Assessment

| Issue | Impact | Mitigation |
|-------|--------|------------|
| {{e.g. ligand X has smiles_used: null (template fallback)}} | {{Bond-order assignments uncertain}} | {{Cross-check via PyMOL inspection; supply validated SMILES on re-run}} |
| {{Numbering offset: PDB X uses ALIAS_A while PDBs Y, Z use ALIAS_B (N-residue offset)}} | {{Same residue can appear non-conserved if numbering not normalized}} | {{Caveat in §4.6; map to UniProt canonical numbering in cross-reference}} |
| {{Coordinate-coincident optional features (IDs A/B, C/D — same ligand atoms in dual roles)}} | {{Inflates feature count; double-counts the same atom}} | {{Treat each coordinate-coincident pair as one ligand atom in scoring}} |
| {{Consensus from N entries / M chemotypes}} | {{Limited chemotype coverage}} | {{Future: incorporate additional structures with diverse chemotypes if available}} |
| {{Excluded ligand validation}} | {{None — built-in negative control}} | {{Use as filter: hits matching only excluded-ligand-like positions should be rejected}} |
| {{No cross-validation against external screening hits}} | {{Predictive performance unmeasured}} | {{Run a retrospective benchmark: score the source ligands + a decoy set against the consensus}} |

---
```

This section is mandatory. Pull issues from the `numbering_offsets` field
and from any per-PDB JSON where `smiles_used` is null. The last two rows
(consensus diversity + lack of cross-validation) are nearly always
applicable.

---

## §8. Deployment Guidance

```markdown
## 8. Deployment Guidance

### 8.1 Workflow for Prospective Screening

\`\`\`
Query (consensus_pharmacophore.json)
    │
    ▼
Conformer generation (ETKDG or OMEGA)
    │
    ▼
Feature annotation (RDKit ChemicalFeatures)
    │
    ▼
Flexible alignment to query (Pharmer / align-it / RDKit)
    │
    ├── Strict mode: require {{n_mandatory}}/{{n_mandatory}} mandatory features matched
    ├── Relaxed mode:     require ≥{{ceil(n_mandatory/2)}} mandatory + ≥3 optional
    └── Scoring mode:     score ≥{{0.75 × n_mandatory}} threshold
    │
    ▼
Post-filter: docking (optional), PAINS, Lipinski
\`\`\`

### 8.2 Tool Compatibility

| Tool              | Input Format                     | Notes                                |
|-------------------|----------------------------------|--------------------------------------|
| **Pharmer**       | `.json` pharmacophore file       | Direct import; supports partial matching |
| **align-it**      | Feature list + SMILES            | Flexible alignment scoring            |
| **RDKit**         | `ChemicalFeatures` + custom code | Full control; sample script available |
| **MOE**           | Pharmacophore editor             | Manual entry of feature positions    |
| **LigandScout**   | `.pml` file                      | Visual inspection; manual query build |

---
```

Fixed structure. Only the threshold numbers in §8.1 are receptor-specific.

---

## §9. Appendix: Excluded Ligand Details

```markdown
## 9. Appendix: Excluded Ligand Details

### 9.1 {{First Excluded Ligand}} — Full Results

| PDB ID (chain) | Filtered Features | Types Retained                    |
|----------------|------------------|-----------------------------------|
| {{pdb}}           | {{n}}                | {{family list or —}}                                 |
| {{...}} | {{...}} | {{...}} |

{{One paragraph closing argument: total raw features, why they consistently fail the interaction filter, what binding site they occupy.}}

### 9.2 {{Next Excluded Ligand}}

(Same format.)

---
```

One sub-subsection per excluded ligand. Sort by # copies (most copies first
— typically cholesterol).

---

## §10. Output File Manifest

```markdown
## 10. Output File Manifest

\`\`\`
pharmacophore/
├── REPORT.md                                            ◄── this report
├── consensus_pharmacophore.json                         ◄── {{n_total}}-feature VS query ({{n_mandatory}} mandatory, {{n_optional}} optional)
│
│   Orthosteric source PDBs — used to derive the consensus
├── {{pdb_a}}/  ({{ligand}}, chain {{X}}, {{n}} copy)
│   ├── {{pdb_a}}{_chain}_all.pharmacophore.json               {{n_all}} features (raw)
│   ├── {{pdb_a}}{_chain}_all.pharmacophore.pml
│   ├── {{pdb_a}}{_chain}_filtered.pharmacophore.json          {{n_filtered}} features ({{1-line role description: anchored on RESA + RESB}})
│   ├── {{pdb_a}}{_chain}_filtered.pharmacophore.pml
│   └── {{pdb_a}}_clean.pdb                                   {{symlink → ../../raw/{{pdb_a}}/{{pdb_a}}_clean.pdb}}
├── {{pdb_b}}/                                                  {{n}} filtered ({{role}})
│
│   Excluded ligands — interaction filter dropped them (see §3.2, §9)
├── {{pdb_c}}/                  {{LIG}}    {{n_chains}}   {{0–N filtered each}}
└── ...
\`\`\`

Each `*_all.pharmacophore.json` contains the unfiltered feature set; each `*_filtered.pharmacophore.json` retains only features whose family has a complementary protein-side feature within the cutoffs in §2. Both formats are paired with a `.pml` script that loads the cleaned receptor + ligand and renders the features as colored pseudo-atoms (HBD blue, HBA red, PosIonizable marine, LumpedHydrophobe forest).

To visualize a per-PDB feature set:
\`\`\`bash
cd pharmacophore/{{pdb_a}}
pymol {{pdb_a}}{_chain}_filtered.pharmacophore.pml
\`\`\`

To regenerate the raw data from scratch (PDB → cleaned monomers + ligands):
\`\`\`bash
python scripts/download_and_extract.py
\`\`\`
```

Sort orthosteric PDBs first (with brief role descriptions), then excluded
PDBs grouped by ligand (with one summary line). Skip PML triples for
zero-feature entries to keep the manifest readable.

---

## Universal style conventions

- **Pull every number from `.report_data.json`.** Don't re-derive counts; the script already did it. If a number isn't in the JSON, that's a bug — flag it and use a `{{compute manually}}` placeholder rather than guessing.
- **Bold the conserved/primary anchors** in tables and prose (e.g., **ASP270**, **GLU296**). Leave secondary residues unbold.
- **Cross-reference file paths** rather than re-quoting JSON. "See `pharmacophore/8U4P/8U4P_chainR_filtered.pharmacophore.json`" beats embedding the JSON in the report.
- **Don't hedge in the executive summary or §6 validation.** Commit to specific findings. Hedging belongs in §7.
- **Keep tables narrow.** If a column doesn't pull weight, drop it. Markdown tables wider than ~110 chars wrap badly in monospace viewers.
- **Numbering aliases** (from `numbering_offsets`) appear in three places: §4.6 caveat, §4.7 ASCII map, and §7 data quality. Use the same canonical form in all three.
- **Emoji and decoration: none.** This is a scientific deliverable.

## Receptor-class context to add

The template is receptor-agnostic. Layer in the right vocabulary based on
what the receptor is (ask the user if unclear — don't guess):

| Class | Pocket descriptors | Ligand categorization |
|---|---|---|
| **GPCR** | TM1–TM7; ECL1/ECL2/ECL3; ICL1/ICL2/ICL3; Ballesteros-Weinstein numbering (Asp^6.58 = pos 58 from most conserved residue in TM6) | Orthosteric / allosteric; agonist / inverse agonist; biased / unbiased |
| **Kinase** | Hinge region; DFG motif (in/out); αC helix (in/out); P-loop; gatekeeper | Type I (ATP-competitive, DFG-in) / Type II (DFG-out) / Type III (allosteric) |
| **Nuclear receptor** | H1–H12; AF-2 surface; LBP (ligand-binding pocket); coactivator-binding groove | Agonist / antagonist; SERM / SARM / etc. |
| **Protease** | Catalytic triad (Ser/His/Asp or Cys/His); oxyanion hole; S1/S2/S3 specificity pockets | Substrate-mimetic / covalent / allosteric |
| **Other enzyme** | Active site residues; substrate channel; allosteric site | Substrate analog / transition-state mimetic / allosteric |

If you're not sure what class the receptor is, **ask the user** before
writing — getting the secondary-structure labels wrong undermines the entire
report.
