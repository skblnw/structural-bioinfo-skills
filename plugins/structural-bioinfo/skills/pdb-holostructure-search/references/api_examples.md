# RCSB API patterns beyond the base case

The bundled `search_holostructures.py` covers the standard "find protein + meaningful ligand" workflow. For niche needs, use these query fragments. They can be inserted into the `nodes` array of the search query the script builds (see `--query-only` to dump the current query as a starting point).

## Endpoints

- Search API: `POST https://search.rcsb.org/rcsbsearch/v2/query`
- Data API (GraphQL): `POST https://data.rcsb.org/graphql`
- Data API (REST): `https://data.rcsb.org/rest/v1/core/{entity}/{pdb_id}/{entity_id}`
- Chemical Component Dictionary: `https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}`

## Search-by-UniProt (no sequence search)

When the user wants only entries explicitly annotated to a single UniProt accession, skip sequence search:

```json
{
  "type": "group", "logical_operator": "and", "label": "nested-attribute",
  "nodes": [
    {"type":"terminal","service":"text","parameters":{
      "attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
      "operator":"exact_match","value":"P00533"}},
    {"type":"terminal","service":"text","parameters":{
      "attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
      "operator":"exact_match","value":"UniProt"}}
  ]
}
```

The `label: "nested-attribute"` matters: without it, accession and database name can match different sub-objects.

## Specific ligand by chem-comp ID

```json
{"type":"terminal","service":"text","parameters":{
  "attribute":"rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
  "operator":"in","value":["ATP","ADP","AMP"]}}
```

## Free vs covalently bound ligand

```json
{"type":"terminal","service":"text","parameters":{
  "attribute":"rcsb_nonpolymer_instance_annotation.type",
  "operator":"exact_match","value":"HAS_NO_COVALENT_LINKAGE"}}
```

Swap to `HAS_COVALENT_LINKAGE` for covalent-inhibitor co-crystals.

## Binding-affinity filter (must be nested)

```json
{
  "type": "group", "logical_operator": "and", "label": "nested-attribute",
  "nodes": [
    {"type":"terminal","service":"text","parameters":{
      "attribute":"rcsb_binding_affinity.type","operator":"exact_match","value":"Kd"}},
    {"type":"terminal","service":"text","parameters":{
      "attribute":"rcsb_binding_affinity.value","operator":"less","value":1000}}
  ]
}
```

Units: nM for IC50/EC50/Ki/Kd; kJ/mol for ΔG/ΔH/-T·ΔS; M⁻¹ for Ka.

## Ligand quality (RSCC, RSR)

Per-instance validation scores live in `rcsb_nonpolymer_instance_validation_score`. Wrap related sub-fields in a `nested-attribute` group:

```json
{
  "type": "group", "logical_operator": "and", "label": "nested-attribute",
  "nodes": [
    {"type":"terminal","service":"text","parameters":{
      "attribute":"rcsb_nonpolymer_instance_validation_score.RSCC",
      "operator":"greater_or_equal","value":0.8}},
    {"type":"terminal","service":"text","parameters":{
      "attribute":"rcsb_nonpolymer_instance_validation_score.RSR",
      "operator":"less_or_equal","value":0.2}}
  ]
}
```

## Group / deduplicate by sequence cluster

```json
"request_options": {
  "group_by": {"aggregation_method":"sequence_identity","similarity_cutoff":90},
  "group_by_return_type":"representatives"
}
```

`similarity_cutoff` accepts 30 / 50 / 70 / 90 / 95 / 100. Useful for "give me one representative per ortholog cluster".

## Counts only

```json
"request_options": {"return_counts": true}
```

## Pulling ligand SMILES / SDF

The Search API never returns chemical-component-level data. Use the Data API:

```graphql
{ chem_comps(comp_ids: ["ATP","STI"]) {
    chem_comp { id name formula formula_weight type }
    rcsb_chem_comp_descriptor { SMILES InChI InChIKey }
} }
```

Or fetch SDF from the file service:

```
https://files.rcsb.org/ligands/download/<COMP_ID>.sdf
https://files.rcsb.org/ligands/view/<COMP_ID>.cif
```

## Gotchas

- **Identity cutoff scale.** Newer docs use 0.0–1.0; some older examples use 0–100. The current API expects 0.0–1.0.
- **`SUBJECT_OF_INVESTIGATION` is sparse on older entries.** It rolled out around 2020. Pre-2018 entries may not carry it; fall back to MW + exclusion list (`--strict-mw` in the bundled script).
- **MW filters live in the chemical schema.** Use `service: "text_chem"` and `chem_comp.formula_weight`, not the `text` service.
- **`nested-attribute` label is mandatory** for any group whose terminals all reference sub-fields of the same nested object (binding affinity, instance validation, reference sequence ids). Skipping it produces silent cross-product matches.
- **Pagination cap is 10,000 per request.** For larger result sets, use `return_all_hits: true` (used by default in the bundled script).
- **Sort key is `sort_by`** in the v2 API. `sort_field`/`sort_key` may still work but are not canonical.
