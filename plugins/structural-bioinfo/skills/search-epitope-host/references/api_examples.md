# API patterns beyond the bundled pipeline

The bundled `search_epitope_host.py` covers exact-match epitope → UniProt parent → structure. Use the patterns below when the user asks for something outside that scope. All examples are runnable with `curl`; the skill stays stdlib-only.

## IEDB IQ-API

Base URL: `https://query-api.iedb.org/`. PostgREST-style filters.

### Exact peptide match (what the skill uses)

```
GET https://query-api.iedb.org/epitope_search?linear_sequence=eq.SIINFEKL
```

### Substring / similar peptide search

PostgREST supports `like` and `ilike`:

```
GET https://query-api.iedb.org/epitope_search?linear_sequence=ilike.*SIINFEK*
```

Use this when the user has an unconfirmed peptide or wants nested epitopes. Cap with `&limit=200` — broad patterns can return thousands.

### Restrict by source organism

```
GET https://query-api.iedb.org/epitope_search?linear_sequence=eq.GILGFVFTL&source_organism_names=cs.{"Influenza A virus"}
```

`cs` is PostgREST's "contains" operator for arrays.

### Selecting only the fields you want

```
&select=structure_id,linear_sequence,curated_source_antigens,source_organism_names,host_organism_names
```

Default response is heavy (T-cell receptor sequences, MHC allele arrays, references, etc.). Always pick a `select` list when scripting.

### Epitope → MHC / T cell / B cell

For receptor or MHC context, follow up with `/epitope_to_tcell`, `/epitope_to_bcell`, `/epitope_to_mhc` using the `structure_id` from the first query.

### Quirks

- The `curated_source_antigens` array can have **mixed source DBs** in one row: `UNIPROT:`, `GENPEPT:`, `NCBI:` prefixes. Filter to what you care about.
- Versioned UniProt accessions (`P01012.2`) — strip the `.N` suffix before calling UniProt.
- Some older entries report a `starting_position` that doesn't match the **current** canonical UniProt sequence (sequence revisions). The bundled script catches this by also doing a substring scan.

## UniProt REST

Base URL: `https://rest.uniprot.org/`.

### Entry as JSON / FASTA / TXT

```
GET https://rest.uniprot.org/uniprotkb/P01012.json
GET https://rest.uniprot.org/uniprotkb/P01012.fasta
GET https://rest.uniprot.org/uniprotkb/P01012.txt
```

### Useful JSON paths

| What | Path |
|---|---|
| Canonical sequence | `sequence.value` |
| Sequence length | `sequence.length` |
| Recommended name | `proteinDescription.recommendedName.fullName.value` |
| Gene symbol | `genes[0].geneName.value` |
| Organism (scientific) | `organism.scientificName` |
| PDB cross-refs | `uniProtKBCrossReferences[?database=="PDB"]` — each has `id`, `properties: [{key:"Method", value:"X-ray"}, {key:"Resolution", value:"1.95 A"}, {key:"Chains", value:"A/B/C/D=2-386"}]` |
| AlphaFold cross-ref | `uniProtKBCrossReferences[?database=="AlphaFoldDB"]` — single entry per accession |

### ID mapping (when the user has a non-UniProt accession)

```
POST https://rest.uniprot.org/idmapping/run -d "from=GeneID&to=UniProtKB&ids=1956"
GET  https://rest.uniprot.org/idmapping/status/{jobId}
GET  https://rest.uniprot.org/idmapping/results/{jobId}
```

Two-step async — poll status until `results` is available. Use when an IEDB row gives a GenPept accession and the user really wants the UniProt-anchored workflow.

### Search by peptide (BLAST-like)

UniProt's REST `search` endpoint accepts a `sequence:` term that does a peptide substring scan in SwissProt (slow, full-text). For real similarity searches, prefer the BLAST endpoint at `https://www.ebi.ac.uk/Tools/sss/ncbiblast/rest/`.

## AlphaFold Database

Base URL: `https://alphafold.ebi.ac.uk/`.

### Prediction metadata

```
GET https://alphafold.ebi.ac.uk/api/prediction/P01012
```

Returns an array; usually one entry for the canonical accession, sometimes additional entries for isoforms (`P00533-2`, `P00533-3`, ...). Each entry includes:

| Field | What |
|---|---|
| `entryId` | e.g., `AF-P01012-F1` (fragment 1 — typically the full sequence; large proteins are split across multiple F-N entries) |
| `uniprotAccession`, `uniprotId` | canonical & mnemonic |
| `uniprotSequence` | full sequence the model was built against |
| `uniprotStart`, `uniprotEnd` | residue range the entry covers (relevant for fragmented predictions) |
| `pdbUrl`, `cifUrl`, `bcifUrl` | downloadable structure files |
| `paeImageUrl`, `paeDocUrl` | predicted-aligned-error products |
| `latestVersion`, `allVersions` | model version history |
| `modelCreatedDate` | ISO date string |

### Fragmented predictions

Proteins over ~2,700 aa are split into overlapping fragments. The API returns multiple `F-N` entries (`AF-Q8WZ42-F1`, `AF-Q8WZ42-F2`, …). When showing structures, list each fragment's range so the user can pick the right one.

### Empty response

For accessions without an AlphaFold model (some TrEMBL entries, viral proteins below the size threshold), the API returns `[]` (not 404). The bundled script treats that as "no AlphaFold model".

## Predicted-structure alternatives

If the user wants more than AlphaFold:

- **ESMFold Atlas** (Meta) — older metagenomic predictions, accessible via `https://api.esmatlas.com/`. Useful for proteins unmodeled by AlphaFold.
- **BFVD** (viral fold database) — `https://bfvd.foldseek.com/` — community-curated viral predictions, often complementary to AlphaFold for short viral antigens.

Neither is wired into the bundled script; trigger them manually only when explicitly requested.

## RCSB structure data (deeper PDB queries)

The UniProt `uniProtKBCrossReferences` block gives you PDB IDs but no ligands, no biological assemblies, no quality scores. For deeper structure mining (apo vs. holo, ligand filtering, resolution windows), pivot to the sister skill `pdb-holostructure-search` — point the user there rather than duplicating its logic.
