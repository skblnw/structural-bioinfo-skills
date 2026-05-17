# RCSB Search API queries issued by this skill

Verbatim JSON payloads for debugging. To capture the exact query without running the API,
pass `--query-only` to `search_pmhc_tcr.py`.

Endpoint: `POST https://search.rcsb.org/rcsbsearch/v2/query`

## Why `seqmotif` and not `sequence` or `text`

For short peptide search we evaluated three RCSB Search services and only one is fit:

| Service | Works for 8–25 aa peptides? | Notes |
|---|---|---|
| `sequence` | ❌ | Backed by BLAST, requires ≥25 aa input — too long for MHC-I/II epitopes. |
| `text` on `entity_poly.pdbx_seq_one_letter_code_can` | ❌ | RCSB does NOT index this attribute for the text service. Returns HTTP 400 *"search is not enabled on [...] attribute"*. |
| `seqmotif` (`pattern_type: "simple"`) | ✅ | Substring-matches any polymer entity sequence with no minimum length. The right tool. |

We initially tried AND-combining a sequence-substring terminal with a UniProt-IN terminal
to require an MHC chain on the same entry — but
`rcsb_polymer_entity_container_identifiers.uniprot_ids` is *also* not text-searchable.
RCSB exposes UniProt-keyed search only via cluster/group services, not via plain text
terminals.

The skill therefore issues a single `seqmotif` query per epitope (no MHC pre-filter) and
the **classifier** (in `search_pmhc_tcr.py`) drops hits that lack an MHC chain after
fetching their metadata. For 8+ aa peptides the over-fetch is tiny (most peptide
substrings of 8+ aa hit fewer than 50 entries even before pre-filtering).

## The query

```json
{
  "return_type": "entry",
  "query": {
    "type": "terminal",
    "service": "seqmotif",
    "parameters": {
      "value": "SIINFEKL",
      "pattern_type": "simple",
      "sequence_type": "protein"
    }
  },
  "request_options": {
    "results_content_type": ["experimental"],
    "return_all_hits": true,
    "sort": [
      { "sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc" }
    ]
  }
}
```

`pattern_type: "simple"` does literal substring matching against every protein polymer
entity in the experimental PDB. `"prosite"` and `"regex"` are also supported if you need
to match an epitope with wildcards (e.g. anchor-residue patterns like `Y..N.K.[LV]`).

## Per-entry metadata — GraphQL Data API

Endpoint: `POST https://data.rcsb.org/graphql`

After the search returns a list of PDB IDs, one batched query collects everything needed
for classification + splitting:

```graphql
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    rcsb_entry_info {
      resolution_combined
      experimental_method
    }
    rcsb_accession_info { deposit_date }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        entity_id
        uniprot_ids
        auth_asym_ids
      }
      entity_poly {
        pdbx_seq_one_letter_code_can
        rcsb_sample_sequence_length
        type
      }
      rcsb_polymer_entity { pdbx_description }
    }
  }
}
```

Batch size is 50 PDB IDs per request (matches the sibling `pdb-holostructure-search`
convention; the Data API tolerates much larger batches, but 50 keeps individual responses
under a megabyte and survives intermittent network blips).

## How `match_type` is decided

After metadata fetch, each classified entry is tagged with one of:

- `free` — the chain that hosts the epitope has a canonical sequence *equal* to the
  epitope. The peptide sits in its own polymer entity. Cleanest for downstream docking
  since the peptide is already isolated.
- `fused` — the chain that hosts the epitope is longer than the epitope (the peptide is
  embedded in a single-chain trimer or fusion construct). The `_pmhc.pdb` derivative
  carries the linker as part of the MHC chain — the peptide is not surgically excised.
- `?` — the classifier saw a peptide-hosting chain but neither equality nor substring
  matched cleanly (rare; usually a non-standard residue collapse).

## Notes on RCSB attribute names

- `entity_poly.pdbx_seq_one_letter_code_can` — canonical sequence with non-standard
  residues mapped to their canonical single-letter codes; the right attribute for
  *post-hoc* epitope substring inspection (we read it from the GraphQL Data API, not the
  Search API).
- `auth_asym_ids` are the chain IDs as they appear in the PDB file column 22 — what the
  splitter needs. `asym_ids` (label) are the mmCIF-internal IDs and differ for entries
  with multiple copies of the same entity.
- `experimental_method` is a string ("X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"); the
  attribute `exptl.method` is the same field accessed differently.
