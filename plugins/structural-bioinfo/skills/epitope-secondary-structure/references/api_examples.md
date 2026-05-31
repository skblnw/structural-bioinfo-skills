# API & method notes for `epitope-secondary-structure`

How the skill turns a UniProt residue range into a per-residue 8-state DSSP assignment across a
parent's experimentally solved structures. All endpoints are public, no-auth.

## 1. Why not just read secondary structure from RCSB?

The RCSB Data API exposes secondary structure as **coarse instance features** only:

```
GET https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{ENTRY}/{ASYM_ID}
  -> rcsb_polymer_instance_feature[].type in {HELIX_P, SHEET, UNASSIGNED_SEC_STRUCT}
     (provenance PROMOTIF; ranges, not per-residue 8-state)
```

That is helix / sheet / other — it cannot reconstruct the 8 DSSP states (H, G, I, E, B, T, S,
coil). So the 8-state assignment is computed locally with `mkdssp`. RCSB is used for the
**residue mapping**, below.

## 2. UniProt → structure residue mapping (SIFTS, via RCSB)

```
GET https://data.rcsb.org/rest/v1/core/entry/{ENTRY}
  -> rcsb_entry_container_identifiers.polymer_entity_ids   # e.g. ["1"]

GET https://data.rcsb.org/rest/v1/core/polymer_entity/{ENTRY}/{ENTITY_ID}
  -> rcsb_polymer_entity_container_identifiers.asym_ids     # label chains, e.g. ["A","B","C","D"]
  -> rcsb_polymer_entity_align[] (provenance SIFTS):
       reference_database_name == "UniProt"
       reference_database_accession == <parent acc>
       aligned_regions[]: {ref_beg_seq_id (UniProt), entity_beg_seq_id (label_seq_id), length}
```

Mapping formula, for a UniProt residue `u`:

```
label_seq_id = u - ref_beg_seq_id + entity_beg_seq_id        (when ref_beg_seq_id <= u < ref_beg_seq_id + length)
```

Worked example (2CH8, BARF1 / P03228): one entity, `asym_ids=[A,B,C,D]`,
`aligned_regions=[{ref_beg_seq_id:21, entity_beg_seq_id:1, length:201}]`. So UniProt 21–221 maps to
label_seq_id 1–201, and the epitope `AFLGERVTL` at UniProt 23–31 maps to label_seq_id 3–11.

## 3. label_seq_id → author numbering (mmCIF `_atom_site`)

DSSP reports the **author** chain + residue number, but SIFTS gives **label** numbering, and the
two differ per chain (e.g. 2CH8 chain A: label 1 = auth 21). Download the coordinate file and read
`_atom_site` — **header-driven**, because column order differs between RCSB and AlphaFold mmCIF:

```
GET https://files.rcsb.org/download/{ENTRY}.cif
```

Per CA atom of model 1 (first altloc), build:
`(label_asym_id, label_seq_id) -> (auth_asym_id, auth_seq_id, ins_code, one_letter, B_factor)`.
The one-letter code cross-checks the mapping against the epitope residue; the B-factor is pLDDT for
AlphaFold models.

## 4. Running mkdssp and parsing the output

`conda install -c conda-forge dssp` provides `mkdssp` (commonly **3.1.2**). It reads PDB or mmCIF:

```
mkdssp -i {ENTRY}.cif -o {ENTRY}.dssp
```

The classic DSSP residue table starts after the `  #  RESIDUE AA STRUCTURE …` header. Fixed
columns (0-based) used here:

| Field | Columns | Notes |
|---|---|---|
| author residue number | `line[5:10]` | integer |
| insertion code | `line[10]` | blank if none |
| author chain id | `line[11]` | |
| amino acid | `line[13]` | `!` marks a chain break (skip) |
| **secondary structure** | `line[16]` | one of `H G I E B T S` or blank → coil `C` |

Then: UniProt `u` → label_seq_id (§2) → `(auth_asym, auth_seq, icode)` (§3) → DSSP `line[16]` (§4).
(This build of mkdssp also appends `CHAIN AUTHCHAIN NUMBER RESNUM` columns giving label+auth side by
side, but the skill relies on the version-stable classic columns + the `_atom_site` map instead.)

**DSSP 3.1.2 cannot parse AlphaFold mmCIF** ("empty protein, or no valid complete residues"); use
DSSP 4 for AlphaFold (and `--use-alphafold`).

## 5. The 8 → 3 collapse

Standard reduction used for the simplified representation:

```
H, G, I  -> H   (helix)
E, B     -> E   (strand)
T, S, C  -> C   (coil; DSSP blank / '-' / DSSP4 'P' all fold here)
```

The 3-state consensus letter is the argmax of the **summed 3-state probability**, which can differ
from collapsing the 8-state argmax when the top 8-state has < 50% mass.
