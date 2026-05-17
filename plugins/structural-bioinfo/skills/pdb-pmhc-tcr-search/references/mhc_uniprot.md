# MHC and accessory chain UniProt accessions

These accessions seed the RCSB Search API filter that asks "does this PDB entry contain at
least one MHC polymer chain?". They are loaded verbatim by `search_pmhc_tcr.py` (see the
`MHC_UNIPROTS` constant). Add or remove entries here if a relevant organism is missing.

## MHC Class I — heavy chain

| Organism | Gene / allele family | UniProt | Notes |
|---|---|---|---|
| Human | HLA-A | P04439 | reference HLA-A*03:01 entry; other -A alleles cross-list to this |
| Human | HLA-B | P01889 | reference HLA-B*07:02 |
| Human | HLA-C | P10321 | reference HLA-C*04:01 |
| Human | HLA-E | P13747 | non-classical Class I |
| Human | HLA-F | P30511 | non-classical |
| Human | HLA-G | P17693 | non-classical |
| Mouse | H-2K | P01901 | H2-K1 |
| Mouse | H-2D | P01899 | H2-D1 |
| Mouse | H-2L | P14430 | H2-L  |
| Mouse | H-2K^b | P01901 | same UniProt as H-2K; allele encoded in description |

## MHC Class I — light chain (β2-microglobulin)

| Organism | UniProt | Notes |
|---|---|---|
| Human | P61769 | β2m |
| Mouse | P01887 | β2m |

## MHC Class II — α chains

| Organism | Locus | UniProt |
|---|---|---|
| Human | HLA-DRA | P01903 |
| Human | HLA-DPA1 | P20036 |
| Human | HLA-DQA1 | P01909 |
| Mouse | H2-Aa | P14434 |
| Mouse | H2-Ea | P14437 |

## MHC Class II — β chains

| Organism | Locus | UniProt |
|---|---|---|
| Human | HLA-DRB1 | P04229 |
| Human | HLA-DPB1 | P04440 |
| Human | HLA-DQB1 | P01920 |
| Mouse | H2-Ab1 | P14483 |
| Mouse | H2-Eb1 | P04230 |

## Sets used by the script

The script groups these into three sets, exposed at module top:

- `MHC_CLASS_I_UNIPROTS`  : Class I heavy chains (human + mouse)
- `MHC_CLASS_II_ALPHA_UNIPROTS`, `MHC_CLASS_II_BETA_UNIPROTS`
- `B2M_UNIPROTS`          : light chain
- `MHC_UNIPROTS`          : union of the four above — used in the search filter

Many PDB entries cross-list to multiple HLA UniProts (alleles share UniProt records); the
classifier matches on the *set* of UniProts on each polymer entity rather than a single ID.
