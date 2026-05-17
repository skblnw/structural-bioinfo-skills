# TCR chain detection heuristics

TCR α and β chains in the PDB are inconsistently annotated:

- Sometimes the polymer entity's `rcsb_polymer_entity_container_identifiers.uniprot_ids` is
  populated with a TRAC/TRBC reference (e.g. P01848, P01850), but only when the structure
  was deposited with that mapping.
- Frequently the UniProt list is empty (engineered constructs, soluble TCRs with
  non-native constant domains), and the only clue is the free-text
  `rcsb_polymer_entity.pdbx_description`.

The skill therefore uses **UniProt accession match OR description regex match** to flag a
polymer entity as `tcr_alpha` or `tcr_beta`.

## UniProt accessions

| Chain | Gene | UniProt | Notes |
|---|---|---|---|
| TCR α constant | TRAC | P01848 | human |
| TCR β constant | TRBC1 | P01850 | human |
| TCR β constant | TRBC2 | A0A075B6X1 | human; less common |
| TCR α constant | Trac | P01849 | mouse |
| TCR β constant | Trbc1 | P01852 | mouse |

V-segment UniProts (TRAV / TRBV) are too numerous and many engineered TCRs use
non-germline constructs that lack a UniProt mapping; the description regex below is the
more robust signal.

## Description regexes (case-insensitive)

A polymer entity is classified as `tcr_alpha` if its `pdbx_description` matches any of:

- `T[- ]cell receptor.*\balpha\b`
- `\bTCR\b.*\balpha\b`
- `\bTCR[\-\s]?α\b`
- `\bTRA[VJC]\b`        (e.g. "TRAV12-2", "TRAC")
- `\balpha[- ]chain\b`  (when co-occurring with "T cell receptor" in title; gated)

And `tcr_beta` for:

- `T[- ]cell receptor.*\bbeta\b`
- `\bTCR\b.*\bbeta\b`
- `\bTCR[\-\s]?β\b`
- `\bTRB[VDJC]\b`
- `\bbeta[- ]chain\b`   (gated as above)

The `alpha-chain` / `beta-chain` patterns are only applied when the description ALSO
contains "T cell receptor", "TCR", or matches the V/J/C-segment patterns — to avoid
misclassifying MHC Class II β-chains (which use "beta chain" in their descriptions) or
hemoglobin α/β as TCR.

## Order of evidence

1. UniProt match (highest precision).
2. Description regex.
3. If neither hits but the entity is the "remaining" short chain in a 5-chain pMHC+TCR
   assembly with two MHC chains and one peptide, the script will *not* guess — instead it
   logs a warning and skips TCR splitting for that entry.

## γδ TCRs

This skill targets αβ TCRs (the overwhelming majority of TCR-pMHC depositions). γδ TCRs
recognise MHC-I-like molecules (CD1, MR1, MICA) and use different chain conventions; they
are *not* classified as TCR-pMHC here and will land in the "unclassified" bucket of the
report. Future work could add `\bTRG[VC]\b` / `\bTRD[VC]\b` and CD1/MR1 UniProts.
