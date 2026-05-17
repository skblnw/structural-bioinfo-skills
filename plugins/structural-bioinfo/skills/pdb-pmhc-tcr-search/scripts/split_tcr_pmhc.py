#!/usr/bin/env python3
"""Split a cached PDB into pMHC-only and TCR-only derivative files.

Given a PDB file (already downloaded by `search_pmhc_tcr.py`) plus the per-component
author chain IDs determined by the GraphQL classifier, write:

  - `<pdb_id>_pmhc.pdb`   : peptide + MHC heavy + β2m (or class-II β) chains
  - `<pdb_id>_tcr.pdb`    : TCR α + TCR β chains (only emitted for tcr_pmhc entries)

Standard junk (water, ions, buffers, cryoprotectants) is dropped via the EXCLUDE_RESIDUES
constant. ATOM records on the kept chains pass through verbatim; HETATM records pass
through only if their residue is not in the exclusion list AND the chain is one of the
kept ones.

Importable as a function or runnable from the CLI for manual splits.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

# Mirror pdb-extractor's exclusion list. Copied (not imported) to keep this skill
# stdlib-only and self-contained; if the upstream list changes, refresh here.
EXCLUDE_RESIDUES: set[str] = {
    "HOH", "DOD", "WAT",
    "NA", "K", "CL", "MG", "CA", "ZN", "FE", "MN", "CU", "CO",
    "CD", "NI", "BR", "I", "F", "LI", "RB", "CS", "SR", "BA",
    "PT", "AU", "HG", "SE",
    "EDO", "GOL", "TRS", "PEG", "PGE", "MPD", "ACT", "CIT",
    "BME", "DMS", "SO4", "PO4", "FMT", "EOH", "IMD",
    "BCT", "BOG", "DMU", "OLC", "PLM", "SCN", "NO3", "NH4",
    "TAM", "HEP", "MES", "PIP", "BIS", "TRI", "BTB",
    "GLY", "SUC", "DTT", "TCE", "DTE", "SGM",
    "UNK", "UNX", "ACE", "NH2", "FOR", "DUM",
    "LDA", "LMT", "OLA", "SDS", "LMU", "LMN", "TGL", "D10",
    "DMX", "OCT", "HEX", "LHG",
}


def _filter_lines(pdb_path: str, keep_chains: set[str]) -> list[str]:
    """Stream the input PDB and keep ATOM/HETATM lines whose chain ID is in keep_chains."""
    if not keep_chains:
        return []
    out: list[str] = []
    last_chain: str | None = None
    with open(pdb_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = line[:6].rstrip()
            if rec in ("ATOM", "HETATM"):
                if len(line) < 22:
                    continue
                chain = line[21]
                if chain not in keep_chains:
                    continue
                if rec == "HETATM":
                    res = line[17:20].strip().upper()
                    if res in EXCLUDE_RESIDUES:
                        continue
                # TER insertion between chains for cleaner output
                if last_chain is not None and chain != last_chain:
                    out.append("TER\n")
                out.append(line if line.endswith("\n") else line + "\n")
                last_chain = chain
            elif rec == "MODEL" or rec == "ENDMDL":
                # Preserve model boundaries verbatim
                out.append(line if line.endswith("\n") else line + "\n")
            elif rec == "END":
                break
    if out and last_chain is not None:
        out.append("TER\n")
    out.append("END\n")
    return out


def split_entry(
    pdb_path: str,
    output_dir: str,
    pdb_id: str,
    peptide_chains: list[str],
    mhc_chains: list[str],
    light_chains: list[str],          # β2m (Class I) or class-II β (Class II)
    tcr_alpha_chains: list[str],
    tcr_beta_chains: list[str],
) -> dict:
    """Write derivative PDB files. Returns a dict of relative output paths."""

    pmhc_chains = set(peptide_chains) | set(mhc_chains) | set(light_chains)
    tcr_chains = set(tcr_alpha_chains) | set(tcr_beta_chains)

    os.makedirs(output_dir, exist_ok=True)
    result: dict = {"pmhc": None, "tcr": None, "pmhc_chains": sorted(pmhc_chains),
                    "tcr_chains": sorted(tcr_chains)}

    if pmhc_chains:
        pmhc_path = os.path.join(output_dir, f"{pdb_id}_pmhc.pdb")
        lines = _filter_lines(pdb_path, pmhc_chains)
        if any(ln.startswith(("ATOM", "HETATM")) for ln in lines):
            with open(pmhc_path, "w", encoding="utf-8") as fh:
                fh.writelines(lines)
            result["pmhc"] = pmhc_path

    if tcr_chains:
        tcr_path = os.path.join(output_dir, f"{pdb_id}_tcr.pdb")
        lines = _filter_lines(pdb_path, tcr_chains)
        if any(ln.startswith(("ATOM", "HETATM")) for ln in lines):
            with open(tcr_path, "w", encoding="utf-8") as fh:
                fh.writelines(lines)
            result["tcr"] = tcr_path

    return result


def _parse_chains(arg: str) -> list[str]:
    if not arg:
        return []
    return [c.strip() for c in arg.replace(";", ",").split(",") if c.strip()]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb-path", required=True, help="Path to the cached full-entry PDB file.")
    p.add_argument("--pdb-id", required=True, help="PDB ID, used in output filenames.")
    p.add_argument("--output-dir", required=True, help="Directory for the _pmhc.pdb / _tcr.pdb files.")
    p.add_argument("--peptide-chains", default="", help="Comma-separated peptide auth chain IDs.")
    p.add_argument("--mhc-chains", default="", help="Comma-separated MHC heavy auth chain IDs.")
    p.add_argument("--light-chains", default="",
                   help="Comma-separated β2m (Class I) or class-II β auth chain IDs.")
    p.add_argument("--tcr-alpha-chains", default="", help="Comma-separated TCR α auth chain IDs.")
    p.add_argument("--tcr-beta-chains", default="", help="Comma-separated TCR β auth chain IDs.")
    p.add_argument("--json", action="store_true", help="Emit result dict as JSON to stdout.")
    args = p.parse_args()

    if not os.path.isfile(args.pdb_path):
        sys.stderr.write(f"[error] --pdb-path not found: {args.pdb_path}\n")
        return 2

    result = split_entry(
        pdb_path=args.pdb_path,
        output_dir=args.output_dir,
        pdb_id=args.pdb_id,
        peptide_chains=_parse_chains(args.peptide_chains),
        mhc_chains=_parse_chains(args.mhc_chains),
        light_chains=_parse_chains(args.light_chains),
        tcr_alpha_chains=_parse_chains(args.tcr_alpha_chains),
        tcr_beta_chains=_parse_chains(args.tcr_beta_chains),
    )

    if args.json:
        json.dump(result, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        if result["pmhc"]:
            sys.stderr.write(f"[info] pMHC → {result['pmhc']}  chains={result['pmhc_chains']}\n")
        if result["tcr"]:
            sys.stderr.write(f"[info] TCR  → {result['tcr']}  chains={result['tcr_chains']}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
