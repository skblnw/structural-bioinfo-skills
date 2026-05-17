#!/usr/bin/env python3
"""
Download PDB/CIF structures from RCSB, resolve protein chains via GraphQL
(UniProt ID), and extract a clean PDB file containing the target monomer
and bound ligands. Supports single PDB ID and CSV batch mode.

Usage:
    # Single PDB ID with explicit ligands
    python download_and_extract.py --pdb-id 3ODU --uniprot P61073 --ligands STI

    # Single PDB ID with auto-ligand detection
    python download_and_extract.py --pdb-id 3ODU --uniprot P61073

    # CSV batch mode
    python download_and_extract.py --csv holostructures.csv --uniprot P61073

    # Single PDB with manual chain and custom output
    python download_and_extract.py --pdb-id 3ODU --chain A --ligands STI --output-dir ./out
"""

import argparse
import csv
import json
import os
import re
import sys
import urllib.request

# ── Auto-ligand filtering ─────────────────────────────────────────────
# When no --ligands specified and no ligand_comp_ids in CSV, keep all
# HETATM residues EXCEPT these common non-ligand artifacts.
EXCLUDE_RESIDUES = {
    # Water
    "HOH", "DOD", "WAT",
    # Monoatomic ions
    "NA", "K", "CL", "MG", "CA", "ZN", "FE", "MN", "CU", "CO",
    "CD", "NI", "BR", "I", "F", "LI", "RB", "CS", "SR", "BA",
    "PT", "AU", "HG", "SE",
    # Buffer / crystallization components
    "EDO", "GOL", "TRS", "PEG", "PGE", "MPD", "ACT", "CIT",
    "BME", "DMS", "SO4", "PO4", "FMT", "EOH", "IMD", "LDA",
    "BCT", "BOG", "DMU", "OLC", "PLM", "SCN", "NO3", "NH4",
    "TAM", "HEP", "MES", "PIP", "BIS", "TRI", "BTB",
    # Cryoprotectants
    "GLY", "SUC", "DTT", "TCE", "DTE", "SGM",
    # Common artifacts / unknown
    "UNK", "UNX", "ACE", "NH2", "FOR", "DUM",
    # Membrane mimetic lipids (exclude unless explicitly requested)
    "LDA", "LMT", "OLA", "SDS", "LMU", "LMN", "TGL", "D10",
    "DMX", "OCT", "HEX", "LHG",
}


def graphql(query):
    """Execute a GraphQL query against RCSB."""
    url = "https://data.rcsb.org/graphql"
    data = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(url, data=data,
                                  headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())


def get_chains_for_uniprot(pdb_id, uniprot_target):
    """Return list of chain IDs matching a UniProt accession via RCSB GraphQL."""
    query = """
    {
      entry(entry_id: "%s") {
        polymer_entities {
          entity_poly {
            pdbx_strand_id
          }
          rcsb_polymer_entity_container_identifiers {
            uniprot_ids
          }
        }
      }
    }
    """ % pdb_id

    result = graphql(query)
    entry = (result.get("data") or {}).get("entry")
    if not entry:
        return []

    chains = []
    for entity in entry.get("polymer_entities") or []:
        uniprot_ids = (
            entity.get("rcsb_polymer_entity_container_identifiers") or {}
        ).get("uniprot_ids") or []
        if uniprot_target in uniprot_ids:
            strand_ids = (
                (entity.get("entity_poly") or {}).get("pdbx_strand_id") or ""
            ).split(",")
            chains.extend(s.strip() for s in strand_ids if s.strip())

    return chains


def cif_to_pdb_lines(cif_path):
    """Parse mmCIF atom_site loop, yield PDB-format ATOM/HETATM lines.

    Uses auth_* fields (author-assigned chain/residue/atom names) so the
    output matches what a legacy PDB file would contain.
    """
    with open(cif_path) as fh:
        lines = fh.readlines()

    # Locate the _atom_site loop and collect column definitions
    col_map = {}
    atom_lines = []
    in_loop = False
    i = 0

    while i < len(lines):
        line = lines[i].rstrip("\n")
        if line.startswith("_atom_site.") and not in_loop:
            col_names = []
            j = i
            while j < len(lines) and lines[j].startswith("_atom_site."):
                col_names.append(lines[j].split(".")[-1].strip())
                j += 1
            col_map = {name: idx for idx, name in enumerate(col_names)}
            in_loop = True
            i = j
            continue
        if in_loop:
            stripped = line.strip()
            if (not stripped or stripped.startswith("_") or
                    stripped.startswith("loop_")):
                break
            if not stripped.startswith("#"):
                atom_lines.append(line)
        i += 1

    def col(name):
        return col_map.get(name)

    # Column indices
    atom_id_idx = col("id")
    group_idx = col("group_PDB")
    type_idx = col("type_symbol")
    alt_loc_idx = col("label_alt_id")
    ins_code_idx = col("pdbx_PDB_ins_code")
    x_idx = col("Cartn_x")
    y_idx = col("Cartn_y")
    z_idx = col("Cartn_z")
    occ_idx = col("occupancy")
    b_idx = col("B_iso_or_equiv")
    elem_idx = col("type_symbol")
    # Author-assigned for PDB-like output
    atm_name_idx = col("auth_atom_id") or col("label_atom_id")
    res_name_idx = col("auth_comp_id")
    chain_idx = col("auth_asym_id")
    res_seq_idx = col("auth_seq_id")

    for line in atom_lines:
        fields = line.split()
        while len(fields) < max(col_map.values()) + 1:
            fields.append("?")

        if x_idx is None or y_idx is None or z_idx is None:
            continue

        def get(idx, default="?"):
            if idx is None or idx >= len(fields):
                return default
            val = fields[idx]
            return val if val not in (".", "?") else ""

        group = get(group_idx, "ATOM")
        atm_id = get(atom_id_idx, "")
        atm_name = get(atm_name_idx, "")
        alt_loc = get(alt_loc_idx, "")
        res_name = get(res_name_idx, "")
        # Truncate to 3 chars for PDB fixed-width columns 18-20.
        # CIF auth_comp_id can be up to 5 chars (e.g. "A1EIE"), which
        # would spill into the chain-ID column (22) and corrupt parsing.
        res_name = res_name[:3]
        chain = get(chain_idx, "")
        res_seq = get(res_seq_idx, "")
        ins_code = get(ins_code_idx, "")
        x = get(x_idx, "0.000")
        y = get(y_idx, "0.000")
        z = get(z_idx, "0.000")
        occ = get(occ_idx, "1.00")
        b = get(b_idx, "0.00")
        elem = get(elem_idx, "")

        # Pad atom name to 4 chars (left-justify for 2-letter elements)
        if len(atm_name) < 4:
            atm_name = (f"{atm_name:<4}" if len(elem) == 2
                        else f" {atm_name:<3}")
        atm_name = atm_name[:4]

        rec = group[:6]
        try:
            serial = int(float(atm_id)) if atm_id else 0
        except (ValueError, TypeError):
            serial = 0
        serial = serial % 100000

        chain = (chain or " ")[0]
        ins_code = (ins_code or " ")[0]
        alt_loc = (alt_loc or " ")[0]

        try:
            res_seq_i = int(float(res_seq))
        except (ValueError, TypeError):
            res_seq_i = 0

        pdb_line = (
            f"{rec:<6}"
            f"{serial:>5} "
            f"{atm_name:<4}"
            f"{alt_loc}"
            f"{res_name:>3}"
            f" {chain}"
            f"{res_seq_i:>4}"
            f"{ins_code}   "
            f"{float(x):>8.3f}"
            f"{float(y):>8.3f}"
            f"{float(z):>8.3f}"
            f"{float(occ):>6.2f}"
            f"{float(b):>6.2f}          "
            f"{elem:>2}"
        )
        # Pad to exactly 80 columns
        pdb_line = pdb_line[:79].ljust(79) + "\n"
        yield pdb_line


def download_coords(pdb_id, output_dir):
    """Download PDB coordinates; fall back to CIF → PDB conversion if needed."""
    pdb_dir = os.path.join(output_dir, pdb_id)
    os.makedirs(pdb_dir, exist_ok=True)
    pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")

    if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 0:
        print(f"    [skip] {pdb_path} already exists", file=sys.stderr)
        return pdb_path

    # Try legacy PDB format first
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"    Trying {pdb_id}.pdb ...", file=sys.stderr)
    try:
        urllib.request.urlretrieve(url, pdb_path)
        if os.path.getsize(pdb_path) > 0:
            return pdb_path
    except Exception:
        pass

    # Fall back to mmCIF, convert on-the-fly
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"    PDB unavailable; downloading {pdb_id}.cif and converting ...",
          file=sys.stderr)
    cif_path = os.path.join(pdb_dir, f"{pdb_id}.cif")
    try:
        urllib.request.urlretrieve(url, cif_path)
        if os.path.getsize(cif_path) == 0:
            raise Exception("empty CIF file")
    except Exception as e:
        print(f"    [error] Download failed: {e}", file=sys.stderr)
        return None

    with open(pdb_path, "w", encoding="utf-8") as f_out:
        for pdb_line in cif_to_pdb_lines(cif_path):
            f_out.write(pdb_line)
        f_out.write("END\n")

    return pdb_path


def extract_monomer_ligand(pdb_path, pdb_id, target_chain, target_comp_ids,
                           output_dir, auto_ligand=False):
    """Extract one chain + matching ligands, write clean PDB. Returns stats."""
    out_path = os.path.join(output_dir, pdb_id, f"{pdb_id}_clean.pdb")
    target_set = set(c.strip().upper()[:3] for c in target_comp_ids)

    atom_count = 0
    kept_ligands = set()

    with open(pdb_path, encoding="utf-8") as f_in, open(out_path, "w", encoding="utf-8") as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                if line[21] == target_chain:
                    f_out.write(line)
                    atom_count += 1
            elif line.startswith("HETATM"):
                res_3 = line[17:20].strip().upper()
                if auto_ligand:
                    # Keep everything except excluded artifacts, target chain only
                    if (line[21] == target_chain and
                            res_3 not in EXCLUDE_RESIDUES):
                        f_out.write(line)
                        atom_count += 1
                        kept_ligands.add(res_3)
                else:
                    # Keep only explicitly requested ligands on target chain
                    res_5 = line[17:22].strip().upper()
                    if (line[21] == target_chain and
                            (res_3 in target_set or res_5 in target_set)):
                        f_out.write(line)
                        atom_count += 1
                        kept_ligands.add(res_3 or res_5)
            elif line.startswith("END"):
                break

        f_out.write("END\n")

    return out_path, atom_count, kept_ligands


def parse_ligand_arg(ligands_str):
    """Split --ligands or CSV ligand_comp_ids into a list of codes."""
    if not ligands_str:
        return []
    # Support commas, semicolons, spaces
    return [c.strip() for c in re.split(r'[;, ]+', ligands_str) if c.strip()]


# ── Main ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Download PDB structures from RCSB, extract target "
                    "monomer + bound ligands, and write clean PDB files."
    )
    parser.add_argument("--pdb-id", help="Single PDB identifier (4 chars)")
    parser.add_argument("--uniprot", help="UniProt accession for target protein")
    parser.add_argument("--csv", help="CSV file with columns: pdb_id, "
                        "ligand_comp_ids (optional)")
    parser.add_argument("--chain", help="Manual chain ID override "
                        "(skip GraphQL resolution)")
    parser.add_argument("--ligands", help="3-letter residue codes "
                        "(comma/semicolon separated)")
    parser.add_argument("--output-dir", default="raw",
                        help="Output root directory (default: raw/)")
    parser.add_argument("--csv-uniprot-col", default="uniprot_ids",
                        help="CSV column containing UniProt IDs (for "
                        "multi-protein entries, default: uniprot_ids)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Build entry list ──────────────────────────────────────────────
    entries = []

    if args.pdb_id:
        entries.append({
            "pdb_id": args.pdb_id,
            "ligand_comp_ids": args.ligands or "",
            "uniprot": args.uniprot,
        })
    elif args.csv:
        with open(args.csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                pdb_id = (row.get("pdb_id") or "").strip()
                if not pdb_id:
                    continue
                # If CSV has uniprot_ids, use first one as target (or
                # use --uniprot as override)
                csv_uniprot = (row.get(args.csv_uniprot_col) or "").strip()
                ligand_str = row.get("ligand_comp_ids", "")
                entries.append({
                    "pdb_id": pdb_id,
                    "ligand_comp_ids": ligand_str,
                    "uniprot": args.uniprot or csv_uniprot.split(";")[0].strip(),
                })
    else:
        parser.error("Either --pdb-id or --csv is required")

    if not entries:
        print("[error] No PDB entries to process", file=sys.stderr)
        sys.exit(1)

    # ── Process entries ───────────────────────────────────────────────
    print(f"Processing {len(entries)} PDB entries\n", file=sys.stderr)

    success_count = 0
    for entry in entries:
        pdb_id = entry["pdb_id"]
        ligand_codes = parse_ligand_arg(entry.get("ligand_comp_ids", ""))
        auto_ligand = not bool(ligand_codes)

        print(f"{pdb_id}:", file=sys.stderr)

        # 1. Resolve chain(s)
        if args.chain:
            target_chain = args.chain
            print(f"    Chain: {target_chain} (manual override)",
                  file=sys.stderr)
        elif entry.get("uniprot"):
            chains = get_chains_for_uniprot(pdb_id, entry["uniprot"])
            if not chains:
                print(f"    [warn] No chains matching {entry['uniprot']} "
                      f"found, skipping", file=sys.stderr)
                continue
            target_chain = chains[0]
            print(f"    Chains matching {entry['uniprot']}: {chains} "
                  f"→ using {target_chain}", file=sys.stderr)
        else:
            print(f"    [error] No --uniprot or --chain specified, "
                  f"skipping", file=sys.stderr)
            continue

        # 2. Download coordinates
        pdb_path = download_coords(pdb_id, args.output_dir)
        if not pdb_path:
            continue

        # 3. Extract monomer + ligands
        out_path, atom_count, kept = extract_monomer_ligand(
            pdb_path, pdb_id, target_chain, ligand_codes,
            args.output_dir, auto_ligand=auto_ligand
        )
        mode = "auto-detect" if auto_ligand else "requested"
        print(f"    Ligands ({mode}): {ligand_codes or 'all non-artifact'}  "
              f"kept: {sorted(kept)}", file=sys.stderr)
        print(f"    → {out_path}  ({atom_count} atoms)", file=sys.stderr)
        print(file=sys.stderr)
        success_count += 1

    print(f"Done. {success_count}/{len(entries)} entries processed.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
