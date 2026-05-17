#!/usr/bin/env python3
"""Search the RCSB PDB for holostructures (protein-ligand complexes) of a target protein and its homologs.

Reads a target (UniProt accession, PDB ID, sequence, or name), runs a sequence-similarity search
against the PDB at a chosen identity cutoff, requires at least one biologically meaningful ligand
(default: RCSB's SUBJECT_OF_INVESTIGATION annotation, fallback: MW > 150 Da + manual exclusion list),
and emits a CSV table with PDB ID, resolution, method, UniProt IDs, ligand chem-comp IDs, names, MWs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import urllib.parse
import urllib.request

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/graphql"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

UNIPROT_ACCESSION_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)
PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
AA_ONLY_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYX*]+$")

# Manual exclusion list — only used when --strict-mw is set (LOI flag handles this otherwise).
# Covers waters, common ions, buffers, cryoprotectants, crystallization additives, and frequently
# non-functional sugars. Keep this conservative; over-excluding hides real cofactors.
EXCLUDED_LIGANDS = {
    # Solvents & cryoprotectants
    "HOH", "DOD", "GOL", "EDO", "PEG", "PG4", "PGE", "PGO", "MPD",
    "DMS", "DMF", "EOH", "MOH", "BME", "ACE", "ACT", "FMT", "ACY",
    # Buffers
    "TRS", "BTB", "EPE", "HEP", "TAM", "MES", "BIS", "MOPS",
    "PIPE", "TAPS", "HEPES", "CIT", "CAC",
    # Ions / small inorganics
    "NA", "K", "MG", "CA", "ZN", "FE", "FE2", "MN", "CU", "CU1",
    "CO", "NI", "CD", "HG", "CS", "RB", "LI", "BA", "SR",
    "CL", "BR", "I", "F", "SO4", "PO4", "NO3", "CO3", "OH",
    "NH4", "BO3", "BO4", "EDT", "IOD", "AZI", "OXY", "PER", "CN",
    # Common crystallization sugars (often non-functional packing artifacts)
    "NAG", "MAN", "BMA", "FUC", "GAL", "GLC", "BGC", "FUL", "XYL", "SUC", "TRE",
    # Unknown / placeholder
    "UNX", "UNL", "UNK", "UNI",
}


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------

def http_post_json(url: str, payload: dict, timeout: int = 60) -> dict:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={"Content-Type": "application/json", "User-Agent": "pdb-holostructure-search/1.0"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8")
        # The Search API returns 204 No Content (empty body) for zero-hit queries —
        # treat that as an empty result set rather than a parse error.
        if not body.strip():
            return {"result_set": [], "total_count": 0}
        return json.loads(body)


def http_get_text(url: str, timeout: int = 30) -> str:
    req = urllib.request.Request(
        url, headers={"User-Agent": "pdb-holostructure-search/1.0"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


# ---------------------------------------------------------------------------
# Target resolution
# ---------------------------------------------------------------------------

def detect_target_type(target: str) -> str:
    """Heuristically classify the user's target string."""
    t = target.strip()
    if t.startswith(">") or "\n" in t:
        return "sequence"
    if " " in t:
        return "name"
    if UNIPROT_ACCESSION_RE.match(t):
        return "uniprot"
    if PDB_ID_RE.match(t):
        return "pdb"
    if len(t) >= 20 and AA_ONLY_RE.match(t.upper()):
        return "sequence"
    return "name"


def parse_fasta(text: str) -> str:
    lines = [ln.strip() for ln in text.splitlines()
             if ln.strip() and not ln.startswith(">")]
    return "".join(lines).upper()


def _uniprot_search(query: str, size: int = 5) -> list[list[str]]:
    """Run a UniProt REST search and return rows of [accession, id, protein_name, organism_name, gene_primary]."""
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,id,protein_name,organism_name,gene_primary",
        "size": size,
    }
    url = UNIPROT_SEARCH_URL + "?" + urllib.parse.urlencode(params)
    text = http_get_text(url)
    rows = text.strip().splitlines()
    if len(rows) <= 1:
        return []
    return [r.split("\t") for r in rows[1:]]


def _is_gene_symbol_like(s: str) -> bool:
    """Heuristic: does this string look like a gene symbol/acronym rather than a descriptive name?

    Gene symbols are short (≤7 chars), contain no spaces, and are mostly alphanumeric.
    Examples: "EGFR", "ATR", "TP53", "BRCA1", "HER2".
    """
    s = s.strip()
    if " " in s:
        return False
    if not (2 <= len(s) <= 7):
        return False
    # Allow letters, digits, hyphens, underscores (e.g. "HER2", "IL-6", "HLA_A")
    return bool(re.match(r"^[A-Za-z0-9][A-Za-z0-9\-_]*$", s))


def resolve_name_to_uniprot(name: str) -> str | None:
    """Resolve a protein or gene name to a UniProt accession.

    Tries progressively narrower queries to avoid the all-too-easy mistake of getting back
    a similarly-named-but-wrong protein (e.g., "thymidine kinase" matching "thymidylate kinase"
    or "ATR" matching "ATR-interacting protein").

    For short, no-space queries that look like gene symbols ("EGFR", "ATR", "TP53"), gene-name
    search is tried *first* because these are overwhelmingly gene-symbol queries and the
    protein-name search is prone to substring matches against unrelated proteins
    ("ATR" → "ATR-interacting protein").  For multi-word queries ("thymidine kinase") the
    protein-name search remains first because those are descriptive names.

    Always logs the candidates so the user can see what was picked and what was passed over —
    and warns when the top hit's protein name doesn't actually contain the query string.
    """
    name = name.strip()

    if _is_gene_symbol_like(name):
        queries = [
            ('gene name',          f'gene:{name} AND reviewed:true'),
            ('exact protein_name', f'protein_name:"{name}" AND reviewed:true'),
            ('general',            f'{name} AND reviewed:true'),
        ]
    else:
        queries = [
            ('exact protein_name', f'protein_name:"{name}" AND reviewed:true'),
            ('gene name',          f'gene:{name} AND reviewed:true'),
            ('general',            f'{name} AND reviewed:true'),
        ]

    rows: list[list[str]] = []
    used_label = None
    for label, q in queries:
        try:
            rows = _uniprot_search(q, size=5)
        except Exception as exc:
            sys.stderr.write(f"[warn] UniProt query ({label}) failed: {exc}\n")
            continue
        if rows:
            used_label = label
            break
    if not rows:
        return None

    # --- Post-filter: for gene-symbol queries, prefer the result whose gene_primary
    # --- matches the query exactly.  This catches cases where gene-name search returned
    # --- multiple results and the top hit by relevance isn't the one the user wants.
    if _is_gene_symbol_like(name):
        gene_exact = [r for r in rows if len(r) >= 5 and r[4].upper() == name.upper()]
        if gene_exact and gene_exact[0][0] != rows[0][0]:
            preferred = gene_exact[0]
            sys.stderr.write(
                f"[info] Preferring gene-exact match {preferred[0]} ({preferred[2][:60] if len(preferred) > 2 else ''}) "
                f"over relevance-top {rows[0][0]} ({rows[0][2][:60] if len(rows[0]) > 2 else ''}).\n"
            )
            # Move the preferred row to the front
            rows = sorted(rows, key=lambda r: (r[0] != preferred[0], rows.index(r)))

    sys.stderr.write(f"[info] UniProt name lookup ('{name}', match strategy: {used_label}):\n")
    for i, r in enumerate(rows[:5]):
        marker = "→" if i == 0 else " "
        gene_info = f" gene={r[4]}" if len(r) > 4 and r[4] else ""
        sys.stderr.write(f"[info]   {marker} {r[0]:<10s}{gene_info:<14s} {r[2] if len(r) > 2 else '':<60s}"
                         f" [{r[3] if len(r) > 3 else '?'}]\n")
    top = rows[0]
    top_pname = top[2].lower() if len(top) > 2 else ""
    top_gene = (top[4] or "").lower() if len(top) > 4 else ""
    if used_label != "exact protein_name" and name.lower() not in top_pname:
        # Suppress warning if the gene name matches (common for gene-symbol queries like "EGFR"
        # where the protein name doesn't literally contain the gene symbol).
        if top_gene != name.lower():
            sys.stderr.write(
                f"[warn] Top match's protein name does not contain '{name}'. "
                "Consider passing an explicit UniProt accession to avoid wrong-protein hits.\n"
            )
    return top[0]


def fetch_uniprot_sequence(accession: str) -> str:
    text = http_get_text(UNIPROT_FASTA_URL.format(acc=accession))
    seq = parse_fasta(text)
    if not seq:
        raise RuntimeError(f"Empty FASTA for UniProt {accession}")
    return seq


def fetch_pdb_protein_sequence(pdb_id: str) -> str:
    """Fetch the longest protein chain's sequence from a PDB entry.

    A PDB entry often contains multiple polymers — e.g., 1ATP has the 350-residue PKA catalytic
    subunit alongside a 20-residue PKI inhibitor peptide. Picking the first polypeptide gives
    the wrong answer for any entry where the catalytic chain isn't entity 1, so we always pick
    the longest protein chain. If the user wants a different chain, they should pass the
    sequence directly.
    """
    query = """
    query($pdb_id: String!) {
      entry(entry_id: $pdb_id) {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers { entity_id }
          entity_poly { pdbx_seq_one_letter_code_can type }
          rcsb_polymer_entity { pdbx_description }
        }
      }
    }
    """
    resp = http_post_json(DATA_URL, {"query": query, "variables": {"pdb_id": pdb_id.upper()}})
    entry = (resp.get("data") or {}).get("entry") or {}
    candidates: list[tuple[int, str, str, str]] = []  # (length, entity_id, description, sequence)
    for ent in entry.get("polymer_entities") or []:
        ep = ent.get("entity_poly") or {}
        if not (ep.get("type") or "").startswith("polypeptide"):
            continue
        seq = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
        if not seq:
            continue
        eid = ((ent.get("rcsb_polymer_entity_container_identifiers") or {})
               .get("entity_id") or "?")
        desc = ((ent.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "")
        candidates.append((len(seq), eid, desc, seq.upper()))
    if not candidates:
        raise RuntimeError(f"No protein chain found in PDB ID {pdb_id}")
    candidates.sort(reverse=True)
    if len(candidates) > 1:
        sys.stderr.write(f"[info] PDB {pdb_id.upper()} has {len(candidates)} protein chains; "
                         f"using longest:\n")
        for length, eid, desc, _ in candidates[:5]:
            marker = "→" if (length, eid) == (candidates[0][0], candidates[0][1]) else " "
            sys.stderr.write(f"[info]   {marker} entity {eid}: {length} aa — {desc[:80]}\n")
    return candidates[0][3]


# ---------------------------------------------------------------------------
# Search query
# ---------------------------------------------------------------------------

def build_search_query(
    sequence: str,
    identity_cutoff: float,
    loi_only: bool,
    method: str | None,
    max_resolution: float | None,
) -> dict:
    nodes: list[dict] = [
        {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "value": sequence,
                "sequence_type": "protein",
                "identity_cutoff": identity_cutoff,
                "evalue_cutoff": 0.1,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                "operator": "greater",
                "value": 0,
            },
        },
    ]

    if loi_only:
        # RCSB's curated SUBJECT_OF_INVESTIGATION flag — high precision but sparse on older
        # entries (annotation programme rolled out ~2020).
        nodes.append({
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_nonpolymer_entity_annotation.type",
                "operator": "exact_match",
                "value": "SUBJECT_OF_INVESTIGATION",
            },
        })
    else:
        # Default: at the server we only require that *some* ligand has MW > 150 Da. We do NOT
        # apply the exclusion list server-side, because `negation:true` on `nonpolymer_comp_id in
        # [excluded]` excludes an entry if it contains *any* listed comp — and since virtually
        # every X-ray entry has waters/ions, that filter wipes out the result set. The exclusion
        # list is therefore applied per-ligand at the metadata stage instead, which yields
        # exactly the right behaviour (drop junk ligands, keep entries that have ≥1 real one).
        nodes.append({
            "type": "terminal",
            "service": "text_chem",
            "parameters": {
                "attribute": "chem_comp.formula_weight",
                "operator": "greater",
                "value": 150,
            },
        })

    if method:
        nodes.append({
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "exptl.method",
                "operator": "exact_match",
                "value": method,
            },
        })
    if max_resolution is not None:
        nodes.append({
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal",
                "value": max_resolution,
            },
        })

    return {
        "return_type": "entry",
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "request_options": {
            "results_content_type": ["experimental"],
            "return_all_hits": True,
            "sort": [
                {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"},
            ],
        },
    }


# ---------------------------------------------------------------------------
# Per-entry metadata via Data API
# ---------------------------------------------------------------------------

DATA_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    rcsb_entry_info {
      resolution_combined
      experimental_method
    }
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers { nonpolymer_comp_id }
      pdbx_entity_nonpoly { name }
      rcsb_nonpolymer_entity_annotation { type }
      nonpolymer_comp { chem_comp { formula_weight type } }
    }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { uniprot_ids }
    }
  }
}
"""


def fetch_entry_metadata(pdb_ids: list[str], loi_only: bool) -> dict[str, dict]:
    """Pull ligand and entity metadata for a list of PDB IDs and apply per-ligand filtering.

    The same filter we used at search time is re-applied here at the per-ligand level so the
    output's ligand columns only contain the meaningful ones (a structure with both ATP and
    sodium will yield only ATP).
    """
    out: dict[str, dict] = {}
    if not pdb_ids:
        return out

    for i in range(0, len(pdb_ids), 50):
        batch = pdb_ids[i:i + 50]
        try:
            resp = http_post_json(DATA_URL, {"query": DATA_QUERY, "variables": {"ids": batch}})
        except Exception as exc:
            sys.stderr.write(f"[warn] Data API batch failed ({batch[0]}..{batch[-1]}): {exc}\n")
            continue
        for entry in (resp.get("data") or {}).get("entries") or []:
            pid = entry.get("rcsb_id")
            if not pid:
                continue
            ligands: list[dict] = []
            for npe in entry.get("nonpolymer_entities") or []:
                cid = ((npe.get("rcsb_nonpolymer_entity_container_identifiers") or {})
                       .get("nonpolymer_comp_id"))
                if not cid:
                    continue
                annotations = npe.get("rcsb_nonpolymer_entity_annotation") or []
                is_loi = any(a.get("type") == "SUBJECT_OF_INVESTIGATION" for a in annotations)
                cc = ((npe.get("nonpolymer_comp") or {}).get("chem_comp") or {})
                mw = cc.get("formula_weight")
                name = (npe.get("pdbx_entity_nonpoly") or {}).get("name") or ""
                if loi_only:
                    if not is_loi:
                        continue
                else:
                    if cid in EXCLUDED_LIGANDS:
                        continue
                    if mw is None or mw <= 150:
                        continue
                ligands.append({"comp_id": cid, "name": name, "mw": mw, "is_loi": is_loi})

            uniprots: list[str] = []
            for pe in entry.get("polymer_entities") or []:
                ids_ = ((pe.get("rcsb_polymer_entity_container_identifiers") or {})
                        .get("uniprot_ids") or [])
                uniprots.extend(ids_)
            uniprots = sorted(set(uniprots))

            res_list = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
            resolution = res_list[0] if res_list else None
            method = (entry.get("rcsb_entry_info") or {}).get("experimental_method") or ""
            title = (entry.get("struct") or {}).get("title") or ""
            out[pid] = {
                "title": title,
                "resolution": resolution,
                "method": method,
                "ligands": ligands,
                "uniprot_ids": uniprots,
            }
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("target", help="UniProt accession, PDB ID, FASTA/sequence, or protein name.")
    p.add_argument("--target-type", choices=["uniprot", "pdb", "sequence", "name", "auto"],
                   default="auto", help="Override auto-detection of TARGET.")
    p.add_argument("--identity", type=float, default=0.90,
                   help="Sequence-identity cutoff for homologs (0.0–1.0). Default 0.90.")
    p.add_argument("--loi-only", action="store_true",
                   help="Use RCSB's curated SUBJECT_OF_INVESTIGATION flag for ligand filtering. "
                        "More precise but undercounts pre-2020 entries that lack the annotation. "
                        "Default is MW>150 Da + manual exclusion list, which works on every entry.")
    p.add_argument("--method",
                   help='Filter by experimental method (e.g., "X-RAY DIFFRACTION").')
    p.add_argument("--max-resolution", type=float, help="Maximum resolution in Å.")
    p.add_argument("--output", "-o", default="-", help="Output CSV path (default: stdout).")
    p.add_argument("--max-results", type=int, help="Cap number of results.")
    p.add_argument("--query-only", action="store_true",
                   help="Print the JSON search query and exit (no API calls run).")
    args = p.parse_args()

    if not (0.0 < args.identity <= 1.0):
        sys.stderr.write("[error] --identity must be in (0, 1].\n")
        return 2

    target_type = args.target_type if args.target_type != "auto" else detect_target_type(args.target)
    sys.stderr.write(f"[info] Target '{args.target[:60]}{'…' if len(args.target) > 60 else ''}' classified as: {target_type}\n")

    try:
        if target_type == "uniprot":
            sequence = fetch_uniprot_sequence(args.target.strip())
        elif target_type == "pdb":
            sequence = fetch_pdb_protein_sequence(args.target.strip())
        elif target_type == "sequence":
            sequence = parse_fasta(args.target) if ">" in args.target else \
                "".join(args.target.split()).upper()
        elif target_type == "name":
            acc = resolve_name_to_uniprot(args.target.strip())
            if not acc:
                sys.stderr.write(f"[error] Could not resolve name '{args.target}' to a UniProt entry.\n")
                return 1
            sequence = fetch_uniprot_sequence(acc)
        else:
            sys.stderr.write(f"[error] Unknown target type: {target_type}\n")
            return 2
    except Exception as exc:
        sys.stderr.write(f"[error] Could not obtain reference sequence: {exc}\n")
        return 1

    if len(sequence) < 20:
        sys.stderr.write(f"[error] Reference sequence too short ({len(sequence)} residues).\n")
        return 1
    sys.stderr.write(f"[info] Reference sequence: {len(sequence)} residues\n")

    query = build_search_query(
        sequence=sequence,
        identity_cutoff=args.identity,
        loi_only=args.loi_only,
        method=args.method,
        max_resolution=args.max_resolution,
    )

    if args.query_only:
        json.dump(query, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    sys.stderr.write(f"[info] Running RCSB search (identity ≥ {args.identity:.0%}, "
                     f"ligand filter = {'SUBJECT_OF_INVESTIGATION (curated)' if args.loi_only else 'MW>150 + exclusion list'})…\n")
    try:
        resp = http_post_json(SEARCH_URL, query)
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")[:500]
        sys.stderr.write(f"[error] Search API HTTP {exc.code}: {body}\n")
        return 1
    except Exception as exc:
        sys.stderr.write(f"[error] Search API request failed: {exc}\n")
        return 1

    hits = resp.get("result_set") or []
    pdb_ids = [h["identifier"] for h in hits]
    if args.max_results:
        pdb_ids = pdb_ids[:args.max_results]
    sys.stderr.write(f"[info] Search returned {len(pdb_ids)} entries.\n")

    if not pdb_ids:
        sys.stderr.write("[info] No matching holostructures found.\n")
        return 0

    sys.stderr.write("[info] Fetching ligand and entity metadata…\n")
    metadata = fetch_entry_metadata(pdb_ids, loi_only=args.loi_only)

    out_fh = sys.stdout if args.output == "-" else open(args.output, "w", newline="", encoding="utf-8")
    try:
        writer = csv.writer(out_fh)
        writer.writerow([
            "pdb_id", "resolution_A", "method", "title",
            "uniprot_ids", "ligand_comp_ids", "ligand_names", "ligand_mw_Da",
        ])
        written = 0
        for pid in pdb_ids:
            m = metadata.get(pid)
            if not m or not m["ligands"]:
                continue
            writer.writerow([
                pid,
                f"{m['resolution']:.2f}" if isinstance(m["resolution"], (int, float)) else "",
                m["method"],
                m["title"],
                ";".join(m["uniprot_ids"]),
                ";".join(l["comp_id"] for l in m["ligands"]),
                ";".join(l["name"] for l in m["ligands"]),
                ";".join(f"{l['mw']:.1f}" if isinstance(l["mw"], (int, float)) else "" for l in m["ligands"]),
            ])
            written += 1
        sys.stderr.write(f"[info] Wrote {written} rows ({len(pdb_ids) - written} entries dropped after metadata filter).\n")
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
