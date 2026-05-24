#!/usr/bin/env python3
"""
extract_report_data.py — bundled with the pharmacophore-report-generator skill.

Walks a project directory of per-PDB pharmacophore JSON files (typically the
output of the pharmacophore-analyzer skill), computes everything the report
needs — per-ligand retention, residue engagement, family retention, file
manifest, numbering-offset detection — and emits a single structured JSON that
the report-writing prompt can consume.

Optionally derives a multi-structure consensus pharmacophore via Kabsch
alignment of the orthosteric set's Cα backbones followed by greedy
single-linkage clustering of features within (2 × tolerance), classifying
clusters as mandatory (sources >= ceil(N/2)) or optional (sources < ceil(N/2)).

Usage:
    python extract_report_data.py \\
        --project-dir <path> \\
        --receptor-name "Receptor Name" \\
        --uniprot-id "PXXXXX" \\
        --output <path>/.report_data.json

See `--help` for all options.
"""
import argparse
import glob
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("numpy is required. Install with: pip install numpy")


# -------- Project layout detection --------

def find_pharmacophore_jsons(project_dir):
    """Find every *.pharmacophore.json in the project, categorize by suffix.

    Returns a list of dicts: {path, pdb_id, chain (or None), variant (filtered/all/unknown)}.
    Handles three common layouts:
      - canonical: pharmacophore/{pdb}/{pdb}{_chainX}_{filtered,all}.pharmacophore.json
      - flat ATR: pharmacophore/{pdb}_{filtered,all}.pharmacophore.json
      - legacy:   results/{filtered,unfiltered}/{pdb}/{pdb}{_chainX}.pharmacophore.json
    """
    paths = glob.glob(os.path.join(project_dir, "**/*.pharmacophore.json"),
                      recursive=True)
    entries = []
    for p in paths:
        rel = os.path.relpath(p, project_dir)
        fname = os.path.basename(p)
        # Strip .pharmacophore.json
        stem = fname[:-len(".pharmacophore.json")]
        # Detect variant
        if stem.endswith("_filtered"):
            variant = "filtered"
            stem = stem[:-len("_filtered")]
        elif stem.endswith("_all"):
            variant = "all"
            stem = stem[:-len("_all")]
        elif "filtered" in rel.lower().split(os.sep):
            variant = "filtered"
        elif "unfiltered" in rel.lower().split(os.sep):
            variant = "all"
        else:
            variant = "unknown"
        # Detect chain: stem like "8U4P_chainR" or just "8U4P"
        m = re.match(r"^([0-9A-Za-z]{4,5})(?:_chain([A-Za-z0-9]+))?(?:_lig\d+)?$",
                     stem)
        if m:
            pdb_id = m.group(1).upper()
            chain = m.group(2)
        else:
            pdb_id = stem.upper()
            chain = None
        entries.append({
            "path": p,
            "rel": rel,
            "stem": stem,
            "pdb_id": pdb_id,
            "chain": chain,
            "variant": variant,
        })
    return entries


def find_consensus_json(project_dir):
    """Look for an existing consensus_pharmacophore.json. Return path or None."""
    candidates = glob.glob(os.path.join(project_dir, "**/consensus_pharmacophore.json"),
                           recursive=True)
    return candidates[0] if candidates else None


def find_clean_pdb(project_dir, pdb_id):
    """Find <pdb_id>_clean.pdb anywhere in the project (raw/ preferred)."""
    paths = glob.glob(os.path.join(project_dir, "**", f"{pdb_id}_clean.pdb"),
                      recursive=True)
    # Resolve symlinks; deduplicate by realpath
    real = list({os.path.realpath(p) for p in paths if os.path.exists(p)})
    if not real:
        return None
    raw_first = sorted(real, key=lambda p: ("raw" not in p, p))
    return raw_first[0]


# -------- Per-PDB JSON loading --------

def load_jsons(entries):
    """Load every JSON, attach as entry['data']."""
    for e in entries:
        with open(e["path"]) as f:
            e["data"] = json.load(f)
    return entries


# -------- Ligand classification --------

def classify_ligands(entries, retention_threshold=0.10,
                       force_orthosteric=None, force_excluded=None):
    """Group filtered/all entries by ligand resname, compute mean retention,
    classify orthosteric (>= threshold) vs excluded.

    `force_orthosteric` / `force_excluded` (sets of resnames) override the
    auto-classification. Useful when bond-order failures depress retention
    for a chemotype that is clearly orthosteric on chemistry grounds (e.g.
    A1E in ATR with 18% retention but unambiguous ATP-pocket geometry).

    Returns:
        {
          "by_ligand": {
              "VH6": {"copies": 2, "all_total": 36, "filtered_total": 25,
                       "mean_retention": 0.69, "classification": "orthosteric",
                       "members": [{pdb, chain, all_count, filtered_count}, ...]},
              ...
          },
          "orthosteric_pdbs": [{"pdb_id": "8U4P", "chain": "R",
                                "ligand": "VH6", "filtered_path": "..."}, ...],
          "excluded_pdbs":    [...],
          "retention_threshold": 0.25,
        }
    """
    # Pair up filtered + all variants per PDB+chain
    pairs = defaultdict(lambda: {"filtered": None, "all": None})
    for e in entries:
        key = (e["pdb_id"], e["chain"])
        if e["variant"] in ("filtered", "all"):
            pairs[key][e["variant"]] = e

    by_ligand = defaultdict(lambda: {"copies": 0, "all_total": 0,
                                      "filtered_total": 0, "members": []})
    for (pdb_id, chain), variants in pairs.items():
        f = variants["filtered"]
        a = variants["all"]
        if f is None:
            continue  # need filtered for classification
        ligand_resname = f["data"]["ligand"]["resname"]
        f_count = len(f["data"]["features"])
        a_count = len(a["data"]["features"]) if a else f["data"]["metadata"].get(
            "n_features_before_filter", f_count)
        bucket = by_ligand[ligand_resname]
        bucket["copies"] += 1
        bucket["all_total"] += a_count
        bucket["filtered_total"] += f_count
        bucket["members"].append({
            "pdb_id": pdb_id, "chain": chain,
            "all_count": a_count, "filtered_count": f_count,
            "filtered_path": f["path"],
            "all_path": a["path"] if a else None,
        })

    force_orthosteric = set(force_orthosteric or [])
    force_excluded = set(force_excluded or [])
    orthosteric_pdbs = []
    excluded_pdbs = []
    for resname, bucket in by_ligand.items():
        retention = (bucket["filtered_total"] / bucket["all_total"]
                     if bucket["all_total"] else 0.0)
        bucket["mean_retention"] = retention
        if resname in force_orthosteric:
            label = "orthosteric"
            bucket["classification_source"] = "user-forced"
        elif resname in force_excluded:
            label = "excluded"
            bucket["classification_source"] = "user-forced"
        else:
            label = ("orthosteric" if retention >= retention_threshold
                      else "excluded")
            bucket["classification_source"] = "auto"
        bucket["classification"] = label
        for m in bucket["members"]:
            row = {"pdb_id": m["pdb_id"], "chain": m["chain"], "ligand": resname,
                   "filtered_path": m["filtered_path"],
                   "all_path": m["all_path"],
                   "filtered_count": m["filtered_count"],
                   "all_count": m["all_count"]}
            (orthosteric_pdbs if label == "orthosteric" else excluded_pdbs).append(row)

    return {
        "by_ligand": dict(by_ligand),
        "orthosteric_pdbs": orthosteric_pdbs,
        "excluded_pdbs": excluded_pdbs,
        "retention_threshold": retention_threshold,
    }


# -------- Residue engagement --------

def aggregate_residue_engagement(entries, orthosteric_pdbs):
    """For every filtered feature in every orthosteric PDB, walk the
    interaction.partner_residue field and aggregate by residue.

    Returns a dict keyed by (residue_name, residue_chain) with counts and
    contributing features.
    """
    ortho_keys = {(o["pdb_id"], o["chain"]): o["ligand"] for o in orthosteric_pdbs}
    by_residue = defaultdict(lambda: {
        "n_donor": 0, "n_acceptor": 0, "n_posion": 0, "n_negion": 0,
        "n_aromatic": 0, "n_hydrophobe": 0, "total": 0,
        "engaged_in": set(),  # set of pdb_ids
        "ligands": set(),
        "features": [],  # detailed list
    })

    for e in entries:
        if e["variant"] != "filtered":
            continue
        key = (e["pdb_id"], e["chain"])
        if key not in ortho_keys:
            continue
        ligand = ortho_keys[key]
        for feat in e["data"]["features"]:
            inter = feat.get("interaction")
            if not inter or not inter.get("partner_residue"):
                continue
            res = inter["partner_residue"]
            chain = inter.get("partner_chain", "?")
            r = by_residue[(res, chain)]
            family = feat["family"]
            r["total"] += 1
            r["engaged_in"].add(e["pdb_id"])
            r["ligands"].add(ligand)
            if family == "Donor":
                r["n_donor"] += 1
            elif family == "Acceptor":
                r["n_acceptor"] += 1
            elif family == "PosIonizable":
                r["n_posion"] += 1
            elif family == "NegIonizable":
                r["n_negion"] += 1
            elif family == "Aromatic":
                r["n_aromatic"] += 1
            elif family == "LumpedHydrophobe":
                r["n_hydrophobe"] += 1
            r["features"].append({
                "pdb_id": e["pdb_id"], "chain": e["chain"],
                "feature_id": feat["id"], "family": family,
                "type": feat.get("type"), "atom": inter.get("partner_atom"),
                "distance_A": inter.get("distance_A"),
                "interaction_type": inter.get("type"),
            })

    # Convert sets to sorted lists for JSON
    out = {}
    for (res, chain), r in by_residue.items():
        r["engaged_in"] = sorted(r["engaged_in"])
        r["ligands"] = sorted(r["ligands"])
        r["residue"] = res
        r["chain"] = chain
        out[f"{res}_{chain}"] = r

    return out


# -------- Family retention --------

def compute_family_retention(entries, orthosteric_pdbs):
    """Sum per-family raw vs filtered features across the orthosteric set."""
    ortho_keys = {(o["pdb_id"], o["chain"]) for o in orthosteric_pdbs}
    families = ["Donor", "PosIonizable", "Acceptor", "NegIonizable",
                "Aromatic", "LumpedHydrophobe"]
    counts = {fam: {"all": 0, "filtered": 0} for fam in families}

    for e in entries:
        key = (e["pdb_id"], e["chain"])
        if key not in ortho_keys:
            continue
        bucket = "filtered" if e["variant"] == "filtered" else (
            "all" if e["variant"] == "all" else None)
        if bucket is None:
            continue
        for feat in e["data"]["features"]:
            fam = feat["family"]
            if fam in counts:
                counts[fam][bucket] += 1

    out = {}
    for fam, c in counts.items():
        retention = c["filtered"] / c["all"] if c["all"] > 0 else None
        if retention is not None:
            inclusion = ("mandatory" if retention >= 0.5 else
                          "optional" if retention > 0 else "exclude")
        else:
            inclusion = "exclude" if c["filtered"] == 0 else "n/a"
        out[fam] = {"all": c["all"], "filtered": c["filtered"],
                     "retention": retention, "inclusion": inclusion}
    return out


# -------- Numbering offset detection --------

def detect_numbering_offsets(residue_engagement):
    """Heuristically detect when the same conserved residue appears under
    different PDB-specific numbering across the orthosteric set.

    Looks for residues with the same residue type (ASP/GLU/etc.) that are
    each engaged exclusively by disjoint PDB sets but with the same partner
    chemistry (likely an N-terminal construct offset).

    Returns a list of suspected aliases:
        [{"alias_a": "ASP262", "pdbs_a": ["8U4P"],
          "alias_b": "ASP270", "pdbs_b": ["8ZPL", "8ZPN"],
          "offset": 8, "residue_type": "ASP"}]
    """
    # Group by residue type
    by_type = defaultdict(list)
    for key, r in residue_engagement.items():
        m = re.match(r"^([A-Z]{3})(\d+)", r["residue"])
        if not m:
            continue
        rtype = m.group(1)
        rnum = int(m.group(2))
        by_type[rtype].append({"name": r["residue"], "num": rnum,
                                "pdbs": set(r["engaged_in"]),
                                "total": r["total"]})

    aliases = []
    for rtype, items in by_type.items():
        if len(items) < 2:
            continue
        # Look for pairs whose PDB sets are disjoint and within ±15 residue offset
        for i, a in enumerate(items):
            for b in items[i + 1:]:
                if a["pdbs"] & b["pdbs"]:
                    continue
                offset = abs(a["num"] - b["num"])
                if 1 <= offset <= 15:
                    # Likely the same residue under different construct numbering
                    aliases.append({
                        "alias_a": a["name"], "pdbs_a": sorted(a["pdbs"]),
                        "alias_b": b["name"], "pdbs_b": sorted(b["pdbs"]),
                        "offset": offset, "residue_type": rtype,
                        "combined_contacts": a["total"] + b["total"],
                    })
    return aliases


# -------- File manifest --------

def build_file_manifest(entries, classification):
    """Return the per-PDB file inventory with role annotation."""
    role_lookup = {}
    for o in classification["orthosteric_pdbs"]:
        role_lookup[(o["pdb_id"], o["chain"])] = "orthosteric_source"
    for e in classification["excluded_pdbs"]:
        role_lookup[(e["pdb_id"], e["chain"])] = "excluded"

    by_pdb = defaultdict(lambda: {"pdb_id": None, "chains": [], "files": []})
    for e in entries:
        key = (e["pdb_id"], e["chain"])
        bucket = by_pdb[e["pdb_id"]]
        bucket["pdb_id"] = e["pdb_id"]
        chain_key = e["chain"] or ""
        if chain_key not in bucket["chains"]:
            bucket["chains"].append(chain_key)
        ligand = e["data"]["ligand"]["resname"]
        f_count = (len(e["data"]["features"])
                    if e["variant"] == "filtered" else None)
        bucket["files"].append({
            "rel": e["rel"], "variant": e["variant"], "chain": e["chain"],
            "ligand": ligand, "n_features": len(e["data"]["features"]),
            "role": role_lookup.get(key),
        })
    return [by_pdb[p] for p in sorted(by_pdb)]


# -------- Consensus derivation: Kabsch alignment --------

def parse_ca_coords(pdb_path):
    """Parse Cα coordinates from a PDB file. Returns dict {(chain, resnum): (x,y,z)}."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") and line[12:16].strip() == "CA"):
                continue
            try:
                chain = line[21]
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords[(chain, resnum)] = (x, y, z)
            except (ValueError, IndexError):
                continue
    return coords


def kabsch_align(p, q):
    """Kabsch algorithm: compute rotation R and translation t that aligns
    points p (source) to q (target). p, q are (N, 3) numpy arrays.

    Returns (R, t) such that aligned_p = (p - centroid_p) @ R + centroid_q
    minimizes RMSD against q.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    cp = p.mean(axis=0)
    cq = q.mean(axis=0)
    H = (p - cp).T @ (q - cq)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = cq - cp @ R
    return R, t


def transform_points(points, R, t):
    return np.asarray(points) @ R + t


def rmsd(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


# -------- Consensus derivation: greedy clustering --------

def derive_consensus(entries, classification, project_dir,
                      cluster_factor=2.0):
    """Derive a multi-structure consensus pharmacophore.

    1. Pick reference PDB (orthosteric structure with most filtered features).
    2. Kabsch-align other orthosteric PDBs' Cα to the reference.
    3. Transform their features into the aligned frame.
    4. Greedy single-linkage cluster features by family within
       (cluster_factor × tolerance).
    5. Classify clusters as mandatory/optional based on source PDB count.

    Returns the consensus dict (same schema as consensus_pharmacophore.json).
    """
    ortho = classification["orthosteric_pdbs"]
    if len(ortho) < 2:
        return None

    # Sort by filtered_count desc → reference is most-feature-rich
    ortho_sorted = sorted(ortho, key=lambda o: -o["filtered_count"])
    reference = ortho_sorted[0]
    others = ortho_sorted[1:]

    # Load reference Cα
    ref_pdb = find_clean_pdb(project_dir, reference["pdb_id"])
    if not ref_pdb:
        sys.stderr.write(f"[warn] could not find clean.pdb for reference "
                         f"{reference['pdb_id']}; skipping consensus\n")
        return None
    ref_ca = parse_ca_coords(ref_pdb)

    # Load reference features (already in reference frame)
    aligned_features = []  # list of (source_member, transformed_position)
    ref_filtered_path = reference["filtered_path"]
    with open(ref_filtered_path) as f:
        ref_data = json.load(f)
    for feat in ref_data["features"]:
        aligned_features.append({
            "pdb": reference["pdb_id"],
            "feature_id": feat["id"],
            "family": feat["family"],
            "type": feat["type"],
            "position": feat["position"],
            "tolerance": feat["tolerance"],
        })

    alignment_summary = [{
        "pdb_id": reference["pdb_id"],
        "role": "reference",
        "common_ca": len(ref_ca),
        "rmsd_before_A": 0.0, "rmsd_after_A": 0.0,
    }]

    for other in others:
        oth_pdb = find_clean_pdb(project_dir, other["pdb_id"])
        if not oth_pdb:
            sys.stderr.write(f"[warn] no clean.pdb for {other['pdb_id']}; "
                             f"skipping in consensus\n")
            continue
        oth_ca = parse_ca_coords(oth_pdb)
        common = sorted(set(ref_ca) & set(oth_ca))
        if len(common) < 10:
            sys.stderr.write(f"[warn] only {len(common)} common Cα between "
                             f"{reference['pdb_id']} and {other['pdb_id']}; "
                             f"skipping in consensus\n")
            continue
        p = np.array([oth_ca[k] for k in common])
        q = np.array([ref_ca[k] for k in common])
        rmsd_before = rmsd(p, q)
        R, t = kabsch_align(p, q)
        p_aligned = transform_points(p, R, t)
        rmsd_after = rmsd(p_aligned, q)
        alignment_summary.append({
            "pdb_id": other["pdb_id"], "role": "aligned",
            "common_ca": len(common),
            "rmsd_before_A": round(rmsd_before, 2),
            "rmsd_after_A": round(rmsd_after, 2),
        })

        # Transform features
        with open(other["filtered_path"]) as f:
            oth_data = json.load(f)
        for feat in oth_data["features"]:
            new_pos = transform_points([feat["position"]], R, t)[0]
            aligned_features.append({
                "pdb": other["pdb_id"],
                "feature_id": feat["id"],
                "family": feat["family"],
                "type": feat["type"],
                "position": new_pos.tolist(),
                "tolerance": feat["tolerance"],
            })

    # Greedy single-linkage cluster by family within cluster_factor * tolerance
    clusters = []
    used = [False] * len(aligned_features)
    for i, fi in enumerate(aligned_features):
        if used[i]:
            continue
        used[i] = True
        members = [fi]
        for j, fj in enumerate(aligned_features):
            if used[j] or fj["family"] != fi["family"]:
                continue
            cutoff = cluster_factor * min(fi["tolerance"], fj["tolerance"])
            d = math.dist(members[-1]["position"], fj["position"])
            if d <= cutoff:
                used[j] = True
                members.append(fj)
        clusters.append(members)

    # Build consensus features (canonical schema: top-level classification,
    # n_sources, members; per-cluster radius_from_centroid)
    n_sources = len({f["pdb"] for f in aligned_features
                     if f["pdb"] in {o["pdb_id"] for o in ortho}})
    threshold = math.ceil(n_sources / 2)
    consensus_features = []
    for cid, members in enumerate(clusters):
        positions = np.array([m["position"] for m in members])
        centroid = positions.mean(axis=0)
        # Radius from centroid: max distance from centroid to any member
        if len(members) > 1:
            radius = float(np.linalg.norm(positions - centroid, axis=1).max())
        else:
            radius = 0.0
        sources = sorted({m["pdb"] for m in members})
        classification_label = ("mandatory" if len(sources) >= threshold
                                  else "optional")
        consensus_features.append({
            "id": cid,
            "family": members[0]["family"],
            "type": members[0]["type"],
            "position": [round(float(c), 3) for c in centroid],
            "tolerance": members[0]["tolerance"],
            "radius_from_centroid": round(radius, 3),
            "classification": classification_label,
            "n_sources": len(sources),
            "n_instances": len(members),
            "members": [{
                "pdb": m["pdb"], "id": m["feature_id"],
                "position": [round(float(c), 3) for c in m["position"]],
            } for m in members],
        })

    n_mandatory = sum(1 for c in consensus_features
                       if c["classification"] == "mandatory")
    n_optional = len(consensus_features) - n_mandatory

    # family_summary: per-family mandatory/optional split
    family_summary = defaultdict(lambda: {"mandatory": 0, "optional": 0})
    for c in consensus_features:
        family_summary[c["family"]][c["classification"]] += 1
    family_summary = dict(family_summary)

    return {
        "method": (f"Kabsch alignment of {n_sources} orthosteric Cα + "
                    f"greedy single-linkage clustering ({cluster_factor}× tolerance)"),
        "reference_pdb": reference["pdb_id"],
        "n_sources": n_sources,
        "alignment": alignment_summary,
        "tolerances_A": {
            "Donor": 1.0, "Acceptor": 1.0, "Aromatic": 1.0,
            "LumpedHydrophobe": 1.5, "PosIonizable": 1.5, "NegIonizable": 1.5,
        },
        "family_summary": family_summary,
        "features": consensus_features,
        "query_summary": {
            "total_features": len(consensus_features),
            "mandatory": n_mandatory,
            "optional": n_optional,
            "screening_strategy": {
                "strict": f"Match all {n_mandatory} mandatory features",
                "relaxed": f"Match >={threshold} mandatory + >=3 optional",
                "scoring": (f"Score = mandatory_matches + 0.5 × optional_matches; "
                             f"threshold >= {0.75 * n_mandatory:.1f} "
                             f"(0.75 × {n_mandatory} mandatory)"),
            },
        },
    }


def write_consensus_json(consensus, out_path, receptor_name, uniprot_id,
                          source_entries, excluded_entries):
    """Write a consensus_pharmacophore.json in the canonical schema."""
    blob = {
        "description": f"{receptor_name} orthosteric pharmacophore consensus for virtual screening",
        "target": (f"{receptor_name} (UniProt {uniprot_id})"
                    if uniprot_id else receptor_name),
        "source_entries": source_entries,
        "excluded_entries": excluded_entries,
        "alignment": consensus["alignment"],
        "method": consensus["method"],
        "tolerances_A": consensus["tolerances_A"],
        "family_summary": consensus["family_summary"],
        "features": consensus["features"],
        "query_summary": consensus["query_summary"],
    }
    with open(out_path, "w") as f:
        json.dump(blob, f, indent=2)


# -------- Loading existing consensus (when not derived) --------

def _feature_classification(feat):
    """Read classification from either top-level or nested-conservation schema."""
    if "classification" in feat:
        return feat["classification"]
    cons = feat.get("conservation") or {}
    return cons.get("classification")


def summarize_existing_consensus(consensus_path):
    """Read a consensus_pharmacophore.json and return its summary stats.

    Robust to two schemas seen in the wild:
      - top-level: feature.{classification, n_sources, members}     (canonical)
      - nested:   feature.conservation.{classification, n_sources}   (legacy)
    """
    if not consensus_path or not os.path.exists(consensus_path):
        return None
    with open(consensus_path) as f:
        d = json.load(f)
    feats = d.get("features", [])
    n_mand = sum(1 for x in feats if _feature_classification(x) == "mandatory")
    n_opt = sum(1 for x in feats if _feature_classification(x) == "optional")
    by_family = Counter(x["family"] for x in feats)
    by_family_mand = Counter(x["family"] for x in feats
                              if _feature_classification(x) == "mandatory")
    return {
        "path": consensus_path,
        "n_total": len(feats), "n_mandatory": n_mand, "n_optional": n_opt,
        "by_family": dict(by_family),
        "by_family_mandatory": dict(by_family_mand),
        "method": d.get("method"),
        "source_pdbs": d.get("source_pdbs") or
                       [s.get("pdb_id") for s in d.get("source_entries", [])],
        "reference_frame": d.get("reference_frame"),
        "query_summary": d.get("query_summary"),
        "source_entries": d.get("source_entries"),
        "excluded_entries": d.get("excluded_entries"),
    }


# -------- Main --------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-dir", required=True,
                    help="Root project directory (must contain per-PDB pharmacophore JSONs)")
    ap.add_argument("--receptor-name", required=True,
                    help="Human-readable receptor name (e.g. 'CXCR4', 'ATR')")
    ap.add_argument("--uniprot-id", default=None,
                    help="UniProt accession of the receptor (e.g. 'P61073'); optional")
    ap.add_argument("--output", default=None,
                    help="Output JSON path (default: <project-dir>/.report_data.json)")
    ap.add_argument("--retention-threshold", type=float, default=0.10,
                    help="Mean retention >= this classifies a ligand as orthosteric "
                         "(default 0.10; lower than CXCR4-style polycations because kinase "
                         "ATP-site binders engage fewer features per ligand)")
    ap.add_argument("--orthosteric-ligands", default="",
                    help="Comma-separated ligand resnames to force-classify as "
                         "orthosteric, overriding retention. Use when bond-order failures "
                         "depress a chemotype's retention (e.g. 'A1E' for ATR).")
    ap.add_argument("--excluded-ligands", default="",
                    help="Comma-separated ligand resnames to force-exclude.")
    ap.add_argument("--cluster-factor", type=float, default=2.0,
                    help="Clustering distance = factor × tolerance (default 2.0)")
    ap.add_argument("--derive-consensus", action="store_true",
                    help="Force derivation of consensus_pharmacophore.json even if one exists")
    ap.add_argument("--no-consensus", action="store_true",
                    help="Skip consensus derivation entirely")
    args = ap.parse_args()

    out_path = args.output or os.path.join(args.project_dir, ".report_data.json")

    print(f"[1/6] Scanning {args.project_dir} for pharmacophore JSONs ...",
          file=sys.stderr)
    entries = find_pharmacophore_jsons(args.project_dir)
    if not entries:
        sys.exit(f"No *.pharmacophore.json files found under {args.project_dir}")
    print(f"      Found {len(entries)} JSON files "
          f"({sum(1 for e in entries if e['variant'] == 'filtered')} filtered, "
          f"{sum(1 for e in entries if e['variant'] == 'all')} all)",
          file=sys.stderr)
    entries = load_jsons(entries)

    force_ortho = [s.strip().upper() for s in args.orthosteric_ligands.split(",") if s.strip()]
    force_excl = [s.strip().upper() for s in args.excluded_ligands.split(",") if s.strip()]
    print(f"[2/6] Classifying ligands (threshold {args.retention_threshold}"
          f"{', force_ortho=' + ','.join(force_ortho) if force_ortho else ''}"
          f"{', force_excl=' + ','.join(force_excl) if force_excl else ''}) ...",
          file=sys.stderr)
    classification = classify_ligands(entries, args.retention_threshold,
                                       force_orthosteric=force_ortho,
                                       force_excluded=force_excl)
    n_ortho = len(classification["orthosteric_pdbs"])
    n_excl = len(classification["excluded_pdbs"])
    print(f"      Orthosteric: {n_ortho} entries; excluded: {n_excl} entries",
          file=sys.stderr)
    for resname, b in classification["by_ligand"].items():
        print(f"        {resname}: {b['copies']} copies, "
              f"retention {b['mean_retention']:.0%} → {b['classification']}",
              file=sys.stderr)

    print(f"[3/6] Aggregating residue engagement across {n_ortho} "
          f"orthosteric entries ...", file=sys.stderr)
    residue_engagement = aggregate_residue_engagement(
        entries, classification["orthosteric_pdbs"])
    print(f"      {len(residue_engagement)} unique residues engaged",
          file=sys.stderr)

    print(f"[4/6] Computing per-family retention ...", file=sys.stderr)
    family_retention = compute_family_retention(
        entries, classification["orthosteric_pdbs"])

    print(f"[5/6] Detecting numbering offsets ...", file=sys.stderr)
    numbering_offsets = detect_numbering_offsets(residue_engagement)
    if numbering_offsets:
        for off in numbering_offsets:
            print(f"      {off['alias_a']} ({','.join(off['pdbs_a'])}) "
                  f"≈ {off['alias_b']} ({','.join(off['pdbs_b'])}) "
                  f"[offset {off['offset']}]", file=sys.stderr)

    file_manifest = build_file_manifest(entries, classification)

    # Consensus
    print(f"[6/6] Consensus pharmacophore ...", file=sys.stderr)
    existing_consensus_path = find_consensus_json(args.project_dir)
    consensus_summary = None
    derived_consensus = None

    if args.no_consensus:
        print(f"      Skipped (--no-consensus)", file=sys.stderr)
        consensus_summary = summarize_existing_consensus(existing_consensus_path)
    elif existing_consensus_path and not args.derive_consensus:
        print(f"      Using existing: {existing_consensus_path}",
              file=sys.stderr)
        consensus_summary = summarize_existing_consensus(existing_consensus_path)
    else:
        if n_ortho < 2:
            print(f"      Skipped (only {n_ortho} orthosteric entry; "
                  f"need >=2 for consensus)", file=sys.stderr)
        else:
            print(f"      Deriving (Kabsch + greedy clustering, "
                  f"factor {args.cluster_factor}) ...", file=sys.stderr)
            derived_consensus = derive_consensus(
                entries, classification, args.project_dir,
                cluster_factor=args.cluster_factor)
            if derived_consensus:
                # Write consensus_pharmacophore.json to the deepest common
                # ancestor of all per-PDB outputs — works for both flat
                # (pharmacophore/{pdb}_filtered.pharmacophore.json) and
                # per-PDB-subdir (pharmacophore/{pdb}/...) layouts.
                ortho_paths = [o["filtered_path"]
                               for o in classification["orthosteric_pdbs"]]
                if len(ortho_paths) == 1:
                    consensus_dir = os.path.dirname(ortho_paths[0])
                else:
                    consensus_dir = os.path.commonpath(
                        [os.path.dirname(p) for p in ortho_paths])
                consensus_out = os.path.join(consensus_dir,
                                             "consensus_pharmacophore.json")
                # Build canonical source_entries / excluded_entries
                source_entries = []
                seen_src = set()
                for o in classification["orthosteric_pdbs"]:
                    key = (o["pdb_id"], o.get("ligand"))
                    if key in seen_src:
                        continue
                    seen_src.add(key)
                    source_entries.append({
                        "pdb_id": o["pdb_id"],
                        "ligand": o["ligand"],
                        "n_filtered_features": o["filtered_count"],
                    })
                excluded_entries = []
                for resname, b in classification["by_ligand"].items():
                    if b["classification"] != "excluded":
                        continue
                    pdbs = sorted({m["pdb_id"] for m in b["members"]})
                    excluded_entries.append({
                        "pdb_ids": pdbs if len(pdbs) > 1 else None,
                        "pdb_id": pdbs[0] if len(pdbs) == 1 else None,
                        "ligand": resname,
                        "reason": (f"Mean retention {b['mean_retention']:.0%} "
                                    f"across {b['copies']} copies; "
                                    f"non-orthosteric"),
                    })
                # strip null fields
                excluded_entries = [{k: v for k, v in e.items() if v is not None}
                                     for e in excluded_entries]
                write_consensus_json(derived_consensus, consensus_out,
                                      args.receptor_name, args.uniprot_id,
                                      source_entries, excluded_entries)
                print(f"      Wrote {consensus_out} "
                      f"({derived_consensus['query_summary']['total_features']} features, "
                      f"{derived_consensus['query_summary']['mandatory']} mandatory)",
                      file=sys.stderr)
                consensus_summary = summarize_existing_consensus(consensus_out)

    # Strip non-serializable bits before dumping
    for ligand in classification["by_ligand"].values():
        for m in ligand.get("members", []):
            m.pop("filtered_path", None)
            m.pop("all_path", None)

    report_data = {
        "receptor": {
            "name": args.receptor_name,
            "uniprot_id": args.uniprot_id,
        },
        "project_dir": os.path.abspath(args.project_dir),
        "n_pdb_entries": len({e["pdb_id"] for e in entries}),
        "ligand_classification": classification,
        "residue_engagement": residue_engagement,
        "feature_retention": family_retention,
        "numbering_offsets": numbering_offsets,
        "file_manifest": file_manifest,
        "consensus_summary": consensus_summary,
        "consensus_derived_this_run": derived_consensus is not None,
    }

    # Strip the per-entry data['data'] from file_manifest (keep only summary)
    with open(out_path, "w") as f:
        json.dump(report_data, f, indent=2, default=lambda o: list(o)
                   if isinstance(o, set) else str(o))
    print(f"\n[done] Report data written to: {out_path}", file=sys.stderr)
    print(f"       File size: {os.path.getsize(out_path) / 1024:.1f} KB",
          file=sys.stderr)


if __name__ == "__main__":
    main()
