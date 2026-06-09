#!/usr/bin/env python3
"""
extract_pharmacophore_data.py — target-agnostic extractor for the pharmacophore-dashboard skill.

Discovers structure-based pharmacophore analysis outputs (RDKit *_filtered/_all
pharmacophore JSON + cleaned PDB receptors), extracts ligands + binding pockets for
3D rendering, computes retention and residue-engagement statistics, and writes a single
JSON data bundle that `build_dashboard.py` turns into a self-contained interactive HTML report.

Input contract (per analyzed ligand copy) — produced by a structure-based pharmacophore
detector such as RDKit BaseFeatures.fdef with an interaction filter:

  <prefix>_filtered.pharmacophore.json   # retained features (have a complementary protein partner)
  <prefix>_all.pharmacophore.json        # all detected features (optional, for retention stats)

Each JSON has: {input_pdb, ligand:{resname,chain,resnum,smiles_used,smiles_source},
  features:[{id,family,type,position:[x,y,z],tolerance,ligand_atom_ids,
             interaction:{partner_residue,partner_chain,partner_atom,distance_A,type}}],
  metadata:{n_features_before_filter,n_features_after_filter}}

An optional consensus query JSON (e.g. *consensus*.json) is auto-detected and embedded.

Usage:
  python3 extract_pharmacophore_data.py PHARM_DIR [--raw-dir DIR] [--out data.json]
      [--consensus FILE] [--config cfg.json] [--pocket-cutoff 5.0]
      [--target-name NAME] [--gene G] [--uniprot U] [--organism O]
"""
import argparse, glob, json, os, re, sys
import numpy as np

FAMILIES = ["Donor", "Acceptor", "PosIonizable", "NegIonizable", "Aromatic", "LumpedHydrophobe"]
# Categories excluded from retention / residue aggregation (negative controls, non-binders).
EXCLUDED_CATS = {"excluded", "control", "negative", "negative control", "inactive", "decoy"}


# ---------------- PDB parsing (fixed columns) ----------------
def parse_pdb(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            name = line[12:16].strip()
            elem = line[76:78].strip()
            if not elem:
                m = re.match(r"\s*([A-Za-z]{1,2})", name)
                elem = (m.group(1)[0] if m else "C")
            atoms.append({"rec": rec, "name": name, "altloc": line[16:17].strip(),
                          "resname": line[17:20].strip(), "chain": line[21:22].strip(),
                          "resseq": line[22:26].strip(), "x": x, "y": y, "z": z,
                          "elem": elem.upper()[:2].strip().capitalize()})
    return atoms


def fmt_pdb_line(i, a):
    name = a["name"]
    nm = ((" " + name) if len(name) < 4 else name).ljust(4)[:4]
    return (f"{a['rec'].ljust(6)}{i:5d} {nm}{a['altloc'][:1] or ' '}{a['resname'][:3]:>3} "
            f"{(a['chain'] or 'A')[:1]}{a['resseq'][:4]:>4}    "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00          {a['elem'][:2]:>2}")


def extract_pocket(atoms, lig_resname, lig_chain, lig_resnum, cutoff=5.0):
    """Return (pdb_text, lig_heavy_atoms, pocket_residue_count). Hydrogens dropped."""
    lig = [a for a in atoms if a["rec"] == "HETATM" and a["resname"] == lig_resname
           and a["chain"] == lig_chain and a["resseq"] == str(lig_resnum)]
    lig_heavy = [a for a in lig if a["elem"] != "H"]
    if not lig_heavy:
        return None, [], 0
    lig_xyz = np.array([[a["x"], a["y"], a["z"]] for a in lig_heavy])
    prot = [a for a in atoms if a["rec"] == "ATOM" and a["elem"] != "H"]
    pocket = []
    if prot:
        prot_xyz = np.array([[a["x"], a["y"], a["z"]] for a in prot])
        d = np.sqrt(((prot_xyz[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(-1)).min(1)
        keep = {(a["chain"], a["resseq"], a["resname"]) for a, dist in zip(prot, d) if dist <= cutoff}
        pocket = [a for a in prot if (a["chain"], a["resseq"], a["resname"]) in keep]
    lines, i = [], 1
    for a in lig_heavy:
        b = dict(a); b["rec"] = "HETATM"; lines.append(fmt_pdb_line(i, b)); i += 1
    for a in pocket:
        b = dict(a); b["rec"] = "ATOM"; lines.append(fmt_pdb_line(i, b)); i += 1
    lines.append("END")
    n_res = len({(a["chain"], a["resseq"]) for a in pocket})
    return "\n".join(lines), lig_heavy, n_res


def find_partner_atom(atoms, resname, resnum, chain, atomname):
    for a in atoms:
        if (a["resname"] == resname and a["resseq"] == str(resnum)
                and a["chain"] == chain and a["name"] == atomname):
            return [round(a["x"], 3), round(a["y"], 3), round(a["z"], 3)]
    for a in atoms:  # fallback: ignore chain
        if a["resname"] == resname and a["resseq"] == str(resnum) and a["name"] == atomname:
            return [round(a["x"], 3), round(a["y"], 3), round(a["z"], 3)]
    return None


def split_residue(token):
    m = re.match(r"([A-Za-z]+)(\d+)", token or "")
    return (m.group(1), m.group(2)) if m else (token, "")


def feature_records(jdata, atoms):
    feats = []
    for f in jdata.get("features", []):
        inter = f.get("interaction") or {}
        pres = inter.get("partner_residue")
        rn, rs = split_residue(pres)
        xyz = find_partner_atom(atoms, rn, rs, inter.get("partner_chain", ""),
                                inter.get("partner_atom", "")) if pres else None
        feats.append({"id": f.get("id"), "family": f.get("family"), "type": f.get("type"),
                      "position": [round(c, 3) for c in f.get("position", [])],
                      "tolerance": f.get("tolerance"),
                      "partner_residue": pres, "partner_resname": rn, "partner_resnum": rs,
                      "partner_chain": inter.get("partner_chain"), "partner_atom": inter.get("partner_atom"),
                      "distance": inter.get("distance_A"), "itype": inter.get("type"),
                      "partner_xyz": xyz})
    return feats


def family_counts(jdata):
    c = {fam: 0 for fam in FAMILIES}
    for f in jdata.get("features", []):
        c[f.get("family")] = c.get(f.get("family"), 0) + 1
    return c


def pdb_id_from_name(fname):
    """Leading token of the filename, e.g. 9L45_filtered... -> 9L45, 8U4P_chainR... -> 8U4P."""
    base = os.path.basename(fname)
    return re.split(r"[._]", base)[0]


def _readable(p):
    """True only for a real, openable file (excludes broken symlinks)."""
    return bool(p) and os.path.isfile(p) and os.access(p, os.R_OK)


def resolve_clean_pdb(input_pdb, bases):
    """Try the JSON's input_pdb against several base dirs, then search by basename.
    Skips broken symlinks (common when a pharmacophore dir links to ../raw/...)."""
    cands = []
    if input_pdb:
        if os.path.isabs(input_pdb):
            cands.append(input_pdb)
        for b in bases:
            cands.append(os.path.normpath(os.path.join(b, input_pdb)))
    for c in cands:
        if _readable(c):
            return c
    # fallback: search for the basename under the bases; prefer real files under a raw/ dir
    if input_pdb:
        bn = os.path.basename(input_pdb)
        hits = []
        for b in bases:
            hits += [h for h in glob.glob(os.path.join(b, "**", bn), recursive=True) if _readable(h)]
        if hits:
            hits.sort(key=lambda h: (0 if (os.sep + "raw" + os.sep) in h else 1, len(h)))
            return hits[0]
    return None


def discover_entries(pharm_dir):
    """Find filtered JSONs, dedupe by (pdb_id,resname,chain,resnum), prefer shallower paths."""
    files = glob.glob(os.path.join(pharm_dir, "**", "*_filtered.pharmacophore.json"), recursive=True)
    files += glob.glob(os.path.join(pharm_dir, "**", "*.filtered.pharmacophore.json"), recursive=True)
    files = sorted(set(files), key=lambda p: (p.count(os.sep), len(p), p))
    seen, chosen = {}, []
    for fp in files:
        try:
            j = json.load(open(fp))
        except Exception:
            continue
        lig = j.get("ligand", {})
        key = (pdb_id_from_name(fp), lig.get("resname"), lig.get("chain"), lig.get("resnum"))
        if key in seen:
            continue
        seen[key] = True
        chosen.append((fp, j))
    return chosen


def classify(n_after, override):
    if override:
        return override
    return "productive" if (n_after or 0) > 0 else "control"


def build_bundle(pharm_dir, raw_dir=None, consensus_path=None, cfg=None,
                 pocket_cutoff=5.0, target_meta=None):
    cfg = cfg or {}
    pharm_dir = os.path.abspath(pharm_dir)
    parent = os.path.dirname(pharm_dir)
    bases = [b for b in [raw_dir, parent, pharm_dir, os.getcwd()] if b]
    overrides = {o.get("key") or o.get("pdb_id"): o for o in cfg.get("entry_overrides", [])}

    entries, retention = [], {f: {"all": 0, "filtered": 0} for f in FAMILIES}
    res_eng, pdb_ids, missing = {}, set(), []

    for fp, fild in discover_entries(pharm_dir):
        pid = pdb_id_from_name(fp); pdb_ids.add(pid)
        ap = re.sub(r"_filtered\.pharmacophore\.json$|\.filtered\.pharmacophore\.json$",
                    "_all.pharmacophore.json", fp)
        alld = json.load(open(ap)) if os.path.isfile(ap) else {"features": []}
        lig = fild.get("ligand", {})
        clean = resolve_clean_pdb(fild.get("input_pdb"), bases)
        if not clean:
            # last resort: any readable *_clean.pdb / *.pdb mentioning pid
            for b in bases:
                hits = [h for h in glob.glob(os.path.join(b, "**", f"{pid}*_clean.pdb"), recursive=True) if _readable(h)] \
                       or [h for h in glob.glob(os.path.join(b, "**", f"{pid}*.pdb"), recursive=True) if _readable(h)]
                if hits:
                    hits.sort(key=lambda h: (0 if (os.sep + "raw" + os.sep) in h else 1, len(h)))
                    clean = hits[0]; break
        atoms = []
        if clean:
            try:
                atoms = parse_pdb(clean)
            except Exception as ex:
                print(f"[extract] WARNING could not parse {clean}: {ex}", file=sys.stderr)
        if not atoms:
            missing.append(pid)
        pocket_pdb, _, n_pocket = (extract_pocket(atoms, lig.get("resname"), lig.get("chain"),
                                                   lig.get("resnum"), pocket_cutoff)
                                   if atoms else (None, [], 0))
        feats = feature_records(fild, atoms)
        fc_all, fc_fil = family_counts(alld), family_counts(fild)
        meta = fild.get("metadata", {})
        n_before = meta.get("n_features_before_filter", sum(fc_all.values()) or sum(fc_fil.values()))
        n_after = meta.get("n_features_after_filter", sum(fc_fil.values()))
        key = re.sub(r"\.pharmacophore\.json$", "", os.path.basename(fp))
        ov = overrides.get(key) or overrides.get(pid) or {}
        cat = classify(n_after, ov.get("category"))
        chain = lig.get("chain", "")
        label = ov.get("label") or (f"{pid}" + (f" · chain {chain}" if chain else ""))
        if cat not in EXCLUDED_CATS:
            for f in FAMILIES:
                retention[f]["all"] += fc_all.get(f, 0)
                retention[f]["filtered"] += fc_fil.get(f, 0)
            for f in feats:
                if f["partner_residue"]:
                    r = res_eng.setdefault(f["partner_residue"], {"residue": f["partner_residue"],
                        **{fam: 0 for fam in FAMILIES}, "total": 0, "structures": set()})
                    r[f["family"]] = r.get(f["family"], 0) + 1
                    r["total"] += 1; r["structures"].add(pid)
        entries.append({"key": key, "label": label, "pdb_id": pid,
                        "chemotype": ov.get("chemotype") or lig.get("resname"),
                        "chemotype_label": ov.get("chemotype_label"),
                        "role": ov.get("role", ""), "category": cat, "note": ov.get("note", ""),
                        "ligand": lig, "n_before": n_before, "n_after": n_after,
                        "family_all": fc_all, "family_filtered": fc_fil,
                        "features": feats, "pocket_pdb": pocket_pdb, "n_pocket_res": n_pocket})

    entries.sort(key=lambda e: (e["category"] in EXCLUDED_CATS, -e["n_after"], e["label"]))
    res_list = sorted(res_eng.values(), key=lambda r: -r["total"])
    for r in res_list:
        r["structures"] = sorted(r["structures"])

    bundle = {"target": target_meta or {}, "entries": entries, "retention": retention,
              "residue_engagement": res_list,
              "provenance": {"pharmacophore_dir": pharm_dir, "n_entries": len(entries),
                             "n_structures": len(pdb_ids), "missing_pdb_for": missing,
                             "pocket_cutoff_A": pocket_cutoff}}
    bundle["target"].setdefault("n_structures", len(pdb_ids))

    # ---- optional consensus query ----
    if consensus_path is None:
        cands = [p for p in glob.glob(os.path.join(pharm_dir, "*.json"))
                 if "consensus" in os.path.basename(p).lower()]
        consensus_path = cands[0] if cands else None
    if consensus_path and os.path.isfile(consensus_path):
        bundle["consensus"] = build_consensus(consensus_path, entries, bases, pocket_cutoff)
        if not bundle["target"] and bundle["consensus"].get("target"):
            bundle["target"].update(bundle["consensus"]["target"])
    return bundle


def build_consensus(path, entries, bases, pocket_cutoff):
    cons = json.load(open(path))
    feats = []
    for f in cons.get("features", []):
        cls = f.get("classification")
        if cls is None and (f.get("required") or f.get("mandatory")):
            cls = "mandatory"
        elif cls is None:
            cls = "optional"
        feats.append({"id": f.get("id"), "family": f.get("family"),
                      "position": [round(c, 3) for c in f.get("position", [])],
                      "tolerance": f.get("tolerance"), "classification": cls,
                      "n_sources": f.get("n_sources", f.get("sources"))})
    # choose a reference structure to overlay the consensus pocket on
    ref_txt = ((cons.get("alignment") or {}).get("reference") or "")
    ref_pid = None
    for e in entries:
        if e["pdb_id"] and e["pdb_id"] in ref_txt:
            ref_pid = e["pdb_id"]; break
    ref_entry = next((e for e in entries if e["pdb_id"] == ref_pid), None) \
        or max((e for e in entries if e["pocket_pdb"]), key=lambda e: e["n_after"], default=None)
    pocket = ref_entry["pocket_pdb"] if ref_entry else None
    summary = cons.get("query_summary") or {}
    if "total_features" not in summary:
        summary = {"total_features": len(feats),
                   "mandatory": sum(1 for f in feats if f["classification"] == "mandatory"),
                   "optional": sum(1 for f in feats if f["classification"] == "optional")}
    return {"summary": summary, "alignment": cons.get("alignment"),
            "source_entries": cons.get("source_entries"), "excluded_entries": cons.get("excluded_entries"),
            "target": cons.get("target"), "features": feats, "pocket_pdb": pocket,
            "frame_pdb": ref_entry["pdb_id"] if ref_entry else None}


def main():
    ap = argparse.ArgumentParser(description="Extract pharmacophore data bundle for the dashboard skill.")
    ap.add_argument("pharm_dir", help="Directory containing *_filtered.pharmacophore.json files")
    ap.add_argument("--raw-dir", help="Base dir to resolve each JSON's input_pdb (default: parent of pharm_dir)")
    ap.add_argument("--consensus", help="Path to a consensus query JSON (default: auto-detect *consensus*.json)")
    ap.add_argument("--config", help="Optional report_config.json (entry_overrides, target, narrative)")
    ap.add_argument("--out", default="pharmviz_data.json", help="Output data bundle path")
    ap.add_argument("--pocket-cutoff", type=float, default=5.0, help="Pocket residue cutoff in Angstrom")
    ap.add_argument("--target-name"); ap.add_argument("--gene")
    ap.add_argument("--uniprot"); ap.add_argument("--organism")
    a = ap.parse_args()
    cfg = json.load(open(a.config)) if a.config and os.path.isfile(a.config) else {}
    tmeta = dict(cfg.get("target", {}))
    for k, v in [("name", a.target_name), ("gene", a.gene), ("uniprot", a.uniprot), ("organism", a.organism)]:
        if v:
            tmeta[k] = v
    if not tmeta.get("name"):
        tmeta["name"] = os.path.basename(os.path.dirname(os.path.abspath(a.pharm_dir))) or "Target"
    bundle = build_bundle(a.pharm_dir, a.raw_dir, a.consensus, cfg, a.pocket_cutoff, tmeta)
    json.dump(bundle, open(a.out, "w"))
    p = bundle["provenance"]
    print(f"[extract] {p['n_entries']} entries from {p['n_structures']} structures -> {a.out} "
          f"({round(os.path.getsize(a.out)/1024)} KB)")
    if p["missing_pdb_for"]:
        print(f"[extract] WARNING: no clean PDB found for: {', '.join(p['missing_pdb_for'])} "
              f"(3D viewer disabled for those entries)")
    if "consensus" in bundle:
        s = bundle["consensus"]["summary"]
        print(f"[extract] consensus query: {s.get('total_features')} features "
              f"(overlay frame: {bundle['consensus'].get('frame_pdb')})")


if __name__ == "__main__":
    main()
