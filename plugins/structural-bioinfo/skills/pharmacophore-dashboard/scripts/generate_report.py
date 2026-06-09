#!/usr/bin/env python3
"""
generate_report.py — one-shot: pharmacophore analysis directory -> interactive HTML dashboard.

Chains extract_pharmacophore_data.py + build_dashboard.py. Use this for the standard case.

Examples:
  # Minimal — auto-discovers everything under the directory:
  python3 generate_report.py /path/to/target/pharmacophore --out TARGET_report.html

  # With receptor PDBs in a sibling 'raw' dir and target metadata:
  python3 generate_report.py atr/pharmacophore --raw-dir atr \
      --target-name "ATR kinase" --gene ATR --uniprot Q13535 --accent "#4f46e5" \
      --out ATR_pharmacophore_report.html

  # With a custom narrative/branding config:
  python3 generate_report.py cxcr4/pharmacophore --config cxcr4_report_config.json \
      --out CXCR4_pharmacophore_report.html

Requires: numpy. The generated HTML loads 3Dmol.js from a CDN at view time (internet needed
for the 3D viewers only; charts/tables/data are fully offline).
"""
import argparse, json, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import extract_pharmacophore_data as extract  # noqa: E402
import build_dashboard as build  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Pharmacophore analysis dir -> interactive HTML dashboard.")
    ap.add_argument("pharm_dir", help="Directory with *_filtered.pharmacophore.json files")
    ap.add_argument("--raw-dir", help="Base dir to resolve each JSON's input_pdb (default: parent of pharm_dir)")
    ap.add_argument("--consensus", help="Consensus query JSON (default: auto-detect *consensus*.json)")
    ap.add_argument("--config", help="report_config.json (entry_overrides, target, narrative, branding)")
    ap.add_argument("--out", default="pharmacophore_report.html")
    ap.add_argument("--data-out", help="Also keep the intermediate data bundle JSON at this path")
    ap.add_argument("--pocket-cutoff", type=float, default=5.0)
    ap.add_argument("--title")
    ap.add_argument("--accent")
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

    bundle = extract.build_bundle(a.pharm_dir, a.raw_dir, a.consensus, cfg, a.pocket_cutoff, tmeta)
    p = bundle["provenance"]
    n_prod = sum(1 for e in bundle["entries"]
                 if e["n_after"] > 0 and e["category"] not in extract.EXCLUDED_CATS)
    print(f"[generate] {p['n_entries']} entries from {p['n_structures']} structures ({n_prod} productive)")
    if p["missing_pdb_for"]:
        print(f"[generate] WARNING no clean PDB for {', '.join(p['missing_pdb_for'])} — "
              f"pass --raw-dir to enable their 3D viewers")
    if "consensus" in bundle:
        print(f"[generate] consensus query: {bundle['consensus']['summary'].get('total_features')} features")

    if a.data_out:
        json.dump(bundle, open(a.data_out, "w"))

    if a.accent:
        cfg["accent"] = a.accent
    title = a.title or f"{tmeta['name']} — Pharmacophore Analysis"
    html = build.build(bundle, cfg, title)
    with open(a.out, "w") as fh:
        fh.write(html)
    print(f"[generate] wrote {a.out}  ({round(len(html)/1024)} KB)")


if __name__ == "__main__":
    main()
