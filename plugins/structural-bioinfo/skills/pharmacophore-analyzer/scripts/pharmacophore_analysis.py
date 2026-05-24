#!/usr/bin/env python3
"""
Compute a 3D pharmacophore from a protein-ligand complex PDB file using RDKit.

Identifies the bound ligand, restores correct bond orders and hydrogens via a
SMILES template (auto-fetched from the RCSB Chemical Component Dictionary by
default), runs the BaseFeatures.fdef feature factory on the ligand, and (by
default) keeps only those features that engage the protein within type-specific
distance cutoffs (LigandScout-style structure-based pharmacophore).

Outputs:
    <prefix>.pharmacophore.json   feature list with family/type/position/tolerance
    <prefix>.pharmacophore.pml    PyMOL scene: ligand sticks + colored spheres

Usage:
    # Auto-detect ligand, auto-fetch SMILES from RCSB
    python pharmacophore_analysis.py 1IEP_clean.pdb --output-prefix 1IEP

    # Explicit ligand + user-supplied SMILES (fully offline)
    python pharmacophore_analysis.py complex.pdb \\
        --ligand-resname STI --smiles "Cc1ccc(...)c1" -o imatinib

    # Skip the interaction filter (every ligand feature)
    python pharmacophore_analysis.py 1IEP_clean.pdb --include-all-features

    # Offline best-effort (no SMILES, accuracy degrades)
    python pharmacophore_analysis.py complex.pdb --no-fetch
"""

import argparse
import json
import math
import os
import sys
import urllib.request
import urllib.error

# ── RDKit bootstrap ──────────────────────────────────────────────────
# RDKit is rarely on the system Python; conda envs are the norm. If we can't
# import it here, search common conda locations for an env that has it and
# re-exec the script under that env's Python. A loop guard prevents infinite
# re-exec when no rdkit-bearing env is installed anywhere.

def _candidate_pythons():
    """Yield candidate Python executables that might have rdkit, ordered by
    preference: env literally named 'rdkit' first, then any sibling env, then
    conda base installs. Roots searched: the user's home conda installs and a
    few common system locations.
    """
    import glob
    home = os.path.expanduser("~")
    roots = [
        f"{home}/miniconda3", f"{home}/anaconda3",
        f"{home}/miniforge3", f"{home}/mambaforge", f"{home}/conda",
        "/opt/conda", "/opt/miniconda3", "/opt/anaconda3",
        "/opt/homebrew/Caskroom/miniconda/base",
    ]
    cp = os.environ.get("CONDA_PREFIX", "")
    if cp:
        # If we're inside .../envs/foo, walk up to the conda root.
        if os.path.basename(os.path.dirname(cp)) == "envs":
            roots.insert(0, os.path.dirname(os.path.dirname(cp)))
        else:
            roots.insert(0, cp)
    seen = set()
    # 1. Prefer an env literally named "rdkit"
    for root in roots:
        py = os.path.join(root, "envs", "rdkit", "bin", "python")
        if py not in seen and os.path.isfile(py):
            seen.add(py); yield py
    # 2. Then any sibling env
    for root in roots:
        for env in sorted(glob.glob(os.path.join(root, "envs", "*"))):
            py = os.path.join(env, "bin", "python")
            if py not in seen and os.path.isfile(py):
                seen.add(py); yield py
    # 3. Then conda base installs themselves
    for root in roots:
        py = os.path.join(root, "bin", "python")
        if py not in seen and os.path.isfile(py):
            seen.add(py); yield py


def _bootstrap_rdkit():
    """Ensure rdkit is importable, re-execing under another Python if needed."""
    try:
        import rdkit  # noqa: F401
        return
    except ImportError:
        pass
    if os.environ.get("PHARM_RDKIT_REEXECED"):
        sys.stderr.write(
            "[error] RDKit is not importable even after auto-discovery.\n"
            "        Install with: conda install -c conda-forge rdkit -n rdkit\n"
            "                  or: pip install rdkit\n"
        )
        sys.exit(2)
    import subprocess
    for py in _candidate_pythons():
        try:
            r = subprocess.run([py, "-c", "import rdkit"],
                               capture_output=True, timeout=10)
        except (subprocess.TimeoutExpired, OSError):
            continue
        if r.returncode == 0:
            sys.stderr.write(
                f"[info] Re-executing under {py} (auto-discovered RDKit env)\n"
            )
            os.environ["PHARM_RDKIT_REEXECED"] = "1"
            os.execv(py, [py] + sys.argv)
    sys.stderr.write(
        "[error] RDKit is not importable and no conda env with rdkit was found.\n"
        "        Searched ~/miniconda3, ~/anaconda3, ~/miniforge3, ~/mambaforge,\n"
        "        /opt/conda, /opt/homebrew/Caskroom/miniconda/base.\n"
        "        Install with: conda install -c conda-forge rdkit -n rdkit\n"
        "                  or: pip install rdkit\n"
    )
    sys.exit(2)


_bootstrap_rdkit()

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures
from rdkit import __version__ as RDKIT_VERSION


# ── Constants ─────────────────────────────────────────────────────────

# Residues to ignore when auto-detecting the ligand (water, ions, buffers,
# crystallization additives, common membrane mimetics). Mirrors pdb-extractor.
EXCLUDE_RESIDUES = {
    "HOH", "DOD", "WAT",
    "NA", "K", "CL", "MG", "CA", "ZN", "FE", "MN", "CU", "CO",
    "CD", "NI", "BR", "I", "F", "LI", "RB", "CS", "SR", "BA",
    "PT", "AU", "HG", "SE",
    "EDO", "GOL", "TRS", "PEG", "PGE", "MPD", "ACT", "CIT",
    "BME", "DMS", "SO4", "PO4", "FMT", "EOH", "IMD", "LDA",
    "BCT", "BOG", "DMU", "OLC", "PLM", "SCN", "NO3", "NH4",
    "TAM", "HEP", "MES", "PIP", "BIS", "TRI", "BTB",
    "SUC", "DTT", "TCE", "DTE", "SGM",
    "UNK", "UNX", "ACE", "NH2", "FOR", "DUM",
    "LMT", "OLA", "SDS", "LMU", "LMN", "TGL", "D10",
    "DMX", "OCT", "HEX", "LHG",
}

# Per-family tolerance radii (Å), LigandScout/MOE convention.
# Directional features get tight radii; hydrophobic and ionizable get larger.
TOL_DEFAULTS = {
    "Donor": 1.0,
    "Acceptor": 1.0,
    "Aromatic": 1.0,
    "LumpedHydrophobe": 1.5,
    "Hydrophobe": 1.5,
    "PosIonizable": 1.5,
    "NegIonizable": 1.5,
    "ZnBinder": 1.0,
}

# PyMOL color names per family.
PML_COLORS = {
    "Donor": "blue",
    "Acceptor": "red",
    "Aromatic": "orange",
    "LumpedHydrophobe": "forest",
    "Hydrophobe": "forest",
    "PosIonizable": "marine",
    "NegIonizable": "firebrick",
    "ZnBinder": "purple",
}

# Short labels shown on the PML pseudoatom.
PML_LABELS = {
    "Donor": "HBD",
    "Acceptor": "HBA",
    "Aromatic": "ARO",
    "LumpedHydrophobe": "HYD",
    "Hydrophobe": "HYD",
    "PosIonizable": "POS",
    "NegIonizable": "NEG",
    "ZnBinder": "ZN",
}

# Per-atom features for protein side. Backbone N is HBD, backbone O is HBA.
# Side-chain mappings are conservative (one or two atoms per residue).
PROT_ATOM_FEATURES = {
    # Side-chain hydroxyls (both donor and acceptor)
    ("SER", "OG"):  ["Donor", "Acceptor"],
    ("THR", "OG1"): ["Donor", "Acceptor"],
    ("TYR", "OH"):  ["Donor", "Acceptor"],
    # Side-chain amides
    ("ASN", "OD1"): ["Acceptor"],
    ("ASN", "ND2"): ["Donor"],
    ("GLN", "OE1"): ["Acceptor"],
    ("GLN", "NE2"): ["Donor"],
    # Side-chain carboxylates
    ("ASP", "OD1"): ["Acceptor", "NegIonizable"],
    ("ASP", "OD2"): ["Acceptor", "NegIonizable"],
    ("GLU", "OE1"): ["Acceptor", "NegIonizable"],
    ("GLU", "OE2"): ["Acceptor", "NegIonizable"],
    # Lysine ammonium
    ("LYS", "NZ"):  ["Donor", "PosIonizable"],
    # Arginine guanidinium
    ("ARG", "NH1"): ["Donor", "PosIonizable"],
    ("ARG", "NH2"): ["Donor", "PosIonizable"],
    ("ARG", "NE"):  ["Donor", "PosIonizable"],
    # Histidine imidazole
    ("HIS", "ND1"): ["Donor", "Acceptor", "PosIonizable"],
    ("HIS", "NE2"): ["Donor", "Acceptor", "PosIonizable"],
    # Tryptophan indole NH
    ("TRP", "NE1"): ["Donor"],
    # Methionine sulfur (weak HBA / hydrophobic)
    ("MET", "SD"):  ["Hydrophobe"],
    # Hydrophobic aliphatic side chains (representative atom)
    ("ALA", "CB"):  ["Hydrophobe"],
    ("VAL", "CB"):  ["Hydrophobe"],
    ("LEU", "CG"):  ["Hydrophobe"],
    ("ILE", "CG1"): ["Hydrophobe"],
    ("PRO", "CG"):  ["Hydrophobe"],
    ("CYS", "SG"):  ["Hydrophobe"],
}

# Backbone atoms get HBD/HBA roles regardless of residue.
BACKBONE_FEATURES = {
    "N": ["Donor"],
    "O": ["Acceptor"],
}

# Aromatic rings to compute centroids for. (Phe/Tyr/Trp/His).
RING_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "TRP": ["CG", "CD1", "NE1", "CE2", "CD2"],  # 5-mem ring; 6-mem also aromatic
}

# Complementary-feature pairs for the interaction filter and the cutoff in Å.
# Each ligand family maps to a list of (protein_family, cutoff_Å) options.
COMPLEMENTS_DEFAULT = {
    "Donor":            [("Acceptor", 3.5)],
    "Acceptor":         [("Donor", 3.5)],
    "PosIonizable":     [("NegIonizable", 5.0)],
    "NegIonizable":     [("PosIonizable", 5.0)],
    "Aromatic":         [("Aromatic", 5.5), ("PosIonizable", 5.5)],
    "LumpedHydrophobe": [("Hydrophobe", 4.5)],
    "Hydrophobe":       [("Hydrophobe", 4.5)],
    "ZnBinder":         [("Zn", 3.0)],
}

# Map family pair → human-readable interaction type (for JSON output).
INTERACTION_LABEL = {
    ("Donor", "Acceptor"): "hbond",
    ("Acceptor", "Donor"): "hbond",
    ("PosIonizable", "NegIonizable"): "salt_bridge",
    ("NegIonizable", "PosIonizable"): "salt_bridge",
    ("Aromatic", "Aromatic"): "pi_stack",
    ("Aromatic", "PosIonizable"): "pi_cation",
    ("LumpedHydrophobe", "Hydrophobe"): "hydrophobic",
    ("Hydrophobe", "Hydrophobe"): "hydrophobic",
    ("ZnBinder", "Zn"): "metal",
}


# ── Helpers ───────────────────────────────────────────────────────────

def log(msg):
    print(msg, file=sys.stderr)


def fetch_smiles_from_rcsb(resname, timeout=10):
    """Pull the idealized SDF for a 3-letter ligand from the RCSB CCD,
    return its canonical SMILES. Returns None on failure.
    """
    url = f"https://files.rcsb.org/ligands/download/{resname}_ideal.sdf"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            sdf_bytes = resp.read()
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        log(f"    [warn] RCSB CCD fetch failed for {resname}: {e}")
        return None
    if not sdf_bytes:
        return None
    template = Chem.MolFromMolBlock(sdf_bytes.decode("utf-8", errors="replace"),
                                    removeHs=False)
    if template is None:
        log(f"    [warn] RCSB CCD SDF for {resname} could not be parsed")
        return None
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(template))
    except Exception as e:
        log(f"    [warn] Could not canonicalize SMILES for {resname}: {e}")
        return None


def get_residue_info(frag):
    """Return the (resname, chain, resnum) of the first heavy atom of a frag,
    or (None, None, None) if PDB residue info is missing.
    """
    if frag.GetNumAtoms() == 0:
        return None, None, None
    info = frag.GetAtomWithIdx(0).GetPDBResidueInfo()
    if info is None:
        return None, None, None
    resname = (info.GetResidueName() or "").strip()
    chain = (info.GetChainId() or "").strip()
    resnum = info.GetResidueNumber()
    return resname, chain, resnum


def is_candidate_ligand(frag, want_resname=None, want_chain=None):
    """Decide whether a fragment is a small-molecule ligand worth analyzing."""
    if frag.GetNumHeavyAtoms() < 6:
        return False
    info = frag.GetAtomWithIdx(0).GetPDBResidueInfo()
    if info is None:
        return False
    resname = (info.GetResidueName() or "").strip().upper()
    chain = (info.GetChainId() or "").strip()
    if want_resname and resname != want_resname.upper():
        return False
    if want_chain and chain != want_chain.upper():
        return False
    if want_resname is None:
        # Auto mode: must be HETATM and not a known artifact.
        if not info.GetIsHeteroAtom():
            return False
        if resname in EXCLUDE_RESIDUES:
            return False
    return True


def identify_ligands(mol, want_resname=None, want_chain=None):
    """Split mol into fragments and return those matching ligand criteria."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    return [f for f in frags
            if is_candidate_ligand(f, want_resname, want_chain)]


def prepare_ligand(frag, smiles=None, sdf_path=None, no_fetch=False):
    """Restore bond orders + hydrogens. Returns (mol, smiles_used, source).

    Priority for the SMILES template:
        1. --smiles
        2. --template-sdf
        3. RCSB CCD fetch (unless --no-fetch)
        4. best-effort SanitizeMol on the raw PDB fragment
    """
    resname, _, _ = get_residue_info(frag)

    template = None
    smiles_used = None
    source = None

    if smiles:
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            log(f"    [error] Could not parse --smiles '{smiles}'")
            sys.exit(1)
        smiles_used = smiles
        source = "user_smiles"
    elif sdf_path:
        template = Chem.MolFromMolFile(sdf_path, removeHs=False)
        if template is None:
            log(f"    [error] Could not parse --template-sdf '{sdf_path}'")
            sys.exit(1)
        smiles_used = Chem.MolToSmiles(Chem.RemoveHs(template))
        source = "user_sdf"
    elif not no_fetch and resname:
        log(f"    Fetching SMILES template for {resname} from RCSB CCD ...")
        smiles_used = fetch_smiles_from_rcsb(resname)
        if smiles_used:
            template = Chem.MolFromSmiles(smiles_used)
            source = "rcsb_ccd"

    if template is None:
        # Best effort: try to sanitize whatever RDKit could read.
        log("    [warn] No SMILES template available — running best-effort. "
            "Aromatic/Donor/Acceptor calls may be unreliable.")
        ligand = Chem.Mol(frag)
        problems = Chem.DetectChemistryProblems(ligand)
        for p in problems:
            log(f"        chemistry problem: {p.GetType()} {p.Message()}")
        try:
            Chem.SanitizeMol(ligand)
        except Exception as e:
            log(f"    [warn] SanitizeMol failed ({e}); proceeding with "
                f"unsanitized mol.")
        ligand = Chem.AddHs(ligand, addCoords=True)
        return ligand, None, "best_effort"

    # Bond-order fix from template, then add hydrogens with 3D coords.
    try:
        ligand = AllChem.AssignBondOrdersFromTemplate(template, frag)
    except Exception as e:
        log(f"    [warn] AssignBondOrdersFromTemplate failed ({e}); "
            "falling back to best-effort.")
        ligand = Chem.Mol(frag)
        try:
            Chem.SanitizeMol(ligand)
        except Exception:
            pass
        ligand = Chem.AddHs(ligand, addCoords=True)
        return ligand, smiles_used, "template_failed_best_effort"

    ligand = Chem.AddHs(ligand, addCoords=True)
    try:
        Chem.SanitizeMol(ligand)
    except Exception as e:
        log(f"    [warn] SanitizeMol after template fix failed: {e}")
    return ligand, smiles_used, source


def compute_ligand_features(ligand, factory):
    """Run the feature factory and return a list of dicts describing each
    feature (skipping the redundant fine-grained `Hydrophobe` family in favor
    of `LumpedHydrophobe`).
    """
    if ligand.GetNumConformers() < 1:
        log("    [error] Ligand has no 3D conformer; cannot compute "
            "feature positions.")
        sys.exit(1)
    feats = factory.GetFeaturesForMol(ligand)
    out = []
    for i, f in enumerate(feats):
        family = f.GetFamily()
        # Drop atom-granular Hydrophobe; LumpedHydrophobe gives one centroid
        # per region and is what pharmacophore tools normally show.
        if family == "Hydrophobe":
            continue
        pos = f.GetPos()
        out.append({
            "id": i,
            "family": family,
            "type": f.GetType(),
            "position": [round(pos.x, 3), round(pos.y, 3), round(pos.z, 3)],
            "tolerance": TOL_DEFAULTS.get(family, 1.0),
            "ligand_atom_ids": list(f.GetAtomIds()),
        })
    # Re-id sequentially after the Hydrophobe drop.
    for i, fdict in enumerate(out):
        fdict["id"] = i
    return out


def parse_pdb_atoms(pdb_path):
    """Parse ATOM/HETATM lines from a PDB textually. Returns a list of dicts.

    We re-parse from text rather than via RDKit because:
      • we only need names + coordinates, not chemistry, and
      • RDKit can fail/sanitize away parts of the protein on real-world PDBs.
    """
    atoms = []
    with open(pdb_path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip() or "A"
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            atoms.append({
                "record": line[:6].strip(),
                "atom_name": atom_name,
                "resname": resname,
                "chain": chain,
                "resnum": resnum,
                "pos": (x, y, z),
            })
    return atoms


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def compute_protein_features(atoms, ligand_positions, neighbor_radius):
    """Build a list of protein-side feature points within `neighbor_radius`
    of any ligand heavy atom. One feature dict per (atom, family).
    """
    # Pre-filter: only protein ATOM records (skip HETATMs handled separately).
    prot_atoms = [a for a in atoms if a["record"] == "ATOM"]

    # Find atoms within neighbor_radius of any ligand position.
    nearby = []
    for a in prot_atoms:
        for lp in ligand_positions:
            if dist(a["pos"], lp) <= neighbor_radius:
                nearby.append(a)
                break

    if not nearby:
        return []

    feats = []

    # Per-atom features (side-chain + backbone).
    for a in nearby:
        key = (a["resname"], a["atom_name"])
        families = list(PROT_ATOM_FEATURES.get(key, []))
        # Backbone N → Donor, O → Acceptor (skip terminal OXT — also Acceptor).
        if a["atom_name"] in BACKBONE_FEATURES and a["record"] == "ATOM":
            families.extend(BACKBONE_FEATURES[a["atom_name"]])
        if a["atom_name"] == "OXT":
            families.append("Acceptor")
        for fam in families:
            feats.append({
                "family": fam,
                "position": a["pos"],
                "residue": f"{a['resname']}{a['resnum']}",
                "chain": a["chain"],
                "atom_name": a["atom_name"],
            })

    # Aromatic ring centroids (compute per (chain, resnum) of Phe/Tyr/Trp/His).
    # Group nearby atoms by residue.
    by_res = {}
    for a in nearby:
        if a["resname"] in RING_ATOMS:
            by_res.setdefault(
                (a["chain"], a["resnum"], a["resname"]), []
            ).append(a)

    for (chain, resnum, resname), atom_list in by_res.items():
        ring_names = set(RING_ATOMS[resname])
        ring_atoms = [a for a in atom_list if a["atom_name"] in ring_names]
        if len(ring_atoms) < 4:
            continue  # incomplete ring
        cx = sum(a["pos"][0] for a in ring_atoms) / len(ring_atoms)
        cy = sum(a["pos"][1] for a in ring_atoms) / len(ring_atoms)
        cz = sum(a["pos"][2] for a in ring_atoms) / len(ring_atoms)
        feats.append({
            "family": "Aromatic",
            "position": (cx, cy, cz),
            "residue": f"{resname}{resnum}",
            "chain": chain,
            "atom_name": "ring",
        })
        # Aromatic side chains are also hydrophobic.
        feats.append({
            "family": "Hydrophobe",
            "position": (cx, cy, cz),
            "residue": f"{resname}{resnum}",
            "chain": chain,
            "atom_name": "ring",
        })

    # Zinc HETATMs as metal partners.
    for a in atoms:
        if a["record"] == "HETATM" and a["resname"] == "ZN":
            for lp in ligand_positions:
                if dist(a["pos"], lp) <= neighbor_radius:
                    feats.append({
                        "family": "Zn",
                        "position": a["pos"],
                        "residue": f"ZN{a['resnum']}",
                        "chain": a["chain"],
                        "atom_name": "ZN",
                    })
                    break

    return feats


def filter_by_interaction(lig_feats, prot_feats, cutoffs):
    """Keep only ligand features that have a complementary protein feature
    within the appropriate cutoff. Records the partner per kept feature.

    `cutoffs` overrides COMPLEMENTS_DEFAULT cutoffs by family-pair.
    """
    kept = []
    for lf in lig_feats:
        family = lf["family"]
        complements = COMPLEMENTS_DEFAULT.get(family, [])
        if not complements:
            continue
        best = None
        for prot_family, default_cutoff in complements:
            cutoff = cutoffs.get((family, prot_family), default_cutoff)
            for pf in prot_feats:
                if pf["family"] != prot_family:
                    continue
                d = dist(lf["position"], pf["position"])
                if d <= cutoff and (best is None or d < best["distance_A"]):
                    label = INTERACTION_LABEL.get((family, prot_family),
                                                  "contact")
                    best = {
                        "partner_residue": pf["residue"],
                        "partner_chain": pf["chain"],
                        "partner_atom": pf["atom_name"],
                        "distance_A": round(d, 3),
                        "type": label,
                    }
        if best is not None:
            kept.append({**lf, "interaction": best})
    return kept


# ── Output writers ────────────────────────────────────────────────────

def write_json(out_path, *, input_pdb, ligand_info, features,
               interaction_filtered, n_before, n_after,
               cutoffs_used):
    payload = {
        "input_pdb": input_pdb,
        "ligand": ligand_info,
        "tolerance_defaults_A": {
            k: v for k, v in TOL_DEFAULTS.items() if k != "Hydrophobe"
        },
        "interaction_cutoffs_A": cutoffs_used,
        "features": features,
        "metadata": {
            "rdkit_version": RDKIT_VERSION,
            "fdef": "BaseFeatures.fdef",
            "interaction_filtered": interaction_filtered,
            "n_features_before_filter": n_before,
            "n_features_after_filter": n_after,
        },
    }
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)


def write_pml(out_path, *, input_pdb, ligand_resname, features):
    pdb_basename = os.path.basename(input_pdb)
    lines = [
        f"# pharmacophore.pml — generated by pharmacophore-analyzer",
        f"# Source PDB: {input_pdb}",
        f"load {pdb_basename}, complex",
        f"hide everything",
        f"show cartoon, polymer",
        f"color grey80, polymer",
    ]
    if ligand_resname:
        lines.extend([
            f"show sticks, resn {ligand_resname}",
            f"color cyan, (resn {ligand_resname} and elem C)",
            f"util.cnc resn {ligand_resname}",
        ])
    lines.append("")
    lines.append("# Pharmacophore features")
    for f in features:
        family = f["family"]
        color = PML_COLORS.get(family, "white")
        label = PML_LABELS.get(family, family[:3].upper())
        x, y, z = f["position"]
        name = f"phar_{family}_{f['id']}"
        lines.append(
            f'pseudoatom {name}, pos=[{x:.3f}, {y:.3f}, {z:.3f}], '
            f'color={color}, label="{label}"'
        )
        lines.append(f"set sphere_scale, {f['tolerance']:.2f}, {name}")
    lines.extend([
        "",
        "show spheres, phar_*",
        "set sphere_transparency, 0.5, phar_*",
        "set label_size, 14",
        "set label_color, black",
        "bg_color white",
    ])
    if ligand_resname:
        lines.append(f"zoom (resn {ligand_resname}), 8")
    else:
        lines.append("zoom phar_*, 5")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ── Main ──────────────────────────────────────────────────────────────

def process_one(frag, *, input_pdb, factory, all_atoms, args, output_prefix):
    """Run the pipeline on a single ligand fragment. Writes JSON and PML."""
    resname, chain, resnum = get_residue_info(frag)
    log(f"\nLigand: resname={resname} chain={chain} resnum={resnum} "
        f"({frag.GetNumHeavyAtoms()} heavy atoms)")

    # 1. Bond orders + Hs
    ligand, smiles_used, smiles_source = prepare_ligand(
        frag,
        smiles=args.smiles,
        sdf_path=args.template_sdf,
        no_fetch=args.no_fetch,
    )
    log(f"    SMILES source: {smiles_source}"
        + (f"   ({smiles_used})" if smiles_used else ""))

    # 2. Ligand features
    lig_feats = compute_ligand_features(ligand, factory)
    log(f"    Ligand features (raw): {len(lig_feats)}")
    for fam in sorted({f['family'] for f in lig_feats}):
        n = sum(1 for f in lig_feats if f['family'] == fam)
        log(f"        {fam}: {n}")

    # 3. Protein-side features within neighbor_radius
    ligand_positions = []
    conf = ligand.GetConformer()
    for atom in ligand.GetAtoms():
        if atom.GetAtomicNum() == 1:  # skip Hs for proximity check
            continue
        p = conf.GetAtomPosition(atom.GetIdx())
        ligand_positions.append((p.x, p.y, p.z))

    prot_feats = compute_protein_features(
        all_atoms, ligand_positions, args.neighbor_radius,
    )
    log(f"    Protein-side features within {args.neighbor_radius} Å: "
        f"{len(prot_feats)}")

    # 4. Interaction filter (unless --include-all-features)
    cutoffs_used = {
        "Donor-Acceptor": args.cutoff_hbond,
        "Acceptor-Donor": args.cutoff_hbond,
        "PosIonizable-NegIonizable": args.cutoff_ionic,
        "NegIonizable-PosIonizable": args.cutoff_ionic,
        "Aromatic-Aromatic": args.cutoff_aromatic,
        "Aromatic-PosIonizable": args.cutoff_aromatic,
        "LumpedHydrophobe-Hydrophobe": args.cutoff_hydrophobic,
        "ZnBinder-Zn": 3.0,
    }
    n_before = len(lig_feats)
    if args.include_all_features:
        # No filter — but still annotate each feature with its nearest partner
        # if one exists, for downstream interpretation.
        cutoff_lookup = {
            ("Donor", "Acceptor"): args.cutoff_hbond,
            ("Acceptor", "Donor"): args.cutoff_hbond,
            ("PosIonizable", "NegIonizable"): args.cutoff_ionic,
            ("NegIonizable", "PosIonizable"): args.cutoff_ionic,
            ("Aromatic", "Aromatic"): args.cutoff_aromatic,
            ("Aromatic", "PosIonizable"): args.cutoff_aromatic,
            ("LumpedHydrophobe", "Hydrophobe"): args.cutoff_hydrophobic,
            ("Hydrophobe", "Hydrophobe"): args.cutoff_hydrophobic,
            ("ZnBinder", "Zn"): 3.0,
        }
        annotated = filter_by_interaction(lig_feats, prot_feats, cutoff_lookup)
        annotated_by_id = {f["id"]: f for f in annotated}
        out_features = []
        for lf in lig_feats:
            if lf["id"] in annotated_by_id:
                out_features.append(annotated_by_id[lf["id"]])
            else:
                out_features.append({**lf, "interaction": None})
        n_after = sum(1 for f in out_features if f.get("interaction"))
        interaction_filtered = False
        log(f"    --include-all-features set: emitting all {len(out_features)} "
            f"features ({n_after} have a protein partner)")
    else:
        cutoff_lookup = {
            ("Donor", "Acceptor"): args.cutoff_hbond,
            ("Acceptor", "Donor"): args.cutoff_hbond,
            ("PosIonizable", "NegIonizable"): args.cutoff_ionic,
            ("NegIonizable", "PosIonizable"): args.cutoff_ionic,
            ("Aromatic", "Aromatic"): args.cutoff_aromatic,
            ("Aromatic", "PosIonizable"): args.cutoff_aromatic,
            ("LumpedHydrophobe", "Hydrophobe"): args.cutoff_hydrophobic,
            ("Hydrophobe", "Hydrophobe"): args.cutoff_hydrophobic,
            ("ZnBinder", "Zn"): 3.0,
        }
        out_features = filter_by_interaction(lig_feats, prot_feats,
                                             cutoff_lookup)
        # Re-id sequentially.
        for i, f in enumerate(out_features):
            f["id"] = i
        n_after = len(out_features)
        interaction_filtered = True
        log(f"    Interaction-filtered features: {n_after} "
            f"(of {n_before})")

    # 5. Write outputs
    json_path = f"{output_prefix}.pharmacophore.json"
    pml_path = f"{output_prefix}.pharmacophore.pml"
    ligand_info = {
        "resname": resname,
        "chain": chain,
        "resnum": resnum,
        "smiles_used": smiles_used,
        "smiles_source": smiles_source,
    }
    write_json(json_path,
               input_pdb=input_pdb,
               ligand_info=ligand_info,
               features=out_features,
               interaction_filtered=interaction_filtered,
               n_before=n_before,
               n_after=n_after,
               cutoffs_used=cutoffs_used)
    write_pml(pml_path,
              input_pdb=input_pdb,
              ligand_resname=resname,
              features=out_features)
    log(f"    → {json_path}")
    log(f"    → {pml_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute a 3D pharmacophore from a protein-ligand PDB "
                    "using RDKit. Writes a JSON feature list and a PyMOL "
                    "PML scene.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  python pharmacophore_analysis.py 1IEP_clean.pdb -o 1IEP\n"
            "  python pharmacophore_analysis.py complex.pdb \\\n"
            "      --ligand-resname STI --smiles 'Cc1...' -o imatinib\n"
            "  python pharmacophore_analysis.py complex.pdb --no-fetch\n"
        ),
    )
    parser.add_argument("pdb_path",
                        help="Path to the input protein-ligand PDB file")
    parser.add_argument("--ligand-resname",
                        help="3-letter ligand residue code (e.g. STI). "
                        "Default: auto-detect.")
    parser.add_argument("--ligand-chain",
                        help="Restrict to a specific chain")
    parser.add_argument("--smiles",
                        help="Canonical SMILES template for bond-order fix")
    parser.add_argument("--template-sdf",
                        help="SDF file with the ligand template")
    parser.add_argument("--no-fetch", action="store_true",
                        help="Skip RCSB CCD SMILES fetch, do best-effort")
    parser.add_argument("--include-all-features", action="store_true",
                        help="Skip the interaction filter; emit every ligand "
                        "feature (still annotates partner if one exists)")
    parser.add_argument("-d", "--output-dir", default="pharmacophore",
                        help="Output directory (auto-created). The input PDB "
                        "is symlinked here so the PML's relative `load` works. "
                        "Default: 'pharmacophore'.")
    parser.add_argument("-o", "--output-prefix",
                        help="Output filename stem within --output-dir. "
                        "Default: PDB stem (e.g. '1IEP_clean'). Pass an "
                        "absolute path to bypass --output-dir entirely.")
    parser.add_argument("--cutoff-hbond", type=float, default=3.5,
                        help="Donor↔Acceptor cutoff in Å (default 3.5)")
    parser.add_argument("--cutoff-hydrophobic", type=float, default=4.5,
                        help="Hydrophobic↔Hydrophobic cutoff in Å (default 4.5)")
    parser.add_argument("--cutoff-aromatic", type=float, default=5.5,
                        help="Aromatic / π-cation cutoff in Å (default 5.5)")
    parser.add_argument("--cutoff-ionic", type=float, default=5.0,
                        help="PosIonizable↔NegIonizable cutoff in Å (default 5.0)")
    parser.add_argument("--neighbor-radius", type=float, default=5.0,
                        help="Pocket-residue inclusion radius in Å (default 5.0)")
    args = parser.parse_args()

    if not os.path.isfile(args.pdb_path):
        log(f"[error] PDB not found: {args.pdb_path}")
        sys.exit(1)

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Symlink the input PDB into output_dir so the PML's relative `load
    # <basename>` resolves regardless of where PyMOL is launched from.
    input_real = os.path.realpath(args.pdb_path)
    input_basename = os.path.basename(args.pdb_path)
    symlink_path = os.path.join(output_dir, input_basename)
    if os.path.lexists(symlink_path):
        try:
            existing_real = os.path.realpath(symlink_path)
        except OSError:
            existing_real = None
        if existing_real == input_real:
            pass  # already correct
        elif os.path.islink(symlink_path):
            os.remove(symlink_path)
            os.symlink(input_real, symlink_path)
            log(f"[info] Updated symlink {symlink_path} -> {input_real}")
        else:
            log(f"[warn] {symlink_path} exists and is not a symlink to the "
                "input PDB; leaving it. The PML will load that file instead.")
    else:
        os.symlink(input_real, symlink_path)
        log(f"[info] Symlinked {symlink_path} -> {input_real}")

    output_prefix_stem = args.output_prefix or os.path.splitext(
        input_basename
    )[0]

    # 1. Read full structure
    log(f"Reading {args.pdb_path} ...")
    mol = Chem.MolFromPDBFile(args.pdb_path, removeHs=False, sanitize=False)
    if mol is None:
        log(f"[error] RDKit failed to parse {args.pdb_path}")
        sys.exit(1)

    # 2. Identify ligand fragment(s)
    candidates = identify_ligands(
        mol, want_resname=args.ligand_resname, want_chain=args.ligand_chain,
    )
    if not candidates:
        log("[error] No ligand fragment matched the criteria. "
            "Try --ligand-resname or check the input PDB.")
        sys.exit(1)
    log(f"Found {len(candidates)} ligand fragment(s)")

    # 3. Build feature factory
    fdef_path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    if not os.path.isfile(fdef_path):
        log(f"[error] BaseFeatures.fdef not found at {fdef_path}")
        sys.exit(1)
    factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)

    # 4. Pre-parse all PDB atoms for protein-side feature inference
    all_atoms = parse_pdb_atoms(args.pdb_path)
    log(f"Parsed {len(all_atoms)} atom records from PDB")

    # 5. Process each ligand copy
    multi = len(candidates) > 1
    for frag in candidates:
        resname, chain, resnum = get_residue_info(frag)
        suffix = f"_chain{chain}_res{resnum}" if multi else ""
        # os.path.join honours absolute prefixes (so users can bypass --output-dir
        # by passing an absolute --output-prefix).
        prefix = os.path.join(output_dir, output_prefix_stem + suffix)
        try:
            process_one(
                frag,
                input_pdb=args.pdb_path,
                factory=factory,
                all_atoms=all_atoms,
                args=args,
                output_prefix=prefix,
            )
        except SystemExit:
            raise
        except Exception as e:
            log(f"    [error] processing {resname} chain {chain}: {e}")
            continue

    log("\nDone.")


if __name__ == "__main__":
    main()
