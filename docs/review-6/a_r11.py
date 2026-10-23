"""Review 6, R11 attacks on master's coordinate_kind / read_structure.
usage: a_r11.py MASTER_DIR FX_DIR. Writes fixtures to FX_DIR, prints the
reader's verdict per case and its cost under tests/data/r11's severity scale."""
import math
import sys
from pathlib import Path

import numpy as np

code, fx = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(code))
import esp_map as em  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

fx.mkdir(parents=True, exist_ok=True)
COST = {"silent wrong geometry": 3, "recorded discard": 2,
        "overrules the 2D header": 2, "BLOCK": 1, "correct": 0}


def verdict(p, record=None):
    mol, probs, n = em.read_structure(str(p), record=record)
    if mol is not None:
        return mol.GetProp("_coords"), mol
    return ("ambiguous" if probs and probs[0][1] == "DIMENSION_AMBIGUOUS"
            else f"BLOCK {probs[0][1]}"), None


def write_mol(path, mol, header="3D", crlf=False):
    blk = Chem.MolToMolBlock(mol)
    lines = blk.split("\n")
    lines[1] = lines[1][:20].ljust(20) + header + lines[1][22:]
    txt = "\n".join(lines)
    if crlf:
        txt = txt.replace("\n", "\r\n")
    Path(path).write_bytes(txt.encode())
    return path


def from_xyz(smiles, xyz):
    """mol with the given heavy-atom (or all-atom) coordinates, in SMILES order."""
    m = Chem.MolFromSmiles(smiles)
    if len(xyz) != m.GetNumAtoms():
        m = Chem.AddHs(m)
    assert len(xyz) == m.GetNumAtoms(), (smiles, len(xyz), m.GetNumAtoms())
    c = Chem.Conformer(m.GetNumAtoms())
    for i, p in enumerate(xyz):
        c.SetAtomPosition(i, [float(v) for v in p])
    c.Set3D(True)
    m.AddConformer(c, assignId=True)
    return m


def mmff(smiles, seed=7):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    props = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
    if props is None:
        AllChem.UFFOptimizeMolecule(m, maxIters=5000)
        ff = "UFF"
    else:
        AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=5000)
        ff = "MMFF94s"
    X = m.GetConformer().GetPositions()
    X = X - X.mean(0)
    dev = np.abs(X @ np.linalg.svd(X)[2][-1]).max()
    return m, ff, dev


def report(name, path, kind, prefer, outcomes, note=""):
    got, mol = verdict(path)
    key = got if got in outcomes else None
    cost = COST[outcomes[key]] if key else 1
    extra = ""
    if mol is not None and mol.GetNumConformers():
        X = mol.GetConformer().GetPositions(); X = X - X.mean(0)
        dev = np.abs(X @ np.linalg.svd(X)[2][-1]).max() if len(X) > 3 else 0.0
        sp3, real, ring4 = em.flat_sp3_problems(mol)
        extra = f" dev={dev:.3f} sp3={sp3} real_angles={real}"
    print(f"{name:44} kind={kind:8} got={got:22} prefer={prefer:9} "
          f"cost={cost} ({outcomes.get(key, got) if key else got}){extra} {note}")


GEO = {"geometry": "correct", "ambiguous": "BLOCK", "drawing": "recorded discard"}
DRW = {"drawing": "correct", "ambiguous": "BLOCK", "geometry": "silent wrong geometry"}
HDR2 = {"drawing": "correct", "ambiguous": "BLOCK", "geometry": "overrules the 2D header"}

print("=== 1. real planar molecules, reference geometries (A) ===")
# tricyanomethanide, D3h: C-C 1.40, C#N 1.16 (crystal structures, [I])
X = [[0, 0, 0]]
for k in range(3):
    t = 2 * math.pi * k / 3
    X += [[1.40 * math.cos(t), 1.40 * math.sin(t), 0], [2.56 * math.cos(t), 2.56 * math.sin(t), 0]]
m = from_xyz("[C-](C#N)(C#N)C#N", X)
print("  RDKit hybridisation of the central C:", m.GetAtomWithIdx(0).GetHybridization())
report("real_tricyanomethanide", write_mol(fx / "real_tricyanomethanide.mol", m), "geometry", "geometry", GEO)
# dicyanamide, C2v: N-C 1.31, C#N 1.16, C-N-C 120 deg, N-C#N linear ([I])
a = math.radians(60)
c1 = [-1.31 * math.sin(a), -1.31 * math.cos(a), 0]
c2 = [1.31 * math.sin(a), -1.31 * math.cos(a), 0]
n1 = [-2.47 * math.sin(a), -2.47 * math.cos(a), 0]
n2 = [2.47 * math.sin(a), -2.47 * math.cos(a), 0]
m = from_xyz("N#C[N-]C#N", [n1, c1, [0, 0, 0], c2, n2])
print("  RDKit hybridisation of the central N:", m.GetAtomWithIdx(2).GetHybridization())
report("real_dicyanamide", write_mol(fx / "real_dicyanamide.mol", m), "geometry", "geometry", GEO)
# disiloxane heavy-atom skeleton: Si-O 1.634, Si-O-Si 144.1 deg (gas ED, [I])
h = math.radians(144.1 / 2)
m = from_xyz("[SiH3]O[SiH3]", [[-1.634 * math.sin(h), -1.634 * math.cos(h), 0], [0, 0, 0],
                               [1.634 * math.sin(h), -1.634 * math.cos(h), 0]])
report("real_disiloxane_heavy_only", write_mol(fx / "real_disiloxane_noH.mol", m), "geometry", "geometry", GEO,
       "(water_heavy_only is in the set at 104.5 deg)")
h = math.radians(110.9 / 2)
m = from_xyz("ClOCl", [[-1.70 * math.sin(h), -1.70 * math.cos(h), 0], [0, 0, 0],
                       [1.70 * math.sin(h), -1.70 * math.cos(h), 0]])
report("real_Cl2O (control, 110.9 deg)", write_mol(fx / "real_Cl2O.mol", m), "geometry", "geometry", GEO)

print("=== 2. real planar molecules from MMFF94s (a geometry any user could supply) ===")
for name, smi in [("squaric_acid", "OC1=C(O)C(=O)C1=O"), ("tropolone", "O=C1C=CC=CC=C1O"),
                  ("vinamidinium", "NC=CC=[NH2+]"), ("allyl_anion", "C=C[CH2-]"),
                  ("acac_anion", "CC(=O)C=C([O-])C"), ("hexapentaene", "C=C=C=C=C=C"),
                  ("malonaldehyde_enol", "O=CC=CO"), ("glyoxal", "O=CC=O"),
                  ("dicyanamide_mmff", "N#C[N-]C#N"), ("tricyanomethanide_mmff", "[C-](C#N)(C#N)C#N")]:
    try:
        m, ff, dev = mmff(smi)
    except Exception as e:
        print(f"{name:44} could not build: {e}")
        continue
    if dev > em.FLAT_TOL:
        print(f"{name:44} {ff} geometry not flat (dev {dev:.3f}); not a flat-file case")
        continue
    p = write_mol(fx / f"real_{name}.mol", m)
    r = Chem.MolFromMolFile(str(p), removeHs=False)
    cr = r.GetConformer()
    lp = em.flat_length_problems(r, [(cr.GetAtomPosition(b.GetBeginAtomIdx())
                                      - cr.GetAtomPosition(b.GetEndAtomIdx())).Length()
                                     for b in r.GetBonds()])
    report(f"real_{name} [{ff}, dev {dev:.3f}]", p, "geometry", "geometry", GEO, f"length_problems={lp}")

print("=== 3. header handling: real C2F4 (z = 0, 'either' C=C) under file variants ===")
src = code / "tests/data/r11/real_C2F4.mol"
base = src.read_text()
report("C2F4 as committed (control)", src, "geometry", "geometry", GEO)
p = fx / "C2F4_crlf.mol"; p.write_bytes(base.replace("\n", "\r\n").encode())
report("C2F4 CRLF line endings", p, "geometry", "geometry", GEO)
sdf = base.rstrip("\n") + "\n$$$$\n"
p = fx / "C2F4_sdf.sdf"; p.write_text(sdf)
report("C2F4 as one-record SDF", p, "geometry", "geometry", GEO)
p = fx / "C2F4_sdf_crlf.sdf"; p.write_bytes(sdf.replace("\n", "\r\n").encode())
report("C2F4 one-record SDF, CRLF", p, "geometry", "geometry", GEO)
p = fx / "C2F4_sdf_props.sdf"
p.write_text(base.rstrip("\n") + "\n> <note>\nline one\n\n> <geometry>\nB3LYP/def2-TZVP\n\n$$$$\n")
report("C2F4 SDF with data fields", p, "geometry", "geometry", GEO)
p = fx / "C2F4_sdf_trailing.sdf"; p.write_text(sdf + "\n# exported by tool X\n")
print("  dims:", em._molfile_header_dims(p), "RDKit records:", len(list(Chem.SDMolSupplier(str(p), removeHs=False))))
report("C2F4 SDF + trailing comment line", p, "geometry", "geometry", GEO)
p = fx / "C2F4_sdf_blank_record.sdf"; p.write_text("$$$$\n" + sdf)
print("  blank-first dims:", em._molfile_header_dims(p), "RDKit records:",
      [m is not None for m in Chem.SDMolSupplier(str(p), removeHs=False)])
g = em.read_structure(str(p))
print("  read_structure (no --record):", g[1][:1] if g[0] is None else g[0].GetProp("_coords"), "n_found", g[2])
p2 = fx / "C2F4_two_records_one_bad.sdf"
bad = base.replace("  2  4  2  3", "  2  4  3  0")    # C#C with two F each: pentavalent C, unparseable
p2.write_text(bad.rstrip("\n") + "\n$$$$\n" + sdf)
print("  bad+good dims:", em._molfile_header_dims(p2), "RDKit:",
      [m is not None for m in Chem.SDMolSupplier(str(p2), removeHs=False)])
for r in (None, 1, 2):
    g = em.read_structure(str(p2), record=r)
    print(f"  bad+good --record {r}:", (g[0].GetProp("_coords") + " title=" + g[0].GetProp("_Name")
                                         if g[0] is not None else g[1][0][1:3]), "n_found", g[2])
v3 = Chem.MolToV3KMolBlock(Chem.MolFromMolFile(str(src), removeHs=False))
lines = v3.split("\n"); lines[1] = lines[1][:20].ljust(20) + "3D" + lines[1][22:]
p = fx / "C2F4_v3000.mol"; p.write_text("\n".join(lines))
print("  v3000 Is3D after RDKit:", Chem.MolFromMolFile(str(p), removeHs=False).GetConformer().Is3D(),
      "dims:", em._molfile_header_dims(p))
report("C2F4 V3000", p, "geometry", "geometry", GEO)
# header line shorter than 22 columns (some writers drop the program/date fields)
lines = base.split("\n"); lines[1] = "  3D"
p = fx / "C2F4_short_header.mol"; p.write_text("\n".join(lines))
print("  short-header dims:", em._molfile_header_dims(p))
report("C2F4 header line '  3D' (short)", p, "geometry", "geometry", GEO)

print("=== 4. 2D header with non-zero z ===")
d = Chem.MolFromSmiles("C[C@H](O)[C@@H](N)C(=O)O")
AllChem.Compute2DCoords(d)
conf = d.GetConformer()
for i, dz in ((0, 0.8), (2, -0.8), (4, 0.8)):
    q = conf.GetAtomPosition(i); conf.SetAtomPosition(i, (q.x, q.y, dz))
p = write_mol(fx / "draw_threonine_2Dheader_zwedge.mol", d, header="2D")
print("  RDKit Is3D:", Chem.MolFromMolFile(str(p)).GetConformer().Is3D(), "dims", em._molfile_header_dims(p))
report("draw threonine, 2D header, wedge z +-0.8", p, "drawing", "drawing", DRW)
for dz in (0.11, 0.3):
    d2 = Chem.MolFromSmiles("Clc1ccccc1"); AllChem.Compute2DCoords(d2)
    cf = d2.GetConformer(); q = cf.GetAtomPosition(0); cf.SetAtomPosition(0, (q.x, q.y, dz))
    p = write_mol(fx / f"draw_chlorobenzene_2Dheader_z{dz}.mol", d2, header="2D")
    report(f"draw chlorobenzene, 2D header, one z={dz}", p, "drawing", "drawing", DRW)
m, ff, dev = mmff("CC(O)CC")
p = write_mol(fx / "real_2butanol_2Dheader.mol", m, header="2D")
report("real 2-butanol, 2D header, non-flat", p, "geometry", "drawing", HDR2)
