"""Fault-inject the corrections proposed in review 3 (fixtures, not patches to the tree)."""
import os, sys, io, contextlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tree"))
import esp_map as em
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
em._quiet_rdkit()
FX = os.path.join(HERE, "fx")

# ---- C1: approved-vs-computed from the COORDINATES that reach PySCF
def smiles_of_spec(spec, charge=0):
    n = len(spec.split(";"))
    g = Chem.MolFromXYZBlock(f"{n}\n\n" + "\n".join(t.strip() for t in spec.split(";")))
    rdDetermineBonds.DetermineBonds(g, charge=charge)
    Chem.AssignStereochemistryFrom3D(g)
    return Chem.MolToSmiles(Chem.RemoveHs(g))

def flat(smi):  # stereo-free canonical form
    m = Chem.MolFromSmiles(smi); Chem.RemoveStereochemistry(m); return Chem.MolToSmiles(m)

print("== C1 coordinate-derived check")
for label, appr_smi, spec_smi in [("healthy CCO", "CCO", "CCO"),
                                  ("healthy acetate", "CC(=O)[O-]", "CC(=O)[O-]"),
                                  ("injected: geometry of CC=O under approval of CCO", "CCO", "CC=O"),
                                  ("injected: E isomer computed for Z approval", "C/C=C\C", "C/C=C/C")]:
    spec, gi = em.from_smiles(spec_smi)
    approved = Chem.MolToSmiles(Chem.RemoveHs(Chem.MolFromSmiles(appr_smi)))
    old = Chem.CanonSmiles(appr_smi) == approved   # what main() compares on the SMILES path (gi['smiles'] is the input string)
    got = smiles_of_spec(spec, gi["charge"])
    new_ok = (got == approved) or (flat(got) == approved and flat(approved) == approved)
    print(f"  {label:48s} shipped check passes={old}  coordinates say {got:22s} proposed check passes={new_ok}")

# ---- C2: PDB read through coordinates + an H-count guard
print("== C2 PDB without CONECT")
pdb = os.path.join(FX, "chlorobenzene.pdb")
m = Chem.MolFromPDBFile(pdb, removeHs=False)
n_explicit_h = sum(a.GetAtomicNum() == 1 for a in m.GetAtoms())
n_after = sum(a.GetAtomicNum() == 1 for a in Chem.AddHs(m).GetAtoms())
print(f"  shipped reader: {Chem.MolToSmiles(Chem.RemoveHs(m))}; file has {n_explicit_h} H, AddHs makes {n_after}"
      f" -> guard 'file has H but AddHs adds more' trips: {n_explicit_h > 0 and n_after > n_explicit_h}")
for ok_file in ("real_3d.mol", "real_3d_noH.mol"):
    mm, _, _ = em.read_structure(os.path.join(HERE, "tree/examples", ok_file))
    h0 = sum(a.GetAtomicNum() == 1 for a in mm.GetAtoms()); h1 = sum(a.GetAtomicNum() == 1 for a in Chem.AddHs(mm).GetAtoms())
    print(f"  control {ok_file:16s}: file H {h0}, after AddHs {h1} -> guard trips: {h0 > 0 and h1 > h0}")
c = m.GetConformer()
xyzb = f"{m.GetNumAtoms()}\n\n" + "\n".join(f"{a.GetSymbol()} {p.x} {p.y} {p.z}" for a, p in zip(m.GetAtoms(), c.GetPositions().tolist() and [c.GetAtomPosition(i) for i in range(m.GetNumAtoms())]))
g = Chem.MolFromXYZBlock(xyzb); rdDetermineBonds.DetermineBonds(g, charge=0)
print(f"  PDB coordinates through the XYZ route: {Chem.MolToSmiles(Chem.RemoveHs(g))}")

# ---- C3: depiction detector that does not confuse 'planar' with 'drawn'
print("== C3 depiction vs planar 3D")
def is_depiction(mol):
    if mol.GetNumConformers() == 0: return True
    c = mol.GetConformer(); P = c.GetPositions()
    if np.abs(P[:, 2]).max() > 1e-4: return False
    L = np.array([np.linalg.norm(P[b.GetBeginAtomIdx()] - P[b.GetEndAtomIdx()]) for b in mol.GetBonds()])
    return L.std() / L.mean() < 0.03      # RDKit/ChemDraw depictions draw every bond at one length
cases = {"planar 3D chlorobenzene (z=0)": os.path.join(FX, "chlorobenzene_3d_planar.mol"),
         "examples/drawn_2d.mol": os.path.join(HERE, "tree/examples/drawn_2d.mol"),
         "examples/real_3d.mol": os.path.join(HERE, "tree/examples/real_3d.mol"),
         "examples/drawn.cdxml": os.path.join(HERE, "tree/examples/drawn.cdxml"),
         "water in xy plane (xyz)": os.path.join(FX, "water_xy.xyz")}
m2 = Chem.MolFromSmiles("Clc1ccccc1"); AllChem.Compute2DCoords(m2); p2 = os.path.join(FX, "rd2d.mol"); Chem.MolToMolFile(m2, p2)
cases["RDKit 2D depiction"] = p2
for k, p in cases.items():
    mm, _, _ = em.read_structure(p)
    shipped = "DISCARDED" in em.from_file(p)[1]["source"]
    print(f"  {k:32s} shipped discards={shipped!s:5s} proposed discards={is_depiction(mm)}")

# ---- C4: --selftest against an injected kernel-sign flip in esp_at_points itself
print("== C4 selftest with the sign flip injected into esp_at_points")
orig = em.esp_at_points
def flipped(mol, dm, pts, mem_budget_gb=0.5):
    nuc = orig(mol, np.zeros_like(dm), pts)       # nuclear term only
    return nuc + (nuc - orig(mol, dm, pts))        # nuclear + electronic (wrong sign)
em.esp_at_points = flipped
buf = io.StringIO()
with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
    rc = em.selftest()
em.esp_at_points = orig
print("  selftest return code with flipped kernel:", rc)
print("  " + "\n  ".join(l for l in buf.getvalue().splitlines() if "kcal/mol" in l or "lone pair negative" in l))
