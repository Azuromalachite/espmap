import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tree"))
import esp_map as em
from rdkit import Chem
from rdkit.Chem import AllChem
em._quiet_rdkit()
FX = os.path.join(HERE, "fx")
def is_depiction(path, mol):
    ext = path.lower().rsplit(".", 1)[-1]
    if ext in ("xyz", "pdb"):
        return False                                   # these formats are never drawings
    if mol.GetNumConformers() == 0:
        return True
    c = mol.GetConformer(); P = c.GetPositions()
    flat = np.abs(P[:, 2]).max() < 1e-4
    if not c.Is3D():
        return True                                    # file says 2D: believe it
    if not flat:
        return False
    # flagged 3D but z == 0: a depiction draws every bond at one length;
    # only decide that way if the molecule HAS more than one bond type
    kinds = {tuple(sorted((b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()))) + (str(b.GetBondType()),) for b in mol.GetBonds()}
    L = np.array([np.linalg.norm(P[b.GetBeginAtomIdx()] - P[b.GetEndAtomIdx()]) for b in mol.GetBonds()])
    return len(kinds) > 1 and L.std() / L.mean() < 0.03
m2 = Chem.MolFromSmiles("Clc1ccccc1"); AllChem.Compute2DCoords(m2)
b = Chem.MolToMolBlock(m2).splitlines(); b[1] = b[1].replace("2D", "3D")
open(os.path.join(FX, "rd2d_flagged3d.mol"), "w").write("\n".join(b) + "\n")
cases = {"planar 3D chlorobenzene (z=0, flag 3D)": os.path.join(FX, "chlorobenzene_3d_planar.mol"),
         "examples/drawn_2d.mol": os.path.join(HERE, "tree/examples/drawn_2d.mol"),
         "examples/real_3d.mol": os.path.join(HERE, "tree/examples/real_3d.mol"),
         "examples/real_3d_noH.mol": os.path.join(HERE, "tree/examples/real_3d_noH.mol"),
         "examples/drawn.cdxml": os.path.join(HERE, "tree/examples/drawn.cdxml"),
         "examples/water.xyz": os.path.join(HERE, "tree/examples/water.xyz"),
         "water in xy plane (xyz)": os.path.join(FX, "water_xy.xyz"),
         "RDKit 2D depiction": os.path.join(FX, "rd2d.mol"),
         "2D depiction mislabelled 3D": os.path.join(FX, "rd2d_flagged3d.mol")}
for k, p in cases.items():
    mm, _, _ = em.read_structure(p)
    shipped = "DISCARDED" in em.from_file(p)[1]["source"]
    print(f"  {k:40s} shipped discards={shipped!s:5s} proposed discards={is_depiction(p, mm)}")
