import os, sys
sys.path.insert(0, os.environ["REPO"]); import esp_map as em
from rdkit import Chem
em._quiet_rdkit()
for s in ["CC(C)O", "[CH3:1][CH:2]([CH3:3])[OH:4]", "CC(C)C(C)C", "[CH3:1][CH:2]([CH3:3])[CH:4]([CH3:5])[CH3:6]",
          "OC(=O)CC(O)(CC(=O)O)C(=O)O", "[OH:1][C:2](=[O:3])[CH2:4][C:5]([OH:6])([CH2:7][C:8]([OH:9])=[O:10])[C:11]([OH:12])=[O:13]"]:
    p = Chem.MolFromSmiles(s)
    v, f = em.check_structure(Chem.RemoveHs(p), maps_from=p, smiles=s)
    print(f"{s:60s} {v:5s}", [(x[0], x[1]) for x in f if "UNSPEC" in x[1]])
