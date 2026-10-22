import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tree"))
import esp_map as em
from rdkit.Chem import AllChem, Descriptors3D
smi = "CC(=O)OCC/C=C\CCCCCCCC/C=C\CCCC"
_, info = em.from_smiles(smi, n_confs=100)
mol = info["mol"]
props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
cids = [c.GetId() for c in mol.GetConformers()]
E = np.array([AllChem.MMFFGetMoleculeForceField(mol, props, confId=c).CalcEnergy() for c in cids])
rg = {c: Descriptors3D.RadiusOfGyration(mol, confId=c) for c in cids}
rows = [json.loads(l) for l in open(os.path.join(HERE, "tree/bench/conformers.jsonl")) if "conf_id" in l]
print("lowest-E conf id now:", cids[int(np.argmin(E))], " recorded: 4")
for r in rows:
    print(f"  conf {r['conf_id']:3d}: recorded rg {r['rg']:.6f}  recomputed {rg[r['conf_id']]:.6f}  dE rec {r['mmff_dE']:.4f} now {E[cids.index(r['conf_id'])]-E.min():.4f}")
