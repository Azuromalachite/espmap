import os, sys, json
REPO = os.environ["REPO"]; sys.path.insert(0, REPO)
import esp_map as em
from rdkit import Chem
em._quiet_rdkit()
r = json.load(open(os.environ["W"] + "/out/r8.json"))
spec = em.atoms_from_geometry(r["geometry"]).replace("; ", ";")
appr = Chem.MolFromSmiles("[NH3+]CC(=O)[O-]")
print("record: validity", r["validity"], "structure_verified", r["input"]["structure_verified"],
      "computed_smiles", r["input"]["computed_smiles"], "chemical_match", r["sites"]["identity"]["chemical_match"])
print("re-verifying the RECORDED (computed) geometry against the approval:", em.verify_computed(spec, appr))
from rdkit.Chem import rdDetermineBonds
rows = [t.split() for t in spec.split(";")]
g = Chem.MolFromXYZBlock(f"{len(rows)}\n\n" + "\n".join(" ".join(x) for x in rows))
rdDetermineBonds.DetermineBonds(g, charge=0)
print("perceived from the recorded geometry:", Chem.MolToSmiles(Chem.RemoveHs(g)))
