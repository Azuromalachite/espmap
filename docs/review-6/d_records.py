"""Review 6, R8/R9 and Arrowfield section B: the forced and invalid paths of
the record and the NPZ, at STO-3G, through the CLI (main()).
usage: d_records.py CODE_DIR OUT_DIR"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

CODE, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(CODE / "tests"))
from record_tree import load_schema, tree_problems  # noqa: E402
SCHEMA = load_schema()

# an ethanol bond table whose coordinates have C1 and O2 swapped: the
# coordinates are a different skeleton (O between the carbons, four-coordinate)
m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
AllChem.EmbedMolecule(m, randomSeed=7)
AllChem.MMFFOptimizeMolecule(m)
c = m.GetConformer()
p1, p2 = c.GetAtomPosition(1), c.GetAtomPosition(2)
c.SetAtomPosition(1, p2)
c.SetAtomPosition(2, p1)
Chem.MolToMolFile(m, str(OUT / "ethanol_swapped.mol"))
# an ethanol whose 3D coordinates are fine but one H sits on O's other side:
# the same skeleton -> verified; used as a control
m2 = Chem.AddHs(Chem.MolFromSmiles("CCO"))
AllChem.EmbedMolecule(m2, randomSeed=7)
Chem.MolToMolFile(m2, str(OUT / "ethanol_ok.mol"))

CASES = {
    "valid_methanol": ["--smiles", "CO"],
    "h2": ["--atoms", "H 0 0 0; H 0 0 0.74"],
    "alh4_water": ["--smiles", "[AlH4-]", "--solvent", "water"],
    "forced_two_waters": ["--atoms", "O 0 0 0; H 0.96 0 0; H -0.24 0.93 0; "
                          "O 3 0 0; H 3.96 0 0; H 2.76 0.93 0", "--force"],
    "forced_charge_override": ["--smiles", "O", "--charge", "2", "--force"],
    "fluoride_gas": ["--smiles", "[F-]"],
    "verify_fail_refused": ["--file", str(OUT / "ethanol_swapped.mol")],
    "verify_fail_forced": ["--file", str(OUT / "ethanol_swapped.mol"), "--force"],
    "ethanol_file_ok": ["--file", str(OUT / "ethanol_ok.mol")],
}
PY = sys.executable
for name, argv in CASES.items():
    js, npz = OUT / f"{name}.json", OUT / f"{name}.npz"
    for f in (js, npz):
        f.unlink(missing_ok=True)
    p = subprocess.run([PY, str(CODE / "esp_map.py"), *argv, "--basis", "sto-3g",
                        "--size", "60", "--out", str(OUT / f"{name}.png"),
                        "--json", str(js), "--npz", str(npz)],
                       capture_output=True, text=True, cwd=OUT)
    print(f"== {name}: exit {p.returncode}")
    if not js.exists():
        print("   no record;", (p.stderr.strip().splitlines() or [""])[-1][:300])
        continue
    r = json.loads(js.read_text())
    inp = r["input"]
    print("   validity", r["validity"])
    print("   gate.forced", (inp["gate"] or {}).get("forced"), "| structure_verified",
          inp["structure_verified"], "| computed_smiles", inp["computed_smiles"],
          "| approved_smiles", inp["approved_smiles"])
    print("   sites_reason:", r.get("sites_reason"))
    print("   electronic problems:", [x[:40] for x in r["results"]["electronic"]["problems"]])
    if r["sites"]:
        print("   partition", r["sites"]["partition"])
    print("   schema problems:", tree_problems(SCHEMA["record"], r))
    if npz.exists():
        z = np.load(npz)
        tri = sorted(k for k in z.files if k.startswith("tri_atom"))
        neg = {k: int((z[k] < 0).sum()) for k in tri}
        print("   npz tri_atom arrays:", tri, "negative entries:", neg,
              "| npz validity:", json.loads(str(z["validity"])))
        print("   npz files == schema npz?", sorted(z.files) == sorted(SCHEMA["npz"]),
              "extra:", sorted(set(z.files) - set(SCHEMA["npz"])),
              "missing:", sorted(set(SCHEMA["npz"]) - set(z.files)))
        print("   presence rule (tri_atom present iff validity.sites):",
              bool(tri) == bool(r["validity"]["sites"]))
