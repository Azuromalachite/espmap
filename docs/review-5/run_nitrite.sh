source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review5-work/env.sh
cd $HOME
$PY - <<'PY' > $W/e_nitrite.out 2>&1
import json, os, subprocess, sys
from rdkit import Chem
W, P3, REPO = os.environ["W"], os.environ["P3"], os.environ["REPO"]
def atoms(path):
    m = Chem.MolFromMolFile(path, removeHs=False)
    return "; ".join(f"{a.GetSymbol()} {p.x:.4f} {p.y:.4f} {p.z:.4f}"
                     for a, p in zip(m.GetAtoms(), m.GetConformer().GetPositions().tolist() and
                                     [m.GetConformer().GetAtomPosition(i) for i in range(m.GetNumAtoms())]))
flat = atoms(P3 + "/bench/fixtures/nitrite.sdf")
tilt = atoms(P3 + "/fixtures/geometries/nitrite.sdf")
for basis in ("sto-3g", "def2-svp"):
    for lab, a in (("planar (bench frame)", flat), ("tilted (P3 input)", tilt)):
        for sp in ("0.26", "0.16"):
            js = f"{W}/out/nit_{basis}_{lab.split()[0]}_{sp}.json"
            subprocess.run([sys.executable, REPO + "/esp_map.py", "--atoms", a, "--charge", "-1",
                            "--solvent", "water", "--basis", basis, "--spacing", sp, "--size", "80",
                            "--out", W + "/out/nit.png", "--json", js], capture_output=True)
            r = json.load(open(js))
            m = {s["label"]: s["nearest_nucleus"]["V_area_mean_kcal_mol"] for s in r["sites"]["reaction_sites"]}
            ar = {s["label"]: s["nearest_nucleus"]["area_A2"] for s in r["sites"]["reaction_sites"]}
            print(f"{basis:8} {lab:22} spacing {sp}: |O1-O2| mean {abs(m['O1'] - m['O2']):.3f} kcal/mol, "
                  f"areas {ar['O1']:.3f}/{ar['O2']:.3f} A2, V_min {r['results']['V_min_kcal_mol']:.3f}", flush=True)
PY
cat $W/e_nitrite.out
