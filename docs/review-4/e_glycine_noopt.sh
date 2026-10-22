source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review4-work/env.sh
cd $W/out
$PY $REPO/esp_map.py --smiles "[NH3+]CC(=O)[O-]" --basis sto-3g --size 80 --out r9.png --json r9.json 2>&1 | grep "V = \|SCF -"
grep -c . r8.log; grep -n "optimis\|Geometry optimization\|converged" r8.log | head -3
$PY - <<'P'
import json
for f in ("r8.json","r9.json"):
    r=json.load(open(f)); print(f, r["settings"]["optimize"], r["analysis_key"][:10], r["results"]["V_min_kcal_mol"], r["results"]["V_max_kcal_mol"], r["geometry"]["sha1"][:10])
P
