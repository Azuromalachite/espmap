source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review4-work/env.sh
mkdir -p $W/out; cd $W/out
E="$PY $REPO/esp_map.py --basis sto-3g --size 80"
WATER="O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692"
echo "######## R1 mapped ethanol"
$E --smiles "[CH3:1][CH2:2][OH:3]" --out r1.png --json r1.json --npz r1.npz 2>&1 | tail -4
echo "######## R2a H2 via --atoms"
$E --atoms "H 0 0 0; H 0 0 0.74" --out r2.png --json r2.json 2>&1 | tail -4; ls r2.* 2>&1
echo "######## R2b H2 via --smiles"
$E --smiles "[H][H]" --out r2b.png --json r2b.json 2>&1 | tail -4; ls r2b.* 2>&1
echo "######## R3 --optimize mapped methanol"
$E --smiles "[CH3:1][OH:2]" --optimize --out r3.png --json r3.json 2>&1 | tail -4
echo "######## R4 clipped water --no-strict"
$E --atoms "$WATER" --margin 2.0 --no-strict --out r4.png --json r4.json --npz r4.npz 2>&1 | tail -6
echo "######## R5 two waters 8 A apart via --atoms"
$E --atoms "$WATER; O 8.0 0.0 0.1173; H 8.0 0.7572 -0.4692; H 8.0 -0.7572 -0.4692" --out r5.png --json r5.json 2>&1 | tail -6
echo "######## R6 sodium acetate --force"
$E --smiles "CC(=O)[O-].[Na+]" --force --out r6.png --json r6.json 2>&1 | tail -5
echo "######## R7 acetate default (H3)"
$E --smiles "CC(=O)[O-]" --out r7.png --json r7.json 2>&1 | grep -i "wrote\|scale\|cmap" | tail -4
