source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review5-work/env.sh
cd $HOME
free -m | head -2
/usr/bin/time -f "peak %M kB, %e s" $PY $W/e_scf.py $REPO nu > $W/e_nu.out 2>&1
/usr/bin/time -f "peak %M kB, %e s" $PY $W/e_scf.py $REPO homo > $W/e_homo.out 2>&1
# determinism: the same record three times (4, 4, 1 threads)
for t in 4a 4b 1; do
  n=${t%[ab]}
  OMP_NUM_THREADS=$n OPENBLAS_NUM_THREADS=$n MKL_NUM_THREADS=$n $PY $REPO/esp_map.py --smiles "CC(=O)[O-]" --solvent water --basis def2-svp --size 80 --out $W/out/det_$t.png --json $W/out/det_$t.json --npz $W/out/det_$t.npz > /dev/null 2>&1
done
$PY - <<'PY' > $W/e_det.out 2>&1
import json, numpy as np, os
W = os.environ["W"] + "/out/"
recs = {t: json.load(open(W + f"det_{t}.json")) for t in ("4a", "4b", "1")}
def flat(d, p=""):
    if isinstance(d, dict):
        for k, v in d.items():
            yield from flat(v, f"{p}.{k}")
    elif isinstance(d, list):
        for i, v in enumerate(d):
            yield from flat(v, f"{p}[{i}]")
    else:
        yield p, d
for a, b in (("4a", "4b"), ("4a", "1")):
    fa, fb = dict(flat(recs[a])), dict(flat(recs[b]))
    diffs = [(k, fa[k], fb.get(k)) for k in fa if fa[k] != fb.get(k)]
    num = [abs(x - y) for k, x, y in diffs if isinstance(x, float) and isinstance(y, float)]
    print(f"{a} vs {b}: {len(diffs)} differing leaves; max |numeric diff| {max(num) if num else 0:.3e}")
    for k, x, y in diffs[:8]:
        print("   ", k, x, y)
    za, zb = np.load(W + f"det_{a}.npz"), np.load(W + f"det_{b}.npz")
    print(f"    npz V max |diff| {np.abs(za['V'] - zb['V']).max():.3e}, faces equal {np.array_equal(za['faces'], zb['faces'])}")
PY
cat $W/e_det.out
