source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review5-work/env.sh
export PYTHONDONTWRITEBYTECODE=1
rm -rf $HOME/p3copy && cp -r $P3 $HOME/p3copy && rm -rf $HOME/p3copy/.git
cd $HOME/p3copy
mkdir -p $W/p3dry
$PY - <<'PY'
import os, subprocess, sys, time
sys.path.insert(0, "fixtures")
import make_records as mk
from pathlib import Path
out = Path(os.environ["W"]) / "p3dry"
for name in ("benzene", "nitrite", "thiocyanate", "enolate"):
    if (out / f"{name}.json").exists():
        continue
    t = time.time()
    p = subprocess.run(["/usr/bin/time", "-f", "%M"] + mk.argv(name, mk.FIXTURES[name], out, "sto-3g"),
                       capture_output=True, text=True)
    print(name, "exit", p.returncode, f"{time.time()-t:.1f}s", "peak kB", p.stderr.strip().splitlines()[-1], flush=True)
    (out / f"{name}.stderr").write_text(p.stderr)
PY
ESPMAP_FIXTURES=$W/p3dry $PY -m pytest -q -p no:cacheprovider tests/test_fixtures.py -k "benzene or nitrite or thiocyanate or enolate" 2>&1 | tail -15
$PY - <<'PY'
import json, os
d = os.environ["W"] + "/p3dry/"
for n in ("benzene", "nitrite", "thiocyanate", "enolate"):
    r = json.load(open(d + n + ".json"))
    s = r["sites"]
    print(n, "| geometry_method", r["key_fields"]["analysis"]["geometry_method"],
          "| key", r["analysis_key"][:10], "| budget", [b["id"] for b in s["budget"]],
          "| code", r["provenance"]["code"], "| conformer_scope", s["conformer_scope"],
          "| labels", [x["label"] for x in s["reaction_sites"]])
    print("    O/C site means nn:", {x["label"]: round(x["nearest_nucleus"]["V_area_mean_kcal_mol"], 3)
                                  for x in s["reaction_sites"]})
PY
