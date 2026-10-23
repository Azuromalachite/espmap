"""Review 6: the C18 'extended' site spread (README/CLAUDE.md 2.69, max 5.16) under the
genus split the bench used vs the PR #73 contact rule (which also calls conf 57 folded).
No SCF: bench/site_partition.jsonl rows. usage: e_extended.py PR_DIR"""
import json, statistics, sys
from pathlib import Path
rows = [json.loads(l) for l in (Path(sys.argv[1]) / "bench/site_partition.jsonl").read_text().splitlines() if l.strip()]
c18 = {int(r["key"].split("|")[-1]): r for r in rows if r["part"] == "c18"}
nn = lambda r: r["sites"]["nearest_nucleus"]
sites = list(nn(next(iter(c18.values()))))
def spread(ids):
    s = [max(nn(c18[i])[k]["V_area_mean"] for i in ids) - min(nn(c18[i])[k]["V_area_mean"] for i in ids) for k in sites]
    return statistics.median(s), max(s)
genus_ext = [i for i, r in c18.items() if r["genus"] == 0]
contact_ext = [i for i in genus_ext if i != 57]          # PR rule: 40, 57, 72 folded (e_select_b.txt)
for name, ids in (("all 10", list(c18)), ("genus 0 (bench, quoted)", genus_ext), ("PR contact rule", contact_ext)):
    m, x = spread(ids)
    print(f"{name:26s} n={len(ids):2d}  median site spread {m:.2f}  max {x:.2f}")
