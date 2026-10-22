import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, "tree/bench")
rows = [json.loads(l) for l in open(f"{B}/results.jsonl") if l.strip()]
ok = {r["key"]: r for r in rows if "error" not in r}
print("rows", len(rows), "ok", len(ok), "errors", len(rows) - len(ok))
print("== density-fit pairs (exact vs df, same basis, 16t)")
for k, r in ok.items():
    if "|exact|16t" in k and k.replace("exact", "df") in ok:
        d = ok[k.replace("exact", "df")]
        print(f"  {k:32s} dVmin={d['V_min']-r['V_min']:+.4f} dVmax={d['V_max']-r['V_max']:+.4f} "
              f"t_scf {r['t_scf']:.0f}->{d['t_scf']:.0f}s ratio {d['t_scf']/r['t_scf']:.2f}  nao {r['nao']}  geom same {d['geom_sha']==r['geom_sha']}")
print("== def2-SVP -> def2-TZVP on V_max / V_min (all molecules with both)")
mols = sorted({r["molecule"] for r in ok.values()})
for m in mols:
    for fit in ("exact", "df"):
        s, t = ok.get(f"{m}|def2-svp|{fit}|16t"), ok.get(f"{m}|def2-tzvp|{fit}|16t")
        if s and t:
            print(f"  {m:14s} {fit:5s} dVmax={t['V_max']-s['V_max']:+.3f}  dVmin={t['V_min']-s['V_min']:+.3f}")
print("== conformers (8 extended)")
c = [json.loads(l) for l in open(f"{B}/conformers.jsonl") if '"conf_id"' in l]
okc = [r for r in c if "error" not in r]
import statistics as st
for f in ("V_min", "V_max", "dipole_D", "area_neg_A2"):
    v = [r[f] for r in okc]
    print(f"  {f:12s} n={len(v)} range={max(v)-min(v):.3f}  /mean={100*(max(v)-min(v))/st.mean(v):.0f}%  /min={100*(max(v)-min(v))/min(v):.0f}%")
f9 = [json.loads(l) for l in open(f"{B}/f9_recheck.jsonl")]
allc = okc + f9
for f in ("V_min", "V_max", "area_neg_A2"):
    v = [r[f] for r in allc]
    print(f"  all 10 {f:12s} range={max(v)-min(v):.3f} /mean={100*(max(v)-min(v))/st.mean(v):.0f}%")
lo = [r for r in okc if r["is_lowest_E"]][0]
print("  folded vs lowest-E: dVmin", [round(r["V_min"]-lo["V_min"], 2) for r in f9], " dVmax", [round(r["V_max"]-lo["V_max"], 2) for r in f9])
