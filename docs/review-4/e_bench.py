import json, statistics, os
R=os.environ["REPO"]+"/bench/"
rows=[json.loads(l) for l in open(R+"site_partition.jsonl") if l.strip()]
print("n rows", len(rows), [r["key"] for r in rows])
c18=[r for r in rows if r["part"]=="c18"]
print("c18 n", len(c18), "genus", [r["genus"] for r in c18], "conf", [r["conf_id"] for r in c18])
for r in rows:
    if r["part"] in ("grid","symmetry"):
        nn=r["sites"]["nearest_nucleus"]; rs=r["sites"]["radius_scaled"]
        print(r["key"], r["basis"], r["solvent"], "V", round(r["V_min"],2), round(r["V_max"],2))
        for s in nn:
            a=nn[s]; b=rs[s]
            print("   ", s, "area %.3f/%.3f mean %s/%s p05 %s vmin %s" % (a["area_A2"], b["area_A2"],
                  round(a.get("V_area_mean",float('nan')),3), round(b.get("V_area_mean",float('nan')),3), round(a.get("V_p05",float('nan')),3), round(a.get("V_min",float('nan')),3)))
# grid area %
grid={}
for r in rows:
    if r["part"]=="grid": grid.setdefault(r["molecule"],{})[r["spacing_bohr"]]=r["sites"]["nearest_nucleus"]
mx=0
for m,v in grid.items():
    for s in v[0.26]:
        a0=v[0.26][s]["area_A2"]; a1=v[0.16][s]["area_A2"]
        if a0>0:
            rel=abs(a1-a0)/a0; 
            if rel>mx: mx=rel; who=(m,s,a0,a1)
print("max area rel change 0.26->0.16", mx, who)
# conformer spread per site (all c18) and extended only
nn=lambda r: r["sites"]["nearest_nucleus"]
sites=list(nn(c18[0]))
sp=[max(nn(r)[s]["V_area_mean"] for r in c18)-min(nn(r)[s]["V_area_mean"] for r in c18) for s in sites]
ext=[r for r in c18 if r["genus"]==0]
sp_e=[max(nn(r)[s]["V_area_mean"] for r in ext)-min(nn(r)[s]["V_area_mean"] for r in ext) for s in sites]
print("c18 conformer site-mean spread median all %.2f max %.2f ; extended-only (n=%d) median %.2f max %.2f"%(statistics.median(sp),max(sp),len(ext),statistics.median(sp_e),max(sp_e)))
# zero-area sites
for r in c18[:1]:
    z=[s for s in nn(r) if nn(r)[s]["area_A2"]<1.0]
    print("tiny sites", [(s, round(nn(r)[s]["area_A2"],3)) for s in nn(r) if nn(r)[s]["area_A2"]<3])
