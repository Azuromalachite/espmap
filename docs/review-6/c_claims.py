"""Review 6: recompute numbers quoted in README.md, CLAUDE.md, docs/roadmap-review.md,
docs/ensemble-spec.md from bench/ and the hostile-set manifests, under their protocols.
usage: c_claims.py MASTER"""
import json, re, statistics, sys
from pathlib import Path
R = Path(sys.argv[1]); B = R / "bench"
J = lambda n: [json.loads(l) for l in (B / n).read_text().splitlines() if l.strip()]


def show(label, value, quoted, where):
    print(f"{label:58s} data {value:>24s} | quoted {quoted:>20s} | {where}")


# ---- timings (README speed table; roadmap F1 table)
rows = {r["key"]: r for r in J("results.jsonl") if "error" not in r}
dup = {}
for r in J("results.jsonl"):
    dup.setdefault(r["key"], []).append(r)
print("duplicate keys in results.jsonl:",
      {k: [round(x.get("t_scf", -1)) for x in v] for k, v in dup.items() if len(v) > 1})
clean = {r["key"]: r for r in J("results_clean_rerun.jsonl")}
for k, q, w in [("chlorobenzene|def2-tzvp|exact|16t", "39 s", "README/roadmap"),
                ("chlorobenzene|def2-tzvp|exact|1t", "150 s on 1 core", "roadmap"),
                ("chlorobenzene|def2-tzvp|df|16t", "21 s", "roadmap"),
                ("decanoic_acid|def2-svp|exact|16t", "66 s", "README/roadmap"),
                ("decanoic_acid|def2-svp|df|16t", "48 s", "roadmap"),
                ("decanoic_acid|def2-tzvp|exact|16t", "211 s", "README/roadmap"),
                ("decanoic_acid|def2-tzvp|df|16t", "97 s", "roadmap"),
                ("BDSF|def2-svp|exact|16t", "76 s", "roadmap"), ("BDSF|def2-svp|df|16t", "59 s", "roadmap"),
                ("BDSF|def2-tzvp|exact|16t", "286 s", "roadmap"), ("BDSF|def2-tzvp|df|16t", "129 s", "roadmap"),
                ("Z3Z13-18OH|def2-svp|df|16t", "103 s", "roadmap"),
                ("Z3Z13-18OH|def2-tzvp|df|16t", "813 s", "roadmap"),
                ("Z3Z13-18OAc|def2-svp|exact|16t", "182 s", "README/roadmap"),
                ("Z3Z13-18OAc|def2-svp|df|16t", "146 s", "roadmap"),
                ("Z3Z13-18OAc|def2-tzvp|exact|16t", "875 s", "README/roadmap/CLAUDE")]:
    ts = [round(x["t_scf"]) for x in dup.get(k, []) if "t_scf" in x]
    show(f"t_scf {k}", str(ts), q, w)
show("t_scf C18 OAc tzvp df (clean rerun)",
     f"{clean['Z3Z13-18OAc|def2-tzvp|df|16t']['t_scf']:.0f}", "1,175 s", "README/roadmap/CLAUDE")
nao = {k.split("|")[0]: r["nao"] for k, r in rows.items() if "|def2-tzvp|" in k}
show("nao (TZVP) per molecule", str(nao), "253/492/566/793/898", "README/roadmap")
rat = {k: rows[k.replace("|exact|", "|df|")]["t_scf"] / r["t_scf"] for k, r in rows.items()
       if "|exact|16t" in k and k.replace("|exact|", "|df|") in rows}
tz = {k.split("|")[0]: round(v, 2) for k, v in rat.items() if "tzvp" in k}
show("DF/exact time ratio at TZVP (uses dict-last row)", str(tz), "0.45-0.53 (<=~600)", "README")
sv = {k.split("|")[0]: round(v, 2) for k, v in rat.items() if "svp|" in k}
show("DF/exact time ratio at SVP", str(sv), "saves 20-28 %", "README")
sh = {}
for k, r in rows.items():
    d = k.replace("|exact|", "|df|")
    if "|exact|16t" in k and d in rows:
        sh[k.split("|")[0] + "|" + k.split("|")[1]] = (rows[d]["V_max"] - r["V_max"],
                                                       rows[d]["V_min"] - r["V_min"],
                                                       rows[d]["e_tot"] - r["e_tot"])
show("DF-exact dV_max (always upward?)", str({k: f"{v[0]:+.4f}" for k, v in sh.items()}), "always upward", "README/CLAUDE")
show("DF-exact |dE| range (Ha)",
     f"{min(abs(v[2]) for v in sh.values()):.1e}..{max(abs(v[2]) for v in sh.values()):.1e}",
     "6e-5 to 3e-4", "roadmap F8")
show("DF-exact max |dV|", f"{max(max(abs(v[0]), abs(v[1])) for v in sh.values()):.3f}",
     "<= 0.018 / <= 0.02", "roadmap/README")
show("C18 OAc TZVP DF-exact dV_max", f"{sh['Z3Z13-18OAc|def2-tzvp'][0]:+.3f}", "+0.006", "roadmap F8")
g = json.loads((B / "gpu_l4.json").read_text())
cpu = rows["decanoic_acid|def2-tzvp|df|16t"]["t_scf"]
show("GPU speedup vs results.jsonl CPU row", f"{cpu / g['gpu']['t_scf_s']:.2f} (cpu {cpu:.1f})",
     "4.3x (cpu 97.4)", "README/roadmap")

# ---- chlorobenzene budget prose endpoints
eb = {r["key"]: r for r in J("error_budget.jsonl")}
s = lambda k: eb[k]["sigma_hole"]
show("basis endpoints svp -> tzvpp", f"{s('basis|def2-svp'):+.2f} -> {s('basis|def2-tzvpp'):+.2f}",
     "+6.35 -> +9.90", "README")
fn = [s("basis|def2-tzvp")] + [s(f"functional|{x}") for x in ("pbe0", "m062x", "hf")]
show("functional range", f"{min(fn):+.2f} to {max(fn):+.2f}", "+9.51 to +10.20", "README")
ge = [s("basis|def2-tzvp"), s("geometry|mmff"), s("geometry|dft")]
show("geometry range", f"{min(ge):+.2f} to {max(ge):+.2f}", "+9.54 to +9.94", "README")
show("cones", str(sorted({round(c['sigma_hole'], 2) for c in eb['basis|def2-tzvp']['cones']})),
     "9 cones, all +9.83", "README")
sc = {(r["basis"], r["solvent"]): r for r in J("solvent_convention.jsonl")}
for bs in ("def2-svp", "def2-tzvp"):
    gas = sc[(bs, None)]["solute_only"]["sigma_hole"]
    w, c = sc[(bs, "water")], sc[(bs, "chloroform")]
    show(f"{bs} water: solute shift; +RF shift",
         f"{w['solute_only']['sigma_hole'] - gas:+.2f} ({w['solute_only']['sigma_hole']:+.2f}); "
         f"{w['with_reaction_field']['sigma_hole'] - gas:+.2f} ({w['with_reaction_field']['sigma_hole']:+.2f})",
         "-3.60 (+6.22); +0.05 (+9.88)" if bs == "def2-tzvp" else "-3.10; +0.57", "README")
    show(f"{bs} chloroform: solute shift; +RF shift",
         f"{c['solute_only']['sigma_hole'] - gas:+.2f} ({c['solute_only']['sigma_hole']:+.2f}); "
         f"{c['with_reaction_field']['sigma_hole'] - gas:+.2f} ({c['with_reaction_field']['sigma_hole']:+.2f})",
         "-2.69 (+7.14); -0.34 (+9.49)" if bs == "def2-tzvp" else "-", "README")
show("gas TZVP sigma-hole: error_budget vs solvent_convention",
     f"{s('basis|def2-tzvp'):.4f} vs {sc[('def2-tzvp', None)]['solute_only']['sigma_hole']:.4f}", "+9.83", "consistency")
print("files in bench/ mentioning bromobenzene:",
      [f.name for f in B.rglob("*") if f.is_file() and "bromo" in f.read_text(errors="ignore").lower()])

# ---- anion geometry
an = {r["key"]: r for r in J("anion_geometry.jsonl")}
show("nitrite MMFF N-O", str(an["nitrite|mmff|water"]["bonds"]), "1.42 (DFT 1.25)", "README")
for sp in ("nitrite", "thiocyanate", "enolate"):
    for so in ("water", "gas"):
        d = {f: an[f"{sp}|mmff|{so}"][f] - an[f"{sp}|dft|{so}"][f] for f in ("V_min", "V_max")}
        print(f"   anion geometry {sp:12s} {so:5s} MMFF-DFT dV_min {d['V_min']:+.2f} dV_max {d['V_max']:+.2f}"
              f"  homo mmff {an[f'{sp}|mmff|{so}'].get('homo_eV', 0):+.2f} dft {an[f'{sp}|dft|{so}'].get('homo_eV', 0):+.2f}")

# ---- conformers
cf = [r for r in J("conformers.jsonl") if "conf_id" in r]
ens = J("conformers.jsonl")[0]["ensemble"]
ok = [r for r in cf if "V_min" in r]
f9 = J("f9_recheck.jsonl")
allr = ok + f9
rng = lambda xs: max(xs) - min(xs)
show("ext/all V_min range", f"{rng([r['V_min'] for r in ok]):.2f} / {rng([r['V_min'] for r in allr]):.2f}",
     "1.10 / 5.20", "roadmap F3")
show("ext/all V_max range", f"{rng([r['V_max'] for r in ok]):.2f} / {rng([r['V_max'] for r in allr]):.2f}",
     "1.18 / 2.48; CLAUDE 2.5", "roadmap/CLAUDE")
dp_e = [r["dipole_D"] for r in ok]; dp_a = [r["dipole_D"] for r in allr]
show("dipole range ext (% mean) / all",
     f"{rng(dp_e):.2f} ({100 * rng(dp_e) / statistics.mean(dp_e):.0f} %) / {rng(dp_a):.2f}",
     "0.77 (43 %) / 0.77", "README/roadmap")
ne = [r["area_neg_A2"] for r in ok]; na = [r["area_neg_A2"] for r in allr]
show("neg area ext / all",
     f"{min(ne):.1f}-{max(ne):.1f} ({100 * rng(ne) / statistics.mean(ne):.0f} %) / "
     f"{min(na):.1f}-{max(na):.1f} ({100 * rng(na) / statistics.mean(na):.0f} %)",
     "68.8-78.5 (13 %) / -117.5 (60 %)", "README/roadmap/CLAUDE")
ta = [r["area_A2"] for r in ok]; taa = [r["area_A2"] for r in allr]
show("total area ext / all", f"{min(ta):.0f}-{max(ta):.0f} / {min(taa):.0f}-{max(taa):.0f}",
     "413-424 / 401-424", "roadmap")
show("folded MMFF dE (conf 40, 72)", str([round(r['mmff_dE'], 1) for r in cf if r['conf_id'] in (40, 72)]),
     "4.8 and 6.4", "roadmap")
rest = [r["V_min"] for r in ok] + [x["V_min"] for x in f9 if x["genus"] != 3]
show("V_min genus-3 folded vs the rest",
     f"{[r['V_min'] for r in f9 if r['genus'] == 3][0]:.2f} vs {min(rest):.1f}..{max(rest):.1f}",
     "-35.70 vs -39.8 to -40.9", "roadmap")
show("ensemble MMFF spread; lowest-E rank", f"{ens['mmff_dE_range']:.1f}; {ens['lowest_E_rg_rank']}",
     "10.5; 87 of 100", "spec/roadmap")
print(f"benzoate: -158.7 to -40.7 spans {158.7 - 40.7:.1f} kcal/mol; README says 117")

# ---- hostile-set counts
r11 = json.loads((R / "tests/data/r11/manifest.json").read_text())["cases"]
st = json.loads((R / "tests/data/stereo/manifest.json").read_text())["cases"]
ens_c = json.loads((R / "tests/data/ensemble/manifest.json").read_text())["cases"]
src = (R / "tests/test_r11_set.py").read_text()
np_ = re.search(r"NOT_PREFERRED_NOW = UNSAFE_NOW \| \{(.*?)\n\}", src, re.S).group(1)
show("R11 committed cases / not preferred", f"{len(r11)} / {len(re.findall(r'[.]mol', np_))}",
     "README 100/seven; roadmap 105/8", "README:182-184")
show("stereo committed cases", str(len(st)), "106", "README:215")
show("ensemble hostile cases", str(len(ens_c)), "29", "spec 7")
