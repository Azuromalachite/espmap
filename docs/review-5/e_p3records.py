"""Review 5: the real P3 fixture records (23dede9), read-only, no SCF.
    python e_p3records.py <p3_copy> <repo>
"""
import json
import sys
from pathlib import Path

P3, REPO = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(P3))
sys.path.insert(0, str(P3 / "tests"))
import numpy as np  # noqa: E402
import esp_map as em  # noqa: E402
from record_tree import load_schema, tree_problems  # noqa: E402

D = P3 / "fixtures" / "records"
schema = load_schema()
names = ["benzene", "nitrite", "thiocyanate", "enolate", "Z3Z13-18OAc_extended", "Z3Z13-18OAc_folded"]


def bad(c):
    raise ValueError(c)


recs = {n: json.loads((D / f"{n}.json").read_text(), parse_constant=bad) for n in names}
print("provenance (code, espmap_version, pyscf, rdkit, numpy, threads):")
for n, r in recs.items():
    p = r["provenance"]
    print(f"  {n:22} {p['code']:28} {p['espmap_version']:7} {p['pyscf']} {p['rdkit']} {p['numpy']} {p['threads']}")
print()
for n, r in recs.items():
    s = r["sites"]
    probs = tree_problems(schema["record"], r)
    z = np.load(D / f"{n}.npz")
    print(f"== {n}: schema problems {probs or 'none'}; npz files == schema: {sorted(z.files) == sorted(schema['npz'])}")
    print(f"   validity {r['validity']}  verified {r['input']['structure_verified']}  gate {r['input']['gate']['verdict']} "
          f"{[f[1] for f in r['input']['gate']['findings']]}")
    ka = r["key_fields"]["analysis"]
    print(f"   key {r['analysis_key'][:10]} recomputes {em.record_key(ka) == r['analysis_key']}; "
          f"geometry_method {ka['geometry_method']}; basis {ka['basis']} solvent {ka['solvent']} df {ka['density_fit']}")
    print(f"   settings.geometry: {r['settings']['geometry']}")
    print(f"   V {r['results']['V_min_kcal_mol']:+.3f} .. {r['results']['V_max_kcal_mol']:+.3f}; genus "
          f"{r['results']['surface']['genus']}; nu {r['results']['gipf']['nu']}; scale {r['scale']['kind']}")
    print(f"   budget {[b['id'] for b in s['budget']]}")
    labs = [(x["label"], x["atoms"]) for x in s["reaction_sites"]]
    amap = {a["index"]: a["atom_map"] for a in s["atoms"]}
    print(f"   sites ({len(labs)}): {[(l, [amap[i] for i in a if amap[i]]) for l, a in labs][:6]} ...")
    from collections import Counter
    c = Counter(i for _, a in labs for i in a)
    over = {i: k for i, k in c.items() if k > 1}
    total = r["results"]["surface_stats"]["area_A2"]
    ssum = sum(x["nearest_nucleus"]["area_A2"] for x in s["reaction_sites"])
    print(f"   atoms in >1 site: {over or 'none'}; sum(site areas)/surface = {ssum / total:.4f}")
    buried = [a["index"] for a in s["atoms"] if a["element"] != "H" and a["nearest_nucleus"]["area_A2"] == 0]
    print(f"   heavy atoms buried (nearest_nucleus): {len(buried)}; identity {s['identity']}")

# the author's symmetry claims
b = {x["label"]: x["nearest_nucleus"]["V_area_mean_kcal_mol"] for x in recs["benzene"]["sites"]["reaction_sites"]}
o = {x["label"]: x["nearest_nucleus"]["V_area_mean_kcal_mol"] for x in recs["nitrite"]["sites"]["reaction_sites"]}
print(f"\nbenzene C-site range {max(b.values()) - min(b.values()):.4f} (claim 0.075); "
      f"nitrite |O1-O2| {abs(o['O1'] - o['O2']):.4f} (claim 0.007)")
bb = {x["label"]: x["radius_scaled"]["V_area_mean_kcal_mol"] for x in recs["benzene"]["sites"]["reaction_sites"]}
oo = {x["label"]: x["radius_scaled"]["V_area_mean_kcal_mol"] for x in recs["nitrite"]["sites"]["reaction_sites"]}
print(f"radius_scaled: benzene {max(bb.values()) - min(bb.values()):.4f}; nitrite {abs(oo['O1'] - oo['O2']):.4f}")
keys = {n: r["analysis_key"][:10] for n, r in recs.items()}
print("analysis_keys:", keys)

# installed vs source version
from importlib import metadata  # noqa: E402
import importlib.util  # noqa: E402
print("\ninstalled espmap metadata version:", metadata.version("espmap"))
d = metadata.distribution("espmap")
print("  dist location:", d._path if hasattr(d, "_path") else "?")
print("  esp_map imported from:", em.__file__)
print("  pyproject (P3):", [ln for ln in (P3 / "pyproject.toml").read_text().splitlines() if ln.startswith("version")])
print("  pyproject (main repo):", [ln for ln in (REPO / "pyproject.toml").read_text().splitlines()
                                   if ln.startswith("version")])
print("  code_identity() here:", em.code_identity())
