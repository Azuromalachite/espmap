"""Review 5: real CLI runs (STO-3G, small molecules), one process each, then
every record checked against docs/record-schema.json and for consistency
between validity, sites_reason, structure_verified, identity and the NPZ.
    python e_records.py <repo> <outdir> <fxdir> [case ...]
"""
import json
import subprocess
import sys
import time
from pathlib import Path

REPO, OUT, FX = (Path(x) for x in sys.argv[1:4])
ONLY = sys.argv[4:]
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))
import numpy as np  # noqa: E402
from record_tree import load_schema, tree_problems  # noqa: E402

(FX / "cco_index.sites.json").write_text(json.dumps(
    {"by": "index", "groups": [{"label": "hydroxyl", "atoms": ["2+H"]},
                               {"label": "CH2", "atoms": ["1+H"]}]}))

CASES = {
    "benzene_charge2": ["--smiles", "c1ccccc1", "--charge", "2"],
    "benzene_spin2": ["--smiles", "c1ccccc1", "--spin", "2"],
    "water_open": ["--smiles", "O", "--margin", "2.0", "--no-strict"],
    "methyl_uks": ["--smiles", "[CH3]", "--spin", "1"],
    "forced_salt": ["--smiles", "[Na+].[O-]C=O", "--force"],
    "nitrite_xyz": ["--file", str(FX / "nitrite_mmff.xyz"), "--charge", "-1"],
    "formate_rf": ["--smiles", "[O-]C=O", "--solvent", "water", "--include-reaction-field"],
    "cco_sites": ["--smiles", "CCO", "--sites", str(FX / "cco_index.sites.json")],
    "alh4": ["--smiles", "[AlH4-]", "--solvent", "water"],
    "water_opt": ["--smiles", "O", "--optimize"],
    "fluoride_gas": ["--smiles", "[F-]"],
    "h2": ["--smiles", "[H][H]"],
    "cco": ["--smiles", "CCO"],
}


def bad_const(c):
    raise ValueError(f"non-JSON constant {c}")


rows = {}
schema = load_schema()
for name, args in CASES.items():
    if ONLY and name not in ONLY:
        continue
    js, npz = OUT / f"{name}.json", OUT / f"{name}.npz"
    cmd = [sys.executable, str(REPO / "esp_map.py"), *args, "--basis", "sto-3g",
           "--size", "80", "--out", str(OUT / f"{name}.png"), "--json", str(js),
           "--npz", str(npz)]
    t = time.time()
    p = subprocess.run(["/usr/bin/time", "-f", "%M", *cmd], capture_output=True, text=True)
    rss = p.stderr.strip().splitlines()[-1]
    (OUT / f"{name}.stderr").write_text(p.stderr)
    (OUT / f"{name}.stdout").write_text(p.stdout)
    print(f"\n=== {name}: {' '.join(args)}\n    exit {p.returncode}, {time.time() - t:.1f} s, "
          f"peak {int(rss) / 1024:.0f} MB" if rss.isdigit() else f"exit {p.returncode}")
    warn = [ln for ln in p.stderr.splitlines() if "WARNING" in ln]
    if warn:
        print("    stderr WARNING lines:", warn[:3])
    if not js.exists():
        print("    no JSON;", p.stderr.strip().splitlines()[-2:])
        continue
    rec = json.loads(js.read_text(), parse_constant=bad_const)
    probs = tree_problems(schema["record"], rec)
    print("    schema problems:", probs or "none")
    inp, s = rec["input"], rec["sites"]
    print(f"    validity {rec['validity']}  sites_reason {rec['sites_reason']!r}")
    print(f"    gate {inp['gate']['verdict'] if inp['gate'] else None} forced "
          f"{inp['gate']['forced'] if inp['gate'] else None}; approved {inp['approved_smiles']!r} "
          f"computed {inp['computed_smiles']!r} verified {inp['structure_verified']} "
          f"charge {inp['charge']} spin {inp['spin']}")
    if s:
        print(f"    identity {s['identity']}  partition {s['partition']}")
        print(f"    budget ids {[b['id'] for b in s['budget']]}")
        print(f"    sites {[(x['label'], x['source']) for x in s['reaction_sites']]}")
    print(f"    geometry_method {rec['key_fields']['analysis']['geometry_method']}  "
          f"analysis_key {rec['analysis_key'][:10]}  optimised_from {rec['geometry']['optimised_from_sha1']}")
    print(f"    V {rec['results']['V_min_kcal_mol']:+.2f} .. {rec['results']['V_max_kcal_mol']:+.2f}; "
          f"gipf nu {rec['results']['gipf']['nu']}; surface genus {rec['results']['surface']['genus']}")
    if npz.exists():
        z = np.load(npz)
        tri = [k for k in z.files if k.startswith("tri_atom_")]
        print(f"    npz validity {json.loads(str(z['validity']))}  tri_atom arrays {tri}")
        extra = sorted(set(z.files) ^ set(schema["npz"]))
        if extra:
            print(f"    npz files differing from schema['npz']: {extra}")
    rows[name] = dict(validity=rec["validity"], schema_problems=probs,
                      analysis_key=rec["analysis_key"])
(OUT / "summary.json").write_text(json.dumps(rows, indent=1))
