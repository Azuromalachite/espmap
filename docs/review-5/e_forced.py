"""Review 5: do validity, sites_reason, structure_verified and identity agree
on the forced paths?
  mapping   verify_mapping fault-injected to fail, --force (CCO, STO-3G)
  glycine   gas-phase zwitterion --optimize --force (STO-3G): N1 end to end
    python e_forced.py <repo> <outdir> mapping|glycine
"""
import json
import sys
from pathlib import Path

REPO, OUT, what = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
sys.path.insert(0, str(REPO))
import esp_map as em  # noqa: E402

js = OUT / f"forced_{what}.json"
if what == "mapping":
    em.verify_mapping = lambda spec, probe, reference=False, tol=1e-3: (False, "INJECTED: atom 1 swapped")
    argv = ["--smiles", "CCO", "--force"]
else:
    argv = ["--smiles", "[NH3+]CC(=O)[O-]", "--optimize", "--force"]
sys.argv = ["esp_map.py", *argv, "--basis", "sto-3g", "--size", "80",
            "--out", str(OUT / f"forced_{what}.png"), "--json", str(js)]
em.main()
r = json.loads(js.read_text())
print("\nRECORD", what)
print("  validity", r["validity"])
print("  sites_reason", r["sites_reason"])
print("  gate", r["input"]["gate"]["verdict"], "forced", r["input"]["gate"]["forced"])
print("  approved_smiles", r["input"]["approved_smiles"], "| computed_smiles", r["input"]["computed_smiles"],
      "| structure_verified", r["input"]["structure_verified"])
print("  geometry_method", r["key_fields"]["analysis"]["geometry_method"],
      "| optimised_from", r["geometry"]["optimised_from_sha1"])
print("  any field recording the mapping result:",
      [k for k in json.dumps(r).split('"') if "mapp" in k.lower()])
print("  V", r["results"]["V_min_kcal_mol"], r["results"]["V_max_kcal_mol"])
