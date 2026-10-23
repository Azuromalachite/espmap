"""Review 6: DF pairs with the clean re-run superseding (as test_claims does),
bromobenzene nao, test counts. usage: c_claims2.py MASTER"""
import json, sys
from pathlib import Path
B = Path(sys.argv[1]) / "bench"
J = lambda n: [json.loads(l) for l in (B / n).read_text().splitlines() if l.strip()]
rows = {r["key"]: r for r in J("results.jsonl") if "error" not in r}
for r in J("results_clean_rerun.jsonl"):
    rows[r["key"]] = r
for k, ex in rows.items():
    d = k.replace("|exact|", "|df|")
    if "|exact|16t" in k and d in rows:
        df = rows[d]
        print(f"{k:38s} dV_max {df['V_max']-ex['V_max']:+.2e} dV_min {df['V_min']-ex['V_min']:+.2e} dE {df['e_tot']-ex['e_tot']:+.2e} Ha  t ratio {df['t_scf']/ex['t_scf']:.2f}")
from pyscf import gto
br = "Br 0 0 0; C 0 0 1.89; C 1.21 0 2.59; C 1.21 0 3.98; C 0 0 4.68; C -1.21 0 3.98; C -1.21 0 2.59; H 2.15 0 2.05; H 2.15 0 4.52; H 0 0 5.76; H -2.15 0 4.52; H -2.15 0 2.05"
for b in ("def2-svp", "def2-tzvp"):
    print("bromobenzene nao", b, gto.M(atom=br, basis=b, verbose=0).nao)
