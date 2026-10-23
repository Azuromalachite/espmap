"""docs/ensemble-spec.md section 8 and review 6's heavy reproductions: the
parts that need an SCF. Written while the shared machine was busy, to run
when it is free. One espmap process per row (the normal CLI, so every row
is an ordinary record), resumable: a row already in
bench/ensemble_dft.jsonl is skipped. One job at a time, 4 threads.

    python bench/ensemble_dft.py --list            # what would run, no SCF
    python bench/ensemble_dft.py genus             # S3: ~7 x 5 min, ~4 GB peak
    python bench/ensemble_dft.py dedup_esp         # dedup at 0.5 vs 0.35 A, small, < 1 GB
    python bench/ensemble_dft.py tco               # S1: diastereomer mix, ~2 x 3 min, ~1 GB
    python bench/ensemble_dft.py bromobenzene      # S25: -0.76 basis shift, ~3 x 3 min, ~1.5 GB

genus        C18 acetate conformers (seed 0xC0FFEE, 100 embedded, as
             bench/conformers.py) whose folded call differs between contact
             rules: 3, 7, 70 (>= 6 bonds folds them, >= 7 does not) and four
             folded only by 5-bond contacts (2, 5, 8, 11). def2-SVP + RI, gas,
             0.26 bohr: the bench conformer protocol. Decides CONTACT_MIN_BONDS
             against surface genus.
dedup_esp    each n-alkane conformer class that a 0.5 A dedup loses and 0.35 A
             keeps, paired with the member that absorbed it: do their site means
             differ by more than the grid row (0.15 kcal/mol)? def2-SVP, gas.
tco          trans-cyclooct-4-enol (R) with each ring helicity, one MMFF
             conformer each, def2-SVP, PCM water: how much a mixed ensemble
             (review 6 S1) moves V_min, V_max and the OH site.
bromobenzene the sigma-hole V_max at def2-SVP, def2-TZVP and def2-TZVPP, gas,
             idealised geometry (--molecule), as error_budget.py measured
             chlorobenzene: the -0.76 kcal/mol (SVP -> TZVPP) in CLAUDE.md and
             README has no bench/ row (S25).
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import esp_map as em  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdMolAlign, rdMolTransforms  # noqa: E402

OUT = HERE / "ensemble_dft.jsonl"
WORK = HERE / "ensemble_dft_work"
SEED = 0xC0FFEE
C18 = "CC(=O)OCC/C=C\\CCCCCCCC/C=C\\CCCC"
GAS_SVP = ["--basis", "def2-svp", "--density-fit", "--spacing", "0.26"]


def done():
    if not OUT.exists():
        return set()
    return {(r["task"], r["key"]) for r in map(json.loads, OUT.read_text().splitlines())}


def write_sdf(path, mol, cid, origin):
    m = Chem.Mol(mol, confId=cid)
    m.SetProp("espmap_geometry_origin", origin)
    w = Chem.SDWriter(str(path))
    w.write(m, confId=cid)
    w.close()


def run_row(task, key, argv, extra=None):
    """One espmap process; appends {task, key, rc, wall, record summary}."""
    if (task, key) in done():
        print(f"{task} {key}: done earlier")
        return
    WORK.mkdir(exist_ok=True)
    js = WORK / f"{task}_{key.replace(' ', '_').replace('/', '-')}.json"
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "4"),
               PYSCF_TMPDIR=str(WORK / "pyscf_tmp"))
    (WORK / "pyscf_tmp").mkdir(exist_ok=True)
    t0 = time.time()
    r = subprocess.run([sys.executable, str(HERE.parent / "esp_map.py"), *argv, "--json", str(js),
                        "--out", str(js.with_suffix(".png"))], capture_output=True, text=True, env=env)
    row = dict(task=task, key=key, rc=r.returncode, wall_s=round(time.time() - t0, 1),
               argv=argv, **(extra or {}))
    if r.returncode == 0:
        rec = json.loads(js.read_text())
        res = rec["results"]
        row.update(analysis_key=rec["analysis_key"], V_min=res["V_min_kcal_mol"],
                   V_max=res["V_max_kcal_mol"], sigma_hole=(res.get("sigma_hole") or {}).get("V_max_kcal_mol")
                   if isinstance(res.get("sigma_hole"), dict) else res.get("sigma_hole"),
                   genus=res["surface"].get("genus"), chi=res["surface"].get("chi"),
                   sites={s["label"]: {rl: (s.get(rl) or {}).get("V_area_mean_kcal_mol")
                                       for rl in em.PARTITION_RULES}
                          for s in (rec.get("sites") or {}).get("reaction_sites") or []},
                   espmap_version=rec["provenance"]["espmap_version"], code=rec["provenance"]["code"])
    else:
        row["stderr_tail"] = "\n".join(r.stderr.strip().splitlines()[-5:])
    with OUT.open("a") as f:
        f.write(json.dumps(row) + "\n")
    print(f"{task} {key}: rc {r.returncode}, {row['wall_s']} s")


def c18():
    mol = Chem.AddHs(Chem.MolFromSmiles(C18))
    em._embed_and_minimise(mol, 100, SEED)
    return mol


def task_genus(listing):
    ids = {3: "folded at >= 6 bonds, not >= 7", 7: "folded at >= 6 bonds, not >= 7",
           70: "folded at >= 6 bonds, not >= 7", 2: "5-bond contacts only",
           5: "5-bond contacts only", 8: "5-bond contacts only", 11: "5-bond contacts only"}
    if listing:
        return [("genus", f"C18 conf {c}", why) for c, why in ids.items()]
    mol = c18()
    WORK.mkdir(exist_ok=True)
    for c, why in ids.items():
        sdf = WORK / f"c18_conf{c}.sdf"
        write_sdf(sdf, mol, c, "etkdg_mmff94s")
        run_row("genus", f"C18 conf {c}", ["--file", str(sdf)] + GAS_SVP, dict(conf_id=c, why=why))


def _chain_class(mol, cid):
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    dm = Chem.GetDistanceMatrix(mol)
    i, j = max(((i, j) for i in heavy for j in heavy), key=lambda p: dm[p])
    path = list(Chem.GetShortestPath(mol, i, j))
    conf = mol.GetConformer(cid)
    s = tuple("a" if abs(d) > 120 else ("+" if d > 0 else "-")
              for d in (rdMolTransforms.GetDihedralDeg(conf, *path[k:k + 4])
                        for k in range(len(path) - 3)))
    flip = {"+": "-", "-": "+", "a": "a"}
    return min({s, s[::-1], tuple(flip[x] for x in s), tuple(flip[x] for x in s[::-1])})


def task_dedup_esp(listing):
    rows = []
    for name, smi in (("n-pentane", "CCCCC"), ("n-hexane", "CCCCCC"), ("n-heptane", "CCCCCCC")):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        *_, conv = em._embed_and_minimise(mol, 100, SEED)
        k35 = em.select_ensemble(mol, conv, len(conv), dedup_rmsd=0.35)["members"]
        k50 = {c for c, _ in em.select_ensemble(mol, conv, len(conv), dedup_rmsd=0.5)["members"]}
        heavy = Chem.RemoveHs(mol)
        for c, _ in k35:
            if c in k50:
                continue
            into = min(k50, key=lambda d: rdMolAlign.GetBestRMS(heavy, heavy, c, d))
            rows.append((name, mol, c, into, _chain_class(mol, c), _chain_class(mol, into)))
    if listing:
        return [("dedup_esp", f"{n} conf {c} (absorbed by {d})", f"{a} vs {b}")
                for n, _, c, d, a, b in rows]
    WORK.mkdir(exist_ok=True)
    for name, mol, c, into, ca, cb in rows:
        for cid in (c, into):
            sdf = WORK / f"{name}_conf{cid}.sdf"
            write_sdf(sdf, mol, cid, "etkdg_mmff94s")
            run_row("dedup_esp", f"{name} conf {cid}", ["--file", str(sdf)] + GAS_SVP,
                    dict(molecule=name, conf_id=cid, pair=[c, into], classes=[list(ca), list(cb)]))


def task_tco(listing):
    smi = "O[C@H]1CC/C=C/CCC1"
    if listing:
        return [("tco", f"(R) helicity {s}", "PCM water def2-SVP") for s in ("+", "-")]
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    *_, conv = em._embed_and_minimise(mol, 30, SEED)
    b = next(x for x in mol.GetBonds() if x.GetBondType() == Chem.BondType.DOUBLE)
    _, quad = em._ring_neighbours(mol, b)
    best = {}
    for c, e in sorted(conv, key=lambda ce: ce[1]):
        s = "+" if rdMolTransforms.GetDihedralDeg(mol.GetConformer(c), *quad) > 0 else "-"
        best.setdefault(s, c)
    WORK.mkdir(exist_ok=True)
    for s, c in best.items():
        sdf = WORK / f"tco_{'p' if s == '+' else 'm'}.sdf"
        write_sdf(sdf, mol, c, "etkdg_mmff94s")
        run_row("tco", f"(R) helicity {s}", ["--file", str(sdf), "--basis", "def2-svp",
                                              "--solvent", "water"], dict(conf_id=c, helicity=s))


def task_bromobenzene(listing):
    rows = [("def2-svp", "def2-SVP"), ("def2-tzvp", "def2-TZVP"), ("def2-tzvpp", "def2-TZVPP")]
    if listing:
        return [("bromobenzene", b, "gas sigma-hole") for b, _ in rows]
    for b, _ in rows:
        run_row("bromobenzene", b, ["--molecule", "bromobenzene", "--basis", b], dict(basis=b))


TASKS = dict(genus=task_genus, dedup_esp=task_dedup_esp, tco=task_tco, bromobenzene=task_bromobenzene)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("task", nargs="?", choices=sorted(TASKS))
    ap.add_argument("--list", action="store_true", help="print the rows each task would run")
    a = ap.parse_args()
    if a.list or not a.task:
        for name, fn in TASKS.items():
            for row in fn(True):
                print(" | ".join(map(str, row)), "| done" if row[:2] in done() else "")
    else:
        TASKS[a.task](False)
