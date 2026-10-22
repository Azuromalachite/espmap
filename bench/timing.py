"""Timing benchmark: what does espmap cost on the molecules we actually want?

Every timing in the repo is from a 6-ring aromatic (141-298 nao). The
roadmap targets C18 pheromones at ~900 nao. This measures, per stage:
embed, SCF, surface, ESP, on a fixed geometry per molecule, with and
without density fitting, so that both cost AND the ESP shift from density
fitting come out of the same run.

    python bench/timing.py                # all runs, cheapest first
    python bench/timing.py --only chloro  # runs whose molecule matches
    python bench/timing.py --list

Each run is a separate process (clean peak memory, a crash loses one row)
and appends one JSON line to bench/results.jsonl. Rows already present
are skipped, so the driver can be stopped and restarted.
"""
import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
RESULTS = HERE / "results.jsonl"

MOLECULES = {
    "chlorobenzene": "Clc1ccccc1",
    "decanoic_acid": "CCCCCCCCCC(=O)O",
    "BDSF": "CCCCCCCCC/C=C\\C(=O)O",                       # cis-2-dodecenoic acid
    "Z3Z13-18OH": "OCC/C=C\\CCCCCCCC/C=C\\CCCC",
    "Z3Z13-18OAc": "CC(=O)OCC/C=C\\CCCCCCCC/C=C\\CCCC",
}

# (molecule, basis, density_fit, threads); threads 0 = all. Cheapest first,
# so an interrupted run still leaves the most informative rows.
ALL = os.cpu_count()
RUNS = [
    ("chlorobenzene", "def2-tzvp", False, 1),   # calibration vs repo's 165 s
    ("chlorobenzene", "def2-tzvp", False, 0),
    ("chlorobenzene", "def2-tzvp", True, 0),
]
for m in ("decanoic_acid", "BDSF"):
    RUNS += [(m, "def2-svp", True, 0), (m, "def2-svp", False, 0),
             (m, "def2-tzvp", True, 0), (m, "def2-tzvp", False, 0)]
RUNS += [
    ("Z3Z13-18OH", "def2-svp", True, 0),
    ("Z3Z13-18OAc", "def2-svp", True, 0),
    ("Z3Z13-18OH", "def2-tzvp", True, 0),
    ("Z3Z13-18OAc", "def2-tzvp", True, 0),
    ("Z3Z13-18OAc", "def2-svp", False, 0),
    ("Z3Z13-18OAc", "def2-tzvp", False, 0),     # the overnight one
]


def key(mol, basis, df, threads):
    return f"{mol}|{basis}|{'df' if df else 'exact'}|{threads or ALL}t"


def one(mol, basis, df, threads):
    """Run a single configuration in THIS process and print one JSON line."""
    from pyscf import lib
    import esp_map as em
    from rdkit import Chem

    lib.num_threads(threads or ALL)
    smi = MOLECULES[mol]
    row = {"key": key(mol, basis, df, threads), "molecule": mol, "smiles": smi,
           "basis": basis, "density_fit": df, "threads": lib.num_threads()}

    # gate before compute (invariant 2), through the one reader (invariant 1)
    probe, probs = em.parse_for_check(smiles=smi, path=None, charge=0)
    verdict, findings = em.check_structure(Chem.RemoveHs(probe), basis, spin=0)
    row["gate"] = verdict
    if verdict == em.BLOCK:
        row["error"] = f"gate BLOCK: {findings}"
        return row

    t = time.perf_counter()
    spec, info = em.from_smiles(smi)
    row["t_embed"] = time.perf_counter() - t
    row["geom_sha"] = hashlib.sha1(spec.encode()).hexdigest()[:10]

    t = time.perf_counter()
    mol_, mf = em.run_scf(spec, basis=basis, density_fit=df)
    row["t_scf"] = time.perf_counter() - t
    row["nao"], row["natm"] = mol_.nao, mol_.natm
    row["e_tot"] = mf.e_tot
    row["scf_cycles"] = getattr(mf, "cycles", None)
    dm = em.total_dm(mf)
    row["dipole_D"] = float(sum(x * x for x in mf.dip_moment(
        unit="Debye", verbose=0)) ** 0.5)

    t = time.perf_counter()
    verts, faces, diag = em.extract_surface(mol_, dm, 0.26)
    row["t_surface"] = time.perf_counter() - t
    row["n_verts"], row["chi"] = diag["n_verts"], diag["chi"]

    t = time.perf_counter()
    V = em.esp_at_points(mol_, dm, verts) * em.HA2KCAL
    row["t_esp"] = time.perf_counter() - t
    row["V_min"], row["V_max"] = float(V.min()), float(V.max())
    return row


def done_keys():
    if not RESULTS.exists():
        return set()
    # rows that errored are retried: an OSError from a full /tmp is not a result
    rows = [json.loads(l) for l in RESULTS.read_text().splitlines() if l.strip()]
    return {r["key"] for r in rows if "error" not in r}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--one", nargs=4, metavar=("MOL", "BASIS", "DF", "THREADS"))
    ap.add_argument("--only", help="substring of molecule name")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.one:
        m, b, df, th = a.one
        try:
            row = one(m, b, df == "1", int(th))
        except Exception as e:                      # record, don't hide
            row = {"key": key(m, b, df == "1", int(th)), "error": repr(e)}
        row["peak_rss_GB"] = resource.getrusage(
            resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2
        print("ROW " + json.dumps(row), flush=True)
        return

    runs = [r for r in RUNS if not a.only or a.only in r[0]]
    if a.list:
        for r in runs:
            print(key(*r))
        return
    done = done_keys()
    for r in runs:
        k = key(*r)
        if k in done:
            print(f"skip  {k}", flush=True)
            continue
        print(f"start {k}  {time.strftime('%H:%M:%S')}", flush=True)
        t = time.perf_counter()
        p = subprocess.run([sys.executable, __file__, "--one", r[0], r[1],
                            "1" if r[2] else "0", str(r[3])],
                           capture_output=True, text=True)
        rows = [l[4:] for l in p.stdout.splitlines() if l.startswith("ROW ")]
        row = json.loads(rows[-1]) if rows else {
            "key": k, "error": f"exit {p.returncode}: {p.stderr[-500:]}"}
        row["t_wall"] = time.perf_counter() - t
        with RESULTS.open("a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"done  {k}  {row['t_wall']:.0f} s"
              + (f"  ERROR {row['error'][:200]}" if "error" in row else
                 f"  nao={row['nao']} scf={row['t_scf']:.0f}s "
                 f"V={row['V_min']:+.2f}..{row['V_max']:+.2f}"), flush=True)


if __name__ == "__main__":
    main()
