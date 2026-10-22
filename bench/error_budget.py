"""The error-budget rows that had no data in the repository (review 3 H8).

Reviews 1 and 2 measured these on the chlorobenzene sigma-hole and wrote the
numbers into prose; nothing in bench/ could reproduce them. This re-measures
each row with the protocol the prose describes, one SCF per row, appending
to error_budget.jsonl (resumable: finished rows are skipped).

  axis        protocol                                          prose -> measured
  basis       def2-SVP / def2-TZVP / def2-TZVPP, B3LYP          +3.55 -> +3.54
  functional  B3LYP / PBE0 / M06-2X / HF at def2-TZVP            0.69 -> 0.69
  grid        spacing 0.26 -> 0.16 bohr at def2-TZVP             0.04 -> 0.04
  region      sigma-hole cone, pr_min x pe_max, 9 combinations   0.00 -> 0.00
  geometry    idealised / SMILES+MMFF94s / B3LYP/def2-SVP opt    0.35 -> 0.40
  stereo      threonine vs allo-threonine, 3 seeds each, SVP     0.9 vs 0.4 -> 1.6-1.9 vs 1.25

Common: B3LYP unless stated, gas phase, iso 0.002, spacing 0.26 bohr,
margin 5.0 bohr; idealised chlorobenzene has C-Cl 1.740 A. The solvent row
is in solvent_convention.jsonl.

    python bench/error_budget.py
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import esp_map as em  # noqa: E402

OUT = HERE / "error_budget.jsonl"
IDEAL = em.monosubstituted_benzene("Cl", 1.740)
THREONINE = {"threonine (2S,3R)": "C[C@@H](O)[C@H](N)C(=O)O",
             "allo-threonine (2S,3S)": "C[C@H](O)[C@H](N)C(=O)O"}
SEEDS = (0xC0FFEE, 1, 2)


def done():
    if not OUT.exists():
        return set()
    return {r["key"] for r in map(json.loads, OUT.read_text().splitlines()) if r}


def i_cl(spec):
    return next(i for i, t in enumerate(spec.split(";")) if t.split()[0] == "Cl")


def measure(spec, basis="def2-tzvp", xc="b3lyp", spacing=0.26, cones=False, sub=True):
    t0 = time.time()
    mol, mf = em.run_scf(spec, basis=basis, xc=xc)
    dm = em.total_dm(mf)
    v, _, diag = em.extract_surface(mol, dm, spacing, iso=em.ISO_BADER, margin=5.0)
    V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
    row = dict(basis=basis, xc=xc, spacing_bohr=spacing, nao=int(mol.nao),
               e_tot_Ha=float(mf.e_tot), n_verts=int(len(v)),
               V_min=float(V.min()), V_max=float(V.max()))
    if sub:
        sh, n = em.sigma_hole(v, V, mol, i_sub=i_cl(spec))
        row.update(sigma_hole=sh, n_cap=n)
    if cones:
        row["cones"] = [dict(pr_min=p, pe_max=e, sigma_hole=s, n_cap=k)
                        for p in (0.0, 0.3, 0.6) for e in (1.0, 1.6, 2.2)
                        for s, k in [em.sigma_hole(v, V, mol, i_sub=i_cl(spec),
                                                   pr_min=p, pe_max=e)]]
    row["t_wall"] = round(time.time() - t0, 1)
    return row


def dft_optimised(spec):
    from pyscf.geomopt.berny_solver import optimize
    _, mf = em.run_scf(spec, basis="def2-svp")
    mol = optimize(mf, maxsteps=50)
    return "; ".join(f"{mol.atom_symbol(i)} " + " ".join(f"{c * em.BOHR:.6f}"
                     for c in mol.atom_coord(i)) for i in range(mol.natm))


def tasks():
    for b in ("def2-svp", "def2-tzvp", "def2-tzvpp"):
        yield f"basis|{b}", "basis", "idealised", lambda b=b: measure(
            IDEAL, basis=b, cones=(b == "def2-tzvp"))
    for xc in ("pbe0", "m062x", "hf"):                 # B3LYP is basis|def2-tzvp
        yield f"functional|{xc}", "functional", "idealised", lambda xc=xc: measure(IDEAL, xc=xc)
    yield "grid|0.16", "grid", "idealised", lambda: measure(IDEAL, spacing=0.16)
    yield "geometry|mmff", "geometry", "SMILES + MMFF94s", lambda: measure(
        em.from_smiles("Clc1ccccc1")[0])
    yield "geometry|dft", "geometry", "B3LYP/def2-SVP optimised", lambda: measure(
        dft_optimised(IDEAL))
    for name, smi in THREONINE.items():
        for s in SEEDS:
            yield (f"stereo|{name}|{s:#x}", "stereo", f"{name}, MMFF seed {s:#x}",
                   lambda smi=smi, s=s: measure(em.from_smiles(smi, seed=s)[0],
                                                basis="def2-svp", sub=False))


def main():
    have = done()
    for key, axis, geometry, run in tasks():
        if key in have:
            continue
        row = dict(key=key, axis=axis, geometry=geometry, **run(),
                   code=em.code_identity()["code"])
        with OUT.open("a") as f:
            f.write(json.dumps(row) + "\n")
        print(key, {k: row[k] for k in ("sigma_hole", "V_min", "V_max") if k in row},
              flush=True)


if __name__ == "__main__":
    main()
