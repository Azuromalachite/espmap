"""Does the chlorobenzene "geometry source" row (0.40 kcal/mol) transfer to small anions?

The row was measured on chlorobenzene (idealised / MMFF / DFT geometry).
For small anions MMFF is far off: an unpublished pilot found nitrite N-O
1.424 A against 1.258 from DFT, and thiocyanate S-C 1.798 against 1.660.
This measures what that does to the surface potential (review 3 H8).

Per species: the MMFF geometry espmap builds from SMILES, and a B3LYP/
def2-TZVP geometry optimised in PCM water from it (pyberny, the --optimize
path). Both evaluated at B3LYP/def2-TZVP, in PCM water (the default
polarised-solute convention) and in the gas phase. Gas-phase anions can be
unbound in this basis; the HOMO energy is recorded so those rows can be
read with that in mind. Resumable: finished rows are skipped.

    python bench/anion_geometry.py
"""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import esp_map as em  # noqa: E402

OUT = HERE / "anion_geometry.jsonl"
SPECIES = {"nitrite": "N(=O)[O-]", "thiocyanate": "[S-]C#N", "enolate": "C=C[O-]"}


def done():
    if not OUT.exists():
        return {}
    rows = [json.loads(l) for l in OUT.read_text().splitlines() if l.strip()]
    return {r["key"]: r for r in rows}


def heavy_bonds(spec):
    """Distances (A) between bonded heavy atoms, by covalent-radius test."""
    from rdkit.Chem import GetPeriodicTable
    pt = GetPeriodicTable()
    at = [(t.split()[0], np.array([float(x) for x in t.split()[1:]])) for t in spec.split(";")]
    out = {}
    for (i, (a, p)), (j, (b, q)) in itertools.combinations(enumerate(at), 2):
        if "H" in (a, b):
            continue
        d = float(np.linalg.norm(p - q))
        if d < 1.25 * (pt.GetRcovalent(a) + pt.GetRcovalent(b)):
            out[f"{a}{i}-{b}{j}"] = round(d, 4)
    return out


def optimise_in_water(spec, charge):
    from pyscf.geomopt.berny_solver import optimize
    _, mf = em.run_scf(spec, basis="def2-tzvp", charge=charge, solvent="water")
    mol = optimize(mf, maxsteps=50)
    return "; ".join(f"{mol.atom_symbol(i)} " + " ".join(f"{c * em.BOHR:.6f}"
                     for c in mol.atom_coord(i)) for i in range(mol.natm))


def measure(spec, charge, solvent):
    t0 = time.time()
    mol, mf = em.run_scf(spec, basis="def2-tzvp", charge=charge, solvent=solvent)
    dm = em.total_dm(mf)
    v, _, _ = em.extract_surface(mol, dm, 0.26, iso=em.ISO_BADER, margin=5.0)
    V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
    homo = float(np.max(mf.mo_energy[mf.mo_occ > 0])) * 27.211386
    return dict(e_tot_Ha=float(mf.e_tot), homo_eV=homo, n_verts=int(len(v)),
                V_min=float(V.min()), V_max=float(V.max()),
                t_wall=round(time.time() - t0, 1))


def main():
    have = done()
    for name, smi in SPECIES.items():
        mmff, info = em.from_smiles(smi)
        charge = info["charge"]
        gkey = f"{name}|dft-geometry"
        if gkey not in have:
            t0 = time.time()
            row = dict(key=gkey, species=name, smiles=smi, charge=charge,
                       method="B3LYP/def2-TZVP, PCM water, pyberny from the MMFF geometry",
                       spec=optimise_in_water(mmff, charge),
                       t_wall=round(time.time() - t0, 1), code=em.code_identity()["code"])
            with OUT.open("a") as f:
                f.write(json.dumps(row) + "\n")
            have[gkey] = row
        geoms = {"mmff": mmff, "dft": have[gkey]["spec"]}
        for g, solvent in itertools.product(geoms, ("water", None)):
            key = f"{name}|{g}|{solvent or 'gas'}"
            if key in have:
                continue
            row = dict(key=key, species=name, smiles=smi, charge=charge, geometry=g,
                       solvent=solvent or "gas", bonds=heavy_bonds(geoms[g]),
                       **measure(geoms[g], charge, solvent),
                       code=em.code_identity()["code"])
            with OUT.open("a") as f:
                f.write(json.dumps(row) + "\n")
            print(key, row["bonds"], round(row["V_min"], 2), round(row["V_max"], 2),
                  flush=True)


if __name__ == "__main__":
    main()
