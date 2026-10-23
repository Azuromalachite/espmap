"""How much do symmetry-equivalent sites disagree from mesh noise alone?

Review 5 R7: the P3 fixture tests asserted |O1 - O2| < 0.1 kcal/mol on
nitrite and a C-site range < 0.5 on benzene, from one or two draws of a
quantity that depends on how the molecule sits in the marching-cubes grid
(0.012-0.119 at STO-3G/def2-SVP over eight draws). This measures the
distribution at the fixtures' own settings: the fixture geometry under
N_ROT random rigid rotations, B3LYP/def2-TZVP in PCM water, 0.26 bohr.
A rotation changes nothing physical, so every spread here is noise, and
tests/test_fixtures.py takes its tolerances from these rows.

    python bench/symmetry_noise.py      # ~20 min, 4 threads; appends, resumable
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
import esp_map as em  # noqa: E402

OUT = HERE / "symmetry_noise.jsonl"
N_ROT, SEED = 10, 20261023
CASES = {"nitrite": (ROOT / "fixtures" / "geometries" / "nitrite.sdf", -1, "O"),
         "benzene": (ROOT / "bench" / "fixtures" / "benzene.sdf", 0, "C")}


def random_rotation(rng):
    """Uniform on SO(3): QR of a Gaussian matrix, sign-fixed."""
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    return q if np.linalg.det(q) > 0 else -q


def spread(path, charge, element, R):
    from rdkit import Chem
    mol = Chem.MolFromMolFile(str(path), removeHs=False)
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    X = mol.GetConformer().GetPositions() @ R.T
    spec = "; ".join(f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(sym, X))
    m, mf = em.run_scf(spec, basis="def2-tzvp", charge=charge, solvent="water")
    dm = em.total_dm(mf)
    verts, faces, diag = em.extract_surface(m, dm, 0.26)
    V = em.esp_at_points(m, dm, verts) * em.HA2KCAL
    rep = em.site_report(verts * em.BOHR, faces, V, m.atom_coords() * em.BOHR, sym)
    means = [s["nearest_nucleus"]["V_area_mean_kcal_mol"] for s in rep["reaction_sites"]
             if sym[s["atoms"][0]] == element]
    return float(max(means) - min(means)), int(diag["n_verts"])


def main():
    done = ({json.loads(l)["key"] for l in OUT.read_text().splitlines()}
            if OUT.exists() else set())
    rng = np.random.default_rng(SEED)
    rots = [random_rotation(rng) for _ in range(N_ROT)]     # same draws on every run
    for name, (path, charge, element) in CASES.items():
        for k, R in enumerate(rots):
            key = f"{name}|rot{k}"
            if key in done:
                continue
            d, nv = spread(path, charge, element, R)
            row = dict(key=key, molecule=name, rotation=k, seed=SEED, element=element,
                       site_mean_spread=d, n_verts=nv, method="B3LYP/def2-TZVP, PCM water",
                       spacing_bohr=0.26, code=em.code_identity()["code"])
            with OUT.open("a") as f:
                f.write(json.dumps(row) + "\n")
            print(f"{name:8s} rot {k}  spread {d:.4f}  n_verts {nv}", flush=True)


if __name__ == "__main__":
    main()
