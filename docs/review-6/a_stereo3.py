"""Review 6: end to end (no SCF) for the noise-read drawings of a_stereo2.
usage: a_stereo3.py MASTER_DIR FX_DIR"""
import sys
from pathlib import Path

import numpy as np

code, fx = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(code))
import esp_map as em  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

for f in sorted(fx.glob("draw3d_threonine_znoise_*_passes.mol")):
    probe, _ = em.parse_for_check(path=str(f))
    v0, f0 = em.check_structure(Chem.RemoveHs(probe), maps_from=probe)
    atoms, info = em.from_file(str(f))
    ok, computed, why = em.verify_computed(atoms, Chem.RemoveHs(probe))
    print(f.name)
    print("  without --require-stereo: verdict", v0, [c for s, c, *_ in f0])
    print("  coords:", probe.GetProp("_coords"), "| source:", info["source"][:70])
    print("  computed (stereo from the generated 3D):", computed, "verified:", ok)

print("=== drawings with z noise below FLAT_TOL that come out 'geometry' ===")
for smi, name in (("CC(N)C(=O)O", "alanine"), ("CC(O)C(N)C(=O)O", "threonine")):
    d = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(d)
    for seed in range(40):
        rng = np.random.default_rng(seed)
        m = Chem.Mol(d)
        c = m.GetConformer()
        z = rng.uniform(-0.09, 0.09, m.GetNumAtoms())
        for i in range(m.GetNumAtoms()):
            q = c.GetAtomPosition(i)
            c.SetAtomPosition(i, (q.x, q.y, float(z[i])))
        c.Set3D(True)
        X = c.GetPositions(); X = X - X.mean(0)
        dev = np.abs(X @ np.linalg.svd(X)[2][-1]).max()
        kind = em.coordinate_kind(m, "mol")
        if kind == "geometry":
            print(f"  {name} seed {seed}: max|z| {np.abs(z).max():.3f} A, best-plane max dev {dev:.3f} A > FLAT_TOL {em.FLAT_TOL} -> geometry")
