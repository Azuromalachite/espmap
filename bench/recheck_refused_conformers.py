"""F9 validation: recompute the two conformers the chi = 2 check refused.

bench/conformers.jsonl records conf 40 (chi = -2) and conf 72 (chi = -4) of
Z3Z13-18OAc as refused with no box-face clipping. If they are closed surfaces
of higher genus (no open or non-manifold edges), the old check was wrong to
refuse them and mesh_topology should pass them. Same embedding, seed, basis
and density fitting as bench/conformers.py, so the geometries are identical.

    python bench/recheck_refused_conformers.py   # appends to bench/f9_recheck.jsonl
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import esp_map as em  # noqa: E402
from conformers import MOLECULES, surface_stats  # noqa: E402

spec0, info = em.from_smiles(MOLECULES["Z3Z13-18OAc"], n_confs=100)
mol = info["mol"]
for cid in (40, 72):
    t = time.perf_counter()
    spec = em._mol_to_spec(mol, cid)
    m, mf = em.run_scf(spec, basis="def2-svp", density_fit=True)
    dm = em.total_dm(mf)
    verts, faces, diag = em.extract_surface(m, dm, 0.26)      # strict: raises if refused
    V = em.esp_at_points(m, dm, verts) * em.HA2KCAL
    row = dict(molecule="Z3Z13-18OAc", conf_id=cid, basis="def2-svp", density_fit=True,
               **{k: diag[k] for k in ("chi", "components", "open_edges",
                                       "nonmanifold_edges", "closed", "genus", "face_max_rho")},
               V_min=float(V.min()), V_max=float(V.max()),
               dipole_D=float(np.linalg.norm(mf.dip_moment(unit="Debye", verbose=0))),
               **surface_stats(verts, faces, V), t_wall=time.perf_counter() - t)
    with open(HERE / "f9_recheck.jsonl", "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
