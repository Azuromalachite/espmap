import os, sys, json, io, contextlib, hashlib, subprocess, ast
import numpy as np
REPO = os.environ["REPO"]; sys.path.insert(0, REPO)
import esp_map as em
from rdkit import Chem
em._quiet_rdkit()
print("== H6 dipole translation invariance, acetate sto-3g")
spec, gi = em.from_smiles("CC(=O)[O-]")
for shift in (0.0, 5.0):
    sp = "; ".join(f"{t.split()[0]} {float(t.split()[1])+shift} {t.split()[2]} {t.split()[3]}" for t in spec.split(";"))
    mol, mf = em.run_scf(sp, basis="sto-3g", charge=-1)
    d = em.dipole_about_nuclear_centre(mol, em.total_dm(mf))
    print(f"  shift {shift} A: |mu| = {d['magnitude_D']:.6f} D  (pyscf origin-0 dip: {np.linalg.norm(mf.dip_moment(verbose=0)):.3f})")
print("== H9 gpu warnings")
buf = io.StringIO()
with contextlib.redirect_stderr(buf):
    try:
        em.run_scf("H 0 0 0; H 0 0 0.74", basis="sto-3g", gpu=True, solvent="water")
    except RuntimeError as e:
        print("  raised:", str(e)[:60])
print("  stderr:", buf.getvalue().strip().replace("\n", " | ")[:300])
print("== selftest with kernel sign flipped inside esp_at_points")
orig = em.esp_at_points
def flipped(mol, dm, pts, mem_budget_gb=0.5):
    V = orig(mol, dm, pts)
    nuc = np.zeros(len(pts))
    for A in range(mol.natm):
        nuc += mol.atom_charge(A) / np.maximum(np.linalg.norm(np.asarray(pts) - mol.atom_coord(A), axis=1), 1e-12)
    return nuc + (nuc - V)            # electronic term added instead of subtracted
em.esp_at_points = flipped
with contextlib.redirect_stdout(io.StringIO()) as o:
    rc = em.selftest()
print("  rc =", rc, [l for l in o.getvalue().splitlines() if "far field" in l])
em.esp_at_points = orig
