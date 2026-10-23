"""Review 6: |e| vs -e on the rho = 0.002 surface of gas-phase anions whose
HOMO is positive at def2-SVP (review 5 R3). usage: c_alie_unbound.py MASTER"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import esp_map as em
from pyscf.dft import numint
for smi in ("C=C[O-]", "N(=O)[O-]", "CC(=O)[O-]"):
    spec, info = em.from_smiles(smi)
    for solv in (None, "water"):
        mol, mf = em.run_scf(spec, basis="def2-svp", charge=-1, solvent=solv)
        homo = mf.mo_energy[mf.mo_occ > 0].max() * em.HA2EV
        verts, faces, _ = em.extract_surface(mol, em.total_dm(mf), 0.3)
        a = em.alie_at_points(mol, mf, verts) * em.HA2EV
        C, e, n = mf.mo_coeff, mf.mo_energy, mf.mo_occ
        occ = n > 0
        ao = numint.eval_ao(mol, verts); w = (ao @ C[:, occ]) ** 2 * n[occ]
        k = (w @ (-e[occ])) / w.sum(1) * em.HA2EV
        print(f"{smi:12s} {solv or 'gas':5s} HOMO {homo:+.2f} eV | min I |e| {a.min():.3f} at v{a.argmin()} ; "
              f"min I -e {k.min():.3f} at v{k.argmin()} ; max |diff| {np.abs(a-k).max():.3f} eV ; "
              f"HOMO weight at the |e| min {(w[a.argmin(), np.argmax(e[occ])]/w[a.argmin()].sum()) if True else 0:.2f}")
