"""Review 6, ALIE: independent dense computation, chunk invariance, UKS
alpha+beta vs total_dm, rho->0, on small systems. usage: c_alie_check.py MASTER"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import esp_map as em
from pyscf.dft import numint

def dense(mol, mf, pts):
    C, e, n = (np.asarray(x, float) for x in (mf.mo_coeff, mf.mo_energy, mf.mo_occ))
    if C.ndim == 2: C, e, n = C[None], e[None], n[None]
    ao = numint.eval_ao(mol, pts)
    num = den = 0
    for s in range(len(C)):
        for i in range(C.shape[2]):
            if n[s, i] <= 0: continue
            p2 = (ao @ C[s][:, i]) ** 2
            num = num + n[s, i] * p2 * (-e[s, i])        # Koopmans form, -e
            den = den + n[s, i] * p2
    return num / den, den

rng = np.random.default_rng(0)
cases = [("H2 RKS sto-3g", "H 0 0 0; H 0 0 0.74", "sto-3g", 0, 0),
         ("water RKS def2-svp", "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", "def2-svp", 0, 0),
         ("OH radical UKS def2-svp", "O 0 0 0; H 0 0 0.97", "def2-svp", 0, 1),
         ("water UKS spin0 forced? no: NH2 radical UKS def2-svp", "N 0 0 0.14; H 0 0.80 -0.49; H 0 -0.80 -0.49", "def2-svp", 0, 1)]
for label, atom, basis, q, s in cases:
    mol, mf = em.run_scf(atom, basis=basis, charge=q, spin=s)
    pts = rng.normal(scale=3.0, size=(5003, 3))          # bohr, around the molecule
    a = em.alie_at_points(mol, mf, pts)
    b = em.alie_at_points(mol, mf, pts, mem_budget_gb=1e-9)   # chunk = 16
    c = em.alie_at_points(mol, mf, pts, mem_budget_gb=1e-5)   # chunk ~ odd size
    d, den = dense(mol, mf, pts)
    rho = numint.eval_rho(mol, numint.eval_ao(mol, pts), em.total_dm(mf), xctype="LDA")
    print(f"{label}: nao={mol.nao} E={mf.e_tot:.8f}")
    print(f"  chunk invariance max|default-16|={np.nanmax(abs(a-b)):.2e} max|default-odd|={np.nanmax(abs(a-c)):.2e}")
    print(f"  vs dense (-e) max abs dev={np.nanmax(abs(a-d)):.2e} Ha; min I={np.nanmin(a)*em.HA2EV:.3f} eV")
    print(f"  sum n|psi|^2 vs rho(total_dm) max rel dev={np.max(abs(den-rho)/np.maximum(rho,1e-300)):.2e}")
    e = np.atleast_2d(mf.mo_energy); n = np.atleast_2d(mf.mo_occ)
    print(f"  occupied energies max (should be <0): {e[n>0].max():+.4f} Ha")
    far = em.alie_at_points(mol, mf, np.array([[0, 0, 40.0], [0, 0, 80.0], [0, 0, 400.0]]))
    print(f"  far field 40/80/400 bohr: {far}")
# unbound anion: |e| vs -e
mol, mf = em.run_scf("F 0 0 0", basis="sto-3g", charge=-1)
e = mf.mo_energy[mf.mo_occ > 0]
print(f"F- gas sto-3g occupied energies {np.round(e,4)}")
pts = rng.normal(scale=3.0, size=(200, 3))
a = em.alie_at_points(mol, mf, pts); d, _ = dense(mol, mf, pts)
print(f"  |e| implementation vs Koopmans -e: max dev {np.nanmax(abs(a-d))*em.HA2EV:.3f} eV; mean I(|e|)={np.nanmean(a)*em.HA2EV:.3f}, mean I(-e)={np.nanmean(d)*em.HA2EV:.3f} eV")
