import os, sys, time, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tree"))
import esp_map as em
from pyscf import lib
print("threads", lib.num_threads())

def sig(spec, **kw):
    t = time.time()
    mol, mf = em.run_scf(spec, **kw)
    dm = em.total_dm(mf)
    v, f, d = em.extract_surface(mol, dm, 0.26)
    V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
    sh, n = em.sigma_hole(v, V, mol, i_sub=6)
    return dict(sh=sh, Vmin=V.min(), Vmax=V.max(), t=time.time() - t), mol, mf, dm, v, V

print("== S7 density fitting on the chlorobenzene sigma-hole (def2-SVP and def2-TZVP)")
spec = em.monosubstituted_benzene("Cl", 1.740)
for basis in ("def2-svp", "def2-tzvp"):
    ex = sig(spec, basis=basis)[0]; df = sig(spec, basis=basis, density_fit=True)[0]
    print(f"  {basis}: exact sh={ex['sh']:+.4f} ({ex['t']:.0f}s)  df sh={df['sh']:+.4f} ({df['t']:.0f}s)  "
          f"diff={df['sh']-ex['sh']:+.4f}  dVmax={df['Vmax']-ex['Vmax']:+.4f}  dVmin={df['Vmin']-ex['Vmin']:+.4f}")

print("== S8 PCM water: size of the solvent reaction-field potential omitted from V (def2-SVP)")
r, mol, mf, dm, v, V = sig(spec, basis="def2-svp", solvent="water")
g, _, _, _, vg, Vg = sig(spec, basis="def2-svp")
print(f"  gas sh={g['sh']:+.3f}   PCM-polarised solute-only sh={r['sh']:+.3f}  (what espmap reports)")
ws = mf.with_solvent
q = None
for k in ("q", "q_sym"):
    try:
        q = np.asarray(ws._intermediates[k]); break
    except Exception as e:
        err = e
if q is not None:
    s = ws.surface["grid_coords"]
    Vrf = np.array([np.sum(q / np.linalg.norm(s - p, axis=1)) for p in v]) * em.HA2KCAL
    co = mol.atom_coords(); X = co[6]; ax = X - co[0]; ax /= np.linalg.norm(ax)
    rel = v - X; pr = rel @ ax; pe = np.linalg.norm(rel - np.outer(pr, ax), axis=1)
    m = (pr > 0.3) & (pe < 1.6)
    i = np.argmax(np.where(m, V, -np.inf))
    print(f"  sum q = {q.sum():+.4f} e (molecule neutral);  reaction-field potential on the surface: "
          f"{Vrf.min():+.2f} .. {Vrf.max():+.2f} kcal/mol;  at the sigma-hole vertex {Vrf[i]:+.2f}")
    print(f"  solute+reaction-field at that vertex: {V[i]+Vrf[i]:+.2f}")
else:
    print("  could not read PCM charges:", repr(err))

print("== S9 dipole_D of an anion depends on where the molecule sits")
spec_a, gi = em.from_smiles("CC(=O)[O-]")
def shift(s, dx):
    return "; ".join(f"{t.split()[0]} {float(t.split()[1])+dx:.6f} {t.split()[2]} {t.split()[3]}" for t in s.split(";"))
for dx in (0.0, 5.0):
    mol, mf = em.run_scf(shift(spec_a, dx), basis="def2-svp", charge=-1)
    print(f"  acetate translated {dx:+.0f} A: e_tot={mf.e_tot:.8f}  |dip|={np.linalg.norm(mf.dip_moment(unit='Debye', verbose=0)):.3f} D")
