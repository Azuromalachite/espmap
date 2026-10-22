"""Review 5, small SCF probes.
  homo   gas-phase HOMO of the P3 anions (and acetate) at def2-SVP: is the
         'unbound anion' condition, which the record does not carry, reachable
         beyond STO-3G?
  nu     gipf nu with one positive triangle of real area vs none.
    python e_scf.py <repo> homo|nu
"""
import sys
from pathlib import Path

REPO = Path(sys.argv[1])
sys.path.insert(0, str(REPO))
import numpy as np  # noqa: E402
import esp_map as em  # noqa: E402

what = sys.argv[2]
if what == "homo":
    for smi in ["CC(=O)[O-]", "N(=O)[O-]", "[S-]C#N", "C=C[O-]", "[F-]", "[OH-]"]:
        spec, gi = em.from_smiles(smi)
        for solv in (None, "water"):
            mol, mf = em.run_scf(spec, basis="def2-svp", charge=gi["charge"], solvent=solv)
            e = np.asarray(mf.mo_energy)
            homo = float(e[np.asarray(mf.mo_occ) > 0].max()) * 27.211386
            w = em.unbound_anion_warning(mf, gi["charge"], solv)
            print(f"{smi:12} {'PCM water' if solv else 'gas':9} def2-SVP nao {mol.nao:3d} "
                  f"HOMO {homo:+6.2f} eV  warning: {'YES' if w else 'no'}", flush=True)
elif what == "nu":
    from skimage import measure
    g = np.linspace(-2, 2, 41)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    v, f, _, _ = measure.marching_cubes(X ** 2 + Y ** 2 + Z ** 2, level=1.0, spacing=(0.1,) * 3)
    f = f.astype(int)
    a = em.triangle_areas(v, f)
    V = -50.0 - 10 * v[:, 2]
    k = int(np.argmax(a))
    for val in (None, 0.01, 0.5, 5.0):
        V2 = V.copy()
        if val is not None:
            V2[f[k]] = val
        g_ = em.gipf_stats(v, f, V2)
        pos = a[V2[f].mean(axis=1) > 0].sum() / a.sum()
        print(f"positive triangle at {val}: positive area fraction {pos:.2e}; nu={g_['nu']} "
              f"sigma2_pos={g_['sigma2_pos_kcal2_mol2']} V_pos_mean={g_['V_pos_mean_kcal_mol']} "
              f"sign_change_ratio={em.sign_change_ratio(V2):.4f}")
