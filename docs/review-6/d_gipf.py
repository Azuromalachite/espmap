"""R14 as fixed: gipf_stats either side of GIPF_MIN_AREA_FRACTION on a
synthetic surface of 1000 equal, separate triangles. usage: d_gipf.py CODE_DIR"""
import sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import esp_map as em
rng = np.random.default_rng(0)
n = 1000
P = np.zeros((3 * n, 3)); faces = np.arange(3 * n).reshape(n, 3)
for t in range(n):
    P[3*t:3*t+3] = [[t, 0, 0], [t + 1, 0, 0], [t, 1, 0]]
for npos in (49, 50, 51):
    V = np.repeat(rng.normal(-10, 3, n), 3)
    V[:3 * npos] = np.repeat(rng.normal(+10, 3, npos), 3)
    g = em.gipf_stats(P, faces, V)
    print(f"positive area {npos/n:.3f}: sigma2_pos {g['sigma2_pos_kcal2_mol2']}, "
          f"sigma2_neg {g['sigma2_neg_kcal2_mol2']:.3f}, sigma2_tot {g['sigma2_tot_kcal2_mol2']:.3f}, "
          f"nu {g['nu']}, reason {g['reason']!r}")
