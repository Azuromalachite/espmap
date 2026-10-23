"""Review 6: spec section 8's E14 margin claim (closest >=5-bond heavy pair
0.58 A beyond the Bondi sum). usage: c_e14.py MASTER"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import esp_map as em
from rdkit import Chem
m = next(iter(Chem.SDMolSupplier(sys.argv[1] + "/tests/data/ensemble/E14_c18_single_gauche_kink.sdf", removeHs=False)))
T = Chem.GetDistanceMatrix(m); X = m.GetConformer().GetPositions()
best = []
for i in range(m.GetNumAtoms()):
    for j in range(i + 1, m.GetNumAtoms()):
        a, b = m.GetAtomWithIdx(i), m.GetAtomWithIdx(j)
        if a.GetAtomicNum() > 1 and b.GetAtomicNum() > 1 and T[i, j] >= 5:
            d = np.linalg.norm(X[i] - X[j]); s = em.BONDI[a.GetSymbol()] + em.BONDI[b.GetSymbol()]
            best.append((d - s, i, j, int(T[i, j]), d))
best.sort()
for x in best[:3]:
    print(f"pair {x[1]}-{x[2]} bonds {x[3]} d {x[4]:.3f} A, beyond Bondi sum by {x[0]:.3f} A")
