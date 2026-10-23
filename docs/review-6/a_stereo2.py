"""Review 6: a flat 3D-header drawing whose z noise RDKit reads as chirality.
usage: a_stereo2.py MASTER_DIR FX_DIR"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np

code, fx = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(code))
import esp_map as em  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

for smi, name in (("CC(N)C(=O)O", "alanine"), ("CC(O)C(N)C(=O)O", "threonine")):
    d = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(d)
    for amp in (0.02, 0.05, 0.09):
        tally = Counter()
        example = None
        for seed in range(40):
            rng = np.random.default_rng(seed)
            m = Chem.Mol(d)
            c = m.GetConformer()
            for i in range(m.GetNumAtoms()):
                q = c.GetAtomPosition(i)
                c.SetAtomPosition(i, (q.x, q.y, float(rng.uniform(-amp, amp))))
            blk = Chem.MolToMolBlock(m).split("\n")
            blk[1] = blk[1][:20].ljust(20) + "3D" + blk[1][22:]
            p = fx / f"draw3d_{name}_znoise.mol"
            p.write_text("\n".join(blk))
            probe, probs = em.parse_for_check(path=str(p))
            if probe is None:
                tally["refused at read: " + probs[0][1]] += 1
                continue
            v, f = em.check_structure(Chem.RemoveHs(probe), maps_from=probe, require_stereo=True)
            key = f"{probe.GetProp('_coords')}/{v}"
            tally[key] += 1
            if v != "BLOCK" and example is None:
                example = (seed, Chem.MolToSmiles(probe))
                keep = fx / f"draw3d_{name}_znoise_{amp}_seed{seed}_passes.mol"
                keep.write_text("\n".join(blk))
        print(f"{name:10} z noise +-{amp}: {dict(tally)}  first passing: {example}")
