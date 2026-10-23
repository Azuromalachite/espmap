"""Review 6, R11 follow-ups: what from_file does with the 2D-header z case,
and record numbering when one SDF record does not parse.
usage: a_r11b.py MASTER_DIR FX_DIR"""
import sys
from pathlib import Path

code, fx = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(code))
import esp_map as em  # noqa: E402
from rdkit import Chem  # noqa: E402

print("=== from_file on the 2D-header drawing with wedge z (no SCF) ===")
for name in ("draw_threonine_2Dheader_zwedge.mol", "draw_chlorobenzene_2Dheader_z0.3.mol"):
    atoms, info = em.from_file(str(fx / name))
    rows = [r.split() for r in atoms.split(";")]
    import numpy as np
    X = np.array([[float(v) for v in r[1:]] for r in rows])
    sym = [r[0] for r in rows]
    heavy = [i for i, s in enumerate(sym) if s != "H"]
    d = [np.linalg.norm(X[i] - X[j]) for k, i in enumerate(heavy) for j in heavy[k + 1:]]
    print(name)
    for k in ("source", "h_source", "geometry", "coords"):
        if k in info:
            print(f"  {k}: {info[k]}")
    print("  info keys:", sorted(info))
    print(f"  shortest heavy-heavy distance {min(d):.3f} A; heavy bond lengths in file ~1.5 (a depiction)")

print("=== record numbering when record 2 of 3 does not parse ===")
base = (code / "tests/data/r11/real_C2F4.mol").read_text().rstrip("\n")
bad = base.replace("  2  4  2  3", "  2  4  3  0")


def titled(block, t):
    lines = block.split("\n"); lines[0] = t
    return "\n".join(lines)


p = fx / "three_records_middle_bad.sdf"
p.write_text("\n$$$$\n".join([titled(base, "pose_1"), titled(bad, "pose_2_bad"),
                               titled(base, "pose_3")]) + "\n$$$$\n")
print("RDKit:", [m.GetProp("_Name") if m is not None else None
                 for m in Chem.SDMolSupplier(str(p), removeHs=False)])
for r in (None, 1, 2, 3):
    m, probs, n = em.read_structure(str(p), record=r)
    print(f"  --record {r}: n_found={n}",
          f"read {m.GetProp('_Name')}" if m is not None else [x[1] for x in probs])
p = fx / "two_records_second_bad.sdf"
p.write_text("\n$$$$\n".join([titled(base, "pose_1"), titled(bad, "pose_2_bad")]) + "\n$$$$\n")
m, probs, n = em.read_structure(str(p))
print(f"two records, second unparseable, no --record: n_found={n}",
      f"read {m.GetProp('_Name')} (no MULTIPLE_RECORDS)" if m is not None else [x[1] for x in probs])
