"""Review 6: supplied_members cost against the number of records (it re-reads the file per record).
usage: e_scale.py PR_DIR FX_DIR"""
import sys, time
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import esp_map as em
src = (Path(sys.argv[1]) / "tests/data/ensemble/E7_poses_one_ligand.sdf").read_text()
for mult in (1, 2, 5, 10):
    p = Path(sys.argv[2]) / f"scale_{20*mult}.sdf"
    p.write_text(src * mult)
    t = time.perf_counter()
    r = em.plan_ensemble(path=p, ensemble="supplied")
    print(f"{len(r['members'])} records: {time.perf_counter()-t:.1f} s", flush=True)
