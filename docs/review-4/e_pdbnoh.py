import os, sys
sys.path.insert(0, os.environ["REPO"]); import esp_map as em
from rdkit import Chem
em._quiet_rdkit()
src = os.environ["REPO"] + "/docs/review-3/fx/chlorobenzene.pdb"
lines = [l for l in open(src) if not (l.startswith(("ATOM", "HETATM")) and l[76:78].strip() == "H")]
p = os.environ["W"] + "/fx/chlorobenzene_noH.pdb"; open(p, "w").writelines(lines)
for f in (src, p):
    m, probs, n = em.read_structure(f)
    if m is None: print(os.path.basename(f), "->", [x[:2] for x in probs]); continue
    v, fi = em.check_structure(Chem.RemoveHs(m), maps_from=m)
    print(os.path.basename(f), "->", Chem.MolToSmiles(Chem.RemoveHs(m)), v, [x[1] for x in fi])
    try:
        spec, gi = em.from_file(f); print("   from_file:", gi["formula"], gi["source"])
    except ValueError as e: print("   from_file refused:", str(e)[:90])
