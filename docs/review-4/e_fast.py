"""Review 4 fast reproductions (no SCF)."""
import os, sys, json, io, contextlib
import numpy as np
REPO = os.environ["REPO"]; sys.path.insert(0, REPO)
import esp_map as em
from rdkit import Chem
from rdkit.Chem import AllChem
em._quiet_rdkit()
FX = os.path.join(REPO, "docs/review-3/fx")
OUT = os.environ["W"] + "/fx"; os.makedirs(OUT, exist_ok=True)

print("== H1: chlorobenzene.pdb without CONECT")
m, probs, n = em.read_structure(os.path.join(FX, "chlorobenzene.pdb"))
print("  read as", Chem.MolToSmiles(Chem.RemoveHs(m)), probs)

print("== H2: coordinates vs approval")
for appr, geo in [("CCO","CCO"),("CC(=O)[O-]","CC(=O)[O-]"),("CCO","CC=O"),("C/C=C\C","C/C=C/C")]:
    spec, gi = em.from_smiles(geo)
    ok, comp, why = em.verify_computed(spec, Chem.MolFromSmiles(appr))
    print(f"  approval {appr:12s} geometry {geo:12s} -> ok={ok} computed={comp} {why}")

print("== H5: coordinate_kind on review-3 fixtures")
for f in sorted(os.listdir(FX)):
    if f.endswith(".pdb"): continue
    mm, pr, _ = em.read_structure(os.path.join(FX, f))
    print(f"  {f:28s}", mm.GetProp('_coords') if mm is not None else pr[0][1])

print("== H11: sign_change_ratio")
for lo in (-0.01*120, -0.04*120, -0.06*120, -0.5*120, -0.01):
    V = np.array([lo, 120.0])
    vl, kind, _ = em.choose_scale(V)
    print(f"  V {lo:+.2f}..+120 ratio {em.sign_change_ratio(V):.3f} -> {kind}")

print("== N1: site_parents with no heavy atom (H2)")
try:
    print(em.site_parents(np.array([[0,0,0],[0,0,0.74]]), ["H","H"]))
except Exception as e:
    print("  RAISES", type(e).__name__, e)

print("== N2: validity of an ungated --atoms input")
print("  validity(gate=None, verified=None, surface ok) ->", em.validity(None, None, {"problems": []}))

print("== N3: --sites names on an unmapped input")
spec, gi = em.from_smiles("CCO")
rows = [t.split() for t in spec.split(";")]
sym = [r[0] for r in rows]; X = np.array([[float(v) for v in r[1:]] for r in rows])
print("  symbols", sym)
for g in ([{"label":"C1-H (meant map 1)","atoms":["1+H"]}], [{"label":"an H","atoms":[5]}], [{"label":"H+H","atoms":["5+H"]}]):
    print("  ", g, "->", em.resolve_site_groups(g, sym, X, None))

print("== N4: analysis_key and basis spelling")
from pyscf import gto
for b in ("def2-tzvp", "def2tzvp", "DEF2_TZVP"):
    k = em.record_key(em.analysis_fields("b3lyp", b, "vacuum", 0.002, 0.26, False))
    mm = gto.M(atom="O 0 0 0; H 0 0 0.96; H 0.93 0 -0.24", basis=b, verbose=0)
    print(f"  basis {b:10s} nao={mm.nao} analysis_key={k[:10]} small-basis warning: {b.lower() in em.SMALL_BASIS}")
for b in ("def2-svp", "def2svp", "def2_svp"):
    print(f"  basis {b:10s} SMALL_BASIS warning fires: {b.lower() in em.SMALL_BASIS}")

print("== N5: interleaved-H SDF: 'index' vs heavy-atom index")
mh = Chem.AddHs(Chem.MolFromSmiles("CCO"))
AllChem.EmbedMolecule(mh, randomSeed=7)
# renumber so each H follows its heavy atom (as many writers do)
order = []
for a in mh.GetAtoms():
    if a.GetAtomicNum() > 1:
        order.append(a.GetIdx()); order += [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum()==1]
mi = Chem.RenumberAtoms(mh, order)
p = os.path.join(OUT, "ethanol_interleaved.sdf"); w = Chem.SDWriter(p); w.write(mi); w.close()
spec, gi = em.from_file(p)
sy = [t.split()[0] for t in spec.split(";")]
print("  spec order", sy)
print("  O is RDKit/record index", sy.index("O"), "but heavy-atom index 2; default label would be", f"O{sy.index('O')}")

print("== N6: per-atom partition delta vs per-site (ethanol-like blob is synthetic; see SCF run)")

print("== N7: _area_quantile is a lower step quantile")
v = np.array([0., 10.]); w = np.array([1., 1.])
print("  p05/p50/p95 of two equal triangles at 0 and 10:", [em._area_quantile(v, w, q) for q in (0.05, 0.5, 0.95)])
v = np.array([0., 10.]); w = np.array([0.05, 0.95])
print("  p05 with weights .05/.95 (cum exactly 0.05):", em._area_quantile(v, w, 0.05))
