"""Fast (no-SCF) reproductions for review 3. Imports esp_map from the copied tree."""
import os, sys, json, subprocess, difflib, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tree"))
import esp_map as em
from rdkit import Chem
from rdkit.Chem import AllChem
em._quiet_rdkit()
D = os.path.join(HERE, "fx"); os.makedirs(D, exist_ok=True)

print("== E1 planar 3D geometry in the xy plane")
m = Chem.AddHs(Chem.MolFromSmiles("Clc1ccccc1"))
p = AllChem.ETKDGv3(); p.randomSeed = 7
AllChem.EmbedMolecule(m, p); AllChem.MMFFOptimizeMolecule(m)
c = m.GetConformer(); X = c.GetPositions()
# rotate into its principal plane: smallest-variance axis -> z
Xc = X - X.mean(0); _, _, vt = np.linalg.svd(Xc); Xr = Xc @ vt.T
for i in range(m.GetNumAtoms()): c.SetAtomPosition(i, Xr[i].tolist())
print("  max|z| after rotation into plane:", float(np.abs(Xr[:, 2]).max()))
xyz = os.path.join(D, "chlorobenzene_dft_like.xyz")
with open(xyz, "w") as f:
    f.write(f"{m.GetNumAtoms()}\nplanar chlorobenzene, 3D geometry\n")
    for a, r in zip(m.GetAtoms(), Xr): f.write(f"{a.GetSymbol()} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}\n")
molf = os.path.join(D, "chlorobenzene_3d_planar.mol")
Chem.MolToMolFile(m, molf)
print("  mol header line:", open(molf).read().splitlines()[1])
for path in (xyz, molf):
    spec, info = em.from_file(path)
    got = np.array([[float(t) for t in s.split()[1:]] for s in spec.split(";")])
    rms = np.sqrt(((got - got.mean(0)) ** 2).sum(1).mean())
    print(f"  {os.path.basename(path)}: source = {info['source']!r}")
    print(f"     ETKDG used: {info['n_embedded']} conformers; out-of-plane max|z| of computed geometry {np.abs(np.linalg.svd(got-got.mean(0))[2][2] @ (got-got.mean(0)).T).max():.3f} A")
# control: same file tilted 10 degrees
th = np.radians(10); R = np.array([[1,0,0],[0,np.cos(th),-np.sin(th)],[0,np.sin(th),np.cos(th)]])
xyz2 = os.path.join(D, "chlorobenzene_tilted.xyz")
with open(xyz2, "w") as f:
    f.write(f"{m.GetNumAtoms()}\ntilted\n")
    for a, r in zip(m.GetAtoms(), Xr @ R.T): f.write(f"{a.GetSymbol()} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}\n")
print("  tilted 10 deg control:", em.from_file(xyz2)[1]["source"])
print("  examples/water.xyz:", em.from_file(os.path.join(HERE, "tree/examples/water.xyz"))[1]["source"])
wxy = os.path.join(D, "water_xy.xyz")
open(wxy, "w").write("3\nwater in xy plane\nO 0.0 0.117 0.0\nH 0.757 -0.469 0.0\nH -0.757 -0.469 0.0\n")
print("  water rotated into xy plane:", em.from_file(wxy)[1]["source"])

print("\n== E2 PDB without CONECT (explicit H, aromatic ring)")
pdb = os.path.join(D, "chlorobenzene.pdb")
lines = []
for i, (a, r) in enumerate(zip(m.GetAtoms(), Xr @ R.T), 1):
    el = a.GetSymbol(); name = (el + str(i))[:4]
    lines.append(f"HETATM{i:5d} {name:<4s} LIG A   1    {r[0]:8.3f}{r[1]:8.3f}{r[2]:8.3f}  1.00  0.00          {el:>2s}")
open(pdb, "w").write("\n".join(lines) + "\nEND\n")
mol, probs = em.parse_for_check(path=pdb)
if mol is None:
    print("  parse problems:", probs)
else:
    v, f = em.check_structure(Chem.RemoveHs(mol))
    print("  read as:", Chem.MolToSmiles(Chem.RemoveHs(mol)), " verdict:", v)
    for s, code, msg, rem in f: print(f"   [{s}] {code}: {msg[:110]} | fix: {rem[:90]}")

print("\n== E3 approved-vs-computed assertion: what each side is derived from")
smi = "CC(N)C(O)C(=O)O"   # threonine, both centres unspecified
appr = Chem.MolToSmiles(Chem.RemoveHs(em.parse_for_check(smiles=smi)[0]))
seen = {}
for seed in (0xC0FFEE, 1, 2, 3, 4, 5, 6, 7):
    spec, gi = em.from_smiles(smi, seed=seed)
    comp = Chem.CanonSmiles(gi["smiles"])
    mm = Chem.Mol(gi["mol"]); 
    # perceive the stereoisomer actually present in the coordinates handed to PySCF
    best = [c.GetId() for c in mm.GetConformers()]
    xyzs = f"{mm.GetNumAtoms()}\n\n" + "\n".join(s.strip() for s in spec.split(";"))
    g = Chem.MolFromXYZBlock(xyzs)
    from rdkit.Chem import rdDetermineBonds
    rdDetermineBonds.DetermineBonds(g, charge=0)
    Chem.AssignStereochemistryFrom3D(g)
    iso = Chem.MolToSmiles(Chem.RemoveHs(g))
    seen.setdefault(iso, []).append(hex(seed))
    print(f"  seed {seed:#x}: approved={appr}  computed(gi)={comp}  assert passes={comp==appr}  coordinates are: {iso}")
print("  distinct stereoisomers computed across 8 seeds:", len(seen))

print("\n== E4 choose_scale knife-edge (verbatim function)")
for lo in (-0.01, +0.01):
    V = np.linspace(lo, 120.0, 5000)
    vlim, kind, notes = em.choose_scale(V, kind="auto")
    t0 = (0 - vlim[0]) / (vlim[1] - vlim[0])
    print(f"  cation surface {lo:+.2f}..+120: kind={kind:10s} vlim=({vlim[0]:+.2f},{vlim[1]:+.2f})  "
          f"fraction of ramp used = {(V.max()-V.min())/(vlim[1]-vlim[0]):.2f}")
print("  default --cmap for a sequential scale stays", "vik", "(main passes a.cmap; choose_scale only adds a note)")
cm = em.build_maps()["vik"]
V = np.array([-170.0, -125.0, -80.0]); vlim, kind, notes = em.choose_scale(V)
for v in V:
    t = (v - vlim[0]) / (vlim[1] - vlim[0]); r, g, b = cm(t)[:3]
    print(f"  anion V={v:+.0f} kcal/mol -> t={t:.2f} rgb=({r:.2f},{g:.2f},{b:.2f}) {'BLUE (reads as positive)' if b > r else 'red/neutral'}")

print("\n== E5 scale_provenance: comparable ignores everything except the flags")
print("  ", em.scale_provenance("diverging", lim=25.0), "<- for ANY basis/xc/solvent/iso/cmap")

print("\n== E6 notebook series refusal with an all-ion series (logic lifted verbatim)")
surfs = [("acetate", {"V": np.array([-170., -80.])}), ("benzoate", {"V": np.array([-158.7, -40.7])})]
one_sign = [n for n, zz in surfs if not (zz["V"].min() < 0 < zz["V"].max())]
refused = bool(one_sign and len(one_sign) < len(surfs))
print(f"  one_sign={one_sign}  refused={refused}  -> rendered on symmetric +/-SHARED_RANGE diverging scale")

print("\n== E7 notebook generator vs shipped notebook")
out = os.path.join(HERE, "rebuilt.ipynb")
subprocess.run([sys.executable, os.path.join(HERE, "tree/notebooks/build_notebook.py"), out], check=True)
a = open(out, encoding="utf-8").read(); b = open(os.path.join(HERE, "tree/notebooks/espmap_colab.ipynb"), encoding="utf-8").read()
print("  identical bytes:", a == b, " sha1", hashlib.sha1(a.encode()).hexdigest()[:10], hashlib.sha1(b.encode()).hexdigest()[:10])
if a != b:
    for l in list(difflib.unified_diff(b.splitlines(), a.splitlines(), lineterm=""))[:30]: print("   ", l)
nb = json.loads(b)
import ast
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] == "code":
        src = "".join(c["source"])
        try: ast.parse(src)
        except SyntaxError as e: print("  cell", i, "SyntaxError", e)
print("  all code cells parse")

print("\n== E8 euler_characteristic callers")
src = open(os.path.join(HERE, "tree/esp_map.py")).read()
print("  occurrences in esp_map.py:", src.count("euler_characteristic("))
