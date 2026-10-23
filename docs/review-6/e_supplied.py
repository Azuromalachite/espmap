"""Review 6, PR #73: supplied-ensemble and generated-ensemble identity attacks.
No SCF. usage: e_supplied.py PR_DIR FX_DIR
Each case builds an SDF the committed hostile set (tests/data/ensemble) does
not have, runs plan_ensemble, and prints ACCEPTED / REFUSED with the codes."""
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1])
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T

import esp_map as em

FX = Path(sys.argv[2])
FX.mkdir(parents=True, exist_ok=True)
em._quiet_rdkit()


def embed(smi, n=1, seed=7):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    AllChem.EmbedMultipleConfs(m, n, p)
    AllChem.MMFFOptimizeMoleculeConfs(m)
    return m


def one(m, cid):
    x = Chem.Mol(m, confId=cid)
    out = Chem.Mol(x)
    out.RemoveAllConformers()
    out.AddConformer(Chem.Conformer(x.GetConformer(cid)), assignId=True)
    return out


def write(name, mols):
    p = FX / name
    w = Chem.SDWriter(str(p))
    for m in mols:
        w.write(m)
    w.close()
    return p


def run(label, **kw):
    try:
        r = em.plan_ensemble(**kw)
        print(f"{label}: ACCEPTED  members={len(r['members'])}  identity={r['identity']['smiles']}")
        return r
    except em.EnsembleRefused as e:
        print(f"{label}: REFUSED   {[(c, m[:110]) for _, c, m, _ in e.findings]}")
    except Exception as e:  # a crash is a finding too
        print(f"{label}: CRASH     {type(e).__name__}: {e}")


def with_coords(template_smiles, src, cid=-1, props=None):
    """A mol from template_smiles (same atom order as src) carrying src's coords."""
    t = Chem.AddHs(Chem.MolFromSmiles(template_smiles))
    assert t.GetNumAtoms() == src.GetNumAtoms(), template_smiles
    c = Chem.Conformer(t.GetNumAtoms())
    pos = src.GetConformer(cid).GetPositions()
    for i in range(t.GetNumAtoms()):
        c.SetAtomPosition(i, pos[i].tolist())
    c.Set3D(True)
    t.AddConformer(c, assignId=True)
    Chem.AssignStereochemistryFrom3D(t)
    return t


print("=== 1. meso / chiral pair (tartaric acid) ===")
meso = embed("O[C@@H](C(=O)O)[C@H](O)C(=O)O", 1)
rr = embed("O[C@@H](C(=O)O)[C@@H](O)C(=O)O", 1)
run("meso then (R,R)", path=write("s1_meso_rr.sdf", [meso, meso, rr]), ensemble="supplied")

print("=== 2. E/Z (pent-2-ene) ===")
E = embed("C/C=C/CC", 1)
Z = embed("C/C=C\\CC", 1)
run("E then Z", path=write("s2_e_z.sdf", [E, E, Z]), ensemble="supplied")

print("=== 3. same nuclei, different charge / electron count: allyl cation, radical, anion ===")
cat = embed("C=C[CH2+]", 2)
rad = with_coords("C=C[CH2]", cat, 1)
ani = with_coords("C=C[CH2-]", cat, 0)
for m, s in ((one(cat, 0), "cation"), (rad, "radical"), (ani, "anion")):
    print(f"   record {s}: formal charge {Chem.GetFormalCharge(m)}, radicals "
          f"{sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())}")
r = run("allyl +1 / radical / -1", path=write("s3_allyl_charge_mix.sdf", [one(cat, 0), rad, ani]),
        ensemble="supplied")
print("   NOTE: plan returns no charge; the members' own charges are",
      [Chem.GetFormalCharge(m) for m in Chem.SDMolSupplier(str(FX / "s3_allyl_charge_mix.sdf"),
                                                             removeHs=False)])
# a neutral analogue: a record that is the triplet (biradical) of a closed-shell molecule
eth = embed("C=C", 2)
trip = with_coords("[CH2][CH2]", eth, 1)
run("ethene then its 1,2-biradical", path=write("s3b_ethene_biradical.sdf", [one(eth, 0), trip]),
    ensemble="supplied")

print("=== 4. isotope labels ===")
d3 = embed("[2H]C([2H])([2H])O", 2)
h3 = with_coords("CO", d3, 1)
run("CD3OH then CH3OH", path=write("s4_cd3_ch3.sdf", [one(d3, 0), h3]), ensemble="supplied")
run("CH3OH then CD3OH", path=write("s4b_ch3_cd3.sdf", [h3, one(d3, 0)]), ensemble="supplied")
# isotope-only stereocentre, both labels: same ESP, different "stereoisomer"
dR = embed("[2H][C@H](C)O", 1)
dS = embed("[2H][C@@H](C)O", 1)
run("(R)-ethanol-1-d then (S)", path=write("s4c_chd.sdf", [dR, dS]), ensemble="supplied")

print("=== 5. H counts that match, on other atoms (4-methylimidazole tautomers) ===")
t1 = embed("Cc1c[nH]cn1", 1)
t2 = embed("Cc1cnc[nH]1", 1)
run("4-Me vs 5-Me imidazole", path=write("s5_imid_tautomers.sdf", [t1, t1, t2]), ensemble="supplied")

print("=== 6. a 2D drawing under a 3D header, with z noise and explicit H ===")
base = embed("CCCO", 2)
flat = Chem.Mol(one(base, 1))
Chem.rdDepictor.Compute2DCoords(flat)
c = flat.GetConformer()
rng = np.random.default_rng(3)
for i in range(flat.GetNumAtoms()):
    p = c.GetAtomPosition(i)
    c.SetAtomPosition(i, (p.x, p.y, float(rng.uniform(-0.3, 0.3))))
c.Set3D(True)
run("record 2 = depiction + z noise", path=write("s6_drawing_3dheader.sdf", [one(base, 0), flat]),
    ensemble="supplied")

print("=== 7. atom maps that disagree between records ===")
mp = embed("[CH3:1][CH2:2][OH:3]", 2)
r1, r2 = one(mp, 0), one(mp, 1)
for a in r2.GetAtoms():                       # swap maps 1 and 3 on record 2
    a.SetAtomMapNum({1: 3, 3: 1}.get(a.GetAtomMapNum(), a.GetAtomMapNum()))
run("maps 1,2,3 then 3,2,1", path=write("s7_maps_swapped.sdf", [r1, r2]), ensemble="supplied")
r3 = one(mp, 1)
for a in r3.GetAtoms():
    a.SetAtomMapNum(5 if a.GetAtomicNum() > 1 else 0)
run("maps then duplicate maps", path=write("s7b_maps_dup.sdf", [r1, r3]), ensemble="supplied")
r4 = one(mp, 1)
for a in r4.GetAtoms():
    a.SetAtomMapNum(0)
run("maps then no maps", path=write("s7c_maps_none.sdf", [r1, r4]), ensemble="supplied")
run("duplicate maps as record 1 (gate sees it)", path=write("s7d_dup_first.sdf", [r3, r1]),
    ensemble="supplied")

print("=== 8. an unparseable record in the middle ===")
good = embed("CCCO", 3)
p8 = write("s8_unparseable_mid.sdf", [one(good, 0), one(good, 1), one(good, 2)])
txt = p8.read_text().split("$$$$\n")
blk = txt[1].split("\n")
# make record 2's first carbon pentavalent: add a bond line C1-C3 is not enough; set a charge
# that breaks valence instead: an M  CHG of -3 on carbon is still parseable, so use a bad bond
for k, line in enumerate(blk):
    if line.startswith("  1  2  1"):
        blk[k] = "  1  2  3  0"                   # C#C with full H count: valence 6
txt[1] = "\n".join(blk)
p8.write_text("$$$$\n".join(txt))
print("   RDKit reads:", [m is not None for m in Chem.SDMolSupplier(str(p8), removeHs=False)])
run("3 records, record 2 unparseable", path=p8, ensemble="supplied")

print("=== 9. diastereotopic methyls swapped by index in one record ===")
s = embed("C[C@H](O)C(C)C", 2)
a1, a2 = one(s, 0), one(s, 1)
m = a2
heavy_h = {i: [n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1]
           for i in (4, 5)}
conf = m.GetConformer()
P = conf.GetPositions()
swap = {4: 5, 5: 4, **dict(zip(heavy_h[4], heavy_h[5])), **dict(zip(heavy_h[5], heavy_h[4]))}
for i, j in swap.items():
    conf.SetAtomPosition(i, P[j].tolist())


def pro(mol):
    """CIP label of C3 when atom 4 is 13C: which diastereotopic methyl atom 4 is."""
    x = Chem.Mol(mol)
    x.GetAtomWithIdx(4).SetIsotope(13)
    Chem.AssignStereochemistryFrom3D(x)
    return x.GetAtomWithIdx(3).GetPropsAsDict().get("_CIPCode", "?")


print(f"   atom 4 is the pro-{pro(a1)} methyl in record 1, pro-{pro(m)} in record 2")
run("record 2 swaps the two methyls' indices", path=write("s9_diastereotopic_swap.sdf", [a1, m]),
    ensemble="supplied")

print("=== 10. stereo elements RDKit does not perceive from 3D ===")


def dihedral(mol, cid, quad):
    return T.GetDihedralDeg(mol.GetConformer(cid), *quad)


# allene axis + a stereocentre: diastereomers (R,Ra) / (R,Sa)
al = embed("C[C@H](O)C=C=CC", 30, seed=11)
heavy = [x.GetIdx() for x in al.GetAtoms()]
sgn = [np.sign(dihedral(al, c.GetId(), (1, 3, 5, 6))) for c in al.GetConformers()]
print(f"   allene axis signs over 30 embeddings of C[C@H](O)C=C=CC: +{sgn.count(1.0)} / -{sgn.count(-1.0)}")
ip = next(i for i, x in enumerate(sgn) if x > 0)
im = next(i for i, x in enumerate(sgn) if x < 0)
run("allene (R,Ra) + (R,Sa)", path=write("s10_allene_diastereomers.sdf",
                                         [one(al, ip), one(al, im)]), ensemble="supplied")
# trans-cyclooctene planar chirality + a centre
tco = embed("O[C@H]1CC/C=C/CCC1", 30, seed=11)
sg = [np.sign(dihedral(tco, c.GetId(), (2, 3, 4, 5))) for c in tco.GetConformers()]
print(f"   TCO helicity signs over 30 embeddings: +{sg.count(1.0)} / -{sg.count(-1.0)}")
ip = next(i for i, x in enumerate(sg) if x > 0)
im = next(i for i, x in enumerate(sg) if x < 0)
run("TCO (R,Rp) + (R,Sp)", path=write("s10b_tco_diastereomers.sdf", [one(tco, ip), one(tco, im)]),
    ensemble="supplied")
# atropisomer + centre: 2'-methyl... a tetra-ortho biaryl with a stereocentre
atr = embed("Cc1cccc(C)c1-c1c(C)cccc1[C@H](C)O", 30, seed=11)
# axis dihedral: C(ortho-Me side) - C - C' - C'(the CH(OH) side)
mt = Chem.MolFromSmarts("[CH3]c1cccc([CH3])c1-c1c([CH3])cccc1[CH](C)O")
mm = atr.GetSubstructMatch(mt)
quad = (mm[1], mm[6], mm[8], mm[14])
sa = [np.sign(dihedral(atr, c.GetId(), quad)) for c in atr.GetConformers()]
print(f"   biaryl axis signs over 30 embeddings: +{sa.count(1.0)} / -{sa.count(-1.0)}")
if 1.0 in sa and -1.0 in sa:
    ip = sa.index(1.0)
    im = sa.index(-1.0)
    run("tetra-ortho biaryl (R,Ra) + (R,Sa)", path=write("s10c_atrop_diastereomers.sdf",
                                                         [one(atr, ip), one(atr, im)]),
        ensemble="supplied")

print("=== 11. GENERATED ensembles from a 3D file: is the file's configuration kept? ===")
thr = embed("C[C@@H](O)[C@H](N)C(=O)O", 1)          # L-threonine (2S,3R)
sdf = write("g1_threonine_3d.sdf", [thr])
xyz = FX / "g1_threonine.xyz"
Chem.MolToXYZFile(thr, str(xyz))
open_smiles = Chem.MolFromSmiles("CC(O)C(N)C(=O)O")


def stereo_mix(plan):
    got = {}
    for mem in plan["members"]:
        ok, comp, why = em.verify_computed(mem["spec"], open_smiles)
        got[comp] = got.get(comp, 0) + 1
    return got


for lab, p in (("3D SDF", sdf), ("XYZ", xyz)):
    probe, _ = em.parse_for_check(path=p)
    print(f"   {lab}: probe SMILES {Chem.MolToSmiles(Chem.RemoveHs(probe))}; "
          f"stereo_left_open={em.stereo_left_open(Chem.RemoveHs(probe))}")
    r = run(f"generated K=10 from threonine {lab}", path=p, ensemble=10)
    if r:
        print(f"   member configurations: {stereo_mix(r)}")

al3 = write("g2_allene_3d.sdf", [one(al, 0)])
r = run("generated K=10 from a 3D allene (R,Ra) SDF", path=al3, ensemble=10)
if r:
    probe, _ = em.parse_for_check(path=al3)
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C=C=CC"))
    signs = []
    for mem in r["members"]:
        X = np.array([[float(v) for v in t.split()[1:]] for t in mem["spec"].split(";")])
        b0, b1, b2 = X[3] - X[1], X[5] - X[3], X[6] - X[5]
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        signs.append(int(np.sign(np.dot(np.cross(n1, n2), b1))))
    print(f"   allene axis signs among members: {signs}")

print("=== 12. GENERATED from SMILES with an unperceived element: TCO ===")
r = run("generated K=10 O[C@H]1CC/C=C/CCC1", smiles="O[C@H]1CC/C=C/CCC1", ensemble=10)
if r:
    signs = []
    for mem in r["members"]:
        X = np.array([[float(v) for v in t.split()[1:]] for t in mem["spec"].split(";")])
        b0, b1, b2 = X[3] - X[2], X[4] - X[3], X[5] - X[4]
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        signs.append(int(np.sign(np.dot(np.cross(n1, n2), b1))))
    print(f"   TCO helicity among members: {signs}")

print("=== 13. tiny molecules ===")
for smi in ("[H][H]", "O", "C"):
    run(f"generated K=3 {smi}", smiles=smi, ensemble=3)
