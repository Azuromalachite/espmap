"""Review 6, PR #73: the cases e_supplied.py built wrongly or could not reach.
usage: e_supplied2.py PR_DIR FX_DIR"""
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1])
from rdkit import Chem
from rdkit.Chem import AllChem

import esp_map as em

FX = Path(sys.argv[2])
em._quiet_rdkit()


def embed(smi, n=1, seed=7):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    AllChem.EmbedMultipleConfs(m, n, p)
    AllChem.MMFFOptimizeMoleculeConfs(m)
    return m


def one(m, cid):
    out = Chem.Mol(m)
    out.RemoveAllConformers()
    out.AddConformer(Chem.Conformer(m.GetConformer(cid)), assignId=True)
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
        print(f"{label}: REFUSED   {[(c, m[:120]) for _, c, m, _ in e.findings]}")


def relabel(src, smiles_same_order):
    """src's coordinates under another Lewis structure with the SAME atom
    order (checked atom by atom)."""
    t = Chem.AddHs(Chem.MolFromSmiles(smiles_same_order))
    assert [a.GetAtomicNum() for a in t.GetAtoms()] == [a.GetAtomicNum() for a in src.GetAtoms()]
    assert [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in t.GetBonds()] == \
           [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in src.GetBonds()]
    c = Chem.Conformer(src.GetConformer())
    c.Set3D(True)
    t.AddConformer(c, assignId=True)
    return t


print("=== 3'. same nuclei, different electron count, non-planar (ethyl +1 / radical / -1) ===")
base = embed("C[CH2]", 3)                       # the radical's geometry, MMFF
recs = [relabel(one(base, 0), "C[CH2+]"), one(base, 1), relabel(one(base, 2), "C[CH2-]")]
for m in recs:
    print(f"   record: charge {Chem.GetFormalCharge(m):+d}, radicals "
          f"{sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())}, "
          f"kind {em.read_structure(write('tmp.sdf', [m]), record=1)[0].GetProp('_coords')}")
run("ethyl cation, radical, anion", path=write("s3c_ethyl_charge_mix.sdf", recs), ensemble="supplied")
run("ethyl radical, cation, anion", path=write("s3d_ethyl_radical_first.sdf",
                                                [recs[1], recs[0], recs[2]]), ensemble="supplied")
# and a protonation-free closed-shell pair with one heavy-atom charge moved:
# methyl nitrite? keep it to an electron-count change on a real species
nb = embed("C[N+](=O)[O-]", 2)
print("   (nitromethane charge-separated vs a nitro radical anion cannot share H counts; skipped)")

print("=== 4'. isotope labels with the atom order preserved ===")
d3 = embed("[2H]C([2H])([2H])O", 2)
h3 = one(d3, 1)
for a in h3.GetAtoms():
    a.SetIsotope(0)
run("CD3OH then CH3OH", path=write("s4_cd3_ch3.sdf", [one(d3, 0), h3]), ensemble="supplied")
run("CH3OH then CD3OH", path=write("s4b_ch3_cd3.sdf", [h3, one(d3, 0)]), ensemble="supplied")

print("=== 9'. diastereotopic methyls swapped by index ===")
s = embed("C[C@H](O)C(C)C", 2)
a1 = one(s, 0)
sw = one(s, 0)                                   # the SAME conformer, methyl labels swapped
h = {i: [n.GetIdx() for n in sw.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1]
     for i in (4, 5)}
P = sw.GetConformer().GetPositions()
perm = {4: 5, 5: 4, **dict(zip(h[4], h[5])), **dict(zip(h[5], h[4]))}
for i, j in perm.items():
    sw.GetConformer().SetAtomPosition(i, P[j].tolist())


def pro(mol):
    x = Chem.Mol(mol)
    x.GetAtomWithIdx(4).SetIsotope(13)
    Chem.AssignStereochemistryFrom3D(x, replaceExistingTags=True)
    from rdkit.Chem import rdCIPLabeler
    rdCIPLabeler.AssignCIPLabels(x)
    a = x.GetAtomWithIdx(3)
    return a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else "none"


print(f"   C3 label with 13C at atom 4: record 1 {pro(a1)}, record 2 {pro(sw)}")
run("record 2 = record 1 with the two diastereotopic methyls' indices swapped",
    path=write("s9_diastereotopic_swap.sdf", [a1, sw]), ensemble="supplied")
