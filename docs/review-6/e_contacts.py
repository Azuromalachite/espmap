"""Review 6: which atom pairs make the C18 'folded' calls (PR #73 contact rule). usage: e_contacts.py PR_DIR"""
import sys
from collections import Counter
sys.path.insert(0, sys.argv[1])
from rdkit import Chem
import esp_map as em
em._quiet_rdkit()
smi = "CC(=O)OCC/C=C\CCCCCCCC/C=C\CCCC"
mol = Chem.AddHs(Chem.MolFromSmiles(smi))
*_, conv = em._embed_and_minimise(mol, em.N_EMBED, 0xC0FFEE)
topo = Chem.GetDistanceMatrix(mol)
dbl = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 2 and b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C"]
print("C=C bonds:", dbl)
def spans_cis(i, j):
    path = Chem.GetShortestPath(mol, i, j)
    return any(a in path and b in path for a, b in dbl)
for cid in (57, 72, 40):
    for i, j, d, lim in em.intramolecular_contacts(mol, cid):
        print(f"conf {cid}: atoms {i}-{j} ({mol.GetAtomWithIdx(i).GetSymbol()}{mol.GetAtomWithIdx(j).GetSymbol()}), "
              f"{int(topo[i, j])} bonds, d {d:.2f} < {lim:.2f}, path crosses a C=C: {spans_cis(i, j)}")
c = Counter()
for cid, _ in conv:
    for i, j, d, lim in em.intramolecular_contacts(mol, cid):
        c[(int(topo[i, j]), spans_cis(i, j))] += 1
print("contacts over 100 conformers by (bonds apart, crosses C=C):", dict(sorted(c.items())))
