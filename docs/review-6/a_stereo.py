"""Review 6, --require-stereo attacks on master (no SCF).
usage: a_stereo.py MASTER_DIR FX_DIR"""
import sys
from pathlib import Path

import numpy as np

code, fx = Path(sys.argv[1]), Path(sys.argv[2])
fx.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(code))
import esp_map as em  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem, rdMolTransforms  # noqa: E402


def gate(smiles=None, path=None):
    probe, probs = em.parse_for_check(smiles=smiles, path=path)
    if probe is None:
        return "PARSE " + ",".join(p[1] for p in probs), [], None
    v, f = em.check_structure(Chem.RemoveHs(probe), maps_from=probe, smiles=smiles,
                              require_stereo=True)
    codes = sorted({c for s, c, *_ in f if s in ("BLOCK", "WARN")})
    return v, codes, probe


def show(name, smiles, expect, why):
    v, codes, probe = gate(smiles=smiles)
    left = em.stereo_left_open(Chem.RemoveHs(probe)) if probe is not None else None
    tag = "ok" if (v == "BLOCK") == (expect == "block") else (
        "UNSAFE-PASS" if expect == "block" else "OVERBLOCK")
    print(f"{tag:11} {name:38} verdict={v:6} want={expect:5} open={left} | {why}")


print("=== A. stereo elements the rule cannot see (should block under the flag) ===")
show("TCO-ol (E)-cyclooct-4-en-1-ol, C1 set", "O[C@H]1CC/C=C/CCC1", "block",
     "planar chirality of the trans ring alkene + C1: two diastereomers (axial/equatorial)")
show("trans-cyclooctene", "C1CCC/C=C/CC1", "block", "planar-chiral; enantiomers only")
show("Curran anilide, o-tBu", "CC(=O)N(C)c1ccccc1C(C)(C)C", "block",
     "Ar-N(amide) axis, a stable atropisomer class [I]; axis is not biaryl")
show("o-tBu anilide + stereocentre", "C[C@H](O)C(=O)N(C)c1ccccc1C(C)(C)C", "block",
     "axis + centre: diastereomers")
show("bridged biaryl, dibenzoxepine 1,11-Me2", "Cc1cccc2c1-c1c(C)cccc1COC2", "block",
     "biaryl axis inside a 7-ring: IsInRing() skips it [I barrier]")
show("1,12-dimethylbenzo[c]phenanthrene", "Cc1ccc2ccc3ccc4cccc(C)c4c3c2c1", "block",
     "[4]helicene made stable by fjord methyls [I]; rule needs >= 5 rings")
show("4-methylcyclohexylidene acetic acid", "CC1CCC(=CC(=O)O)CC1", "block",
     "axially chiral alkylidenecycloalkane (resolved 1909) [I]")
show("4-Me-cyclohexylidene + centre", "C[C@H](O)C=C1CCC(C)CC1", "block", "axis + centre")
show("P-stereogenic phosphine", "CCP(C)c1ccccc1", "block", "P inversion ~30+ kcal/mol [I]")
show("sulfoxide", "CS(=O)c1ccccc1", "block", "S stereocentre, stable")
show("sulfonium", "CC[S+](C)c1ccccc1", "block", "stable S+ centre")
show("Troger's base", "Cc1ccc2c(c1)CN1Cc3cc(C)ccc3N(C2)C1", "block", "bridgehead N pair, chiral")
show("cyclic allene 1,2-cyclononadiene", "C1=C=CCCCCCC1", "block", "axially chiral cyclic allene")
print("=== B. elements that should NOT block ===")
show("acyclic tertiary amine", "CCN(C)CCC", "pass", "N inverts fast")
show("2,2',6-trifluorobiphenyl", "Fc1cccc(F)c1-c1ccccc1F", "pass", "ortho score 3 of F only; F small [I]")
show("2,2'-diiodo... control 2,2'-dimethylbiphenyl", "Cc1ccccc1-c1ccccc1C", "pass", "rotates")
show("CHD label only", "[2H][C@@H](O)C", "pass", "isotope-only centre (spec says excluded)")
show("CHD unspecified", "[2H]C(O)C", "pass", "isotope-only centre unspecified")

print("=== C. TCO-ol: which diastereomer does the embedder make? ===")
smi = "O[C@H]1CC/C=C/CCC1"
mol = Chem.AddHs(Chem.MolFromSmiles(smi))
cids = AllChem.EmbedMultipleConfs(mol, 30, randomSeed=0xC0FFEE)
res = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s", maxIters=5000)
ring = [a.GetIdx() for a in mol.GetAtoms()]
# ring atoms in SMILES order: O0 C1 C2 C3 C4 C5 C6 C7 C8 ; trans C4=C5
approval = Chem.MolFromSmiles(smi)
signs = []
for cid, (ok, e) in zip(cids, res):
    conf = mol.GetConformer(cid)
    dih = rdMolTransforms.GetDihedralDeg(conf, 3, 4, 5, 6)
    spec = em._mol_to_spec(mol, cid)
    v, computed, why = em.verify_computed(spec, approval)
    signs.append((np.sign(dih), round(e, 2), v, computed))
pos = [s for s in signs if s[0] > 0]
neg = [s for s in signs if s[0] < 0]
print(f"30 ETKDG+MMFF conformers: dihedral C3-C4=C5-C6 > 0: {len(pos)}, < 0: {len(neg)}")
print("  verify_computed accepts all:", all(s[2] for s in signs), "; computed SMILES set:", {s[3] for s in signs})
if pos and neg:
    print(f"  lowest MMFF energy, + helicity {min(s[1] for s in pos):.2f}, - helicity {min(s[1] for s in neg):.2f} kcal/mol")
best, ff, spread, n_emb, n_conv = em._embed_and_minimise(Chem.AddHs(Chem.MolFromSmiles(smi)), 20, 0xC0FFEE)
print("  (from_smiles keeps the lowest of 20; its helicity is decided by which embeddings converge)")
# ensure +/- are not one diastereomer in two conformations: compare the best of each by RMSD
if pos and neg:
    ip = cids[[i for i, s in enumerate(signs) if s[0] > 0][0]]
    ineg = cids[[i for i, s in enumerate(signs) if s[0] < 0][0]]
    h = Chem.RemoveHs(mol)
    print(f"  heavy-atom best RMSD between a + and a - conformer: {AllChem.GetBestRMS(h, h, ip, ineg):.2f} A")

print("=== D. 3D-header drawing with z noise: chirality read from noise? ===")
d = Chem.MolFromSmiles("CC(O)C(N)C(=O)O")
AllChem.Compute2DCoords(d)
rng = np.random.default_rng(3)
for amp in (0.0, 0.001, 0.01, 0.05):
    m = Chem.Mol(d)
    c = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        q = c.GetAtomPosition(i)
        c.SetAtomPosition(i, (q.x, q.y, float(rng.uniform(-amp, amp))))
    blk = Chem.MolToMolBlock(m).split("\n")
    blk[1] = blk[1][:20].ljust(20) + "3D" + blk[1][22:]
    p = fx / f"draw3d_threonine_znoise_{amp}.mol"
    p.write_text("\n".join(blk))
    v, codes, probe = gate(path=str(p))
    tags = [str(a.GetChiralTag()) for a in probe.GetAtoms() if a.GetIdx() in (1, 3)] if probe is not None else None
    print(f"  z noise +-{amp}: coords={probe.GetProp('_coords') if probe is not None else None} "
          f"chiral tags on C1,C3={tags} verdict={v} codes={codes}")
