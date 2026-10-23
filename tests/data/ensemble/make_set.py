"""Hostile set for ensembles and comparison records (docs/ensemble-spec.md
section 7), written BEFORE any ensemble code. Run from the repository root:

    python tests/data/ensemble/make_set.py

It rewrites tests/data/ensemble/*.sdf and manifest.json deterministically.

Each case is a recipe (what to run) and an expectation:
  - outcome "refuse": the run or comparison must stop with a stated reason,
    before any SCF where the spec says so (`before_compute`);
  - outcome "accept": it must run, and `checks` lists what the record must
    show (keys named as in the spec's section 4 and 6).
`worst` is the outcome the case exists to prevent. `cost` is "fast" when the
check needs no SCF (refusals, fixture identity), "slow" otherwise.
"""
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

OUT = Path(__file__).resolve().parent
for p in OUT.glob("*.sdf"):
    p.unlink()
CASES = []
SEED = 0xC0FFEE
C18 = "CC(=O)OCC/C=C\\CCCCCCCC/C=C\\CCCC"


def case(cid, name, recipe, expect, worst, cost, why):
    assert expect["outcome"] in ("refuse", "accept"), cid
    CASES.append(dict(id=cid, name=name, recipe=recipe, expect=expect, worst=worst,
                      cost=cost, why=why))


def confs(smiles, n, seed=SEED):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=n, params=p))
    AllChem.MMFFOptimizeMoleculeConfs(m, mmffVariant="MMFF94s", maxIters=2000)
    return m, cids


def write_sdf(name, mols):
    w = Chem.SDWriter(str(OUT / name))
    for m, cid, title in mols:
        m = Chem.Mol(m)
        m.SetProp("_Name", title)
        w.write(m, confId=cid)
    w.close()
    return name


# --- supplied-ensemble files ---------------------------------------------------
IBU_S = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"
IBU_R = "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"
m, cids = confs(IBU_S, 20)
f_poses = write_sdf("E7_poses_one_ligand.sdf", [(m, c, f"pose {i + 1}") for i, c in enumerate(cids)])

a, ac = confs("Oc1ccccn1", 3)            # 2-hydroxypyridine
b, bc = confs("O=c1cccc[nH]1", 3)        # 2-pyridone, its tautomer
f_taut = write_sdf("E8_ligand_and_tautomer.sdf",
                   [(a, c, "2-hydroxypyridine") for c in ac] + [(b, c, "2-pyridone") for c in bc])

m, cids = confs(IBU_S, 4)
perm = list(reversed(range(m.GetNumAtoms())))
mp = Chem.RenumberAtoms(m, perm)
f_perm = write_sdf("E9_permuted_atom_order.sdf",
                   [(m, c, f"conf {i + 1}") for i, c in enumerate(cids[:3])]
                   + [(mp, cids[3], "conf 4, atoms renumbered")])

s, sc = confs(IBU_S, 5)
r, rc = confs(IBU_R, 1)
f_inv = write_sdf("E15_one_centre_inverted.sdf",
                  [(s, c, f"(S) conf {i + 1}") for i, c in enumerate(sc)] + [(r, rc[0], "(R)")])

p1, c1 = confs("CCCO", 3)                # propan-1-ol
p2, c2 = confs("CC(C)O", 2)              # propan-2-ol: same formula, other molecule
f_iso = write_sdf("E17_constitutional_isomer_mixed.sdf",
                  [(p1, c, "propan-1-ol") for c in c1] + [(p2, c, "propan-2-ol") for c in c2])

# --- E14: octadecyl acetate, all-anti except one gauche+/gauche- kink ------------
# Saturated, so no cis double bond confounds the kink (a first attempt on the
# Z,Z-diene put an unrelated 3.8 A contact beside its cis alkene).
SAT = "CCCCCCCCCCCCCCCCCCOC(C)=O"
m = Chem.AddHs(Chem.MolFromSmiles(SAT))
p = AllChem.ETKDGv3()
p.randomSeed = SEED
AllChem.EmbedMolecule(m, p)
conf = m.GetConformer()
chain = list(range(0, 18)) + [18, 19]          # C1..C18, O, C(=O)
for i in range(len(chain) - 3):
    rdMolTransforms.SetDihedralDeg(conf, *chain[i:i + 4], 180.0)
k = 7                                          # one g+/g- pair mid-chain
rdMolTransforms.SetDihedralDeg(conf, *chain[k:k + 4], 65.0)
rdMolTransforms.SetDihedralDeg(conf, *chain[k + 1:k + 5], -65.0)
w = Chem.SDWriter(str(OUT / "E14_c18_single_gauche_kink.sdf"))
m.SetProp("_Name", "octadecyl acetate, all-anti except one g+g- kink")
m.SetProp("kink_atoms", " ".join(map(str, chain[k:k + 5])))
m.SetProp("chain_atoms", " ".join(map(str, chain)))
w.write(m)
w.close()

# --- ensemble cases --------------------------------------------------------------
case("E1", "benzene, K = 10", {"smiles": "c1ccccc1", "ensemble": 10},
     {"outcome": "accept", "checks": {"k_requested": 10, "k_obtained": 1}},
     "ten copies averaged as if independent", "slow", "a rigid molecule has one conformer")
case("E2", "cyclohexane", {"smiles": "C1CCCCC1", "ensemble": 10},
     {"outcome": "accept", "checks": {"max_members": 2, "no_duplicates": True}},
     "duplicate chairs counted twice", "slow", "the chair, and a twist-boat only if distinct")
case("E3", "18-crown-6 (macrocycle)", {"smiles": "C1COCCOCCOCCOCCOCCO1", "ensemble": 10},
     {"outcome": "accept", "checks": {"folded_by_genus_alone": 0}},
     "every member called folded because the ring opens a hole", "slow",
     "ring genus is a diagnostic, not the folded rule")
case("E4", "C18 acetate, K = 10", {"smiles": C18, "ensemble": 10},
     {"outcome": "accept", "checks": {"folded_count_min": 1,
                                      "stats_over": ["all", "extended"]}},
     "folded members averaged in without a flag", "slow",
     "the measured site spread is 5.1 with folded members, 2.69 without")
case("E5", "E4 with the folded members' surfaces failing integrity",
     {"smiles": C18, "ensemble": 10, "inject": "surface_failure_on_folded"},
     {"outcome": "accept", "checks": {"excluded_min": 1, "validity.bias": False}},
     "a converged-looking ensemble of only extended members", "slow",
     "the chi = 2 check once excluded exactly the folded C18 members (F9)")
case("E6", "n-butane: gauche+ and gauche- (mirror-related, achiral molecule)",
     {"smiles": "CCCC", "ensemble": 10},
     {"outcome": "accept", "checks": {"max_members": 2, "no_duplicates": True}},
     "symmetry copies inflating n", "slow",
     "anti and one gauche; the mirror gauche is the same conformer for an achiral molecule")
case("E7", "supplied: 20 conformers ('poses') of one ligand",
     {"file": f_poses, "ensemble": "supplied"},
     {"outcome": "accept", "checks": {"k_obtained_max": 20, "protocol.source": "supplied"}},
     "refused, or computed as pose 1 only", "slow", "a docking run is a legitimate ensemble")
case("E8", "supplied: a ligand mixed with its tautomer", {"file": f_taut, "ensemble": "supplied"},
     {"outcome": "refuse", "before_compute": True},
     "a mixed 'ensemble'", "fast", "one ensemble, one molecule")
case("E9", "supplied: one record's atom order permuted", {"file": f_perm, "ensemble": "supplied"},
     {"outcome": "accept", "checks": {"member_atom_order": "verified_mapping"},
      "or": {"outcome": "refuse", "before_compute": True}},
     "site rows from the wrong atoms", "slow",
     "reorder by verified mapping, or refuse; never trust file order across records")
case("E10", "members under two different bases",
     {"members_from": [{"smiles": C18, "basis": "def2-svp"}, {"smiles": C18, "basis": "sto-3g"}]},
     {"outcome": "refuse", "before_compute": False},
     "a mixed-level ensemble", "slow", "every member shares one analysis_key")
case("E11", "--weight dft, one member 5 kcal/mol below the rest",
     {"smiles": "CCCCCCO", "ensemble": 6, "optimize": True, "weight": "dft"},
     {"outcome": "accept", "checks": {"kish_n_reported": True}},
     "a weighted mean presented as an ensemble of K", "slow", "Kish n makes dominance visible")
case("E12", "MMFF energies offered as weights",
     {"smiles": "CCCCCCO", "ensemble": 6, "weight": "mmff"},
     {"outcome": "refuse", "before_compute": True},
     "a near-single-member average from force-field energies", "fast",
     "MMFF orders members, never weights them")
case("E13", "threonine with no stereo, generated", {"smiles": "CC(O)C(N)C(=O)O", "ensemble": 10},
     {"outcome": "refuse", "before_compute": True},
     "an average over diastereomers", "fast", "generated ensembles need stereo fully specified")
case("E14", "octadecyl acetate with a single g+g- kink", {"file": "E14_c18_single_gauche_kink.sdf",
                                                   "classify_only": True},
     {"outcome": "accept", "checks": {"folded": False}},
     "an ordinary kink counted as folding", "fast",
     "an extended chain with one ordinary kink; whether the spec's contact rule "
     "(>= 5 bonds, Bondi + 0.5 A) passes it is what this case tests")
case("E15", "supplied: one centre inverted in one member", {"file": f_inv, "ensemble": "supplied"},
     {"outcome": "refuse", "before_compute": True},
     "a mixed-stereoisomer ensemble", "fast", "members must share one 3D stereo completion")
case("E16", "--weight dft on single points of MMFF geometries",
     {"smiles": "CCCCCCO", "ensemble": 6, "weight": "dft"},
     {"outcome": "refuse", "before_compute": True},
     "weights carrying force-field strain", "fast", "weighting needs optimised members")
case("E17", "supplied: a constitutional isomer mixed in", {"file": f_iso, "ensemble": "supplied"},
     {"outcome": "refuse", "before_compute": True},
     "a record whose identity and members disagree", "fast",
     "one InChIKey per ensemble, checked on every member")

# --- comparison cases ------------------------------------------------------------
ACET = "CC(=O)[O-]"
case("C1", "A compared with itself",
     {"A": {"smiles": "CCCCCCO", "ensemble": 6}, "B": "same_as_A",
      "sites": {"by": "smarts", "pairs": [{"name": "O", "A": "[OX2H]", "B": "[OX2H]"}]}},
     {"outcome": "accept", "checks": {"delta_mean_abs_max": 1e-9, "resolved": False}},
     "resolved", "slow", "no difference is no difference")
case("C2", "A and B at different bases",
     {"A": {"smiles": "CCO", "basis": "def2-svp"}, "B": {"smiles": "CCCO", "basis": "sto-3g"},
      "sites": {"by": "smarts", "pairs": [{"name": "O", "A": "[OX2H]", "B": "[OX2H]"}]}},
     {"outcome": "refuse"}, "a cross-level difference", "slow", "analysis_key fields must match")
case("C3", "SMARTS matching two non-equivalent atoms in B",
     {"A": {"smiles": "CCO"}, "B": {"smiles": "CC(O)CCO"},
      "sites": {"by": "smarts", "pairs": [{"name": "OH", "A": "[OX2H]", "B": "[OX2H]"}]}},
     {"outcome": "refuse"}, "an arbitrary match", "slow",
     "butane-1,3-diol has two different OH groups")
case("C4", "a difference smaller than the partition sensitivity",
     {"A": {"smiles": "CC(C)CO"}, "B": {"smiles": "CCCCO"},
      "sites": {"by": "smarts", "pairs": [{"name": "OH", "A": "[OX2H]", "B": "[OX2H]"}]},
      "inject": "delta_below_partition"},
     {"outcome": "accept", "checks": {"resolved": False}}, "resolved", "slow",
     "the partition term is part of U")
case("C5", "mean difference whose sign flips across 30 % of member pairs",
     {"A": {"smiles": "CCCCCCO", "ensemble": 10}, "B": {"smiles": "CCCCCCCO", "ensemble": 10},
      "sites": {"by": "smarts", "pairs": [{"name": "OH", "A": "[OX2H]", "B": "[OX2H]"}]},
      "inject": "sign_consistency_0.7"},
     {"outcome": "accept", "checks": {"resolved": False}}, "resolved on means alone", "slow",
     "sign_consistency below s_min")
case("C6", "two single-geometry C18 records",
     {"A": {"smiles": C18}, "B": {"smiles": "CC(=O)OCCC/C=C\\CCCCCCC/C=C\\CCCC"},
      "sites": {"by": "smarts", "pairs": [{"name": "C=O", "A": "[CX3]=O", "B": "[CX3]=O"}]}},
     {"outcome": "accept", "checks": {"flags": ["CONFORMERS_NOT_SAMPLED"]}},
     "resolved with no conformer caveat", "slow", "single geometries of flexible molecules")
case("C7", "a non-converged ensemble vs a converged one",
     {"A": {"smiles": C18, "ensemble": 3}, "B": {"smiles": C18, "ensemble": 10},
      "sites": {"by": "smarts", "pairs": [{"name": "C=O", "A": "[CX3]=O", "B": "[CX3]=O"}]}},
     {"outcome": "accept", "checks": {"flags_any": ["NOT_CONVERGED"]}}, "silent", "slow",
     "non-convergence travels with the comparison")
case("C8", "correspondence by atom index",
     {"A": {"smiles": "CCO"}, "B": {"smiles": "CCCO"},
      "sites": {"by": "index", "pairs": [{"name": "O", "A": [2], "B": [3]}]}},
     {"outcome": "refuse"}, "index-based pairing", "fast", "indices carry no chemistry (R5)")
case("C9", "a map pair whose atoms differ in element",
     {"A": {"smiles": "[CH3:1][OH:2]"}, "B": {"smiles": "[CH3:1][NH2:2]"},
      "sites": {"by": "map", "pairs": [{"name": "X", "A": [2], "B": [2]}]}},
     {"outcome": "refuse"}, "a non-homologous 'pair'", "fast", "a map pair asserts homology")
case("C10", "SMARTS matching acetate's two equivalent oxygens",
     {"A": {"smiles": ACET, "solvent": "water"}, "B": {"smiles": "CCC(=O)[O-]", "solvent": "water"},
      "sites": {"by": "smarts", "pairs": [{"name": "O", "A": "[OX1]", "B": "[OX1]"}]}},
     {"outcome": "accept", "checks": {"site_atoms_A": 2, "site_atoms_B": 2}},
     "refused, or one arbitrary oxygen", "slow", "symmetry-equivalent matches are one union site")
case("C11", "equal analysis_key, different espmap versions",
     {"A": {"smiles": "CCO"}, "B": {"smiles": "CCCO"}, "inject": "B_version_differs",
      "sites": {"by": "smarts", "pairs": [{"name": "OH", "A": "[OX2H]", "B": "[OX2H]"}]}},
     {"outcome": "accept", "checks": {"flags": ["VERSIONS_DIFFER"], "versions_recorded": 2}},
     "silent", "slow", "the key says comparable; the flag lets a reader check")
case("C12", "an ensemble vs a single record of the same flexible molecule",
     {"A": {"smiles": C18, "ensemble": 10}, "B": {"smiles": C18},
      "sites": {"by": "smarts", "pairs": [{"name": "C=O", "A": "[CX3]=O", "B": "[CX3]=O"}]}},
     {"outcome": "accept", "checks": {"flags": ["CONFORMERS_NOT_SAMPLED"], "single_side": "B"}},
     "refused, or the single side treated as sampled", "slow", "mixed sources are allowed, flagged")

(OUT / "manifest.json").write_text(json.dumps(
    {"about": __doc__.split("\n\n")[0], "spec": "docs/ensemble-spec.md section 7",
     "cases": CASES}, indent=1) + "\n", encoding="utf-8")
print(len(CASES), "cases")
