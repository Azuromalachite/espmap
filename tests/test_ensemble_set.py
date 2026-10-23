"""Ensemble and comparison hostile set (docs/ensemble-spec.md section 7),
committed before any ensemble code.

The fast tests check the manifest and that every fixture file is what its
case says (a tautomer mix really is two tautomers, an inverted member really
inverts one centre, the kinked C18 really has no contact under the spec's
rule). The behaviour tests are xfail(strict=True) until the ensemble code
exists: a PR implementing a case removes it from NOT_YET, and may not edit an
expectation to pass.
"""
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import esp_map as em

_DIR = Path(__file__).resolve().parent / "data" / "ensemble"
_CASES = json.loads((_DIR / "manifest.json").read_text(encoding="utf-8"))["cases"]
_BY_ID = {c["id"]: c for c in _CASES}

# every case is unimplemented; implementation PRs remove entries
NOT_YET = {c["id"] for c in _CASES}


def _records(name):
    return [m for m in Chem.SDMolSupplier(str(_DIR / name), removeHs=False)]


def _key(m, stereo=True):
    """Identity as the spec checks it: canonical SMILES (isomeric, with stereo
    from 3D). Not the InChIKey: standard InChI merges mobile-H tautomers, so
    2-hydroxypyridine and 2-pyridone share one key (found by test_E8)."""
    m = Chem.RemoveHs(Chem.Mol(m))
    if stereo:
        Chem.AssignStereochemistryFrom3D(m)
    return Chem.MolToSmiles(m, isomericSmiles=stereo)


def test_the_manifest_is_complete_and_consistent():
    ids = [c["id"] for c in _CASES]
    assert len(ids) == len(set(ids)) == 29
    assert {i for i in ids if i[0] == "E"} == {f"E{i}" for i in range(1, 18)}
    assert {i for i in ids if i[0] == "C"} == {f"C{i}" for i in range(1, 13)}
    files = {c["recipe"]["file"] for c in _CASES if "file" in c["recipe"]}
    assert files == {p.name for p in _DIR.glob("*.sdf")}
    for c in _CASES:
        assert c["expect"]["outcome"] in ("refuse", "accept") and c["worst"], c["id"]
        assert c["cost"] in ("fast", "slow"), c["id"]


def test_E7_is_one_ligand_in_one_configuration():
    keys = {_key(m) for m in _records(_BY_ID["E7"]["recipe"]["file"])}
    assert len(_records(_BY_ID["E7"]["recipe"]["file"])) == 20 and len(keys) == 1


def test_E8_really_mixes_two_tautomers():
    ms = _records(_BY_ID["E8"]["recipe"]["file"])
    assert len({_key(m) for m in ms}) == 2
    # and the standard InChIKey cannot tell them apart: identity must not use it
    assert len({Chem.MolToInchiKey(m) for m in ms}) == 1
    assert len({rdMolDescriptors.CalcMolFormula(m) for m in ms}) == 1


def test_E9_is_one_molecule_with_one_record_reordered():
    ms = _records(_BY_ID["E9"]["recipe"]["file"])
    assert len({_key(m) for m in ms}) == 1
    orders = {tuple(a.GetSymbol() for a in m.GetAtoms()) for m in ms}
    assert len(orders) == 2


def test_E15_inverts_exactly_one_member():
    ms = _records(_BY_ID["E15"]["recipe"]["file"])
    keys = [_key(m) for m in ms]
    assert len(set(keys)) == 2 and sorted(keys.count(k) for k in set(keys)) == [1, 5]
    assert len({_key(m, stereo=False) for m in ms}) == 1     # same connectivity


def test_E17_mixes_constitutional_isomers():
    ms = _records(_BY_ID["E17"]["recipe"]["file"])
    assert len({_key(m, stereo=False) for m in ms}) == 2
    assert len({rdMolDescriptors.CalcMolFormula(m) for m in ms}) == 1


def test_E14_has_one_gauche_kink_and_is_otherwise_all_anti():
    """The fixture claim only. Whether the folded rule passes it is the
    behaviour test (E14), not this one."""
    from rdkit.Chem import rdMolTransforms
    m = _records(_BY_ID["E14"]["recipe"]["file"])[0]
    c = m.GetConformer()
    chain = [int(x) for x in m.GetProp("chain_atoms").split()]
    kink = [int(x) for x in m.GetProp("kink_atoms").split()]
    gauche = []
    for i in range(len(chain) - 3):
        d = rdMolTransforms.GetDihedralDeg(c, *chain[i:i + 4])
        if abs(abs(d) - 180) > 5:
            gauche.append((chain[i:i + 4], d))
    assert len(gauche) == 2
    (q1, d1), (q2, d2) = gauche
    assert q1 == kink[0:4] and q2 == kink[1:5]
    assert 55 < abs(d1) < 75 and 55 < abs(d2) < 75 and np.sign(d1) != np.sign(d2)


def _params(cost):
    return [pytest.param(c, id=c["id"],
                         marks=[pytest.mark.xfail(strict=True, reason="ensembles not implemented")]
                         if c["id"] in NOT_YET else [])
            for c in _CASES if c["cost"] == cost]


def _run(case):
    """Execute a case against the implementation (docs/ensemble-spec.md section 9
    steps 3-5). Written against the spec's names; adapt the call, never the
    expectation."""
    if not hasattr(em, "run_ensemble") or not hasattr(em, "compare_records"):
        pytest.fail("ensembles and comparison records are not implemented yet")
    raise NotImplementedError("wire the case recipe to run_ensemble / compare_records")


@pytest.mark.parametrize("case", _params("fast"))
def test_fast_case(case):
    _run(case)


@pytest.mark.slow
@pytest.mark.parametrize("case", _params("slow"))
def test_slow_case(case):
    _run(case)
