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
from rdkit.Chem import rdMolDescriptors, rdMolTransforms

import esp_map as em

_DIR = Path(__file__).resolve().parent / "data" / "ensemble"
_CASES = json.loads((_DIR / "manifest.json").read_text(encoding="utf-8"))["cases"]
_BY_ID = {c["id"]: c for c in _CASES}

# implementation PRs remove entries. Step 3 (selection and the supplied
# identity checks) implements the refusals before compute and E14's rule, and
# with review 6's fixes E18, E20-E22, E24 and E25 (E23 refused already).
NOT_YET = {c["id"] for c in _CASES} - {"E8", "E9", "E12", "E13", "E14", "E15", "E16",
                                       "E17", "E18", "E20", "E21", "E22", "E23",
                                       "E24", "E25", "E26", "E27", "E28"}


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
    assert len(ids) == len(set(ids)) == 40
    assert {i for i in ids if i[0] == "E"} == {f"E{i}" for i in range(1, 29)}
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


def _signs(name, quad):
    return [1 if rdMolTransforms.GetDihedralDeg(m.GetConformer(), *quad) > 0 else -1
            for m in _records(name)]


def test_E19_E20_E21_fix_one_centre_and_vary_what_rdkit_misses():
    for cid, quad in (("E20", (3, 4, 5, 6)), ("E21", (1, 3, 5, 6))):
        name = _BY_ID[cid]["recipe"]["file"]
        assert len({_key(m) for m in _records(name)}) == 1, cid     # RDKit: one stereoisomer
        assert set(_signs(name, quad)) == {1, -1}, cid               # the torsion: two
    assert _signs(_BY_ID["E19"]["recipe"]["file"], (1, 3, 5, 6)) == [1]


def test_E22_E23_are_one_skeleton_with_three_electron_counts():
    for cid in ("E22", "E23"):
        ms = _records(_BY_ID[cid]["recipe"]["file"])
        assert sorted(Chem.GetFormalCharge(m) for m in ms) == [-1, 0, 1]
        assert len({tuple(a.GetSymbol() for a in m.GetAtoms()) for m in ms}) == 1
    first = _records(_BY_ID["E23"]["recipe"]["file"])[0]
    assert sum(a.GetNumRadicalElectrons() for a in first.GetAtoms()) == 1


def test_E24_maps_differ_on_the_same_atoms():
    ms = _records(_BY_ID["E24"]["recipe"]["file"])
    maps = [tuple(a.GetAtomMapNum() for a in m.GetAtoms() if a.GetAtomicNum() > 1) for m in ms]
    assert maps == [(1, 2, 3), (3, 2, 1)]


def test_E25_record_2_is_unreadable_and_the_others_are_not():
    assert [m is None for m in _records(_BY_ID["E25"]["recipe"]["file"])] == [False, True, False]


def test_E26_swaps_labels_only():
    a, b = _records(_BY_ID["E26"]["recipe"]["file"])
    assert _key(a) == _key(b)                          # one molecule, one configuration
    pa, pb = a.GetConformer().GetPositions(), b.GetConformer().GetPositions()
    assert not np.allclose(pa, pb)                     # the labels moved
    assert sorted(map(tuple, np.round(pa, 4))) == sorted(map(tuple, np.round(pb, 4)))  # the points did not


def test_E27_E28_are_two_isotopologues_in_both_orders():
    for cid in ("E27", "E28"):
        ms = _records(_BY_ID[cid]["recipe"]["file"])
        assert sorted(sum(a.GetIsotope() == 2 for a in m.GetAtoms()) for m in ms) == [0, 3]
        assert len({_key(m, stereo=False) for m in ms}) == 2         # SMILES shows the labels


def _params(cost):
    return [pytest.param(c, id=c["id"],
                         marks=[pytest.mark.xfail(strict=True, reason="ensembles not implemented")]
                         if c["id"] in NOT_YET else [])
            for c in _CASES if c["cost"] == cost]


def _mol_block_from_spec(spec, recipe):
    """The file's first record carrying spec's coordinates (same atom order)."""
    m = Chem.Mol(_records(recipe["file"])[0])
    conf = m.GetConformer()
    for i, row in enumerate(spec.split(";")):
        conf.SetAtomPosition(i, [float(x) for x in row.split()[1:]])
    return Chem.MolToMolBlock(m)


def _refused_before_compute(expect):
    return any(e.get("outcome") == "refuse" and e.get("before_compute")
               for e in (expect, expect.get("or", {})))


def _plan_kwargs(recipe):
    kw = dict(ensemble=recipe["ensemble"], weight=recipe.get("weight"),
              optimize=recipe.get("optimize", False))
    if "file" in recipe:
        kw["path"] = _DIR / recipe["file"]
    else:
        kw["smiles"] = recipe["smiles"]
    return kw


# the protocol the slow cases compute under: the C18 measurements' (gas,
# def2-SVP + RI), cheap enough for the shared machine one member at a time
PROTOCOL = ["--basis", "def2-svp", "--density-fit"]


def _get(rec, dotted):
    for k in dotted.split("."):
        rec = rec[k]
    return rec


def _check(rec, checks):
    """Each manifest check, by the name the manifest uses."""
    ms = rec["members"]
    for name, want in checks.items():
        if name == "max_members" or name == "k_obtained_max":
            assert rec["k_obtained"] <= want, name
        elif name == "no_duplicates":
            assert len({m["geometry_sha1"] for m in ms}) == len(ms) == rec["k_obtained"]
        elif name == "folded_count_min":
            assert sum(m["folded"] for m in ms) >= want
        elif name == "stats_over":
            assert set(want) <= set(rec["results"]["V_min_kcal_mol"])
        elif name == "excluded_min":
            assert len(rec["excluded"]) >= want
        elif name == "member_atom_order":                  # "verified_mapping"
            assert all(m["validity"]["input"] for m in ms)
        elif name == "configurations":
            assert len({json.dumps(m["configuration"], sort_keys=True) for m in ms}) == want
        elif name == "kish_n_reported":
            assert (rec["kish_n"] is not None) == want
        elif name == "folded_by_genus_alone":
            assert sum(m["folded"] and not m["contacts"] for m in ms) == want
        else:                                              # a record field, dotted
            assert _get(rec, name) == want, name


def _run(case, tmp_path=None):
    """Execute a case against the implementation (docs/ensemble-spec.md section 9
    steps 3-5). Written against the spec's names; adapt the call, never the
    expectation."""
    r, e = case["recipe"], case["expect"]
    if r.get("classify_only"):
        got = em.shape_descriptors(_records(r["file"])[0])
        assert {k: got[k] for k in e["checks"]} == e["checks"]
        return
    if r.get("plan_only"):
        # an expectation decided before compute (E26: relabelling)
        plan = em.plan_ensemble(**_plan_kwargs(r))
        for name, want in e["checks"].items():
            if name == "relabelled_records":
                assert [m["record"] for m in plan["members"] if m["relabelled"]] == want
            elif name == "pro_labels_agree":
                mols = [Chem.MolFromMolBlock(_mol_block_from_spec(m["spec"], r)) for m in plan["members"]]
                labels = {json.dumps(sorted(em.prochiral_labels(x).items(), key=str)) for x in mols}
                assert (len(labels) == 1) == want
        return
    if "ensemble" in r and _refused_before_compute(e):
        # plan_ensemble is everything before the first SCF. A case whose
        # expectation also allows "accept" (E9's and E19's "or") passes by
        # refusing, or goes on to the runner.
        if e["outcome"] == "refuse":
            with pytest.raises(em.EnsembleRefused):
                em.plan_ensemble(**_plan_kwargs(r))
            return
        try:
            em.plan_ensemble(**_plan_kwargs(r))
        except em.EnsembleRefused:
            return
    if "A" in r:
        if not hasattr(em, "compare_records"):
            pytest.fail("comparison records are not implemented yet")
        raise NotImplementedError("wire the comparison recipe to compare_records (step 5)")
    if r.get("inject") == "surface_failure_on_folded":
        # E5: a folded member's surface fails integrity (a 0.5 bohr box clips
        # it, and espmap refuses a clipped surface), as the chi = 2 check once
        # did to exactly the folded C18 members (F9)
        plan = em.plan_ensemble(**_plan_kwargs(r))
        folded = {em.member_label(i) for i, m in enumerate(plan["members"]) if m["folded"]}
        assert folded, "E5 needs a folded member in the plan"

        def runner(argv, cwd, env, timeout):
            label = Path(argv[argv.index("--file") + 1]).stem
            extra = ["--margin", "0.5"] if label in folded else []
            return em._run_member_process(argv + extra, cwd, env, timeout)
        rec = em.run_ensemble(plan, tmp_path, PROTOCOL, runner=runner)
        _check(rec, e["checks"])
        return
    if r.get("inject"):
        raise NotImplementedError(f"injection {r['inject']!r} is not wired yet")
    if "members_from" in r:
        # E10: one member per entry, each under its own basis; must be refused
        plan = em.plan_ensemble(smiles=r["members_from"][0]["smiles"],
                                ensemble=len(r["members_from"]))
        bases = [m["basis"] for m in r["members_from"]]

        def runner(argv, cwd, env, timeout):
            i = int(Path(argv[argv.index("--file") + 1]).stem.split("_")[1])
            return em._run_member_process(argv + ["--basis", bases[i], "--density-fit"],
                                          cwd, env, timeout)
        with pytest.raises(em.EnsembleRefused):
            em.run_ensemble(plan, tmp_path, runner=runner)
        return
    plan = em.plan_ensemble(**_plan_kwargs(r))
    rec = em.run_ensemble(plan, tmp_path, PROTOCOL + (["--optimize"] if r.get("optimize") else []))
    assert e["outcome"] == "accept"
    _check(rec, e.get("checks", {}))


@pytest.mark.parametrize("case", _params("fast"))
def test_fast_case(case, tmp_path):
    _run(case, tmp_path)


@pytest.mark.slow
@pytest.mark.parametrize("case", _params("slow"))
def test_slow_case(case, tmp_path):
    _run(case, tmp_path)
