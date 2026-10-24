"""espmap compare and espmap.compare/1 (docs/ensemble-spec.md section 6),
without an SCF: real fixture records where the chemistry matters (site
recombination), and synthetic records with known numbers for the rule."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

import esp_map as em

_FIX = Path(__file__).resolve().parents[1] / "fixtures" / "records"


def _fixture(name):
    return json.loads((_FIX / f"{name}.json").read_text(encoding="utf-8"))


# --- recombination: the site value IS the record's reaction-site value ----------
@pytest.mark.parametrize("name", ["Z3Z13-18OAc_extended", "Z3Z13-18OAc_folded", "benzene",
                                  "enolate", "nitrite"])
def test_a_recombined_site_equals_the_records_reaction_site(name):
    """Arrowfield's check: triangles are owned per atom, so the area-weighted
    mean over the union of the atoms' rows is the reaction site's own mean,
    hydrogens included, under both rules."""
    rec = _fixture(name)
    for rs in rec["sites"]["reaction_sites"]:
        for rule in em.PARTITION_RULES:
            want = (rs.get(rule) or {}).get("V_area_mean_kcal_mol")
            got = em.site_value(rec, rs["atoms"], rule)
            assert (got is None) == (want is None), (rs["label"], rule)
            if want is not None:
                assert got == pytest.approx(want, abs=1e-9), (rs["label"], rule)


def test_map_plus_h_names_the_same_atoms_as_the_records_c3_h_site():
    rec = _fixture("Z3Z13-18OAc_extended")
    mol = em.record_mol(rec)
    c3h = next(r for r in rec["sites"]["reaction_sites"] if r["label"] == "C3-H")
    atoms, why = em.resolve_site(mol, ["3+H"], "map")
    assert why is None and atoms == sorted(c3h["atoms"])


def test_record_mol_is_the_approved_molecule_in_record_order():
    for name in ("Z3Z13-18OAc_extended", "enolate", "nitrite", "benzene"):
        rec = _fixture(name)
        mol = em.record_mol(rec)
        assert [a.GetSymbol() for a in mol.GetAtoms()] == rec["geometry"]["symbols"]
        appr = Chem.MolFromSmiles(rec["input"]["approved_smiles"])
        for a in appr.GetAtoms():
            a.SetAtomMapNum(0)
        m = Chem.Mol(mol)
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        # connectivity, charges and bond orders; stereo is not carried (SMARTS
        # sites are matched without it)
        got = Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=False)
        assert got == Chem.MolToSmiles(appr, isomericSmiles=False), name


# --- correspondence -----------------------------------------------------------
def _mol(smiles):
    return Chem.AddHs(Chem.MolFromSmiles(smiles))


def _codes(fn, *a):
    with pytest.raises(em.CompareRefused) as e:
        fn(*a)
    return {c for _, c, _, _ in e.value.findings}


def test_sites_by_index_are_refused():                                     # C8
    corr = {"by": "index", "pairs": [{"name": "O", "A": [2], "B": [3]}]}
    assert _codes(em.check_correspondence, _mol("CCO"), _mol("CCCO"), corr) == {"COMPARE_BY_INDEX"}


def test_a_map_pair_of_two_elements_is_refused():                         # C9
    corr = {"by": "map", "pairs": [{"name": "X", "A": [2], "B": [2]}]}
    assert _codes(em.check_correspondence, _mol("[CH3:1][OH:2]"), _mol("[CH3:1][NH2:2]"),
                  corr) == {"COMPARE_MAP_ELEMENT"}


def test_smarts_matching_two_different_atoms_is_refused():                # C3
    corr = {"by": "smarts", "pairs": [{"name": "OH", "A": "[OX2H]", "B": "[OX2H]"}]}
    assert _codes(em.check_correspondence, _mol("CCO"), _mol("CC(O)CCO"), corr) == {
        "COMPARE_SITE_UNRESOLVED"}


def test_symmetry_equivalent_matches_are_one_site():                      # C10
    corr = {"by": "smarts", "pairs": [{"name": "O", "A": "[OX1]", "B": "[OX1]"}]}
    (p,) = em.check_correspondence(_mol("CC(=O)[O-]"), _mol("CCC(=O)[O-]"), corr)
    assert len(p["atoms_A"]) == 2 and len(p["atoms_B"]) == 2


def test_include_h_adds_the_sites_hydrogens():
    corr = {"by": "smarts", "pairs": [{"name": "OH", "A": "[OX2H]", "B": "[OX2H]",
                                       "include_h": True}]}
    (p,) = em.check_correspondence(_mol("CCO"), _mol("CCCO"), corr)
    assert len(p["atoms_A"]) == 2 and len(p["atoms_B"]) == 2


# --- synthetic records with known numbers ---------------------------------------
def _record(smiles, site_means, vmin=-20.0, key="k" * 40, version="0.3.0", valid=True, prov=None):
    """A record of smiles with every atom given area 1 and, under each rule,
    the mean site_means[rule] (same for every atom)."""
    m = _mol(smiles)
    return dict(
        schema="espmap.record/1", analysis_key=key,
        key_fields=dict(analysis=dict(xc="b3lyp", basis="def2svp" if key == "k" * 40 else "sto3g")),
        validity=dict(input=True, surface=True, electronic=True, sites=valid),
        geometry=dict(symbols=[a.GetSymbol() for a in m.GetAtoms()], bonds=em.record_bonds(
            Chem.MolFromSmiles(smiles)), sha1="0" * 40),
        input=dict(approved_smiles=smiles, computed_smiles=smiles),
        provenance=dict(espmap_version=version, **(prov or {})),
        results=dict(V_min_kcal_mol=vmin, V_max_kcal_mol=-vmin),
        sites=dict(atoms=[dict(index=a.GetIdx(), atom_map=a.GetAtomMapNum() or None,
                               element=a.GetSymbol(),
                               **{r: dict(area_A2=1.0, V_area_mean_kcal_mol=site_means[r])
                                  for r in em.PARTITION_RULES}) for a in m.GetAtoms()],
                   reaction_sites=[]))


def _ensemble(tmp, name, recs, folded=None, converged=True, unbiased=True, excluded=()):
    d = tmp / name
    (d / "conformers").mkdir(parents=True)
    members = []
    for i, r in enumerate(recs):
        f = d / "conformers" / f"conf_{i:02d}.json"
        f.write_text(json.dumps(r))
        members.append(dict(conformer=f"conf_{i:02d}", record=f"conformers/conf_{i:02d}.json",
                            record_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
                            geometry_sha1=r["geometry"]["sha1"],
                            folded=bool(folded and folded[i])))
    e = dict(schema="espmap.ensemble/1", subject=name, analysis_key=recs[0]["analysis_key"],
             members=members, excluded=list(excluded),
             convergence=dict(converged=converged),
             validity=dict(unbiased=unbiased, convergence=converged))
    (d / "ensemble.json").write_text(json.dumps(e))
    return d / "ensemble.json"


def _single(tmp, name, rec):
    p = tmp / f"{name}.json"
    p.write_text(json.dumps(rec))
    return p


NN, RS = em.PARTITION_RULES
CORR = {"by": "smarts", "pairs": [{"name": "OH", "A": "[OX2H]", "B": "[OX2H]", "include_h": True}]}


def test_a_record_compared_with_itself_is_unresolved(tmp_path):                 # C1
    recs = [_record("CCCCCCO", {NN: v, RS: v + 1}) for v in (1.0, 2.0, 3.0)]
    a = _ensemble(tmp_path, "a", recs)
    c = em.compare_records(a, a, CORR)
    s = c["sites"][0]["by_rule"][NN]
    assert abs(s["delta_mean"]) < 1e-9 and s["resolved"] is False
    assert c["flags"] == []


def test_the_rule_resolves_a_clear_difference_and_names_itself(tmp_path):
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: v, RS: v}) for v in (1.0, 1.2, 1.4)])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}) for v in (5.0, 5.2, 5.4)])
    c = em.compare_records(a, b, CORR)
    s = c["sites"][0]["by_rule"][NN]
    sd = np.std([1.0, 1.2, 1.4], ddof=1)
    assert s["delta_mean"] == pytest.approx(4.0) and s["U"] == pytest.approx(2 * np.hypot(sd, sd))
    assert s["sign_consistency"] == 1.0 and s["resolved"] is True
    assert s["resolved_under"] == c["rule"]["id"] and c["rule"]["status"] == "provisional"
    assert c["uncertainty_not_included"] == ["method", "basis", "solvent", "geometry_source"]


def test_a_delta_inside_the_partition_is_unresolved(tmp_path):                  # C4
    a = _single(tmp_path, "a", _record("CC(C)CO", {NN: 1.0, RS: 1.0}))
    b = _single(tmp_path, "b", _record("CCCCO", {NN: 2.0, RS: 6.0}))     # partition 4
    c = em.compare_records(a, b, CORR)
    s = c["sites"][0]
    assert s["partition"] == pytest.approx(4.0)
    assert s["by_rule"][NN]["delta_mean"] == pytest.approx(1.0) and s["by_rule"][NN]["resolved"] is False
    assert "CONFORMERS_NOT_SAMPLED" in c["flags"]
    assert s["by_rule"][NN]["conformational"]["A"]["measured"] == "not_measured"


def test_low_sign_consistency_is_unresolved(tmp_path):                          # C5
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: v, RS: v}) for v in (0.0, 0.1, 9.0)])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}) for v in (8.0, 8.1, 1.0)])
    s = em.compare_records(a, b, CORR)["sites"][0]["by_rule"][NN]
    assert s["sign_consistency"] < 0.9 and s["resolved"] is False


def test_different_levels_of_theory_are_refused(tmp_path):                      # C2
    a = _single(tmp_path, "a", _record("CCO", {NN: 1.0, RS: 1.0}))
    b = _single(tmp_path, "b", _record("CCCO", {NN: 1.0, RS: 1.0}, key="s" * 40))
    assert _codes(em.compare_records, a, b, CORR) == {"COMPARE_ANALYSIS_DIFFERS"}


def test_flags_versions_convergence_and_bias(tmp_path):                         # C7, C11
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: 1.0, RS: 1.0})] * 3, converged=False,
                  unbiased=False, excluded=[dict(conformer="conf_03", folded=True)])
    # R7-7 (the user's triage of review 7): v0.2.1 vs v0.3.0 is now refused
    # (K1); two releases under one contract still compare, with the flag
    b = _single(tmp_path, "b", _record("CCCCCCCO", {NN: 1.0, RS: 1.0}, version="0.3.1",
                                       prov=dict(code="v0.3.1", measurement_contract="espmap.numeric/1")))
    c = em.compare_records(a, b, CORR)
    assert set(c["flags"]) == {"CONFORMERS_NOT_SAMPLED", "VERSIONS_DIFFER", "NOT_CONVERGED",
                               "SELECTION_BIASED"}
    assert c["inputs"]["A"]["excluded"] == 1 and c["inputs"]["A"]["excluded_folded"] == [True]
    assert {c["inputs"]["A"]["espmap_version"], c["inputs"]["B"]["espmap_version"]} == {"0.3.0", "0.3.1"}


def test_per_member_rows_and_descriptive_partitions_are_in_the_record(tmp_path):
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: v, RS: v + d}) for v, d in
                                  ((1.0, 0.5), (2.0, 2.0))], folded=[False, True])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}) for v in (3.0, 4.0)])
    s = em.compare_records(a, b, CORR)["sites"][0]
    assert [r["folded"] for r in s["members"]["A"]] == [False, True]
    assert [r[NN] for r in s["members"]["A"]] == [1.0, 2.0]
    assert s["partition_descriptive"]["max_member"] == pytest.approx(2.0)
    assert s["partition_descriptive"]["between_side_means"] == pytest.approx(1.25)


def test_global_extrema_are_compared_like_a_site(tmp_path):
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: 1, RS: 1}, vmin=v) for v in (-20, -21)])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: 1, RS: 1}, vmin=v) for v in (-30, -31)])
    g = em.compare_records(a, b, CORR)["global_extrema"]["V_min_kcal_mol"]
    assert g["delta_mean"] == pytest.approx(-10.0) and g["partition"] == "not applicable"
    assert [r["value"] for r in g["members"]["B"]] == [-30, -31]


def test_a_changed_member_record_is_refused(tmp_path):
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: 1, RS: 1})] * 2)
    f = a.parent / "conformers" / "conf_01.json"
    f.write_text(f.read_text().replace('"0.3.0"', '"9.9.9"'))
    assert _codes(em.compare_records, a, a, CORR) == {"COMPARE_MEMBER_CHANGED"}


def test_an_invalid_site_table_is_refused(tmp_path):
    a = _single(tmp_path, "a", _record("CCO", {NN: 1, RS: 1}, valid=False))
    assert _codes(em.compare_records, a, a, CORR) == {"COMPARE_INVALID"}


def test_the_command_line(tmp_path):
    import subprocess
    import sys
    a = _single(tmp_path, "a", _record("CCO", {NN: 1, RS: 1}))
    b = _single(tmp_path, "b", _record("CCCO", {NN: 3, RS: 3}))
    corr = tmp_path / "corr.json"
    corr.write_text(json.dumps(CORR))
    r = subprocess.run([sys.executable, em.__file__, "compare", str(a), str(b), "--sites", str(corr),
                        "--out", str(tmp_path / "c.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    c = json.loads((tmp_path / "c.json").read_text())
    assert c["schema"] == "espmap.compare/1" and c["sites"][0]["name"] == "OH"


# --- the user's decision (2026-10-24): a single-geometry side is never resolved -------
@pytest.mark.parametrize("sides", ["single_single", "single_ensemble", "ensemble_single"])
def test_a_single_geometry_side_is_never_resolved(tmp_path, sides):
    """Conformational uncertainty is not measured on a single geometry, so the
    comparison may report delta and partition but never cross into
    'resolved', however large delta is."""
    a_rec = [_record("CCCCCCO", {NN: v, RS: v}) for v in (1.0, 1.1, 1.2)]
    b_rec = [_record("CCCCCCCO", {NN: v, RS: v}) for v in (30.0, 30.1, 30.2)]
    a = (_single(tmp_path, "a", a_rec[0]) if sides.startswith("single")
         else _ensemble(tmp_path, "a", a_rec))
    b = (_single(tmp_path, "b", b_rec[0]) if sides.endswith("single")
         else _ensemble(tmp_path, "b", b_rec))
    c = em.compare_records(a, b, CORR)
    for rule in em.PARTITION_RULES:
        s = c["sites"][0]["by_rule"][rule]
        assert abs(s["delta_mean"]) > 25 and s["resolved"] is False
    for q in ("V_min_kcal_mol", "V_max_kcal_mol"):
        assert c["global_extrema"][q]["resolved"] is False
    assert "CONFORMERS_NOT_SAMPLED" in c["flags"]
    single = [s for s in "AB" if c["inputs"][s]["kind"] == "single"]
    for s in single:
        assert c["sites"][0]["by_rule"][NN]["conformational"][s]["measured"] == "not_measured"



def test_a_one_member_ensemble_is_a_single_geometry(tmp_path):
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: 1.0, RS: 1.0})])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}) for v in (30.0, 30.1, 30.2)])
    c = em.compare_records(a, b, CORR)
    s = c["sites"][0]["by_rule"][NN]
    assert s["resolved"] is False and s["conformational"]["A"]["measured"] == "not_measured"
    assert "CONFORMERS_NOT_SAMPLED" in c["flags"] and "side A" in s["not_resolvable"]


# --- review 7: counterexamples before the fixes ------------------------------------
def _no_resolution(s):
    return all(s["by_rule"][rl]["resolved"] is not True for rl in em.PARTITION_RULES)


def test_r7_e8_a_site_missing_on_one_member_is_not_evaluated_on_n_minus_1(tmp_path):
    """Review 7 R7-E8: the folded member buries the site (None under both
    rules). The comparison may not quietly use the other two."""
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: v, RS: v}) for v in (1.0, 1.2)]
                  + [_record("CCCCCCO", {NN: None, RS: None})], folded=[False, False, True])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}) for v in (5.0, 5.2, 5.4)])
    try:
        c = em.compare_records(a, b, CORR)
    except em.CompareRefused as e:
        assert {code for _, code, _, _ in e.value.findings} == {"COMPARE_SITE_INCOMPLETE"}
        return
    s = c["sites"][0]
    assert _no_resolution(s)
    assert s.get("evaluable") is False and "conf_02" in json.dumps(s.get("not_evaluable"))


def test_r7_e9_an_unconverged_ensemble_never_resolves(tmp_path):
    """Review 7 R7-E9: a huge delta, but side A's convergence is false."""
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: v, RS: v}, vmin=-v)
                                  for v in (1.0, 1.2, 1.4)], converged=False)
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}, vmin=-v)
                                  for v in (50.0, 50.2, 50.4)])
    c = em.compare_records(a, b, CORR)
    s = c["sites"][0]
    assert _no_resolution(s)
    assert "convergence" in s["by_rule"][NN].get("not_resolvable", "")
    assert all(g["resolved"] is not True for g in c["global_extrema"].values())


def test_r7_e10_convergence_not_judgeable_never_resolves(tmp_path):
    """Review 7 R7-E10: K = 2 on each side, so convergence is null."""
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: v, RS: v}) for v in (1.0, 1.2)],
                  converged=None)
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}) for v in (50.0, 50.2)],
                  converged=None)
    s = em.compare_records(a, b, CORR)["sites"][0]
    assert _no_resolution(s)
    assert "convergence" in s["by_rule"][NN].get("not_resolvable", "")


def test_r7_e12_a_parent_claiming_another_analysis_key_is_refused(tmp_path):
    """Review 7 R7-E12: both parents claim key x, their members were computed
    under key k. The parent's claim is what the key check reads."""
    paths = []
    for name, smi in (("a", "CCCCCCO"), ("b", "CCCCCCCO")):
        p = _ensemble(tmp_path, name, [_record(smi, {NN: v, RS: v}) for v in (1.0, 2.0, 3.0)])
        e = json.loads(p.read_text())
        e["analysis_key"] = "x" * 40
        p.write_text(json.dumps(e))
        paths.append(p)
    assert _codes(em.compare_records, *paths, CORR)


@pytest.mark.parametrize("doc", [
    "record2",                                    # a future record, otherwise readable
    "foreign",                                    # not an espmap file at all
])
def test_r7_e13_an_unknown_schema_is_refused_safely(tmp_path, doc):
    """Review 7 R7-E13: refuse with COMPARE_SCHEMA_UNSUPPORTED; never
    interpret the file, never crash with a KeyError."""
    good = _single(tmp_path, "good", _record("CCO", {NN: 1.0, RS: 1.0}))
    rec = (dict(_record("CCCO", {NN: 1.0, RS: 1.0}), schema="espmap.record/2")
           if doc == "record2" else {"hello": 1})
    bad = _single(tmp_path, "bad", rec)
    for args in ((good, bad), (bad, good)):
        assert _codes(em.compare_records, *args, CORR) == {"COMPARE_SCHEMA_UNSUPPORTED"}


# --- review 7 R7-10: the record validates against docs/compare-schema.json -------
def _validator():
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((Path(__file__).resolve().parents[1] / "docs" / "compare-schema.json")
                        .read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(schema)


def test_r7_10_comparison_records_validate_against_the_schema(tmp_path):
    v = _validator()
    ens = lambda n, smi, vals, **kw: _ensemble(tmp_path, n, [_record(smi, {NN: x, RS: x}) for x in vals], **kw)
    cases = [
        (ens("a", "CCCCCCO", (1.0, 1.2, 1.4)), ens("b", "CCCCCCCO", (5.0, 5.2, 5.4))),     # resolved
        (_single(tmp_path, "s1", _record("CCO", {NN: 1.0, RS: 2.0})),
         _single(tmp_path, "s2", _record("CCCO", {NN: 3.0, RS: 3.5}))),                   # singles
        (ens("c", "CCCCCCO", (1.0, 1.2), converged=None),
         ens("d", "CCCCCCCO", (9.0, 9.2), converged=None)),                               # null
        (_ensemble(tmp_path, "e", [_record("CCCCCCO", {NN: 1.0, RS: 1.0}),
                                   _record("CCCCCCO", {NN: None, RS: None})]),
         ens("f", "CCCCCCCO", (2.0, 2.1, 2.2))),                                          # buried
    ]
    seen = set()
    for a, b in cases:
        c = json.loads(json.dumps(em.compare_records(a, b, CORR)))       # as written to disk
        errors = sorted(v.iter_errors(c), key=lambda e: list(e.path))
        assert not errors, [f"{list(e.path)}: {e.message}" for e in errors]
        seen.add(c["sites"][0]["by_rule"][NN]["resolved"])
    assert seen == {True, False}


def test_r7_10_the_schema_rejects_a_resolved_site_that_was_not_evaluable(tmp_path):
    v = _validator()
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: x, RS: x}) for x in (1.0, 1.2, 1.4)])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: x, RS: x}) for x in (5.0, 5.2, 5.4)])
    c = json.loads(json.dumps(em.compare_records(a, b, CORR)))
    assert not list(v.iter_errors(c))
    c["sites"][0]["evaluable"] = False                     # a tampered record
    assert list(v.iter_errors(c))


# --- review 7 R7-7: a measurement contract, not a version warning ----------------
# Same analysis_key, different espmap: only a declared contract may make the
# numbers comparable. Counterexamples first; strict xfails until the fix.
R77 = pytest.mark.xfail(strict=True, reason="R7-7: a version difference only warns")
CONTRACT_CODES = {"COMPARE_CONTRACT_UNKNOWN", "COMPARE_CONTRACT_DIFFERS"}


def _sides(tmp, a, b):
    """Two single records, A hexanol and B heptanol, provenance as given."""
    return (_single(tmp, "a", _record("CCCCCCO", {NN: 1.0, RS: 1.0}, **a)),
            _single(tmp, "b", _record("CCCCCCCO", {NN: 2.0, RS: 2.0}, **b)))


def _refused_by_contract(tmp, a, b):
    try:
        em.compare_records(*_sides(tmp, a, b), CORR)
    except em.CompareRefused as e:
        return {c for _, c, _, _ in e.findings} <= CONTRACT_CODES
    return False


def test_r7_7_k1_a_v021_record_does_not_compare_with_a_v030_one(tmp_path):
    """Even under one analysis_key (a key collision, or an edited record):
    v0.2.x declares no contract."""
    assert _refused_by_contract(tmp_path, dict(version="0.2.1", prov=dict(code="v0.2.1")),
                                dict(version="0.3.0", prov=dict(code="v0.3.0")))


def test_r7_7_k2_a_dev_build_does_not_compare_with_a_release(tmp_path):
    assert _refused_by_contract(tmp_path, dict(version="0.3.0.dev0", prov=dict(code="v0.2.1-86-g976d343")),
                                dict(version="0.3.0", prov=dict(code="v0.3.0")))


def test_r7_7_k4_one_dev_version_two_commits_do_not_compare(tmp_path):
    """The same version string from two commits is two programs."""
    assert _refused_by_contract(tmp_path, dict(version="0.3.1.dev0", prov=dict(code="v0.3.0-3-gaaaaaaa")),
                                dict(version="0.3.1.dev0", prov=dict(code="v0.3.0-5-gbbbbbbb")))


def test_r7_7_k6_two_contracts_do_not_compare(tmp_path):
    assert _refused_by_contract(
        tmp_path, dict(version="0.3.1", prov=dict(code="v0.3.1", measurement_contract="espmap.numeric/1")),
        dict(version="0.4.0", prov=dict(code="v0.4.0", measurement_contract="espmap.numeric/2")))


def test_r7_7_k7_an_unknown_contract_is_refused(tmp_path):
    assert _refused_by_contract(
        tmp_path, dict(version="0.3.1", prov=dict(code="v0.3.1", measurement_contract="espmap.numeric/9")),
        dict(version="0.3.1", prov=dict(code="v0.3.1-2-gccccccc", measurement_contract="espmap.numeric/9")))


def test_r7_7_k8_a_contract_that_contradicts_its_release_is_refused(tmp_path):
    """v0.3.0 is numeric/1 by declaration; a v0.3.0 record claiming another
    contract was edited or forged."""
    assert _refused_by_contract(
        tmp_path, dict(version="0.3.0", prov=dict(code="v0.3.0", measurement_contract="espmap.numeric/2")),
        dict(version="0.3.0", prov=dict(code="v0.3.0")))


def test_r7_7_k9_an_ensemble_of_two_programs_is_refused(tmp_path):
    """Members computed by two espmap builds: the ensemble has no one contract."""
    a = _ensemble(tmp_path, "a", [_record("CCCCCCO", {NN: v, RS: v}, version=ver, prov=dict(code=code))
                                  for v, ver, code in ((1.0, "0.3.0", "v0.3.0"), (1.1, "0.3.0", "v0.3.0"),
                                                       (1.2, "0.3.0.dev0", "v0.2.1-86-g976d343"))])
    b = _ensemble(tmp_path, "b", [_record("CCCCCCCO", {NN: v, RS: v}) for v in (5.0, 5.1, 5.2)])
    with pytest.raises(em.CompareRefused):
        em.compare_records(a, b, CORR)


def test_r7_7_k3_one_build_compares_with_itself(tmp_path):
    """A dev build is comparable with records from the very same commit."""
    a, b = _sides(tmp_path, dict(version="0.3.1.dev0", prov=dict(code="v0.3.0-3-gaaaaaaa")),
                  dict(version="0.3.1.dev0", prov=dict(code="v0.3.0-3-gaaaaaaa")))
    assert "VERSIONS_DIFFER" not in em.compare_records(a, b, CORR)["flags"]


def test_r7_7_k5_two_releases_under_one_contract_compare_with_a_flag(tmp_path):
    """v0.3.0 (numeric/1 by declaration) and a later release that writes the
    same contract: comparable, and the record says the versions differ."""
    a, b = _sides(tmp_path, dict(version="0.3.0", prov=dict(code="v0.3.0")),
                  dict(version="0.3.1", prov=dict(code="v0.3.1", measurement_contract="espmap.numeric/1")))
    assert "VERSIONS_DIFFER" in em.compare_records(a, b, CORR)["flags"]
