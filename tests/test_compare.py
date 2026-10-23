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
def _record(smiles, site_means, vmin=-20.0, key="k" * 40, version="0.3.0", valid=True):
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
        provenance=dict(espmap_version=version),
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
                            folded=bool(folded and folded[i])))
    e = dict(schema="espmap.ensemble/1", subject=name, analysis_key=recs[0]["analysis_key"],
             members=members, excluded=list(excluded),
             convergence=dict(converged=converged), validity=dict(unbiased=unbiased))
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
    b = _single(tmp_path, "b", _record("CCCCCCCO", {NN: 1.0, RS: 1.0}, version="0.2.1"))
    c = em.compare_records(a, b, CORR)
    assert set(c["flags"]) == {"CONFORMERS_NOT_SAMPLED", "VERSIONS_DIFFER", "NOT_CONVERGED",
                               "SELECTION_BIASED"}
    assert c["inputs"]["A"]["excluded"] == 1 and c["inputs"]["A"]["excluded_folded"] == [True]
    assert {c["inputs"]["A"]["espmap_version"], c["inputs"]["B"]["espmap_version"]} == {"0.3.0", "0.2.1"}


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
@pytest.mark.xfail(strict=True, reason="sd 0 on a single side lets resolved be true")
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
