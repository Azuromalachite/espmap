"""Ensemble selection and the supplied-ensemble identity checks
(docs/ensemble-spec.md sections 1-3, step 3). Fast: embedding and force
fields only, no SCF. The hostile set's own cases are in test_ensemble_set.py;
these pin WHY each refusal happens, so a refusal for the wrong reason fails.
"""
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

import esp_map as em

_DIR = Path(__file__).resolve().parent / "data" / "ensemble"


def _refusal(**kw):
    with pytest.raises(em.EnsembleRefused) as e:
        em.plan_ensemble(**kw)
    return [(c, m) for _, c, m, _ in e.value.findings]


def _supplied(name):
    return _refusal(path=_DIR / name, ensemble="supplied")


def test_a_tautomer_mix_is_refused_on_each_tautomer_record():
    got = _supplied("E8_ligand_and_tautomer.sdf")
    assert {c for c, _ in got} == {"MEMBER_IDENTITY"}
    assert [m.split(" is ")[0] for _, m in got] == ["record 4", "record 5", "record 6"]


def test_a_reordered_record_is_refused_not_trusted_by_file_order():
    got = _supplied("E9_permuted_atom_order.sdf")
    assert got == [("MEMBER_ATOM_ORDER", got[0][1])] and got[0][1].startswith("record 4")


def test_one_inverted_centre_is_a_stereo_refusal_on_that_record_only():
    got = _supplied("E15_one_centre_inverted.sdf")
    assert [c for c, _ in got] == ["MEMBER_STEREO"] and got[0][1].startswith("record 6")


def test_constitutional_isomers_are_refused():
    assert {c for c, _ in _supplied("E17_constitutional_isomer_mixed.sdf")} == {"MEMBER_IDENTITY"}


def test_twenty_poses_of_one_ligand_are_one_supplied_ensemble():
    p = em.plan_ensemble(path=_DIR / "E7_poses_one_ligand.sdf", ensemble="supplied")
    assert p["source"] == "supplied" and len(p["members"]) == 20
    assert p["protocol"]["k_obtained"] == 20
    assert p["identity"]["smiles"] == "CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1"
    # the coordinates are the file's, record by record
    first = next(iter(Chem.SDMolSupplier(str(_DIR / "E7_poses_one_ligand.sdf"), removeHs=False)))
    assert p["members"][0]["spec"] == em._mol_to_spec(first)


def _write(path, mols):
    w = Chem.SDWriter(str(path))
    for m in mols:
        w.write(m)
    w.close()
    return path


def _poses(n):
    return list(Chem.SDMolSupplier(str(_DIR / "E7_poses_one_ligand.sdf"), removeHs=False))[:n]


def test_a_member_without_hydrogens_is_refused(tmp_path):
    ms = _poses(3)
    ms[1] = Chem.RemoveHs(ms[1])
    got = _refusal(path=_write(tmp_path / "noh.sdf", ms), ensemble="supplied")
    assert [c for c, _ in got] == ["MEMBER_MISSING_H"]


def test_a_drawn_member_is_refused(tmp_path):
    ms = _poses(3)
    flat = Chem.Mol(ms[2])
    Chem.rdDepictor.Compute2DCoords(flat)
    got = _refusal(path=_write(tmp_path / "drawn.sdf", ms[:2] + [flat]), ensemble="supplied")
    assert got == [("MEMBER_NOT_3D", "record 3 is not a 3D geometry")]


def test_a_single_record_is_not_a_supplied_ensemble(tmp_path):
    got = _refusal(path=_write(tmp_path / "one.sdf", _poses(1)), ensemble="supplied")
    assert [c for c, _ in got] == ["ENSEMBLE_ONE_RECORD"]


@pytest.mark.parametrize("kw, code", [
    (dict(smiles="CCO", ensemble="supplied"), "ENSEMBLE_SOURCE"),
    (dict(smiles="CCO", ensemble=0), "ENSEMBLE_K"),
    (dict(smiles="CCO", ensemble=True), "ENSEMBLE_K"),
    (dict(smiles="CCO", ensemble=101), "ENSEMBLE_K"),
    (dict(smiles="CCCCCCO", ensemble=6, weight="mmff"), "WEIGHT_MMFF"),
    (dict(smiles="CCCCCCO", ensemble=6, weight="dft"), "WEIGHT_NEEDS_OPTIMIZE"),
    (dict(smiles="CCCCCCO", ensemble=6, weight="boltzmann"), "WEIGHT_UNKNOWN"),
    (dict(smiles="CC(O)C(N)C(=O)O", ensemble=10), "STEREO_REQUIRED"),
    (dict(smiles="C/C=C/C(O)CC", ensemble=4), "STEREO_REQUIRED"),
])
def test_refused_before_embedding(kw, code, monkeypatch):
    def no_embedding(*a, **k):
        raise AssertionError("embedded before refusing")
    monkeypatch.setattr(em, "_embed_and_minimise", no_embedding)
    assert code in {c for c, _ in _refusal(**kw)}


def test_dft_weights_with_optimised_members_pass_the_preflight():
    assert em._weight_problems("dft", optimize=True) == []


def test_the_embedder_returns_every_converged_conformer():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCO"))
    best, ff, spread, n_emb, n_conv, conv = em._embed_and_minimise(mol, 12, 7)
    assert len(conv) == n_conv <= n_emb == mol.GetNumConformers() == 12
    assert best == min(conv, key=lambda ce: ce[1])[0]
    energies = [e for _, e in conv]
    assert spread == pytest.approx(max(energies) - min(energies))


@pytest.mark.parametrize("smiles, k, obtained", [
    ("c1ccccc1", 10, 1),       # E1: rigid
    ("C1CCCCC1", 10, 2),       # E2: the chair and a twist-boat, distinct at 0.35 A
    ("CCCC", 10, 2),           # E6: anti and gauche; gauche- is gauche+'s mirror
])
def test_selection_never_pads(smiles, k, obtained):
    p = em.plan_ensemble(smiles=smiles, ensemble=k)
    pr = p["protocol"]
    assert (pr["k_requested"], pr["k_obtained"], len(p["members"])) == (k, obtained, obtained)
    assert pr["n_unique"] == obtained


def test_mirror_images_merge_only_for_an_achiral_molecule():
    assert em._is_achiral(Chem.MolFromSmiles("CCCC"))
    assert em._is_achiral(Chem.MolFromSmiles("O[C@@H](C(=O)O)[C@H](O)C(=O)O"))   # meso
    assert not em._is_achiral(Chem.MolFromSmiles("C[C@H](O)CC"))
    assert em._is_achiral(Chem.MolFromSmiles("[CH3:1]C(=O)[OH:2]"))   # maps are labels


def test_mirror_images_are_one_conformer_of_an_achiral_molecule(monkeypatch):
    """n-pentane has four conformers up to reflection (aa, ag, gg, g+g-) and
    six without it. A tight threshold isolates the mirror rule from the
    RMSD one: at 0.5 A butane's g+ and g- are already within the threshold
    (4 heavy atoms), and pentane loses a real conformer (section 8 measures
    the threshold)."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCC"))
    *_, conv = em._embed_and_minimise(mol, 100, 0xC0FFEE)
    assert em.select_ensemble(mol, conv, 50, dedup_rmsd=0.1)["n_unique"] == 4
    monkeypatch.setattr(em, "_is_achiral", lambda m: False)
    assert em.select_ensemble(mol, conv, 50, dedup_rmsd=0.1)["n_unique"] == 6


def test_selection_is_deterministic_and_starts_at_the_lowest_energy():
    a = em.plan_ensemble(smiles="CCCCCCO", ensemble=6)
    b = em.plan_ensemble(smiles="CCCCCCO", ensemble=6)
    assert [m["spec"] for m in a["members"]] == [m["spec"] for m in b["members"]]
    e = [m["mmff_energy"] for m in a["members"]]
    assert e[0] == min(e)
    assert len({m["embed_conf_id"] for m in a["members"]}) == len(e) == 6


def test_selection_spans_the_rg_range():
    """The extremes of Rg among the distinct conformers are both chosen."""
    p = em.plan_ensemble(smiles="CCCCCCO", ensemble=6)
    rg = sorted(m["rg"] for m in p["members"])
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCCCO"))
    *_, conv = em._embed_and_minimise(mol, em.N_EMBED, 0xC0FFEE)
    every = em.select_ensemble(mol, conv, em.N_EMBED)["members"]
    all_rg = [em.shape_descriptors(mol, c)["rg"] for c, _ in every]
    assert rg[0] == pytest.approx(min(all_rg)) and rg[-1] == pytest.approx(max(all_rg))


def test_contacts_follow_the_stated_rule():
    p = em.plan_ensemble(smiles="CCCCCCO", ensemble=6)
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCCCO"))
    topo = Chem.GetDistanceMatrix(mol)
    assert p["folded_count"] == sum(bool(m["contacts"]) for m in p["members"])
    for m in p["members"]:
        for i, j, d in m["contacts"]:
            assert topo[i, j] >= em.CONTACT_MIN_BONDS
            limit = sum(em.BONDI[mol.GetAtomWithIdx(x).GetSymbol()] for x in (i, j))
            assert d < limit + em.CONTACT_SLACK


def test_a_single_gauche_kink_is_not_folded():
    m = next(iter(Chem.SDMolSupplier(str(_DIR / "E14_c18_single_gauche_kink.sdf"),
                                     removeHs=False)))
    s = em.shape_descriptors(m)
    assert s["folded"] is False and s["contacts"] == []
    assert 0 < s["asphericity"] <= 1 and s["extent"] > 2 * s["rg"]


def test_shape_descriptors_of_a_rod_and_a_ball():
    rod = Chem.AddHs(Chem.MolFromSmiles("C#CC#CC#C"))
    em._embed_and_minimise(rod, 1, 1)
    assert em.shape_descriptors(rod)["asphericity"] > 0.95
    ball = Chem.AddHs(Chem.MolFromSmiles("C12C3C4C1C5C2C3C45"))      # cubane
    em._embed_and_minimise(ball, 1, 1)
    assert em.shape_descriptors(ball)["asphericity"] < 0.02


# --- review 6: committed as strict xfails before the fixes (c92a448) ----------
C18 = "CC(=O)OCC/C=C\CCCCCCCC/C=C\CCCC"
BENCH_C18_IDS = [4, 28, 34, 40, 53, 54, 57, 67, 72, 82]      # bench/conformers.jsonl conf_id


def test_selection_reproduces_the_measured_c18_members():
    """The conformer rows (5.1 / 2.69 / 5.2) were measured on these ten; the
    docstring cites bench/conformers.py select(), so the selection must be it."""
    p = em.plan_ensemble(smiles=C18, ensemble=10)
    assert sorted(m["embed_conf_id"] for m in p["members"]) == BENCH_C18_IDS


@pytest.mark.xfail(strict=True, reason="review 6 S4, deferred to section 8: re-spreading "
                   "changes the measured C18 members (test above)")
def test_a_collision_with_the_lowest_energy_member_re_spreads():
    """1-octanol, K = 3: the lowest-energy conformer is the most extended, so
    the Rg-spread pick at the top collides with it. The third member must
    come from the middle of the range, not be the second most compact."""
    p = em.plan_ensemble(smiles="CCCCCCCCO", ensemble=3)
    rg = sorted(m["rg"] for m in p["members"])
    lo, hi = rg[0], rg[-1]
    assert lo + 0.25 * (hi - lo) <= rg[1] <= lo + 0.75 * (hi - lo)


def test_a_capped_symmetry_search_never_keeps_a_symmetry_copy():
    """(tBuOCH2)3CH has 1296 heavy-atom automorphisms, past the cap of 1000."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)(C)OCC(COC(C)(C)C)COC(C)(C)C"))
    *_, conv = em._embed_and_minimise(mol, 30, 0xC0FFEE)
    exact = em.select_ensemble(mol, conv, 30, max_matches=10 ** 6)["n_unique"]
    try:
        capped = em.select_ensemble(mol, conv, 30)
    except ValueError:
        return                                   # a refusal is a correct outcome
    assert capped["n_unique"] == exact


def test_records_keep_their_file_numbers():
    path = _DIR / "E25_record_2_unreadable.sdf"
    m, probs, n = em.read_structure(path, record=2)
    assert m is None and n == 3 and "record 2" in probs[0][2]
    m, probs, n = em.read_structure(path, record=3)
    assert m is not None and m.GetProp("_Name") == "conf 3"


def test_a_generated_ensemble_keeps_the_configuration_its_3d_file_fixes():
    from rdkit.Chem import rdMolTransforms
    path = _DIR / "E19_allene_one_axis_3d.sdf"
    try:
        p = em.plan_ensemble(path=path, ensemble=10)
    except em.EnsembleRefused:
        return                                   # refusing is allowed (E19's "or")
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C=C=CC"))
    signs = set()
    for m in p["members"]:
        rows = [r.split() for r in m["spec"].split(";")]
        c = Chem.Conformer(len(rows))
        for i, r in enumerate(rows):
            c.SetAtomPosition(i, [float(x) for x in r[1:]])
        signs.add(rdMolTransforms.GetDihedralDeg(c, 1, 3, 5, 6) > 0)
    assert signs == {True}


def test_the_gate_leaves_a_small_ring_trans_bond_open_and_nothing_else():
    """S1, gate half: no SMILES fixes trans-cyclooctene's helicity."""
    def left_open(smiles):
        return [k for k, _ in em.stereo_left_open(Chem.MolFromSmiles(smiles))]
    assert any("trans double bond" in k for k in left_open("C1CC/C=C/CCC1"))
    assert any("trans double bond" in k for k in left_open("O[C@H]1CC/C=C/CCC1"))
    assert left_open("C1CC/C=C\\CCC1") == []                 # cis-cyclooctene
    assert left_open("C1CCC/C=C/CCCC1") == []                 # a ten-ring: no helicity
    assert left_open("c1ccccc1") == [] and left_open("C1=CCCCC1") == []


def test_configuration_signs_tell_the_two_helicities_and_axes_apart():
    for name in ("E20_tco_two_helicities.sdf", "E21_allene_two_axes.sdf"):
        ms = list(Chem.SDMolSupplier(str(_DIR / name), removeHs=False))
        signs = [em.configuration_signs(Chem.RemoveHs(m)) for m in ms]
        assert len(signs[0]) == 1 and signs[0] != signs[-1], name
    assert em.configuration_signs(Chem.RemoveHs(_poses(1)[0])) == {}   # ibuprofen: none


def test_every_supplied_record_is_gated_whatever_its_place():
    """S2: the radical is refused by the gate in second place as in first."""
    for name in ("E22_one_skeleton_three_charges.sdf", "E23_radical_first.sdf"):
        codes = {c for c, _ in _refusal(path=_DIR / name, ensemble="supplied")}
        assert "OPEN_SHELL" in codes, name
    codes = {c for c, _ in _refusal(path=_DIR / "E22_one_skeleton_three_charges.sdf",
                                    ensemble="supplied")}
    assert "MEMBER_CHARGE" in codes


def test_the_plan_carries_the_charge_and_spin_every_member_shares():
    p = em.plan_ensemble(path=_DIR / "E7_poses_one_ligand.sdf", ensemble="supplied")
    assert (p["identity"]["charge"], p["identity"]["n_radical"]) == (0, 0)


def test_a_generated_ensemble_from_3d_drops_the_other_configuration():
    p = em.plan_ensemble(path=_DIR / "E19_allene_one_axis_3d.sdf", ensemble=10)
    assert p["protocol"]["n_other_configuration"] > 0


def test_an_unreadable_record_is_named_in_a_supplied_ensemble():
    got = _supplied("E25_record_2_unreadable.sdf")
    assert [c for c, _ in got] == ["RECORD_UNREADABLE"] and got[0][1].startswith("record 2 ")


def test_two_records_with_one_unreadable_are_still_two_records():
    """Invariant 9: the unreadable second record no longer hides the first's
    company (it was read as a single-record file)."""
    m, probs, n = em.read_structure(_DIR / "E25_record_2_unreadable.sdf")
    assert m is None and n == 3 and probs[0][1] == "MULTIPLE_RECORDS"


@pytest.mark.parametrize("smiles", ["Cl/C1=C/CCCCCC1", "O/C1=C/CCCCCC1",
                                    "C/C1=C/CCCCCC1", "C/C1=C(C)/CCCCCC1",
                                    "C1=C\\CCCCC1"])
def test_a_specified_cis_ring_bond_written_as_E_is_not_left_open(smiles):
    """1-chloro- and 1-hydroxy-cis-cyclooctene: the substituent outranks the
    ring carbon, so the CIS ring bond is labelled E. Fully specified, it
    must pass --require-stereo."""
    assert em.stereo_left_open(Chem.MolFromSmiles(smiles)) == []


@pytest.mark.parametrize("smiles", ["Cl/C1=C\\CCCCCC1", "C/C1=C\\CCCCCC1", "C/C1=C(C)\\CCCCCC1",
                                    "C1CC/C=C/CCC1", "O[C@H]1CC/C=C/CCC1"])
def test_a_ring_trans_bond_is_left_open_whatever_its_cip_label(smiles):
    """The trans rings, labelled E or Z: 1-chloro-trans-cyclooctene is Z."""
    assert any("trans double bond" in k for k, _ in em.stereo_left_open(Chem.MolFromSmiles(smiles)))


def test_the_ring_reading_agrees_with_an_embedded_geometry():
    """Truth from 3D, as in Arrowfield's probe: the in-ring C-C=C-C torsion."""
    from rdkit.Chem import rdMolTransforms
    for smiles in ("Cl/C1=C/CCCCCC1", "Cl/C1=C\\CCCCCC1", "C/C1=C/CCCCCC1",
                   "C/C1=C\\CCCCCC1", "O/C1=C/CCCCCC1", "C1CC/C=C/CCC1", "C/C1=C(C)\\CCCCCC1",
                   "C/C1=C(C)/CCCCCC1"):
        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        em._embed_and_minimise(m, 1, 7)
        b = next(x for x in m.GetBonds() if x.GetBondType() == Chem.BondType.DOUBLE)
        _, quad = em._ring_neighbours(m, b)
        trans_3d = abs(rdMolTransforms.GetDihedralDeg(m.GetConformer(), *quad)) > 90
        assert bool(em.small_ring_trans_bonds(Chem.MolFromSmiles(smiles))) == trans_3d, smiles


# --- review 6 S12-S14 ------------------------------------------------------------
def test_a_molecule_without_heavy_atoms_does_not_crash():
    try:
        p = em.plan_ensemble(smiles="[H][H]", ensemble=1)
    except em.EnsembleRefused:
        return
    assert p["members"][0]["rg"] is None and p["members"][0]["folded"] is False


def test_a_supplied_file_is_parsed_once(monkeypatch):
    real = Chem.SDMolSupplier
    calls = []

    def counting(*a, **k):
        calls.append(a)
        return real(*a, **k)
    monkeypatch.setattr(Chem, "SDMolSupplier", counting)
    em.plan_ensemble(path=_DIR / "E7_poses_one_ligand.sdf", ensemble="supplied")
    assert len(calls) <= 2                   # the gate's probe read, and one for the members


def test_the_stereo_backstop_holds_when_verify_computed_accepts_some_completion(monkeypatch):
    """S14 (M1): verify_computed accepts SOME completion when the approval
    leaves stereo open; the members must still agree with record 1."""
    real = em.verify_computed

    def connectivity_only(spec, approved, **kw):
        ok, computed, why = real(spec, approved, **kw)
        return (True, computed, "") if "stereo" in why else (ok, computed, why)
    monkeypatch.setattr(em, "verify_computed", connectivity_only)
    assert [c for c, _ in _supplied("E15_one_centre_inverted.sdf")] == ["MEMBER_STEREO"]


def test_a_generated_member_that_fails_identity_is_refused(monkeypatch):
    """S14 (M3): the per-member identity check on generated members."""
    monkeypatch.setattr(em, "verify_computed", lambda spec, approved, **kw: (False, "CC", "x"))
    assert "MEMBER_IDENTITY" in {c for c, _ in _refusal(smiles="CCO", ensemble=2)}


def test_the_contact_rule_is_the_spec_s_numbers():
    """S14 (M7, M8): docs/ensemble-spec.md section 3 states >= 6 bonds apart and
    Bondi + 0.5 A (measured, section 8.2); a test that reads the module's
    constants checks nothing."""
    assert (em.CONTACT_MIN_BONDS, em.CONTACT_SLACK) == (6, 0.5)
    # a 5-bond pair is never a contact, however close: pentan-1-ol's C1...O
    m5 = Chem.AddHs(Chem.MolFromSmiles("CCCCCO"))
    em._embed_and_minimise(m5, 1, 1)
    c5 = m5.GetConformer()
    o5 = next(a.GetIdx() for a in m5.GetAtoms() if a.GetSymbol() == "O")
    p5 = c5.GetPositions()
    c5.SetAtomPosition(0, (p5[o5] + (p5[0] - p5[o5]) / np.linalg.norm(p5[0] - p5[o5]) * 3.0).tolist())
    assert (0, o5) not in {(i, j) for i, j, *_ in em.intramolecular_contacts(m5)}
    # and on a geometry: hexan-1-ol's C1...O pair is 6 bonds apart
    m = Chem.AddHs(Chem.MolFromSmiles("CCCCCCO"))
    em._embed_and_minimise(m, 1, 1)
    conf = m.GetConformer()
    o = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O")
    pos = conf.GetPositions()
    v = pos[0] - pos[o]
    for d, contact in ((3.70, True), (3.74, False)):     # limit 1.70 + 1.52 + 0.5 = 3.72
        conf.SetAtomPosition(0, (pos[o] + v / np.linalg.norm(v) * d).tolist())
        got = {(i, j) for i, j, *_ in em.intramolecular_contacts(m)}
        assert ((0, o) in got) is contact, d



# --- review 6 S10, S11: guards (what must still be accepted) -------------------------
def test_an_achiral_molecules_mirror_conformers_are_not_refused(tmp_path):
    """S10 checks prochiral labels only in a chiral molecule: a conformer and
    its mirror image are both valid members of an achiral one, and the mirror
    swaps enantiotopic labels."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)CO"))
    em._embed_and_minimise(mol, 2, 3)
    mirror = Chem.Mol(mol, confId=0)
    c = mirror.GetConformer(0)
    for i, p in enumerate(c.GetPositions()):
        c.SetAtomPosition(i, (p[0], p[1], -p[2]))
    path = _write(tmp_path / "mirror.sdf", [Chem.Mol(mol, confId=0), mirror])
    assert em.plan_ensemble(path=path, ensemble="supplied")["protocol"]["k_obtained"] == 2


def test_two_unswapped_conformers_of_a_chiral_molecule_are_accepted(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C(C)C"))
    *_, conv = em._embed_and_minimise(mol, 4, 11)
    path = _write(tmp_path / "ok.sdf", [Chem.Mol(mol, confId=c) for c, _ in conv[:2]])
    assert em.plan_ensemble(path=path, ensemble="supplied")["protocol"]["k_obtained"] == 2


def test_two_conformers_of_one_isotopologue_are_accepted(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("[2H]C([2H])([2H])O"))
    *_, conv = em._embed_and_minimise(mol, 4, 5)
    path = _write(tmp_path / "d3.sdf", [Chem.Mol(mol, confId=c) for c, _ in conv[:2]])
    assert em.plan_ensemble(path=path, ensemble="supplied")["protocol"]["k_obtained"] == 2


def test_tert_butyl_rotation_is_not_a_label_swap(tmp_path):
    """Three equivalent groups rotate into each other: no fixed pro-label."""
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C(C)(C)C"))
    *_, conv = em._embed_and_minimise(mol, 6, 2)
    path = _write(tmp_path / "tbu.sdf", [Chem.Mol(mol, confId=c) for c, _ in conv[:3]])
    assert em.plan_ensemble(path=path, ensemble="supplied")["protocol"]["k_obtained"] == 3
