"""Fast tests: no SCF, no DFT. These run in seconds and gate every push.

Every case here is a defect that actually shipped at some point, lifted from
docs/adversarial-review-1.md and docs/adversarial-review-2.md. A test that was
never red is not evidence, so each one names the finding it guards.
"""
import inspect
import sys

import numpy as np
import pytest

import esp_map as em


# --------------------------------------------------------------- colormaps
def _lum(c):
    return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]


def test_every_colormap_is_oriented_correctly():
    """Review 1: three sequential maps shipped inverted, one of them
    introduced by a 'fix' made on reasoning rather than measurement."""
    for name, cmap in sorted(em.build_maps().items()):
        lo, hi = np.array(cmap(0.0)[:3]), np.array(cmap(1.0)[:3])
        if name in em.SEQUENTIAL:
            assert _lum(lo) < _lum(hi), f"{name}: V_min must be the darkest end"
        else:
            assert lo[0] > lo[2], f"{name}: negative end must be warm"
            assert hi[2] > hi[0], f"{name}: positive end must be cool"


def test_diverging_maps_have_low_chroma_at_zero():
    """The whole colormap argument in one assertion: a diverging scale must be
    least saturated where the data is least meaningful."""
    from colorspacious import cspace_convert
    good = {"vik", "RdBu", "PuOr", "pindr", "broc", "cork"}
    for name in good & set(em.build_maps()):
        rgb = np.array([em.build_maps()[name](0.5)[:3]])
        j = cspace_convert(np.clip(rgb, 0, 1), "sRGB1", "CAM02-UCS")[0]
        assert np.hypot(j[1], j[2]) < 6.0, f"{name}: chroma at V=0 too high"


# ------------------------------------------------------------------ scales
def _anion():
    return np.linspace(-162.0, -44.0, 200)


def _neutral():
    return np.linspace(-12.0, 24.0, 200)


def test_sign_changing_field_gets_a_symmetric_diverging_scale():
    vlim, kind, _ = em.choose_scale(_neutral(), kind="auto")
    assert kind == "diverging"
    assert vlim[0] == pytest.approx(-vlim[1])


def test_ion_switches_to_a_sequential_scale():
    """Review 1: 'symmetric about zero, always' is wrong for an ion, whose
    surface never changes sign."""
    vlim, kind, notes = em.choose_scale(_anion(), kind="auto")
    assert kind == "sequential"
    assert vlim == pytest.approx((-162.0, -44.0))
    assert any("does not change sign" in n for n in notes)


@pytest.mark.parametrize("ratio, expect", [
    (0.01, "sequential"), (0.04, "sequential"), (0.05, "sequential"),
    (0.06, "diverging"), (0.5, "diverging"),
])
@pytest.mark.parametrize("sign", [+1, -1])      # cation with a negative patch, and the mirror
def test_a_sliver_of_the_other_sign_does_not_make_a_diverging_scale(ratio, expect, sign):
    """Review 3 H11: -0.01..+120 was diverging and +0.01..+120 sequential, so
    0.02 kcal/mol decided the picture. Now the minority sign must exceed 5 %
    of the larger extremum, and the notes print the ratio."""
    V = sign * np.array([-ratio * 120.0, 3.0, 120.0])
    assert em.sign_change_ratio(V) == pytest.approx(ratio)
    _, kind, notes = em.choose_scale(V, kind="auto")
    assert kind == expect
    if expect == "sequential":
        assert any(f"{ratio:.3f}" in n for n in notes)


def test_sign_change_ratio_is_zero_when_the_field_has_one_sign():
    assert em.sign_change_ratio(_anion()) == 0.0
    assert em.sign_change_ratio(np.array([0.0, 5.0, 120.0])) == 0.0


def test_default_colormap_follows_the_scale_kind():
    """Review 3, H3: --cmap defaulted to vik on every scale, so an anion's
    sequential scale was drawn diverging and its -80 kcal/mol end came out
    dark blue, the colour every caption reads as positive."""
    assert em.resolve_cmap("sequential") in em.SEQUENTIAL
    assert em.resolve_cmap("diverging") not in em.SEQUENTIAL
    _, kind, _ = em.choose_scale(_anion())
    assert em.resolve_cmap(kind) in em.SEQUENTIAL


@pytest.mark.parametrize("kind, cmap", [
    ("sequential", "vik"),       # diverging map on an ion: the H3 defect
    ("sequential", "RdBu"),
    ("diverging", "lajolla"),    # sequential map on a sign-changing field hides the sign
])
def test_mismatched_colormap_is_refused_unless_forced(kind, cmap):
    with pytest.raises(ValueError, match="does not fit"):
        em.resolve_cmap(kind, cmap)
    assert em.resolve_cmap(kind, cmap, force=True) == cmap


@pytest.mark.parametrize("kind, cmap", [("sequential", "batlow"), ("diverging", "PuOr")])
def test_matching_explicit_colormap_is_kept(kind, cmap):
    assert em.resolve_cmap(kind, cmap) == cmap


_MAPS = {"vik", "RdBu", "lajolla", "batlow", "viridis"}


@pytest.mark.parametrize("scale, cmap, force, refused", [
    ("sequential", "vik", False, True),    # kind known: refuse before the SCF
    ("diverging", "lajolla", False, True),
    ("sequential", "vik", True, False),    # --force-cmap
    ("auto", "vik", False, False),         # kind unknown until the SCF: decide later
    ("auto", "nonsense", False, True),     # unknown map: refuse before the SCF
    ("sequential", None, False, False),
])
def test_cmap_preflight_refuses_what_it_can_before_the_scf(scale, cmap, force, refused):
    """#11 follow-up: resolve_cmap ran after the SCF, so `--cmap vik` on an
    ion did the whole calculation (15 min on a C18) and then wrote nothing."""
    err = em.cmap_preflight(scale, cmap, force, available=_MAPS)
    assert bool(err) is refused


def test_mismatch_found_after_the_scf_falls_back_and_says_so():
    name, warning = em.cmap_after_scf("sequential", "vik")
    assert name == "lajolla" and "only known after the SCF" in warning
    assert em.cmap_after_scf("sequential", "vik", force=True) == ("vik", None)
    assert em.cmap_after_scf("diverging", None) == ("vik", None)


def test_every_caption_uses_one_reading_rule():
    """#11 follow-up: the notebook's static image kept 'warm = negative' on
    ions after the 3D colourbar was fixed."""
    from pathlib import Path
    assert em.scale_caption("diverging") == "warm = negative"
    assert "warm" not in em.scale_caption("sequential")
    nb = (Path(__file__).resolve().parents[1] / "notebooks" / "build_notebook.py").read_text(
        encoding="utf-8")
    display = nb[nb.index("5 · Display"):nb.index("6 · Package")]
    assert "em.scale_caption(kind)" in display and "warm = negative" not in display


def _diag(**over):
    d = dict(closed=True, open_edges=0, nonmanifold_edges=0, components=1, genus=0,
             chi=2, face_max_rho=1e-6, margin=5.0, strict=True, problems=[])
    d.update(over)
    return d


def test_surface_report_tells_clipped_from_folded():
    """Review 3, H4: the JSON kept only chi, and chi = -2 is both a clipped
    surface run with --no-strict (V_min 23 % wrong) and a valid genus-2 fold."""
    clipped = em.surface_report(_diag(closed=False, open_edges=266, genus=None, chi=-2,
                                      face_max_rho=4e-3, margin=2.0, strict=False,
                                      problems=["... CLIPPED ...", "... OPEN ..."]), 0.002)
    folded = em.surface_report(_diag(genus=2, chi=-2), 0.002)
    assert clipped["chi"] == folded["chi"] == -2          # the old JSON's only field
    assert clipped["valid"] is False and clipped["problems"] and clipped["strict"] is False
    assert folded["valid"] is True and folded["problems"] == [] and folded["genus"] == 2
    assert em.surface_report(_diag(), 0.002)["valid"] is True


def test_surface_report_states_the_surface_it_is_on():
    """Spec section 6: the record itself says every area is on the density
    isosurface, not a solvent-accessible surface."""
    d = em.surface_report(_diag(), 0.001)["definition"]
    assert d == {"property": "electron_density", "isovalue_e_bohr3": 0.001}


@pytest.mark.parametrize("gate, verified, problems, expect", [
    (None, None, [], False),                 # no gate ran: never a pass (review 4 N3)
    ({"verdict": "INFO", "forced": False}, True, [], True),
    ({"verdict": "BLOCK", "forced": True}, True, [], False),       # --force
    ({"verdict": "INFO", "forced": False}, False, [], False),      # H2 mismatch
])
def test_validity_input(gate, verified, problems, expect):
    v = em.validity(gate, verified, {"problems": problems}, electronic={"problems": []})
    assert v["input"] is expect and v["surface"] is True and v["sites"] is False


def test_validity_surface_is_false_when_the_surface_has_problems():
    v = em.validity(None, None, em.surface_report(_diag(problems=["CLIPPED"]), 0.002),
                    electronic={"problems": []})
    assert v["surface"] is False


def test_gate_report_marks_a_forced_rejection():
    """Review 3, H4: a structure the gate REJECTED and --force computed was
    written to the JSON as approved_smiles, with no verdict."""
    findings = [(em.BLOCK, "MULTIPLE_FRAGMENTS", "msg", "fix")]
    g = em.gate_report(em.BLOCK, findings, forced=True)
    assert g == {"verdict": "BLOCK", "forced": True,
                 "findings": [["BLOCK", "MULTIPLE_FRAGMENTS"]]}
    # --force on a structure the gate accepted overrides nothing
    assert em.gate_report(em.WARN, [(em.WARN, "X", "m", "f")], forced=True)["forced"] is False


def test_lim_does_not_reimpose_symmetry_on_a_sequential_scale():
    """Review 2, G6: --lim 200 on an anion gave -200..+200, more than half the
    ramp unreachable -- the exact failure the sequential branch exists for."""
    vlim, kind, notes = em.choose_scale(_anion(), lim=200.0, kind="auto")
    assert kind == "sequential"
    assert vlim[1] < 0.0, "upper bound must stay inside the data's sign"
    assert any("does not apply here" in n for n in notes)


def test_explicit_sequential_bounds_are_honoured():
    vlim, kind, _ = em.choose_scale(_anion(), vmin=-170.0, vmax=-40.0, kind="auto")
    assert vlim == pytest.approx((-170.0, -40.0))


def test_inverted_sequential_bounds_are_rejected():
    with pytest.raises(SystemExit):
        em.choose_scale(_anion(), vmin=-40.0, vmax=-170.0, kind="auto")


# ---------------------------------------------------------- structure gate
def _verdict(smiles, **kw):
    from rdkit import Chem
    mol, problems = em.parse_for_check(smiles=smiles)
    if mol is None:
        return em.BLOCK, problems
    return em.check_structure(Chem.RemoveHs(mol), **kw)


def test_clean_structure_is_accepted_with_nothing_to_flag():
    verdict, findings = _verdict("Clc1ccccc1")
    assert verdict == em.INFO and findings == []


@pytest.mark.parametrize("smiles,code", [
    ("*c1ccccc1", "WILDCARD_ATOM"),
    ("C[CH2]", "OPEN_SHELL"),
    ("CC(=O)[O-].[Na+]", "MULTIPLE_FRAGMENTS"),
])
def test_uncomputable_structures_are_blocked(smiles, code):
    verdict, findings = _verdict(smiles)
    assert verdict == em.BLOCK
    assert any(c == code for _, c, _, _ in findings)


def test_valence_error_reports_a_verdict_not_a_traceback():
    """Review 1: this used to raise a raw RDKit traceback out of main()."""
    verdict, findings = _verdict("C[N](C)(C)(C)C")
    assert verdict == em.BLOCK
    assert any("Valence" in c for _, c, _, _ in findings)


def test_radical_is_allowed_once_spin_is_given():
    """Review 2, G1: the gate and the solver must agree about what is
    computable. With --spin the pipeline is unrestricted, so this is a note."""
    verdict, findings = _verdict("C[CH2]", spin=1)
    assert verdict != em.BLOCK
    assert any(c == "OPEN_SHELL_UNRESTRICTED" for _, c, _, _ in findings)


def test_mismatched_spin_is_flagged():
    """Was an INFO note; review 5 R2 made it a WARN that reaches the record."""
    _, findings = _verdict("C[CH2]", spin=3)
    assert [s for s, c, _, _ in findings if c == "SPIN_OVERRIDE"] == [em.WARN]


@pytest.mark.parametrize("smiles,n_centres", [("CC(N)C(=O)O", 1),
                                              ("CC(N)C(O)C(=O)O", 2)])
def test_unspecified_stereocentres_are_reported(smiles, n_centres):
    _, findings = _verdict(smiles)
    hit = [f for f in findings if f[1] == "UNSPECIFIED_TETRAHEDRAL"]
    assert hit, "an undrawn stereocentre must be reported"
    # one centre is a note (enantiomers are mirror-identical); two are a
    # warning (diastereomers are different molecules)
    assert hit[0][0] == (em.INFO if n_centres == 1 else em.WARN)


def test_undrawn_double_bond_geometry_is_warned():
    _, findings = _verdict("CC=CC")
    assert any(c == "UNSPECIFIED_DOUBLE_BOND" for _, c, _, _ in findings)


def test_specified_stereo_is_not_flagged():
    """The check must be able to pass: a fully specified structure is clean."""
    _, findings = _verdict("C/C=C/C")
    assert not any(c.startswith("UNSPECIFIED") for _, c, _, _ in findings)


def test_charge_is_noted_so_the_scale_can_switch():
    _, findings = _verdict("c1ccccc1C(=O)[O-]")
    assert any(c == "CHARGED" for _, c, _, _ in findings)


# ------------------------------------------------------------------- input
def test_2d_depiction_is_detected_and_its_coordinates_discarded(tmp_path):
    """The single most dangerous input: a drawn structure is a flat picture,
    and feeding it to an SCF gives a planar molecule and a wrong map."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles("Clc1ccccc1")
    AllChem.Compute2DCoords(m)
    p = tmp_path / "drawn.mol"
    p.write_text(Chem.MolToMolBlock(m))
    _, info = em.from_file(str(p))
    assert "DISCARDED" in info["source"]


def test_real_3d_coordinates_are_used_as_is(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("Clc1ccccc1"))
    ps = AllChem.ETKDGv3(); ps.randomSeed = 42
    AllChem.EmbedMolecule(m, ps); AllChem.MMFFOptimizeMolecule(m)
    p = tmp_path / "real.mol"
    p.write_text(Chem.MolToMolBlock(m))
    _, info = em.from_file(str(p))
    assert "as-is" in info["source"]


_ROOT = __import__("pathlib").Path(__file__).resolve().parent.parent


def _flat_nitrite_mol(tmp_path, header="3D", length=1.258):
    """A nitrite in the xy plane (115 deg), Kekule: N=O and N-O(-). One
    element pair. At 1.258 A it is the DFT geometry; at 1.5 A a drawing."""
    import math
    from rdkit import Chem
    m = Chem.RWMol(Chem.MolFromSmiles("N(=O)[O-]"))
    c = Chem.Conformer(3)
    h = math.radians(115 / 2)
    c.SetAtomPosition(0, (0.0, 0.0, 0.0))
    c.SetAtomPosition(1, (length * math.sin(h), length * math.cos(h), 0.0))
    c.SetAtomPosition(2, (-length * math.sin(h), length * math.cos(h), 0.0))
    c.Set3D(header == "3D")
    m.AddConformer(c)
    p = tmp_path / "nitrite.mol"
    p.write_text(Chem.MolToMolBlock(m))
    return p


@pytest.mark.parametrize("path, expect", [
    # review 3 H5: real planar geometries were thrown away as "drawings"
    ("docs/review-3/fx/chlorobenzene_3d_planar.mol", "as-is"),
    ("docs/review-3/fx/water_xy.xyz", "as-is"),
    ("docs/review-3/fx/chlorobenzene_dft_like.xyz", "as-is"),
    # drawings must still be caught, whatever the header claims
    ("docs/review-3/fx/rd2d.mol", "DISCARDED"),
    ("docs/review-3/fx/rd2d_flagged3d.mol", "DISCARDED"),
    ("examples/drawn_2d.mol", "DISCARDED"),
    ("examples/drawn.cdxml", "DISCARDED"),
    # unchanged behaviour on ordinary 3D input
    ("examples/real_3d.mol", "as-is"),
    ("examples/real_3d_noH.mol", "as-is"),
    ("examples/water.xyz", "as-is"),
])
def test_planar_geometry_is_kept_and_drawings_are_discarded(path, expect):
    _, info = em.from_file(str(_ROOT / path))
    assert expect in info["source"], (path, info["source"])


def test_symmetric_flat_ion_at_a_drawing_length_is_refused_not_guessed(tmp_path):
    """Kekule nitrite has two bond ORDERS but one element pair. Drawn at a
    drawing program's 1.5 A it looks like a symmetric molecule: refused.
    Counting orders once called the DFT geometry a drawing and MMFF replaced
    it with N-O 1.42 A."""
    p = _flat_nitrite_mol(tmp_path, length=1.5)
    mol, probs = em.parse_for_check(path=str(p), charge=-1)
    assert mol is None and probs[0][1] == "DIMENSION_AMBIGUOUS"
    with pytest.raises(ValueError):
        em.from_file(str(p))


def test_symmetric_flat_ion_at_real_lengths_is_kept(tmp_path):
    """Review 5 R11: the DFT nitrite (1.258 A, within 10 % of N-O and N=O)
    is a geometry; v0.2.0 refused it as DIMENSION_AMBIGUOUS."""
    _, info = em.from_file(str(_flat_nitrite_mol(tmp_path)))
    assert "as-is" in info["source"], info["source"]


def test_symmetric_flat_ion_with_a_2d_header_is_a_drawing(tmp_path):
    _, info = em.from_file(str(_flat_nitrite_mol(tmp_path, header="2D")))
    assert "DISCARDED" in info["source"]


def _multi_sdf(tmp_path, smiles=("CCO", "CCCl", "c1ccccc1")):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    p = tmp_path / "poses.sdf"
    w = Chem.SDWriter(str(p))
    for s in smiles:
        m = Chem.AddHs(Chem.MolFromSmiles(s))
        AllChem.EmbedMolecule(m, randomSeed=1)
        w.write(m)
    w.close()
    return p


def test_multi_record_file_is_refused_without_record(tmp_path):
    """One geometry in, one record out: 20 docking poses used to become an
    analysis of pose 1 with a warning nobody reads in a pipeline."""
    mol, probs = em.parse_for_check(path=str(_multi_sdf(tmp_path)))
    assert mol is None and probs[0][1] == "MULTIPLE_RECORDS"


def test_record_picks_one_record(tmp_path):
    from rdkit import Chem
    p = _multi_sdf(tmp_path)
    mol, probs = em.parse_for_check(path=str(p), record=2)
    assert not probs and Chem.MolToSmiles(Chem.RemoveHs(mol)) == "CCCl"
    _, info = em.from_file(str(p), record=2)
    assert info["record"] == {"index": 2, "of": 3}


@pytest.mark.parametrize("record", [0, 4])
def test_record_out_of_range_is_refused(tmp_path, record):
    mol, probs = em.parse_for_check(path=str(_multi_sdf(tmp_path)), record=record)
    assert mol is None and probs[0][1] == "RECORD_OUT_OF_RANGE"


def test_cli_refuses_multi_record_before_compute_and_accepts_record(tmp_path, monkeypatch):
    p = str(_multi_sdf(tmp_path))
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--file", p])
    with pytest.raises(SystemExit):
        em.main()
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--file", p, "--record", "3"])
    with pytest.raises(_ReachedCompute):
        em.main()


def _spec_xyz(spec):
    return np.array([[float(x) for x in t.split()[1:]] for t in spec.split(";")])


def test_supplied_hydrogens_are_frozen_when_others_are_added(tmp_path):
    """A docking pose often keeps its polar hydrogens and drops the rest. Only
    the added hydrogens may move: relaxing the supplied O-H in isolation
    replaces an orientation set by the pocket with the gas-phase optimum."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    m = Chem.AddHs(Chem.MolFromSmiles("CCCO"))
    ps = AllChem.ETKDGv3(); ps.randomSeed = 42
    AllChem.EmbedMolecule(m, ps); AllChem.MMFFOptimizeMolecule(m)
    o = next(a for a in m.GetAtoms() if a.GetSymbol() == "O")
    h_o = next(n.GetIdx() for n in o.GetNeighbors() if n.GetSymbol() == "H")
    c1 = next(n.GetIdx() for n in o.GetNeighbors() if n.GetSymbol() == "C")
    c2 = next(n.GetIdx() for n in m.GetAtomWithIdx(c1).GetNeighbors()
              if n.GetSymbol() == "C")
    # an O-H rotamer the isolated molecule would not choose
    rdMolTransforms.SetDihedralDeg(m.GetConformer(), c2, c1, o.GetIdx(), h_o, 120.0)
    rw = Chem.RWMol(m)
    for i in sorted((a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "H"
                     and a.GetNeighbors()[0].GetSymbol() == "C"), reverse=True):
        rw.RemoveAtom(i)
    pose = rw.GetMol()
    p = tmp_path / "pose.mol"
    p.write_text(Chem.MolToMolBlock(pose))

    spec, info = em.from_file(str(p))
    xyz = _spec_xyz(spec)
    supplied = Chem.MolFromMolFile(str(p), removeHs=False).GetConformer().GetPositions()
    assert len(xyz) == m.GetNumAtoms()
    # every supplied atom, the O-H hydrogen included, is where the file put it
    assert np.abs(xyz[:len(supplied)] - supplied).max() < 1e-5
    assert "7 hydrogen(s) added" in info["source"]


def test_xyz_connectivity_is_perceived(tmp_path):
    """Review 2, G3: xyz carries no bonds, so every xyz file was rejected as a
    salt with a remedy about counter-ions."""
    p = tmp_path / "w.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.117\nH 0.0 0.757 -0.469\n"
                 "H 0.0 -0.757 -0.469\n")
    mol, problems, _ = em.read_structure(str(p))
    assert mol is not None and not problems
    assert mol.GetNumBonds() == 2


# chlorobenzene, all 12 atoms, docking-style HETATM records, no CONECT
_CLPH_PDB = """\
HETATM    1 Cl1  LIG A   1      -3.060  -0.000  -0.000  1.00  0.00          Cl
HETATM    2 C2   LIG A   1      -1.339  -0.000  -0.000  1.00  0.00           C
HETATM    3 C3   LIG A   1      -0.647   1.192   0.210  1.00  0.00           C
HETATM    4 C4   LIG A   1       0.749   1.191   0.210  1.00  0.00           C
HETATM    5 C5   LIG A   1       1.446   0.000  -0.000  1.00  0.00           C
HETATM    6 C6   LIG A   1       0.749  -1.191  -0.210  1.00  0.00           C
HETATM    7 C7   LIG A   1      -0.647  -1.192  -0.210  1.00  0.00           C
HETATM    8 H8   LIG A   1      -1.186   2.121   0.374  1.00  0.00           H
HETATM    9 H9   LIG A   1       1.293   2.118   0.373  1.00  0.00           H
HETATM   10 H10  LIG A   1       2.533   0.000  -0.000  1.00  0.00           H
HETATM   11 H11  LIG A   1       1.293  -2.118  -0.373  1.00  0.00           H
HETATM   12 H12  LIG A   1      -1.186  -2.121  -0.374  1.00  0.00           H
END
"""


def test_pdb_without_bond_records_is_read_as_the_right_molecule(tmp_path):
    """Review 3, H1 (BLOCKER): proximity bonding made every bond single, so
    chlorobenzene was read as chlorocyclohexane (ClC1CCCCC1), AddHs added six
    hydrogens to a flat ring, and the gate reported nothing to flag."""
    from rdkit import Chem
    p = tmp_path / "clph.pdb"
    p.write_text(_CLPH_PDB)
    mol, problems, _ = em.read_structure(str(p))
    assert not problems
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == "Clc1ccccc1"
    _, info = em.from_file(str(p))
    assert info["formula"] == "C6H5Cl"
    assert "hydrogen(s) added" not in info["source"]


def test_heavy_atom_only_pdb_is_refused_not_misread(tmp_path):
    """Bond orders cannot be perceived without the hydrogens. Refusing is
    correct; guessing reproduces H1."""
    p = tmp_path / "heavy.pdb"
    p.write_text("".join(l + "\n" for l in _CLPH_PDB.splitlines()
                         if l.startswith("HETATM") and not l.endswith("H")) + "END\n")
    mol, problems, _ = em.read_structure(str(p))
    assert mol is None and problems[0][1] == "PDB_CONNECTIVITY"


def test_h_count_guard_catches_a_reader_that_guesses_wrong(tmp_path, monkeypatch):
    """Independent of the reader: replay the OLD PDB reader (proximity bonds,
    all single) and check from_file refuses its output by hydrogen count."""
    from rdkit import Chem
    p = tmp_path / "clph.pdb"
    p.write_text(_CLPH_PDB)

    def old_reader(path, charge=0, record=None):
        m = Chem.MolFromPDBFile(str(path), removeHs=False)
        m.SetProp("_bonds", "perceived")
        m.SetProp("_coords", "geometry")
        return m, [], 1
    monkeypatch.setattr(em, "read_structure", old_reader)
    with pytest.raises(ValueError, match="H_COUNT_MISMATCH"):
        em.from_file(str(p))


def test_file_bond_orders_with_only_polar_h_are_accepted(tmp_path):
    """A docking pose may keep only its polar hydrogens. MOL/SDF carry their
    own bond orders, so adding the rest is correct and must not be refused."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    full = Chem.AddHs(Chem.MolFromSmiles("OCc1ccccc1"))
    AllChem.EmbedMolecule(full, randomSeed=1)
    Chem.Kekulize(full, clearAromaticFlags=True)
    rw = Chem.RWMol(full)
    for a in sorted((a.GetIdx() for a in full.GetAtoms() if a.GetSymbol() == "H"
                     and a.GetNeighbors()[0].GetSymbol() != "O"), reverse=True):
        rw.RemoveAtom(a)
    p = tmp_path / "polarH.mol"
    Chem.MolToMolFile(rw.GetMol(), str(p))
    _, info = em.from_file(str(p))
    assert info["formula"] == "C7H8O" and "7 hydrogen(s) added" in info["source"]


def test_charged_xyz_needs_and_accepts_the_charge(tmp_path):
    """Review 2, G10: DetermineBonds hard-coded charge=0 while the error
    message told the user to set the charge explicitly."""
    p = tmp_path / "ac.xyz"
    p.write_text("7\nacetate\nC 0.000 0.000 0.000\nC 1.520 0.000 0.000\n"
                 "O 2.150 1.070 0.000\nO 2.150 -1.070 0.000\n"
                 "H -0.390 1.020 0.000\nH -0.390 -0.510 0.880\n"
                 "H -0.390 -0.510 -0.880\n")
    assert em.read_structure(str(p), charge=0)[0] is None    # must fail
    mol, problems, _ = em.read_structure(str(p), charge=-1)  # must succeed
    assert mol is not None and not problems


# ------------------------------------------------------- approved = computed
# Review 3, H2: the assertion compared the input with itself and never read
# the coordinates PySCF receives. These check the coordinates.
@pytest.mark.parametrize("smiles", [
    "CCO", "CC(=O)[O-]", r"C/C=C\C", "Clc1ccccc1",
    "[N-]=C=S",           # perceived from geometry as N#C[S-]: bond orders are not checked
    "CS(C)=O",            # likewise C[S+](C)[O-]
    "CC(O)C(N)C(=O)O",    # stereo left open: any completion is the approved structure
    "C[CH2]",             # radical
    "[2H]OC",             # isotope: carried by the approval, invisible in coordinates
])
def test_computed_structure_matches_its_own_approval(smiles):
    from rdkit import Chem
    spec, _ = em.from_smiles(smiles)
    ok, computed, reason = em.verify_computed(spec, Chem.MolFromSmiles(smiles))
    assert ok, reason


@pytest.mark.parametrize("approved, computed", [
    ("CCO", "CC=O"),                        # another molecule (H count)
    ("COC", "CCO"),                         # constitutional isomer, same formula
    ("Oc1ncccc1", "O=c1[nH]cccc1"),         # tautomer
    (r"C/C=C\C", "C/C=C/C"),                # E computed for a Z approval
    ("C[C@@H](O)CC", "C[C@H](O)CC"),        # the other enantiomer
])
def test_computing_a_different_structure_is_detected(approved, computed):
    from rdkit import Chem
    spec, _ = em.from_smiles(computed)
    ok, _, reason = em.verify_computed(spec, Chem.MolFromSmiles(approved))
    assert not ok and reason


class _ReachedCompute(Exception):
    pass


def _stop_before_scf(monkeypatch):
    def stop(*a, **k):
        raise _ReachedCompute
    monkeypatch.setattr(em, "estimate_cost", stop)


def test_optimize_without_pyberny_is_refused_before_the_scf(monkeypatch):
    """--optimize imported pyscf.geomopt.berny_solver only after the first SCF,
    and pyberny was in no dependency list: the run computed, then crashed."""
    _stop_before_scf(monkeypatch)
    monkeypatch.setitem(sys.modules, "berny", None)        # import berny -> ImportError
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CCO", "--optimize"])
    with pytest.raises(SystemExit, match="pyberny"):
        em.main()


def test_optimize_with_pyberny_reaches_the_compute(monkeypatch):
    import types
    _stop_before_scf(monkeypatch)
    monkeypatch.setitem(sys.modules, "berny", types.ModuleType("berny"))
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CCO", "--optimize"])
    with pytest.raises(_ReachedCompute):
        em.main()


def test_main_refuses_coordinates_that_are_not_the_approved_structure(monkeypatch):
    """Fault injection at the boundary the old check could not see: the gate
    approves CCO, the geometry step hands PySCF acetaldehyde. Only the
    coordinates are swapped; the info dict still describes CCO, so a check
    that reads anything but the coordinates passes."""
    real = em.from_smiles
    monkeypatch.setattr(em, "from_smiles",
                        lambda s, **k: (real("CC=O", **k)[0], real(s, **k)[1]))
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CCO"])
    with pytest.raises(SystemExit, match="not the one about to be computed"):
        em.main()


def test_main_refuses_a_contradicting_cmap_before_the_scf(tmp_path, monkeypatch):
    """--scale fixes the kind up front, so --cmap vik on it is refused before
    any compute (cmap_preflight wired into main), with no image written."""
    out = tmp_path / "ac.png"
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CC(=O)[O-]",
                                      "--scale", "sequential", "--cmap", "vik",
                                      "--out", str(out)])
    with pytest.raises(SystemExit, match="does not fit"):
        em.main()
    assert not out.exists()


def test_main_lets_a_healthy_structure_through_to_the_scf(monkeypatch):
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CC(=O)[O-]"])
    with pytest.raises(_ReachedCompute):
        em.main()


def test_a_forced_unverifiable_structure_warns_instead_of_refusing(monkeypatch, capsys):
    """A salt is BLOCKed by the gate; --force computes it anyway, and distance
    bonding ties Na+ to the carboxylate, so the identity cannot be confirmed.
    The user has already overridden the gate: say so, do not refuse."""
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CC(=O)[O-].[Na+]",
                                      "--force"])
    with pytest.raises(_ReachedCompute):
        em.main()
    assert "could not be verified" in capsys.readouterr().err


def test_smiles_geometry_is_deterministic():
    assert em.from_smiles("CCO")[0] == em.from_smiles("CCO")[0]


def test_a_different_seed_gives_a_different_conformer():
    """The determinism test must be able to fail."""
    assert em.from_smiles("CCCCO", seed=1)[0] != em.from_smiles("CCCCO", seed=2)[0]


# ------------------------------------------------------------------- mesh
def test_euler_characteristic_of_a_closed_surface():
    """chi = 2 is the clipped-isosurface detector. A tetrahedron is the
    smallest closed genus-0 mesh."""
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    assert em.euler_characteristic(faces) == 2


def test_euler_characteristic_detects_an_open_surface():
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3]])   # one face removed
    assert em.euler_characteristic(faces) != 2


def _torus_faces(nu=12, nv=8):
    """A closed triangulated torus: genus 1, chi = 0."""
    idx = lambda i, j: (i % nu) * nv + (j % nv)
    f = []
    for i in range(nu):
        for j in range(nv):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            f += [[a, b, c], [a, c, d]]
    return np.array(f)


_TET = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])


def test_closed_higher_genus_surface_passes_with_a_note():
    """F9: two compact C18 conformers were refused at chi = -2 / -4 with
    nothing clipped. A closed surface with handles is valid."""
    topo = em.mesh_topology(_torus_faces())
    assert topo == dict(chi=0, components=1, open_edges=0, nonmanifold_edges=0,
                        closed=True, genus=1)
    problems, notes = em.surface_problems(topo, 1e-6, em.ISO_BADER)
    assert problems == [] and "genus 1" in notes[0]


def test_two_closed_pieces_pass_with_a_warning():
    """Genus is summed over pieces (chi = 2C - 2g): two separate tori give
    C = 2, g = 2. Closed, so not refused, but the extra piece is flagged."""
    t = _torus_faces()
    topo = em.mesh_topology(np.vstack([t, t + t.max() + 1]))
    assert (topo["components"], topo["chi"], topo["genus"]) == (2, 0, 2)
    problems, notes = em.surface_problems(topo, 1e-6, em.ISO_BADER)
    assert problems == [] and any("2 separate closed surfaces" in n for n in notes)


def test_open_mesh_is_still_refused():
    """The reason the check exists: a clipped surface is open."""
    topo = em.mesh_topology(_TET[:3])
    assert topo["open_edges"] == 3 and not topo["closed"] and topo["genus"] is None
    problems, _ = em.surface_problems(topo, 1e-6, em.ISO_BADER)
    assert any("OPEN" in p and "not a single closed" in p for p in problems)


def test_density_on_the_box_face_is_refused_even_if_closed():
    problems, _ = em.surface_problems(em.mesh_topology(_TET), 5e-3, em.ISO_BADER)
    assert any("CLIPPED" in p for p in problems)


def test_sphere_like_surface_passes_silently():
    topo = em.mesh_topology(_TET)
    assert (topo["closed"], topo["genus"], topo["components"]) == (True, 0, 1)
    assert em.surface_problems(topo, 1e-6, em.ISO_BADER) == ([], [])


def test_nonmanifold_edge_is_refused():
    fan = np.array([[0, 1, 2], [0, 1, 3], [0, 1, 4]])     # one edge, three triangles
    problems, _ = em.surface_problems(em.mesh_topology(fan), 1e-6, em.ISO_BADER)
    assert any("more than two" in p for p in problems)


# ------------------------------------------------------------------- cost
def test_cost_estimate_needs_no_scf_and_knows_the_spin():
    """Review 2, G8: estimate_cost omitted spin, so the guard meant to prevent
    an expensive mistake was the first thing to crash on a radical."""
    spec, info = em.from_smiles("C[CH2]")
    natm, nao = em.estimate_cost(spec, "def2-svp", charge=0, spin=1)
    assert natm == 7 and nao > 0


def test_gpu_without_gpu4pyscf_fails_before_any_work():
    """--gpu on a machine without GPU4PySCF must say what is missing, not
    build a molecule and die inside to_gpu() with an opaque error."""
    pytest.importorskip("pyscf")
    try:
        import gpu4pyscf  # noqa: F401
        pytest.skip("gpu4pyscf is installed here")
    except ImportError:
        pass
    with pytest.raises(RuntimeError, match="gpu4pyscf"):
        em.run_scf("H 0 0 0; H 0 0 0.74", basis="sto-3g", gpu=True)


@pytest.mark.parametrize("spin, solvent, df, expect", [
    (0, None, True, []),                                  # the validated path
    (0, "vacuum", True, []),
    (0, "water", True, ["PCM"]),
    (2, None, True, ["open shell"]),
    (0, None, False, ["exact integrals"]),
    (1, "water", False, ["open shell", "PCM", "exact integrals"]),
])
def test_gpu_untested_names_each_unvalidated_part(spin, solvent, df, expect):
    got = em.gpu_untested(spin, solvent, df)
    assert len(got) == len(expect) and all(e in g for e, g in zip(expect, got))


class _Stop(Exception):
    pass


def test_run_scf_warns_on_an_untested_gpu_path(monkeypatch, capsys):
    """Review 3 H9 validation: with gpu4pyscf stubbed, run_scf(gpu=True,
    solvent="water") must warn. The notebook's GPU = auto put every GPU
    runtime with a solvent on this path silently."""
    import types
    from pyscf import gto
    monkeypatch.setitem(sys.modules, "gpu4pyscf", types.ModuleType("gpu4pyscf"))

    def stop(*a, **k):
        raise _Stop
    monkeypatch.setattr(gto, "M", stop)
    with pytest.raises(_Stop):
        em.run_scf("H 0 0 0; H 0 0 0.74", basis="sto-3g", gpu=True,
                   solvent="water", density_fit=True)
    err = capsys.readouterr().err
    assert "WARNING: --gpu with a PCM solvent (water) is untested" in err
    with pytest.raises(_Stop):
        em.run_scf("H 0 0 0; H 0 0 0.74", basis="sto-3g", gpu=True, density_fit=True)
    assert "WARNING" not in capsys.readouterr().err


def test_notebook_warns_on_an_untested_gpu_path():
    from pathlib import Path
    nb = (Path(__file__).resolve().parents[1] / "notebooks" / "build_notebook.py").read_text(
        encoding="utf-8")
    run = nb[nb.index("4 · Run"):]
    assert run.index("em.gpu_untested(SPIN, SOLVENT, DENSITY_FIT)") < run.index("cmd = [")


@pytest.mark.parametrize("kind, lim, vmin, vmax, comparable", [
    ("diverging", None, None, None, False),     # auto-scaled to this molecule
    ("diverging", 50.0, None, None, True),      # --lim fixes a diverging scale
    ("diverging", None, -60.0, 60.0, False),    # --vmin/--vmax ignored here
    ("sequential", None, None, None, False),    # the defect: was reported True
    ("sequential", 50.0, None, None, False),    # --lim ignored for sequential
    ("sequential", None, -170.0, None, False),  # one end auto, still per-molecule
    ("sequential", None, -170.0, -30.0, True),  # both ends fixed
])
def test_scale_is_comparable_only_when_the_user_fixed_it(kind, lim, vmin, vmax,
                                                         comparable):
    """An auto-scaled sequential map (one anion at -170..-80, another at
    -130..-30) was written to --json as comparable. Panels are comparable
    only when the limits came from the user, through the flag that applies
    to that scale kind."""
    p = em.scale_provenance(kind, lim, vmin, vmax)
    assert p == {"source": "user" if comparable else "auto",
                 "colour_scale_fixed": comparable}


_A = dict(xc="b3lyp", basis="def2-svp", solvent="vacuum", iso=0.002,
          spacing=0.2, density_fit=False, geometry_method="etkdg_mmff94s")
_R = dict(cmap="vik", kind="diverging", vlim=(-25.0, 25.0), elev=20.0, azim=-60.0,
          size=900, ambient=0.35)


def _keys(a=None, r=None):
    return (em.record_key(em.analysis_fields(**{**_A, **(a or {})})),
            em.record_key(em.render_fields(**{**_R, **(r or {})})))


def test_cmap_changes_only_the_render_key_and_basis_only_the_analysis_key():
    """Spec 8b. Review 3 proposed one compare_key over method AND colour, which
    would call two identical calculations incomparable over a colormap; and
    'comparable' used to depend on the colour flags alone, across bases."""
    a0, r0 = _keys()
    a1, r1 = _keys(r={"cmap": "RdBu"})
    a2, r2 = _keys(a={"basis": "def2-tzvp"})
    assert a1 == a0 and r1 != r0
    assert a2 != a0 and r2 == r0


@pytest.mark.parametrize("change", [
    {"xc": "B3LYP"}, {"basis": "def2-SVP"}, {"solvent": None}, {"solvent": "Vacuum"},
])
def test_one_calculation_gets_one_analysis_key(change):
    """PySCF is case-insensitive: def2-SVP and def2-svp are the same run."""
    assert _keys(a=change)[0] == _keys()[0]


@pytest.mark.parametrize("change", [
    {"solvent": "water"}, {"iso": 0.001}, {"spacing": 0.16}, {"density_fit": True},
    {"xc": "pbe0"}, {"geometry_method": "dft:b3lyp/def2svp"},
    {"geometry_method": "undeclared"},
])
def test_every_method_axis_changes_the_analysis_key(change):
    assert _keys(a=change)[0] != _keys()[0]


_EX = __import__("pathlib").Path(__file__).resolve().parents[1] / "examples"


@pytest.mark.parametrize("make, expect", [
    (lambda: em.from_smiles("CCO")[1], "etkdg_mmff94s"),
    (lambda: em.from_file(str(_EX / "real_3d.mol"))[1], "undeclared"),
    (lambda: em.from_file(str(_EX / "real_3d_noH.mol"))[1], "undeclared+mmff94s_h"),
    (lambda: em.from_file(str(_EX / "drawn_2d.mol"))[1], "etkdg_mmff94s"),
])
def test_geometry_method_names_how_the_coordinates_were_made(make, expect):
    """Review 4 N2: the key omitted the geometry source, the largest measured
    axis on small anions (MMFF vs DFT up to 6.5 kcal/mol)."""
    info = make()
    assert em.geometry_method(info["source"], info["forcefield"]) == expect


def test_geometry_method_for_cli_only_sources():
    assert em.geometry_method(atoms=True) == "undeclared"
    assert em.geometry_method(molecule="chlorobenzene") == "idealised"
    opt = ("B3LYP", "def2-TZVP", "Water")
    assert em.geometry_method("3D generated from SMILES by ETKDGv3",
                              optimize=opt) == "dft:b3lyp/def2tzvp/water"
    assert em.geometry_method(atoms=True, optimize=("b3lyp", "sto-3g", "vacuum")) \
        == "dft:b3lyp/sto3g"
    with pytest.raises(ValueError):
        em.geometry_method("from a source nobody wrote down")


def test_review_glycine_pair_no_longer_shares_a_key():
    """r8 (--optimize) and r9 (no --optimize) shared analysis_key 62f86d4d6b
    across a 25-29 kcal/mol difference in V_min/V_max."""
    import json
    out = __import__("pathlib").Path(__file__).resolve().parents[1] / "docs" / "review-4" / "out"
    keys = []
    for name in ("r8.json", "r9.json"):
        r = json.loads((out / name).read_text())
        s = r["settings"]
        f = dict(r["key_fields"]["analysis"], geometry_method=em.geometry_method(
            "3D generated from SMILES by ETKDGv3", "MMFF94s",
            optimize=(s["xc"], s["basis"], s["solvent"]) if s["optimize"] else None))
        keys.append(em.record_key(f))
    assert keys[0] != keys[1]


def test_different_molecules_made_the_same_way_share_a_key():
    a = em.from_smiles("CCO")[1]
    b = em.from_smiles("c1ccccc1Cl")[1]
    assert em.geometry_method(a["source"], a["forcefield"]) == \
        em.geometry_method(b["source"], b["forcefield"])


def test_auto_scaled_panels_get_different_render_keys():
    assert _keys(r={"vlim": (-30.0, 30.0)})[1] != _keys()[1]


def test_solvent_esp_is_recorded_with_solvent_only():
    assert em.analysis_fields(**{**_A, "solvent": "water"})["solvent_esp"] == em.SOLVENT_ESP[False]
    assert em.analysis_fields(**_A)["solvent_esp"] is None


def test_reaction_field_convention_is_a_different_analysis_key():
    """Review 3 H7: the two conventions give solvent shifts of opposite sign
    on the chlorobenzene sigma-hole. Numbers under one must never be compared
    with numbers under the other."""
    water = {**_A, "solvent": "water"}
    a = em.analysis_fields(**water)
    b = em.analysis_fields(**water, reaction_field=True)
    assert b["solvent_esp"] == em.SOLVENT_ESP[True]
    assert em.record_key(a) != em.record_key(b)
    assert em.analysis_fields(**_A, reaction_field=True)["solvent_esp"] is None


class _FakePCM:
    def __init__(self, q, centres, zeta):
        import types
        self.with_solvent = types.SimpleNamespace(
            _intermediates={"q": np.asarray(q, float)},
            surface={"grid_coords": np.asarray(centres, float),
                     "charge_exp": np.asarray(zeta, float)})


def test_reaction_field_is_coulombic_far_away_and_finite_on_a_charge():
    """Smeared charges q erf(zeta r)/r: a point charge at long range, and
    2 zeta q / sqrt(pi) on top of one, where a bare 1/r would blow up."""
    from math import erf, pi, sqrt
    centres = [[0, 0, 0], [3.0, 0, 0]]
    mf = _FakePCM([0.4, -0.1], centres, [2.0, 1.5])
    far = np.array([[60.0, 0, 0]])
    assert em.reaction_field_at_points(mf, far)[0] == pytest.approx(0.4 / 60 - 0.1 / 57, rel=1e-9)
    on = em.reaction_field_at_points(mf, np.array([[0.0, 0, 0]]))[0]
    assert on == pytest.approx(2 * 2.0 * 0.4 / sqrt(pi) - 0.1 * erf(1.5 * 3.0) / 3.0)


class _FakeSCF:
    def __init__(self, e, occ):
        self.mo_energy, self.mo_occ = np.asarray(e, float), np.asarray(occ, float)


@pytest.mark.parametrize("e, occ, charge, solvent, warns", [
    ([-0.5, 0.02, 0.3], [2, 2, 0], -1, None, True),        # unbound HOMO
    ([-0.5, -0.02, 0.3], [2, 2, 0], -1, None, False),      # bound
    ([-0.5, 0.02, 0.3], [2, 2, 0], -1, "water", False),    # continuum binds it
    ([-0.5, 0.02, 0.3], [2, 2, 0], 0, None, False),        # not an anion
    ([[-0.5, 0.02], [-0.5, 0.3]], [[1, 1], [1, 0]], -1, None, True),   # UKS
])
def test_unbound_gas_phase_anion_is_warned(e, occ, charge, solvent, warns):
    w = em.unbound_anion_warning(_FakeSCF(e, occ), charge, solvent)
    assert bool(w) is warns
    if w:
        assert "--solvent" in w


def test_reaction_field_without_solvent_is_refused_before_compute(monkeypatch):
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CCO",
                                      "--include-reaction-field"])
    with pytest.raises(SystemExit, match="needs --solvent"):
        em.main()


def _fake_dm(mol, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(mol.nao, mol.nao))
    return x @ x.T / mol.nao


def test_ion_dipole_does_not_depend_on_where_the_molecule_sits():
    """Review 3 H6: acetate gave 6.42 D as embedded and 30.35 D translated by
    5 A. About the centre of nuclear charge it must not move. Any fixed AO
    density matrix will do: translation carries the basis with the atoms."""
    from pyscf import gto
    spec, _ = em.from_smiles("CC(=O)[O-]")
    shifted = "; ".join(f"{t.split()[0]} " + " ".join(
        f"{float(x) + 5.0:.10f}" for x in t.split()[1:]) for t in spec.split(";"))
    m0 = gto.M(atom=spec, basis="sto-3g", charge=-1, verbose=0)
    m1 = gto.M(atom=shifted, basis="sto-3g", charge=-1, verbose=0)
    dm = _fake_dm(m0)
    d0, d1 = em.dipole_about_nuclear_centre(m0, dm), em.dipole_about_nuclear_centre(m1, dm)
    assert np.allclose(d0["vector_D"], d1["vector_D"], atol=1e-6)
    assert np.allclose(np.array(d1["origin"]["xyz_angstrom"]) -
                       np.array(d0["origin"]["xyz_angstrom"]), 5.0, atol=1e-8)


def test_neutral_dipole_matches_pyscf():
    """For a neutral molecule the origin does not matter, so the value must
    equal PySCF's own dip_moment on the same density."""
    from pyscf import gto, scf
    m = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692",
              basis="sto-3g", verbose=0)
    mf = scf.RHF(m)
    dm = _fake_dm(m)
    dm *= m.nelectron / np.trace(dm @ m.intor("int1e_ovlp"))      # neutral
    ref = mf.dip_moment(mol=m, dm=dm, unit="Debye", verbose=0)
    got = em.dipole_about_nuclear_centre(m, dm)
    assert np.allclose(got["vector_D"], ref, atol=1e-5)


def test_sigma_hole_record_names_the_axis_atom():
    """The JSON's i_bond is the atom sigma_hole actually used, and a missing
    sigma-hole is null with a reason, never silently absent."""
    from pyscf import gto
    spec = em.monosubstituted_benzene(*em.MOLECULES["chlorobenzene"])
    mol = gto.M(atom=spec, basis="sto-3g", verbose=0)
    co = mol.atom_coords()
    axis = co[6] - co[0]
    verts = co[6] + 3.0 * axis / np.linalg.norm(axis) + np.zeros((1, 3))
    rec, why = em.sigma_hole_record(verts, np.array([7.5]), mol, 6)
    assert why is None and rec == {"V_max_kcal_mol": 7.5, "n_cap_vertices": 1,
                                   "i_sub": 6, "i_bond": 0}
    assert em.sigma_hole_record(verts, np.array([7.5]), mol, None,
                                "fluorine does not form one") == (
        None, "fluorine does not form one")


def test_scale_provenance_agrees_with_choose_scale():
    """The flags scale_provenance trusts are the ones choose_scale obeys."""
    anion = np.array([-170.0, -120.0, -80.0])
    vlim, kind, _ = em.choose_scale(anion, vmin=-200.0, vmax=0.0)
    assert kind == "sequential" and vlim == (-200.0, 0.0)
    vlim, kind, _ = em.choose_scale(anion, lim=50.0)          # ignored
    assert kind == "sequential" and vlim == (-170.0, -80.0)
    neutral = np.array([-40.0, 10.0, 55.0])
    vlim, kind, _ = em.choose_scale(neutral, lim=60.0)
    assert kind == "diverging" and vlim == (-60.0, 60.0)


def _rgb(s):
    return np.array([int(x) for x in s[4:-1].split(",")]) / 255.0


def test_3d_view_colours_run_the_same_way_as_the_image():
    """The interactive view must not re-derive colour orientation: the same
    trap as the three maps that shipped inverted (review 1)."""
    for name, cmap in sorted(em.build_maps().items()):
        cs = em.plotly_colorscale(cmap)
        ts = [t for t, _ in cs]
        assert ts[0] == 0.0 and ts[-1] == 1.0 and ts == sorted(ts), name
        lo, hi = _rgb(cs[0][1]), _rgb(cs[-1][1])
        assert np.allclose(lo, cmap(0.0)[:3], atol=1 / 255), name
        assert np.allclose(hi, cmap(1.0)[:3], atol=1 / 255), name
        if name not in em.SEQUENTIAL:
            assert lo[0] > lo[2] and hi[2] > hi[0], f"{name}: warm must be negative"


def test_3d_view_maps_vlim_and_has_no_specular_highlight():
    pytest.importorskip("plotly")
    verts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    V = np.array([-30.0, -5.0, 10.0, 40.0])
    fig = em.esp_figure(verts, faces, V, ["O", "H"], np.array([[0, 0, 0], [0, 0, 1.8]]),
                        em.build_maps()["vik"], (-40.0, 40.0))
    surf = fig.data[0]
    assert (surf.cmin, surf.cmax) == (-40.0, 40.0)    # V_min side gets t=0, warm
    assert list(surf.intensity) == list(V)
    assert surf.lighting.specular == 0.0
    assert len(fig.data[1].x) == 3                     # one O-H bond: 2 points + gap
    assert "warm = negative" in surf.colorbar.title.text


def test_3d_colourbar_does_not_say_warm_is_negative_on_a_sequential_scale():
    """Review 3, H3: the HTML colourbar said 'warm = negative' on an ion too."""
    pytest.importorskip("plotly")
    verts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    V = np.array([-170.0, -150.0, -100.0, -80.0])
    fig = em.esp_figure(verts, faces, V, ["O"], np.zeros((1, 3)),
                        em.build_maps()["lajolla"], (-170.0, -80.0), kind="sequential")
    title = fig.data[0].colorbar.title.text
    assert "warm" not in title and "more negative" in title


# ------------------------------------------------------------------- sites
# per-site-spec P2: the triangle partition and the area-weighted statistics.
# Meshes come from marching cubes on a model density (atom-centred
# exponentials), the same construction extract_surface uses, so the
# triangles are as uneven as real ones.
def _blob_mesh(coords_A, scale=None, iso=0.05, h=0.12):
    """(verts_A, faces) of the iso-surface of sum_i exp(-2 |r - R_i| / s_i)."""
    from skimage import measure
    X = np.asarray(coords_A, float)
    s = np.ones(len(X)) if scale is None else np.asarray(scale, float)
    lo = X.min(0) - 3.5
    n = np.ceil((X.max(0) + 3.5 - lo) / h).astype(int)
    g = np.stack(np.meshgrid(*[lo[k] + h * np.arange(n[k]) for k in range(3)],
                             indexing="ij"), -1)
    rho = sum(np.exp(-2.0 * np.linalg.norm(g - X[i], axis=-1) / s[i]) for i in range(len(X)))
    v, f, _, _ = measure.marching_cubes(rho, level=iso)
    return lo + v * h, f.astype(int)


_CH2O = (["C", "O", "H", "H"],
         [[0.0, 0.0, 0.0], [1.21, 0.0, 0.0], [-0.55, 0.94, 0.0], [-0.55, -0.94, 0.0]])


def _field(P):
    """A smooth signed model potential: negative near +x (the 'O' end)."""
    return -30.0 * P[:, 0] + 5.0 * P[:, 1] ** 2


def _three_triangles():
    """Areas 1, 3, 1 with V = -4, 8, 4 on every vertex."""
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0],
                  [5, 0, 0], [8, 0, 0], [5, 2, 0],
                  [20, 0, 0], [21, 0, 0], [20, 2, 0]], float)
    faces = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
    V = np.repeat([-4.0, 8.0, 4.0], 3)
    return P, faces, V


def test_gipf_variances_are_taken_about_each_signs_own_mean():
    """Review 4 N11: V_area_sd is one SD about the overall mean, cited as
    GIPF. GIPF's sigma^2_tot is sigma_+^2 + sigma_-^2, each about its own
    sign's mean, with nu = sigma_+^2 sigma_-^2 / (sigma^2_tot)^2."""
    P, faces, V = _three_triangles()
    g = em.gipf_stats(P, faces, V)
    assert g["V_pos_mean_kcal_mol"] == pytest.approx((3 * 8 + 1 * 4) / 4)     # 7
    assert g["V_neg_mean_kcal_mol"] == pytest.approx(-4.0)
    assert g["sigma2_pos_kcal2_mol2"] == pytest.approx((3 * 1 + 1 * 9) / 4)  # 3
    assert g["sigma2_neg_kcal2_mol2"] == pytest.approx(0.0)
    assert g["sigma2_tot_kcal2_mol2"] == pytest.approx(3.0)
    assert g["nu"] == pytest.approx(0.0)
    sd = em.area_stats(P, faces, V)["V_area_sd_kcal_mol"]
    assert sd ** 2 != pytest.approx(g["sigma2_tot_kcal2_mol2"])  # not the same quantity


def test_gipf_balance_is_a_quarter_when_both_signs_vary_equally():
    P, faces, _ = _three_triangles()
    P = np.vstack([P, P[3:6] + [0, 10, 0]])
    faces = np.vstack([faces, [[9, 10, 11]]])        # areas 1, 3, 1, 3
    V = np.repeat([-4.0, 8.0, 4.0, -8.0], 3)        # mirror images in sign
    g = em.gipf_stats(P, faces, V)
    assert g["sigma2_pos_kcal2_mol2"] == pytest.approx(g["sigma2_neg_kcal2_mol2"])
    assert g["nu"] == pytest.approx(0.25)


def test_gipf_on_a_one_signed_surface_leaves_the_missing_sign_null():
    """An anion's surface has no positive part: its mean and variance are
    undefined, not zero, and nu with them."""
    P, faces, V = _three_triangles()
    g = em.gipf_stats(P, faces, V - 20.0)
    assert g["V_pos_mean_kcal_mol"] is None and g["sigma2_pos_kcal2_mol2"] is None
    assert g["nu"] is None
    assert g["sigma2_tot_kcal2_mol2"] == pytest.approx(g["sigma2_neg_kcal2_mol2"])


def test_area_stats_are_area_weighted_not_vertex_averaged():
    """Two triangles of areas 1 and 3: a vertex average would weight them by
    vertex count. The spec forbids it (marching-cubes vertex density is a
    mesh artefact, not chemistry)."""
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0],           # area 1
                  [5, 0, 0], [8, 0, 0], [5, 2, 0]], float)   # area 3
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    V = np.array([-4.0, -4.0, -4.0, 8.0, 8.0, 8.0])
    s = em.area_stats(P, faces, V)
    assert s["area_A2"] == pytest.approx(4.0)
    assert s["area_neg_A2"] == pytest.approx(1.0) and s["area_pos_A2"] == pytest.approx(3.0)
    assert s["V_area_mean_kcal_mol"] == pytest.approx((1 * -4 + 3 * 8) / 4)
    assert s["V_min_kcal_mol"] == -4.0 and s["V_max_kcal_mol"] == 8.0
    assert s["V_p05_kcal_mol"] == -4.0 and s["V_p95_kcal_mol"] == 8.0


@pytest.mark.parametrize("rule", ["nearest_nucleus", "radius_scaled"])
def test_site_areas_conserve_the_surface_exactly(rule):
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym)
    total = em.area_stats(P, faces, _field(P))["area_A2"]
    atoms = sum(a[rule]["area_A2"] for a in rep["atoms"])
    sites = sum(s[rule]["area_A2"] for s in rep["reaction_sites"])
    assert atoms == pytest.approx(total, rel=1e-12)
    assert sites == pytest.approx(total, rel=1e-12)


def test_hydrogens_fold_into_their_parent_site():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym)
    by_index = {a["index"]: a for a in rep["atoms"]}
    c_site = next(s for s in rep["reaction_sites"] if s["atoms"] == [0, 2, 3])
    expect = sum(by_index[i]["nearest_nucleus"]["area_A2"] for i in (0, 2, 3))
    assert c_site["nearest_nucleus"]["area_A2"] == pytest.approx(expect)
    assert [s["label"] for s in rep["reaction_sites"]] == ["C0", "O1"]
    assert by_index[2]["parent"] == 0 and by_index[0]["parent"] is None


def test_site_extrema_are_the_global_extrema():
    """The same vertices, so exactly equal: a partition that loses triangles,
    or reads V from the wrong index, breaks this."""
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    V = _field(P)
    rep = em.site_report(P, faces, V, X, sym)
    for rule in em.PARTITION_RULES:
        mins = [s[rule]["V_min_kcal_mol"] for s in rep["reaction_sites"]]
        maxs = [s[rule]["V_max_kcal_mol"] for s in rep["reaction_sites"]]
        assert min(mins) == V[np.unique(faces)].min()
        assert max(maxs) == V[np.unique(faces)].max()


def test_a_buried_atom_is_reported_with_zero_area():
    """'Not on the surface' must stay distinguishable from 'not computed'."""
    # an atom at the centre of an octahedron: any surface point more than
    # ~1.2 A out is nearer one of the six than the centre, by either rule
    sym = ["C"] + ["O"] * 6
    X = [[0.0, 0.0, 0.0]] + [[s * 1.4 if k == j else 0.0 for k in range(3)]
                             for j in range(3) for s in (1, -1)]
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym)
    c_atom = rep["atoms"][0]
    c_site = next(s for s in rep["reaction_sites"] if s["label"] == "C0")
    for entry in (c_atom, c_site):
        assert entry["nearest_nucleus"]["area_A2"] == 0.0
        assert entry["nearest_nucleus"]["V_area_mean_kcal_mol"] is None
        assert entry["delta"]["V_area_mean_kcal_mol"] is None
    assert len(rep["atoms"]) == 7


def test_radius_scaled_rule_can_differ_from_nearest_nucleus():
    """A centroid 1.5 A from H and 1.6 A from C belongs to H by distance but
    to C by d / r_Bondi (1.25 vs 0.94)."""
    sym = ["C", "H"]
    X = np.array([[0.0, 0.0, 0.0], [3.1, 0.0, 0.0]])
    tri = np.array([[1.6, 0.0, 0.0], [1.6, 0.01, 0.0], [1.6, 0.0, 0.01]])
    owner = {r: em.partition_surface(tri, np.array([[0, 1, 2]]), X, sym, r)[0]
             for r in em.PARTITION_RULES}
    assert owner == {"nearest_nucleus": 1, "radius_scaled": 0}


def test_radius_scaled_is_null_for_an_element_without_a_radius():
    sym = ["Fe", "O"]
    X = [[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]]
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym)
    assert rep["partition"]["compared"] is None and rep["partition"]["radii"] is None
    assert all(a["radius_scaled"] is None and a["delta"] is None for a in rep["atoms"])
    assert rep["atoms"][0]["nearest_nucleus"]["area_A2"] > 0


def test_site_block_is_withheld_from_an_invalid_record():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    for bad in (dict(input=False, surface=True), dict(input=True, surface=False)):
        block, reason = em.sites_block(dict(bad, sites=False), P, faces, _field(P), X, sym)
        assert block is None and reason
    block, reason = em.sites_block(dict(input=True, surface=True, sites=False),
                                   P, faces, _field(P), X, sym)
    assert block["schema"] == em.SITES_SCHEMA and reason is None


# ------------------------------------------------- review 5 R8: say which check
@pytest.mark.parametrize("gate, verified, mapped, says", [
    (dict(forced=True), True, True, "gate"),
    (None, None, None, "no structure gate ran"),
    (dict(forced=False), False, True, "structure not verified"),
    (dict(forced=False), True, False, "atom order not verified"),
    (dict(forced=True), False, False, "gate"),          # every cause is named
])
def test_input_reason_names_the_failing_check(gate, verified, mapped, says):
    """R8: a forced atom-order failure read 'forced past the gate, or
    structure not verified', neither of which had happened."""
    reason = em.input_invalid_reason(gate, verified, mapped)
    assert says in reason
    if gate and gate["forced"] and verified is False:
        assert "structure not verified" in reason and "atom order" in reason


def test_sites_block_passes_the_input_reason_through():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    block, reason = em.sites_block(dict(input=False, surface=True, sites=False),
                                   P, faces, _field(P), X, sym,
                                   input_reason="input not valid: atom order not verified")
    assert block is None and "atom order not verified" in reason


def test_a_failed_identity_check_records_what_was_perceived():
    """R8: verify_computed returned None for the computed SMILES on a
    connectivity mismatch, so a forced record could not say what it ran."""
    from rdkit import Chem
    spec, _ = em.from_smiles("COC")
    ok, computed, why = em.verify_computed(spec, Chem.MolFromSmiles("CCO"))
    assert ok is False and computed == "COC" and "connectivity" in why


# ------------------------------------------------ review 5 R9: NPZ partition
def test_npz_partition_follows_the_site_table():
    """R9: the NPZ wrote tri_atom_* whenever input and surface were valid, so
    H2 (no heavy atom, no site table) still shipped a partition."""
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    assert em.tri_atom_arrays(False, P, faces, X, sym) == {}
    arrays = em.tri_atom_arrays(True, P, faces, X, sym)
    assert set(arrays) == {"tri_atom_nearest_nucleus", "tri_atom_radius_scaled"}
    assert all(len(v) == len(faces) and v.min() >= 0 for v in arrays.values())


def test_npz_omits_a_rule_that_cannot_run():
    """R9: absence had three encodings; an unavailable rule was an array of
    -1. It is now omitted, as the schema's npz list already allows."""
    sym = ["Fe", "O"]
    X = [[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]]
    P, faces = _blob_mesh(X)
    assert set(em.tri_atom_arrays(True, P, faces, X, sym)) == {"tri_atom_nearest_nucleus"}


# ------------------------------------------- review 4 N5: one record contract
def test_partition_names_its_primary_and_compared_rules():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym)
    assert rep["partition"] == {"primary": "nearest_nucleus", "compared": "radius_scaled",
                                "unit": "triangle", "radii": "Bondi; Mantina 2009 for B"}


def test_a_site_block_that_fails_its_own_checks_is_withheld(monkeypatch):
    """validity.sites was just 'a block exists'; conservation and nullability
    were asserted only in tests. A triangle lost from one atom must now
    withhold the block with a reason."""
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    good = em.site_report(P, faces, _field(P), X, sym)
    assert em.site_checks(good, P, faces, _field(P)) == []
    real = em.site_report

    def lossy(*a, **k):
        rep = real(*a, **k)
        rep["atoms"][1]["nearest_nucleus"]["area_A2"] *= 0.9
        return rep
    monkeypatch.setattr(em, "site_report", lossy)
    block, reason = em.sites_block(dict(input=True, surface=True, sites=False),
                                   P, faces, _field(P), X, sym)
    assert block is None and "conserve" in reason


def test_site_checks_catch_a_value_on_a_zero_area_atom():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym)
    rep["atoms"][0]["nearest_nucleus"]["area_A2"] = 0.0
    assert any("null" in p for p in em.site_checks(rep, P, faces, _field(P)))


def test_geometry_provenance_counts_every_hydrogen_source():
    """Only the keys that occurred were written ({'generated_relaxed': 6}),
    so a consumer had to guess that an absent key meant zero."""
    ctx = em.site_context(solvent_esp=None, heavy_atoms="generated",
                          h_source=[None, "generated_relaxed", "generated_relaxed"],
                          relaxation="ETKDGv3 + MMFF94s, all atoms, isolated molecule",
                          n_rot=0, surface=dict(genus=0, components=1),
                          identity={}, budget=[])
    gp = ctx["geometry_provenance"]
    assert gp["hydrogens"] == {"supplied": 0, "generated": 0, "generated_relaxed": 2,
                               "optimised": 0}
    assert gp["relaxation"].startswith("ETKDGv3") and gp["caveat"]
    assert ctx["conformer_scope"].startswith("single conformer; 0 rotatable")


def _budget_ids(**kw):
    return [r["id"] for r in em.budget_rows(**kw)]


def test_budget_rows_follow_the_run_settings():
    gas = _budget_ids(basis="def2-SVP", solvent=None, charge=0, n_rot=0, geometry_method="etkdg_mmff94s")
    assert {"basis", "solvent water", "solvent chloroform", "partition C18",
            "partition small anions", "grid mean"} <= set(gas)
    assert "conformer C18" not in gas and "anion geometry" not in gas
    wet = _budget_ids(basis="def2-tzvpp", solvent="water", charge=-1, n_rot=3, geometry_method="etkdg_mmff94s")
    assert "basis" not in wet and not any(i.startswith("solvent") for i in wet)
    # N9: any spelling PySCF resolves to def2-TZVPP drops the basis row too
    assert "basis" not in _budget_ids(basis="DEF2_TZVPP", solvent=None, charge=0, n_rot=0, geometry_method="etkdg_mmff94s")
    assert {"conformer C18", "anion geometry"} <= set(wet)
    # --atoms: the rotatable-bond count is unknown, so the conformer row stays
    assert "conformer C18" in _budget_ids(basis="sto-3g", solvent=None, charge=0, n_rot=None, geometry_method="undeclared")


def test_every_budget_row_names_a_bench_source():
    for r in em.budget_rows(basis="def2-svp", solvent=None, charge=-1, n_rot=None, geometry_method="etkdg_mmff94s"):
        assert set(r) == {"id", "axis", "statistic", "kcal_mol", "measured_on", "source"}
        assert r["source"].startswith("bench/") and isinstance(r["kcal_mol"], float)


def test_the_site_block_matches_the_committed_schema():
    """The sites subtree, both directions, without an SCF; the slow suite
    checks the whole record from main()."""
    from record_tree import load_schema, tree_problems
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    hs = [None, None, "supplied", "supplied"]
    block, _ = em.sites_block(dict(input=True, surface=True, sites=False),
                              P, faces, _field(P), X, sym, atom_maps={0: 1, 1: 2}, h_source=hs)
    block.update(em.site_context(
        solvent_esp=None, heavy_atoms="supplied_frozen", h_source=hs,
        relaxation="none", n_rot=None, surface=dict(genus=0, components=1),
        identity=dict(chemical_smiles="C=O", mapped_smiles="[CH2:1]=[O:2]",
                      chemical_match=True, mapped_match=True),
        budget=em.budget_rows(basis="sto-3g", solvent=None, charge=0, n_rot=None, geometry_method="undeclared")))
    assert tree_problems(load_schema()["record"]["sites"], block, "sites") == []


def test_the_schema_comparison_runs_in_both_directions():
    from record_tree import tree_problems
    schema = {"a": "int", "b": {"$nullable": True, "c": ["str?"]}}
    assert tree_problems(schema, {"a": 1, "b": None}) == []
    assert tree_problems(schema, {"a": 1, "b": {"c": ["x", None]}}) == []
    assert tree_problems(schema, {"b": None}) == [
        "record.a: declared by the schema, absent from the record"]
    assert tree_problems(schema, {"a": 1, "b": None, "z": 0}) == [
        "record.z: written by the code, absent from the schema"]
    assert tree_problems(schema, {"a": True, "b": None}) == ["record.a: bool, expected int"]


def test_record_schema_names_its_version():
    from record_tree import load_schema
    s = load_schema()
    assert s["schema"] == em.RECORD_SCHEMA == "espmap.record/1"
    assert s["record"]["schema"] == "str" and s["record"]["sites"]["schema"] == "str"


# ----------------------------- review 5: the contract before the P3 freeze
def test_default_site_labels_use_the_atom_map_on_a_mapped_input():
    """R5: nitrite's 'O1' was atom index 1 = map 2, beside C18 'C3-H' = map 3.
    A mapped input now labels default sites by map, an unmapped one by
    index, and source says which."""
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    mapped = em.site_report(P, faces, _field(P), X, sym, atom_maps={0: 1, 1: 2})
    assert [s["label"] for s in mapped["reaction_sites"]] == ["C1", "O2"]
    assert {s["source"] for s in mapped["reaction_sites"]} == {"default:map"}
    plain = em.site_report(P, faces, _field(P), X, sym)
    assert [s["label"] for s in plain["reaction_sites"]] == ["C0", "O1"]
    assert {s["source"] for s in plain["reaction_sites"]} == {"default:index"}


def test_overlapping_sites_say_so():
    """Review 5: the C18 'acetyl' and 'carbonyl O' groups share map 32, so
    the site areas summed to 1.041 of the surface with nothing saying so."""
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym,
                         groups=[("carbonyl", [0, 1, 2, 3]), ("O", [1]), ("CH2", [0, 2, 3])])
    ov = {s["label"]: s["overlaps"] for s in rep["reaction_sites"]}
    assert ov == {"carbonyl": ["O", "CH2"], "O": ["carbonyl"], "CH2": ["carbonyl"]}
    assert rep["sites_tile_surface"] is False
    default = em.site_report(P, faces, _field(P), X, sym)
    assert default["sites_tile_surface"] is True
    assert all(s["overlaps"] == [] for s in default["reaction_sites"])


def test_caller_sites_that_miss_atoms_do_not_tile_the_surface():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym, groups=[("O", [1])])
    assert rep["sites_tile_surface"] is False


def test_gipf_needs_both_signs_above_a_minimum_area():
    """R14: one triangle at +0.01 on an all-negative surface turned nu from
    null into 0.0. The same knife-edge review 3 H11 removed from the scale."""
    P, faces, V = _three_triangles()
    P[6:9] = P[6] + 0.1 * (P[6:9] - P[6])          # third triangle: area 1 -> 0.01
    V = V - 20.0                                   # all negative
    V[6:9] = 0.01                                  # one small positive triangle
    a = em.triangle_areas(P, faces)
    assert a[2] / a.sum() < em.GIPF_MIN_AREA_FRACTION
    g = em.gipf_stats(P, faces, V)
    assert g["nu"] is None and g["sigma2_pos_kcal2_mol2"] is None
    assert "positive" in g["reason"] and g["min_area_fraction"] == em.GIPF_MIN_AREA_FRACTION
    both = em.gipf_stats(*_three_triangles())
    assert both["reason"] is None


def test_record_bonds_are_the_approved_structure_in_record_order():
    """A1: a consumer could rebuild bonds only from a mapped SMILES. The
    record now carries them: one Kekule structure, hydrogens included."""
    from rdkit import Chem
    probe = Chem.MolFromSmiles("c1ccccc1O")
    bonds = em.record_bonds(probe)
    heavy = [b for b in bonds if b[0] < 7 and b[1] < 7]
    assert len(bonds) == 13 and len(heavy) == 7
    assert sorted(b[2] for b in heavy if b[0] < 6 and b[1] < 6) == [1, 1, 1, 2, 2, 2]
    assert all(b[0] < b[1] for b in bonds)
    h = Chem.AddHs(probe)
    assert {(b[0], b[1]) for b in bonds} == {
        tuple(sorted((x.GetBeginAtomIdx(), x.GetEndAtomIdx()))) for x in h.GetBonds()}


def test_schema_vocabularies_are_enforced():
    """R10: spec and tests used verdict 'ACCEPT', which the code never
    writes. Closed vocabularies are now enums the comparison checks."""
    from record_tree import load_schema, tree_problems
    gate = load_schema()["record"]["input"]["gate"]
    ok = {"verdict": "INFO", "forced": False, "findings": []}
    assert tree_problems(gate, ok) == []
    assert tree_problems(gate, dict(ok, verdict="ACCEPT")) == [
        "record.verdict: 'ACCEPT', expected one of INFO|WARN|BLOCK"]
    assert tree_problems({"x": "enum:a|b?"}, {"x": None}) == []


def test_open_surface_genus_is_allowed_to_be_null():
    """R4: mesh_topology writes genus None on an open surface; the schema
    said int, so a --no-strict record failed its own contract."""
    from record_tree import load_schema, tree_problems
    surf = load_schema()["record"]["results"]["surface"]
    assert tree_problems(surf["genus"], None) == []


def test_h_source_records_which_hydrogens_the_file_supplied(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(m, randomSeed=3)
    rw = Chem.RWMol(m)
    for i in sorted((a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "H"
                     and a.GetNeighbors()[0].GetSymbol() == "C"), reverse=True):
        rw.RemoveAtom(i)
    p = tmp_path / "pose.mol"
    p.write_text(Chem.MolToMolBlock(rw.GetMol()))
    _, info = em.from_file(str(p))
    hs = [s for s in info["h_source"] if s is not None]
    assert hs.count("supplied") == 1 and hs.count("generated_relaxed") == 5
    _, info = em.from_smiles("CCO")
    assert {s for s in info["h_source"] if s is not None} == {"generated_relaxed"}


# ---------------------------------------------------------- atom identity
# per-site-spec P1. Maps are metadata: they never change what is computed,
# but a site table keyed by them must point at the right atoms.
def _map_codes(smiles=None, path=None):
    from rdkit import Chem
    mol, _ = em.parse_for_check(smiles=smiles, path=path)
    _, findings = em.check_structure(Chem.RemoveHs(mol), maps_from=mol, smiles=smiles)
    return {code for _, code, *_ in findings}


@pytest.mark.parametrize("smiles, code", [
    ("[CH3:1]CO", "PARTIAL_ATOM_MAP"),               # a table with holes
    ("[CH3:1][CH2:1][OH:2]", "DUPLICATE_ATOM_MAP"),
    ("[CH3:1][CH2:2][O:3][H:4]", "MAPPED_HYDROGEN"),  # RDKit's parser drops this map silently
])
def test_ambiguous_atom_maps_are_blocked(smiles, code):
    assert code in _map_codes(smiles)
    from rdkit import Chem
    mol, _ = em.parse_for_check(smiles=smiles)
    verdict, _ = em.check_structure(Chem.RemoveHs(mol), maps_from=mol, smiles=smiles)
    assert verdict == em.BLOCK


@pytest.mark.parametrize("smiles", ["[CH3:1][CH2:2][OH:3]", "CCO"])
def test_full_or_absent_maps_pass(smiles):
    assert not _map_codes(smiles) & {"PARTIAL_ATOM_MAP", "DUPLICATE_ATOM_MAP",
                                     "MAPPED_HYDROGEN"}


def test_a_mapped_hydrogen_in_a_file_is_blocked(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("[CH3:1][CH2:2][OH:3]"))
    AllChem.EmbedMolecule(m, randomSeed=1)
    m.GetAtomWithIdx(m.GetNumAtoms() - 1).SetAtomMapNum(9)
    p = tmp_path / "h.sdf"
    w = Chem.SDWriter(str(p)); w.write(m); w.close()
    assert "MAPPED_HYDROGEN" in _map_codes(path=str(p))


@pytest.mark.parametrize("smiles", ["[CH3:1][CH:2]([CH3:3])[OH:4]",
                                    "[CH3:1][C:2]([CH3:3])([CH3:4])[OH:5]"])
def test_a_mapped_molecule_with_equivalent_groups_is_not_refused(smiles):
    """Maps make isopropyl's two methyls distinct labels, and RDKit then sees
    a stereocentre the chemistry does not have: verify_computed refused a
    correctly mapped isopropanol as 'stereochemistry differs'. Chemical
    identity is compared with maps stripped; verify_mapping owns the maps."""
    from rdkit import Chem
    spec, _ = em.from_smiles(smiles)
    ok, computed, reason = em.verify_computed(spec, Chem.MolFromSmiles(smiles))
    assert ok, reason
    assert ":" not in computed                        # chemical identity, no maps


def test_atom_order_must_match_the_approval_by_index():
    """verify_computed searches for SOME match, so it accepts a geometry whose
    atoms were reordered. A site table keyed by index cannot: spec atom i
    must be approval atom i."""
    from rdkit import Chem
    probe = Chem.MolFromSmiles("[CH3:1][CH:2]([CH3:3])[OH:4]")
    spec, _ = em.from_smiles(Chem.MolToSmiles(probe))
    assert em.verify_mapping(spec, probe)[0]
    rows = spec.split(";")
    rows[0], rows[2] = rows[2], rows[0]               # swap the two methyl carbons
    swapped = ";".join(rows)
    assert em.verify_computed(swapped, Chem.RemoveHs(probe))[0]   # same molecule
    ok, reason = em.verify_mapping(swapped, probe)
    assert not ok and reason


def test_supplied_coordinates_must_stay_where_the_file_put_them(tmp_path):
    """A pose: map 1 must sit where the adapter wrote atom 1. Swap the two
    methyl GROUPS consistently: chemistry and connectivity are unchanged, so
    only the comparison with the file catches it."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("[CH3:1][CH:2]([CH3:3])[OH:4]"))
    AllChem.EmbedMolecule(m, randomSeed=2)
    p = tmp_path / "pose.sdf"
    w = Chem.SDWriter(str(p)); w.write(m); w.close()
    probe, _ = em.parse_for_check(path=str(p))
    spec, info = em.from_file(str(p))
    assert em.verify_mapping(spec, probe, reference=True)[0]
    sym, X = zip(*[(t.split()[0], [float(v) for v in t.split()[1:]]) for t in spec.split(";")])
    X = np.array(X)
    nb = lambda i: [n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1]
    perm = list(range(len(sym)))
    for a, b in [(0, 2)] + list(zip(nb(0), nb(2))):
        perm[a], perm[b] = b, a
    swapped = "; ".join(f"{sym[i]} {X[perm[i]][0]:.6f} {X[perm[i]][1]:.6f} {X[perm[i]][2]:.6f}"
                        for i in range(len(sym)))
    assert em.verify_mapping(swapped, probe)[0]       # internally consistent
    ok, reason = em.verify_mapping(swapped, probe, reference=True)
    assert not ok and "file" in reason


def test_main_refuses_coordinates_whose_atom_order_moved(monkeypatch):
    from rdkit import Chem
    real = em.from_smiles

    def reordered(s, **k):
        spec, info = real(s, **k)
        rows = spec.split(";")
        rows[0], rows[2] = rows[2], rows[0]
        return ";".join(rows), info
    monkeypatch.setattr(em, "from_smiles", reordered)
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[CH3:1][CH:2]([CH3:3])[OH:4]"])
    with pytest.raises(SystemExit, match="atom order"):
        em.main()


# ------------------------------------------------------------ caller sites
# per-site-spec section 1: Arrowfield names its own sites ("C3-H",
# "carbonyl O"); espmap only checks that every named atom exists.
def _ethanol_mapped():
    from rdkit import Chem
    spec, info = em.from_smiles("[CH3:1][CH2:2][OH:3]")
    sym, X = zip(*[(t.split()[0], [float(v) for v in t.split()[1:]]) for t in spec.split(";")])
    maps = {a.GetIdx(): a.GetAtomMapNum() for a in info["mol"].GetAtoms() if a.GetAtomMapNum()}
    return list(sym), np.array(X), maps


def _ethanol_plain():
    spec, _ = em.from_smiles("CCO")
    sym, X = zip(*[(t.split()[0], [float(v) for v in t.split()[1:]]) for t in spec.split(";")])
    return list(sym), np.array(X)


def test_caller_groups_resolve_maps_and_hydrogens():
    sym, X, maps = _ethanol_mapped()
    groups = [{"label": "C1-H", "atoms": ["1+H"]}, {"label": "hydroxyl", "atoms": ["3+H"]},
              {"label": "C2", "atoms": [2]}]
    by, res = em.resolve_site_groups({"by": "map", "groups": groups}, sym, X, maps)
    got = dict(res)
    parent = em.site_parents(X, sym)
    assert by == "map"
    assert got["C1-H"] == sorted([0] + [i for i, p in enumerate(parent) if p == 0])
    assert len(got["C1-H"]) == 4 and got["C2"] == [1]
    assert got["hydroxyl"] == sorted([2] + [i for i, p in enumerate(parent) if p == 2])


def test_unmapped_input_names_atoms_by_index():
    sym, X = _ethanol_plain()
    by, res = em.resolve_site_groups({"by": "index", "groups": [{"label": "O", "atoms": [2]}]},
                                     sym, X, None)
    assert by == "index" and dict(res) == {"O": [2]}


@pytest.mark.parametrize("groups, words", [
    ([{"label": "x", "atoms": [7]}], "7"),                         # no such map
    ([{"label": "x", "atoms": []}], "no atoms"),
    ([{"label": "a", "atoms": [1]}, {"label": "a", "atoms": [2]}], "twice"),
    ([{"label": "x", "atoms": ["1+X"]}], r"1\+X"),                 # malformed token
    ([{"atoms": [1]}], "label"),
])
def test_bad_caller_groups_are_refused_with_a_reason(groups, words):
    sym, X, maps = _ethanol_mapped()
    with pytest.raises(ValueError, match=words):
        em.resolve_site_groups({"by": "map", "groups": groups}, sym, X, maps)


# Review 4 N6: a name meant map numbers on mapped input and atom indices on
# unmapped input, silently. Maps are 1-based and indices 0-based, so a lost
# map always landed on a plausible NEIGHBOUR.
@pytest.mark.parametrize("spec, words", [
    ([{"label": "C1-H", "atoms": ["1+H"]}], "by"),                 # how atoms are counted, unsaid
    ({"groups": [{"label": "C1-H", "atoms": ["1+H"]}]}, "by"),
    ({"by": "name", "groups": []}, "by"),
])
def test_a_sites_file_must_say_how_it_counts_atoms(spec, words):
    sym, X, maps = _ethanol_mapped()
    with pytest.raises(ValueError, match=words):
        em.resolve_site_groups(spec, sym, X, maps)


def test_map_names_on_an_unmapped_input_are_refused():
    """['1+H'] written for map 1 resolved to C index 1 when the maps were lost."""
    sym, X = _ethanol_plain()
    with pytest.raises(ValueError, match="no atom maps"):
        em.resolve_site_groups({"by": "map", "groups": [{"label": "C1-H", "atoms": ["1+H"]}]},
                               sym, X, None)


@pytest.mark.parametrize("atoms", [[5], ["5+H"]])
def test_an_index_naming_a_hydrogen_is_refused(atoms):
    sym, X = _ethanol_plain()
    assert sym[5] == "H"
    with pytest.raises(ValueError, match="hydrogen"):
        em.resolve_site_groups({"by": "index", "groups": [{"label": "h", "atoms": atoms}]},
                               sym, X, None)


def test_an_index_is_the_file_order_index_hydrogens_included(tmp_path):
    """An SDF that interleaves its hydrogens: O is index 7, not heavy atom 2."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(m, randomSeed=1)
    order = [0, 3, 4, 5, 1, 6, 7, 2, 8]                   # C H H H C H H O H
    m2 = Chem.RenumberAtoms(m, order)
    p = tmp_path / "interleaved.sdf"
    w = Chem.SDWriter(str(p)); w.write(m2); w.close()
    spec, info = em.from_file(str(p))
    sym = [t.split()[0] for t in spec.split(";")]
    X = np.array([[float(v) for v in t.split()[1:4]] for t in spec.split(";")])
    _, res = em.resolve_site_groups({"by": "index", "groups": [{"label": "O", "atoms": [7]}]},
                                    sym, X, None)
    assert sym[7] == "O" and dict(res) == {"O": [7]}


def test_sites_record_how_they_were_named():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    rep = em.site_report(P, faces, _field(P), X, sym, groups=[("O", [1])], source="caller:index")
    assert rep["reaction_sites"][0]["source"] == "caller:index"


def test_overlapping_caller_groups_each_get_their_own_values():
    sym, X = _CH2O
    P, faces = _blob_mesh(X)
    groups = [("CO", [0, 1]), ("O", [1]), ("C-H", [0, 2, 3])]
    rep = em.site_report(P, faces, _field(P), X, sym, groups=groups)
    s = {g["label"]: g for g in rep["reaction_sites"]}
    assert set(s) == {"CO", "O", "C-H"} and all(g["source"] == "caller" for g in s.values())
    by_atom = {a["index"]: a["nearest_nucleus"]["area_A2"] for a in rep["atoms"]}
    assert s["CO"]["nearest_nucleus"]["area_A2"] == pytest.approx(by_atom[0] + by_atom[1])
    assert s["O"]["nearest_nucleus"]["area_A2"] == pytest.approx(by_atom[1])


def test_main_refuses_bad_groups_before_the_scf(tmp_path, monkeypatch):
    import json
    p = tmp_path / "g.json"
    p.write_text(json.dumps({"by": "map", "groups": [{"label": "C9-H", "atoms": ["9+H"]}]}))
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[CH3:1][CH2:2][OH:3]",
                                      "--sites", str(p)])
    with pytest.raises(SystemExit, match="--sites"):
        em.main()


def test_main_lets_good_groups_through_to_the_scf(tmp_path, monkeypatch):
    import json
    p = tmp_path / "g.json"
    p.write_text(json.dumps({"by": "map", "groups": [{"label": "C1-H", "atoms": ["1+H"]}]}))
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[CH3:1][CH2:2][OH:3]",
                                      "--sites", str(p)])
    with pytest.raises(_ReachedCompute):
        em.main()

# ------------------------------------------------------------- provenance (H14)
def _water(charge=0, spin=0):
    from pyscf import gto
    return gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692",
                 basis="sto-3g", charge=charge, spin=spin, verbose=0)


def test_geometry_record_hashes_the_nuclei_only():
    """Review 3 H14: the JSON had no coordinates. The sha1 covers symbols and
    positions, so an ion and a radical on the same nuclei share it."""
    g0 = em.geometry_record(_water(), 0, 0)
    g1 = em.geometry_record(_water(charge=1, spin=1), 1, 1)
    assert g0["symbols"] == ["O", "H", "H"] and len(g0["xyz_angstrom"]) == 3
    assert g0["sha1"] == g1["sha1"] and (g1["charge"], g1["spin"]) == (1, 1)
    assert g0["xyz_angstrom"][1][1] == pytest.approx(0.7572, abs=1e-12)
    assert "not charge or spin" in g0["sha1_covers"]


def test_geometry_round_trips_through_atoms_exactly():
    """The record is how to recompute a run: --atoms built from it must give
    the same nuclei to float precision."""
    from pyscf import gto
    m = _water()
    g = em.geometry_record(m, 0, 0)
    m2 = gto.M(atom=em.atoms_from_geometry(g), basis="sto-3g", verbose=0)
    assert np.abs(m2.atom_coords() - m.atom_coords()).max() < 1e-12
    assert em.geometry_record(m2, 0, 0)["sha1"] == g["sha1"]


def test_code_identity_is_unknown_outside_a_repo_that_tracks_this_file(tmp_path,
                                                                        monkeypatch):
    """`git rev-parse HEAD` reported whatever repository happened to enclose
    the module, and a clean HEAD for an edited tree. Untracked -> "unknown"."""
    fake = tmp_path / "esp_map.py"
    fake.write_text("")
    monkeypatch.setattr(em, "__file__", str(fake))
    assert em.code_identity()["code"] == "unknown"


def test_provenance_records_what_moves_results():
    p = em.provenance()
    for k in ("espmap_version", "code", "python", "numpy", "scipy", "scikit_image",
              "pyscf", "rdkit", "threads"):
        assert p[k] not in (None, ""), k
    assert isinstance(p["threads"], int) and p["threads"] >= 1


def test_espmap_version_is_the_module_that_ran(monkeypatch):
    """Review 5 R18: espmap_version read the installed metadata, so a stale
    0.1.0 editable install labelled every P3 record made by 0.1.2 code."""
    from importlib import metadata
    monkeypatch.setattr(metadata, "version", lambda name: "0.1.0")
    ident = em.code_identity()
    assert ident["espmap_version"] == em.__version__
    assert ident["installed_distribution_version"] == "0.1.0"


def test_the_version_has_one_source():
    """pyproject.toml cannot drift from __version__ (test_claims checks
    CITATION.cff against it)."""
    from pathlib import Path
    try:
        import tomllib
    except ImportError:                     # Python 3.10
        tomllib = pytest.importorskip("tomli")
    root = Path(__file__).resolve().parents[1]
    proj = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    assert "version" not in proj["project"] and "version" in proj["project"]["dynamic"]
    assert proj["tool"]["setuptools"]["dynamic"]["version"] == {"attr": "esp_map.__version__"}


def test_notebook_pins_rdkit_and_stops_on_a_different_geometry():
    """Review 3 H14: RDKit was unpinned, so the GPU check printed DIFFERENT on
    the geometry hash and went on comparing energies of different molecules."""
    from pathlib import Path
    nb = (Path(__file__).resolve().parents[1] / "notebooks" / "build_notebook.py").read_text(
        encoding="utf-8")
    assert "'rdkit==" in nb and "'pyscf==" in nb
    gpu = nb[nb.index("GPU vs CPU check"):]
    assert gpu.index("raise RuntimeError(f\"geometry hash") < gpu.index("em.run_scf(")


def test_code_identity_matches_git_describe_inside_this_repo():
    """The positive case: where git can read this checkout and tracks the
    module, the identity is `git describe --dirty`, not "unknown"."""
    import os
    import subprocess
    here = os.path.dirname(os.path.abspath(em.__file__))
    try:
        d = subprocess.run(["git", "describe", "--always", "--dirty", "--abbrev=12"],
                           cwd=here, capture_output=True, text=True, timeout=10)
    except OSError:
        pytest.skip("no git")
    if d.returncode != 0:
        pytest.skip("git cannot read this checkout")
    assert em.code_identity()["code"] == d.stdout.strip()


# ------------------------------------------------------------- notebook (H10)
# Review 3 H10. The notebook's logic is tested from the GENERATED .ipynb: the
# code a user runs, which must also be what build_notebook.py produces.
from pathlib import Path  # noqa: E402
_NB_DIR = Path(__file__).resolve().parents[1] / "notebooks"


def _nb_code():
    import json
    nb = json.loads((_NB_DIR / "espmap_colab.ipynb").read_text(encoding="utf-8"))
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


def _nb_function(name, **env):
    """Exec one top-level function out of the notebook's code cells."""
    for src in _nb_code():
        lines = src.splitlines()
        for i, line in enumerate(lines):
            if line.startswith(f"def {name}("):
                body = [line]
                for nxt in lines[i + 1:]:
                    if nxt.strip() and not nxt[0].isspace():
                        break
                    body.append(nxt)
                ns = dict(env)
                exec("\n".join(body), ns)
                return ns[name]
    raise AssertionError(f"no function {name} in the notebook")


def test_committed_notebook_is_what_the_builder_writes(tmp_path):
    """CLAUDE.md: never hand-edit the .ipynb. A stale or hand-edited notebook
    would also make the tests below test the wrong code."""
    import subprocess
    out = tmp_path / "nb.ipynb"
    subprocess.run([sys.executable, str(_NB_DIR / "build_notebook.py"), str(out)],
                   check=True, capture_output=True)
    norm = lambda p: p.read_text(encoding="utf-8").replace("\r\n", "\n")
    assert norm(out) == norm(_NB_DIR / "espmap_colab.ipynb")


def test_install_marker_follows_the_branch(tmp_path):
    """The marker ignored BRANCH: after one install, changing it neither
    cloned nor reinstalled, and the old code kept running silently."""
    write = _nb_function("write_install_marker")
    current = _nb_function("install_marker_current")
    m = tmp_path / "ESPMAP_READY"
    assert not current(str(m), "master")                  # nothing installed yet
    write(str(m), "master", "abc1234")
    assert current(str(m), "master")
    assert not current(str(m), "v0.2.0")                  # branch changed: reinstall
    assert "abc1234" in m.read_text()


def test_cache_key_covers_colour_range_and_code_and_trusts_only_finished_runs():
    run = next(s for s in _nb_code() if "4 · Run" in s)
    settings = run[run.index("settings = dict("):run.index("tag = ")]
    assert "lim=LIM" in settings and "espmap=ESPMAP_COMMIT" in settings
    # settings.json is written last, so only a finished run is reused
    assert 'os.path.isfile(os.path.join(JOB, "settings.json"))' in run
    assert "if os.path.isfile(npz)" not in run


def test_an_all_ion_series_gets_a_shared_sequential_scale():
    """An all-ion series passed the old check and was drawn on a symmetric
    +/-60 diverging scale: the failure the notebook's own text says it refuses."""
    plan = _nb_function("series_scale", em=em)
    acetate, benzoate = np.linspace(-170, -80, 50), np.linspace(-158, -41, 50)
    kind, vlim = plan([acetate, benzoate], 60)
    assert kind == "sequential" and vlim == pytest.approx((-170.0, -41.0))
    kind, vlim = plan([np.linspace(-30, 20, 50), np.linspace(-12, 40, 50)], 60)
    assert kind == "diverging" and vlim == (-60, 60)
    with pytest.raises(ValueError, match="separate"):
        plan([acetate, np.linspace(-30, 20, 50)], 60)    # mixed: still refused


@pytest.mark.parametrize("phrase", [
    "ESP extrema move little with conformation",   # 5.2 kcal/mol once folded
    "about ±3.5 kcal/mol in absolute values",      # one signed, molecule-specific shift
    "decanoic acid ~1 min",                         # a 16-thread laptop figure
])
def test_notebook_drops_review_3_stale_claims(phrase):
    text = (_NB_DIR / "build_notebook.py").read_text(encoding="utf-8")
    assert phrase not in text


# ------------------------------------------------- --optimize identity (review 4 N1)
class _SCFStop(Exception):
    pass


def _optimize_harness(monkeypatch, optimised_geometry, stop_on_second_scf=True):
    """Run main() with --optimize and no real SCF: run_scf builds the Mole
    only, and berny's optimize returns `optimised_geometry` (or the input)."""
    import types
    from pyscf import gto
    calls = []

    def fake_run_scf(spec, basis="sto-3g", charge=0, spin=0, **k):
        calls.append(spec)
        if len(calls) > 1 and stop_on_second_scf:
            raise _SCFStop
        mol = gto.M(atom=spec, basis="sto-3g", charge=charge, spin=spin, verbose=0)
        return mol, types.SimpleNamespace(mol=mol)

    def fake_optimize(mf, maxsteps=50):
        if optimised_geometry is None:
            return mf.mol
        return gto.M(atom=em.atoms_from_geometry(optimised_geometry), basis="sto-3g",
                     charge=mf.mol.charge, verbose=0)

    # stand in for pyscf's berny_solver itself, so neither pyberny nor a real
    # optimisation is needed (CI's [test] extra does not install pyberny)
    monkeypatch.setitem(sys.modules, "berny", types.ModuleType("berny"))
    monkeypatch.setitem(sys.modules, "pyscf.geomopt.berny_solver",
                        types.SimpleNamespace(optimize=fake_optimize))
    monkeypatch.setattr(em, "run_scf", fake_run_scf)
    return calls


def _neutral_glycine():
    import json
    from pathlib import Path
    p = Path(__file__).resolve().parents[1] / "docs" / "review-4" / "out" / "r8.json"
    return json.loads(p.read_text())["geometry"]


def test_optimisation_that_moves_a_proton_is_refused(monkeypatch):
    """Review 4 N1: the zwitterion's proton moves in a gas-phase optimisation,
    and the run computed neutral glycine under structure_verified: true."""
    _optimize_harness(monkeypatch, _neutral_glycine())
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[NH3+]CC(=O)[O-]",
                                      "--optimize", "--basis", "sto-3g"])
    with pytest.raises(SystemExit, match="after --optimize"):
        em.main()


def test_forced_past_a_failed_reverification_it_warns(monkeypatch, capsys):
    calls = _optimize_harness(monkeypatch, _neutral_glycine())
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[NH3+]CC(=O)[O-]",
                                      "--optimize", "--basis", "sto-3g", "--force"])
    with pytest.raises(_SCFStop):
        em.main()
    assert len(calls) == 2
    assert "after --optimize" in capsys.readouterr().err


def test_an_optimisation_that_keeps_the_structure_proceeds(monkeypatch):
    calls = _optimize_harness(monkeypatch, None)       # optimizer returns the input
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[CH3:1][OH:2]",
                                      "--optimize", "--basis", "sto-3g"])
    with pytest.raises(_SCFStop):                       # reached the second SCF
        em.main()
    assert len(calls) == 2


def test_optimizing_a_supplied_pose_is_refused_before_compute(monkeypatch):
    """Optimising a kept file geometry discards exactly what it was for."""
    import types
    from pathlib import Path
    _stop_before_scf(monkeypatch)
    monkeypatch.setitem(sys.modules, "berny", types.ModuleType("berny"))
    f = Path(__file__).resolve().parents[1] / "examples" / "real_3d.mol"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--file", str(f), "--optimize"])
    with pytest.raises(SystemExit, match="would move the coordinates"):
        em.main()
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--file", str(f), "--optimize",
                                      "--force"])
    with pytest.raises(_ReachedCompute):
        em.main()


# --------------------------------------------- --atoms is gated (review 4 N3)
_WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
_TWO_WATERS = _WATER + "; O 8 0 0.1173; H 8 0.7572 -0.4692; H 8 -0.7572 -0.4692"
_NITRITE = "N 0 0 0.45; O 0 1.06 -0.2; O 0 -1.06 -0.2"


def _cli(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["esp_map.py", *args])
    em.main()


def test_two_molecules_in_atoms_are_refused_before_compute(monkeypatch):
    """Review 4 N3: two waters 8 A apart came out validity.input true with a
    site table; the gate (which --atoms skipped) refuses them."""
    _stop_before_scf(monkeypatch)
    with pytest.raises(SystemExit):
        _cli(monkeypatch, "--atoms", _TWO_WATERS)


@pytest.mark.parametrize("args", [
    ("--atoms", _WATER),
    ("--atoms", _NITRITE, "--charge", "-1"),
    ("--atoms", _TWO_WATERS, "--force"),
    ("--molecule", "chlorobenzene"),
    (),                                            # the default --molecule
])
def test_gated_coordinates_reach_compute(monkeypatch, args):
    _stop_before_scf(monkeypatch)
    with pytest.raises(_ReachedCompute):
        _cli(monkeypatch, *args)


@pytest.mark.parametrize("atoms, charge", [
    ("O 0 0 0; H 0 0", None),                      # syntax
    ("Xx 0 0 0; H 0 0 1", None),                   # not an element
    (_NITRITE, None),                              # anion without --charge
])
def test_unreadable_or_uncharged_atoms_are_refused(monkeypatch, atoms, charge):
    _stop_before_scf(monkeypatch)
    extra = [] if charge is None else ["--charge", charge]
    with pytest.raises(SystemExit):
        _cli(monkeypatch, "--atoms", atoms, *extra)


def test_atoms_can_be_checked_without_computing(monkeypatch):
    with pytest.raises(SystemExit) as e:
        _cli(monkeypatch, "--atoms", _WATER, "--check")
    assert e.value.code == 0
    with pytest.raises(SystemExit) as e:          # rejected: a message, never 0
        _cli(monkeypatch, "--atoms", _TWO_WATERS, "--check")
    assert e.value.code not in (0, None)


# ------------------------------------------ no heavy atom (review 4 N7)
def test_site_parents_without_a_heavy_atom_is_all_none():
    """H2 passed the gate, then crashed after the SCF: argmin of an empty list."""
    assert em.site_parents([[0, 0, 0], [0, 0, 1.4]], ["H", "H"]) == [None, None]


def test_sites_block_names_why_a_hydrogen_only_molecule_has_no_sites():
    valid = {"input": True, "surface": True, "sites": False}
    P = np.array([[0.0, 0, 2], [1.0, 0, 2], [0, 1.0, 2]])
    sites, why = em.sites_block(valid, P, np.array([[0, 1, 2]]), np.zeros(3),
                                np.array([[0, 0, 0], [0, 0, 0.74]]), ["H", "H"])
    assert sites is None and "no heavy atom" in why


# ------------------------------------------- basis spellings (review 4 N9)
_SPELLINGS = ["def2-tzvp", "def2tzvp", "DEF2_TZVP", "Def2-TZVP"]


def test_pyscf_resolves_the_spellings_to_one_basis():
    """The premise: these are one calculation, so they must get one key."""
    from pyscf import gto
    nao = {gto.M(atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587", basis=b,
                 verbose=0).nao for b in _SPELLINGS}
    assert len(nao) == 1


def test_one_basis_gets_one_analysis_key_however_it_is_spelled():
    """Review 4 N9: def2-tzvp, def2tzvp and DEF2_TZVP gave three keys."""
    keys = {_keys(a={"basis": b})[0] for b in _SPELLINGS}
    assert len(keys) == 1
    assert em.analysis_fields(**{**_A, "basis": "DEF2_TZVP"})["basis"] == "def2tzvp"


@pytest.mark.parametrize("name", ["def2svp", "DEF2-SVP", "def2_svp", "sto3g", "STO-3G",
                                  "6-31G*"])
def test_double_zeta_is_recognised_however_it_is_spelled(name):
    """--basis def2svp skipped the +3.5 kcal/mol double-zeta warning."""
    assert em.basis_key(name) in em.SMALL_BASIS
    assert em.basis_key("def2-tzvp") not in em.SMALL_BASIS


@pytest.mark.parametrize("name", ["cc-pVDZ", "aug-cc-pVDZ", "6-31G(d)", "6-31+G*",
                                  "6-31G**", "6-31G(d,p)", "def2-SVPD", "pcseg-1"])
def test_review5_double_zeta_names_warn(name):
    """Review 5 R13: none of these warned."""
    assert em.basis_key(name) in em.SMALL_BASIS


@pytest.mark.parametrize("name", ["def2-tzvp", "cc-pVTZ", "6-311G**", "pcseg-2"])
def test_triple_zeta_does_not_warn(name):
    assert em.basis_key(name) not in em.SMALL_BASIS


def test_every_double_zeta_name_is_a_basis_pyscf_runs():
    from pyscf import gto
    for b in em.SMALL_BASIS_NAMES:
        assert gto.M(atom="O 0 0 0; H 0 0.76 0.5; H 0 -0.76 0.5", basis=b).nao > 0, b


@pytest.mark.parametrize("a, b", [("6-31G(d)", "6-31g*"), ("6-31G(d,p)", "6-31G**"),
                                  ("6-31+G(d,p)", "6-31+g**")])
def test_pople_polarisation_spellings_share_a_key(a, b):
    """Review 5 R13: 6-31G(d) and 6-31g* are one basis (same nao), two keys."""
    from pyscf import gto
    water = "O 0 0 0; H 0 0.76 0.5; H 0 -0.76 0.5"
    assert gto.M(atom=water, basis=a).nao == gto.M(atom=water, basis=b).nao
    assert _keys(a={"basis": a})[0] == _keys(a={"basis": b})[0]


@pytest.mark.parametrize("xc", ["b3lyp", "B3LYP", "b3lypg", "HYB_GGA_XC_B3LYP"])
def test_one_functional_gets_one_key_however_it_is_spelled(xc):
    """Review 5 R6: these resolve to one functional (one energy to 1e-8 Ha)
    and got three keys; the key now hashes what PySCF resolves."""
    assert em.analysis_fields(**{**_A, "xc": xc})["xc"] == "b3lyp"
    assert _keys(a={"xc": xc})[0] == _keys()[0]


def test_b3lyp_with_vwn5_is_a_different_key():
    assert _keys(a={"xc": "b3lyp5"})[0] != _keys()[0]
    assert em.xc_key("M06-2X") == em.xc_key("m062x") == "m062x"


def test_an_unnamed_functional_is_keyed_by_its_terms():
    k = em.xc_key("0.5*HF + 0.5*B88, LYP")
    assert k.startswith("libxc:hyb(0.5,") and k == em.xc_key("0.5*HF+0.5*B88,LYP")


@pytest.mark.parametrize("label", [
    "etkdg_mmff94s", "etkdg_uff", "mmff94s", "mmff94", "uff", "idealised", "undeclared",
    "undeclared+mmff94s_h", "undeclared+added_h", "dft:b3lyp/def2tzvp/water",
    "dft:b3lyp/def2tzvp/water+mmff94s_h", "dft:pbe0/631g*", "dft:b3lyp5/sto3g",
])
def test_every_geometry_method_label_can_be_declared(label):
    """Review 5 R6: a record re-run through --atoms must be able to say how its
    coordinates were made, or it is keyed "undeclared" and the key changes."""
    assert em.declared_geometry_method(label) == label
    assert em.geometry_method(atoms=True, declared=label) == label


def test_spellings_in_a_declared_dft_label_are_normalised():
    assert em.declared_geometry_method("dft:B3LYPG/6-31G(d)/Water") == "dft:b3lyp/631g*/water"
    assert em.declared_geometry_method("dft:nonsense/sto-3g") is None


def test_molecule_and_the_same_atoms_declared_idealised_share_a_label():
    """Review 5 R6: --molecule and its own atom string through --atoms had
    two keys (idealised_builtin / atoms_string)."""
    assert em.geometry_method(molecule="chlorobenzene") == \
        em.geometry_method(atoms=True, declared="idealised")


def test_reproduce_argv_reaches_compute_with_the_records_method(monkeypatch):
    rec = dict(geometry=dict(symbols=["O", "H", "H"],
                             xyz_angstrom=[[0, 0, 0.1173], [0, 0.7572, -0.4692],
                                           [0, -0.7572, -0.4692]]),
               input=dict(charge=0, spin=0),
               settings=dict(xc="b3lyp", basis="sto-3g", solvent="water", iso=0.002,
                             spacing_bohr=0.26, margin_bohr=5.0, density_fit=True,
                             include_reaction_field=True),
               key_fields=dict(analysis=dict(geometry_method="etkdg_mmff94s")))
    argv = em.reproduce_argv(rec)
    assert argv[argv.index("--geometry-origin") + 1] == "etkdg_mmff94s"
    assert "--density-fit" in argv and "--include-reaction-field" in argv
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", *argv])
    with pytest.raises(_ReachedCompute):
        em.main()


# --------------------------------------- declared geometry origin (review 5 R1)
@pytest.mark.parametrize("text, label", [
    ("dft:b3lyp/def2-tzvp/water", "dft:b3lyp/def2tzvp/water"),
    ("DFT:B3LYP/def2-SVP", "dft:b3lyp/def2svp"),
    ("B3LYP/def2-TZVP, PCM water, Cartesian BFGS gtol 3e-4 Ha/A", "dft:b3lyp/def2tzvp/water"),
    ("PBE0/def2-SVP", "dft:pbe0/def2svp"),
    ("MMFF94s conformer 4 of 100 (from_smiles seed 0xC0FFEE)", "mmff94s"),
    ("mmff94s", "mmff94s"), ("UFF", "uff"), ("etkdg_mmff94s", "etkdg_mmff94s"),
    ("ETKDG_UFF", "etkdg_uff"),
    ("drawn by hand", None), ("", None), (None, None),
])
def test_declared_geometry_origin_is_normalised(text, label):
    """The fixture SDFs' own free text and the label vocabulary both read."""
    assert em.declared_geometry_method(text) == label


def _nitrite_files(tmp_path):
    """One nitrite geometry, written twice: an SDF that declares a DFT origin
    and an XYZ that cannot declare anything."""
    import math
    from rdkit import Chem
    m = Chem.MolFromSmiles("N(=O)[O-]")
    c = Chem.Conformer(3)
    t = math.radians(35)                       # bench/fixtures/nitrite.sdf, tilted rigidly
    for i, (x, y) in enumerate([(0.0037, 0.4427), (-1.0666, -0.2126), (1.0628, -0.2301)]):
        c.SetAtomPosition(i, (x, y * math.cos(t), y * math.sin(t)))
    c.Set3D(True)
    m.AddConformer(c)
    sdf = tmp_path / "nitrite.sdf"
    m.SetProp("geometry", "B3LYP/def2-TZVP, PCM water, Cartesian BFGS")
    w = Chem.SDWriter(str(sdf)); w.write(m); w.close()
    xyz = tmp_path / "nitrite.xyz"
    Chem.MolToXYZFile(m, str(xyz))
    return sdf, xyz


def test_a_declared_dft_file_and_an_undeclared_file_get_different_keys(tmp_path):
    """Review 5 R1 validation: a DFT-geometry nitrite SDF and an XYZ of the
    same atoms were both file_frozen and shared analysis_key 85d0b5cbef."""
    sdf, xyz = _nitrite_files(tmp_path)
    labels = []
    for f in (sdf, xyz):
        info = em.from_file(str(f), charge=-1)[1]
        labels.append(em.geometry_method(info["source"], info["forcefield"],
                                         declared=info["declared_geometry"]))
    assert labels == ["dft:b3lyp/def2tzvp/water", "undeclared"]
    keys = {em.record_key(em.analysis_fields(**{**_A, "geometry_method": g})) for g in labels}
    assert len(keys) == 2


def test_a_dft_file_and_optimize_at_that_level_share_a_label():
    assert (em.geometry_method("file 3D coordinates used as-is",
                               declared="B3LYP/def2-TZVP, PCM water")
            == em.geometry_method(optimize=("b3lyp", "def2-tzvp", "water")))


@pytest.mark.parametrize("method, has_anion_row", [
    ("etkdg_mmff94s", True), ("undeclared", True), ("mmff94s", True),
    ("undeclared", True), ("dft:b3lyp/def2tzvp/water", False),
    ("dft:b3lyp/def2tzvp/water+mmff94s_h", False),
])
def test_the_anion_geometry_row_follows_the_geometry_method(method, has_anion_row):
    """Review 5 R1: the DFT-geometry anion fixtures carried the MMFF -> DFT
    6.5 kcal/mol row, an uncertainty they do not have."""
    ids = _budget_ids(basis="def2-tzvp", solvent="water", charge=-1, n_rot=0,
                      geometry_method=method)
    assert ("anion geometry" in ids) is has_anion_row


@pytest.mark.parametrize("basis, has_basis_row", [
    ("def2-svp", True), ("STO-3G", True), ("def2-tzvp", False), ("def2-tzvpp", False),
])
def test_the_basis_row_is_for_double_zeta_runs(basis, has_basis_row):
    """It measures def2-SVP -> def2-TZVPP; review 5 found it on def2-TZVP runs."""
    ids = _budget_ids(basis=basis, solvent=None, charge=0, n_rot=0,
                      geometry_method="etkdg_mmff94s")
    assert ("basis" in ids) is has_basis_row


@pytest.mark.parametrize("args, match", [
    (["--smiles", "CCO", "--geometry-origin", "mmff94s"], "needs --file or --atoms"),
    (["--atoms", "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692",
      "--geometry-origin", "a guess"], "not a recognised origin"),
])
def test_geometry_origin_flag_is_checked_before_compute(monkeypatch, args, match):
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", *args])
    with pytest.raises(SystemExit, match=match):
        em.main()


def test_atoms_with_a_declared_origin_reach_compute(monkeypatch):
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--atoms",
                                      "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692",
                                      "--geometry-origin", "dft:b3lyp/def2-tzvp"])
    with pytest.raises(_ReachedCompute):
        em.main()


def test_a_flag_that_contradicts_the_file_is_refused_before_compute(tmp_path, monkeypatch):
    """The file's own declaration wins; a disagreeing --geometry-origin is an
    error, not a silent override. An agreeing one is accepted."""
    sdf, _ = _nitrite_files(tmp_path)
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--file", str(sdf),
                                      "--geometry-origin", "mmff94s"])
    with pytest.raises(SystemExit, match="contradicts the file"):
        em.main()
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--file", str(sdf),
                                      "--geometry-origin", "dft:b3lyp/def2-tzvp/water"])
    with pytest.raises(_ReachedCompute):
        em.main()


def test_the_explicit_property_is_read(tmp_path):
    from rdkit import Chem
    sdf, _ = _nitrite_files(tmp_path)
    m = next(iter(Chem.SDMolSupplier(str(sdf), removeHs=False)))
    m.ClearProp("geometry")
    m.SetProp("espmap_geometry_origin", "mmff94s")
    p = tmp_path / "declared.sdf"
    w = Chem.SDWriter(str(p)); w.write(m); w.close()
    info = em.from_file(str(p), charge=-1)[1]
    assert em.geometry_method(info["source"], declared=info["declared_geometry"]) == "mmff94s"


def test_a_declared_etkdg_file_shares_the_smiles_label():
    """The C18 fixtures are ETKDG + MMFF94s conformers, what --smiles makes;
    declared so, they get the same label, hence the same key (review 5 R6)."""
    info = em.from_smiles("CCO")[1]
    assert (em.geometry_method("file 3D coordinates used as-is", declared="etkdg_mmff94s")
            == em.geometry_method(info["source"], info["forcefield"]))


# ------------------------------------- the gate sees --charge/--spin (review 5 R2)
def _codes(findings):
    return {c for _, c, _, _ in findings}


@pytest.mark.parametrize("smiles, charge, blocked", [
    ("c1ccccc1", 2, True),              # the review's case: a dication called benzene
    ("CC(=O)[O-]", 0, True),            # acetate computed as neutral
    ("CC(=O)[O-]", -1, False),          # agrees with the drawing
    ("c1ccccc1", None, False),          # no override given
    ("c1ccccc1", 0, False),
])
def test_a_charge_that_is_not_the_drawn_one_is_blocked(smiles, charge, blocked):
    verdict, findings = _verdict(smiles, charge=charge)
    assert ("CHARGE_OVERRIDE" in _codes(findings)) is blocked
    if blocked:
        assert verdict == em.BLOCK


@pytest.mark.parametrize("smiles, spin, code, sev", [
    ("O=O", 2, "SPIN_OVERRIDE", "WARN"),                 # a deliberate triplet: runs, recorded
    ("c1ccccc1", 2, "SPIN_OVERRIDE", "WARN"),            # the review's case
    ("[CH3]", 1, "OPEN_SHELL_UNRESTRICTED", "INFO"),     # agrees with the drawing
    ("[CH3]", 3, "SPIN_OVERRIDE", "WARN"),               # disagrees: was an INFO note
])
def test_a_spin_that_is_not_the_drawn_one_is_recorded(smiles, spin, code, sev):
    verdict, findings = _verdict(smiles, spin=spin)
    hit = [s for s, c, _, _ in findings if c == code]
    assert hit and hit[0] == getattr(em, sev)
    assert verdict != em.BLOCK


def test_cli_refuses_a_charge_override_before_compute(monkeypatch):
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "c1ccccc1", "--charge", "2"])
    with pytest.raises(SystemExit, match="rejected"):
        em.main()
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "c1ccccc1", "--charge", "2",
                                      "--force"])
    with pytest.raises(_ReachedCompute):
        em.main()


def test_atoms_charge_is_not_an_override(monkeypatch):
    """--atoms perceives its bonds with --charge, so the two always agree."""
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--atoms", _NITRITE, "--charge", "-1"])
    with pytest.raises(_ReachedCompute):
        em.main()


# ------------------------------------- unbound anion in the record (review 5 R3)
class _MO:
    def __init__(self, e, occ):
        self.mo_energy, self.mo_occ = np.asarray(e, float), np.asarray(occ, float)


@pytest.mark.parametrize("e, occ, charge, solvent, bad", [
    ([-0.5, 0.0063, 0.3], [2, 2, 0], -1, None, True),     # acetate, def2-SVP gas: +0.17 eV
    ([-0.5, -0.20, 0.3], [2, 2, 0], -1, "water", False),   # the same in PCM water
    ([-0.5, 0.0063, 0.3], [2, 2, 0], 0, None, False),      # not an anion
    ([[-0.5, 0.01], [-0.5, 0.3]], [[1, 1], [1, 0]], -1, None, True),   # UKS
])
def test_unbound_anion_is_an_electronic_problem(e, occ, charge, solvent, bad):
    rep = em.electronic_report(_MO(e, occ), charge, solvent)
    assert bool(rep["problems"]) is bad
    if bad:
        assert rep["problems"][0].startswith("UNBOUND_ANION_GAS: ")
    v = em.validity({"forced": False}, True, {"problems": []}, electronic=rep)
    assert v["electronic"] is (not bad)


def test_frontier_orbitals_in_eV_and_lumo_absent_without_virtuals():
    homo, lumo = em.frontier_orbitals_eV(_MO([-0.5, -0.2, 0.1], [2, 2, 0]))
    assert homo == pytest.approx(-0.2 * em.HA2EV) and lumo == pytest.approx(0.1 * em.HA2EV)
    assert em.frontier_orbitals_eV(_MO([-0.5, -0.2], [2, 2]))[1] is None


def test_sites_are_withheld_when_the_electronic_structure_is_not_valid():
    valid = {"input": True, "surface": True, "electronic": False, "sites": False}
    P = np.array([[0.0, 0, 2], [1.0, 0, 2], [0, 1.0, 2]])
    sites, why = em.sites_block(valid, P, np.array([[0, 1, 2]]), np.zeros(3),
                                np.array([[0, 0, 0], [0, 0, 1.2]]), ["C", "O"])
    assert sites is None and "electronic" in why


# ------------------------------------------------ P4: ALIE, without an SCF
def _fake_mf(C, e, n):
    from types import SimpleNamespace
    return SimpleNamespace(mo_coeff=np.asarray(C, float), mo_energy=np.asarray(e, float),
                           mo_occ=np.asarray(n, float))


def _h2_and_points():
    from pyscf import gto
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="6-31g", verbose=0)
    pts = np.array([[0.0, 0.0, 0.7], [1.5, 0.3, -0.4], [0.2, -2.0, 2.5]])
    return mol, pts


def test_alie_one_occupied_orbital_is_its_energy_everywhere():
    """docs/p4-spec.md acceptance test 1: exact for ANY coefficients, so a
    wrong weight or a missing |e| cannot hide."""
    mol, pts = _h2_and_points()
    rng = np.random.default_rng(1)
    C = rng.normal(size=(mol.nao, mol.nao))
    I = em.alie_at_points(mol, _fake_mf(C, [-0.6, -0.2, 0.3, 0.9], [2, 0, 0, 0]), pts)
    assert np.allclose(I, 0.6, atol=1e-12)


def test_alie_is_the_density_weighted_mean_of_occupied_energies():
    from pyscf.dft import numint
    mol, pts = _h2_and_points()
    rng = np.random.default_rng(2)
    C = rng.normal(size=(mol.nao, mol.nao))
    e, n = [-0.9, -0.4, 0.2, 0.7], [2, 2, 0, 0]          # virtuals must not count
    ao = numint.eval_ao(mol, pts)
    w = 2 * (ao @ C[:, :2]) ** 2
    want = (w @ np.abs(e[:2])) / w.sum(axis=1)
    assert np.allclose(em.alie_at_points(mol, _fake_mf(C, e, n), pts), want, atol=1e-12)


def test_alie_uks_with_equal_spins_equals_rks():
    """Acceptance test 2 without an SCF: UKS sums alpha and beta with n = 1,
    so identical spin channels must reproduce RKS (n = 2) exactly."""
    mol, pts = _h2_and_points()
    rng = np.random.default_rng(3)
    C = rng.normal(size=(mol.nao, mol.nao))
    e = [-0.9, -0.4, 0.2, 0.7]
    rks = em.alie_at_points(mol, _fake_mf(C, e, [2, 2, 0, 0]), pts)
    uks = em.alie_at_points(mol, _fake_mf([C, C], [e, e], [[1, 1, 0, 0], [1, 1, 0, 0]]), pts)
    assert np.allclose(uks, rks, atol=1e-12)


def test_alie_open_shell_counts_each_spin_with_its_own_energy():
    mol, pts = _h2_and_points()
    rng = np.random.default_rng(4)
    C = rng.normal(size=(mol.nao, mol.nao))
    # alpha: two occupied; beta: one, with different energies
    mf = _fake_mf([C, C], [[-0.8, -0.3, 0.1, 0.5], [-0.7, 0.0, 0.2, 0.6]],
                  [[1, 1, 0, 0], [1, 0, 0, 0]])
    from pyscf.dft import numint
    ao = numint.eval_ao(mol, pts)
    wa, wb = (ao @ C[:, :2]) ** 2, (ao @ C[:, :1]) ** 2
    want = (wa @ [0.8, 0.3] + wb @ [0.7]) / (wa.sum(axis=1) + wb.sum(axis=1))
    assert np.allclose(em.alie_at_points(mol, mf, pts), want, atol=1e-12)


def test_alie_is_nan_where_no_orbital_reaches():
    mol, _ = _h2_and_points()
    C = np.eye(mol.nao)
    I = em.alie_at_points(mol, _fake_mf(C, [-0.6, -0.2, 0.3, 0.9], [2, 0, 0, 0]),
                          np.array([[0.0, 0.0, 200.0]]))
    assert np.isnan(I).all()


# review 6 S9: |e| turned an unbound occupied orbital into a plausible ALIE
@pytest.mark.parametrize("e_occ", [[0.05], [-0.5, 0.02], [0.0]])
def test_alie_refuses_an_unbound_occupied_orbital(e_occ):
    mol, pts = _h2_and_points()
    C = np.eye(mol.nao)
    k = len(e_occ)
    e = list(e_occ) + [0.6] * (mol.nao - k)
    n = [2] * k + [0] * (mol.nao - k)
    with pytest.raises(em.ALIEUnbound):
        em.alie_at_points(mol, _fake_mf(C, e, n), pts)


@pytest.mark.slow
def test_alie_refuses_the_gas_phase_enolate():
    """Review 6's reproduction: the gas-phase enolate's HOMO is above zero at
    B3LYP/def2-SVP (+1.38 eV), and |e| gave a plausible 2.19 eV minimum."""
    spec, _ = em.from_smiles("C=C[O-]")
    mol, mf = em.run_scf(spec, basis="def2-svp", charge=-1)
    with pytest.raises(em.ALIEUnbound):
        em.alie_at_points(mol, mf, np.array([[0.0, 0.0, 4.0]]))


# ---------------------------------------- review 6 LOW findings (S18-S23)
def _separate_triangles(values):
    """One unit right triangle per value, no shared vertices: equal areas."""
    P, faces, V = [], [], []
    for k, v in enumerate(values):
        base = np.array([3.0 * k, 0.0, 0.0])
        P += [base, base + [1, 0, 0], base + [0, 1, 0]]
        faces.append([3 * k, 3 * k + 1, 3 * k + 2])
        V += [v, v, v]
    return np.array(P), np.array(faces), np.array(V, float)


def test_sigma2_tot_is_null_when_a_sign_is_below_the_floor():
    vals = [5.0 + 0.1 * i for i in range(49)] + [-10.0 - 0.1 * i for i in range(951)]
    g = em.gipf_stats(*_separate_triangles(vals))           # 4.9 % positive
    assert g["sigma2_pos_kcal2_mol2"] is None and g["sigma2_tot_kcal2_mol2"] is None


def test_sigma2_tot_is_the_one_variance_when_the_other_sign_has_no_area():
    vals = [-10.0 - 0.1 * i for i in range(100)]
    g = em.gipf_stats(*_separate_triangles(vals))
    assert g["sigma2_tot_kcal2_mol2"] == pytest.approx(g["sigma2_neg_kcal2_mol2"])


def test_the_geometry_label_follows_what_b3lyp_means():
    import subprocess
    import sys
    code = ("import pyscf.__config__ as c; c.B3LYP_WITH_VWN5 = True\n"
            "import esp_map as em\n"
            "print(em.xc_key('b3lyp'), em._dft_label('b3lyp', 'def2-svp', None))")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(Path(em.__file__).parent)).stdout.split()
    assert out == ["b3lyp5", "dft:b3lyp5/def2svp"]


def test_the_geometry_label_is_unchanged_for_plain_b3lyp():
    assert em._dft_label("B3LYP", "def2-SVP", "water") == "dft:b3lyp/def2svp/water"


def test_a_connectivity_failure_does_not_also_report_the_atom_order():
    from rdkit import Chem
    probe = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    spec, _ = em.from_smiles("CCO")
    rows = [r.split() for r in spec.split(";")]
    rows[1][1:], rows[2][1:] = rows[2][1:], rows[1][1:]       # swap C2 and O
    swapped = "; ".join(" ".join(r) for r in rows)
    verified, _, mapped_ok = em.check_identity(swapped, probe, "CCO", False, True, "test")
    assert verified is False and mapped_ok is None
    assert "atom order" not in em.input_invalid_reason(None, verified, mapped_ok)


@pytest.mark.parametrize("basis", ["sto-6g", "pcseg-0", "6-31++g", "6-31g(2d)", "6-31g(2df,p)",
                                   "def2-svp", "cc-pvdz", "6-31g*"])
def test_double_zeta_and_smaller_bases_are_small(basis):
    assert em.basis_is_small(basis)


@pytest.mark.parametrize("basis", ["def2-tzvp", "def2-tzvpp", "cc-pvtz", "6-311g**", "def2-qzvp",
                                   "pcseg-2", "aug-cc-pvtz"])
def test_triple_zeta_and_larger_bases_are_not_small(basis):
    assert not (em.basis_key(basis) in em.SMALL_BASIS)
    if hasattr(em, "basis_is_small"):
        assert not em.basis_is_small(basis)


def test_alie_chunking_and_spins_against_a_direct_reference():
    """S18: the identities cannot see a dropped beta channel when alpha and
    beta agree, beta given alpha's energies, or a chunk skipped. Distinct
    spins, 60 points and a budget that forces several chunks can."""
    from pyscf.dft import numint
    mol, _ = _h2_and_points()
    rng = np.random.default_rng(7)
    pts = rng.normal(scale=1.5, size=(60, 3))
    Ca, Cb = rng.normal(size=(2, mol.nao, mol.nao))
    ea, eb = [-0.9, -0.5, 0.2, 0.6], [-0.7, -0.3, 0.1, 0.4]
    na, nb = [1, 1, 0, 0], [1, 0, 0, 0]
    mf = _fake_mf([Ca, Cb], [ea, eb], [na, nb])
    ao = numint.eval_ao(mol, pts)
    wa, wb = (ao @ Ca[:, :2]) ** 2, (ao @ Cb[:, :1]) ** 2
    want = (wa @ [0.9, 0.5] + wb @ [0.7]) / (wa.sum(1) + wb.sum(1))
    tiny = 16 * 8 * 2 * mol.nao / 1e9                          # 16 points per chunk
    assert np.allclose(em.alie_at_points(mol, mf, pts, mem_budget_gb=tiny), want, atol=1e-12)


class _GaussMol:
    """Two atoms, an analytic density; enough for extract_surface."""
    nelectron = 2

    def atom_coords(self):
        return np.array([[0.0, 0.0, 0.0], [2.6, 0.4, -0.3]])


def _gauss_density(mol, dm, origin, n, spacing, block=20000):
    ax = [origin[i] + np.arange(n[i]) * spacing for i in range(3)]
    P = np.stack(np.meshgrid(*ax, indexing="ij"), axis=-1)
    return sum(np.exp(-2.0 * np.linalg.norm(P - c, axis=-1))
               for c in mol.atom_coords())


@pytest.mark.xfail(strict=True, reason="R7-6: margin moves the lattice, not the key")
def test_r7_e11_margin_changes_the_key_or_nothing(monkeypatch):
    """Review 7 R7-E11. Measured on b11d7c4 (chlorobenzene, butanol, def2-SVP,
    0.26 bohr): margins 5.0 / 5.13 / 5.20 give one analysis_key but V_min
    moves 0.13 kcal/mol and a site mean 0.20. Either margin enters the key, or
    the surface does not depend on it."""
    monkeypatch.setattr(em, "density_on_grid", _gauss_density)
    surf = {}
    for m in (5.0, 5.13, 5.20):
        v, f, _ = em.extract_surface(_GaussMol(), None, 0.26, iso=0.002, margin=m)
        surf[m] = {tuple(x) for x in np.round(v, 9)}
    identical = surf[5.0] == surf[5.13] == surf[5.20]
    keyed = "margin" in inspect.signature(em.analysis_fields).parameters
    assert identical or keyed
