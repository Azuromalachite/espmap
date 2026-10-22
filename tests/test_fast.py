"""Fast tests: no SCF, no DFT. These run in seconds and gate every push.

Every case here is a defect that actually shipped at some point, lifted from
docs/adversarial-review-1.md and docs/adversarial-review-2.md. A test that was
never red is not evidence, so each one names the finding it guards.
"""
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
    (None, None, [], True),                                        # --atoms
    ({"verdict": "ACCEPT", "forced": False}, True, [], True),
    ({"verdict": "BLOCK", "forced": True}, True, [], False),       # --force
    ({"verdict": "ACCEPT", "forced": False}, False, [], False),    # H2 mismatch
])
def test_validity_input(gate, verified, problems, expect):
    v = em.validity(gate, verified, {"problems": problems})
    assert v["input"] is expect and v["surface"] is True and v["sites"] is False


def test_validity_surface_is_false_when_the_surface_has_problems():
    v = em.validity(None, None, em.surface_report(_diag(problems=["CLIPPED"]), 0.002))
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
    _, findings = _verdict("C[CH2]", spin=3)
    msg = " ".join(m for _, _, m, _ in findings)
    assert "does not match" in msg


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


def _flat_nitrite_mol(tmp_path, header="3D"):
    """A DFT-like nitrite in the xy plane (N-O 1.258 A, 115 deg), Kekule:
    N=O and N-O(-). One element pair, so lengths cannot decide."""
    import math
    from rdkit import Chem
    m = Chem.RWMol(Chem.MolFromSmiles("N(=O)[O-]"))
    c = Chem.Conformer(3)
    h = math.radians(115 / 2)
    c.SetAtomPosition(0, (0.0, 0.0, 0.0))
    c.SetAtomPosition(1, (1.258 * math.sin(h), 1.258 * math.cos(h), 0.0))
    c.SetAtomPosition(2, (-1.258 * math.sin(h), 1.258 * math.cos(h), 0.0))
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


def test_symmetric_flat_ion_with_a_3d_header_is_refused_not_guessed(tmp_path):
    """Kekule nitrite has two bond ORDERS but one element pair, and a DFT
    geometry gives both bonds 1.258 A. Counting orders called it a drawing and
    MMFF replaced it with N-O 1.42 A; lengths cannot tell the two apart."""
    p = _flat_nitrite_mol(tmp_path)
    mol, probs = em.parse_for_check(path=str(p), charge=-1)
    assert mol is None and probs[0][1] == "DIMENSION_AMBIGUOUS"
    with pytest.raises(ValueError):
        em.from_file(str(p))


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
    {"xc": "pbe0"}, {"geometry_method": "dft_opt:b3lyp/def2-svp"},
    {"geometry_method": "file_frozen"},
])
def test_every_method_axis_changes_the_analysis_key(change):
    assert _keys(a=change)[0] != _keys()[0]


_EX = __import__("pathlib").Path(__file__).resolve().parents[1] / "examples"


@pytest.mark.parametrize("make, expect", [
    (lambda: em.from_smiles("CCO")[1], "etkdg_mmff94s"),
    (lambda: em.from_file(str(_EX / "real_3d.mol"))[1], "file_frozen"),
    (lambda: em.from_file(str(_EX / "real_3d_noH.mol"))[1], "file_frozen+mmff94s_h"),
    (lambda: em.from_file(str(_EX / "drawn_2d.mol"))[1], "etkdg_mmff94s"),
])
def test_geometry_method_names_how_the_coordinates_were_made(make, expect):
    """Review 4 N2: the key omitted the geometry source, the largest measured
    axis on small anions (MMFF vs DFT up to 6.5 kcal/mol)."""
    info = make()
    assert em.geometry_method(info["source"], info["forcefield"]) == expect


def test_geometry_method_for_cli_only_sources():
    assert em.geometry_method(atoms=True) == "atoms_string"
    assert em.geometry_method(molecule="chlorobenzene") == "idealised_builtin"
    opt = ("B3LYP", "def2-TZVP", "Water")
    assert em.geometry_method("3D generated from SMILES by ETKDGv3",
                              optimize=opt) == "dft_opt:b3lyp/def2-tzvp/water"
    assert em.geometry_method(atoms=True, optimize=("b3lyp", "sto-3g", "vacuum")) \
        == "dft_opt:b3lyp/sto-3g"
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
    assert rep["partition"]["radius_scaled"] is None
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


def test_caller_groups_resolve_maps_and_hydrogens():
    sym, X, maps = _ethanol_mapped()
    groups = [{"label": "C1-H", "atoms": ["1+H"]}, {"label": "hydroxyl", "atoms": ["3+H"]},
              {"label": "C2", "atoms": [2]}]
    got = dict(em.resolve_site_groups(groups, sym, X, maps))
    parent = em.site_parents(X, sym)
    assert got["C1-H"] == sorted([0] + [i for i, p in enumerate(parent) if p == 0])
    assert len(got["C1-H"]) == 4 and got["C2"] == [1]
    assert got["hydroxyl"] == sorted([2] + [i for i, p in enumerate(parent) if p == 2])


def test_unmapped_input_names_atoms_by_index():
    spec, _ = em.from_smiles("CCO")
    sym, X = zip(*[(t.split()[0], [float(v) for v in t.split()[1:]]) for t in spec.split(";")])
    got = dict(em.resolve_site_groups([{"label": "O", "atoms": [2]}], list(sym), np.array(X), None))
    assert got == {"O": [2]}


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
        em.resolve_site_groups(groups, sym, X, maps)


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
    p.write_text(json.dumps([{"label": "C9-H", "atoms": ["9+H"]}]))
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[CH3:1][CH2:2][OH:3]",
                                      "--sites", str(p)])
    with pytest.raises(SystemExit, match="--sites"):
        em.main()


def test_main_lets_good_groups_through_to_the_scf(tmp_path, monkeypatch):
    import json
    p = tmp_path / "g.json"
    p.write_text(json.dumps([{"label": "C1-H", "atoms": ["1+H"]}]))
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
