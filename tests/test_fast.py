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
                                      problems=["... CLIPPED ...", "... OPEN ..."]))
    folded = em.surface_report(_diag(genus=2, chi=-2))
    assert clipped["chi"] == folded["chi"] == -2          # the old JSON's only field
    assert clipped["valid"] is False and clipped["problems"] and clipped["strict"] is False
    assert folded["valid"] is True and folded["problems"] == [] and folded["genus"] == 2
    assert em.surface_report(_diag())["valid"] is True


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

    def old_reader(path, charge=0):
        m = Chem.MolFromPDBFile(str(path), removeHs=False)
        m.SetProp("_bonds", "perceived")
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


# A colour flag must not cost a DFT run: PR #11 refused a mismatched or
# unknown --cmap only after the SCF, with nothing written.
@pytest.mark.parametrize("argv", [
    ["--scale", "sequential", "--cmap", "vik"],     # kind known up front: mismatch
    ["--scale", "diverging", "--cmap", "lajolla"],
    ["--cmap", "jet"],                              # not a map this program has
])
def test_a_colour_flag_that_cannot_work_is_refused_before_the_scf(monkeypatch, argv):
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CCO", *argv])
    with pytest.raises(SystemExit, match="does not fit|must be one of"):
        em.main()


@pytest.mark.parametrize("argv", [
    ["--cmap", "vik"],                              # auto scale: kind not known yet
    ["--scale", "diverging", "--cmap", "broc"],
    ["--scale", "sequential", "--cmap", "vik", "--force-cmap"],
])
def test_a_colour_flag_that_may_work_reaches_the_scf(monkeypatch, argv):
    _stop_before_scf(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CCO", *argv])
    with pytest.raises(_ReachedCompute):
        em.main()


def test_scale_caption_matches_the_selftest_orientation():
    assert em.scale_caption("diverging") == "warm = negative"
    assert em.scale_caption("sequential") == "darker = more negative"


def test_notebook_image_caption_follows_the_scale_kind():
    """H3's caption half, in the notebook: the static image was titled 'warm =
    negative' on sequential scales too. A source check, because the notebook
    cell runs only after a real SCF."""
    from pathlib import Path
    src = (Path(__file__).resolve().parents[1] / "notebooks" / "build_notebook.py").read_text(
        encoding="utf-8")
    title = next(line for line in src.splitlines() if "plt.title(f\"{NAME}" in line)
    assert "warm = negative" not in title and "em.scale_caption(kind)" in title


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
                 "comparable": comparable}


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
