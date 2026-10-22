"""Fast tests: no SCF, no DFT. These run in seconds and gate every push.

Every case here is a defect that actually shipped at some point, lifted from
docs/adversarial-review-1.md and docs/adversarial-review-2.md. A test that was
never red is not evidence, so each one names the finding it guards.
"""
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


def test_xyz_connectivity_is_perceived(tmp_path):
    """Review 2, G3: xyz carries no bonds, so every xyz file was rejected as a
    salt with a remedy about counter-ions."""
    p = tmp_path / "w.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.117\nH 0.0 0.757 -0.469\n"
                 "H 0.0 -0.757 -0.469\n")
    mol, problems, _ = em.read_structure(str(p))
    assert mol is not None and not problems
    assert mol.GetNumBonds() == 2


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
