"""Slow tests: these run real SCF calculations. Minutes, not seconds.

Run with `pytest -m slow` or `pytest --runslow`. They are the only tests that
can catch a physics error, so they gate a release even though they do not gate
every push.
"""
import numpy as np
import pytest

import esp_map as em

pytestmark = pytest.mark.slow

WATER = "O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692"


# ----------------------------------------------------------------- physics
def test_selftest_passes():
    assert em.selftest(basis="def2-svp") == 0


def test_esp_sign_and_magnitude_on_water():
    """Reference windows from the published range of water ESP extrema."""
    mol, mf = em.run_scf(WATER, basis="def2-svp")
    dm = em.total_dm(mf)
    b = 1.0 / em.BOHR
    pts = np.array([[0.0, 0.0, 1.6173],     # lone-pair lobe
                    [0.0, 1.546, -1.055],   # beyond an O-H bond
                    [0.0, 0.0, 25.0]]) * b  # far field
    v = em.esp_at_points(mol, dm, pts) * em.HA2KCAL
    assert -70 < v[0] < -30, "lone pair must be strongly negative"
    assert 35 < v[1] < 80, "donor hydrogen must be positive"
    assert abs(v[2]) < 2, "a neutral molecule's potential must decay to zero"


def test_a_flipped_coulomb_kernel_is_caught():
    """Fault injection. int1e_grids returns the POSITIVE kernel, so the
    electronic term is subtracted; getting it backwards gives a smooth,
    structured, entirely wrong map. The far field is the sharpest detector."""
    mol, mf = em.run_scf(WATER, basis="def2-svp")
    dm = em.total_dm(mf)
    pts = np.array([[0.0, 0.0, 25.0]]) / em.BOHR
    v_ok = em.esp_at_points(mol, dm, pts)[0] * em.HA2KCAL

    nuc = sum(mol.atom_charge(a) /
              max(np.linalg.norm(pts[0] - mol.atom_coord(a)), 1e-12)
              for a in range(mol.natm))
    ints = mol.intor("int1e_grids", grids=pts)
    v_broken = (nuc + np.einsum("pij,ij->p", ints, dm)[0]) * em.HA2KCAL

    assert abs(v_ok) < 2, "correct sign decays"
    assert v_broken > 100, "flipped sign leaves an undecayed monopole"


# ------------------------------------------------------- surface integrity
def test_isosurface_is_closed_with_adequate_margin():
    mol, mf = em.run_scf(WATER, basis="def2-svp")
    _, _, diag = em.extract_surface(mol, em.total_dm(mf), 0.26, margin=5.0)
    assert diag["chi"] == 2
    assert diag["face_max_rho"] < em.ISO_BADER


def test_clipped_isosurface_is_refused():
    """A clipped surface changed a measured V_max by 26% while every other
    number in the run stayed plausible."""
    spec = em.monosubstituted_benzene("Cl", 1.740)
    mol, mf = em.run_scf(spec, basis="def2-svp")
    dm = em.total_dm(mf)
    with pytest.raises(RuntimeError, match="CLIPPED|not a single closed"):
        em.extract_surface(mol, dm, 0.26, margin=3.0, strict=True)
    # and it still renders when the user explicitly opts out
    _, _, diag = em.extract_surface(mol, dm, 0.26, margin=3.0, strict=False)
    assert diag["chi"] != 2


# -------------------------------------------------------------- sigma-hole
def test_sigma_hole_of_chlorobenzene_is_in_the_literature_range():
    spec = em.monosubstituted_benzene("Cl", 1.740)
    mol, mf = em.run_scf(spec, basis="def2-tzvp")
    dm = em.total_dm(mf)
    v, _, _ = em.extract_surface(mol, dm, 0.26, margin=5.0)
    V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
    sh, n = em.sigma_hole(v, V, mol, i_sub=6)
    assert 8.0 < sh < 12.0, f"def2-TZVP gas-phase value was +9.83, got {sh}"
    assert n > 50


def test_sigma_hole_axis_comes_from_the_bond_not_a_ring():
    """Review 2, G2: the axis was (X - mean of atoms 0..5). On CC(=O)OCCCl
    that is 40 degrees off and flipped the reported sign."""
    spec, gi = em.from_smiles("CC(=O)OCCCl")
    mol, mf = em.run_scf(spec, basis="def2-svp", charge=gi["charge"])
    dm = em.total_dm(mf)
    v, _, _ = em.extract_surface(mol, dm, 0.30, margin=5.0)
    V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
    sym = [mol.atom_symbol(i) for i in range(mol.natm)]
    i_cl = sym.index("Cl")
    sh, _ = em.sigma_hole(v, V, mol, i_sub=i_cl)
    assert sh > 0.5, f"a chloroalkyl sigma-hole is positive, got {sh}"


def test_bromine_sigma_hole_exceeds_chlorine():
    """The chemistry this tool exists to show: the sigma-hole grows down the
    group. If this inverts, the axis or the surface is wrong."""
    out = {}
    for name, (sym, d) in [("Cl", ("Cl", 1.740)), ("Br", ("Br", 1.890))]:
        mol, mf = em.run_scf(em.monosubstituted_benzene(sym, d), basis="def2-svp")
        dm = em.total_dm(mf)
        v, _, _ = em.extract_surface(mol, dm, 0.26, margin=5.0)
        V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
        out[name] = em.sigma_hole(v, V, mol, i_sub=6)[0]
    assert out["Br"] > out["Cl"] + 2.0, out


# ----------------------------------------------------------------- solvent
def test_solvent_reduces_the_sigma_hole():
    """Measured: gas +9.83 -> water +6.22 at def2-TZVP, an axis as large as
    the basis-set error. Guards the PCM wiring, not the physics."""
    spec = em.monosubstituted_benzene("Cl", 1.740)
    vals = {}
    for solv in [None, "water"]:
        mol, mf = em.run_scf(spec, basis="def2-svp", solvent=solv)
        dm = em.total_dm(mf)
        v, _, _ = em.extract_surface(mol, dm, 0.30, margin=5.0)
        V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
        vals[solv] = em.sigma_hole(v, V, mol, i_sub=6)[0]
    assert vals["water"] < vals[None] - 1.0, vals


# --------------------------------------------------------------- open shell
def test_radical_uses_an_unrestricted_wavefunction():
    """Review 2, G1: RKS accepted spin=2 and returned a converged, meaningless
    number. <S^2> is the invariant that distinguishes the two."""
    methyl = ("C 0 0 0; H 0 0 1.09; H 1.03 0 -0.36; H -0.51 0.89 -0.36")
    mol, mf = em.run_scf(methyl, basis="def2-svp", spin=1)
    assert mf.spin_square()[0] == pytest.approx(0.75, abs=0.05)
    dm = em.total_dm(mf)
    assert dm.ndim == 2, "consumers need the total (alpha+beta) density"
    s = mol.intor("int1e_ovlp")
    assert np.einsum("ij,ji->", dm, s) == pytest.approx(mol.nelectron, abs=1e-6)


# --------------------------------------------------------------- end to end
def test_render_writes_an_image(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    mol, mf = em.run_scf(WATER, basis="def2-svp")
    dm = em.total_dm(mf)
    v, f, _ = em.extract_surface(mol, dm, 0.30, margin=5.0)
    V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
    img = em.render(v, f, V, em.build_maps()["vik"], (-60, 60), size=80, ss=1)
    assert img.shape == (80, 80, 3)
    assert 0.0 <= img.min() and img.max() <= 1.0
    assert img.std() > 0.01, "a flat image means the surface was not painted"


def test_a_late_colormap_refusal_still_writes_the_data(tmp_path, monkeypatch):
    """With --scale auto the kind is known only after the SCF. A mismatched
    --cmap still refuses the picture, but the JSON and NPZ hold valid numbers
    and must be written first, so the compute is not thrown away."""
    import json
    import sys
    out, js, npz = tmp_path / "ac.png", tmp_path / "ac.json", tmp_path / "ac.npz"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CC(=O)[O-]",
                                      "--basis", "sto-3g", "--cmap", "vik",
                                      "--out", str(out), "--json", str(js),
                                      "--npz", str(npz)])
    with pytest.raises(SystemExit, match="does not fit"):
        em.main()
    assert not out.exists(), "no image in a colormap the scale contradicts"
    rec = json.loads(js.read_text())
    assert rec["scale"]["kind"] == "sequential"
    assert rec["scale"]["cmap"] is None and "does not fit" in rec["scale"]["cmap_refused"]
    assert rec["results"]["V_max"] < 0
    assert np.load(npz)["V"].max() < 0
