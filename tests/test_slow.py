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
    """Reference windows asserted by --selftest (wide by design; not taken from a cited source)."""
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


@pytest.mark.parametrize("extra, requested, substituted", [
    ([], None, False),                  # default flags: the map follows the scale
    (["--cmap", "vik"], "vik", True),   # diverging map asked for, scale known only after V
])
def test_an_anion_is_drawn_sequential_end_to_end(tmp_path, monkeypatch, extra,
                                                  requested, substituted):
    """Review 3, H3, through main() where it shipped: acetate never changes
    sign, so its scale is sequential and so is the map actually drawn, with
    or without a contradicting --cmap. sto-3g still gives an all-negative
    surface (-181 to -86 kcal/mol)."""
    import json
    import sys
    out, js = tmp_path / "ac.png", tmp_path / "ac.json"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CC(=O)[O-]",
                                      "--basis", "sto-3g", "--size", "120",
                                      "--out", str(out), "--json", str(js), *extra])
    em.main()
    scale = json.loads(js.read_text())["scale"]
    assert scale["kind"] == "sequential"
    assert scale["sign_change_ratio"] == 0.0 and scale["sign_change_threshold"] == 0.05
    assert scale["cmap"] in em.SEQUENTIAL
    assert scale["cmap_requested"] == requested
    assert scale["cmap_substituted"] is substituted
    assert out.exists()


def test_json_record_carries_sigma_hole_keys_and_validity(tmp_path, monkeypatch):
    """Review 3 H6 through main(): chlorobenzene's sigma-hole V_max was printed
    and never written; V_min/V_max had no unit; the record had no keys."""
    import json
    import sys
    js, npz = tmp_path / "cb.json", tmp_path / "cb.npz"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "Clc1ccccc1",
                                      "--basis", "sto-3g", "--size", "120",
                                      "--out", str(tmp_path / "cb.png"),
                                      "--json", str(js), "--npz", str(npz)])
    em.main()
    rec = json.loads(js.read_text())
    sh = rec["results"]["sigma_hole"]
    assert sh["V_max_kcal_mol"] is not None and sh["n_cap_vertices"] > 0
    assert rec["results"]["sigma_hole_reason"] is None
    assert "V_min_kcal_mol" in rec["results"] and "V_min" not in rec["results"]
    # sites became True with P2: a valid record now carries a site block
    assert rec["validity"] == {"input": True, "surface": True, "electronic": True, "sites": True}
    assert rec["sites"]["schema"] == em.SITES_SCHEMA and rec["sites_reason"] is None
    assert rec["analysis_key"] == em.record_key(rec["key_fields"]["analysis"])
    assert rec["render_key"] == em.record_key(rec["key_fields"]["render"])
    assert rec["results"]["dipole"]["origin"]["kind"] == "centre_of_nuclear_charge"
    z = np.load(npz)
    assert str(z["analysis_key"]) == rec["analysis_key"] and "bohr" in str(z["units"])


def test_reaction_field_reads_pyscf_charges_correctly():
    """The reaction field uses PySCF's internal PCM charges, so check the
    reading against PySCF's own solvation energy, E = q . v_solute / 2, and
    the physics: around an anion the surface charge is +|Q|(eps-1)/eps and
    its potential on the molecular surface is positive."""
    spec, _ = em.from_smiles("CC(=O)[O-]")
    mol, mf = em.run_scf(spec, basis="sto-3g", charge=-1, solvent="water")
    ws = mf.with_solvent
    q, v = ws._intermediates["q"], ws._intermediates["v_grids"]
    assert 0.5 * float(np.dot(q, v)) == pytest.approx(mf.scf_summary["e_solvent"], rel=1e-6)
    eps = ws.eps
    assert q.sum() == pytest.approx((eps - 1) / eps, abs=0.02)
    verts, _, _ = em.extract_surface(mol, em.total_dm(mf), 0.4, margin=5.0)
    assert em.reaction_field_at_points(mf, verts).min() > 0


def _run_json(tmp_path, monkeypatch, argv):
    import json
    import sys
    js = tmp_path / "rec.json"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", *argv, "--basis", "sto-3g",
                                      "--size", "80", "--out", str(tmp_path / "x.png"),
                                      "--json", str(js)])
    em.main()
    return json.loads(js.read_text())


def test_a_gas_phase_unbound_anion_withholds_its_sites(tmp_path, monkeypatch):
    """Review 5 R3 through main(), on the record P2 used to compute: nitrite
    in the gas phase at STO-3G has HOMO ~ +7 eV, so its ESP is an artefact of
    the basis. The record says so, withholds the site table, and the NPZ
    carries no partition, while the record still matches its schema."""
    from record_tree import load_schema, tree_problems
    h = np.radians(116.4) / 2
    atoms = (f"N 0 0 0; O {-1.2577 * np.sin(h):.6f} {-1.2577 * np.cos(h):.6f} 0; "
             f"O {1.2577 * np.sin(h):.6f} {-1.2577 * np.cos(h):.6f} 0")
    npz = tmp_path / "g.npz"
    rec = _run_json(tmp_path, monkeypatch, ["--atoms", atoms, "--charge", "-1",
                                            "--npz", str(npz)])
    el = rec["results"]["electronic"]
    assert el["homo_eV"] > 0
    assert any(p.startswith("UNBOUND_ANION_GAS") for p in el["problems"])
    assert rec["validity"]["electronic"] is False and rec["validity"]["sites"] is False
    assert rec["validity"]["input"] is True and rec["validity"]["surface"] is True
    assert rec["sites"] is None and "electronic" in rec["sites_reason"]
    assert not any(k.startswith("tri_atom_") for k in np.load(npz).files)
    assert tree_problems(load_schema()["record"], rec) == []


def test_sites_end_to_end_on_a_symmetric_nitrite(tmp_path, monkeypatch):
    """per-site-spec P2 through main(): a C2v nitrite built symmetric, so its
    two O sites must agree; the site areas must add up to the surface."""
    h = np.radians(116.4) / 2
    atoms = (f"N 0 0 0; O {-1.2577 * np.sin(h):.6f} {-1.2577 * np.cos(h):.6f} 0; "
             f"O {1.2577 * np.sin(h):.6f} {-1.2577 * np.cos(h):.6f} 0")
    # in PCM water, as the P3 nitrite fixture is: in the gas phase at STO-3G
    # its HOMO is +7 eV, an unbound anion whose record now says so (review 5 R3)
    rec = _run_json(tmp_path, monkeypatch, ["--atoms", atoms, "--charge", "-1",
                                            "--solvent", "water"])
    assert rec["validity"]["sites"] is True and rec["sites"]["schema"] == em.SITES_SCHEMA
    total = rec["results"]["surface_stats"]["area_A2"]
    g = rec["results"]["gipf"]                      # an anion: no positive part
    assert g["sigma2_pos_kcal2_mol2"] is None and g["sigma2_neg_kcal2_mol2"] > 0
    for rule in em.PARTITION_RULES:
        s = {x["label"]: x[rule] for x in rec["sites"]["reaction_sites"]}
        assert sum(v["area_A2"] for v in s.values()) == pytest.approx(total, rel=1e-9)
        assert s["O1"]["area_A2"] == pytest.approx(s["O2"]["area_A2"], rel=0.01)
        assert s["O1"]["V_area_mean_kcal_mol"] == pytest.approx(
            s["O2"]["V_area_mean_kcal_mol"], abs=0.5)
        assert min(v["V_min_kcal_mol"] for v in s.values()) == pytest.approx(
            rec["results"]["V_min_kcal_mol"])


def test_a_generated_polar_hydrogen_is_flagged_on_its_site(tmp_path, monkeypatch):
    rec = _run_json(tmp_path, monkeypatch, ["--smiles", "CCO"])
    flags = {s["label"]: s["h_generated_polar"] for s in rec["sites"]["reaction_sites"]}
    assert flags == {"C0": False, "C1": False, "O2": True}
    assert rec["sites"]["geometry_provenance"]["hydrogens"] == {
        "supplied": 0, "generated": 0, "generated_relaxed": 6, "optimised": 0}


def test_record_matches_the_committed_schema_in_both_directions(tmp_path, monkeypatch):
    """Review 4 N5: the spec and the record had drifted (budget, compared,
    radii, relaxation). One mapped chloromethane run fills nearly every
    optional branch: sigma-hole, atom maps, sites, NPZ."""
    from record_tree import load_schema, tree_problems
    npz = tmp_path / "m.npz"
    rec = _run_json(tmp_path, monkeypatch, ["--smiles", "[CH3:1][Cl:2]", "--npz", str(npz)])
    schema = load_schema()
    assert rec["schema"] == em.RECORD_SCHEMA
    assert tree_problems(schema["record"], rec) == []
    assert sorted(np.load(npz).files) == sorted(schema["npz"])
    # A1: bonds in record order, each one a real contact in the geometry
    xyz = np.array(rec["geometry"]["xyz_angstrom"])
    assert len(rec["geometry"]["bonds"]) == 4
    assert all(np.linalg.norm(xyz[i] - xyz[j]) < 2.0 for i, j, _ in rec["geometry"]["bonds"])


@pytest.mark.parametrize("argv", [
    ["--smiles", "O", "--margin", "2.0", "--no-strict"],                    # open surface
    ["--atoms", "O 0 0 0; H 0.96 0 0; H -0.24 0.93 0; "
                "O 3 0 0; H 3.96 0 0; H 2.76 0.93 0", "--force"],              # forced BLOCK
    ["--smiles", "[CH3]", "--spin", "1"],                                    # UKS
], ids=["open_surface", "forced_block", "uks"])
def test_other_record_branches_match_the_schema(tmp_path, monkeypatch, argv):
    """R4: the one-record schema test never wrote an open surface, whose
    genus is null. These are the branches review 5 ran by hand."""
    from record_tree import load_schema, tree_problems
    rec = _run_json(tmp_path, monkeypatch, argv)
    assert tree_problems(load_schema()["record"], rec) == []

def test_json_geometry_reproduces_the_run_through_atoms(tmp_path, monkeypatch):
    """Review 3 H14 validation: re-running a JSON's recorded geometry through
    --atoms reproduces V_min/V_max to 1e-6. Review 5 R6: and the analysis_key,
    which the test never compared (etkdg_mmff94s vs atoms_string)."""
    import json
    import sys
    recs = []
    for i, src in enumerate([["--smiles", "CCO", "--basis", "sto-3g"], None]):
        js = tmp_path / f"r{i}.json"
        if src is None:
            src = em.reproduce_argv(recs[0])
        monkeypatch.setattr(sys, "argv", ["esp_map.py", *src,
                                          "--size", "80", "--out", str(tmp_path / f"r{i}.png"),
                                          "--json", str(js)])
        em.main()
        recs.append(json.loads(js.read_text()))
    a, b = (r["results"] for r in recs)
    assert abs(a["V_min_kcal_mol"] - b["V_min_kcal_mol"]) < 1e-6
    assert abs(a["V_max_kcal_mol"] - b["V_max_kcal_mol"]) < 1e-6
    assert recs[0]["geometry"]["sha1"] == recs[1]["geometry"]["sha1"]
    assert recs[0]["analysis_key"] == recs[1]["analysis_key"]
    assert recs[1]["key_fields"]["analysis"]["geometry_method"] == "etkdg_mmff94s"
    assert recs[0]["provenance"]["code"]


def test_optimize_keeps_atom_maps_and_records_the_optimisation(tmp_path, monkeypatch):
    """Review 4 N1: --optimize dropped the atom maps of a mapped input and
    left settings.geometry saying 'generated by ETKDGv3'."""
    import json
    import sys
    pytest.importorskip("berny")
    js = tmp_path / "m.json"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "[CH3:1][OH:2]",
                                      "--optimize", "--basis", "sto-3g", "--size", "80",
                                      "--out", str(tmp_path / "m.png"), "--json", str(js)])
    em.main()
    rec = json.loads(js.read_text())
    maps = [a["atom_map"] for a in rec["sites"]["atoms"] if a["atom_map"]]
    assert sorted(maps) == [1, 2]
    assert rec["input"]["structure_verified"] is True
    assert "optimised" in rec["settings"]["geometry"]
    g = rec["geometry"]
    assert g["optimised_from_sha1"] and g["optimised_from_sha1"] != g["sha1"]


def test_forced_ungated_input_is_not_recorded_as_valid(tmp_path, monkeypatch):
    """Review 4 N3 validation: the two-water --atoms record must not say
    validity.input true, and must name the reason."""
    import json
    import sys
    js = tmp_path / "w2.json"
    atoms = ("O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692; "
             "O 8 0 0.1173; H 8 0.7572 -0.4692; H 8 -0.7572 -0.4692")
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--atoms", atoms, "--force",
                                      "--basis", "sto-3g", "--size", "80",
                                      "--out", str(tmp_path / "w2.png"), "--json", str(js)])
    em.main()
    rec = json.loads(js.read_text())
    assert rec["validity"]["input"] is False and rec["validity"]["sites"] is False
    assert rec["input"]["gate"]["forced"] is True
    assert ["BLOCK", "MULTIPLE_FRAGMENTS"] in rec["input"]["gate"]["findings"]


@pytest.mark.parametrize("src", [["--atoms", "H 0 0 0; H 0 0 0.74"],
                                 ["--smiles", "[H][H]"]])
def test_hydrogen_molecule_writes_its_record(tmp_path, monkeypatch, src):
    """Review 4 N7 validation: H2 wrote the PNG, then crashed with no JSON."""
    import json
    import sys
    js = tmp_path / "h2.json"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", *src, "--basis", "sto-3g",
                                      "--size", "80", "--out", str(tmp_path / "h2.png"),
                                      "--json", str(js), "--npz", str(tmp_path / "h2.npz")])
    em.main()
    rec = json.loads(js.read_text())
    assert rec["sites"] is None and "no heavy atom" in rec["sites_reason"]
    assert rec["validity"]["sites"] is False and rec["validity"]["input"] is True
    # review 5 R9: no site table, so no partition in the NPZ either
    assert not any(k.startswith("tri_atom_") for k in np.load(tmp_path / "h2.npz").files)


def _npz_run(tmp_path, monkeypatch, *args):
    import json
    import sys
    js, npz = tmp_path / "r.json", tmp_path / "r.npz"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", *args, "--basis", "sto-3g",
                                      "--size", "80", "--out", str(tmp_path / "r.png"),
                                      "--json", str(js), "--npz", str(npz)])
    em.main()
    return json.loads(js.read_text()), np.load(npz)


def test_npz_of_an_invalid_record_carries_its_validity_and_no_partition(tmp_path,
                                                                        monkeypatch):
    """Review 4 N8 validation: a clipped surface forced through with
    --no-strict gave a JSON with validity.surface false and no sites, while
    its NPZ still carried tri_atom_* and no validity at all."""
    import json
    rec, z = _npz_run(tmp_path, monkeypatch, "--smiles", "O", "--margin", "2.0",
                      "--no-strict")
    v = json.loads(str(z["validity"]))
    assert v == rec["validity"] and v["surface"] is False
    assert json.loads(str(z["surface_problems"])) == rec["results"]["surface"]["problems"]
    assert not [k for k in z.files if k.startswith("tri_atom_")]


def test_npz_of_a_valid_record_matches_its_json(tmp_path, monkeypatch):
    import json
    rec, z = _npz_run(tmp_path, monkeypatch, "--smiles", "[CH3:1][OH:2]")
    assert json.loads(str(z["validity"])) == rec["validity"]
    assert rec["validity"]["sites"] is True
    assert [k for k in z.files if k.startswith("tri_atom_")]
    assert str(z["geometry_sha1"]) == rec["geometry"]["sha1"]
    assert sorted(int(m) for m in z["atom_maps"] if m) == [1, 2]


def test_double_zeta_warning_fires_for_an_unhyphenated_name(capsys):
    """Review 4 N9 validation: --basis def2svp ran without the warning."""
    em.run_scf("H 0 0 0; H 0 0 0.74", basis="def2svp")
    assert "double-zeta" in capsys.readouterr().err


def test_forced_charge_override_is_recorded_as_invalid_input(tmp_path, monkeypatch):
    """Review 5 R2 validation: benzene --charge 2 was recorded as verified
    benzene with validity all true; forced, it is now validity.input false
    with the CHARGE_OVERRIDE finding."""
    import json
    import sys
    js = tmp_path / "bz2.json"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "c1ccccc1", "--charge", "2",
                                      "--force", "--basis", "sto-3g", "--size", "80",
                                      "--out", str(tmp_path / "bz2.png"), "--json", str(js)])
    em.main()
    rec = json.loads(js.read_text())
    assert ["BLOCK", "CHARGE_OVERRIDE"] in rec["input"]["gate"]["findings"]
    assert rec["input"]["gate"]["forced"] is True
    assert rec["validity"]["input"] is False and rec["sites"] is None


def test_deliberate_triplet_runs_with_the_spin_finding_recorded(tmp_path, monkeypatch):
    import json
    import sys
    js = tmp_path / "o2.json"
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "O=O", "--spin", "2",
                                      "--basis", "sto-3g", "--size", "80",
                                      "--out", str(tmp_path / "o2.png"), "--json", str(js)])
    em.main()
    rec = json.loads(js.read_text())
    assert ["WARN", "SPIN_OVERRIDE"] in rec["input"]["gate"]["findings"]
    assert rec["validity"]["input"] is True


@pytest.mark.parametrize("solvent, valid", [(None, False), ("water", True)])
def test_gas_phase_unbound_anion_record_says_so(tmp_path, monkeypatch, solvent, valid):
    """Review 5 R3 validation: acetate at def2-SVP in the gas phase (HOMO
    +0.17 eV) was recorded fully valid with a site table."""
    import json
    import sys
    js = tmp_path / "ac.json"
    extra = ["--solvent", solvent] if solvent else []
    monkeypatch.setattr(sys, "argv", ["esp_map.py", "--smiles", "CC(=O)[O-]", "--basis",
                                      "def2-svp", *extra, "--size", "80",
                                      "--out", str(tmp_path / "ac.png"), "--json", str(js)])
    em.main()
    rec = json.loads(js.read_text())
    el = rec["results"]["electronic"]
    assert rec["validity"]["electronic"] is valid and (el["homo_eV"] > 0) is (not valid)
    assert (rec["sites"] is not None) is valid
    if not valid:
        assert el["problems"][0].startswith("UNBOUND_ANION_GAS")


def test_alie_rks_and_uks_agree_on_a_closed_shell():
    """docs/p4-spec.md acceptance test 2 with real SCFs: water as UKS with
    spin 0 must give RKS's ALIE (same orbitals, alpha + beta)."""
    from pyscf import dft
    water = "O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692"
    mol, rks = em.run_scf(water, basis="def2-svp")
    uks = dft.UKS(mol)
    uks.xc, uks.conv_tol = "b3lyp", rks.conv_tol
    uks.grids.level = rks.grids.level
    uks.kernel(dm0=np.array([rks.make_rdm1() / 2] * 2))
    pts = np.array([[0.0, 0.0, 3.06], [0.0, 2.9, -2.0], [1.5, 1.0, 1.0]])
    a = em.alie_at_points(mol, rks, pts)
    b = em.alie_at_points(mol, uks, pts)
    assert np.all(a > 0) and np.allclose(a, b, atol=1e-6)
