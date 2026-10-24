"""The P3 fixture records: real records, checked against the frozen contract
and against the files they were computed from. No SCF: reads fixtures/.

ESPMAP_FIXTURES=<dir> points the tests at another set (a --basis sto-3g dry
run of fixtures/make_records.py); the settings test is then skipped, since a
dry run overrides the basis on purpose.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

import esp_map as em
from record_tree import load_schema, tree_problems

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fixtures"))
import make_records as mk  # noqa: E402

DIR = Path(os.environ.get("ESPMAP_FIXTURES", ROOT / "fixtures" / "records"))
DRY = "ESPMAP_FIXTURES" in os.environ
NAMES = list(mk.FIXTURES)
# until the first real run is committed; the P3 PR commits all six records
pytestmark = pytest.mark.skipif(not (DIR / "benzene.json").exists(),
                                reason=f"no fixture records in {DIR} yet")
C18 = [n for n in NAMES if n.startswith("Z3Z13")]


def _rec(name):
    return json.loads((DIR / f"{name}.json").read_text())


def _npz(name):
    return np.load(DIR / f"{name}.npz")


@pytest.mark.parametrize("name", NAMES)
def test_fixture_matches_the_record_schema(name):
    rec, schema = _rec(name), load_schema()
    assert rec["schema"] == em.RECORD_SCHEMA
    assert tree_problems(schema["record"], rec) == []
    assert sorted(_npz(name).files) == sorted(schema["npz"])


@pytest.mark.parametrize("name", NAMES)
def test_fixture_is_valid_and_verified(name):
    rec = _rec(name)
    # every validity field, electronic included (review 5 R3): no fixture is an
    # anion left unbound by its basis
    assert rec["validity"] == {"input": True, "surface": True, "electronic": True,
                               "sites": True}
    assert rec["results"]["electronic"]["problems"] == []
    assert rec["input"]["structure_verified"] is True
    assert rec["sites"]["identity"]["mapped_match"] is True
    assert rec["sites_reason"] is None


@pytest.mark.parametrize("name", NAMES)
def test_fixture_keys_recompute(name):
    rec = _rec(name)
    assert rec["analysis_key"] == em.record_key(rec["key_fields"]["analysis"])
    assert rec["render_key"] == em.record_key(rec["key_fields"]["render"])


@pytest.mark.parametrize("name", NAMES)
def test_fixture_geometry_is_its_committed_file_unmoved(name):
    """Every atom supplied, nothing relaxed: the record's coordinates are the
    SDF's, and its sha1 recomputes."""
    from rdkit import Chem
    rec, g = _rec(name), _rec(name)["geometry"]
    mol = Chem.MolFromMolFile(str(mk.FIXTURES[name]["file"]), removeHs=False)
    xyz = mol.GetConformer().GetPositions()
    assert g["symbols"] == [a.GetSymbol() for a in mol.GetAtoms()]
    assert np.abs(np.array(g["xyz_angstrom"]) - xyz).max() < 1e-4
    canon = "\n".join(f"{s} {x:.8f} {y:.8f} {z:.8f}"
                      for s, (x, y, z) in zip(g["symbols"], g["xyz_angstrom"]))
    assert g["sha1"] == em._sha1(canon)
    gp = rec["sites"]["geometry_provenance"]
    assert gp["relaxation"] == "none" and gp["heavy_atoms"] == "supplied_frozen"
    assert gp["hydrogens"]["supplied"] == g["symbols"].count("H")


@pytest.mark.parametrize("name", NAMES)
def test_fixture_sites_conserve_the_surface(name):
    rec = _rec(name)
    total = rec["results"]["surface_stats"]["area_A2"]
    for rule in em.PARTITION_RULES:
        atoms = sum(a[rule]["area_A2"] for a in rec["sites"]["atoms"])
        assert atoms == pytest.approx(total, rel=1e-9)
        if rec["sites"]["reaction_sites"][0]["source"] == "default":
            sites = sum(s[rule]["area_A2"] for s in rec["sites"]["reaction_sites"])
            assert sites == pytest.approx(total, rel=1e-9)


@pytest.mark.parametrize("name", NAMES)
def test_fixture_npz_is_the_same_record(name):
    """The NPZ's partition reproduces the JSON's atom areas: a consumer that
    regroups sites from the NPZ starts from the same surface."""
    rec, z = _rec(name), _npz(name)
    assert str(z["analysis_key"]) == rec["analysis_key"]
    assert json.loads(str(z["validity"])) == rec["validity"]
    assert str(z["geometry_sha1"]) == rec["geometry"]["sha1"]
    a = em.triangle_areas(z["verts"] * em.BOHR, z["faces"])
    for rule in em.PARTITION_RULES:
        owner = z[f"tri_atom_{rule}"]
        for at in rec["sites"]["atoms"]:
            assert a[owner == at["index"]].sum() == pytest.approx(
                at[rule]["area_A2"], rel=1e-9, abs=1e-12)


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.skipif(DRY, reason="a dry run's basis override selects the basis row")
def test_fixture_budget_is_the_one_its_settings_select(name):
    """Explicit row ids per fixture, not budget_rows() recomputed: review 5 R1
    found the old test checking the producer against itself, so it passed
    while DFT geometries carried the MMFF->DFT anion row.
    - small four: DFT geometries (no anion-geometry row), def2-TZVP (no
      double-zeta basis row), PCM water (no gas->solvent rows), rigid;
    - C18: def2-SVP (basis row), gas (both solvent rows), flexible."""
    always = {"functional", "grid", "partition small anions", "partition C18", "grid mean"}
    expect = (always | {"basis", "solvent water", "solvent chloroform", "conformer C18"}
              if name in C18 else always)
    assert {r["id"] for r in _rec(name)["sites"]["budget"]} == expect


@pytest.mark.skipif(DRY, reason="a dry run overrides the basis on purpose")
@pytest.mark.parametrize("name", NAMES)
def test_fixture_settings_are_the_declared_protocol(name):
    fx, s = mk.FIXTURES[name], _rec(name)["settings"]
    assert em.basis_key(s["basis"]) == em.basis_key(fx["basis"])
    assert s["solvent"] == (fx["solvent"] or "vacuum")
    assert s["density_fit"] is fx["density_fit"]


@pytest.mark.skipif(DRY, reason="a dry run is not a committed set")
def test_fixtures_come_from_one_clean_named_revision():
    """Review 5 R17: nothing stopped a partial regeneration mixing revisions,
    or a record whose code was 'unknown'."""
    prov = [_rec(n)["provenance"] for n in NAMES]
    codes = {p["code"] for p in prov}
    assert len(codes) == 1, codes
    code = codes.pop()
    assert code != "unknown" and not code.endswith("-dirty"), code
    for k in ("espmap_version", "pyscf", "rdkit", "numpy", "threads"):
        assert len({p[k] for p in prov}) == 1, k


def version_agrees_with_code(version, code):
    """Review 5 R18: the version a record names must fit the code that made
    it. An exact release tag vX.Y.Z means version == X.Y.Z. A build after a
    tag, vX.Y.Z-N-g<sha>, names the release being built, so version > X.Y.Z:
    the fixtures are made before the v0.2.0 tag exists, because the tag
    needs them merged. Unknown or dirty code is refused elsewhere."""
    import re
    from packaging.version import Version
    m = re.fullmatch(r"v(\d+\.\d+\.\d+)(-\d+-g[0-9a-f]+)?", code)
    if m is None:
        return False
    return (Version(version) == Version(m[1]) if m[2] is None
            else Version(version) > Version(m[1]))


@pytest.mark.parametrize("version, code, ok", [
    ("0.2.0", "v0.2.0", True), ("0.2.0", "v0.1.2-90-g1a2b3c4d5e6f", True),
    ("0.1.0", "v0.1.2-82-g878332291704", False),       # the stale-metadata records
    ("0.1.2", "v0.1.2-82-g878332291704", False), ("0.1.2", "v0.2.0", False),
    ("0.2.0", "unknown", False),
])
def test_version_rule(version, code, ok):
    assert version_agrees_with_code(version, code) is ok


@pytest.mark.skipif(DRY, reason="a dry run is not a committed set")
def test_fixture_version_agrees_with_its_code():
    prov = _rec(NAMES[0])["provenance"]
    assert version_agrees_with_code(prov["espmap_version"], prov["code"]), prov


# V_min, V_max and the mean of the site holding atom 0 (nearest-nucleus), as
# computed. The records reproduced to <= 1e-10 kcal/mol across checkouts and
# thread counts (#44, review 5), so 1e-6 catches a changed result, not noise.
PINNED = {
    "benzene": (-18.410642400668962, 21.500566795728, 2.127480708211835),
    "nitrite": (-167.75227445962423, -132.13490018495082, -140.77754117328),
    "thiocyanate": (-141.26273840198374, -94.82782124016789, -118.32279758131904),
    "enolate": (-166.67195269173814, -78.27554378692678, -114.16333654057361),
    "Z3Z13-18OAc_extended": (-40.737152865309504, 19.87057637821997, -1.376987102617478),
    "Z3Z13-18OAc_folded": (-35.702106899216, 22.20914819815874, 3.685052850206189),
}


# The contract PINNED was computed under (review 7 R7-7). A change that moves
# a pinned number under an unchanged analysis_key is a new contract: bump
# em.MEASUREMENT_CONTRACT and this together, never PINNED alone.
PINNED_UNDER = "espmap.numeric/1"


def test_pinned_values_name_their_contract():
    assert em.MEASUREMENT_CONTRACT == PINNED_UNDER
    assert PINNED_UNDER in em.MEASUREMENT_CONTRACTS


@pytest.mark.skipif(DRY, reason="a dry run uses another basis")
@pytest.mark.parametrize("name", NAMES)
def test_fixture_values_are_pinned(name):
    rec = _rec(name)
    site = next(s for s in rec["sites"]["reaction_sites"] if 0 in s["atoms"])
    got = (rec["results"]["V_min_kcal_mol"], rec["results"]["V_max_kcal_mol"],
           site["nearest_nucleus"]["V_area_mean_kcal_mol"])
    assert got == pytest.approx(PINNED[name], abs=1e-6)


def _site_means(name, prefix):
    return {s["label"]: s["nearest_nucleus"]["V_area_mean_kcal_mol"]
            for s in _rec(name)["sites"]["reaction_sites"] if s["label"].startswith(prefix)}


def _noise_tolerance(molecule):
    """1.5x the largest symmetric-site spread mesh noise produced over 10
    random rigid rotations at the fixtures' own settings (review 5 R7:
    bench/symmetry_noise.py). A rotation changes nothing physical, so a
    spread beyond the noise means the partition broke symmetry."""
    rows = [json.loads(l) for l in (ROOT / "bench" / "symmetry_noise.jsonl").read_text()
            .splitlines() if l.strip()]
    spreads = [r["site_mean_spread"] for r in rows if r["molecule"] == molecule]
    assert len(spreads) == 10, molecule
    return 1.5 * max(spreads)


def test_symmetry_tolerances_come_from_the_noise_bench():
    # nitrite O sites 0.010-0.043, benzene C sites 0.041-0.081 kcal/mol
    assert _noise_tolerance("nitrite") == pytest.approx(0.0645, abs=1e-3)
    assert _noise_tolerance("benzene") == pytest.approx(0.1215, abs=1e-3)


@pytest.mark.skipif(DRY, reason="the tolerance is calibrated at the fixtures' own basis")
def test_benzene_six_carbon_sites_agree():
    """Equivalent by symmetry; the independent check on the partition (N12)."""
    m = _site_means("benzene", "C")
    assert len(m) == 6 and max(m.values()) - min(m.values()) < _noise_tolerance("benzene")


@pytest.mark.skipif(DRY, reason="the tolerance is calibrated at the fixtures' own basis")
def test_nitrite_two_oxygen_sites_agree():
    a, b = _site_means("nitrite", "O").values()      # labels come from the maps (R5)
    assert abs(a - b) < _noise_tolerance("nitrite")


@pytest.mark.parametrize("name", [n for n in NAMES if n not in C18])
def test_default_site_labels_are_atom_maps(name):
    """Review 5 R5: nitrite's 'O1' was atom index 1 = map 2. On a mapped
    input a default site's number is its heavy atom's map."""
    rec = _rec(name)
    maps = {a["index"]: a["atom_map"] for a in rec["sites"]["atoms"]}
    for s in rec["sites"]["reaction_sites"]:
        heavy = [i for i in s["atoms"] if rec["geometry"]["symbols"][i] != "H"]
        assert s["source"] == "default:map"
        assert s["label"] == f"{rec['geometry']['symbols'][heavy[0]]}{maps[heavy[0]]}"


@pytest.mark.parametrize("name", C18)
def test_c18_sites_are_named_by_chain_carbon(name):
    rec = _rec(name)
    labels = [s["label"] for s in rec["sites"]["reaction_sites"]]
    assert labels == [g["label"] for g in mk.C18_SITES["groups"]]
    assert all(s["source"] == "caller:map" for s in rec["sites"]["reaction_sites"])
    by_map = {a["atom_map"]: a["element"] for a in rec["sites"]["atoms"] if a["atom_map"]}
    assert all(by_map[n] == "C" for n in range(1, 19))
    assert (by_map[21], by_map[22], by_map[31], by_map[32]) == ("C", "C", "O", "O")


def test_c18_extended_and_folded_differ_in_topology():
    assert _rec("Z3Z13-18OAc_extended")["results"]["surface"]["genus"] == 0
    assert _rec("Z3Z13-18OAc_folded")["results"]["surface"]["genus"] > 0


@pytest.mark.parametrize("name", [n for n, f in mk.FIXTURES.items() if "source" in f])
def test_fixture_input_has_its_sources_internal_geometry(name):
    """A renumbered or rotated input is the bench geometry: every interatomic
    distance agrees, so the ESP is the same molecule's."""
    from rdkit import Chem
    fx = mk.FIXTURES[name]
    X, Y = (Chem.MolFromMolFile(str(p), removeHs=False).GetConformer().GetPositions()
            for p in (fx["source"], fx["file"]))
    d = lambda P: np.linalg.norm(P[:, None] - P[None], axis=2)  # noqa: E731
    assert np.abs(d(X) - d(Y)).max() < 1e-3        # SDF stores 4 decimals


def test_chain_maps_number_the_chain_from_the_ester_oxygen():
    from rdkit import Chem
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)OCC/C=C\\CCCCCCCC/C=C\\CCCC"))
    maps = mk.chain_maps(mol)
    sym = {m: mol.GetAtomWithIdx(i).GetSymbol() for i, m in maps.items()}
    assert sym[31] == sym[32] == "O" and sym[1] == "C"
    c1 = next(i for i, m in maps.items() if m == 1)
    assert any(maps.get(n.GetIdx()) == 31 for n in mol.GetAtomWithIdx(c1).GetNeighbors())
    c3, c4 = (next(i for i, m in maps.items() if m == k) for k in (3, 4))
    assert mol.GetBondBetweenAtoms(c3, c4).GetBondTypeAsDouble() == 2.0
