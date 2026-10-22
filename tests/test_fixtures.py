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
    assert rec["validity"] == {"input": True, "surface": True, "sites": True}
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
def test_fixture_budget_is_the_one_its_settings_select(name):
    rec = _rec(name)
    s, n_rot = rec["settings"], None
    scope = rec["sites"]["conformer_scope"]
    if scope.startswith("single conformer; "):
        n_rot = int(scope.split("; ")[1].split()[0])
    solvent = None if s["solvent"] == "vacuum" else s["solvent"]
    assert rec["sites"]["budget"] == em.budget_rows(
        basis=s["basis"], solvent=solvent, charge=rec["geometry"]["charge"], n_rot=n_rot)


@pytest.mark.skipif(DRY, reason="a dry run overrides the basis on purpose")
@pytest.mark.parametrize("name", NAMES)
def test_fixture_settings_are_the_declared_protocol(name):
    fx, s = mk.FIXTURES[name], _rec(name)["settings"]
    assert em.basis_key(s["basis"]) == em.basis_key(fx["basis"])
    assert s["solvent"] == (fx["solvent"] or "vacuum")
    assert s["density_fit"] is fx["density_fit"]


def _site_means(name, prefix):
    return {s["label"]: s["nearest_nucleus"]["V_area_mean_kcal_mol"]
            for s in _rec(name)["sites"]["reaction_sites"] if s["label"].startswith(prefix)}


def test_benzene_six_carbon_sites_agree():
    """Equivalent by symmetry; the independent check on the partition (N12).
    The bench rows spread 0.19 kcal/mol at 0.26 bohr in the gas phase."""
    m = _site_means("benzene", "C")
    assert len(m) == 6 and max(m.values()) - min(m.values()) < 0.5


def test_nitrite_two_oxygen_sites_agree():
    m = _site_means("nitrite", "O")
    assert len(m) == 2 and abs(m["O1"] - m["O2"]) < 0.1


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
