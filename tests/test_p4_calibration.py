"""P4 calibration set (docs/p4-spec.md section 6), committed before ALIE.

The fast tests check the manifest itself. The slow test runs each case with
--alie once ALIE exists. Until then every case is xfail(strict=True), so the
ALIE PR must turn each expectation green or challenge it in the open; it may
not edit an expectation to pass.
"""
import json
import re
from pathlib import Path

import pytest
from rdkit import Chem

import esp_map as em

_DIR = Path(__file__).resolve().parent / "data" / "p4_calibration"
_CASES = json.loads((_DIR / "manifest.json").read_text(encoding="utf-8"))["cases"]
_TYPES = {"lower", "not_clearly_lower", "report_only"}

# every case fails until ALIE is implemented; the ALIE PR removes entries
NOT_YET = {c["name"] for c in _CASES if c["expect"]["type"] != "report_only"}


def test_the_calibration_manifest_is_consistent():
    names = [c["name"] for c in _CASES]
    assert len(names) == len(set(names))
    on_disk = {p.name for p in _DIR.glob("*.sdf")}
    assert {c["file"] for c in _CASES} == on_disk
    for c in _CASES:
        m = next(iter(Chem.SDMolSupplier(str(_DIR / c["file"]), removeHs=False)))
        ref = Chem.MolFromSmiles(c["smiles"])
        for site, idx in c["sites"].items():
            for i in idx:
                assert i < ref.GetNumAtoms(), (c["name"], site)
                assert m.GetAtomWithIdx(i).GetSymbol() != "H", (c["name"], site)
                assert m.GetAtomWithIdx(i).GetSymbol() == ref.GetAtomWithIdx(i).GetSymbol()
        e = c["expect"]
        assert e["type"] in _TYPES, c["name"]
        if e["type"] != "report_only":
            assert set(e["a"]) | set(e["b"]) <= set(c["sites"]), c["name"]
        assert c["sources"], c["name"]
        for s in c["sources"]:
            assert re.fullmatch(r"10\.\d{4,9}/\S+", s["doi"]), s
            assert s["evidence"] in {"metadata", "abstract", "full_text"}, s


def _site_min(rec, indices):
    """The lowest per-atom ALIE minimum over `indices` (nearest-nucleus rule).
    Written against docs/p4-spec.md section 3; adapt if the record/2 keys
    change in review, never the expectation."""
    atoms = rec["sites"]["atoms"]
    return min(atoms[i]["alie"]["nearest_nucleus"]["min_eV"] for i in indices)


@pytest.mark.slow
@pytest.mark.parametrize("case", [
    pytest.param(c, id=c["name"],
                 marks=[pytest.mark.xfail(strict=True, reason="ALIE not implemented (P4)")]
                 if c["name"] in NOT_YET else [])
    for c in _CASES if c["expect"]["type"] != "report_only"])
def test_alie_meets_the_calibration_expectation(case, tmp_path, monkeypatch):
    import sys
    if not hasattr(em, "alie_at_points"):
        pytest.fail("ALIE is not implemented yet (docs/p4-spec.md step 4)")
    s, js = case["settings"], tmp_path / "r.json"
    monkeypatch.setattr(sys, "argv", [
        "esp_map.py", "--file", str(_DIR / case["file"]), "--alie",
        "--xc", s["xc"], "--basis", s["basis"], "--solvent", s["solvent"],
        "--size", "80", "--out", str(tmp_path / "r.png"), "--json", str(js)])
    em.main()
    rec = json.loads(js.read_text())
    e, sites = case["expect"], case["sites"]
    a = _site_min(rec, [i for k in e["a"] for i in sites[k]])
    b = _site_min(rec, [i for k in e["b"] for i in sites[k]])
    if e["type"] == "lower":
        assert a < b, (case["name"], a, b, case["why"])
    else:
        assert not a < b - e["tol_eV"], (case["name"], a, b, case["why"])
