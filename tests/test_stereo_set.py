"""--require-stereo against a hostile input set built before any rule
(tests/data/stereo/make_set.py, manifest.json).

Each case runs the way main() runs it: parse_for_check, then check_structure
on the hydrogen-free molecule with require_stereo=True. The verdict is "block"
when the gate returns BLOCK, else "pass".
  - safe: the verdict is in "allowed". Passing an unspecified stereo element
    computes an arbitrary isomer under a flag that promised not to.
  - preferred: the verdict is the one a good rule gives.
The xfail sets list where the current code fails. They are strict: a rule
that fixes a case must remove it from the set, so the sets only shrink.
"""
import inspect
import json
from pathlib import Path

import pytest

import esp_map as em

_DIR = Path(__file__).resolve().parent / "data" / "stereo"
_CASES = json.loads((_DIR / "manifest.json").read_text(encoding="utf-8"))["cases"]
# until the rule lands, check_structure has no such argument: master is
# measured as it stands, with nothing that can require stereo
_HAS_FLAG = "require_stereo" in inspect.signature(em.check_structure).parameters

# what master (v0.2.1, no --require-stereo) gets wrong; a rule removes entries.
# Every unspecified case passes: the gate only warns or notes about stereo.
UNSAFE_NOW = {
    # tetrahedral centres left open in SMILES, partly or wholly
    "smi_ala_unspec", "smi_threonine_half", "smi_threonine_none", "smi_tartaric_unspec",
    "smi_trihydroxyglutaric_C3_unspec", "smi_cholesterol_none", "smi_glucose_open_partial",
    # ring cis/trans, ring junctions and bridgeheads
    "smi_14dimethylcyclohexane_unspec", "smi_decalin_unspec", "smi_camphor_unspec",
    # S, P and N centres
    "smi_sulfoxide_unspec", "smi_phosphine_unspec", "smi_ammonium_unspec",
    "smi_troger_unspec", "smi_aziridine_C_unspec",
    # E/Z: alkenes, large rings, C=N, N=N, an enolate
    "smi_2butene_unspec", "smi_hexadiene_half", "smi_cyclooctene_unspec",
    "smi_cyclodecene_unspec", "smi_oxime_unspec", "smi_nmethylimine_unspec",
    "smi_azobenzene_unspec", "smi_enolate_unspec",
    # axial and helical: allene, even cumulene, exocyclic alkene, spiro, and
    # three that RDKit cannot see at all (two atropisomers and a helicene)
    "smi_allene_unspec", "smi_butatriene_unspec", "smi_ethylidene_cyclohexane",
    "smi_spiro_unspec", "smi_binol", "smi_nitro_biphenyl", "smi_helicene6",
    # drawings: no wedge, a wavy bond, an 'either' double bond, a 3D-flagged flat drawing
    "draw_ala_plain", "draw_threonine_plain", "draw_ala_wavy", "draw_sulfoxide_plain",
    "draw_14dmch_plain", "draw_2butene_either", "draw_2butene_stripped", "draw3d_ala_plain",
    # --record 2 picks the unwedged drawing from a two-record SDF
    "sdf_record2_drawing",
}
NOT_PREFERRED_NOW = UNSAFE_NOW | set()


def _gate(case, require_stereo):
    from rdkit import Chem
    inp = case["input"]
    kw = ({"smiles": inp} if case["route"] == "smiles" else
          {"path": str(_DIR / inp), "record": case.get("record")} if case["route"] == "file" else
          {"atoms": (_DIR / inp).read_text(encoding="utf-8").strip()})
    probe, probs = em.parse_for_check(**kw)
    assert probe is not None, probs
    extra = {"require_stereo": True} if (require_stereo and _HAS_FLAG) else {}
    verdict, findings = em.check_structure(Chem.RemoveHs(probe), maps_from=probe,
                                           smiles=kw.get("smiles"), **extra)
    return verdict, findings


def _verdict(case):
    return "block" if _gate(case, True)[0] == em.BLOCK else "pass"


def _params(xfail_set, cases):
    return [pytest.param(c, id=c["name"],
                         marks=[pytest.mark.xfail(strict=True, reason="--require-stereo")]
                         if c["name"] in xfail_set else [])
            for c in cases]


def test_the_manifest_is_consistent():
    names = [c["name"] for c in _CASES]
    assert len(names) == len(set(names))
    files = {c["input"] for c in _CASES if c["route"] != "smiles"}
    on_disk = {p.name for p in _DIR.iterdir()} - {"manifest.json", "make_set.py"}
    assert files == on_disk
    assert UNSAFE_NOW <= set(names) and NOT_PREFERRED_NOW <= set(names)
    for c in _CASES:
        assert c["prefer"] is None or c["prefer"] in c["allowed"], c["name"]
        # "allowed" is exactly the verdicts whose outcome is correct or a BLOCK
        safe = {v for v, o in c["outcomes"].items() if o in ("correct", "BLOCK")}
        assert set(c["allowed"]) == safe, c["name"]


@pytest.mark.parametrize("case", _CASES, ids=[c["name"] for c in _CASES])
def test_without_the_flag_no_case_is_blocked(case):
    """The set isolates stereo: nothing else in a case may BLOCK, and stereo
    alone never blocks unless the flag asks it to."""
    verdict, findings = _gate(case, False)
    assert verdict != em.BLOCK, [f for f in findings if f[0] == em.BLOCK]


@pytest.mark.parametrize("case", _params(UNSAFE_NOW, _CASES))
def test_stereo_is_required_safely(case):
    assert _verdict(case) in case["allowed"], case["why"]


@pytest.mark.parametrize("case", _params(NOT_PREFERRED_NOW,
                                         [c for c in _CASES if c["prefer"]]))
def test_stereo_requirement_gives_the_preferred_verdict(case):
    assert _verdict(case) == case["prefer"], case["why"]
