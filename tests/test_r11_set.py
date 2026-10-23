"""Review 5 R11: the drawing-vs-geometry rule against a hostile input set built
before any rule (tests/data/r11/make_set.py, manifest.json).

Two checks per input file, both on coordinate_kind as read_structure applies it:
  - safe: the verdict is in "allowed". A drawing computed as-is, or a real
    geometry discarded, fails it.
  - preferred: the verdict is the one a good rule gives.
The xfail sets list where the current rule fails. They are strict: a rule
that fixes a case must remove it from the set, so the sets only shrink.
"""
import json
from pathlib import Path

import pytest

import esp_map as em

_DIR = Path(__file__).resolve().parent / "data" / "r11"
_CASES = json.loads((_DIR / "manifest.json").read_text(encoding="utf-8"))["cases"]

# what the current rule gets wrong; a new rule removes entries. v0.2.0 had 4
# unsafe and 25 not preferred; the rule written against this set has none
# unsafe and these two not preferred, both BLOCKs:
UNSAFE_NOW = set()
NOT_PREFERRED_NOW = UNSAFE_NOW | {
    # one element pair, uniform, at 1.5 A: the same evidence as draw_CO2_1.5
    # and draw_nitrate_1.5, which prefer "ambiguous"
    "draw_benzene_1.5.mol",
    # Kekule O=O at 1.278 A is 11 % from its reference: outside PLAUSIBLE_TOL
    "real_ozone.mol",
}


def _verdict(path):
    mol, probs, _ = em.read_structure(str(path))
    if mol is not None:
        return mol.GetProp("_coords")
    assert probs and probs[0][1] == "DIMENSION_AMBIGUOUS", probs
    return "ambiguous"


def test_the_header_is_read_from_the_file_per_record(tmp_path):
    """RDKit calls 3D-headed C2F4 2D (every z = 0 and an "either" double
    bond), which discarded a real geometry. The reader takes the dimension
    from line 2 of each record instead."""
    from rdkit import Chem
    real = (_DIR / "real_C2F4.mol").read_text(encoding="utf-8")
    assert not Chem.MolFromMolBlock(real).GetConformer().Is3D()  # the RDKit quirk
    drawn = (_DIR / "real_chlorobenzene_2Dheader.mol").read_text(encoding="utf-8")
    p = tmp_path / "two.sdf"
    p.write_text(drawn.rstrip() + "\n$$$$\n" + real.rstrip() + "\n$$$$\n", encoding="utf-8")
    assert em._molfile_header_dims(p) == ["2D", "3D"]
    assert em.read_structure(str(p), record=1)[0].GetProp("_coords") == "drawing"
    assert em.read_structure(str(p), record=2)[0].GetProp("_coords") == "geometry"


def _params(xfail_set, cases):
    return [pytest.param(c, id=c["file"],
                         marks=[pytest.mark.xfail(strict=True, reason="review 5 R11")]
                         if c["file"] in xfail_set else [])
            for c in cases]


def test_the_manifest_lists_every_file_once():
    files = [c["file"] for c in _CASES]
    assert len(files) == len(set(files))
    assert set(files) == {p.name for p in _DIR.glob("*.mol")}
    assert UNSAFE_NOW <= set(files) and NOT_PREFERRED_NOW <= set(files)
    for c in _CASES:
        assert c["prefer"] is None or c["prefer"] in c["allowed"], c["file"]
        # "allowed" is exactly the verdicts whose outcome is correct or a BLOCK
        safe = {v for v, o in c["outcomes"].items() if o in ("correct", "BLOCK")}
        assert set(c["allowed"]) == safe, c["file"]


@pytest.mark.parametrize("case", _params(UNSAFE_NOW, _CASES))
def test_flat_input_is_classified_safely(case):
    assert _verdict(_DIR / case["file"]) in case["allowed"], case["why"]


@pytest.mark.parametrize("case", _params(NOT_PREFERRED_NOW,
                                         [c for c in _CASES if c["prefer"]]))
def test_flat_input_gets_the_preferred_verdict(case):
    assert _verdict(_DIR / case["file"]) == case["prefer"], case["why"]
