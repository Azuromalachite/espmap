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

# what the rule on master (v0.2.0) gets wrong; a new rule removes entries
UNSAFE_NOW = {
    "draw_cyclohexanol_stretched_1.5.mol",   # non-uniform drawing -> computed flat
    "draw_chlorobenzene_zjitter_1.5.mol",    # z noise in a drawing -> computed flat
    "real_C6F6.mol", "real_C2F4.mol",        # real, uniform within 3 % -> discarded
    # added after the first candidate rule (section K of make_set.py):
    # coplanar drawings off z = 0, and z noise -> computed flat
    "draw_chlorobenzene_tilted_1.5.mol", "draw_cyclohexanol_tilted_1.5.mol",
    "draw_thf_tilted_1.5.mol", "draw_chlorobenzene_at_z5_1.5.mol",
    "draw_cyclohexanol_at_z5_1.5.mol", "draw_thf_at_z5_1.5.mol",
    "draw_thf_znoise_0.002_1.45.mol", "draw_thf_znoise_0.02_1.45.mol",
    "draw_thf_znoise_0.05_1.45.mol",
    "real_isoxazole_noH.mol", "real_anisole_noH.mol",   # real skeletons discarded
}
NOT_PREFERRED_NOW = UNSAFE_NOW | {
    "draw_hexane_1.0.mol", "draw_hexane_1.4.mol", "draw_hexane_1.45.mol",
    "draw_hexane_1.5.mol", "draw_hexane_1.54.mol", "draw_benzene_1.5.mol",
    "real_benzene_skeleton.mol", "real_anti_butane_skeleton.mol", "dft_nitrite.mol",
    "real_CO2.mol", "real_azide.mol", "real_nitrate.mol", "real_carbonate.mol",
    "real_BF3.mol", "real_SO3.mol", "real_ozone.mol", "real_nitronium.mol",
    "real_water.mol", "real_HCN_heavy.mol", "real_N2.mol", "real_CO.mol",
    "draw_hexane_squash_0.8_1.5.mol", "draw_hexane_squash_1.25_1.5.mol",
    "draw_cyclobutane_1.45.mol", "real_dimethylether_noH.mol", "real_cyclopropane_noH.mol",
}


def _verdict(path):
    mol, probs, _ = em.read_structure(str(path))
    if mol is not None:
        return mol.GetProp("_coords")
    assert probs and probs[0][1] == "DIMENSION_AMBIGUOUS", probs
    return "ambiguous"


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
