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
import json
from pathlib import Path

import pytest

import esp_map as em

_DIR = Path(__file__).resolve().parent / "data" / "stereo"
_CASES = json.loads((_DIR / "manifest.json").read_text(encoding="utf-8"))["cases"]
# what the current rule gets wrong; a rule removes entries. Master (v0.2.1)
# had no --require-stereo and passed all 39 unspecified cases silently; the
# rule written against this set (stereo_left_open) leaves none.
UNSAFE_NOW = {
    # review 6 S7, added before any rule change: an unwedged drawing under a 3D
    # header with z noise keeps the chirality read from the noise, which
    # --require-stereo then treats as specified
    "draw3d_ala_znoise_0.05_s1", "draw3d_ala_znoise_0.05_s2", "draw3d_ala_znoise_0.05_s3",
    "draw3d_threonine_znoise_0.02_s1", "draw3d_threonine_znoise_0.02_s2",
    "draw3d_threonine_znoise_0.02_s3", "draw3d_threonine_znoise_0.05_s1",
    "draw3d_threonine_znoise_0.05_s3",
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
    verdict, findings = em.check_structure(Chem.RemoveHs(probe), maps_from=probe,
                                           smiles=kw.get("smiles"),
                                           require_stereo=require_stereo)
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


# Off the frozen set: the rule's own thresholds attacked before review (flat
# polycyclics that must not read as helicenes, the ortho-count boundary,
# cumulenes that are not stereogenic). Expected verdicts, not set entries.
@pytest.mark.parametrize("smiles, verdict", [
    ("c1ccc2c(c1)ccc1ccccc12", "pass"),                                  # phenanthrene
    ("c1ccc2c(c1)ccc1c3ccccc3ccc21", "pass"),                           # chrysene
    ("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1", "pass"),                       # pentacene
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61", "pass"),                 # coronene
    ("c1ccc2c(c1)ccc1ccc3ccccc3c12", "pass"),                           # [4]helicene
    ("c1ccc2c(c1)ccc1ccc3ccc4ccccc4c3c12", "block"),                    # [5]helicene
    ("c1ccc2c(c1)ccc1ccc3ccc4ccc5ccc6ccccc6c5c4c3c12", "block"),        # [7]helicene
    ("Cc1ccccc1-c1ccccc1C", "pass"),                                    # two ortho: rotates
    ("Cc1cccc(Cl)c1-c1ccccc1Br", "block"),                              # tri-ortho, chiral
    ("Cc1cccc(C)c1-c1c(C)cccc1Cl", "pass"),                             # one end symmetric
    ("c1ccc2c(-c3cccc4ccccc34)cccc2c1", "block"),                       # 1,1'-binaphthyl
    ("c1ccc(-c2ccccn2)cc1", "pass"),                                    # 2-phenylpyridine
    ("c1ccccc1C=C=Cc1ccccc1", "block"),                                 # 1,3-diphenylallene
    ("CC(C)=C=C(C)C", "pass"),                                          # tetramethylallene
    ("CC(C)=C=O", "pass"),                                              # a ketene
    ("CN=C=NC", "pass"),                                                # a carbodiimide
    ("CC=C=C=C=CC", "block"),                                           # an odd cumulene
])
def test_the_rule_on_its_own_boundaries(smiles, verdict):
    from rdkit import Chem
    got = "block" if em.stereo_left_open(Chem.MolFromSmiles(smiles)) else "pass"
    assert got == verdict, smiles


def test_the_cli_flag_blocks_before_compute(monkeypatch):
    """--check --require-stereo: an open centre exits non-zero with
    STEREO_REQUIRED; the same SMILES passes without the flag, and a fixed
    one passes with it."""
    import sys
    calls = []

    def run(*args):
        monkeypatch.setattr(sys, "argv", ["esp_map.py", *args])
        with pytest.raises(SystemExit) as e:
            em.main()
        calls.append(args)
        return e.value.code

    assert run("--smiles", "CC(N)C(=O)O", "--check", "--require-stereo") not in (0, None)
    assert run("--smiles", "CC(N)C(=O)O", "--check") == 0
    assert run("--smiles", "C[C@H](N)C(=O)O", "--check", "--require-stereo") == 0
    assert run("--atoms", (_DIR / "atoms_ala.atoms").read_text(encoding="utf-8").strip(),
               "--check", "--require-stereo") == 0


@pytest.mark.parametrize("smiles", ["CC=C=C=CC", "C/C=C=C=C/C", "C/C=C=C=C\C", "CC=C=CC"])
def test_a_cumulene_names_the_only_remedy_that_works(smiles):
    """RDKit drops / and \ on a cumulene (both butatriene forms parse to
    STEREONONE), so the remedy must not offer them: only 3D coordinates fix a
    stereogenic allene or cumulene. Held-out probe of #69 at c3e137b."""
    from rdkit import Chem
    v, findings = em.check_structure(Chem.MolFromSmiles(smiles), smiles=smiles,
                                     require_stereo=True)
    assert v == em.BLOCK
    req = next(f for f in findings if f[1] == "STEREO_REQUIRED")
    assert "cumulene" in req[3] and "3D" in req[3], req[3]
