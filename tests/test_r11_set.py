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

# what the current rule gets wrong; a new rule removes entries. Against v0.2.0
# the 100 cases had 15 unsafe and 40 not preferred; the rule now has none
# unsafe and these not preferred, every one a BLOCK:
UNSAFE_NOW = set()
NOT_PREFERRED_NOW = UNSAFE_NOW | {
    # Kekule O=O at 1.278 A is 11 % from its reference: outside PLAUSIBLE_TOL
    "real_ozone.mol",
    # aromatic N-O is 10.0 % long against the covalent-radius reference; only
    # just over, so refused rather than discarded (DRAWN_BOND_FACTOR)
    "real_isoxazole.mol", "real_isoxazole_noH.mol",
    # a flat sp3 four-membered ring is refused: oxetane is nearly planar
    "draw_cyclobutane_1.45.mol",
    # a flat sp3-sp3 ring 0.05 A off its plane is refused, not called a
    # drawing, so a nearly flat real ring (1,3-cyclohexadiene) is not discarded
    "draw_thf_znoise_0.05_1.45.mol",
    # C-C up to +6.6 %, past the 5 % for sp3 carbon, so not a geometry; with
    # one element pair an implausible file is refused, not called a drawing
    "draw_hexane_squash_1.25_1.5.mol",
    # three C-C at 1.502 A, one pair: a drawing length; ring angles are the
    # same drawn or real
    "real_cyclopropane_noH.mol",
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


# Off the frozen set: classes an independent probe of #60's first head found
# (squashed chains, tilted, moved and noisy drawings, flat cyclobutane), and
# real flat rings whose lone CH2 is planar in reality too.

def _flat_file(tmp_path, name, smiles, length=None, sy=1.0, tilt=0.0, dz=0.0,
               noise=0.0, real=False, flatten=True):
    import math
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem
    if real:      # MMFF94s heavy atoms, rotated into the xy plane
        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
        AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=5000)
        m = Chem.RemoveHs(m)
        X = m.GetConformer().GetPositions()
        X -= X.mean(axis=0)
        X = X @ np.linalg.svd(X)[2].T
        if flatten:
            X[:, 2] = 0.0
    else:         # an RDKit depiction at a uniform bond length
        m = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(m, canonOrient=True)
        X = m.GetConformer().GetPositions() * length / 1.5
        X[:, 1] *= sy
        X[:, 2] = 0.0
    t = math.radians(tilt)
    X = X @ np.array([[1, 0, 0], [0, math.cos(t), -math.sin(t)], [0, math.sin(t), math.cos(t)]]).T
    X[:, 2] += dz + noise * np.where(np.arange(len(X)) % 2, 1, -1)
    c = m.GetConformer()
    for i, p in enumerate(X):
        c.SetAtomPosition(i, p.tolist())
    c.Set3D(True)
    p = tmp_path / f"{name}.mol"
    p.write_text(Chem.MolToMolBlock(m), encoding="utf-8")
    return p


@pytest.mark.parametrize("name, kw", [
    ("hexane squashed", dict(smiles="CCCCCC", length=1.5, sy=0.8)),
    ("hexane stretched", dict(smiles="CCCCCC", length=1.5, sy=1.25)),
    ("ether squashed", dict(smiles="CCOCC", length=1.45, sy=0.8)),
    ("ether stretched", dict(smiles="CCOCC", length=1.45, sy=1.25)),
    ("cyclobutane", dict(smiles="C1CCC1", length=1.45)),
    ("chlorobenzene tilted", dict(smiles="Clc1ccccc1", length=1.4, tilt=35)),
    ("chlorobenzene at z = 5", dict(smiles="Clc1ccccc1", length=1.45, dz=5.0)),
    ("chlorobenzene z noise", dict(smiles="Clc1ccccc1", length=1.5, noise=0.01)),
    ("thf z noise", dict(smiles="C1CCOC1", length=1.45, noise=0.05)),
    ("cyclohexanol tilted", dict(smiles="OC1CCCCC1", length=1.5, tilt=20)),
    ("cyclopentadiene", dict(smiles="C1=CCC=C1", length=1.45)),
])
def test_an_off_set_drawing_is_never_computed_flat(tmp_path, name, kw):
    assert _verdict(_flat_file(tmp_path, "d", **kw)) in ("drawing", "ambiguous"), name


@pytest.mark.parametrize("name, smiles", [
    ("fluorene", "c1ccc2c(c1)Cc1ccccc1-2"), ("cyclopentadiene", "C1=CCC=C1"),
    ("1,4-cyclohexadiene", "C1=CCC=CC1"), ("thiophene", "c1ccsc1"),
])
def test_a_real_flat_ring_with_a_lone_ch2_is_kept(tmp_path, name, smiles):
    """A flat ring puckers only with a bond between two sp3 ring atoms."""
    for tilt in (0, 35):
        assert _verdict(_flat_file(tmp_path, "r", smiles, real=True, tilt=tilt)) == "geometry", name


@pytest.mark.parametrize("name, smiles", [
    ("uracil", "O=c1cc[nH]c(=O)[nH]1"), ("cytosine", "Nc1cc[nH]c(=O)n1"),
    ("adenine", "Nc1ncnc2[nH]cnc12"), ("guanine", "Nc1nc2[nH]cnc2c(=O)[nH]1"),
])
def test_a_real_flat_nucleobase_is_kept(tmp_path, name, smiles):
    """RDKit calls a lactam ring's localised C-C aromatic (uracil C4-C5
    1.477 A, +5.6 %); a 5 % aromatic C-C class refused real uracil and
    cytosine (re-probe of #60 at 92492ed)."""
    for tilt in (0, 35):
        assert _verdict(_flat_file(tmp_path, "b", smiles, real=True, tilt=tilt)) == "geometry", name


def test_a_nearly_flat_real_ring_is_refused_not_discarded(tmp_path):
    """H-free 1,3-cyclohexadiene (MMFF94s) lies within FLAT_TOL of a plane
    and has a CH2-CH2 ring bond; it was discarded as a drawing until a ring
    measurably off its plane was refused instead (re-probe of #60)."""
    p = _flat_file(tmp_path, "chd", "C1=CC=CCC1", real=True, flatten=False)
    assert _verdict(p) == "ambiguous"


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
