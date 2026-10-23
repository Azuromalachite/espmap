"""Build the R11 hostile input set: 80 planar and near-planar structure files,
mostly exactly z = 0 with deliberate near-planar controls (a pyramidal
aniline N, one atom at z = 1e-5, a drawing with +-1e-4 A z noise), that a
drawing-vs-geometry rule must classify, each with what the gate should do.

Written BEFORE any rule (review 5 R11; two rules written first both failed on
counterexamples found afterwards). Run from the repository root:

    python tests/data/r11/make_set.py

It rewrites tests/data/r11/*.mol and manifest.json deterministically.

Verdicts are coordinate_kind's: "geometry" (coordinates used as-is),
"drawing" (discarded, a 3D conformer generated, recorded as DISCARDED),
"ambiguous" (a DIMENSION_AMBIGUOUS BLOCK). Each case lists
  - allowed: verdicts that are SAFE. A drawing must never be "geometry" (it
    would be computed flat, silently); a real geometry must never be
    "drawing" (it would be replaced). "ambiguous" is always safe: it refuses
    with a remedy.
  - prefer: the verdict a good rule gives. A BLOCK on a real geometry is
    safe but costs the user a re-save, which is what R11 complained about.
    null where the coordinates fit both readings equally (a hexagon at the
    aromatic C-C length is a drawing and benzene's skeleton).

tests/test_r11_set.py checks every case: "safe" (verdict in allowed) and
"preferred" (verdict == prefer), each with an xfail list of what the current
rule gets wrong, which a new rule must shrink.
"""
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

OUT = Path(__file__).resolve().parent
CASES = []
DRAW = ["drawing", "ambiguous"]          # a drawing: never computed as-is
REAL = ["geometry", "ambiguous"]         # a geometry: never discarded
EITHER = ["geometry", "drawing", "ambiguous"]


def _write(name, mol, flag3d, allowed, prefer, kind, source, why):
    mol = Chem.Mol(mol)
    mol.GetConformer().Set3D(flag3d)
    mol.SetProp("_Name", name)
    (OUT / f"{name}.mol").write_text(Chem.MolToMolBlock(mol), encoding="utf-8")
    assert prefer is None or prefer in allowed, name
    CASES.append(dict(file=f"{name}.mol", kind=kind, allowed=allowed, prefer=prefer,
                      outcomes=_outcomes(kind, flag3d),
                      header="3D" if flag3d else "2D", source=source, why=why))


# What each verdict does to the user, worst first, so a rule can trade misses
# honestly: a silent wrong geometry is worse than a recorded discard, which
# is worse than a BLOCK with a remedy.
SEVERITY = {"silent wrong geometry": 3, "recorded discard": 2,
            "overrules the file's 2D header": 2, "BLOCK": 1, "correct": 0}


def _outcomes(kind, flag3d):
    if not flag3d:                      # the header says 2D: the file's own word
        return {"drawing": "correct", "ambiguous": "BLOCK",
                "geometry": "overrules the file's 2D header"}
    if kind == "drawing":
        return {"drawing": "correct", "ambiguous": "BLOCK",
                "geometry": "silent wrong geometry"}
    if kind == "geometry":
        return {"geometry": "correct", "ambiguous": "BLOCK",
                "drawing": "recorded discard"}
    return {"geometry": "correct", "drawing": "correct", "ambiguous": "BLOCK"}


def _set_xy(mol, xy, z=None):
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y) in enumerate(xy):
        conf.SetAtomPosition(i, (float(x), float(y), 0.0 if z is None else float(z[i])))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def _depiction(smiles, length, hs=False, stretch=None, jitter=None):
    """RDKit 2D depiction (every bond 1.5 A) rescaled to a uniform `length`;
    stretch=(factor) scales x only (a hand-tweaked, non-uniform drawing);
    jitter puts +-jitter on z (a converter that wrote noise into z)."""
    m = Chem.MolFromSmiles(smiles)
    if hs:
        m = Chem.AddHs(m)
    AllChem.Compute2DCoords(m, canonOrient=True)
    c = m.GetConformer()
    X = np.array([[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y]
                  for i in range(m.GetNumAtoms())])
    X *= length / 1.5
    if stretch:
        X[:, 0] *= stretch
    z = None
    if jitter:
        z = [jitter * (1 if i % 2 else -1) for i in range(m.GetNumAtoms())]
    return _set_xy(m, X, z)


def _plane(mol, tol=0.05):
    """Rotate a conformer into the xy plane and set z = 0; refuse a molecule
    that is not planar to tol (A)."""
    c = mol.GetConformer()
    X = np.array([list(c.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    X -= X.mean(axis=0)
    _, _, vt = np.linalg.svd(X)
    Y = X @ vt.T
    assert np.abs(Y[:, 2]).max() < tol, (Chem.MolToSmiles(mol), np.abs(Y[:, 2]).max())
    return _set_xy(mol, Y[:, :2])


def _mmff(smiles, keep_h=True):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=2000)
    m = _plane(m)
    return m if keep_h else Chem.RemoveHs(m)


def _built(smiles, xy):
    """A small molecule at given in-plane coordinates (atom order of smiles)."""
    ps = Chem.SmilesParserParams()
    ps.removeHs = False                 # keep [H] atoms written in the SMILES
    m = Chem.MolFromSmiles(smiles, ps)
    assert m.GetNumAtoms() == len(xy), smiles
    return _set_xy(m, xy)


def _polar(r, deg):
    return (r * math.cos(math.radians(deg)), r * math.sin(math.radians(deg)))


def _fixture(name):
    sup = Chem.SDMolSupplier(str(OUT.parents[2] / "bench" / "fixtures" / f"{name}.sdf"),
                             removeHs=False)
    return _plane(next(iter(sup)))


# --- A. drawings at editor scales, 3D flag, z = 0 (noise control: A's last) ---
SCALES = (0.825, 1.0, 1.2, 1.35, 1.4, 1.45, 1.5, 1.54, 1.7, 2.0)
for L in SCALES:
    _write(f"draw_chlorobenzene_{L}", _depiction("Clc1ccccc1", L), True, DRAW, "drawing",
           "drawing", "RDKit depiction rescaled", "a drawing at a uniform editor scale")
for smi, nm in (("OC1CCCCC1", "cyclohexanol"), ("C1CCOC1", "thf"),
                ("CC(O)CN", "aminopropanol"), ("CC(C)(C)O", "tbutanol"),
                ("CCCCCC", "hexane")):
    for L in (1.0, 1.4, 1.45, 1.5, 1.54):
        _write(f"draw_{nm}_{L}", _depiction(smi, L), True, DRAW, "drawing", "drawing",
               "RDKit depiction rescaled",
               "sp3 atoms drawn flat: no real coordinates put these in a plane")
_write("draw_benzene_1.5", _depiction("c1ccccc1", 1.5), True, DRAW, "drawing", "drawing",
       "RDKit depiction", "one element pair at a drawing length")
_write("draw_chlorobenzene_withH_1.5", _depiction("Clc1ccccc1", 1.5, hs=True), True, DRAW,
       "drawing", "drawing", "RDKit depiction with explicit H",
       "some editors draw H; C-H at 1.5 A is no geometry")
_write("draw_cyclohexanol_stretched_1.5", _depiction("OC1CCCCC1", 1.5, stretch=1.12), True,
       DRAW, "drawing", "drawing", "RDKit depiction, x stretched 12 %",
       "a hand-tweaked drawing: bonds no longer uniform, still flat sp3")
_write("draw_chlorobenzene_zjitter_1.5", _depiction("Clc1ccccc1", 1.5, jitter=1e-4), True,
       DRAW, "drawing", "drawing", "RDKit depiction, z = +-1e-4",
       "a converter wrote noise into z; the x,y are still a drawing")
# --- B. drawings with an honest 2D header --------------------------------------
for L in (0.825, 1.5):
    _write(f"draw2d_cyclohexanol_{L}", _depiction("OC1CCCCC1", L), False, DRAW,
           "drawing", "drawing", "RDKit depiction, 2D header", "the header says 2D")
# --- C. indistinguishable: a drawing that is also a fair planar skeleton --------
_write("draw_benzene_1.39", _depiction("c1ccccc1", 1.39), True, EITHER, None,
       "either", "RDKit depiction at 1.39 A",
       "a regular hexagon at the aromatic C-C length IS benzene's heavy skeleton")
_write("draw_pyridine_1.4", _depiction("c1ccncc1", 1.4), True, EITHER, None,
       "either", "RDKit depiction at 1.40 A",
       "a planar aromatic at near-real lengths; computing it as-is is harmless")
# --- D. real planar geometries with hydrogens (MMFF94s, rotated into xy) --------
for smi, nm in (("c1ccccc1", "benzene"), ("c1ccncc1", "pyridine"), ("Clc1ccccc1", "chlorobenzene"),
                ("C=O", "formaldehyde"), ("C=C", "ethylene"), ("c1cc[nH]c1", "pyrrole"),
                ("c1ccoc1", "furan"), ("c1ccc2ccccc2c1", "naphthalene"), ("O=CN", "formamide")):
    _write(f"real_{nm}", _mmff(smi), True, REAL, "geometry", "geometry",
           "RDKit ETKDG + MMFF94s, rotated into the xy plane", "a real planar molecule")
# --- E. real planar geometries without hydrogens -------------------------------
for smi, nm in (("Fc1c(F)c(F)c(F)c(F)c1F", "C6F6"), ("FC(F)=C(F)F", "C2F4"),
                ("Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl", "C6Cl6")):
    _write(f"real_{nm}", _mmff(smi), True, REAL, "geometry", "geometry",
           "RDKit ETKDG + MMFF94s, rotated into the xy plane", "real, no H to break uniformity")
_write("real_benzene_skeleton", _mmff("c1ccccc1", keep_h=False), True, REAL, "geometry",
       "geometry", "MMFF94s benzene, H removed", "a real heavy-atom skeleton, one element pair")
_write("real_anti_butane_skeleton", _built("CCCC", [(0, 0), _polar(1.53, 0),
       (1.53 - 1.53 * math.cos(math.radians(112)), 1.53 * math.sin(math.radians(112))),
       (1.53 - 1.53 * math.cos(math.radians(112)) + 1.53, 1.53 * math.sin(math.radians(112)))]),
       True, REAL, "geometry", "geometry", "built: C-C 1.53 A, 112 deg, anti",
       "anti-butane's carbons are coplanar in a real geometry")
# --- F. DFT geometries from bench/fixtures (B3LYP/def2-TZVP, PCM water) ---------
for nm in ("nitrite", "thiocyanate", "benzene", "enolate"):
    _write(f"dft_{nm}", _fixture(nm), True, REAL, "geometry", "geometry",
           f"bench/fixtures/{nm}.sdf, rotated into the xy plane", "a real DFT geometry")
# --- G. small inorganic and linear species at typical literature lengths --------
SMALL = [
    ("O=C=O", "CO2", [(-1.16, 0), (0, 0), (1.16, 0)]),
    ("[N-]=[N+]=[N-]", "azide", [(-1.18, 0), (0, 0), (1.18, 0)]),
    ("[O-][N+](=O)[O-]", "nitrate", [_polar(1.25, 90), (0, 0), _polar(1.25, 210), _polar(1.25, 330)]),
    ("[O-]C(=O)[O-]", "carbonate", [_polar(1.29, 90), (0, 0), _polar(1.29, 210), _polar(1.29, 330)]),
    ("FB(F)F", "BF3", [_polar(1.31, 90), (0, 0), _polar(1.31, 210), _polar(1.31, 330)]),
    ("O=S(=O)=O", "SO3", [_polar(1.42, 90), (0, 0), _polar(1.42, 210), _polar(1.42, 330)]),
    ("[O-][O+]=O", "ozone", [_polar(1.278, 90 + 58.4), (0, 0), _polar(1.278, 90 - 58.4)]),
    ("O=[N+]=O", "nitronium", [(-1.15, 0), (0, 0), (1.15, 0)]),
    ("O", "water_heavy_only", [(0, 0)]),
    ("[H]O[H]", "water", [_polar(0.958, 90 + 52.25), (0, 0), _polar(0.958, 90 - 52.25)]),
    ("C#N", "HCN_heavy", [(0, 0), (1.156, 0)]),
    ("[H]C#N", "HCN", [(-1.066, 0), (0, 0), (1.156, 0)]),
    ("N#N", "N2", [(0, 0), (1.098, 0)]),
    ("[C-]#[O+]", "CO", [(0, 0), (1.128, 0)]),
]
for smi, nm, xy in SMALL:
    _write(f"real_{nm}", _built(smi, xy), True,
           REAL, "geometry", "geometry", "built at a typical literature bond length (rounded)",
           "a real small molecule; one atom, a diatomic or a symmetric triatomic is planar")
# --- H. single-pair species at a drawing length --------------------------------
# draw_CO2_1.5 preferred "ambiguous" until the user corrected it: the graph
# asserts C=O, and no real C=O is near 1.5 A, so the length is evidence of a
# drawing. Nitrate's Kekule N-O single bond can be 1.5 A, so it stays ambiguous.
_write("draw_CO2_1.5", _built("O=C=O", [(-1.5, 0), (0, 0), (1.5, 0)]), True, DRAW,
       "drawing", "drawing", "built at 1.5 A",
       "one pair at a drawing length, but C=O at 1.5 A is +21 % from its reference "
       "(real CO2 1.16 A): no real CO2 has it")
_write("draw_nitrate_1.5", _built("[O-][N+](=O)[O-]",
       [_polar(1.5, 90), (0, 0), _polar(1.5, 210), _polar(1.5, 330)]), True, DRAW, "ambiguous",
       "drawing", "built at 1.5 A", "one pair at a drawing length: undecidable")
# --- I. nearly planar real geometries (some z != 0) ----------------------------
m = Chem.AddHs(Chem.MolFromSmiles("Nc1ccccc1"))
AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=2000)
_write("real_aniline_pyramidal_N", m, True, REAL, "geometry", "geometry",
       "RDKit ETKDG + MMFF94s, not flattened", "NH2 is pyramidal, so some z != 0")
m = _mmff("Clc1ccccc1")
c = m.GetConformer()
for i in range(m.GetNumAtoms()):
    p = c.GetAtomPosition(i)
    c.SetAtomPosition(i, (p.x, p.y, 1e-5 if i == 0 else 0.0))
_write("real_chlorobenzene_znoise", m, True, REAL, "geometry", "geometry",
       "MMFF94s chlorobenzene with one z = 1e-5", "a real planar geometry printed with noise")
# --- J. a real planar geometry under a 2D header ---------------------------------
_write("real_chlorobenzene_2Dheader", _mmff("Clc1ccccc1"), False, DRAW, "drawing",
       "geometry", "MMFF94s chlorobenzene saved with a 2D header",
       "the header is the file's own statement; the reader cannot overrule it")


# --- K. attacks on the first candidate rule (PR #60 at cbc4c20) ---------------
# Added after that rule, from an independent probe: each class beat it or
# master. They join the set as cases like any other.
def _xyz(mol):
    c = mol.GetConformer()
    return np.array([list(c.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])


def _set_xyz(mol, X):
    return _set_xy(mol, X[:, :2], X[:, 2])


def _tilted(mol, ax=35.0, ay=20.0):
    a, b = math.radians(ax), math.radians(ay)
    rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    return _set_xyz(mol, _xyz(mol) @ (ry @ rx).T)


def _mmff_heavy(smiles):
    """A real geometry whose heavy atoms are planar (H removed, then
    rotated into xy): what a heavy-atom-only file of it looks like."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=5000)
    return _plane(Chem.RemoveHs(m))


_noise = np.random.default_rng(20261023)
for smi, nm in (("Clc1ccccc1", "chlorobenzene"), ("OC1CCCCC1", "cyclohexanol"), ("C1CCOC1", "thf")):
    _write(f"draw_{nm}_tilted_1.5", _tilted(_depiction(smi, 1.5)), True, DRAW, "drawing",
           "drawing", "RDKit depiction, rotated 35 deg about x and 20 deg about y",
           "a coplanar drawing out of the xy plane: z != 0 but still flat")
    m = _depiction(smi, 1.5)
    X = _xyz(m)
    X[:, 2] = 5.0
    _write(f"draw_{nm}_at_z5_1.5", _set_xyz(m, X), True, DRAW, "drawing", "drawing",
           "RDKit depiction translated to z = 5", "every z equal but not zero: still flat")
for zn in (0.002, 0.02, 0.05):
    m = _depiction("C1CCOC1", 1.45)
    X = _xyz(m)
    X[:, 2] = _noise.uniform(-zn, zn, len(X))
    _write(f"draw_thf_znoise_{zn}_1.45", _set_xyz(m, X), True, DRAW, "drawing", "drawing",
           f"RDKit depiction at 1.45 A, z uniform in +-{zn} A",
           "a flat sp3 ring with z noise; real puckering is ~0.2 A or more")
for smi, nm, L in (("CCCCCC", "hexane", 1.5), ("CCOCC", "diethylether", 1.45)):
    for sy in (0.8, 1.25):
        m = _depiction(smi, L)
        X = _xyz(m)
        X[:, 1] *= sy
        _write(f"draw_{nm}_squash_{sy}_{L}", _set_xyz(m, X), True, DRAW, "drawing", "drawing",
               f"RDKit depiction at {L} A, y scaled by {sy}",
               "a hand-squashed flat chain: sp3 angles neither 120 nor real")
for smi, nm in (("C1CCC1", "cyclobutane"), ("OC1CCC1", "cyclobutanol")):
    _write(f"draw_{nm}_1.45", _depiction(smi, 1.45), True, DRAW, "drawing", "drawing",
           "RDKit depiction at 1.45 A", "a flat sp3 4-ring; a real one is puckered")
for smi, nm in (("o1nccc1", "isoxazole"), ("COc1ccccc1", "anisole"), ("COC", "dimethylether"),
                ("C1CC1", "cyclopropane")):
    _write(f"real_{nm}_noH", _mmff_heavy(smi), True, REAL, "geometry", "geometry",
           "RDKit ETKDG + MMFF94s, H removed, heavy atoms rotated into xy",
           "a real heavy-atom skeleton that is planar")
_write("real_isoxazole", _mmff("o1nccc1"), True, REAL, "geometry", "geometry",
       "RDKit ETKDG + MMFF94s, rotated into the xy plane",
       "a real heteroaromatic whose N-O sits near a length tolerance")

# --- L. nucleobases and a nearly planar real ring -------------------------------
# Added after the probe of PR #60 at 92492ed refused real uracil: nucleobases
# are an obvious espmap input, and a real one from a DFT run must be kept.
def _mmff_raw(smiles, keep_h=True, seed=11):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=5000)
    if not keep_h:
        m = Chem.RemoveHs(m)
    X = _xyz(m)
    X -= X.mean(axis=0)
    X = X @ np.linalg.svd(X)[2].T              # best-fit plane -> xy, z kept
    return _set_xyz(m, X), float(np.abs(X[:, 2]).max())


for smi, nm in (("O=c1cc[nH]c(=O)[nH]1", "uracil"), ("Nc1cc[nH]c(=O)n1", "cytosine"),
                ("Nc1nc2[nH]cnc2c(=O)[nH]1", "guanine")):
    m, dev = _mmff_raw(smi)
    X = _xyz(m)
    X[:, 2] = 0.0                               # a planar DFT geometry is written flat
    _write(f"real_{nm}", _set_xyz(m, X), True, REAL, "geometry", "geometry",
           f"RDKit ETKDG + MMFF94s with H, rotated into xy, then z = 0 "
           f"(max |z| was {dev:.2f} A)", "a real nucleobase, planar, with H")
m, dev = _mmff_raw("Cc1c[nH]c(=O)[nH]c1=O")
_write("real_thymine", m, True, REAL, "geometry", "geometry",
       f"RDKit ETKDG + MMFF94s with H, rotated, z kept (max |z| {dev:.2f} A: methyl H)",
       "a real nucleobase whose methyl H leave the plane")
m, dev = _mmff_raw("C1=CC=CCC1", keep_h=False)
_write("real_13chd_noH", m, True, REAL, "geometry", "geometry",
       f"RDKit ETKDG + MMFF94s, H removed, rotated, z kept (max |z| {dev:.3f} A)",
       "a real ring that is nearly, not exactly, planar: within 0.1 A of its plane")

# --- M. review 6: classes the R11 rule at v0.2.1 gets wrong -------------------
# Added as cases before any rule change (review 6 S8, S15-S17), several
# variants per class so a fix cannot special-case one file.

def _write_text(name, text, flag3d, allowed, prefer, kind, source, why, ext="sdf"):
    (OUT / f"{name}.{ext}").write_text(text, encoding="utf-8")
    assert prefer is None or prefer in allowed, name
    CASES.append(dict(file=f"{name}.{ext}", kind=kind, allowed=allowed, prefer=prefer,
                      outcomes=_outcomes(kind, flag3d),
                      header="3D" if flag3d else "2D", source=source, why=why))


def _block(mol, flag3d=True):
    mol = Chem.Mol(mol)
    mol.GetConformer().Set3D(flag3d)
    return Chem.MolToMolBlock(mol)


# S15: text after the last $$$$ (a converter's trailer, a hand-added comment)
# must not shift which header line belongs to which record.
for nm, mol, flag3d, kind, allowed, prefer, tail in (
        ("real_C2F4_trailing_comment", _mmff("FC(F)=C(F)F"), True, "geometry", REAL, "geometry",
         "converted by hand from a Gaussian log\n"),
        ("real_C6F6_trailing_lines", _mmff("Fc1c(F)c(F)c(F)c(F)c1F"), True, "geometry", REAL, "geometry",
         "\n\nend of file\n"),
        ("real_C2F4_trailing_header_like", _mmff("FC(F)=C(F)F"), True, "geometry", REAL, "geometry",
         "note\n  -OEChem-10232614502D\n"),
        ("draw_chlorobenzene_2Dheader_trailing", _depiction("Clc1ccccc1", 1.5), False, "drawing", DRAW,
         "drawing", "exported by an editor\n")):
    _write_text(nm, _block(mol, flag3d) + "$$$$\n" + tail, flag3d, allowed, prefer, kind,
                "one SDF record, then text after the last $$$$",
                "trailing text after the only record: the header that counts is that record's own")

# S16: drawings with z noise up to 0.09 A. The largest |z| stays under the
# flatness tolerance, but a least-squares plane can tilt to fit the noise and
# report a larger deviation. Every seed is kept, not only those that fail.
_noise16 = np.random.default_rng(20261024)
for smi, nm, L in (("OC1CCCCC1", "cyclohexanol", 1.5), ("CC(O)CN", "aminopropanol", 1.45),
                   ("CCCCCC", "hexane", 1.5)):
    for zn in (0.08, 0.09):
        for seed in (1, 2):
            m = _depiction(smi, L)
            X = _xyz(m)
            X[:, 2] = _noise16.uniform(-zn, zn, len(X))
            _write(f"draw_{nm}_znoise_{zn}_s{seed}_{L}", _set_xyz(m, X), True, DRAW, "drawing",
                   "drawing", f"RDKit depiction at {L} A, z uniform in +-{zn} A",
                   "a flat sp3 drawing whose noise stays inside the flatness tolerance; "
                   "real puckering is ~0.2 A or more")
# review 6's own reproductions (docs/review-6/a_stereo3.py): small molecules,
# where a few atoms let the fitted plane tilt furthest
for smi, nm, seed in (("CC(N)C(=O)O", "alanine", 12), ("CC(N)C(=O)O", "alanine", 35),
                      ("CC(O)C(N)C(=O)O", "threonine", 11)):
    m = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(m)
    X = _xyz(m)
    X[:, 2] = np.random.default_rng(seed).uniform(-0.09, 0.09, len(X))
    _write(f"draw_{nm}_znoise_0.09_review6_s{seed}", _set_xyz(m, X), True, DRAW, "drawing",
           "drawing", f"RDKit depiction (1.5 A), z uniform in +-0.09 A, numpy seed {seed}",
           "review 6 S16: max |z| under the flatness tolerance, best-plane deviation over it")

# S17: real species at bond lengths and angles near the rule's thresholds.
# Literature-typical values, rounded; heavy atoms only where H are absent anyway
# or where a heavy-atom-only file is what a user would have.
_write("real_dicyanamide", _built("[N-](C#N)C#N", [(0, 0), _polar(1.31, 30), (
       _polar(1.31, 30)[0] + 1.16 * math.cos(math.radians(30 - 6)),
       _polar(1.31, 30)[1] + 1.16 * math.sin(math.radians(30 - 6))), _polar(1.31, 150), (
       _polar(1.31, 150)[0] + 1.16 * math.cos(math.radians(150 + 6)),
       _polar(1.31, 150)[1] + 1.16 * math.sin(math.radians(150 + 6)))]), True, REAL, "geometry",
       "geometry", "built: N-C 1.31 A, C#N 1.16 A, C-N-C 120 deg, N-C-N 174 deg",
       "real dicyanamide: its central N-C is 11 % short of a single-bond reference because it is "
       "delocalised, not because it is drawn")
_write("real_cyanamide_noH", _built("NC#N", [(0, 0), (1.35, 0), (1.35 + 1.16, 0)]), True, REAL,
       "geometry", "geometry", "built: N-C 1.35 A, C#N 1.16 A, linear",
       "real cyanamide's heavy atoms: a short, delocalised N-C single bond")
for ang in (144.0, 150.0):
    _write(f"real_disiloxane_noH_{int(ang)}", _built("[SiH3]O[SiH3]", [
           _polar(1.63, 90 + ang / 2), (0, 0), _polar(1.63, 90 - ang / 2)]), True, REAL,
           "geometry", "geometry", f"built: Si-O 1.63 A, Si-O-Si {ang:.0f} deg (heavy atoms only)",
           "Si-O-Si is wide in reality (about 142-150 deg); an sp3 O angle window for C-O-C "
           "does not apply to O between silicons")
_write("real_trisilylamine_noH", _built("[SiH3]N([SiH3])[SiH3]", [
       _polar(1.73, 90), (0, 0), _polar(1.73, 210), _polar(1.73, 330)]), True, REAL, "geometry",
       "geometry", "built: Si-N 1.73 A, planar N (heavy atoms only)",
       "N(SiH3)3 is planar at N in reality: a flat N between heavy atoms is not a drawing")

# Held-out probe of the S7/S8/S15-S17 fix (s16_3d_one_atom_up_0.35): a drawing
# with one atom lifted out of its plane is no longer flat, so it skipped every
# drawing test and was computed as given. Its x,y are still a drawing's: every
# bond at one drawing length when projected. Real non-flat skeletons with
# bonds near 1.5 A (heavy atoms only) guard the other side.
for smi, nm, idx, dz in (("CCOC(C)=O", "methylpropanoate", 3, 0.35),
                         ("CCOC(C)=O", "methylpropanoate", 3, 0.8),
                         ("Clc1ccccc1", "chlorobenzene", 0, 0.5),
                         ("CC(O)CN", "aminopropanol", 1, 1.0)):
    m = _depiction(smi, 1.5)
    X = _xyz(m)
    X[idx, 2] = dz
    _write(f"draw_{nm}_atom{idx}_up_{dz}_1.5", _set_xyz(m, X), True, DRAW, "drawing", "drawing",
           f"RDKit depiction at 1.5 A, atom {idx} lifted {dz} A out of the plane",
           "an edited drawing: not flat any more, but every bond projects to one drawing length")
for smi, nm in (("C1CCNCC1", "piperidine"), ("CCNCC", "diethylamine"), ("CC(C)(C)N", "tbutylamine")):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=5000)
    _write(f"real_{nm}_noH_3d", Chem.RemoveHs(m), True, REAL, "geometry", "geometry",
           "RDKit ETKDG + MMFF94s, H removed, not flattened",
           "a real non-flat skeleton whose C-C and C-N lengths sit near 1.5 A")

# S8: a 2D header is the file's own statement, also when z is not 0.
m = _depiction("C[C@@H](O)[C@H](N)C(=O)O", 1.5)
X = _xyz(m)
for i, a in enumerate(m.GetAtoms()):
    if a.GetSymbol() in ("O", "N") and a.GetDegree() == 1:
        X[i, 2] = 0.8 if i % 2 else -0.8                   # wedge ends lifted, as some editors do
_write("draw2d_threonine_wedge_z_1.5", _set_xyz(m, X), False, DRAW, "drawing", "drawing",
       "RDKit depiction at 1.5 A, 2D header, terminal O/N at z +-0.8",
       "a 2D drawing whose editor wrote wedge depth into z; the header says 2D")
for smi, nm in (("Clc1ccccc1", "chlorobenzene"), ("CC(O)CN", "aminopropanol")):
    m = _depiction(smi, 1.5)
    X = _xyz(m)
    X[0, 2] = 0.3
    _write(f"draw2d_{nm}_one_z_1.5", _set_xyz(m, X), False, DRAW, "drawing", "drawing",
           "RDKit depiction at 1.5 A, 2D header, atom 0 at z = 0.3",
           "one atom nudged out of the plane in an editor; the header still says 2D")
m = Chem.AddHs(Chem.MolFromSmiles("CC(O)CC"))
AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=2000)
_write("real_2butanol_2Dheader", m, False, DRAW, "drawing", "geometry",
       "RDKit ETKDG + MMFF94s 2-butanol saved with a 2D header",
       "the header is the file's own statement; the reader cannot overrule it (as "
       "real_chlorobenzene_2Dheader, but not flat)")

(OUT / "manifest.json").write_text(json.dumps(
    {"about": __doc__.split("\n\n")[0], "severity": SEVERITY, "cases": CASES}, indent=1) + "\n", encoding="utf-8")
print(len(CASES), "cases")
