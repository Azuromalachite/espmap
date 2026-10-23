"""Build the --require-stereo hostile input set: structures that a rule
requiring fully specified stereochemistry must sort into "pass" (every real
stereo element fixed, or none) and "block" (a real stereo element whose
configuration the input does not fix), each with what the gate should do.

Written BEFORE any rule (the R11 lesson: rules written first failed on
counterexamples found afterwards). Run from the repository root:

    python tests/data/stereo/make_set.py

It rewrites tests/data/stereo/* and manifest.json deterministically.

Every case runs the way main() runs it: parse_for_check (the SMILES, the file
or the --atoms string), then check_structure on the hydrogen-free molecule,
with require_stereo=True. The verdict is "block" when the gate returns BLOCK,
else "pass". Without the flag no case may BLOCK, so the set isolates stereo.

Each case lists
  - kind: "unspecified" (a real stereo element the input does not fix),
    "specified" (none, or all fixed), or "either" (stereo the ESP cannot see,
    or chemistry where counting it is a judgement call).
  - allowed: the verdicts that are SAFE. "pass" on an unspecified case
    computes an arbitrary isomer under a flag that promised not to, with
    nothing in the record to say so. "block" is always safe: it refuses
    with a remedy.
  - prefer: what a good rule says. A block of a fully specified structure is
    safe but costs the user a redraw. null where neither is better.
  - route: smiles, file (MOL/SDF/XYZ/PDB) or atoms (an --atoms string).

The stereo of a 3D input is the stereo of its coordinates, so every 3D file
and --atoms string is specified. A 2D drawing is specified by its wedges
and by how its double bonds are drawn, except where the file marks a bond as
"either" (crossed double bond, wavy single bond).
"""
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

OUT = Path(__file__).resolve().parent
CASES = []
SEED = 0xC0FFEE
BLOCK_ONLY = ["block"]
BOTH = ["pass", "block"]

# what each verdict does to the user, worst first
SEVERITY = {"silent arbitrary isomer": 3, "BLOCK": 1, "correct": 0}


def _outcomes(kind):
    if kind == "unspecified":
        return {"block": "correct", "pass": "silent arbitrary isomer"}
    if kind == "specified":
        return {"pass": "correct", "block": "BLOCK"}
    return {"pass": "correct", "block": "BLOCK"}       # either: no wrong isomer


def _case(name, route, inp, kind, prefer, source, why, record=None):
    allowed = BLOCK_ONLY if kind == "unspecified" else BOTH
    assert prefer is None or prefer in allowed, name
    c = dict(name=name, route=route, input=inp, kind=kind, allowed=allowed,
             prefer=prefer, outcomes=_outcomes(kind), source=source, why=why)
    if record is not None:
        c["record"] = record
    CASES.append(c)


def smi(name, s, kind, prefer, why):
    assert Chem.MolFromSmiles(s) is not None, s
    _case(name, "smiles", s, kind, prefer, "SMILES", why)


def _embed(s, seed=SEED):
    m = Chem.AddHs(Chem.MolFromSmiles(s))
    assert AllChem.EmbedMolecule(m, randomSeed=seed) == 0, s
    AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=2000)
    return m


def _round(m, nd=4):
    """Round coordinates so the files are byte-stable across platforms."""
    c = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        p = c.GetAtomPosition(i)
        c.SetAtomPosition(i, (round(p.x, nd), round(p.y, nd), round(p.z, nd)))
    return m


def _write(fname, text):
    (OUT / fname).write_text(text, encoding="utf-8", newline="\n")


def file3d(name, s, fmt, why):
    """A 3D geometry: its coordinates fix every stereo element."""
    m = _round(_embed(s))
    if fmt == "mol":
        _write(f"{name}.mol", Chem.MolToMolBlock(m))
    elif fmt == "sdf":
        _write(f"{name}.sdf", Chem.MolToMolBlock(m) + "$$$$\n")
    elif fmt == "xyz":
        _write(f"{name}.xyz", Chem.MolToXYZBlock(m))
    elif fmt == "pdb":
        _write(f"{name}.pdb", Chem.MolToPDBBlock(m))
    _case(name, "file", f"{name}.{fmt}", "specified", "pass",
          f"RDKit ETKDG + MMFF94s of {s}, as .{fmt}", why)


def atoms3d(name, s, why):
    m = _round(_embed(s))
    c = m.GetConformer()
    spec = "; ".join(f"{a.GetSymbol()} {c.GetAtomPosition(i).x:.4f} "
                     f"{c.GetAtomPosition(i).y:.4f} {c.GetAtomPosition(i).z:.4f}"
                     for i, a in enumerate(m.GetAtoms()))
    _write(f"{name}.atoms", spec + "\n")
    _case(name, "atoms", f"{name}.atoms", "specified", "pass",
          f"--atoms string from ETKDG + MMFF94s of {s}", why)


def drawing(name, s, kind, prefer, why, wedge=True, strip=False, either_bond=None,
            wavy_from=None, flag3d=False, hs=False):
    """A 2D drawing as an editor would export it.
    strip: drop all stereo before drawing (no wedges; RDKit marks unspecified
    double bonds 'either'). either_bond / wavy_from: mark one bond 'either'
    (stereo field 3 on a double bond, 4 on a single bond from an atom).
    flag3d: the 2D coordinates under a 3D header (z = 0)."""
    m = Chem.MolFromSmiles(s)
    if strip:
        Chem.RemoveStereochemistry(m)
    if hs:
        m = Chem.AddHs(m)
    AllChem.Compute2DCoords(m, canonOrient=True)
    if wedge and not flag3d:
        Chem.WedgeMolBonds(m, m.GetConformer())
    m.GetConformer().Set3D(flag3d)
    block = Chem.MolToMolBlock(m)
    if either_bond is not None or wavy_from is not None:
        lines = block.split("\n")
        n_at, n_b = int(lines[3][:3]), int(lines[3][3:6])
        for k in range(4 + n_at, 4 + n_at + n_b):
            a1, a2, order = int(lines[k][:3]), int(lines[k][3:6]), int(lines[k][6:9])
            pair = {a1 - 1, a2 - 1}
            if either_bond is not None and pair == set(either_bond) and order == 2:
                lines[k] = lines[k][:9] + "  3" + lines[k][12:]
            if wavy_from is not None and a1 - 1 == wavy_from and order == 1:
                lines[k] = lines[k][:9] + "  4" + lines[k][12:]
                wavy_from = None
        block = "\n".join(lines)
    _write(f"{name}.mol", block)
    _case(name, "file", f"{name}.mol", kind, prefer,
          f"RDKit 2D depiction of {s}" + (", stereo stripped" if strip else "")
          + (", 3D header, z = 0" if flag3d else ""), why)


S, U, E = "specified", "unspecified", "either"

# --- A. no stereo at all (controls) -------------------------------------------
smi("smi_benzene", "c1ccccc1", S, "pass", "no stereo element")
smi("smi_ethanol", "CCO", S, "pass", "no stereo element")
smi("smi_3pentanol", "CCC(O)CC", S, "pass", "two identical ethyls: not a stereocentre")
smi("smi_isobutylene", "CC(C)=C", S, "pass", "two identical methyls on one end: no E/Z")
smi("smi_propene", "CC=C", S, "pass", "a terminal =CH2: no E/Z")
smi("smi_tetramethylethylene", "CC(C)=C(C)C", S, "pass", "identical pairs on both ends: no E/Z")
smi("smi_cyclohexene", "C1CCC=CC1", S, "pass", "a double bond in a six-ring is always cis")
smi("smi_cycloheptene", "C1CCCC=CC1", S, "pass", "a seven-ring double bond is cis in practice")
smi("smi_butadiene", "C=CC=C", S, "pass", "s-cis/s-trans are conformers, not stereoisomers")
smi("smi_nmethylacetamide", "CNC(C)=O", S, "pass", "amide rotamers are conformers, not stereoisomers")
smi("smi_triethylamine_like", "CCN(C)CCC", S, "pass", "a tertiary amine inverts in microseconds: not a stereocentre")
smi("smi_cyclohexanone_4tbu_sym", "CC(C)(C)C1CCC(=O)CC1", S, "pass", "C4 sees two identical ring arms: no stereo")
smi("smi_stereo_mark_on_nothing", "C[C@H](C)O", S, "pass", "an @ on an atom with two methyls marks nothing real")

# --- B. tetrahedral centres, SMILES ---------------------------------------------
smi("smi_ala_S", "C[C@H](N)C(=O)O", S, "pass", "one centre, specified")
smi("smi_ala_unspec", "CC(N)C(=O)O", U, "block", "one centre unspecified; the enantiomers' ESP is mirrored, but the flag promises a fixed isomer")
smi("smi_ala_explicitH", "[H][C@](C)(N)C(=O)O", S, "pass", "specified with an explicit H in the SMILES")
smi("smi_threonine_full", "C[C@@H](O)[C@H](N)C(=O)O", S, "pass", "two centres, both specified")
smi("smi_threonine_half", "C[C@@H](O)C(N)C(=O)O", U, "block", "one of two centres specified: diastereomers differ")
smi("smi_threonine_none", "CC(O)C(N)C(=O)O", U, "block", "two centres unspecified: four stereoisomers, two diastereomer pairs")
smi("smi_tartaric_meso", "O[C@@H](C(=O)O)[C@H](O)C(=O)O", S, "pass", "meso tartaric acid: both centres specified")
smi("smi_tartaric_unspec", "OC(C(=O)O)C(O)C(=O)O", U, "block", "meso or (R,R)/(S,S): different compounds")
smi("smi_trihydroxyglutaric_C3_unspec", "O[C@@H](C(=O)O)C(O)[C@H](O)C(=O)O",
    U, "block", "C2 and C4 fixed meso-like; C3 is pseudo-asymmetric (r/s) and unspecified")
smi("smi_trihydroxyglutaric_full", "O[C@@H](C(=O)O)[C@@H](O)[C@H](O)C(=O)O",
    S, "pass", "every centre specified, C3 included")
smi("smi_cholesterol_full",
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
    S, "pass", "eight centres, all specified")
smi("smi_cholesterol_none", "CC(CCCC(C)C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
    U, "block", "eight centres unspecified: 256 stereoisomers")
smi("smi_glucose_open_partial", "OC[C@@H](O)[C@@H](O)C(O)[C@H](O)C=O",
    U, "block", "one of four centres left unspecified")

# --- C. ring stereo without a CIP centre, and ring junctions --------------------
smi("smi_14dimethylcyclohexane_unspec", "CC1CCC(C)CC1", U, "block", "cis and trans are different compounds, though neither carbon is a CIP centre")
smi("smi_14dimethylcyclohexane_trans", "C[C@H]1CC[C@@H](C)CC1", S, "pass", "ring stereo specified")
smi("smi_12dimethylcyclohexane_cis", "C[C@@H]1CCCC[C@@H]1C", S, "pass", "cis, meso: specified")
smi("smi_decalin_unspec", "C1CCC2CCCCC2C1", U, "block", "cis- or trans-decalin")
smi("smi_decalin_trans", "C1CC[C@H]2CCCC[C@@H]2C1", S, "pass", "ring junction specified")

# --- D. heteroatom centres --------------------------------------------------------
smi("smi_sulfoxide_unspec", "CS(=O)CC", U, "block", "a sulfoxide S is a configurationally stable stereocentre")
smi("smi_sulfoxide_spec", "C[S@](=O)CC", S, "pass", "sulfoxide specified")
smi("smi_phosphine_unspec", "CP(CC)c1ccccc1", U, "block", "a phosphine P(III) inverts only above ~30 kcal/mol: stable")
smi("smi_phosphine_spec", "C[P@](CC)c1ccccc1", S, "pass", "phosphine specified")
smi("smi_ammonium_unspec", "C[N+](CC)(CCC)CCCC", U, "block", "a quaternary N with four different groups cannot invert")
smi("smi_ammonium_spec", "C[N@+](CC)(CCC)CCCC", S, "pass", "quaternary N specified")
smi("smi_troger_unspec", "Cc1ccc2c(c1)CN1Cc3cc(C)ccc3N(C1)C2", U, "block", "Troeger's base: bridgehead N cannot invert")
smi("smi_aziridine_N", "CN1CC1(C)C", E, None,
    "the N is the only candidate; an N-alkyl aziridine inverts at roughly 17-19 kcal/mol, "
    "slowly but at room temperature: whether it counts is a judgement call")
smi("smi_aziridine_C_unspec", "CN1CC1C", U, "block", "the ring carbon is a real stereocentre, unspecified")

# --- E. double bonds ---------------------------------------------------------------
smi("smi_2butene_unspec", "CC=CC", U, "block", "E and Z are different compounds")
smi("smi_2butene_E", "C/C=C/C", S, "pass", "E specified")
smi("smi_2butene_Z", "C/C=C\\C", S, "pass", "Z specified")
smi("smi_hexadiene_half", "C/C=C/C=CC", U, "block", "one of two E/Z bonds specified")
smi("smi_hexadiene_full", "C/C=C/C=C/C", S, "pass", "both specified")
smi("smi_cyclooctene_unspec", "C1CCCC=CCC1", U, "block", "an eight-ring double bond can be E: trans-cyclooctene is isolable")
smi("smi_cyclooctene_Z", "C1CCC/C=C\\CC1", S, "pass", "cis-cyclooctene specified")
smi("smi_cyclodecene_unspec", "C1CCCCC=CCCC1", U, "block", "a ten-ring double bond: E and Z both stable")
smi("smi_oxime_unspec", "CC(=NO)CC", U, "block", "oxime E/Z isomers are separable")
smi("smi_oxime_E", "CC/C(C)=N/O", S, "pass", "oxime specified")
smi("smi_nmethylimine_unspec", "CC(=NC)CC", U, "block", "an N-alkyl ketimine has E/Z isomers")
smi("smi_azobenzene_unspec", "c1ccc(cc1)N=Nc1ccccc1", U, "block", "cis- and trans-azobenzene")
smi("smi_azobenzene_E", "c1ccc(cc1)/N=N/c1ccccc1", S, "pass", "trans-azobenzene specified")
smi("smi_NH_imine", "CC(=N)CC", E, None, "an N-H imine's E/Z interconvert quickly (tautomerism): whether it counts is a judgement call")
smi("smi_enolate_unspec", "CC=C([O-])CC", U, "block", "E and Z enolates are different species")
smi("smi_enolate_Z", "C/C=C(/[O-])CC", S, "pass", "enolate specified")

# --- F. axial stereo the graph-level tools may not see ------------------------------
smi("smi_allene_unspec", "CC=C=CC", U, "block", "1,3-dimethylallene is axially chiral: two enantiomers")
smi("smi_allene_achiral", "CC(C)=C=CC", S, "pass", "two identical methyls on one end: no axial chirality")
smi("smi_butatriene_unspec", "CC=C=C=CC", U, "block", "an even cumulene has cis/trans isomers")
smi("smi_binol", "Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12", U, "block",
    "BINOL's atropisomers are configurationally stable (barrier ~37 kcal/mol); SMILES cannot fix them")
smi("smi_biphenyl_ortho1", "Cc1ccccc1-c1ccccc1", S, "pass", "one ortho methyl: rotation is free at room temperature")
smi("smi_nitro_biphenyl", "Cc1cccc([N+](=O)[O-])c1-c1c(C)cccc1[N+](=O)[O-]", U, "block",
    "6,6'-dimethyl-2,2'-dinitrobiphenyl: the classic resolved atropisomer; SMILES cannot fix it")
smi("smi_tetramethylbiphenyl", "Cc1cccc(C)c1-c1c(C)cccc1C", S, "pass",
    "four ortho methyls, but each ring is symmetric: no atropisomers")
smi("smi_helicene6", "c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c12", U, "block",
    "[6]helicene: P and M helices, barrier ~36 kcal/mol; SMILES cannot fix the helix")
smi("smi_ethylidene_cyclohexane", "CC=C1CCC(C)CC1", U, "block",
    "1-ethylidene-4-methylcyclohexane is axially chiral, with no stereocentre and no E/Z")
smi("smi_camphor_unspec", "CC1(C)C2CCC1(C)C(=O)C2", U, "block", "camphor's two bridgeheads: (+) and (-) camphor")
smi("smi_camphor_spec", "C[C@@]12CC[C@@H](C1(C)C)CC2=O", S, "pass", "(1R,4R)-camphor specified")
smi("smi_spiro_unspec", "CC1CC2(C1)CC(C)C2", U, "block",
    "2,6-dimethylspiro[3.3]heptane is axially chiral: two enantiomers")
smi("smi_isotope_unspec", "CC([2H])O", E, "pass", "a centre made only by an isotope: the ESP cannot see it, espmap ignores isotopes")
smi("smi_isotope_spec", "C[C@H]([2H])O", E, "pass", "isotopic centre specified")

# --- G. 2D drawings (MOL files from an editor) ------------------------------------
drawing("draw_ala_wedged", "C[C@H](N)C(=O)O", S, "pass", "wedges fix the centre")
drawing("draw_ala_plain", "C[C@H](N)C(=O)O", U, "block", "no wedge: the centre is not fixed", strip=True)
drawing("draw_threonine_wedged", "C[C@@H](O)[C@H](N)C(=O)O", S, "pass", "both centres wedged")
drawing("draw_threonine_plain", "C[C@@H](O)[C@H](N)C(=O)O", U, "block", "no wedges on two centres", strip=True)
drawing("draw_2butene_E", "C/C=C/C", S, "pass", "a 2D double bond is drawn E or Z: the drawing fixes it")
drawing("draw_2butene_either", "C/C=C/C", U, "block", "a crossed ('either') double bond says E/Z is unknown",
        either_bond=(1, 2))
drawing("draw_2butene_stripped", "CC=CC", U, "block", "drawn with no E/Z: RDKit writes the bond as 'either'", strip=True)
drawing("draw_ala_wavy", "C[C@H](N)C(=O)O", U, "block", "a wavy bond at the centre says its configuration is unknown",
        strip=True, wavy_from=1)
drawing("draw_sulfoxide_plain", "C[S@](=O)CC", U, "block", "a sulfoxide drawn without a wedge", strip=True)
drawing("draw_cyclooctene_Z", "C1CCC/C=C\\CC1", S, "pass", "an eight-ring double bond drawn cis")
drawing("draw_14dmch_plain", "C[C@H]1CC[C@@H](C)CC1", U, "block", "ring stereo drawn without wedges", strip=True)
drawing("draw3d_ala_plain", "C[C@H](N)C(=O)O", U, "block",
        "a drawing under a 3D header with every z = 0 and no wedges: the flat centre fixes nothing",
        strip=True, flag3d=True)
drawing("draw_cyclohexene_crossed", "C1CCC=CC1", S, "pass",
        "a crossed bond in a six-ring: no E/Z exists, so the 'either' mark says nothing",
        either_bond=(3, 4))
drawing("draw_ethanol_wavy", "CCO", S, "pass", "a wavy bond from a CH2: no stereocentre there",
        wavy_from=1)
drawing("draw_ala_wedged_H", "C[C@H](N)C(=O)O", S, "pass",
        "explicit hydrogens drawn, the wedge on the C-H bond", hs=True)
drawing("draw3d_crotonaldehyde_E", "C/C=C/C=O", S, "pass",
        "a drawing under a 3D header: the drawn double bond still reads E. Two element "
        "pairs, so the R11 rule calls it a drawing (2-butene alone is undecidable there)",
        flag3d=True)

# --- H. 3D geometries: the coordinates fix every element ---------------------------
for s, nm in (("C[C@H](N)C(=O)O", "ala"), ("C[C@@H](O)[C@H](N)C(=O)O", "threonine"),
              ("C/C=C/C", "2butene_E"), ("C[S@](=O)CC", "sulfoxide"),
              ("C[C@H]1CC[C@@H](C)CC1", "14dmch_trans"), ("C1CCC/C=C\\CC1", "cyclooctene_Z")):
    file3d(f"geo_{nm}_mol", s, "mol", "3D coordinates fix every stereo element")
for s, nm in (("C[C@H](N)C(=O)O", "ala"), ("C/C=C/C", "2butene_E"),
              ("C[C@@H]1CCCC[C@H]1C", "12dmch_trans")):
    file3d(f"geo_{nm}_xyz", s, "xyz", "an XYZ file's coordinates fix every stereo element")
    file3d(f"geo_{nm}_pdb", s, "pdb", "a PDB file's coordinates fix every stereo element")
    atoms3d(f"atoms_{nm}", s, "--atoms coordinates fix every stereo element")
file3d("geo_cholesterol_sdf",
       "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
       "sdf", "eight centres fixed by 3D coordinates")
file3d("geo_allene_mol", "CC=C=CC", "mol",
       "a 3D allene: the coordinates fix the axis, which no SMILES can")
file3d("geo_binol_mol", "Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12", "mol",
       "a 3D BINOL: the coordinates fix the atropisomer")

# --- I. one SDF, two records: --record picks a 3D geometry or a plain drawing ------
_geo = _round(_embed("C[C@H](N)C(=O)O"))
_draw = Chem.MolFromSmiles("CC(N)C(=O)O")
AllChem.Compute2DCoords(_draw, canonOrient=True)
_write("two_records.sdf", Chem.MolToMolBlock(_geo) + "$$$$\n" + Chem.MolToMolBlock(_draw) + "$$$$\n")
_case("sdf_record1_geometry", "file", "two_records.sdf", S, "pass",
      "record 1: ETKDG + MMFF94s alanine; record 2: a plain 2D drawing",
      "--record 1 picks the 3D geometry, which fixes the centre", record=1)
_case("sdf_record2_drawing", "file", "two_records.sdf", U, "block",
      "record 1: ETKDG + MMFF94s alanine; record 2: a plain 2D drawing",
      "--record 2 picks the drawing, whose centre has no wedge", record=2)

# --- J. review 6 S7: a flat drawing under a 3D header, with z noise -----------------
# A converter or a hand edit writes small z values into an unwedged drawing.
# The chirality then read from 3D is the noise's, not the author's: two seeds
# of one drawing gave two diastereomers, both verified. Every seed is kept.
import random  # noqa: E402

_rng = random.Random(20261024)
for s, nm in (("C[C@H](N)C(=O)O", "ala"), ("C[C@@H](O)[C@H](N)C(=O)O", "threonine")):
    for zn in (0.02, 0.05):
        for seed in (1, 2, 3):
            m = Chem.MolFromSmiles(s)
            Chem.RemoveStereochemistry(m)
            AllChem.Compute2DCoords(m, canonOrient=True)
            c = m.GetConformer()
            for i in range(m.GetNumAtoms()):
                p = c.GetAtomPosition(i)
                c.SetAtomPosition(i, (p.x, p.y, _rng.uniform(-zn, zn)))
            c.Set3D(True)
            name = f"draw3d_{nm}_znoise_{zn}_s{seed}"
            _write(f"{name}.mol", Chem.MolToMolBlock(m))
            _case(name, "file", f"{name}.mol", U, "block",
                  f"RDKit 2D depiction of {s}, stereo stripped, 3D header, z uniform in +-{zn} A",
                  "an unwedged drawing with z noise: any chirality read from 3D is the noise's")

(OUT / "manifest.json").write_text(json.dumps(
    {"about": __doc__.split("\n\n")[0], "severity": SEVERITY, "cases": CASES},
    indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
print(len(CASES), "cases")
