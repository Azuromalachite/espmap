#!/usr/bin/env python3
"""
esp_map.py — electrostatic-potential maps for organic molecules,
                with a colormap that is not a dataviz disaster.

    pip install pyscf rdkit scikit-image matplotlib cmcrameri colorspacious
    python esp_map.py --selftest
    python esp_map.py --file drawn.cdxml --check          # is the drawing unambiguous?
    python esp_map.py --smiles "Clc1ccccc1" --lim 25 --out esp.png
    python esp_map.py --file drawn.cdxml --lim 25 --out esp.png
    python esp_map.py --molecule bromobenzene --cmap vik --lim 25 --out esp.png

FIXED IN v4 (adversarial review 2, see esp-adversarial-review-2.md)
    --spin built a RESTRICTED solver: RKS on O with spin=2 returned a
      converged, meaningless -74.998688 Ha against UKS -75.000759. Now UKS
      when spin != 0, with <S^2> reported, so radicals are actually computable.
    sigma-hole axis was hard-coded to a benzene ring (atoms 0..5). Axis errors
      of 25-40 deg on other halides; CC(=O)OCCCl reported -0.05 where the
      correct axis gives +1.84. Now taken from the substituent's own bond.
    .xyz carried no connectivity, so EVERY xyz file was rejected as a salt.
      Bonds are now perceived and sanitised.
    One file reader and one conformer embedder, where there were two of each:
      the xyz repair had landed in the gate's copy and not the geometry path.
    --lim on a SEQUENTIAL scale re-imposed (-lim, +lim) -- the exact symmetric
      range that branch exists to avoid. --vmin/--vmax added instead.
    The gate parsed one object while the SCF computed another; the approved
      structure is now asserted against the computed one. (Review 3, H2: that
      assertion compared the input with itself. verify_computed now checks
      the coordinates PySCF receives.)
    Fluorine no longer reports a sigma-hole, because it does not form one.

NEW IN v3 — the structure check (--check)

    A drawing can be perfectly legible to a chemist and still be ambiguous to
    a machine, and the pipeline will happily compute a flawless ESP map of a
    molecule you did not mean. Every structure is therefore gated BEFORE any
    geometry or SCF work:

      REJECTED  valence/charge errors, R-groups and query atoms, radicals
                without --spin (open shells need UKS), and disconnected fragments -- a
                salt has no single ESP surface and the counter-ion dominates
                it (sodium acetate forced through reaches +113 kcal/mol,
                against +61 for acetic acid's own proton).
      WARNED    undrawn E/Z geometry; two or more stereocentres with no
                wedges; elements outside the tested set.
      NOTED     one unspecified stereocentre (the two enantiomers have
                mirror-image ESP -- identical values, reflected picture),
                reachable tautomers, net charge, isotopes.
      --require-stereo makes every stereo element the input leaves open a
                REJECTION: centres, E/Z, 'either' marks, stereogenic
                cumulenes, hindered biaryl axes and helicenes. 3D
                coordinates fix every element.

    The stereochemistry wording is measured, not assumed: threonine's two
    diastereomers differ in mean ESP extrema by 1.6 (V_min) and 1.9 (V_max)
    kcal/mol at def2-SVP, against a spread of up to 1.25 between MMFF seeds
    of one diastereomer (bench/error_budget.jsonl; the earlier "0.9 against
    0.4" had no data behind it and did not reproduce, review 3 H8).
    Stereochemistry moves the numbers about as much as conformer choice does,
    and the SHAPE of the surface more.

NEW IN v3 — structure input from a drawing program

    --file reads .mol .sdf .cdxml .xyz .pdb, i.e. what ChemDraft, ChemDraw,
    Ketcher or MarvinSketch export. The load-bearing distinction is 2D vs 3D:
    a drawn structure exported as MOL/SDF/CDXML is a FLAT PICTURE, every z
    zero. Feeding it to an SCF would give a planar molecule and a confidently
    wrong map, so a 2D file has its coordinates discarded, a 3D conformer
    generated, and says so loudly. A file with real 3D coordinates is used
    as-is; missing hydrogens are added and relaxed with the heavy atoms
    frozen, which reproduces the all-atom result exactly (+6.61 both ways).

NEW IN v3 — SMILES input
    --smiles builds 3D coordinates with RDKit ETKDGv3 + MMFF94s. Measured
    against two reference geometries for the chlorobenzene sigma-hole:

        idealised bond lengths   +9.83 kcal/mol   (C-Cl 1.7400 A)
        SMILES + MMFF94s         +9.94            (C-Cl 1.7209 A)
        B3LYP/def2-SVP optimised +9.54

    a 0.40 kcal/mol spread, an order of magnitude below the basis-set error
    (bench/error_budget.jsonl, def2-TZVP). A force-field geometry is adequate
    for this quantity on this molecule. It is NOT for small anions: MMFF
    puts nitrite's N-O at 1.42 A against 1.25 and moves V_min/V_max by up to
    6.5 kcal/mol (bench/anion_geometry.jsonl). Use --optimize there.

    ETKDG is stochastic, so --seed is fixed by default: without it the same
    SMILES gives a different conformer every run. Conformer effect on the ESP
    extrema, lowest- vs highest-energy MMFF conformer (review 2; prose, with
    no bench/ row, unlike the C18 figures below):

        paracetamol (2 rotatable bonds)  V_max +61.13 -> +62.33
        glycerol    (4 rotatable bonds)  V_max +60.33 -> +61.26

    about 1 kcal/mol for these small molecules, and ~1.1 across 8 EXTENDED
    conformers of a C18 acetate. That does not survive folding: the two
    compact C18 conformers widened V_min to 5.2 and V_max to 2.5 kcal/mol
    and changed the negative surface area by 60 % (bench/conformers.jsonl,
    bench/f9_recheck.jsonl). Extrema are robust only while the molecule does
    not fold onto itself; the SHAPE of the surface is never robust.

    CHARGED SPECIES CHANGE THE COLORMAP PROBLEM. Benzoate's surface runs
    -158.7 to -40.7 kcal/mol -- no sign change anywhere. A diverging map
    centred on zero puts that entire molecule at one end of the ramp and
    hides 118 kcal/mol of real structure. v3 detects this and switches to a
    sequential scale over the actual range. The "symmetric about zero, always"
    rule in esp-colormaps.md is correct for NEUTRAL molecules only.

CHANGES IN v2 (all driven by the adversarial review, see esp-adversarial-review.md)
    F1  default basis raised def2-SVP -> def2-TZVP. At double zeta the
        chlorobenzene sigma-hole is wrong by +3.5 kcal/mol, which was enough
        to invert a published claim. --basis is still free, but a
        double-zeta basis now prints a warning naming the measured error.
    F3  the "density integral ~= nelec" check is NOT a validation of this
        pipeline and no longer pretends to be: it varies -5%..+7% with grid
        spacing while V_max moves 0.5%. Reported as a diagnostic, never gating.
    F4  the isosurface is now checked to be CLOSED (Euler characteristic 2)
        and strictly interior to the grid box. A clipped surface silently
        changed V_max by 26% in the review. Hard failure by default.
    F5  int1e_grids chunk is sized from nao to a memory budget; the old fixed
        chunk asked for 9.7 GB at 900 basis functions.
    F6  --lim no longer silently auto-scales, which is the one thing the
        accompanying document tells you never to do.
    F7  --selftest exercises the ESP sign and magnitude against known values.
        The v1 sign error was found by hand; now it is an assertion.
    F8  per-vertex normal orientation replaces a global flip whose test
        statistic was numerical noise.

WHAT IT DOES
    1. DFT single point (PySCF)
    2. electron density on a regular grid
    3. marching cubes at rho = 0.002 e/bohr^3 (Bader et al. 1987; encloses
       98.5-99.1 % of the electrons here, bench/isosurface_enclosed.jsonl)
    4. ESP evaluated ONLY at the isosurface vertices, not on a full grid
       -- ~10^4 points instead of ~10^6, so this costs seconds, and the
          potential is evaluated exactly where it is displayed rather than
          interpolated onto the surface from a cube file. An established
          approach (Bulat et al. 2010; Lu & Chen 2012), see docs/prior-art.md
    5. z-buffered render with deliberately weak shading, so the colormap's
       lightness carries the data instead of fighting the lighting

COLOR CONVENTION
    Every colormap is oriented warm/red = NEGATIVE V = electron-rich, the
    entrenched chemistry reading. Fixing perception should not silently flip
    the sign convention too.

WHAT THIS IS AND IS NOT
    An instrument for making legible pictures and for order-of-magnitude
    reasoning. Quantitative V_max values carry, at minimum, a basis-set
    uncertainty; gas phase unless --solvent is given; geometries are idealised unless you
    supply your own. Do not quote a number from it against a margin smaller
    than the convergence tests reported in the review.
"""
import argparse
import os
import sys

import numpy as np

# the one place the version is written: pyproject.toml reads it from here
# (review 5 R18), and provenance.espmap_version is this module's own
__version__ = "0.2.1"

BOHR = 0.52917721092
HA2KCAL = 627.5095
AU2DEBYE = 2.541746473
# esp_at_points sums the solute's nuclei and (PCM-polarised) electrons only;
# the potential of the PCM surface charges is not added (review 3 H7).
SOLVENT_ESP = {False: "polarised solute only (reaction field excluded)",
               True: "polarised solute + PCM reaction field"}
ISO_BADER = 0.002

# basis sets at or below double zeta; penalty measured in the review
def basis_key(name):
    """One spelling per basis set. PySCF resolves def2-TZVP, def2tzvp and
    DEF2_TZVP to the same basis; review 4 N9 found them giving three
    analysis_keys, and --basis def2svp skipping the double-zeta warning.
    Review 5 R13: Pople polarisation is spelled two ways, 6-31G(d) = 6-31G*
    and 6-31G(d,p) = 6-31G** (same nao), so those are folded too."""
    k = "".join(c for c in name.lower() if c not in "-_ ")
    return k.replace("(d,p)", "**").replace("(d)", "*")


# Review 5 R13: the list stopped at six names, so cc-pVDZ, 6-31G** or
# pcseg-1 ran without the double-zeta warning or the basis budget row.
# tests/test_fast.py checks every name here resolves in PySCF (def2-SV(P),
# listed before, does not, so no run could carry it).
SMALL_BASIS_NAMES = ("sto-3g", "3-21g", "6-31g", "6-31g*", "6-31g**", "6-31+g*",
                     "6-31+g**", "6-31++g**", "def2-svp", "def2-svpd", "cc-pvdz",
                     "aug-cc-pvdz", "pcseg-1", "aug-pcseg-1")
SMALL_BASIS = {basis_key(b) for b in SMALL_BASIS_NAMES}

# dielectric constants for --solvent. Measured effect on the chlorobenzene
# sigma-hole: gas +9.83 -> chloroform +7.14 -> water +6.22 kcal/mol, i.e. the
# solvent axis is CO-DOMINANT with the basis-set axis (+3.54). Do not quote a
# condensed-phase number from a gas-phase run. That is the polarised-solute
# convention; with the reaction field added the shifts are +0.05 (water) and
# -0.34 (chloroform), bench/solvent_convention.jsonl (review 3 H7).
SOLVENTS = {"vacuum": None, "water": 78.3553, "dmso": 46.826, "methanol": 32.613,
            "acetonitrile": 35.688, "chloroform": 4.7113, "toluene": 2.3741,
            "hexane": 1.8819}

# nao above which an SCF stops being an interactive operation (measured:
# 253 nao -> 165 s, 298 -> 224 s on one core)
NAO_WARN, NAO_REFUSE = 450, 1200


# ------------------------------------------------------------------- geometry
def from_smiles(smiles, n_confs=20, seed=0xC0FFEE, max_iters=2000):
    """SMILES -> (PySCF atom string, info dict) via ETKDGv3 + MMFF94s.

    seed is REQUIRED in effect: ETKDG is stochastic and an unseeded run makes
    the whole pipeline non-deterministic for the same input.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    # one embedding implementation, shared with from_file: v3 had two copies
    best, ff, spread, n_emb, n_conv, _ = _embed_and_minimise(mol, n_confs, seed,
                                                             max_iters=max_iters)
    info = dict(forcefield=ff, n_embedded=n_emb, n_converged=n_conv,
                e_spread=spread,
                charge=Chem.GetFormalCharge(mol),
                n_radical=int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())),
                formula=Chem.rdMolDescriptors.CalcMolFormula(mol),
                n_rot=int(AllChem.CalcNumRotatableBonds(mol)), seed=seed,
                source="3D generated from SMILES by ETKDGv3", n_found=1,
                smiles=smiles, mol=mol, h_source=_h_source(mol, 0, True))
    return _mol_to_spec(mol, best), info


def _h_source(mol, n_supplied, relaxed):
    """Per atom: None for heavy atoms; for hydrogens 'supplied' (index below
    n_supplied: AddHs appends), else 'generated_relaxed' or 'generated'. The
    site block flags a polar H that was not supplied (spec section 4)."""
    return [None if a.GetAtomicNum() > 1 else
            "supplied" if a.GetIdx() < n_supplied else
            "generated_relaxed" if relaxed else "generated"
            for a in mol.GetAtoms()]


def conformer_note(n_rot, e_spread):
    """The warning printed for a generated conformer of a flexible molecule.

    A function so tests/test_claims.py can check its numbers against the
    measurements in bench/: this text once kept a claim (~1 kcal/mol, from
    glycerol) for a day after the C18 data had disproved it.
    """
    return (f"note: {n_rot} rotatable bond(s), force-field conformer spread "
            f"{e_spread:.1f} kcal/mol. This map is ONE conformer (the "
            f"lowest-energy one). Measured: among extended conformers the ESP "
            f"extrema move ~1 kcal/mol, but a molecule that can fold onto "
            f"itself can move them much more (5.2 kcal/mol on V_min for a C18 "
            f"acetate) and the surface shape changes either way. For a long "
            f"flexible molecule, sample shape-diverse conformers before quoting "
            f"a number.")


def _mol_to_spec(mol, conf_id=-1):
    conf = mol.GetConformer(conf_id)
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    return "; ".join(
        f"{s} {conf.GetAtomPosition(i).x:.6f} {conf.GetAtomPosition(i).y:.6f} "
        f"{conf.GetAtomPosition(i).z:.6f}" for i, s in enumerate(syms))


def _embed_and_minimise(mol, n_confs, seed, max_iters=2000):
    """ETKDGv3 + MMFF94s over n_confs, embedded into mol. Returns (best_conf_id,
    forcefield, spread, n_embedded, n_converged, converged), where converged
    is [(conf_id, energy)] for every conformer the force field converged, in
    embedding order: a single map uses the best one, an ensemble
    (select_ensemble) chooses among all of them. The only embedder
    (invariant 1)."""
    from rdkit.Chem import AllChem
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.useSmallRingTorsions = True
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=p)
    if len(cids) == 0:
        raise RuntimeError("RDKit could not embed any conformer")
    if AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s") is None:
        res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=max_iters)
        ff = "UFF (MMFF94s has no parameters for this molecule)"
    else:
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s",
                                                maxIters=max_iters)
        ff = "MMFF94s"
    ok = np.array([r[0] for r in res])
    en = np.where(ok == 0, np.array([r[1] for r in res], float), np.inf)
    if not np.isfinite(en).any():
        raise RuntimeError("no conformer converged under the force field")
    fin = en[np.isfinite(en)]
    converged = [(int(c), float(e)) for c, e in zip(cids, en) if np.isfinite(e)]
    return (int(cids[int(np.argmin(en))]), ff, float(fin.max() - fin.min()),
            len(cids), int((ok == 0).sum()), converged)


def coordinate_kind(mol, fmt, tol=0.03):
    """What the coordinates in a file are: "none", "geometry", "drawing", or
    "ambiguous" (which the reader turns into a BLOCK).

    Review 3 H5: the old rule ("every z is zero" = drawing) also discarded
    real planar geometries -- a DFT-optimised chlorobenzene in the xy plane,
    water from an .xyz -- and so depended on how many decimals a file printed.
      - .xyz and .pdb are never drawings;
      - MOL/SDF/CDXML whose header says 2D are drawings. The header is read
        from the file (_header_dim), not from RDKit's flag, which turns a
        3D header into 2D when every z is 0 and a bond carries a stereo mark
        (real C2F4 with an "either" double bond was discarded as a drawing);
      - a 3D file whose atoms all lie within FLAT_TOL of one plane is flat,
        wherever the plane is: a drawing tilted, moved to z = 5 or given z
        noise is still a drawing. A non-flat file is a geometry. A flat one:
      - an sp3 atom lying flat where no real molecule can put it
        (flat_sp3_problems) makes it a DRAWING; a flat sp3 four-membered
        ring (real cyclobutane puckers, oxetane nearly does not) makes it
        AMBIGUOUS;
      - every bond the same length (within tol) at a drawing program's
        standard (DRAWING_BOND_LENGTHS) makes it a DRAWING with two or more
        element pairs. With one pair a symmetric molecule could look the
        same, so it is AMBIGUOUS unless the length fits none of the bonds
        drawn: benzene and CO2 at 1.5 A are DRAWINGS (aromatic C-C +7.3 %,
        C=O +21 %), nitrate is AMBIGUOUS (its N-O single bond fits 1.5 A,
        +9.5 %). An acyclic two-coordinate sp3 atom at a
        real angle is evidence against a drawing and skips this test
        (anti-butane's carbon skeleton, 1.53 A, 112 deg);
      - otherwise it is a GEOMETRY if every bond is within tolerance of its
        reference length (flat_length_problems);
      - if not, a DRAWING when its bonds are uniform across two or more
        element pairs (a drawing at an unusual scale) and some bond is
        DRAWN_BOND_FACTOR tolerances off, else AMBIGUOUS: a real H-free
        isoxazole (N-O just over) is refused, not discarded.
    Review 5 R11: master discarded real C6F6 and C2F4, and computed a
    stretched drawing and a z-jittered one flat. Two rules written before a
    test set each opened a silent hole (a chlorobenzene drawing at 1.4 A;
    flat cyclohexanol and THF), and the first rule written against the
    frozen tests/data/r11 set opened more off it (squashed chains, tilted
    drawings, flat cyclobutane), found by an independent probe.
    Residual, by construction: a drawing within tolerance of a real
    geometry is computed as one (hexane drawn at 1.45 A with y stretched
    15 % has C-C 1.51 A and 113 deg).
    """
    import numpy as np
    if mol.GetNumConformers() == 0:
        return "none"
    if fmt in ("xyz", "pdb"):
        return "geometry"
    c = mol.GetConformer()
    header = mol.GetProp("_header_dim") if mol.HasProp("_header_dim") else None
    # Review 6 S8: a 2D header is the file's own statement, also when some z
    # is not 0 (an editor writing wedge depth, a nudged atom).
    if header == "2D" or not (c.Is3D() or header == "3D"):
        return "drawing"
    X = np.array([list(c.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    X -= X.mean(axis=0)
    # Flat when one plane holds every atom within FLAT_TOL. The least-squares
    # plane can tilt to fit z noise and leave a larger maximum deviation than
    # the z = const plane does (review 6 S16: a drawing with |z| <= 0.08 A
    # was 0.101 A off its fitted plane); both are planes, so take the smaller.
    # The best z = const plane sits at the midrange: half the z span.
    dev = (min(np.abs(X @ np.linalg.svd(X)[2][-1]).max(), np.ptp(X[:, 2]) / 2)
           if len(X) > 3 else 0.0)
    if dev > FLAT_TOL:
        # An edited drawing (one atom lifted out of the page) is not flat, but
        # its bonds still project to one drawing-program length; no real
        # non-flat molecule's do (held-out probe of the review 6 fix).
        return "drawing" if _projects_to_a_drawing(mol, X) else "geometry"
    if mol.GetNumBonds() == 0:
        return "geometry"
    sp3_problems, real_angles, ring4 = flat_sp3_problems(mol)
    if sp3_problems:
        # a nearly flat real ring (H-free 1,3-cyclohexadiene, 0.09 A off its
        # plane) is refused, not discarded; a drawing is flat to ~0
        rings_only = all(p.startswith("flat ") for p in sp3_problems)
        return "ambiguous" if rings_only and dev > RING_OFF_PLANE else "drawing"
    if ring4:
        return "ambiguous"
    pairs, lengths = set(), []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        pairs.add(tuple(sorted((b.GetBeginAtom().GetSymbol(),
                                b.GetEndAtom().GetSymbol()))))
        lengths.append((c.GetAtomPosition(i) - c.GetAtomPosition(j)).Length())
    mean = sum(lengths) / len(lengths)
    uniform = max(abs(x - mean) for x in lengths) <= tol * mean
    at_drawing_length = any(abs(mean - d) <= tol * d for d in DRAWING_BOND_LENGTHS)
    if uniform and at_drawing_length and not real_angles:
        if len(pairs) > 1:
            return "drawing"
        # one pair: a symmetric molecule could look the same, unless the
        # length fits none of the bonds drawn (benzene, CO2 at 1.5 A; nitrate's
        # N-O single bond fits 1.5 A, so it stays undecided)
        if len(flat_length_problems(mol, lengths, aromatic_cc=True)) == len(lengths):
            return "drawing"
        return "ambiguous"
    if not flat_length_problems(mol, lengths):
        return "geometry"
    if uniform and len(pairs) > 1 and flat_length_problems(mol, lengths, DRAWN_BOND_FACTOR):
        return "drawing"
    return "ambiguous"


def _projects_to_a_drawing(mol, X):
    """True when a non-flat file's bonds, projected onto its best-fit plane or
    onto z = 0, are all within PROJECTED_TOL of one drawing-program length
    (DRAWING_BOND_LENGTHS) across two or more element pairs: a drawing with an
    atom lifted out of the page keeps its in-plane lengths. A real non-flat
    heavy-atom skeleton (piperidine, diethylamine: C-C 1.53, C-N 1.47) projects
    to lengths that vary far more than this."""
    import numpy as np
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    pairs = {tuple(sorted((mol.GetAtomWithIdx(i).GetSymbol(), mol.GetAtomWithIdx(j).GetSymbol())))
             for i, j in bonds}
    if len(bonds) < 2 or len(pairs) < 2:
        return False
    D = np.array([X[j] - X[i] for i, j in bonds])
    normal = np.linalg.svd(X)[2][-1]
    for L in (np.sqrt(np.maximum((D ** 2).sum(1) - (D @ normal) ** 2, 0.0)),
              np.linalg.norm(D[:, :2], axis=1)):
        if any(np.all(np.abs(L / d - 1) <= PROJECTED_TOL) for d in DRAWING_BOND_LENGTHS):
            return True
    return False


# A non-flat file is still a drawing when every bond projects to within this
# fraction of one drawing length (_projects_to_a_drawing). A drawing program
# writes them equal to the last decimal; real skeletons differ by 4 % or more
# (C-C 1.53 vs C-N 1.47) before projection shortens them unevenly.
PROJECTED_TOL = 0.01
# A file is flat when every atom lies within this of its best-fit plane (A).
# A real non-planar molecule leaves it by far more (a CH2's hydrogens by
# ~0.9, ammonia's N by ~0.3); a converter's z noise stays well inside it.
FLAT_TOL = 0.1
# A flat file whose only sp3 problem is a ring is refused rather than called
# a drawing when it sits further than this off its plane (A).
RING_OFF_PLANE = 0.01
# Standard bond lengths drawing programs write (A): RDKit and most editors
# 1.5, Marvin 1.54, ChemDraw 0.825 (its 14.4 pt default in MOL units).
DRAWING_BOND_LENGTHS = (1.5, 1.54, 0.825)
# A bond in a flat file may sit this far from its reference length and still
# be real. Measured on the R11 set's real geometries: worst aromatic +4.0
# (furan C-O); worst non-aromatic -9.8 (enolate C-O), -9.2 carbonate, -8.8
# nitrate, where one Kekule form is compared with a delocalised bond. Ozone,
# +11.3, falls outside and is refused rather than computed. A chlorobenzene
# drawing at 1.4 A has C-Cl -21 %.
PLAUSIBLE_TOL = 0.10
# A bond to an sp3 carbon or a localised C=C has no Kekule ambiguity, and
# real ones sit within 4 % of the reference (C-C 1.50-1.54, C-O 1.43, C-N
# 1.47, C-F up to 1.38, C-Cl 1.78; C=C 1.33, the enolate's 1.37 +3.6 %). A
# squashed hexane drawn at 1.5 A has C-C +5.3 %; cyclopentadiene drawn at
# 1.45 A has C=C +9.6 %. Aromatic C-C joins only in the one-pair test (all
# carbon; worst real +2.3 %, furan; benzene drawn at 1.5 A is +7.3 %).
TIGHT_BOND_TOL = 0.05
# single, aromatic, double, triple: fraction of the single-bond length
_ORDER_FACTOR = {1.0: 1.0, 1.5: 0.92, 2.0: 0.87, 3.0: 0.79}
# An implausible flat file is called a drawing only when some bond is this
# many tolerances off; nearer, it is refused. Measured: every drawing that
# reaches that step is 1.58x or more (phenylacetate drawn at 1.4 A), real
# H-free isoxazole 1.00x (aromatic N-O +10.0 %).
DRAWN_BOND_FACTOR = 1.3
# Angles a real acyclic two-coordinate sp3 atom takes (deg): H2S 92, water
# 104.5, chains ~112-114. A drawing's zigzag is 120.
SP3_REAL_ANGLE = (90.0, 118.0)


def bond_reference_length(bond):
    """A rough expected length (A) for an RDKit bond: the sum of the two
    atoms' covalent radii (RDKit's table, Cordero et al. 2008) scaled for
    bond order. Used only to ask whether a flat file could be a real
    geometry (coordinate_kind)."""
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    r = (pt.GetRcovalent(bond.GetBeginAtom().GetAtomicNum())
         + pt.GetRcovalent(bond.GetEndAtom().GetAtomicNum()))
    return r * _ORDER_FACTOR.get(bond.GetBondTypeAsDouble(), 1.0)


def flat_length_problems(mol, lengths, factor=1.0, aromatic_cc=False):
    """Bonds of a flat file further from their reference length than a real
    bond goes (TIGHT_BOND_TOL for a bond to an sp3 carbon or a C=C, else
    PLAUSIBLE_TOL, each times factor); empty if every bond could be real.

    aromatic_cc puts aromatic C-C in the tight class too. Only the one-pair
    test asks for it, where the file is all carbon: in a heteroaromatic
    lactam RDKit's aromatic C-C is localised (uracil's C4-C5 1.477 A, +5.6 %,
    and C5=C6 1.335 A), and a 5 % class refused real uracil and cytosine."""
    from rdkit import Chem
    sp3 = Chem.HybridizationType.SP3
    out = []
    for b, L in zip(mol.GetBonds(), lengths):
        ref = bond_reference_length(b)
        ends = (b.GetBeginAtom(), b.GetEndAtom())
        tight = (any(a.GetSymbol() == "C" and a.GetHybridization() == sp3 for a in ends)
                 or (b.GetBondTypeAsDouble() in ((1.5, 2.0) if aromatic_cc else (2.0,))
                     and ends[0].GetSymbol() == ends[1].GetSymbol() == "C"))
        tol = (TIGHT_BOND_TOL if tight else PLAUSIBLE_TOL) * factor
        if abs(L / ref - 1) > tol:
            out.append(f"{ends[0].GetSymbol()}-{ends[1].GetSymbol()} {L:.3f} A vs {ref:.2f}")
    return out


def flat_sp3_problems(mol):
    """(problems, real_angles, ring4) for the sp3 atoms of a flat file.

    problems name sp3 atoms that no real molecule puts in a plane:
    - three or more neighbours: a tetrahedral centre is never flat
      (t-butanol, 1-amino-2-propanol's C2, cyclohexanol's C1);
    - in a ring of five to eight with a bond between two sp3 ring atoms:
      the eclipsed bond puckers it (cyclopentane, THF, cyclohexane). A lone
      sp3 atom between sp2 ones stays flat (fluorene's and
      cyclopentadiene's CH2, 1,4-cyclohexadiene);
    - acyclic, two neighbours, outside SP3_REAL_ANGLE: a drawn hexane at
      120 deg, or one squashed to 130.
    real_angles counts acyclic two-coordinate sp3 atoms inside it: evidence
    the coordinates were computed, not drawn (anti-butane's skeleton,
    water). Ring angles are fixed by the ring in a drawing too, so they
    count for neither. ring4 lists flat four-membered rings with such a
    bond: cyclobutane puckers but oxetane barely does, so they are refused.
    Hybridisation is RDKit's, which calls conjugated N and O (amides,
    anilines, pyrrole, furan) sp2."""
    import math
    from rdkit import Chem
    sp3 = Chem.HybridizationType.SP3
    out, real, ring4 = [], 0, []
    c = mol.GetConformer()
    for ring in mol.GetRingInfo().BondRings():
        if any(mol.GetBondWithIdx(k).GetBeginAtom().GetHybridization() == sp3
               and mol.GetBondWithIdx(k).GetEndAtom().GetHybridization() == sp3
               for k in ring):
            if len(ring) == 4:
                ring4.append(ring)
            elif 5 <= len(ring) <= 8:
                out.append(f"flat {len(ring)}-ring with an sp3-sp3 bond")
    for a in mol.GetAtoms():
        if a.GetHybridization() != sp3:
            continue
        # Review 6 S17: O and N bonded to silicon are wide or planar in
        # reality (Si-O-Si 142-150 deg, N(SiH3)3 planar), so their angle and
        # flatness say nothing about a drawing; the bond lengths still decide.
        if a.GetSymbol() in ("O", "N") and any(n.GetSymbol() == "Si" for n in a.GetNeighbors()):
            continue
        nb = [n.GetIdx() for n in a.GetNeighbors()]
        tag = f"sp3 {a.GetSymbol()}{a.GetIdx()}"
        if len(nb) >= 3:
            out.append(f"{tag} with {len(nb)} neighbours lies flat")
        elif len(nb) == 2 and not a.IsInRing():
            p, q, r = (c.GetAtomPosition(i) for i in (nb[0], a.GetIdx(), nb[1]))
            u, v = p - q, r - q
            ang = math.degrees(math.acos(max(-1.0, min(1.0, u.DotProduct(v)
                                                       / (u.Length() * v.Length())))))
            if SP3_REAL_ANGLE[0] <= ang <= SP3_REAL_ANGLE[1]:
                real += 1
            else:
                out.append(f"{tag} bent {ang:.0f} deg")
    return out, real, ring4


def from_file(path, n_confs=20, seed=0xC0FFEE, charge=0, record=None):
    """Read .mol/.sdf/.cdxml/.xyz/.pdb and return (PySCF atom string, info).

    The load-bearing distinction is a DRAWING vs a GEOMETRY. A structure drawn
    in ChemDraft, ChemDraw, Ketcher or MarvinSketch and exported as
    MOL/SDF/CDXML is a picture: feeding it to an SCF would produce a confidently
    wrong map, so its coordinates are discarded and a conformer generated.
    read_structure decides which (coordinate_kind) and blocks when it cannot.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    # ONE reader, shared with the gate. v3 had a second copy here, so a fix
    # applied to parse_for_check (XYZ bond perception) silently did not reach
    # the path that actually builds the geometry.
    mol, problems, n_found = read_structure(path, charge=charge, record=record)
    if mol is None:
        raise ValueError(problems[0][2] if problems else
                         f"RDKit could not read any molecule from {path}")

    has_conf = mol.GetNumConformers() > 0
    drawing = mol.GetProp("_coords") == "drawing"
    # how the file's coordinates were made, if it says (review 5 R1); read
    # before AddHs, which does not carry SDF data fields over
    declared = next((mol.GetProp(k) for k in ("espmap_geometry_origin", "geometry")
                     if mol.HasProp(k)), None)
    n_h_before = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")
    perceived = mol.HasProp("_bonds") and mol.GetProp("_bonds") == "perceived"

    if has_conf and not drawing:
        n_supplied = mol.GetNumAtoms()
        mol = Chem.AddHs(mol, addCoords=True)
        n_h_added = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H") - n_h_before
        # Independent of the reader: if bond orders were INFERRED from the
        # coordinates and the file already had hydrogens, adding more means
        # the inferred bonds disagree with the file (review 3 H1: 5 H in the
        # file, 11 after AddHs). Files that carry their own bond orders
        # (MOL/SDF/CDXML) may legitimately omit some H, e.g. a docking pose
        # that kept only its polar hydrogens, so they are not refused here.
        if perceived and n_h_before and n_h_added:
            raise ValueError(
                f"H_COUNT_MISMATCH: the file has {n_h_before} hydrogen(s) but the "
                f"bonds perceived from its coordinates need {n_h_added} more. The "
                f"perceived structure is not the file's molecule; supply MOL/SDF "
                f"with explicit bond orders.")
        source = "file 3D coordinates used as-is"
        ff, spread, n_emb, n_conv = "none", 0.0, 1, 1
        relaxed = False
        if n_h_added:
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            if props is not None:
                field = AllChem.MMFFGetMoleculeForceField(mol, props)
                # AddHs appends, so indices below n_supplied are the file's
                # atoms. Freeze all of them, supplied H included: a pose's
                # O-H orientation was set by its environment, and relaxing it
                # here would replace that with the isolated-molecule optimum.
                for i in range(n_supplied):
                    field.AddFixedPoint(i)
                field.Minimize(maxIts=500)
                ff = "MMFF94s, added hydrogens only (supplied atoms frozen)"
                relaxed = True
            source += f"; {n_h_added} hydrogen(s) added"
        best = -1
        h_source = _h_source(mol, n_supplied, relaxed)
    else:
        mol = Chem.AddHs(mol)
        best, ff, spread, n_emb, n_conv, _ = _embed_and_minimise(mol, n_confs, seed)
        h_source = _h_source(mol, 0, True)
        source = ("2D depiction: drawn coordinates DISCARDED, 3D generated "
                  "by ETKDGv3" if drawing
                  else "no coordinates in file, 3D generated by ETKDGv3")

    info = dict(forcefield=ff, n_embedded=n_emb, n_converged=n_conv,
                e_spread=spread, charge=Chem.GetFormalCharge(mol),
                n_radical=int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())),
                formula=Chem.rdMolDescriptors.CalcMolFormula(mol),
                n_rot=int(AllChem.CalcNumRotatableBonds(mol)), seed=seed,
                source=source, n_found=n_found,
                record=dict(index=record or 1, of=n_found),
                smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)), mol=mol,
                h_source=h_source,
                declared_geometry=None if drawing or not has_conf else declared)
    return _mol_to_spec(mol, best), info


# ------------------------------------------------------- admissibility gate
# Elements this pipeline has actually been exercised on. Wider is not refused,
# only flagged: def2 covers H-Rn, but closed-shell RKS on a transition metal is
# usually wrong for reasons that have nothing to do with the basis.
TESTED_ELEMENTS = set("H B C N O F Si P S Cl Br I Se".split())

BLOCK, WARN, INFO = "BLOCK", "WARN", "INFO"


def atom_map_findings(mol, smiles=None):
    """Atom maps are all or nothing on heavy atoms, unique, and never on a
    hydrogen (per-site-spec section 1, v1). A site table keyed by maps with
    holes or repeats is worse than no maps. mol must still carry its
    hydrogens (RemoveHs drops mapped ones), and for SMILES the string is
    re-read with hydrogens kept: RDKit's parser merges [H:4] into an implicit
    hydrogen and loses the map without a word."""
    from rdkit import Chem
    out = []
    atoms = list(mol.GetAtoms())
    if smiles:
        ps = Chem.SmilesParserParams()
        ps.removeHs = False
        raw = Chem.MolFromSmiles(smiles, ps)
        if raw is not None:
            atoms = list(raw.GetAtoms())
    h_maps = [a.GetAtomMapNum() for a in atoms if a.GetAtomicNum() == 1 and a.GetAtomMapNum()]
    heavy = [a.GetAtomMapNum() for a in atoms if a.GetAtomicNum() > 1]
    mapped = [m for m in heavy if m]
    if h_maps:
        out.append((BLOCK, "MAPPED_HYDROGEN",
                    f"atom map(s) {h_maps} on hydrogen; in this version a hydrogen "
                    f"belongs to its heavy atom's site and has no identity of its own",
                    "remove the map numbers from the hydrogens"))
    if mapped and len(mapped) < len(heavy):
        out.append((BLOCK, "PARTIAL_ATOM_MAP",
                    f"{len(mapped)} of {len(heavy)} heavy atoms carry an atom map",
                    "map every heavy atom or none"))
    dup = sorted({m for m in mapped if mapped.count(m) > 1})
    if dup:
        out.append((BLOCK, "DUPLICATE_ATOM_MAP", f"atom map(s) {dup} used more than once",
                    "renumber the atom maps so each is unique"))
    return out


# --require-stereo. A biaryl axis is treated as a stable atropisomer when its
# ortho positions, across both rings, score at least this: 1 per substituent,
# 2 per peri ring fusion (a naphthalene C8a: its peri C-H crowds the axis
# like a substituent, on a rigid ring). Tri- and tetra-ortho biaryls, BINOL
# and 1,1'-binaphthyl qualify; 2,2'-dimethylbiphenyl rotates at room
# temperature and does not. This approximates, from the graph alone, the
# barrier classes of LaPlante et al. 2011 (docs/prior-art.md): axes below
# about 20 kcal/mol interconvert within a minute. Borderline axes racemise
# over hours, and a refusal is the safe side.
ATROPISOMER_MIN_ORTHO = 3
# A helicene is treated as configurationally stable from five ortho-fused
# rings; [5]helicene, like 1,1'-binaphthyl, racemises over hours at room
# temperature, [6] does not.
HELICENE_MIN_RINGS = 5


def coordinates_fix_stereo(mol):
    """True when mol carries a real 3D geometry, whose coordinates fix every
    stereo element: a 3D MOL/SDF that coordinate_kind did not call a drawing,
    an XYZ or PDB file, an --atoms string. A SMILES has no conformer; a 2D
    drawing, or one under a 3D header with z = 0, fixes only what its wedges
    and drawn double bonds say."""
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        return False
    return not (mol.HasProp("_coords") and mol.GetProp("_coords") == "drawing")


def atropisomer_axes(mol):
    """Biaryl single bonds that are stable stereogenic axes: not in a ring,
    between two aromatic atoms, with ATROPISOMER_MIN_ORTHO or more substituted
    ortho positions in all, and each end unsymmetric (its two ortho atoms not
    equivalent; four ortho methyls on two symmetric rings are not chiral).
    RDKit's stereo perception does not see these, and no SMILES can fix one."""
    from rdkit import Chem
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    ri = mol.GetRingInfo()
    out = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing():
            continue
        ends = (b.GetBeginAtom(), b.GetEndAtom())
        if not all(a.GetIsAromatic() and a.GetDegree() == 3 for a in ends):
            continue
        n_sub, chiral = 0, True
        for a, partner in (ends, ends[::-1]):
            ortho = [n for n in a.GetNeighbors() if n.GetIdx() != partner.GetIdx()]
            if ranks[ortho[0].GetIdx()] == ranks[ortho[1].GetIdx()]:
                chiral = False
            for o in ortho:
                ring = next((set(r) for r in ri.AtomRings()
                             if a.GetIdx() in r and o.GetIdx() in r), set())
                if any(n.GetAtomicNum() > 1 and n.GetIdx() != a.GetIdx()
                       and n.GetIdx() not in ring for n in o.GetNeighbors()):
                    n_sub += 2 if ri.NumAtomRings(o.GetIdx()) > 1 else 1
        if chiral and n_sub >= ATROPISOMER_MIN_ORTHO:
            out.append(b.GetIdx())
    return out


# Rings of up to 8 atoms: trans-cyclooctene's helicity is configurationally
# stable; trans-cyclononene racemises at room temperature, below the ~20
# kcal/mol barrier that ATROPISOMER_MIN_ORTHO's class follows (laplante2011).
# Below 8 atoms RDKit does not perceive ring double-bond stereo, so a trans
# mark there is dropped, not flagged; such a molecule is not isolable anyway.
SMALL_TRANS_RING = 8


def _ring_neighbours(mol, b):
    """For a ring bond a=b, the ring atoms next to each end in the bond's
    smallest ring: the C-C=C-C torsion inside the ring."""
    ri = mol.GetRingInfo()
    ring = min((r for r in ri.BondRings() if b.GetIdx() in r), key=len)
    atoms = {mol.GetBondWithIdx(k).GetBeginAtomIdx() for k in ring} | \
        {mol.GetBondWithIdx(k).GetEndAtomIdx() for k in ring}
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    ni = next(n.GetIdx() for n in b.GetBeginAtom().GetNeighbors()
              if n.GetIdx() != j and n.GetIdx() in atoms)
    nj = next(n.GetIdx() for n in b.GetEndAtom().GetNeighbors()
              if n.GetIdx() != i and n.GetIdx() in atoms)
    return len(ring), (ni, i, j, nj)


def small_ring_trans_bonds(mol):
    """Double bonds marked E/trans whose smallest ring has SMALL_TRANS_RING or
    fewer atoms (review 6 S1). The E label fixes the bond; the ring then
    twists one way or the other around it (trans-cyclooctene's planar
    chirality), and neither RDKit's stereo perception nor a SMILES sees that.
    A trans double bond read from 3D coordinates counts too.

    Whether the RING is trans is read relative to known atoms, never from the
    E/Z label: E/Z is a CIP label, and once a substituent outranks the ring
    carbon (1-chloro-cis-cyclooctene, Cl/C1=C/CCCCCC1) "E" means the ring is
    cis (Arrowfield's review of #73). FindPotentialStereo gives cis or trans
    relative to controllingAtoms[0] and [2]; with exactly one of the two in
    the ring, the ring's relation is the opposite."""
    from rdkit import Chem
    out = []
    for si in Chem.FindPotentialStereo(mol):
        if si.type != Chem.StereoType.Bond_Double or str(si.specified) != "Specified":
            continue
        b = mol.GetBondWithIdx(si.centeredOn)
        if not b.IsInRing():
            continue
        size, quad = _ring_neighbours(mol, b)
        if size > SMALL_TRANS_RING:
            continue
        ring = set(quad)
        ref = (si.controllingAtoms[0], si.controllingAtoms[2])
        n_in = sum(a in ring for a in ref)
        trans = si.descriptor == Chem.StereoDescriptor.Bond_Trans
        if (trans if n_in != 1 else not trans):
            out.append(b.GetIdx())
    return out


def configuration_signs(mol, conf_id=-1):
    """The configuration of each stereo element RDKit does not perceive, read
    from one conformer: {element: sign}. Review 6 S1: an ensemble mixed
    trans-cyclooctene helicities and allene axes because every stereo check
    relied on RDKit's model, which has neither.
      - a trans double bond in a ring of SMALL_TRANS_RING or fewer atoms
        (the in-ring C-C=C-C torsion is beyond 90 degrees): its sign;
      - a stereogenic cumulene (cumulene_chains): the torsion between the
        highest-ranked heavy substituent at each end, its sign for an axis,
        cis or trans for a planar one;
      - an atropisomeric biaryl axis (atropisomer_axes): the torsion between
        the highest-ranked ortho atoms, its sign.
    mol is heavy atoms only, in one fixed order, so two conformers' results
    compare element by element."""
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms
    conf = mol.GetConformer(conf_id)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    out = {}
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE or not b.IsInRing():
            continue
        size, quad = _ring_neighbours(mol, b)
        d = rdMolTransforms.GetDihedralDeg(conf, *quad)
        if size <= SMALL_TRANS_RING and abs(d) > 90:
            out[f"ring trans bond {b.GetIdx()}"] = 1 if d > 0 else -1
    for bonds, stereogenic in cumulene_chains(mol):
        if not stereogenic:
            continue
        atoms = [x for k in bonds for x in (mol.GetBondWithIdx(k).GetBeginAtomIdx(),
                                            mol.GetBondWithIdx(k).GetEndAtomIdx())]
        ends = [a for a in set(atoms) if atoms.count(a) == 1]
        chain = set(atoms)
        subs = []
        for e in sorted(ends):
            heavy = [n.GetIdx() for n in mol.GetAtomWithIdx(e).GetNeighbors()
                     if n.GetIdx() not in chain]
            subs.append(max(heavy, key=lambda x: (ranks[x], -x)) if heavy else None)
        if None in subs:
            continue
        (e0, e1), (s0, s1) = sorted(ends), subs
        d = rdMolTransforms.GetDihedralDeg(conf, s0, e0, e1, s1)
        axial = len(bonds) % 2 == 0
        out[f"cumulene {bonds}"] = ((1 if d > 0 else -1) if axial else
                                    ("trans" if abs(d) > 90 else "cis"))
    for k in atropisomer_axes(mol):
        b = mol.GetBondWithIdx(k)
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        oi = max((n.GetIdx() for n in b.GetBeginAtom().GetNeighbors() if n.GetIdx() != j),
                 key=lambda x: ranks[x])
        oj = max((n.GetIdx() for n in b.GetEndAtom().GetNeighbors() if n.GetIdx() != i),
                 key=lambda x: ranks[x])
        out[f"biaryl axis {k}"] = 1 if rdMolTransforms.GetDihedralDeg(conf, oi, i, j, oj) > 0 else -1
    return out


def helicene_lengths(mol):
    """The ring count of each helicene: a chain of aromatic six-rings, each
    fused to the next, where every middle ring is angularly fused (its two
    fusion bonds one bond apart) and every turn goes the same way (consecutive
    rings share the inner atom). A zigzag phenacene alternates its turns and
    stays flat. Returns the longest such chain from each starting ring."""
    ri = mol.GetRingInfo()
    rings = [r for r in ri.AtomRings()
             if len(r) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    if len(rings) < 4 or len(rings) > 40:
        return []

    def fusion(i, j):
        s = set(rings[i]) & set(rings[j])
        return frozenset(s) if len(s) == 2 else None

    nbr = {i: [j for j in range(len(rings)) if j != i and fusion(i, j)]
           for i in range(len(rings))}

    def turn(i, f_in, f_out):
        """(u, v) for ring i entered through fusion bond f_in and left through
        f_out when the two are angular: u in f_in and v in f_out are joined by
        a ring edge (the inner side). None when linear or otherwise, or when
        the inner atoms belong to a third ring: coronene's outer rings wind
        round a central ring that holds their inner atoms, and stay flat."""
        r = rings[i]
        for k in range(6):
            for u, v in ((r[k], r[(k + 1) % 6]), (r[(k + 1) % 6], r[k])):
                if u in f_in - f_out and v in f_out - f_in:
                    if ri.NumAtomRings(u) > 2 or ri.NumAtomRings(v) > 2:
                        return None
                    return u, v
        return None

    best = []

    def walk(path, last_v):
        best.append(len(path))
        i = path[-1]
        for j in nbr[i]:
            if j in path:
                continue
            if len(path) == 1:
                walk(path + [j], None)
                continue
            t = turn(i, fusion(path[-2], i), fusion(i, j))
            # a helix winds when each turn starts where the previous one ended
            if t is None or (last_v is not None and t[0] != last_v):
                continue
            walk(path + [j], t[1])

    for s in range(len(rings)):
        walk([s], None)
    return best


def cumulene_chains(mol):
    """[(bond indices, stereogenic)] for each run of cumulated double bonds
    (an allene C=C=C, a butatriene C=C=C=C). RDKit reports each of their
    double bonds as a potential E/Z bond, the achiral ones too. A cumulene is
    stereogenic (axial for an allene, cis/trans for a butatriene) only when
    both terminal atoms carry two different substituents: 3-methylbuta-1,2-
    diene's (CH3)2C= end makes it achiral."""
    from rdkit import Chem
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    dbl = Chem.BondType.DOUBLE

    def central(a):
        return a.GetDegree() == 2 and all(b.GetBondType() == dbl for b in a.GetBonds())

    seen, out = set(), []
    for c in mol.GetAtoms():
        if not central(c) or c.GetIdx() in seen:
            continue
        chain_atoms, bonds, ends = {c.GetIdx()}, set(), []
        for b0 in c.GetBonds():                     # walk out each way
            prev, b = c, b0
            while True:
                bonds.add(b.GetIdx())
                nxt = b.GetOtherAtom(prev)
                chain_atoms.add(nxt.GetIdx())
                if not central(nxt):
                    ends.append(nxt)
                    break
                b = next(x for x in nxt.GetBonds() if x.GetIdx() != b.GetIdx())
                prev = nxt
        seen |= chain_atoms

        def unsymmetric(t):
            subs = [n.GetIdx() for n in t.GetNeighbors() if n.GetIdx() not in chain_atoms]
            if t.GetTotalNumHs() + len(subs) != 2:
                return False
            if len(subs) == 1:
                return True                         # a group and an H
            return len(subs) == 2 and ranks[subs[0]] != ranks[subs[1]]

        out.append((sorted(bonds), all(unsymmetric(t) for t in ends)))
    return out


def stereo_left_open(mol):
    """What the input does not fix, for --require-stereo: a list of
    (code, message). Empty when every stereo element is fixed, or there is
    none. A real 3D geometry fixes everything (coordinates_fix_stereo).
    Otherwise:
      - every potential stereo element RDKit finds that is Unspecified, or
        marked Unknown (a crossed double bond, a wavy single bond), except
        one that exists only because of an isotope label: the ESP cannot see
        an isotope, and the gate already says isotopes are ignored;
      - each stereogenic cumulene, judged as a whole (cumulene_chains)
        instead of bond by bond;
      - each atropisomeric biaryl axis (atropisomer_axes);
      - each helicene of HELICENE_MIN_RINGS or more rings."""
    from rdkit import Chem
    if coordinates_fix_stereo(mol):
        return []
    plain = Chem.Mol(mol)
    for a in plain.GetAtoms():
        a.SetIsotope(0)
    without_isotopes = {(str(e.type), e.centeredOn) for e in Chem.FindPotentialStereo(plain)}
    chains = cumulene_chains(mol)
    in_cumulene = {k for bonds, _ in chains for k in bonds}
    out = [(f"cumulene at bonds {bonds}", "an axis or cis/trans no SMILES can fix")
           for bonds, stereogenic in chains if stereogenic]
    for e in Chem.FindPotentialStereo(mol):
        if str(e.specified) == "Specified" or (str(e.type), e.centeredOn) not in without_isotopes:
            continue
        if "Bond" in str(e.type) and e.centeredOn in in_cumulene:
            continue                                # judged as a whole above
        what = ("atom" if "Atom" in str(e.type) else "bond")
        out.append((f"{str(e.type)} {what} {e.centeredOn}",
                    "marked unknown" if str(e.specified) == "Unknown" else "unspecified"))
    for k in atropisomer_axes(mol):
        out.append((f"biaryl axis at bond {k}", "an atropisomer no SMILES can fix"))
    for k in small_ring_trans_bonds(mol):
        out.append((f"trans double bond {k} in a ring of {SMALL_TRANS_RING} or fewer atoms",
                    "planar chirality (the ring's helicity) no SMILES can fix"))
    longest = max(helicene_lengths(mol), default=0)
    if longest >= HELICENE_MIN_RINGS:
        out.append((f"[{longest}]helicene", "a helix no SMILES can fix"))
    return out


def check_structure(mol, basis="def2-tzvp", spin=0, maps_from=None, smiles=None,
                    charge=None, require_stereo=False):
    """Can this drawing be translated into SMILES unambiguously, and is the
    result a molecule this pipeline can compute?

    A drawn structure can be perfectly legible to a chemist and still be
    ambiguous to a machine. Each finding below names what to change IN THE
    DRAWING, because that is where the fix belongs.

    maps_from: the parsed molecule WITH its hydrogens, for the atom-map rules
    (callers pass RemoveHs(mol) here, which drops mapped hydrogens).

    charge and spin are what the SCF will be given (--charge, --spin). The
    SCF receives nuclei, charge and spin; review 5 R2 found the gate checking
    only the drawing, so benzene --charge 2 was recorded as verified benzene.

    require_stereo (--require-stereo) makes every stereo element the input
    leaves open a BLOCK (stereo_left_open); without it stereo only warns or
    notes, as before. Written against tests/data/stereo, built first.

    Returns (verdict, findings) where findings is a list of
    (severity, code, message, remedy).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    out = atom_map_findings(maps_from if maps_from is not None else mol, smiles)

    # --- 1. does it survive sanitisation at all ---
    probs = Chem.DetectChemistryProblems(mol)
    for p in probs:
        out.append((BLOCK, p.GetType(), p.Message(),
                    "fix the valence or charge on that atom in the drawing"))
    if probs:
        return BLOCK, out

    # --- 2. things that are not a single closed-shell molecule ---
    wild = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if wild:
        out.append((BLOCK, "WILDCARD_ATOM",
                    f"{len(wild)} wildcard/R-group atom(s) at index {wild}",
                    "replace R / * / query atoms with explicit atoms; a generic "
                    "structure has no electron density"))
    formal = Chem.GetFormalCharge(mol)
    if charge is not None and charge != formal:
        out.append((BLOCK, "CHARGE_OVERRIDE",
                    f"--charge {charge:+d} but the structure's formal charge is "
                    f"{formal:+d}: the species computed would not be the one drawn",
                    "draw the charge on the atom that carries it, or drop --charge"))
    rad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    if rad and spin == 0:
        out.append((BLOCK, "OPEN_SHELL",
                    f"{rad} unpaired electron(s): this is a radical, and no "
                    f"--spin was given, so it would be computed with a "
                    f"RESTRICTED wavefunction that cannot represent it",
                    f"pass --spin {rad} to use unrestricted DFT (UKS), or "
                    f"complete the valences / add the missing charge"))
    elif rad and spin == rad:
        out.append((INFO, "OPEN_SHELL_UNRESTRICTED",
                    f"{rad} unpaired electron(s); --spin {spin} selected, so "
                    f"UKS will be used", "none needed"))
    elif spin and spin != rad:
        # a deliberate triplet (O2 drawn O=O, --spin 2) is legitimate, so this
        # warns rather than refuses -- and it is written to the record
        out.append((WARN, "SPIN_OVERRIDE",
                    f"--spin {spin} but the structure as drawn has {rad} unpaired "
                    f"electron(s); it will be computed with UKS at spin {spin}, "
                    f"which is not the electronic state the drawing shows",
                    "draw the radical centres, or check the multiplicity"))
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) > 1:
        names = ", ".join(Chem.MolToSmiles(f) for f in frags[:4])
        out.append((BLOCK, "MULTIPLE_FRAGMENTS",
                    f"{len(frags)} disconnected fragments: {names}",
                    "draw one molecule. A salt or solvate has no single ESP "
                    "surface, and the counter-ion dominates the potential"))
    exotic = sorted({a.GetSymbol() for a in mol.GetAtoms()}
                    - TESTED_ELEMENTS - {"*"})
    if exotic:
        out.append((WARN, "UNTESTED_ELEMENT",
                    f"element(s) outside the tested set: {', '.join(exotic)}",
                    "def2 basis sets cover these, but closed-shell RKS and the "
                    "rho=0.002 surface convention were not validated on them"))

    # --- 3. ambiguity that SMILES would silently carry through ---
    si = Chem.FindPotentialStereo(mol)
    unspec_atoms = [e.centeredOn for e in si
                    if str(e.type) == "Atom_Tetrahedral"
                    and str(e.specified) == "Unspecified"]
    unspec_bonds = [e.centeredOn for e in si
                    if str(e.type) == "Bond_Double"
                    and str(e.specified) == "Unspecified"]
    if unspec_atoms:
        n = len(unspec_atoms)
        if n == 1:
            msg = ("1 stereocentre drawn without a wedge or hash. The two "
                   "enantiomers have mirror-image ESP: every value identical, "
                   "the picture reflected.")
            sev = INFO
        else:
            msg = (f"{n} stereocentres drawn without wedges or hashes: "
                   f"{2 ** n} stereoisomers, and they are diastereomers, i.e. "
                   f"genuinely different molecules. Measured on threonine, the "
                   f"ESP EXTREMA differ by ~1.6-1.9 kcal/mol, about as much as "
                   f"conformer choice moves them (up to 1.25), and the SHAPE of "
                   f"the surface and where the features sit on it differ too.")
            sev = WARN
        out.append((sev, "UNSPECIFIED_TETRAHEDRAL",
                    msg + f" Atom index: {unspec_atoms}",
                    "add wedge/hash bonds in the drawing, or accept that the "
                    "map shows one arbitrary isomer"))
    if unspec_bonds:
        out.append((WARN, "UNSPECIFIED_DOUBLE_BOND",
                    f"{len(unspec_bonds)} double bond(s) with undrawn E/Z "
                    f"geometry at bond index {unspec_bonds}. E and Z are "
                    f"different compounds with different shapes.",
                    "draw the substituents cis or trans explicitly"))
    if require_stereo:
        left = stereo_left_open(mol)
        if left:
            out.append((BLOCK, "STEREO_REQUIRED",
                        f"--require-stereo, but {len(left)} stereo element(s) are not "
                        f"fixed by the input: "
                        + "; ".join(f"{w} ({why})" for w, why in left),
                        "give every centre a wedge or @/@@, every ordinary double bond "
                        "a drawn or /\\ geometry; a stereogenic allene or cumulene, an "
                        "atropisomer or a helix can only be fixed by 3D coordinates "
                        "(MOL/SDF with z, XYZ), as SMILES marks on them are dropped"))

    # --- 4. tautomers: the drawn form is a choice, often an unconscious one ---
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        te = rdMolStandardize.TautomerEnumerator()
        te.SetMaxTautomers(32)
        n_taut = len(te.Enumerate(Chem.Mol(mol)))
        if n_taut > 1:
            out.append((INFO, "TAUTOMERS",
                        f"{n_taut} tautomers reachable from the drawn form; the "
                        f"ESP is computed for exactly the one you drew",
                        "confirm the drawn tautomer is the one you mean"))
    except Exception:
        pass

    # --- 5. things that change the OUTPUT, not the input ---
    q = Chem.GetFormalCharge(mol)
    if q != 0:
        out.append((INFO, "CHARGED",
                    f"net charge {q:+d}: the surface will not change sign, so a "
                    f"diverging colormap centred on zero is the wrong instrument",
                    "none needed; a sequential scale is selected automatically"))
    iso = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsotope()]
    if iso:
        out.append((INFO, "ISOTOPES",
                    f"{len(iso)} isotopic label(s): ignored, the electronic "
                    f"structure does not depend on nuclear mass",
                    "none needed"))

    # --- 6. round trip: does SMILES actually carry this structure? ---
    try:
        smi = Chem.MolToSmiles(mol)
        back = Chem.MolFromSmiles(smi)
        if back is None:
            out.append((BLOCK, "ROUNDTRIP_FAIL",
                        f"the generated SMILES {smi!r} does not parse back",
                        "simplify the drawing; this is an RDKit round-trip bug "
                        "or an exotic bonding pattern"))
        elif Chem.MolToSmiles(back) != smi:
            out.append((WARN, "ROUNDTRIP_DRIFT",
                        f"SMILES is not stable under a round trip: {smi} -> "
                        f"{Chem.MolToSmiles(back)}",
                        "the structure carries information SMILES cannot hold"))
    except Exception as e:
        out.append((BLOCK, "ROUNDTRIP_ERROR", str(e), "simplify the drawing"))

    # --- 7. cost ---
    try:
        heavy = mol.GetNumHeavyAtoms()
        if heavy > 25:
            out.append((INFO, "SIZE",
                        f"{heavy} heavy atoms; the SCF will dominate the runtime",
                        "use --basis def2-svp for a first look"))
    except Exception:
        pass

    sev = BLOCK if any(s == BLOCK for s, *_ in out) else (
        WARN if any(s == WARN for s, *_ in out) else INFO)
    return sev, out


def print_check(verdict, findings, label=""):
    icon = {BLOCK: "BLOCK", WARN: "WARN ", INFO: "note "}
    head = {BLOCK: "REJECTED", WARN: "ACCEPTED WITH WARNINGS", INFO: "ACCEPTED"}
    print(f"--- structure check{(' — ' + label) if label else ''} ---")
    if not findings:
        print("  nothing to flag: unambiguous, closed-shell, single molecule")
    for sev, code, msg, remedy in findings:
        print(f"  [{icon[sev]}] {code}")
        print(f"          {msg}")
        print(f"          fix: {remedy}")
    print(f"verdict: {head[verdict]}")
    return verdict


def _quiet_rdkit():
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")


def parse_for_check(smiles=None, path=None, charge=0, record=None, atoms=None):
    """Parse to an RDKit mol WITHOUT generating 3D coordinates.

    The admissibility gate must run before embedding: a radical or an R-group
    structure should be rejected in milliseconds, not after twenty conformers
    have been generated and a force field has complained about them.
    Returns (mol or None, list of BLOCK findings explaining a parse failure).

    atoms is an --atoms string ("Sym x y z; ..."). Its bonds are perceived
    from the coordinates the way an XYZ file's are, so it passes the same
    gate. Review 4 N3: --atoms skipped the gate, and two waters 8 A apart
    came out validity.input: true with a site table.
    """
    from rdkit import Chem
    _quiet_rdkit()
    if atoms is not None:
        try:
            m = _mol_from_atoms(atoms)
        except ValueError as e:
            return None, [(BLOCK, "ATOMS_SYNTAX", str(e),
                           "write --atoms as 'Sym x y z; Sym x y z; ...' in Angstrom")]
        m = _perceive_bonds(m, charge, "ATOMS")
        return (None, m) if isinstance(m, list) else (m, [])
    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol, []
        raw = Chem.MolFromSmiles(smiles, sanitize=False)
        probs = Chem.DetectChemistryProblems(raw) if raw is not None else []
        if probs:
            return None, [(BLOCK, p.GetType(), p.Message(),
                           "fix the valence or charge on that atom in the drawing")
                          for p in probs]
        return None, [(BLOCK, "UNPARSEABLE", f"RDKit cannot read {smiles!r}",
                       "check the SMILES syntax")]
    mol, problems, _ = read_structure(path, charge=charge, record=record)
    return mol, problems


def _mol_from_atoms(spec):
    """An --atoms string as an RDKit mol: atoms and a conformer, no bonds."""
    from rdkit import Chem
    rows = [t.split() for t in spec.split(";") if t.strip()]
    if not rows or any(len(r) != 4 for r in rows):
        raise ValueError(f"cannot read {spec[:60]!r}: every atom needs a symbol "
                         f"and three coordinates")
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(rows))
    for i, (sym, *xyz) in enumerate(rows):
        try:
            rw.AddAtom(Chem.Atom(sym.capitalize()))
            conf.SetAtomPosition(i, [float(v) for v in xyz])
        except (RuntimeError, ValueError):
            raise ValueError(f"atom {i}: {' '.join([sym, *xyz])!r} is not 'Sym x y z'") from None
    conf.Set3D(True)
    rw.AddConformer(conf, assignId=True)
    return rw.GetMol()


def _perceive_bonds(m, charge, fmt):
    """Coordinates only -> mol with perceived connectivity and bond orders.

    Rebuilds the molecule from elements and positions, so no bond the reader
    guessed survives, then runs RDKit's xyz2mol. Marks the result
    _bonds=perceived so from_file can tell orders read from a file apart from
    orders inferred here. Returns the mol, or a list holding one BLOCK.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    rw = Chem.RWMol()
    for a in m.GetAtoms():
        rw.AddAtom(Chem.Atom(a.GetAtomicNum()))
    conf = Chem.Conformer(m.GetNumAtoms())
    src = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        conf.SetAtomPosition(i, src.GetAtomPosition(i))
    conf.Set3D(True)
    rw.AddConformer(conf, assignId=True)
    mol = rw.GetMol()
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
        Chem.SanitizeMol(mol)              # AddHs needs implicit valences
    except Exception as e:
        return [(BLOCK, f"{fmt}_CONNECTIVITY",
                 f"could not perceive bonds from {fmt} coordinates at charge "
                 f"{charge:+d}: {e}",
                 "include every hydrogen in the file and pass --charge if the "
                 "species is not neutral, or supply a format that carries bond "
                 "orders (MOL/SDF/CDXML)")]
    mol.SetProp("_bonds", "perceived")
    return mol


def _molfile_records(path):
    """The MOL blocks of a MOL/SDF file, in file order: the text before each
    $$$$ that holds a counts line (V2000 or V3000). Text that is not a record
    -- a comment or a converter's trailer after the last $$$$ -- is not one,
    so it cannot shift which header belongs to which record (review 6 S15: a
    C2F4 SDF with one trailing line was read under no header and discarded)."""
    import re
    from pathlib import Path
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
    except OSError:
        return []
    out = []
    for k, rec in enumerate(re.split(r"^\$\$\$\$[^\n]*$", text, flags=re.M)):
        if k:
            rec = rec[1:] if rec.startswith("\n") else rec
        if re.search(r"^.{33,39}V[23]000\s*$", rec, flags=re.M):
            out.append(rec)
    return out


def _molfile_header_dims(path):
    """The dimension code on line 2 of each MOL/SDF record ("2D", "3D" or
    ""), in file order. coordinate_kind needs what the file says: RDKit's
    conformer flag turns a 3D header into 2D when every z is 0 and a bond
    carries a stereo mark."""
    out = []
    for rec in _molfile_records(path):
        lines = rec.split("\n")
        out.append(lines[1][20:22].strip().upper() if len(lines) > 1 else "")
    return out


def _as_drawing(m, block):
    """A record judged a drawing, re-read as the 2D drawing it is: z set to
    0 and the header to 2D, so its stereo comes from its wedges and drawn
    double bonds only. Read under a 3D header, RDKit takes chirality from z,
    and z noise in an unwedged drawing then fixes a configuration nobody drew
    (review 6 S7: two seeds of one threonine drawing gave two diastereomers,
    both passing --require-stereo). The record's properties are kept."""
    from rdkit import Chem
    lines = block.split("\n")
    if len(lines) < 4 or "V2000" not in lines[3]:
        return m                                    # V3000: left as RDKit read it
    n_at = int(lines[3][:3])
    lines[1] = lines[1][:20] + "2D" + lines[1][22:] if len(lines[1]) >= 22 else lines[1]
    for k in range(4, 4 + n_at):
        lines[k] = lines[k][:20] + f"{0.0:10.4f}" + lines[k][30:]
    d = Chem.MolFromMolBlock("\n".join(lines), removeHs=False)
    if d is None or d.GetNumAtoms() != m.GetNumAtoms():
        return m
    d.GetConformer().Set3D(False)
    for p in m.GetPropNames(includePrivate=True, includeComputed=False):
        if not d.HasProp(p):
            d.SetProp(p, m.GetProp(p))
    return d


def read_structure(path, charge=0, record=None):
    """The single structure reader. Returns (mol or None, problems, n_found).

    record (1-based) picks one record from a multi-record file; without it a
    file holding more than one record is a BLOCK. record="all" returns
    ([(mol or None, problems)] for every record in file order, [], n_found)
    from one parse: a supplied ensemble re-read the file once per record,
    quadratic in its length (review 6 S13).

    charge is needed only for XYZ, which carries coordinates and no bonds:
    bond perception has to be told the net charge or it fails. v4.0 hard-coded
    0 while telling the user to \"set the charge explicitly\", naming a remedy
    the code did not implement."""
    from rdkit import Chem
    _quiet_rdkit()
    ext = str(path).lower().rsplit(".", 1)[-1]
    blocks = []                                     # MOL blocks, aligned with mols
    try:
        if ext == "cdxml":
            mols = [m for m in Chem.MolsFromCDXMLFile(str(path)) if m is not None]
        elif ext in ("sdf", "mol"):
            mols = (list(Chem.SDMolSupplier(str(path), removeHs=False))
                    if ext == "sdf" else
                    [Chem.MolFromMolFile(str(path), removeHs=False)])
            recs = _molfile_records(path)
            if len(recs) == len(mols):
                blocks = recs
                for m, rec in zip(mols, recs):
                    if m is not None:
                        lines = rec.split("\n")
                        m.SetProp("_header_dim", lines[1][20:22].strip().upper()
                                  if len(lines) > 1 else "")
        elif ext in ("xyz", "pdb"):
            # Neither format carries trustworthy bond ORDERS. XYZ has no bonds
            # at all (v3 rejected every xyz as "MULTIPLE_FRAGMENTS"). PDB has
            # proximity bonds, all single unless CONECT says otherwise: review 3
            # H1 read chlorobenzene.pdb as chlorocyclohexane, AddHs filled the
            # "missing" valences, and every check passed. Both now go through
            # one route: keep the coordinates, discard any bonds, perceive
            # connectivity and orders from geometry.
            if ext == "xyz":
                m = Chem.MolFromXYZFile(str(path))
            else:
                m = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False,
                                        proximityBonding=False)
            if m is not None:
                m = _perceive_bonds(m, charge, ext.upper())
                if isinstance(m, list):              # a BLOCK finding
                    return None, m, 0
            mols = [m]
        else:
            return None, [(BLOCK, "UNSUPPORTED_FORMAT",
                           f"'.{ext}' is not a format this reads",
                           "export as MOL, SDF, CDXML, XYZ or PDB")], 0
    except Exception as e:
        return None, [(BLOCK, "READ_ERROR", str(e), "check the file is intact")], 0
    # Review 6 S6: an unreadable record keeps its place. Dropping it renumbered
    # the rest (--record 2 read the third record) and let a two-record file
    # whose second record failed pass as one record.
    if all(m is None for m in mols):
        return None, [(BLOCK, "UNPARSEABLE",
                       f"RDKit read no valid molecule from {path}",
                       "re-export from the drawing program; check for R-groups "
                       "or query atoms the exporter could not represent")], 0
    # One geometry in, one record out. A multi-record file (a docking run's
    # poses, a CDXML page of fragments) used to warn and silently take the
    # first; in a pipeline 20 poses became an analysis of pose 1.
    n = len(mols)
    if record is None and n > 1:
        return None, [(BLOCK, "MULTIPLE_RECORDS",
                       f"{path} holds {n} records; one run computes one",
                       f"pass --record N (1..{n}) to choose one, or split the file")], n

    def finish(m, r):
        if m is None:
            return None, [(BLOCK, "RECORD_UNREADABLE",
                           f"record {r} of {path} could not be read",
                           f"fix or remove record {r}; the others keep "
                           f"their numbers")]
        m.SetProp("_format", ext)
        kind = coordinate_kind(m, ext)
        if kind == "ambiguous":
            return None, [(BLOCK, "DIMENSION_AMBIGUOUS",
                           f"{path} is flagged 3D but flat (every atom within "
                           f"{FLAT_TOL} A of one plane), and its bonds and angles "
                           f"cannot tell a drawing from a real planar molecule",
                           "write a 2D header if this is a drawing; if it is a "
                           "geometry, save it as .xyz")]
        if kind == "drawing" and blocks and m.GetConformer().Is3D():
            m = _as_drawing(m, blocks[r - 1])
            m.SetProp("_format", ext)
        m.SetProp("_coords", kind)
        return m, []

    if record == "all":
        return [finish(m, r) for r, m in enumerate(mols, start=1)], [], n
    if record is not None and not 1 <= record <= n:
        return None, [(BLOCK, "RECORD_OUT_OF_RANGE",
                       f"--record {record}, but {path} holds {n} record(s)",
                       f"pass --record between 1 and {n}")], n
    m, problems = finish(mols[(record or 1) - 1], record or 1)
    return m, problems, n


def _skeleton(m):
    """Elements and connectivity only: every bond single, no charge, no
    radical, no isotope. What two Lewis structures of one molecule share."""
    from rdkit import Chem
    rw = Chem.RWMol()
    for a in m.GetAtoms():
        b = Chem.Atom(a.GetAtomicNum())
        b.SetNoImplicit(True)
        rw.AddAtom(b)
    for bd in m.GetBonds():
        rw.AddBond(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx(), Chem.BondType.SINGLE)
    s = rw.GetMol()
    s.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(s)
    return s


def verify_computed(spec, approved, max_stereo=1024):
    """Is the geometry PySCF will receive the structure the gate approved?
    Returns (ok, computed_smiles, reason).

    Review 3 H2: the old assertion compared the input string with itself (and,
    for files, the reader with itself), and never read `spec`. This reads only
    `spec`. The SCF sees nuclei, charge and spin, never bond orders, so what is
    checked is what the coordinates fix:
      - connectivity, hydrogens included: perceived from distances and matched
        to the approval's skeleton. Catches another molecule, an isomer, a
        tautomer, and a reader that guessed bond orders and so added the wrong
        hydrogens;
      - stereo: assigned from 3D on the matched approval, and accepted if it is
        one completion of what the approval specified. Catches E for Z and the
        wrong enantiomer; a centre the user left open may come out either way.
    Bond orders and charge placement are taken from the approval, so another
    resonance form ([N-]=C=S vs N#C[S-]) is not a false refusal. Isotope
    positions cannot be checked: coordinates carry elements only.

    Atom maps are stripped first: this is CHEMICAL identity. Maps make
    graph-equivalent groups distinct labels (isopropyl's two methyls), RDKit
    then perceives a stereocentre the chemistry does not have, and a correctly
    mapped isopropanol was refused. Which atom is which is verify_mapping's job.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Chem.EnumerateStereoisomers import (EnumerateStereoisomers,
                                                   StereoEnumerationOptions)
    approved = Chem.Mol(approved)
    for a in approved.GetAtoms():
        a.SetAtomMapNum(0)
    rows = [t.split() for t in spec.split(";")]
    g = Chem.MolFromXYZBlock(f"{len(rows)}\n\n" + "\n".join(" ".join(r) for r in rows))
    rdDetermineBonds.DetermineConnectivity(g)
    tmpl = Chem.AddHs(approved)
    if tmpl.GetNumAtoms() != g.GetNumAtoms():
        return False, None, (f"{g.GetNumAtoms()} atoms in the geometry, "
                             f"{tmpl.GetNumAtoms()} in the approved structure")
    sk_t, sk_g = _skeleton(tmpl), _skeleton(g)
    match = (sk_g.GetSubstructMatch(sk_t)
             if sk_t.GetNumBonds() == sk_g.GetNumBonds() else ())
    if not match:
        # review 5 R8: say what WAS computed, so a forced record names it
        return False, _perceived_smiles(g, Chem.GetFormalCharge(approved)), (
            "the connectivity perceived from the coordinates is not the "
            "approved structure's")
    conf = Chem.Conformer(tmpl.GetNumAtoms())
    for i, j in enumerate(match):
        conf.SetAtomPosition(i, [float(x) for x in rows[j][1:4]])
    conf.Set3D(True)
    tmpl.RemoveAllConformers()
    tmpl.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(tmpl)
    computed = Chem.MolToSmiles(Chem.RemoveHs(tmpl))
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True,
                                    maxIsomers=max_stereo + 1)
    allowed = {Chem.MolToSmiles(m) for m in EnumerateStereoisomers(approved, options=opts)}
    if len(allowed) > max_stereo:
        return True, computed, (f"more than {max_stereo} stereo completions; "
                                f"connectivity verified, stereo not")
    if computed not in allowed:
        return False, computed, "the stereochemistry in the coordinates differs"
    return True, computed, ""


def _perceived_smiles(g, charge):
    """The molecule a set of coordinates is, as far as the coordinates say:
    bond orders perceived at the approved total charge, or connectivity
    alone when orders cannot be assigned. For the record of a forced run."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    try:
        m = Chem.Mol(g)
        rdDetermineBonds.DetermineBondOrders(m, charge=charge)
        return Chem.MolToSmiles(Chem.RemoveHs(m))
    except Exception:
        return Chem.MolToSmiles(Chem.RemoveHs(_skeleton(g), sanitize=False))


def input_invalid_reason(gate, verified, mapped_ok):
    """Why validity.input is false, naming every failing check (review 5 R8:
    a forced atom-order failure was reported as 'forced past the gate, or
    structure not verified', neither of which had happened)."""
    causes = []
    if gate is None:
        causes.append("no structure gate ran")
    elif gate["forced"]:
        causes.append("the gate's BLOCK was forced")
    if verified is False:
        causes.append("structure not verified (input.computed_smiles is what was computed)")
    if mapped_ok is False:
        # the record has no separate field for this until a schema change
        # (spec section 6, "Deferred"); this reason is where it is recorded
        causes.append("atom order not verified (verify_mapping)")
    return "input not valid: " + ("; ".join(causes) or "unknown")


def check_identity(spec, probe, approved, kept, force, stage):
    """Run both identity checks on `spec` and refuse on a mismatch.

    Returns (verified, computed_smiles, mapped_ok). With force a mismatch
    warns instead, and the caller records it (validity.input false). Called
    before the SCF and again after --optimize: review 4 N1 found the
    optimiser replacing the coordinates after the only check had run, so a
    gas-phase glycine zwitterion was computed as neutral glycine under a
    record that said structure_verified: true.
    """
    verified, computed, why = verify_computed(spec, _rdkit_chem().RemoveHs(probe))
    if not verified:
        msg = (f"{stage}: the structure that passed the check is not the one "
               f"about to be computed: {why}\n  approved:  {approved}\n  "
               f"computing: {computed or '(no match)'}")
        if not force:
            raise SystemExit(msg + "\nrefusing to continue (--force computes it anyway)")
        print(f"WARNING: {msg}\nthe computed structure could not be verified; "
              f"continuing because of --force", file=sys.stderr)
    elif why:
        print(f"note: {why}", file=sys.stderr)
    # Index by index, so per-site rows and atom maps point at the atoms the
    # input named; against the file's own coordinates when they were kept.
    mapped_ok, why_m = verify_mapping(spec, probe, reference=kept)
    if not mapped_ok and verified:
        msg = f"{stage}: atom order of the computed geometry is not the input's: {why_m}"
        if not force:
            raise SystemExit(msg + "\nrefusing to continue (--force computes it anyway)")
        print(f"WARNING: {msg}; continuing because of --force", file=sys.stderr)
    return verified, computed, mapped_ok


def _rdkit_chem():
    from rdkit import Chem
    return Chem


def verify_mapping(spec, probe, reference=False, tol=1e-3):
    """Is spec atom i the approval's atom i? Returns (ok, reason).

    verify_computed searches for SOME match between spec and the approval, so
    it accepts a geometry whose atoms were reordered, and a site table keyed
    by index (and by atom map) would then point at the wrong atoms. This
    checks the correspondence itself (per-site-spec section 1, `mapped`):
      - element and bonded neighbours of every atom, by index, hydrogens
        included, against AddHs(probe), whose order is the spec's order;
      - with reference=True (a file whose 3D coordinates were kept, e.g. a
        docking pose), every atom the file supplied must still be where the
        file put it: map 5 points at the coordinates the adapter wrote.
    A consistent swap of graph-equivalent groups passes the first test and
    only the second can catch it. Bond orders are not compared, as in
    verify_computed: the SCF never sees them.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    rows = [t.split() for t in spec.split(";")]
    tmpl = Chem.AddHs(probe)
    if tmpl.GetNumAtoms() != len(rows):
        return False, f"{len(rows)} atoms in the geometry, {tmpl.GetNumAtoms()} in the approval"
    for i, a in enumerate(tmpl.GetAtoms()):
        if a.GetSymbol() != rows[i][0]:
            return False, (f"atom order differs: atom {i} is {rows[i][0]} in the geometry, "
                           f"{a.GetSymbol()} in the approval")
    g = Chem.MolFromXYZBlock(f"{len(rows)}\n\n" + "\n".join(" ".join(r) for r in rows))
    rdDetermineBonds.DetermineConnectivity(g)
    for i in range(len(rows)):
        got = sorted(n.GetIdx() for n in g.GetAtomWithIdx(i).GetNeighbors())
        want = sorted(n.GetIdx() for n in tmpl.GetAtomWithIdx(i).GetNeighbors())
        if got != want:
            return False, (f"atom order differs: atom {i} ({rows[i][0]}) is bonded to "
                           f"{got} in the geometry, {want} in the approval")
    if reference and probe.GetNumConformers():
        P = probe.GetConformer().GetPositions()
        X = np.array([[float(v) for v in r[1:4]] for r in rows[:len(P)]])
        moved = np.linalg.norm(X - P, axis=1)
        if moved.max() > tol:
            i = int(moved.argmax())
            return False, (f"atom {i} ({rows[i][0]}) is {moved[i]:.3f} A from where the "
                           f"file put it; the file's atom order is not the computed one")
    return True, ""


# ---------------------------------------------------------------- ensembles
# docs/ensemble-spec.md sections 1-3, written against tests/data/ensemble.
# Everything here runs BEFORE any SCF: an ensemble that mixes molecules,
# tautomers or stereoisomers is refused here, not averaged later.
CONTACT_MIN_BONDS = 5          # a 1,5 pair of a gauche kink is 4 bonds apart
CONTACT_SLACK = 0.5            # A over the Bondi sum
DEDUP_RMSD = 0.5               # A, heavy atoms; provisional until section 8
N_EMBED = 100
SELECTION_RULES = ("rg_spread",)


class EnsembleRefused(ValueError):
    """An ensemble that may not be computed. .findings holds the gate-style
    (severity, code, message, remedy) tuples, every one a BLOCK."""

    def __init__(self, findings):
        self.findings = findings
        super().__init__("; ".join(f"{c}: {m}" for _, c, m, _ in findings))


def _heavy_xyz(mol, conf_id=-1):
    conf = mol.GetConformer(conf_id)
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return idx, np.array([list(conf.GetAtomPosition(i)) for i in idx])


def intramolecular_contacts(mol, conf_id=-1):
    """Heavy-atom pairs at least CONTACT_MIN_BONDS bonds apart and closer than
    their Bondi radii plus CONTACT_SLACK (spec section 3). Returns
    [(i, j, distance, limit)]. A pair with an atom BONDI lacks is not judged."""
    from rdkit import Chem
    topo = Chem.GetDistanceMatrix(mol)
    idx, X = _heavy_xyz(mol, conf_id)
    out = []
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            if topo[i, j] < CONTACT_MIN_BONDS:
                continue
            ri = BONDI.get(mol.GetAtomWithIdx(i).GetSymbol())
            rj = BONDI.get(mol.GetAtomWithIdx(j).GetSymbol())
            if ri is None or rj is None:
                continue
            d = float(np.linalg.norm(X[a] - X[b]))
            if d < ri + rj + CONTACT_SLACK:
                out.append((i, j, d, ri + rj + CONTACT_SLACK))
    return out


def shape_descriptors(mol, conf_id=-1):
    """Per-member shape (spec section 3), heavy atoms, unweighted, in A:
    radius of gyration, maximum extent, asphericity normalised to [0, 1]
    (0 a sphere, 1 a rod: (l1 - (l2 + l3)/2) / (l1 + l2 + l3) of the gyration
    tensor's eigenvalues, l1 largest), the contacts, and folded, which is
    'has a contact'. Surface genus is added per member once its surface
    exists; it is a diagnostic, never part of the folded rule. With no heavy
    atom (H2) there is no shape: null descriptors, not folded (review 6 S12)."""
    _, X = _heavy_xyz(mol, conf_id)
    if len(X) == 0:
        return dict(rg=None, extent=None, asphericity=None, contacts=[], folded=False)
    Y = X - X.mean(axis=0)
    lam = np.sort(np.linalg.eigvalsh(Y.T @ Y / len(Y)))[::-1]
    D = np.linalg.norm(X[:, None] - X[None], axis=-1)
    contacts = intramolecular_contacts(mol, conf_id)
    return dict(rg=float(np.sqrt(lam.sum())), extent=float(D.max()),
                asphericity=float((lam[0] - 0.5 * (lam[1] + lam[2])) / lam.sum())
                if lam.sum() > 0 else 0.0,
                contacts=[[int(i), int(j), round(d, 3)] for i, j, d, _ in contacts],
                folded=bool(contacts))


def _is_achiral(mol):
    """True when the molecule is its own mirror image as specified: inverting
    every tetrahedral centre gives the same canonical SMILES. Then a
    conformer and its mirror image are the same conformer (gauche+ and
    gauche- butane); for a chiral molecule the mirror image is the other
    enantiomer and is never compared."""
    from rdkit import Chem
    m = Chem.RemoveHs(Chem.Mol(mol))
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    mirror = Chem.Mol(m)
    flip = {Chem.ChiralType.CHI_TETRAHEDRAL_CW: Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW: Chem.ChiralType.CHI_TETRAHEDRAL_CW}
    for a in mirror.GetAtoms():
        a.SetChiralTag(flip.get(a.GetChiralTag(), a.GetChiralTag()))
    return Chem.MolToSmiles(m) == Chem.MolToSmiles(mirror)


def _rg_spread(rg, energy, k):
    """The lowest-energy conformer, plus k-1 spread evenly across the Rg
    range: the question is how much SHAPE moves the ESP, and the k lowest
    minima of a long chain are often near-duplicates.

    Exactly bench/conformers.py select(), the rule the measured conformer
    rows (5.1 / 2.69 / 5.2) used, collisions included: when an Rg pick is
    the lowest-energy conformer, the next most compact unpicked one tops up.
    On the C18 acetate that happens (its lowest energy sits at rank 87 of
    100), so re-spreading instead would change the members those rows were
    measured on. Review 6 S4 found the top-up crowding the compact end
    (1-octanol, K = 3); changing it is a protocol change for section 8, with
    the rows re-measured, not a fix here."""
    best = int(np.argmin(energy))
    order = np.argsort(rg)
    picks = [best]
    for q in np.linspace(0, len(order) - 1, max(k - 1, 0)):
        i = int(order[int(round(q))])
        if i not in picks:
            picks.append(i)
    for i in order:                         # top up if Rg ties collapsed picks
        if len(picks) >= k:
            break
        if int(i) not in picks:
            picks.append(int(i))
    return picks[:k]


SYMMETRY_MATCH_CAP = 10000


def select_ensemble(mol, converged, k, dedup_rmsd=DEDUP_RMSD, rule="rg_spread",
                    max_matches=SYMMETRY_MATCH_CAP):
    """Choose up to k members from the conformers embedded in mol (spec
    section 2). converged is _embed_and_minimise's [(conf_id, energy)].

    1. Deduplicate by symmetry-aware heavy-atom RMSD below dedup_rmsd, in
       order of increasing energy, so the lower-energy copy is kept. The
       symmetry matches are found once and passed to GetBestRMS. For an
       achiral molecule the mirror image of a conformer is compared too.
       A molecule with max_matches or more automorphisms is refused
       (ValueError): a truncated list kept symmetry copies as distinct
       conformers (review 6 S5; three tBu groups give 1296).
    2. Choose by `rule`; only "rg_spread" exists. It ranks by RDKit's
       mass-weighted all-atom radius of gyration, the one bench/conformers.py
       ranked by (review 6 S4); shape_descriptors' heavy-atom rg is the
       per-member descriptor, not the selection key.
    Never pads: k_obtained may be smaller than k (benzene gives 1). The force
    field's energies order and choose; they are returned, never used as
    weights. Returns a dict: members [(conf_id, energy)] in selection order,
    k_requested, k_obtained, n_converged, n_unique, and the parameters."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    if rule not in SELECTION_RULES:
        raise ValueError(f"unknown selection rule {rule!r}; one of {SELECTION_RULES}")
    if not converged:
        raise ValueError("no converged conformer to select from")
    heavy = Chem.RemoveHs(mol)              # its own copy: GetBestRMS moves the probe
    matches = heavy.GetSubstructMatches(heavy, uniquify=False, useChirality=False,
                                        maxMatches=max_matches)
    if len(matches) >= max_matches:
        raise ValueError(f"{max_matches} or more symmetry matches; deduplication "
                         f"would keep symmetry copies")
    maps = [list(zip(m, range(heavy.GetNumAtoms()))) for m in matches]
    achiral = _is_achiral(mol)
    mirror = None
    if achiral:
        mirror = Chem.Mol(heavy)
        for c in mirror.GetConformers():
            for i in range(mirror.GetNumAtoms()):
                p = c.GetAtomPosition(i)
                c.SetAtomPosition(i, (p.x, p.y, -p.z))
    ordered = sorted(converged, key=lambda ce: (ce[1], ce[0]))
    kept = []
    for cid, e in ordered:
        dup = False
        for kid, _ in kept:
            for probe in (heavy, mirror) if achiral else (heavy,):
                if rdMolAlign.GetBestRMS(probe, heavy, cid, kid, map=maps) < dedup_rmsd:
                    dup = True
                    break
            if dup:
                break
        if not dup:
            kept.append((cid, e))
    from rdkit.Chem import Descriptors3D
    rg = np.array([Descriptors3D.RadiusOfGyration(mol, confId=cid) for cid, _ in kept])
    energy = np.array([e for _, e in kept])
    picks = _rg_spread(rg, energy, k)
    return dict(members=[kept[i] for i in picks], k_requested=int(k),
                k_obtained=len(picks), n_converged=len(converged), n_unique=len(kept),
                rule=rule, dedup_rmsd=dedup_rmsd, mirror_images_merged=achiral,
                n_symmetry_matches=len(maps))


def _ensemble_block(code, message, remedy):
    return (BLOCK, code, message, remedy)


def _weight_problems(weight, optimize):
    """MMFF energies never weight (section 2); DFT energies weight only
    DFT-optimised members, because single points on force-field geometries
    carry the force field's strain."""
    if weight in (None, "none"):
        return []
    if weight == "mmff":
        return [_ensemble_block(
            "WEIGHT_MMFF", "force-field energies may order members, never weight them",
            "drop the weights (statistics are unweighted), or use --weight dft "
            "with --optimize")]
    if weight == "dft":
        return [] if optimize else [_ensemble_block(
            "WEIGHT_NEEDS_OPTIMIZE",
            "DFT single points on force-field geometries are not conformer energies",
            "add --optimize so every member is a DFT minimum, or drop --weight")]
    return [_ensemble_block("WEIGHT_UNKNOWN", f"no weighting called {weight!r}",
                            "use --weight dft (with --optimize), or none")]


def _gate(probe, smiles, basis, spin, charge, require_stereo):
    from rdkit import Chem
    verdict, findings = check_structure(Chem.RemoveHs(probe), basis, spin=spin,
                                        maps_from=probe, smiles=smiles, charge=charge,
                                        require_stereo=require_stereo)
    return verdict, findings


def _identity(mol):
    """Canonical isomeric SMILES (the check) and InChIKey (a label only:
    standard InChI merges mobile-H tautomers, test_E8)."""
    from rdkit import Chem
    m = Chem.RemoveHs(Chem.Mol(mol))
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return dict(smiles=Chem.MolToSmiles(m), inchikey=Chem.MolToInchiKey(m))


def supplied_members(path, charge=None, basis="def2-tzvp", spin=0):
    """Read every record of a multi-record file through read_structure (the
    one reader) and check that they are conformers of ONE molecule in ONE
    configuration with ONE atom order (spec section 1). Record 1 is the
    approval. Returns (probe, [mol]); raises EnsembleRefused with every
    problem found, not only the first.

    Per member: a real 3D geometry with every hydrogen supplied (adding
    hydrogens would be generating geometry), verify_computed against the
    approval (another molecule, isomer or tautomer fails here), the same
    stereo completion as record 1 assigned from each member's own 3D
    coordinates (verify_computed alone accepts SOME completion), and
    verify_mapping index by index: a reordered record is refused, never
    trusted by file order.

    Review 6: every record passes the gate (S2: only record 1 was gated, so
    the ethyl cation, radical and anion made one ensemble, and the outcome
    depended on which came first), has record 1's formal charge, radical
    count and heavy-atom maps (the SCF gets one charge and spin; a site named
    by map must name the same atom in every member), and keeps its own record
    number (S6). No trans ring bond, allene axis or biaryl axis may differ
    from record 1's (S1: RDKit's stereo model sees none of them)."""
    from rdkit import Chem
    q = charge or 0
    records, probs, n_found = read_structure(path, charge=q, record="all")
    if not records:
        raise EnsembleRefused(probs)
    first, probs = records[0]
    if first is None:
        raise EnsembleRefused(probs or [_ensemble_block(
            "MEMBER_UNREADABLE", "record 1 could not be read", "fix record 1")])
    if n_found < 2:
        raise EnsembleRefused([_ensemble_block(
            "ENSEMBLE_ONE_RECORD", f"{path} holds one record",
            "supply a multi-record SDF, or use --ensemble K to generate conformers")])
    members, out = [first], []
    for r, (m, p) in enumerate(records[1:], start=2):
        if m is None:
            out += [(s, c, msg if msg.startswith(f"record {r} ") else f"record {r}: {msg}", fix)
                    for s, c, msg, fix in p] or [
                _ensemble_block("MEMBER_UNREADABLE", f"record {r} could not be read",
                                f"fix or remove record {r}")]
        members.append(m)
    approval = Chem.RemoveHs(first)
    ref_smiles = None

    def state(m):
        return (Chem.GetFormalCharge(m), sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()),
                tuple(a.GetAtomMapNum() for a in m.GetAtoms() if a.GetAtomicNum() > 1))

    ref_state, ref_signs = state(first), None
    for r, m in enumerate(members, start=1):
        if m is None:
            continue
        if r > 1:
            verdict, found = check_structure(Chem.RemoveHs(m), basis, spin=spin, maps_from=m,
                                             charge=charge)
            blocks = [(s, c, f"record {r}: {msg}", fix) for s, c, msg, fix in found if s == BLOCK]
            if blocks:
                out += blocks
                continue
            got = state(m)
            for code, what, a, b in (("MEMBER_CHARGE", "formal charge", got[0], ref_state[0]),
                                     ("MEMBER_SPIN", "unpaired electrons", got[1], ref_state[1]),
                                     ("MEMBER_MAPS", "atom maps", got[2], ref_state[2])):
                if a != b:
                    out.append(_ensemble_block(
                        code, f"record {r}: {what} {a}, record 1 {b}",
                        f"remove record {r}: one ensemble is one species, its atoms "
                        f"named the same way in every record"))
            if got != ref_state:
                continue
        if not coordinates_fix_stereo(m):
            out.append(_ensemble_block(
                "MEMBER_NOT_3D", f"record {r} is not a 3D geometry",
                f"supply real 3D coordinates for record {r}, or remove it"))
            continue
        if Chem.AddHs(m).GetNumAtoms() != m.GetNumAtoms():
            out.append(_ensemble_block(
                "MEMBER_MISSING_H", f"record {r} lacks hydrogens",
                f"supply every hydrogen in record {r}: adding them would generate geometry"))
            continue
        spec = _mol_to_spec(m)
        ok, computed, why = verify_computed(spec, approval)
        if r == 1:
            ref_smiles = computed
        if not ok:
            code = "MEMBER_STEREO" if "stereo" in why else "MEMBER_IDENTITY"
            out.append(_ensemble_block(
                code, f"record {r} is not record 1's molecule: {why} "
                      f"(record 1 {ref_smiles}, record {r} {computed or 'no match'})",
                f"remove record {r}: one ensemble is one molecule, one tautomer, "
                f"one protonation state and one stereoisomer"))
            continue
        if computed != ref_smiles:
            out.append(_ensemble_block(
                "MEMBER_STEREO", f"record {r} is another stereoisomer ({computed}, "
                                 f"record 1 {ref_smiles})",
                f"remove record {r}, or compute each stereoisomer as its own ensemble"))
            continue
        mapped, why_m = verify_mapping(spec, first)
        if not mapped:
            out.append(_ensemble_block(
                "MEMBER_ATOM_ORDER", f"record {r}: {why_m}",
                f"write record {r} with record 1's atom order"))
            continue
        signs = configuration_signs(Chem.RemoveHs(m))
        if ref_signs is None:
            ref_signs = signs
        elif signs != ref_signs:
            diff = sorted(k for k in signs if signs[k] != ref_signs.get(k))
            out.append(_ensemble_block(
                "MEMBER_STEREO", f"record {r} is another stereoisomer: {', '.join(diff)} "
                                 f"differ from record 1 (RDKit does not see these)",
                f"remove record {r}, or compute each stereoisomer as its own ensemble"))
    if out:
        raise EnsembleRefused(out)
    return first, members


def plan_ensemble(smiles=None, path=None, ensemble=None, charge=None, spin=0,
                  basis="def2-tzvp", weight=None, optimize=False, n_embed=N_EMBED,
                  seed=0xC0FFEE, dedup_rmsd=DEDUP_RMSD, rule="rg_spread", max_iters=2000):
    """Everything an ensemble needs before its first SCF: parse, gate, then
    either generate and select (ensemble = K) or read and check the supplied
    members (ensemble = "supplied"). Raises EnsembleRefused before any
    compute; with no --force, because an ensemble record has no
    validity.input to carry a forced member.

    Returns a dict: source, identity (the approved canonical isomeric SMILES
    and its InChIKey), gate findings, protocol, and members, each with its
    spec (the coordinates PySCF will receive), its shape descriptors, and
    for a generated ensemble its embedding conformer id and force-field
    energy (ordering only, never a weight)."""
    from rdkit import Chem
    problems = _weight_problems(weight, optimize)
    supplied = ensemble == "supplied"
    if not supplied and not (isinstance(ensemble, int) and not isinstance(ensemble, bool)
                             and ensemble >= 1):
        problems.append(_ensemble_block(
            "ENSEMBLE_K", f"--ensemble must be a positive K or 'supplied', not {ensemble!r}",
            "give the number of conformers to compute"))
    if supplied and not path:
        problems.append(_ensemble_block(
            "ENSEMBLE_SOURCE", "a supplied ensemble is read from a multi-record file",
            "pass the conformers as one SDF, or use --ensemble K with --smiles"))
    if not supplied and isinstance(ensemble, int) and ensemble > n_embed:
        problems.append(_ensemble_block(
            "ENSEMBLE_K", f"K = {ensemble} is more than the {n_embed} conformers embedded",
            "lower K, or embed more conformers"))
    if problems:
        raise EnsembleRefused(problems)
    q = charge or 0
    probe, probs = parse_for_check(smiles=smiles, path=path, charge=q,
                                   record=1 if supplied else None)
    if probe is None:
        raise EnsembleRefused(probs)
    # one stereoisomer per ensemble: a generated one must have nothing open,
    # or ETKDG picks a configuration per embedding (E13)
    verdict, findings = _gate(probe, smiles, basis, spin, charge, not supplied)
    if verdict == BLOCK:
        raise EnsembleRefused([f for f in findings if f[0] == BLOCK])
    if supplied:
        first, mols = supplied_members(path, charge=charge, basis=basis, spin=spin)
        members = [dict(record=r, spec=_mol_to_spec(m), **shape_descriptors(m))
                   for r, m in enumerate(mols, start=1)]
        protocol = dict(source="supplied", path=str(path), k_requested=len(mols),
                        k_obtained=len(mols))
        identity = _identity(first)
    else:
        mol = Chem.AddHs(Chem.RemoveHs(probe))
        mol.RemoveAllConformers()
        _, ff, spread, n_emb, n_conv, converged = _embed_and_minimise(
            mol, n_embed, seed, max_iters=max_iters)
        # Review 6 S1: the embedder works from the graph, which carries no
        # trans-ring helicity, allene axis or biaryl axis. A 3D input fixes
        # them: keep only the conformers in the file's configuration. From a
        # SMILES, --require-stereo has already refused any such element.
        heavy = Chem.RemoveHs(mol)
        ref = (configuration_signs(Chem.RemoveHs(probe))
               if coordinates_fix_stereo(probe) else None)
        n_other = 0
        if ref:
            keep = [(c, e) for c, e in converged if configuration_signs(heavy, c) == ref]
            n_other = len(converged) - len(keep)
            if not keep:
                raise EnsembleRefused([_ensemble_block(
                    "ENSEMBLE_CONFIGURATION",
                    f"no embedded conformer has the input's configuration {ref}",
                    "supply the conformers (--ensemble supplied)")])
            converged = keep
        try:
            sel = select_ensemble(mol, converged, ensemble, dedup_rmsd=dedup_rmsd, rule=rule)
        except ValueError as e:
            raise EnsembleRefused([_ensemble_block(
                "SYMMETRY_TOO_LARGE", str(e),
                "supply the conformers (--ensemble supplied)")]) from None
        signs = [configuration_signs(heavy, c) for c, _ in sel["members"]]
        if any(s != signs[0] for s in signs):
            raise EnsembleRefused([_ensemble_block(
                "MEMBER_STEREO", f"the selected members differ in {signs}",
                "fix every stereo element in the input")])
        approval = Chem.RemoveHs(probe)
        members, bad = [], []
        for cid, e in sel["members"]:
            spec = _mol_to_spec(mol, cid)
            ok, computed, why = verify_computed(spec, approval)
            mapped, why_m = verify_mapping(spec, approval)
            if not (ok and mapped):
                bad.append(_ensemble_block(
                    "MEMBER_IDENTITY", f"embedded conformer {cid}: {why or why_m}",
                    "report this: the embedder changed the approved structure"))
            members.append(dict(embed_conf_id=cid, mmff_energy=e, spec=spec,
                                **shape_descriptors(mol, cid)))
        if bad:
            raise EnsembleRefused(bad)
        protocol = dict(source="generated", n_embed=n_embed, n_embedded=n_emb,
                        seed=seed, forcefield=ff,
                        forcefield_spread=spread, n_other_configuration=n_other,
                        **{k: v for k, v in sel.items() if k != "members"})
        identity = _identity(probe)
    identity.update(charge=Chem.GetFormalCharge(probe),
                    n_radical=int(sum(a.GetNumRadicalElectrons() for a in probe.GetAtoms())))
    return dict(identity=identity, source=protocol["source"], protocol=protocol,
                gate=findings, members=members,
                folded_count=sum(m["folded"] for m in members),
                weight=weight or "none")


def estimate_cost(spec, basis, charge=0, spin=0):
    """nao without running an SCF. gto.M is cheap; the SCF is not.

    spin is required: gto.M refuses an electron count inconsistent with it,
    so omitting it made the cost guard crash on every open-shell molecule.
    """
    from pyscf import gto
    m = gto.M(atom=spec, basis=basis, charge=charge, spin=spin,
              unit="Angstrom", verbose=0)
    return m.natm, m.nao


# ---------------------------------------------------------------- electronic
def total_dm(mf):
    """Total (alpha+beta) density matrix, whichever solver produced it.

    UKS make_rdm1() returns a (2, nao, nao) stack. Every consumer here wants
    the TOTAL density: rho on the grid and the ESP electronic term are both
    functions of rho_alpha + rho_beta.
    """
    dm = np.asarray(mf.make_rdm1())
    return dm[0] + dm[1] if dm.ndim == 3 else dm


def gpu_untested(spin=0, solvent=None, density_fit=False):
    """The parts of a --gpu run outside what bench/gpu_l4.json validated
    (closed shell, gas phase, density fitting). Empty means validated."""
    out = []
    if spin:
        out.append(f"an open shell (spin={spin}, UKS)")
    if solvent not in (None, "vacuum"):
        out.append(f"a PCM solvent ({solvent})")
    if not density_fit:
        out.append("exact integrals (no --density-fit)")
    return out


def run_scf(atom, basis="def2-tzvp", xc="b3lyp", charge=0, spin=0, verbose=0,
            solvent=None, density_fit=False, gpu=False):
    """Restricted for a closed shell, UNRESTRICTED when spin != 0.

    v3 built an RKS object whatever the spin. PySCF does not refuse: on
    O with spin=2 (nelec 5,3) RKS returned -74.998688 Ha against UKS
    -75.000759, i.e. a converged, plausible, physically meaningless number
    from a restricted wavefunction describing an open shell.

    density_fit is OFF by default and the error budget was measured without
    it. Measured since (bench/results.jsonl): it moves V_min/V_max by <= 0.02
    kcal/mol across 7 pairs up to 898 nao, always upward on V_max.

    gpu runs the SCF ONLY on GPU4PySCF and hands back a CPU object, so
    total_dm and everything downstream are unchanged. Validated once
    (bench/gpu_l4.json, at commits 2313965 and 183ccd7): RKS + density_fit,
    gas phase, decanoic acid, where V_min/V_max matched the CPU to <0.001
    kcal/mol. UKS, PCM and runs without density fitting are untested on the
    GPU, and say so (review 3 H9: the contract was written and not enforced).
    """
    if gpu:
        for why in gpu_untested(spin, solvent, density_fit):
            print(f"WARNING: --gpu with {why} is untested; validated only for "
                  f"closed-shell, gas-phase, density-fitted runs (bench/gpu_l4.json). "
                  f"Check one result against a CPU run before relying on it.",
                  file=sys.stderr)
        try:
            import gpu4pyscf  # noqa: F401
        except ImportError:
            raise RuntimeError("gpu=True needs gpu4pyscf (pip install "
                               "gpu4pyscf-cuda12x) and an NVIDIA GPU") from None
    from pyscf import dft, gto
    mol = gto.M(atom=atom, basis=basis, charge=charge, spin=spin,
                unit="Angstrom", verbose=verbose)
    mf = dft.UKS(mol) if spin else dft.RKS(mol)
    if density_fit:
        mf = mf.density_fit()
    mf.xc = xc
    if solvent is not None:
        eps = SOLVENTS[solvent] if isinstance(solvent, str) else float(solvent)
        if eps is not None:
            mf = mf.PCM()                 # polarises the solute density
            mf.with_solvent.eps = eps
    mf.grids.level = 3
    if gpu:
        mf = mf.to_gpu()
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge ({basis}/{xc}"
                           f"{', spin=%d' % spin if spin else ''})")
    if gpu:
        mf = mf.to_cpu()    # np.asarray in total_dm refuses CuPy arrays
    if basis_key(basis) in SMALL_BASIS:
        print(f"WARNING: {basis} is double-zeta or smaller. Measured error on "
              f"the chlorobenzene sigma-hole V_max is +3.5 kcal/mol relative to "
              f"def2-TZVPP.\n         Adequate for a picture; do not quote a "
              f"number from it.", file=sys.stderr)
    return mol, mf


def density_on_grid(mol, dm, origin, n, spacing, block=20000):
    from pyscf.dft import numint
    ax = [origin[i] + np.arange(n[i]) * spacing for i in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    rho = np.empty(len(pts))
    for s in range(0, len(pts), block):
        ao = numint.eval_ao(mol, pts[s:s + block], deriv=0)
        rho[s:s + block] = numint.eval_rho(mol, ao, dm, xctype="LDA")
    return rho.reshape(n)


def esp_at_points(mol, dm, pts, mem_budget_gb=0.5):
    """V(r) in Hartree/e at arbitrary points (bohr).

        V(r) = sum_A Z_A/|r-R_A|  -  integral rho(r')/|r-r'| d3r'

    int1e_grids returns the POSITIVE Coulomb kernel <u|1/|r-r'||v>, so the
    electronic term is SUBTRACTED. Getting this backwards produces a map that
    looks entirely plausible and is sign-inverted; --selftest asserts it.

    F5: chunk is sized from nao, not fixed. The allocation is
    (chunk, nao, nao) float64, quadratic in basis size.
    """
    pts = np.asarray(pts, float)
    V = np.zeros(len(pts))
    for A in range(mol.natm):
        d = np.linalg.norm(pts - mol.atom_coord(A), axis=1)
        V += mol.atom_charge(A) / np.maximum(d, 1e-12)

    chunk = int(mem_budget_gb * 1e9 / (8 * mol.nao ** 2))
    chunk = max(16, min(chunk, 2000))
    for s in range(0, len(pts), chunk):
        ints = mol.intor("int1e_grids", grids=pts[s:s + chunk])
        V[s:s + chunk] -= np.einsum("pij,ij->p", ints, dm)
    return V


class ALIEUnbound(ValueError):
    """An occupied orbital at or above zero energy: the system is unbound in
    this model (a gas-phase anion), and no ionisation energy exists to
    average. Code ALIE_UNBOUND."""
    code = "ALIE_UNBOUND"


def alie_at_points(mol, mf, pts, mem_budget_gb=0.5):
    """Average local ionisation energy (Hartree) at arbitrary points (bohr).

        I(r) = sum_i n_i |psi_i(r)|^2 (-e_i)  /  sum_i n_i |psi_i(r)|^2

    over occupied orbitals (Sjoberg et al. 1990; docs/p4-spec.md section 1),
    evaluated at the surface vertices like esp_at_points, never on a grid.
    RKS weights each spatial orbital by its occupation (2); UKS sums alpha and
    beta separately with occupation 1, so the denominator is rho_a + rho_b,
    as total_dm gives (invariant 5). Orbital energies are Kohn-Sham, so I
    depends on the functional far more than V does.

    -e_i, not |e_i| as Sjoberg 1990 writes it: the two agree for a bound
    system, and |e| turned a gas-phase anion's positive HOMO into a
    plausible ALIE (enolate, def2-SVP: 2.19 eV, review 6 S9). Any occupied
    orbital at or above zero raises ALIEUnbound instead.

    The one-orbital identity -- a single occupied orbital gives I = -e_1
    everywhere, whatever the basis -- is asserted in --selftest.
    """
    from pyscf.dft import numint
    pts = np.asarray(pts, float)
    C, e, n = (np.asarray(x, dtype=float) for x in (mf.mo_coeff, mf.mo_energy, mf.mo_occ))
    if C.ndim == 2:                     # RKS: one spin channel
        C, e, n = C[None], e[None], n[None]
    top = max(float(es[ns > 0].max()) for es, ns in zip(e, n) if (ns > 0).any())
    if top >= 0:
        raise ALIEUnbound(f"highest occupied orbital at {top * HA2EV:+.2f} eV: unbound "
                          f"in this model; compute it in solvent, or report no ALIE")
    chunk = max(16, min(int(mem_budget_gb * 1e9 / (8 * 2 * mol.nao)), 20000))
    out = np.empty(len(pts))
    for s in range(0, len(pts), chunk):
        ao = numint.eval_ao(mol, pts[s:s + chunk], deriv=0)
        num = np.zeros(len(ao))
        den = np.zeros(len(ao))
        for Cs, es, ns in zip(C, e, n):
            occ = ns > 0
            w = (ao @ Cs[:, occ]) ** 2 * ns[occ]    # n_i |psi_i|^2 at each point
            num += w @ -es[occ]
            den += w.sum(axis=1)
        # NaN where no occupied orbital reaches the point (den underflows to 0,
        # e.g. 25 A out); never on a rho = 0.002 vertex
        out[s:s + chunk] = np.divide(num, den, out=np.full(len(num), np.nan), where=den > 0)
    return out


def reaction_field_at_points(mf, pts):
    """Potential (Hartree/e) of the PCM surface charges at points (bohr).

    esp_at_points gives the potential of the PCM-POLARISED SOLUTE; the
    solvent's own response, the apparent surface charges on the cavity, is a
    second term. Review 3 H7: on the chlorobenzene sigma-hole the two
    conventions give solvent shifts of opposite sign, and only one was stated.
    PySCF's charges are Gaussian-smeared (exponent charge_exp; the
    switching/Gaussian scheme of Lange & Herbert 2010 and York & Karplus 1999,
    as cited in pyscf/solvent/pcm.py), so each contributes q erf(zeta r) / r:
    finite at a vertex that sits near the cavity, where a bare point charge
    would blow up.
    """
    from scipy.special import erf
    ws = mf.with_solvent
    q = np.asarray(ws._intermediates["q"]).ravel()
    centres = ws.surface["grid_coords"]
    zeta = ws.surface["charge_exp"]
    pts = np.asarray(pts, float)
    V = np.zeros(len(pts))
    for s in range(0, len(pts), 2000):
        d = np.linalg.norm(pts[s:s + 2000, None, :] - centres[None], axis=2)
        with np.errstate(invalid="ignore", divide="ignore"):
            k = np.where(d > 1e-12, erf(zeta * d) / d, 2.0 * zeta / np.sqrt(np.pi))
        V[s:s + 2000] = k @ q
    return V


HA2EV = 27.211386


def frontier_orbitals_eV(mf):
    """(HOMO, LUMO) in eV; for UKS the highest occupied and lowest virtual
    over both spins. LUMO is None when the basis has no virtual orbital."""
    e = np.atleast_2d(np.asarray(mf.mo_energy, dtype=float))
    occ = np.atleast_2d(np.asarray(mf.mo_occ, dtype=float))
    virt = e[occ == 0]
    return (float(e[occ > 0].max()) * HA2EV,
            float(virt.min()) * HA2EV if virt.size else None)


def unbound_anion_warning(mf, charge, solvent):
    """A gas-phase anion whose HOMO is above zero has an electron the basis
    is holding in place: in a complete basis it would leave. Its surface
    potential then depends on how diffuse the basis is rather than on the
    chemistry. A continuum binds it, so name --solvent as the remedy."""
    if charge >= 0 or solvent:
        return None
    homo = frontier_orbitals_eV(mf)[0]
    if homo <= 0:
        return None
    return (f"gas-phase anion with HOMO {homo:+.2f} eV: the outer "
            f"electron is unbound and held only by the basis, so the map tracks "
            f"the basis set, not the molecule. Compute it in a continuum "
            f"(--solvent water) or as the neutral/protonated form.")


def electronic_report(mf, charge, solvent):
    """results.electronic: frontier orbitals and the problems that make the
    electronic structure, and so the map, not the molecule's. Review 5 R3:
    the unbound-anion condition was a stderr warning, and the record of a
    gas-phase acetate at def2-SVP (HOMO +0.17 eV) said validity all true."""
    homo, lumo = frontier_orbitals_eV(mf)
    warn = unbound_anion_warning(mf, charge, solvent)
    # each problem leads with a stable code, so a consumer branches on the
    # code, not the prose (as gate findings do)
    return dict(homo_eV=homo, lumo_eV=lumo,
                problems=[f"UNBOUND_ANION_GAS: {warn}"] if warn else [])


# ------------------------------------------------------------------- surface
def euler_characteristic(faces):
    e = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [0, 2]]])
    n_edges = len(np.unique(np.sort(e, axis=1), axis=0))
    return len(np.unique(faces)) - n_edges + len(faces)


def mesh_topology(faces):
    """What the mesh IS, measured directly, rather than inferred from chi alone.

    chi = 2 was the clipping detector, but chi conflates two things: an OPEN
    mesh (clipped by the box: edges used by one triangle) and a CLOSED mesh of
    higher genus (a handle, e.g. a folded chain touching itself at the
    isovalue). Two of ten C18 acetate conformers (bench/conformers.jsonl, the
    two most compact) were refused at chi = -2 and -4 although nothing was
    clipped. Counting edge use separates the cases: every edge of a closed
    manifold is shared by exactly two triangles.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    e = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [0, 2]]]), axis=1)
    _, use = np.unique(e, axis=0, return_counts=True)
    used = np.unique(faces)
    remap = np.full(faces.max() + 1, -1)
    remap[used] = np.arange(len(used))
    f = remap[faces]
    a = np.concatenate([f[:, 0], f[:, 1], f[:, 0]])
    b = np.concatenate([f[:, 1], f[:, 2], f[:, 2]])
    adj = coo_matrix((np.ones(len(a)), (a, b)), shape=(len(used),) * 2)
    n_comp = int(connected_components(adj, directed=False)[0])
    chi = int(len(used) - len(use) + len(faces))
    open_e, nonman = int((use == 1).sum()), int((use > 2).sum())
    closed = open_e == 0 and nonman == 0
    # closed orientable surfaces: chi = 2 * components - 2 * total genus
    genus = (2 * n_comp - chi) // 2 if closed else None
    return dict(chi=chi, components=n_comp, open_edges=open_e,
                nonmanifold_edges=nonman, closed=closed, genus=genus)


def surface_problems(topo, face_max_rho, iso):
    """(problems, notes): problems refuse the surface, notes only inform."""
    problems, notes = [], []
    if face_max_rho > iso:
        problems.append(f"density on a grid-box face is {face_max_rho:.2e} > iso={iso}; "
                        f"the isosurface is CLIPPED (raise --margin)")
    if topo["open_edges"]:
        problems.append(f"the mesh is OPEN ({topo['open_edges']} edges belong to one "
                        f"triangle only): not a single closed surface")
    if topo["nonmanifold_edges"]:
        problems.append(f"{topo['nonmanifold_edges']} edges are shared by more than two "
                        f"triangles: not a single closed surface")
    if topo["closed"] and topo["genus"]:
        notes.append(f"closed surface of genus {topo['genus']} (chi = {topo['chi']}): the "
                     f"molecule folds back on itself at the isovalue. Valid; the surface "
                     f"has {topo['genus']} handle(s), so the 3D view may show a tunnel.")
    if topo["closed"] and topo["components"] > 1:
        notes.append(f"WARNING: {topo['components']} separate closed surfaces (an enclosed "
                     f"void, or fragments the gate should have caught). Areas and extrema "
                     f"include all of them.")
    return problems, notes


def extract_surface(mol, dm, spacing, iso=ISO_BADER, margin=5.0, strict=True):
    """F4: a clipped isosurface is the silent failure of this pipeline.

    At margin 3.5 bohr the review measured the chlorobenzene sigma-hole at
    +4.69 instead of +6.35 kcal/mol -- a 26% error -- while every other number
    in the run stayed plausible. Two cheap, independent detectors: density on
    the grid-box faces, and the topology of the mesh (open edges). chi = 2 is
    no longer required: a folded molecule can give a closed surface of higher
    genus, which is valid (see mesh_topology).
    """
    from skimage import measure
    co = mol.atom_coords()
    lo = co.min(0) - margin
    n = np.maximum(np.ceil((co.max(0) + margin - lo) / spacing).astype(int), 8)
    rho = density_on_grid(mol, dm, lo, tuple(n), spacing)

    face_max = max(rho[0].max(), rho[-1].max(), rho[:, 0].max(),
                   rho[:, -1].max(), rho[:, :, 0].max(), rho[:, :, -1].max())
    verts_i, faces, _, _ = measure.marching_cubes(rho, level=iso)
    verts = lo + verts_i * spacing
    faces = faces.astype(int)
    topo = mesh_topology(faces)

    problems, notes = surface_problems(topo, face_max, iso)
    diag = dict(grid=tuple(int(x) for x in n), n_verts=len(verts),
                face_max_rho=float(face_max), **topo,
                rho_integral=float(rho.sum() * spacing ** 3),
                nelec=int(mol.nelectron), margin=margin, spacing=spacing,
                strict=bool(strict), problems=list(problems))
    for n_ in notes:
        print(("" if n_.startswith("WARNING") else "note: ") + n_, file=sys.stderr)
    if problems:
        msg = "isosurface integrity failed:\n  - " + "\n  - ".join(problems)
        if strict:
            raise RuntimeError(msg + "\n(pass --no-strict to render anyway)")
        print("WARNING: " + msg, file=sys.stderr)
    return verts, faces, diag


# -------------------------------------------------------------------- render
def look_at(elev_deg, azim_deg):
    e, a = np.radians(elev_deg), np.radians(azim_deg)
    fwd = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    up0 = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(fwd, up0)) > 0.999:
        up0 = np.array([0.0, 1.0, 0.0])
    right = np.cross(up0, fwd)
    right /= np.linalg.norm(right)
    return np.stack([right, np.cross(fwd, right), fwd])


def render(verts, faces, scalars, cmap, vlim, size=700, elev=20, azim=172,
           ambient=0.62, diffuse=0.38, ss=3, pad=1.10,
           light_dir=(-0.35, 0.45, 0.82), bg=(1, 1, 1)):
    """Orthographic, Lambertian, no specular.

    ambient/diffuse default to a 0.62..1.0 shading span so geometry reads
    without swamping the colormap's own lightness. Specular is absent by
    design: a white highlight at an arbitrary location reads as a potential
    extremum.
    """
    M = look_at(elev, azim)
    P = verts @ M.T
    c = 0.5 * (P[:, :2].min(0) + P[:, :2].max(0))
    half = pad * 0.5 * (P[:, :2].max(0) - P[:, :2].min(0)).max()
    S = size * ss
    xy = (P[:, :2] - c) / (2 * half) * S + S / 2
    z = P[:, 2]

    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    fn = np.cross(v1 - v0, v2 - v0)
    vn = np.zeros_like(verts)
    for k in range(3):
        np.add.at(vn, faces[:, k], fn)
    nr = np.linalg.norm(vn, axis=1, keepdims=True)
    nr[nr == 0] = 1
    vn_s = (vn / nr) @ M.T
    # F8: orient PER VERTEX toward the camera. The v1 global flip keyed on
    # mean(vn_z), which on a closed surface is ~1e-3 of the per-vertex scale,
    # i.e. numerical noise. (Measured immaterial for closed meshes; it is the
    # open-mesh case that made it worth removing.)
    vn_s = np.where((vn_s[:, 2] > 0)[:, None], -vn_s, vn_s)

    t_v = np.clip((scalars - vlim[0]) / (vlim[1] - vlim[0]), 0, 1)
    L = np.array(light_dir, float)
    L /= np.linalg.norm(L)
    ndl = np.clip(-(vn_s @ np.array([L[0], L[1], -abs(L[2])])), 0, 1)
    sh_v = ambient + diffuse * ndl
    rgb_v = cmap(t_v)[:, :3]

    img = np.ones((S, S, 3)) * np.array(bg, float)
    zb = np.full((S, S), np.inf)

    e1 = xy[faces[:, 1]] - xy[faces[:, 0]]
    e2 = xy[faces[:, 2]] - xy[faces[:, 0]]
    fnz = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]   # numpy>=2 has no 2-D cross
    for fi in np.where(np.abs(fnz) > 1e-12)[0]:
        i0, i1, i2 = faces[fi]
        x0, y0 = xy[i0]
        x1, y1 = xy[i1]
        x2, y2 = xy[i2]
        xmin = max(int(np.floor(min(x0, x1, x2))), 0)
        xmax = min(int(np.ceil(max(x0, x1, x2))), S - 1)
        ymin = max(int(np.floor(min(y0, y1, y2))), 0)
        ymax = min(int(np.ceil(max(y0, y1, y2))), S - 1)
        if xmin > xmax or ymin > ymax:
            continue
        den = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        if abs(den) < 1e-12:
            continue
        PX, PY = np.meshgrid(np.arange(xmin, xmax + 1), np.arange(ymin, ymax + 1))
        w0 = ((y1 - y2) * (PX - x2) + (x2 - x1) * (PY - y2)) / den
        w1 = ((y2 - y0) * (PX - x2) + (x0 - x2) * (PY - y2)) / den
        w2 = 1 - w0 - w1
        m = (w0 >= -1e-9) & (w1 >= -1e-9) & (w2 >= -1e-9)
        if not m.any():
            continue
        zz = w0 * z[i0] + w1 * z[i1] + w2 * z[i2]
        sl = (slice(ymin, ymax + 1), slice(xmin, xmax + 1))
        hit = m & (zz < zb[sl])
        if not hit.any():
            continue
        col = (w0[..., None] * rgb_v[i0] + w1[..., None] * rgb_v[i1]
               + w2[..., None] * rgb_v[i2])
        sh = (w0 * sh_v[i0] + w1 * sh_v[i1] + w2 * sh_v[i2])[..., None]
        img[sl][hit] = np.clip(col * sh, 0, 1)[hit]
        zb[sl][hit] = zz[hit]
    return np.flipud(img.reshape(size, ss, size, ss, 3).mean(axis=(1, 3)))


# ------------------------------------------------------------------ colormaps
def build_maps():
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap as LSC
    m = {
        "rainbow": plt.get_cmap("jet").reversed(),
        "rgb3": LSC.from_list("rgb3", ["#e00000", "#00c000", "#0000e0"]),
        "RdBu": plt.get_cmap("RdBu"),
        "PuOr": plt.get_cmap("PuOr"),
        "pindr": LSC.from_list("pindr", ["#c23a35", "#ea8b70", "#ddd8cd",
                                         "#7f8ddc", "#2340c8"]),
    }
    try:
        from cmcrameri import cm as crm
        m["vik"] = crm.vik.reversed()
        m["broc"] = crm.broc.reversed()
        m["cork"] = crm.cork.reversed()
        # sequential, for ions: no sign change to straddle, so a diverging map
        # would spend half its ramp on values the molecule never takes
        # NOT reversed: lajolla's t=0 is already its dark end, so V_min (most
        # electron-rich) renders darkest, keeping the intense = electron-rich
        # reading the diverging maps use. Verified against benzoate, where
        # V_min sits on the carboxylate O and V_max on a ring H.
        m["lajolla"] = crm.lajolla
        m["batlow"] = crm.batlow
    except ImportError:
        print("note: cmcrameri absent; vik/broc/cork/lajolla/batlow unavailable",
              file=sys.stderr)
    m["viridis"] = plt.get_cmap("viridis")
    return m


SEQUENTIAL = {"lajolla", "batlow", "viridis"}


# The minority sign must reach this fraction of the larger extremum before a
# surface counts as changing sign (review 3 H11).
SIGN_CHANGE_MIN = 0.05


def sign_change_ratio(V):
    """min(|V_min|, |V_max|) / max(...) if V spans zero, else 0.

    Review 3 H11: `V.min() < 0 < V.max()` was a knife-edge. A cation surface
    from -0.01 to +120 got a diverging +-120 scale using half the ramp, and
    from +0.01 to +120 a sequential one using all of it: 0.02 kcal/mol, far
    below the 3.5 kcal/mol basis error, decided the picture.
    """
    lo, hi = float(np.min(V)), float(np.max(V))
    if not lo < 0.0 < hi:
        return 0.0
    return min(-lo, hi) / max(-lo, hi)


def choose_scale(V, lim=None, kind="auto", vmin=None, vmax=None):
    """Decide the colour scale from the data, not from habit.

    A diverging scale symmetric about zero is right only if the surface
    actually CHANGES SIGN. For an ion it does not: benzoate runs -158.7 to
    -40.7 kcal/mol, so centring on zero maps the whole molecule to one end
    and hides 118 kcal/mol of structure. "Changes sign" means the minority
    sign reaches SIGN_CHANGE_MIN of the larger extremum (sign_change_ratio).

    Returns (vlim, kind, notes).
    """
    data_min, data_max = float(V.min()), float(V.max())
    ratio = sign_change_ratio(V)
    spans_zero = ratio > SIGN_CHANGE_MIN
    sign_txt = (f"minority sign {ratio:.3f} of the larger extremum, threshold "
                f"{SIGN_CHANGE_MIN}")
    notes = []
    if kind == "auto":
        kind = "diverging" if spans_zero else "sequential"
        if not spans_zero:
            notes.append(f"surface does not change sign appreciably ({data_min:+.1f} to "
                         f"{data_max:+.1f} kcal/mol; {sign_txt}): using a SEQUENTIAL "
                         f"scale over the actual range. A diverging scale centred on "
                         f"zero would put the whole molecule at one end of the ramp.")
    if kind == "diverging":
        if not spans_zero:
            notes.append(f"WARNING: diverging scale requested but V barely or never "
                         f"changes sign ({data_min:+.1f} to {data_max:+.1f}; {sign_txt}); "
                         f"nearly every vertex falls on one side of the midpoint colour.")
        hi = abs(lim) if lim else max(abs(data_min), abs(data_max))
        vlim = (-hi, hi)
        asym = (max(abs(data_min), abs(data_max))
                / max(min(abs(data_min), abs(data_max)), 1e-9))
        if lim is None and asym > 2.5:
            notes.append(f"range is asymmetric ({data_min:+.1f} / {data_max:+.1f}, ratio "
                         f"{asym:.1f}): symmetric limits leave part of the ramp "
                         f"unused. That is the correct trade for a signed quantity.")
    else:
        # v3 bug: --lim on a sequential scale re-imposed (-lim, +lim), i.e.
        # exactly the symmetric range this branch exists to avoid. For an
        # anion running -162..-44 that put half the ramp outside the data.
        if lim is not None:
            notes.append(f"--lim is a SYMMETRIC range for a diverging scale and "
                         f"does not apply here; ignoring it. Use --vmin/--vmax "
                         f"to set an explicit sequential range.")
        lo = data_min if vmin is None else vmin
        hi = data_max if vmax is None else vmax
        if lo >= hi:
            raise SystemExit(f"--vmin {lo} must be below --vmax {hi}")
        vlim = (lo, hi)
    return vlim, kind, notes


DEFAULT_CMAP = {"diverging": "vik", "sequential": "lajolla"}


def resolve_cmap(kind, requested=None, force=False):
    """The colormap actually used, chosen AFTER the scale kind is known.

    Review 3 H3: --cmap defaulted to 'vik' whatever the scale, so an anion's
    sequential scale was drawn with a diverging map: acetate's -80 kcal/mol end
    came out dark BLUE, the colour every caption says means positive, and the
    map's white neutral sat at -125. A note was printed and the image shipped.
    The reverse mismatch is as bad: a sequential map on a sign-changing field
    hides the sign. Both now refuse unless force (--force-cmap) is given.
    """
    if requested is None:
        return DEFAULT_CMAP[kind]
    if not force and (requested in SEQUENTIAL) != (kind == "sequential"):
        want = "a sequential map (lajolla, batlow, viridis)" if kind == "sequential" \
            else "a diverging map (vik, RdBu, PuOr, broc, cork, pindr)"
        raise ValueError(f"--cmap {requested} does not fit a {kind} scale: use {want}, "
                         f"leave --cmap out to choose automatically, or pass "
                         f"--force-cmap to draw it anyway")
    return requested


def scale_caption(kind):
    """How to read the colours, for every caption: one source, so the 3D
    colourbar and the notebook image cannot disagree (the image kept saying
    'warm = negative' on ions after the colourbar was fixed)."""
    return "warm = negative" if kind == "diverging" else "darker = more negative"


def cmap_after_scf(kind, requested=None, force=False):
    """(cmap, warning or None) once the scale kind is known from the data.

    With --scale auto the kind is only known after the SCF, which can take 15
    minutes on a C18 molecule. Refusing then would throw that work away with
    no image, JSON or NPZ written, so a mismatch falls back to the map that
    fits, says so, and the JSON records both. Where the kind IS known up front
    (--scale given), main refuses before the SCF instead (cmap_preflight).
    """
    try:
        return resolve_cmap(kind, requested, force), None
    except ValueError as e:
        fallback = DEFAULT_CMAP[kind]
        return fallback, (f"WARNING: {e}. The surface turned out {kind}, which was "
                          f"only known after the SCF, so it is drawn with "
                          f"{fallback} instead (recorded in --json).")


def cmap_preflight(scale, requested=None, force=False, available=None):
    """Every colour-map check that does not need the SCF, run before it.
    Returns an error message, or None."""
    if requested is not None and available is not None and requested not in available:
        return f"--cmap must be one of {sorted(available)}"
    if requested is not None and scale in ("diverging", "sequential"):
        try:
            resolve_cmap(scale, requested, force)
        except ValueError as e:
            return str(e)
    return None


# ------------------------------------------------------------------ molecules
MOLECULES = {"bromobenzene": ("Br", 1.890), "chlorobenzene": ("Cl", 1.740),
             "iodobenzene": ("I", 2.095), "fluorobenzene": ("F", 1.354),
             "benzene": ("H", 1.080)}


def surface_definition(iso):
    """Every area and value in a record is on this surface. It is the
    electron-density isosurface, NOT a solvent-accessible surface."""
    return dict(property="electron_density", isovalue_e_bohr3=float(iso))


def surface_report(diag, iso):
    """The surface-integrity block of --json.

    Review 3 H4: the JSON kept only chi, and chi = -2 is what a clipped surface
    run with --no-strict gives AND what a valid folded genus-2 surface gives.
    A downstream reader could not tell a 23 %-wrong V_min from a good one.
    """
    return dict(definition=surface_definition(iso),
                closed=bool(diag["closed"]), open_edges=int(diag["open_edges"]),
                nonmanifold_edges=int(diag["nonmanifold_edges"]),
                components=int(diag["components"]), genus=diag["genus"],
                chi=int(diag["chi"]), face_max_rho=float(diag["face_max_rho"]),
                margin_bohr=float(diag["margin"]), strict=bool(diag["strict"]),
                problems=list(diag["problems"]), valid=not diag["problems"])


def gate_report(verdict, findings, forced):
    """The structure-gate block of --json: what the gate said, and whether it
    was overridden. Review 3 H4: a structure REJECTED by the gate and computed
    with --force was written out as `approved_smiles`."""
    return dict(verdict=verdict, forced=bool(forced and verdict == BLOCK),
                findings=[[sev, code] for sev, code, *_ in findings])


def scale_provenance(kind, lim=None, vmin=None, vmax=None):
    """Where the colour limits came from, and whether the user fixed them.

    Fixed means another molecule run with the same flags maps a colour to the
    same kcal/mol. It says nothing about the numbers: review 3 H6 found it
    named "comparable" across def2-SVP vs TZVP and gas vs PCM; whether numbers
    compare is analysis_key, whether pictures do is render_key. The first
    --json marked every sequential scale comparable; an auto-scaled one (acetate
    -170..-80 next to another anion at -130..-30) is not. --lim only fixes a
    DIVERGING scale; --vmin/--vmax only fix a SEQUENTIAL one, and only as a
    pair (choose_scale ignores the others for that kind).
    """
    if kind == "diverging":
        fixed = lim is not None
    else:
        fixed = vmin is not None and vmax is not None
    return dict(source="user" if fixed else "auto", colour_scale_fixed=fixed)


def plotly_colorscale(cmap, n=33):
    """matplotlib colormap -> plotly colorscale, sampled t = 0 .. 1 in order.

    Built from the SAME cmap object render() uses, so the orientation that
    --selftest asserts (warm = negative on diverging maps, V_min darkest on
    sequential ones) carries over instead of being re-derived for plotly.
    """
    out = []
    for t in np.linspace(0.0, 1.0, n):
        r, g, b = (int(round(255 * c)) for c in cmap(float(t))[:3])
        out.append([float(t), f"rgb({r},{g},{b})"])
    return out


# covalent radii (A) for drawing bonds in the 3D view; TESTED_ELEMENTS only
_RCOV = {"H": .31, "B": .84, "C": .76, "N": .71, "O": .66, "F": .57, "Si": 1.11,
         "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39, "Se": 1.20}
_ATOM_RGB = {"H": "#f2f2f2", "C": "#555555", "N": "#3050f8", "O": "#ff0d0d",
             "F": "#90e050", "Cl": "#1ff01f", "Br": "#a62929", "I": "#940094",
             "S": "#d9c21f", "P": "#ff8000"}


def esp_figure(verts, faces, V, symbols, coords, cmap, vlim, title="",
               kind="diverging"):
    """Interactive 3D surface (plotly): rotate, zoom, hover to read V.

    verts and coords in BOHR (as computed); shown in Angstrom. Lighting has
    no specular term for the same reason render() has none: a white
    highlight at an arbitrary place reads as a potential extremum. The
    surface starts opaque so its colours are exact; a button makes it
    see-through to show which atom a feature sits on.
    """
    import plotly.graph_objects as go
    P, X = np.asarray(verts) * BOHR, np.asarray(coords) * BOHR
    surf = go.Mesh3d(
        x=P[:, 0], y=P[:, 1], z=P[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        intensity=V, intensitymode="vertex",
        colorscale=plotly_colorscale(cmap), cmin=vlim[0], cmax=vlim[1],
        colorbar=dict(title="V (kcal/mol)<br>" + scale_caption(kind), len=0.75),
        lighting=dict(ambient=0.62, diffuse=0.38, specular=0.0,
                      roughness=1.0, fresnel=0.0),
        hovertemplate="V = %{intensity:+.1f} kcal/mol<extra></extra>",
        name="ESP surface", opacity=1.0)
    bx, by, bz = [], [], []
    for a in range(len(symbols)):
        for b in range(a + 1, len(symbols)):
            cut = 1.15 * (_RCOV.get(symbols[a], .8) + _RCOV.get(symbols[b], .8))
            if np.linalg.norm(X[a] - X[b]) < cut:
                for u in (X[a], X[b]):
                    bx.append(u[0]); by.append(u[1]); bz.append(u[2])
                bx.append(None); by.append(None); bz.append(None)
    bonds = go.Scatter3d(x=bx, y=by, z=bz, mode="lines", hoverinfo="skip",
                         line=dict(color="#444444", width=5), name="bonds")
    atoms = go.Scatter3d(
        x=X[:, 0], y=X[:, 1], z=X[:, 2], mode="markers", name="atoms",
        text=symbols, hovertemplate="%{text}<extra></extra>",
        marker=dict(size=[3 if s == "H" else 5 for s in symbols],
                    color=[_ATOM_RGB.get(s, "#ff1493") for s in symbols],
                    line=dict(color="#222222", width=1)))
    fig = go.Figure([surf, bonds, atoms])
    fig.update_layout(
        title=title, showlegend=False, margin=dict(l=0, r=0, t=40, b=0),
        scene=dict(aspectmode="data", xaxis_visible=False,
                   yaxis_visible=False, zaxis_visible=False),
        updatemenus=[dict(type="buttons", direction="right", x=0.0, y=1.0,
                          xanchor="left", yanchor="top", buttons=[
            dict(label="opaque", method="restyle", args=[{"opacity": 1.0}, [0]]),
            dict(label="see-through", method="restyle",
                 args=[{"opacity": 0.45}, [0]])])])
    return fig


def monosubstituted_benzene(sub,c_sub, rcc=1.390, rch=1.080):
    """Idealised planar C6H5X. NOT an optimised geometry.

    Measured dV_max/dR ~ -14 kcal/mol per Angstrom for the chlorobenzene
    sigma-hole, so a 0.01 A geometry error costs ~0.15 kcal/mol. The built-in
    C-Cl 1.740 is within 0.001 A of the microwave value, which is luck rather
    than method -- pass --atoms with an optimised geometry for real work.
    """
    out = [("C", (rcc * np.cos(np.radians(60 * i)),
                  rcc * np.sin(np.radians(60 * i)), 0.0)) for i in range(6)]
    out.append((sub, (rcc + c_sub, 0.0, 0.0)))
    out += [("H", ((rcc + rch) * np.cos(np.radians(60 * i)),
                   (rcc + rch) * np.sin(np.radians(60 * i)), 0.0))
            for i in range(1, 6)]
    return "; ".join(f"{s} {a:.4f} {b:.4f} {c:.4f}" for s, (a, b, c) in out)


def sigma_axis_partner(coords, i_sub):
    """The atom the sigma-hole axis runs from: the substituent's nearest atom.
    One function so the JSON's i_bond is the atom sigma_hole actually used."""
    d = np.linalg.norm(coords - coords[i_sub], axis=1)
    d[i_sub] = np.inf
    return int(d.argmin())


def sigma_hole_record(verts, V, mol, i_sub, reason=None):
    """results.sigma_hole for the JSON (review 3 H6: the number every
    error-budget row is about was printed and never written). Returns
    (record or None, reason or None)."""
    if i_sub is None:
        return None, reason or "no Cl/Br/I"
    i_bond = sigma_axis_partner(mol.atom_coords(), i_sub)
    sh, ns = sigma_hole(verts, V, mol, i_sub=i_sub, i_bond=i_bond)
    return dict(V_max_kcal_mol=None if np.isnan(sh) else sh, n_cap_vertices=ns,
                i_sub=int(i_sub), i_bond=i_bond), None


def dipole_about_nuclear_centre(mol, dm):
    """Dipole (Debye) about the centre of nuclear charge, and that centre (A).

    Review 3 H6: PySCF's dip_moment takes the origin at (0,0,0), and for a
    charged species the dipole depends on the origin -- acetate gave 6.42 D as
    embedded and 30.35 D translated by 5 A, same energy. About the centre of
    nuclear charge the nuclear term vanishes and the value moves with the
    molecule. For a neutral molecule it equals the origin-free dipole.
    """
    Z = mol.atom_charges().astype(float)
    centre = Z @ mol.atom_coords() / Z.sum()
    with mol.with_common_orig(centre):
        r = mol.intor_symmetric("int1e_r", comp=3)
    mu_au = -np.einsum("xij,ji->x", r, np.asarray(dm))
    mu = mu_au * AU2DEBYE
    return dict(vector_D=[float(x) for x in mu], magnitude_D=float(np.linalg.norm(mu)),
                origin=dict(kind="centre_of_nuclear_charge",
                            xyz_angstrom=[float(x) for x in centre * BOHR]))


def record_key(fields):
    """sha1 of the canonical JSON of fields: equal keys, equal fields."""
    import hashlib
    import json
    return hashlib.sha1(json.dumps(fields, sort_keys=True).encode()).hexdigest()


def geometry_record(mol, charge, spin):
    """The coordinates actually computed, in PySCF atom order (review 3 H14:
    the JSON had none, and ETKDG output depends on the RDKit version, so a
    seed is not a substitute). Full float precision, so feeding the record
    back through --atoms reproduces V; the sha1 is taken at %.8f so it does
    not depend on how a float happens to be printed.

    The sha1 covers the nuclei only (symbols and positions), not charge or
    spin: an anion and a radical on the same nuclei share it, and differ in
    the calculation. No atom maps here; the site block carries those."""
    syms = [mol.atom_symbol(i) for i in range(mol.natm)]
    xyz = mol.atom_coords() * BOHR
    canon = "\n".join(f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, (x, y, z) in zip(syms, xyz))
    return dict(symbols=syms, xyz_angstrom=[[float(c) for c in r] for r in xyz],
                sha1=_sha1(canon), sha1_covers=GEOMETRY_SHA1_COVERS,
                charge=int(charge), spin=int(spin))


GEOMETRY_SHA1_COVERS = ("symbols and xyz_angstrom at %.8f, one 'Sym x y z' line "
                        "per atom, newline-joined; not charge or spin")


def atoms_from_geometry(geom):
    """--atoms string for a JSON geometry block: the way to recompute a record."""
    return "; ".join(f"{s} {x!r} {y!r} {z!r}"
                     for s, (x, y, z) in zip(geom["symbols"], geom["xyz_angstrom"]))


def reproduce_argv(rec):
    """Command-line arguments that recompute a JSON record from its own
    coordinates under the same analysis_key (review 5 R6): the atoms, the
    settings the key hashes, and the record's geometry_method declared with
    --geometry-origin, since --atoms alone is keyed "undeclared"."""
    s, inp = rec["settings"], rec["input"]
    argv = ["--atoms", atoms_from_geometry(rec["geometry"]),
            "--charge", str(inp["charge"]), "--spin", str(inp["spin"]),
            "--xc", s["xc"], "--basis", s["basis"], "--solvent", s["solvent"],
            "--iso", repr(s["iso"]), "--spacing", repr(s["spacing_bohr"]),
            "--margin", repr(s["margin_bohr"]),
            "--geometry-origin", rec["key_fields"]["analysis"]["geometry_method"]]
    if s["density_fit"]:
        argv.append("--density-fit")
    if s["include_reaction_field"]:
        argv.append("--include-reaction-field")
    return argv


def _sha1(text):
    import hashlib
    return hashlib.sha1(text.encode()).hexdigest()


def code_identity():
    """Which espmap produced a number. Review 3 H14: `git rev-parse HEAD` was
    null outside git, reported a clean HEAD for an edited tree, and inside
    someone else's repository would report THEIR commit. Now the commit is
    taken only if this very file is tracked there, with --dirty, and
    "unknown" otherwise -- never a clean commit by default.

    Review 5 R18: espmap_version was the INSTALLED metadata, which names
    whatever espmap distribution the interpreter holds, not the file that
    ran: every P3 record said 0.1.0 (a stale editable install) from 0.1.2
    code. It is now __version__ of this module; the installed metadata is
    kept, under its own name, as installed_distribution_version."""
    import subprocess
    from importlib import metadata
    try:
        installed = metadata.version("espmap")
    except metadata.PackageNotFoundError:
        installed = None
    here = os.path.dirname(os.path.abspath(__file__))
    code = "unknown"
    try:
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch",
                                  os.path.basename(__file__)], cwd=here,
                                 capture_output=True, text=True, timeout=10)
        if tracked.returncode == 0:
            d = subprocess.run(["git", "describe", "--always", "--dirty", "--abbrev=12"],
                               cwd=here, capture_output=True, text=True, timeout=10)
            code = d.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        pass
    return dict(espmap_version=__version__, code=code,
                installed_distribution_version=installed)


def provenance():
    """Versions and threads that move results, for --json (review 3 H14)."""
    import platform
    from importlib import metadata

    def ver(dist):
        try:
            return metadata.version(dist)
        except metadata.PackageNotFoundError:
            return None
    from pyscf import lib
    blas = None
    try:
        from threadpoolctl import threadpool_info
        blas = [dict(api=p.get("internal_api"), version=p.get("version"),
                     threads=p.get("num_threads")) for p in threadpool_info()
                if p.get("user_api") == "blas"]
    except ImportError:
        pass
    return dict(**code_identity(), python=platform.python_version(),
                numpy=np.__version__, scipy=ver("scipy"),
                scikit_image=ver("scikit-image"), pyscf=ver("pyscf"),
                rdkit=ver("rdkit"), matplotlib=ver("matplotlib"),
                threads=int(lib.num_threads()), blas=blas)


def analysis_fields(xc, basis, solvent, iso, spacing, density_fit, geometry_method,
                    reaction_field=False):
    """What makes two records' NUMBERS comparable (spec section 5). Colour
    settings are deliberately absent: review 3's single compare_key would
    have called two identical calculations incomparable over a colormap.

    Normalised first, since PySCF is case-insensitive: def2-TZVP and
    def2-tzvp are one calculation and must get one key. Site and partition
    schema versions are NOT here: they cannot change a global V_min, and the
    site block carries its own schema.

    geometry_method (see geometry_method()) is how the coordinates were
    made. Review 4 N2: it was absent, yet on small anions it is the largest
    measured axis (MMFF vs DFT up to 6.5 kcal/mol, above the basis row)."""
    solv = None if solvent is None or solvent.lower() == "vacuum" else solvent.lower()
    return dict(xc=xc_key(xc), basis=basis_key(basis), solvent=solv,
                solvent_esp=SOLVENT_ESP[bool(reaction_field)] if solv else None,
                spacing_bohr=float(spacing), density_fit=bool(density_fit),
                geometry_method=geometry_method,
                surface_definition=surface_definition(iso))


# display name -> the PySCF spelling that defines it. "b3lyp" is pinned to
# b3lypg (VWN-RPA, PySCF's default and Gaussian's B3LYP), so a site that sets
# B3LYP_WITH_VWN5 records its "b3lyp" as "b3lyp5", not under our name.
_XC_NAMES = (("hf", "hf"), ("b3lyp", "b3lypg"), ("b3lyp5", "b3lyp5"),
             ("pbe0", "pbe0"), ("pbe", "pbe"), ("blyp", "blyp"), ("bp86", "bp86"),
             ("tpss", "tpss"), ("m062x", "m06-2x"), ("m06", "m06"),
             ("camb3lyp", "cam-b3lyp"), ("wb97x", "wb97x"))


def _xc_resolved(xc):
    """The functional PySCF will run for xc, in a hashable canonical form:
    the hybrid coefficients and the sorted (libxc id, weight) terms."""
    from pyscf.dft import libxc
    hyb, fns = libxc.parse_xc(xc)
    r = lambda v: round(float(v), 12)  # noqa: E731
    return (tuple(r(h) for h in hyb), tuple(sorted((int(i), r(w)) for i, w in fns)))


def xc_name(xc):
    """The display name of the functional xc resolves to, or None."""
    want = _xc_resolved(xc)
    for name, spelling in _XC_NAMES:
        try:
            if _xc_resolved(spelling) == want:
                return name
        except Exception:               # a name this PySCF build lacks
            continue
    return None


def xc_key(xc):
    """One key per functional. Review 5 R6: the key hashed xc.lower(), so
    b3lyp, b3lypg and HYB_GGA_XC_B3LYP (one functional, one energy) got three
    keys, and B3LYP_WITH_VWN5 changed what "b3lyp" meant under one key.
    A functional without a display name is keyed by its resolved libxc form."""
    name = xc_name(xc)
    if name:
        return name
    (a, b, w), fns = _xc_resolved(xc)
    return f"libxc:hyb({a:g},{b:g},{w:g})" + "".join(f"+{i}*{c:g}" for i, c in fns)


def _dft_label(xc, basis, solvent):
    """dft:<xc>/<basis_key>[/<solvent>], or None if PySCF cannot read xc. A
    display name is kept as written, so a label reads the same everywhere."""
    solv = None if solvent in (None, "vacuum") else solvent.lower()
    if xc.lower() not in {n for n, _ in _XC_NAMES}:
        try:
            xc = xc_name(xc) or xc
        except Exception:
            return None
    return f"dft:{xc.lower()}/{basis_key(basis)}" + (f"/{solv}" if solv else "")


_FORCE_FIELDS = {"mmff94s": "mmff94s", "mmff94": "mmff94", "uff": "uff"}
# every label geometry_method() emits besides dft:..., each with an optional
# hydrogen suffix for a file whose hydrogens were added here
_BASE_LABELS = (*_FORCE_FIELDS, "etkdg_mmff94s", "etkdg_uff", "idealised", "undeclared")
_H_SUFFIXES = ("+mmff94s_h", "+added_h")


def declared_geometry_method(text):
    """Normalise a declared geometry origin to the geometry_method vocabulary,
    or None when it cannot be read.

    Accepts every label geometry_method() emits ("dft:b3lyp/def2tzvp/water",
    "mmff94s", "etkdg_mmff94s", "undeclared+mmff94s_h", "idealised"), so a
    record re-run from its own coordinates keeps its key (reproduce_argv,
    review 5 R6), and the free text bench/fixture_geometries.py writes
    ("B3LYP/def2-TZVP, PCM water, Cartesian BFGS ..."; "MMFF94s conformer 4
    of 100 ...")."""
    import re
    t = (text or "").strip()
    if not t:
        return None
    low = t.lower()
    suffix = next((s for s in _H_SUFFIXES if low.endswith(s)), "")
    core = low[:len(low) - len(suffix)]
    if re.fullmatch(r"dft:[a-z0-9_\-]+/[a-z0-9()*+,_\-]+(/[a-z0-9_\-]+)?", core):
        xc, rest = core[4:].split("/", 1)
        basis, _, solv = rest.partition("/")
        lab = _dft_label(xc, basis, solv or None)
        return lab + suffix if lab else None
    if core in _BASE_LABELS:
        return core + suffix
    first = re.split(r"[\s,;(]", low, 1)[0]
    if first in _FORCE_FIELDS:
        return _FORCE_FIELDS[first]
    m = re.match(r"([a-z0-9\-]+)/([a-z0-9()*+_\-]+)(?:\s*,\s*pcm\s+([a-z0-9_\-]+))?", low)
    if m and not m.group(1).startswith(tuple(_FORCE_FIELDS)):
        return _dft_label(m.group(1), m.group(2), m.group(3))
    return None


def geometry_method(source=None, forcefield=None, atoms=False, molecule=None,
                    optimize=None, declared=None):
    """How the computed coordinates were made, as one normalised label.

    source/forcefield are from_smiles/from_file's info["source"] and
    info["forcefield"]; atoms and molecule mark --atoms and --molecule
    runs; optimize is (xc, basis, solvent) when --optimize replaced the
    coordinates, which then decides the label on its own; declared is a
    file's or --geometry-origin's declared origin (declared_geometry_method).

    Review 5 R1: every kept file geometry was "file_frozen", so a DFT-geometry
    fixture shared analysis_key with an MMFF file of the same molecule, on
    the axis worth up to 6.5 kcal/mol for anions. A file now says what made
    it, or is "undeclared". A file declared as DFT and an --optimize run at
    the same level share "dft:...".

    Review 5 R6: the label names the method, not the input channel. It was
    "idealised_builtin" for --molecule and "atoms_string:undeclared" for the
    same atoms through --atoms, one calculation under two keys; a record
    re-run through --atoms declares its method with --geometry-origin."""
    if optimize is not None:
        return _dft_label(*optimize)
    # text that cannot be read counts as undeclared (main warns, and checks
    # --geometry-origin before any compute)
    label = declared_geometry_method(declared) if declared else None
    if atoms:
        return label or "undeclared"
    if molecule is not None:
        return "idealised"
    if source is None:
        raise ValueError("geometry_method needs a source, --atoms or --molecule")
    if source.startswith("file 3D coordinates used as-is"):
        base = label or "undeclared"
        if "hydrogen(s) added" not in source or base.endswith(_H_SUFFIXES):
            return base
        return base + ("+mmff94s_h" if (forcefield or "").startswith("MMFF94s")
                       else "+added_h")
    if "ETKDG" in source:
        # _embed_and_minimise falls back to UFF where MMFF94s has no parameters
        return "etkdg_uff" if (forcefield or "").startswith("UFF") else "etkdg_mmff94s"
    raise ValueError(f"unrecognised geometry source: {source!r}")


def validity(gate, verified, surface, *, electronic):
    """One place a consumer reads validity from (spec section 6), never from
    stderr, chi or a filename. sites is False until a site block exists.

    input is true only if a gate RAN and passed (not forced) and the
    computed structure was not found to differ. Review 4 N3: absence of a
    gate used to read as a pass, so ungated --atoms records said true.
    electronic is false when results.electronic has problems (review 5 R3)."""
    return dict(input=gate is not None and not gate["forced"] and verified is not False,
                surface=not surface["problems"], electronic=not electronic["problems"],
                sites=False)


def render_fields(cmap, kind, vlim, elev, azim, size, ambient):
    """What makes two PICTURES share one visual scale. vlim is the resolved
    range, so two auto-scaled panels match only if their ranges do."""
    return dict(cmap=cmap, kind=kind, vmin=round(float(vlim[0]), 6),
                vmax=round(float(vlim[1]), 6), elev=float(elev), azim=float(azim),
                size=int(size), ambient=float(ambient))


# --------------------------------------------------------------------- sites
# per-site-spec P2 (docs/per-site-spec.md). Atom maps and caller-defined
# reaction sites arrive with P1; here sites are the default groups (each heavy
# atom with its hydrogens) and atoms carry whatever map the input had.
SITES_SCHEMA = "espmap.sites/1"
# The whole --json record (review 4 N5). Its key tree is docs/record-schema.json,
# checked in both directions by the slow suite; outside analysis_key, like
# SITES_SCHEMA, so a schema bump never makes two calculations incomparable.
RECORD_SCHEMA = "espmap.record/1"
PARTITION_RULES = ("nearest_nucleus", "radius_scaled")
# Bondi van der Waals radii (A), with Mantina et al. 2009 for B, which Bondi
# does not list. Every TESTED_ELEMENTS entry has one; an element outside gets
# no radius-scaled partition rather than a guessed radius.
BONDI = {"H": 1.20, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
         "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
         "Se": 1.90}
_POLAR = {"N", "O", "S"}


def triangle_areas(P, faces):
    tri = P[faces]
    return 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)


def _area_quantile(v, w, q):
    o = np.argsort(v, kind="stable")
    c = np.cumsum(w[o]) / w.sum()
    return float(v[o][min(int(np.searchsorted(c, q)), len(v) - 1)])


def area_stats(P, faces, V):
    """Area-weighted statistics of V over a set of triangles (P in A).

    Weighted by triangle area, with V averaged over each triangle's three
    vertices, as bench/conformers.py measured. Never over vertices:
    marching-cubes vertex density follows grid orientation and curvature,
    not chemistry. Extrema are over the triangles' vertices, where V is
    computed exactly, and are NOT robust on a small set: lead with the mean
    and percentiles. One helper for the whole surface and for every site.
    V_area_sd is about the overall mean; GIPF's variances, each about its
    own sign's mean, are gipf_stats.
    """
    keys = ("V_area_mean_kcal_mol", "V_area_sd_kcal_mol", "V_p05_kcal_mol",
            "V_p50_kcal_mol", "V_p95_kcal_mol", "V_min_kcal_mol", "V_min_at_A",
            "V_max_kcal_mol", "V_max_at_A")
    a = triangle_areas(P, faces) if len(faces) else np.zeros(0)
    if a.sum() == 0:
        return dict(area_A2=0.0, area_pos_A2=0.0, area_neg_A2=0.0, **dict.fromkeys(keys))
    vf = V[faces].mean(axis=1)
    mean = float((a * vf).sum() / a.sum())
    idx = np.unique(faces)
    i_min, i_max = idx[np.argmin(V[idx])], idx[np.argmax(V[idx])]
    return dict(area_A2=float(a.sum()), area_pos_A2=float(a[vf > 0].sum()),
                area_neg_A2=float(a[vf < 0].sum()),
                V_area_mean_kcal_mol=mean,
                V_area_sd_kcal_mol=float(np.sqrt((a * (vf - mean) ** 2).sum() / a.sum())),
                V_p05_kcal_mol=_area_quantile(vf, a, 0.05),
                V_p50_kcal_mol=_area_quantile(vf, a, 0.50),
                V_p95_kcal_mol=_area_quantile(vf, a, 0.95),
                V_min_kcal_mol=float(V[i_min]), V_min_at_A=[float(x) for x in P[i_min]],
                V_max_kcal_mol=float(V[i_max]), V_max_at_A=[float(x) for x in P[i_max]])


# A fraction of the surface AREA (review 5 R14 suggested 5 %); not the same
# quantity as SIGN_CHANGE_MIN, which is a ratio of potential values.
GIPF_MIN_AREA_FRACTION = 0.05


def gipf_stats(P, faces, V):
    """GIPF's global surface statistics (Murray et al. 1994, murray1994gipf):
    the mean of each sign, sigma_+^2 and sigma_-^2 each about its own sign's
    mean, sigma^2_tot = sigma_+^2 + sigma_-^2, and the balance
    nu = sigma_+^2 sigma_-^2 / (sigma^2_tot)^2, at most 0.25.

    Review 4 N11: area_stats' V_area_sd is one SD about the overall mean and
    was cited as GIPF's variance; it is not, and stays as what it is. GIPF
    averaged over surface points; here, as everywhere in espmap, over
    triangles weighted by area (V averaged per triangle), so compare with
    Multiwfn only when its weighting is the same. A sign covering less than
    GIPF_MIN_AREA_FRACTION of the surface has a null mean and variance
    (undefined, not zero), and nu is then null, with the reason: review 5
    R14 found one triangle at +0.01 kcal/mol turning nu from null into 0.0,
    the knife-edge review 3 H11 removed from the scale. surface_stats holds
    the signed areas that decide it.
    """
    a = triangle_areas(P, faces)
    vf = V[faces].mean(axis=1)
    total = a.sum()
    small = []

    def part(m, name):
        if total == 0 or a[m].sum() / total < GIPF_MIN_AREA_FRACTION:
            small.append(name)
            return None, None
        mu = float((a[m] * vf[m]).sum() / a[m].sum())
        return mu, float((a[m] * (vf[m] - mu) ** 2).sum() / a[m].sum())
    mp, sp = part(vf > 0, "positive")
    mn, sn = part(vf < 0, "negative")
    tot = None if sp is None and sn is None else (sp or 0.0) + (sn or 0.0)
    nu = (sp * sn / tot ** 2 if sp is not None and sn is not None and tot > 0 else None)
    reason = (None if not small else
              f"{' and '.join(small)} area below {GIPF_MIN_AREA_FRACTION:.0%} of the surface: "
              f"that sign's statistics and nu are undefined")
    return dict(V_pos_mean_kcal_mol=mp, V_neg_mean_kcal_mol=mn,
                sigma2_pos_kcal2_mol2=sp, sigma2_neg_kcal2_mol2=sn,
                sigma2_tot_kcal2_mol2=tot, nu=nu, reason=reason,
                min_area_fraction=GIPF_MIN_AREA_FRACTION,
                weighting="triangle area, V averaged per triangle")


def record_bonds(probe):
    """The approved structure's bonds, hydrogens included, in record atom
    order: [[i, j, order], ...] with i < j and one Kekule structure (orders
    1-3). Review 5 A1: a consumer could rebuild bonds only from a mapped
    SMILES. Valid only when verify_mapping confirmed the order (main writes
    null otherwise). The SCF never sees bond orders (invariant 3), so they
    are the drawing's, not a computed property."""
    from rdkit import Chem
    m = Chem.AddHs(probe)
    Chem.Kekulize(m, clearAromaticFlags=True)
    return sorted([min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   int(b.GetBondTypeAsDouble())] for b in m.GetBonds())


def partition_surface(P, faces, coords, symbols, rule):
    """The atom owning each triangle, by its centroid: nearest nucleus, or
    nearest by d / r_Bondi. Triangles, not vertices, so atom areas sum to the
    surface area exactly. None if rule needs a radius an element lacks."""
    X = np.asarray(coords, float)
    d = np.linalg.norm(P[faces].mean(axis=1)[:, None, :] - X[None], axis=2)
    if rule == "radius_scaled":
        if any(s not in BONDI for s in symbols):
            return None
        d = d / np.array([BONDI[s] for s in symbols])[None, :]
    return np.argmin(d, axis=1)


def site_parents(coords, symbols):
    """Each hydrogen's heavy atom (the nearest one); None for heavy atoms,
    and None for every atom of a molecule with no heavy atom (H2): review 4
    N7, argmin over an empty list crashed the run after the SCF."""
    X = np.asarray(coords, float)
    heavy = [i for i, s in enumerate(symbols) if s != "H"]
    if not heavy:
        return [None] * len(symbols)
    return [None if s != "H" else
            heavy[int(np.argmin(np.linalg.norm(X[heavy] - X[i], axis=1)))]
            for i, s in enumerate(symbols)]


_DELTA = ("area_A2", "V_area_mean_kcal_mol", "V_p05_kcal_mol", "V_p95_kcal_mol",
          "V_min_kcal_mol", "V_max_kcal_mol")


def _delta(nn, rs):
    if rs is None:
        return None
    return {k: None if nn[k] is None or rs[k] is None else rs[k] - nn[k] for k in _DELTA}


def load_site_groups(path):
    """--sites file: {"by": "map" | "index", "groups": [{"label", "atoms"}, ...]}."""
    import json
    with open(path) as f:
        return json.load(f)


def resolve_site_groups(spec, symbols, coords, atom_maps):
    """Caller-named reaction sites -> (by, [(label, sorted atom indices)]).

    The file says how it counts atoms (review 4 N6): "by": "map" names atoms
    by atom-map number and is refused on an input with no maps; "by": "index"
    names them by RDKit atom index in FILE order, hydrogens included, and
    refuses an index that is a hydrogen (name its heavy atom with "N+H").
    Before, a name meant a map on mapped input and an index otherwise,
    silently: maps are 1-based and indices 0-based, so a lost map always
    landed on a plausible neighbour. "3+H" adds that atom's hydrogens (a C-H
    site is the unit a desaturation or oxidation question is about). espmap
    only checks that every named atom exists: it does not know what an
    enolate is. Raises ValueError naming the fault. Groups may overlap.
    """
    if not isinstance(spec, dict) or spec.get("by") not in ("map", "index"):
        raise ValueError('a --sites file says how it counts atoms: {"by": "map" | "index", '
                         '"groups": [{"label": ..., "atoms": [...]}, ...]}')
    by, groups = spec["by"], spec.get("groups")
    if not isinstance(groups, list):
        raise ValueError('"groups" must be a list of {label, atoms}')
    if by == "map":
        if not atom_maps:
            raise ValueError('"by": "map", but the input has no atom maps (were they lost '
                             'on the way?); write the input with maps, or use "by": "index"')
        by_name = {m: i for i, m in atom_maps.items()}
    else:
        by_name = {i: i for i, s in enumerate(symbols)}
    naming = "atom map" if by == "map" else "atom index"
    parent = site_parents(coords, symbols)
    out, seen = [], set()
    for g in groups:
        label = g.get("label") if isinstance(g, dict) else None
        if not isinstance(label, str) or not label:
            raise ValueError(f"group {g!r} has no label")
        if label in seen:
            raise ValueError(f"label {label!r} is used twice")
        seen.add(label)
        members = set()
        for tok in g.get("atoms") or []:
            name, with_h = tok, False
            if isinstance(tok, str):
                head, plus, tail = tok.partition("+")
                if not head.isdigit() or (plus and tail != "H"):
                    raise ValueError(f"{label!r}: {tok!r} is not an atom ('3' or '3+H')")
                name, with_h = int(head), bool(plus)
            if not isinstance(name, int) or name not in by_name:
                raise ValueError(f"{label!r}: no {naming} {tok} in this molecule")
            i = by_name[name]
            if symbols[i] == "H":
                raise ValueError(f"{label!r}: {naming} {tok} is a hydrogen; a site names "
                                 f"heavy atoms, and 'N+H' adds a heavy atom's hydrogens")
            members.add(i)
            if with_h:
                members.update(j for j, p in enumerate(parent) if p == i)
        if not members:
            raise ValueError(f"{label!r} names no atoms")
        out.append((label, sorted(members)))
    return by, out


def site_report(P, faces, V, coords, symbols, atom_maps=None, h_source=None, groups=None,
                source="caller"):
    """atoms and reaction_sites under both partition rules (spec sections
    1-3). P and coords in A. Reaction sites are the caller's groups
    [(label, atom indices)], recorded with `source` ("caller:map" or
    "caller:index", how the --sites file named them), or by default each
    heavy atom with its hydrogens ("default"). Every atom is listed, a buried
    one with area 0 and null potentials, so 'not on the surface' never reads
    as 'not computed'."""
    parent = site_parents(coords, symbols)
    owner = {r: partition_surface(P, faces, coords, symbols, r) for r in PARTITION_RULES}

    def per_rule(members):
        out = {}
        for r in PARTITION_RULES:
            out[r] = (None if owner[r] is None else
                      area_stats(P, faces[np.isin(owner[r], members)], V))
        out["delta"] = _delta(out["nearest_nucleus"], out["radius_scaled"])
        return out

    atoms = []
    for i, s in enumerate(symbols):
        atoms.append(dict(index=i, atom_map=(atom_maps or {}).get(i), element=s,
                          parent=parent[i],
                          h_source=(h_source[i] if h_source and s == "H" else None),
                          **per_rule([i])))
    if groups is None:
        # review 5 R5: labelled by atom map on a mapped input, as caller sites
        # are, so nitrite's O sites read O2/O3 (its maps), not O1/O2 (indices)
        # (a partial map, BLOCKed by the gate unless forced, falls back to indices)
        by_map = bool(atom_maps) and all(h in atom_maps for h, s in enumerate(symbols)
                                         if s != "H")
        source = "default:map" if by_map else "default:index"
        groups = [(f"{symbols[h]}{atom_maps[h] if by_map else h}",
                   [h] + [i for i, p in enumerate(parent) if p == h])
                  for h, s in enumerate(symbols) if s != "H"]
    sites = []
    for label, members in groups:
        # a generated hydrogen in the group, on a polar atom of the group
        polar_gen = h_source is not None and any(
            symbols[i] == "H" and h_source[i] not in (None, "supplied")
            and parent[i] in members and symbols[parent[i]] in _POLAR for i in members)
        sites.append(dict(label=label, source=source, atoms=list(members),
                          h_generated_polar=bool(polar_gen), **per_rule(list(members))))
    # Caller groups may overlap or leave atoms out; then site areas do not
    # sum to the surface (the C18 fixtures summed to 1.041). Say so.
    for s in sites:
        s["overlaps"] = [t["label"] for t in sites
                         if t is not s and set(t["atoms"]) & set(s["atoms"])]
    counts = np.bincount([i for s in sites for i in s["atoms"]], minlength=len(symbols))
    compared = owner["radius_scaled"] is not None
    return dict(schema=SITES_SCHEMA,
                partition=dict(primary="nearest_nucleus",
                               compared="radius_scaled" if compared else None,
                               unit="triangle",
                               radii="Bondi; Mantina 2009 for B" if compared else None),
                atoms=atoms, reaction_sites=sites,
                sites_tile_surface=bool((counts == 1).all()))


def site_checks(rep, P, faces, V):
    """What validity.sites now rests on (review 4 N5): the atoms' areas add
    up to the surface under each rule, every atom is listed once, and a
    zero-area entry has null potentials while a non-zero one has values.
    Returns the problems found; empty means the block passed."""
    out = []
    total = float(triangle_areas(P, faces).sum())
    if [a["index"] for a in rep["atoms"]] != list(range(len(rep["atoms"]))):
        out.append("atoms are not listed once each, in index order")
    for r in PARTITION_RULES:
        rows = [a[r] for a in rep["atoms"] if a[r] is not None]
        if not rows:
            continue
        s = sum(x["area_A2"] for x in rows)
        if abs(s - total) > 1e-9 * max(total, 1.0):
            out.append(f"{r}: atom areas {s:.6f} do not conserve the surface {total:.6f} A^2")
    for e in rep["atoms"] + rep["reaction_sites"]:
        for r in PARTITION_RULES:
            x = e[r]
            if x is not None and (x["area_A2"] == 0) != (x["V_area_mean_kcal_mol"] is None):
                name = e.get("label", e.get("index"))
                out.append(f"{name} {r}: potentials must be null exactly when the area is 0")
    return out


def tri_atom_arrays(has_sites, P, faces, coords, symbols):
    """The NPZ's per-triangle partition, {"tri_atom_<rule>": owner index}.

    Written exactly when the record has a site table: review 4 N8 withheld it
    from invalid input or surface, and review 5 R9 found it still written for
    H2 (no heavy atom) and for a block that failed its own checks. A rule
    that cannot run (an element without a Bondi radius) is omitted rather
    than filled with -1: absence has one encoding, a missing array."""
    if not has_sites:
        return {}
    out = {}
    for r in PARTITION_RULES:
        t = partition_surface(P, faces, coords, symbols, r)
        if t is not None:
            out[f"tri_atom_{r}"] = t
    return out


def sites_block(valid, P, faces, V, coords, symbols, input_reason=None, **kw):
    """(site report, None), or (None, reason) when the record is not valid:
    per-site values near a clipped box face, or for a forced-through
    structure, would be worse than the global ones (spec section 3). A block
    that fails its own checks is withheld the same way. input_reason is
    input_invalid_reason(...), which names the failing check (review 5 R8)."""
    if not valid["input"]:
        return None, input_reason or "input not valid"
    if not valid["surface"]:
        return None, "surface not valid (see results.surface.problems)"
    if not valid.get("electronic", True):
        return None, "electronic structure not valid (see results.electronic.problems)"
    if all(s == "H" for s in symbols):
        return None, ("no heavy atom: sites are heavy atoms with their hydrogens, "
                      "and this molecule has none")
    rep = site_report(P, faces, V, coords, symbols, **kw)
    problems = site_checks(rep, P, faces, V)
    if problems:
        return None, "site block failed its own checks: " + "; ".join(problems)
    return rep, None


H_SOURCES = ("supplied", "generated", "generated_relaxed", "optimised")


def site_context(*, solvent_esp, heavy_atoms, h_source, relaxation, n_rot, surface,
                 identity, budget):
    """The run-level fields main() adds to a site block (spec sections 4-6).
    Every hydrogen source is counted, zeros included: an absent key read as
    zero only if the consumer guessed so (review 4 N5)."""
    n_h = [x for x in h_source if x is not None]
    return dict(
        solvent_esp=solvent_esp,
        geometry_provenance=dict(
            heavy_atoms=heavy_atoms,
            hydrogens={k: n_h.count(k) for k in H_SOURCES},
            relaxation=relaxation,
            caveat="hydrogens not supplied were placed for the isolated molecule, "
                   "without any protein or other environment"),
        conformer_scope=(f"single conformer; {n_rot} rotatable bond(s)"
                         if n_rot is not None else "single geometry, as given"),
        topology=dict(genus=surface["genus"], components=surface["components"]),
        identity=identity, budget=budget)


# Measured sensitivities a record points to (spec section 5). Each value is
# recomputed from its bench/ file by tests/test_claims.py, under the same id
# the README table uses. A row is evidence about the method on the molecule
# named, not an error bar on this run's sites.
BUDGET = (
    dict(id="basis", axis="basis def2-SVP -> def2-TZVPP", statistic="shift, sigma-hole V_max",
         kcal_mol=3.54, measured_on="chlorobenzene", source="bench/error_budget.jsonl"),
    dict(id="solvent water", axis="solvent gas -> PCM water, polarised solute only",
         statistic="shift, sigma-hole V_max", kcal_mol=-3.60, measured_on="chlorobenzene",
         source="bench/solvent_convention.jsonl"),
    dict(id="solvent chloroform", axis="solvent gas -> PCM chloroform, polarised solute only",
         statistic="shift, sigma-hole V_max", kcal_mol=-2.69, measured_on="chlorobenzene",
         source="bench/solvent_convention.jsonl"),
    dict(id="functional", axis="functional B3LYP / PBE0 / M06-2X / HF",
         statistic="spread, sigma-hole V_max", kcal_mol=0.69, measured_on="chlorobenzene",
         source="bench/error_budget.jsonl"),
    dict(id="grid", axis="grid spacing 0.26 -> 0.16 bohr", statistic="shift, sigma-hole V_max",
         kcal_mol=0.04, measured_on="chlorobenzene", source="bench/error_budget.jsonl"),
    dict(id="anion geometry", axis="geometry MMFF -> DFT, PCM water",
         statistic="max |shift|, V_min and V_max", kcal_mol=6.5,
         measured_on="nitrite, thiocyanate, enolate", source="bench/anion_geometry.jsonl"),
    dict(id="partition small anions", axis="partition rule, radius_scaled - nearest_nucleus",
         statistic="max |delta|, site area-weighted mean", kcal_mol=3.0,
         measured_on="nitrite, thiocyanate, enolate", source="bench/site_partition.jsonl"),
    dict(id="partition C18", axis="partition rule, radius_scaled - nearest_nucleus",
         statistic="max |delta|, site area-weighted mean", kcal_mol=3.95,
         measured_on="Z3Z13-18:OAc, 10 conformers", source="bench/site_partition.jsonl"),
    dict(id="conformer C18", axis="conformer, extended and folded",
         statistic="median site range, site area-weighted mean", kcal_mol=5.1,
         measured_on="Z3Z13-18:OAc, 10 conformers", source="bench/site_partition.jsonl"),
    dict(id="grid mean", axis="grid spacing 0.26 -> 0.16 bohr",
         statistic="max |shift|, site area-weighted mean", kcal_mol=0.15,
         measured_on="nitrite, thiocyanate, enolate, chlorobenzene",
         source="bench/site_partition.jsonl"),
)


def budget_rows(*, basis, solvent, charge, n_rot, geometry_method):
    """The BUDGET rows that apply to a run with these settings: the basis
    row for a double-zeta run (it measures def2-SVP -> def2-TZVPP), solvent
    rows for a gas-phase run, the anion-geometry row for an anion whose
    geometry is not a DFT one, and the conformer row unless the molecule is
    known to be rigid.

    Review 5 R1: the anion row ("MMFF -> DFT geometry, up to 6.5") was
    attached on charge alone, so the DFT-geometry P3 fixtures carried an
    uncertainty they do not have; the basis row was attached to def2-TZVP
    runs, which it does not describe."""
    keep = {"basis": basis_key(basis) in SMALL_BASIS,
            "solvent water": solvent is None, "solvent chloroform": solvent is None,
            "anion geometry": charge < 0 and not geometry_method.startswith("dft:"),
            "conformer C18": n_rot is None or n_rot > 0}
    return [dict(r) for r in BUDGET if keep.get(r["id"], True)]


def sigma_hole(verts, V, mol, i_sub=6, pr_min=0.3, pe_max=1.6, i_bond=None):
    """V_max on the outer cap of the substituent, along the X-to-neighbour axis.

    The sigma-hole and its surface V_s,max are established (Clark et al.
    2007; Murray, Lane & Politzer 2009; docs/prior-art.md). What is espmap's
    own is this implementation: the cone around the bond axis, below.

    v3 took the axis as (X - mean of atoms 0..5), hard-coding "the first six
    atoms are the benzene ring". For any other molecule that reference point
    is meaningless: measured axis errors of 25-40 degrees on simple alkyl
    halides, and on CC(=O)OCCCl it reported -0.05 kcal/mol where the correct
    axis gives +1.84 -- a sign flip on the reported quantity.

    The axis is now the bond vector from the substituent's single heavy
    neighbour to the substituent, which is what a sigma-hole is defined along.
    i_bond may be supplied; otherwise the nearest atom is used and that
    choice is reported.

    The cone thresholds were a review concern that tested NEGATIVE: V_max was
    identical (+9.84) across pr_min in {0, 0.3, 0.6} x pe_max in {1.0, 1.6,
    2.2}, because the maximum is interior to even the tightest cone. Returned
    with its vertex count so an empty or boundary-hugging selection is visible.
    """
    co = mol.atom_coords()
    X = co[i_sub]
    if i_bond is None:
        i_bond = sigma_axis_partner(co, i_sub)
    axis = X - co[i_bond]
    axis /= np.linalg.norm(axis)
    rel = verts - X
    pr = rel @ axis
    pe = np.linalg.norm(rel - np.outer(pr, axis), axis=1)
    m = (pr > pr_min) & (pe < pe_max)
    if not m.any():
        return float("nan"), 0
    return float(V[m].max()), int(m.sum())


# ------------------------------------------------------------------- selftest
def selftest(basis="def2-svp"):
    """F7: the v1 sign error was caught by hand. Now it is an assertion.

    Water, B3LYP. The windows are deliberately wide enough to pass at any sane
    basis and narrow enough that a kernel-sign flip or a units error fails by
    orders of magnitude.
    """
    water = "O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692"
    mol, mf = run_scf(water, basis=basis)
    dm = total_dm(mf)
    B = 1.0 / BOHR
    probes = {
        "lone pair, +z 1.5 A from O": (np.array([0.0, 0.0, 1.6173]), (-70.0, -30.0)),
        "1.0 A beyond an O-H bond": (np.array([0.0, 1.546, -1.055]), (35.0, 80.0)),
        "far field, 25 A": (np.array([0.0, 0.0, 25.0]), (-2.0, 2.0)),
    }
    pts = np.array([p for p, _ in probes.values()]) * B
    V = esp_at_points(mol, dm, pts) * HA2KCAL
    ok = True
    print(f"--- selftest, B3LYP/{basis}, water ---")
    for (label, (_, (lo, hi))), v in zip(probes.items(), V):
        good = lo <= v <= hi
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {label:30s} {v:+8.2f} kcal/mol   "
              f"expect [{lo:+.0f}, {hi:+.0f}]")
    good = bool(V[0] < 0 < V[1])
    ok &= good
    print(f"  {'PASS' if good else 'FAIL'}  lone pair negative AND donor H positive")

    # ALIE (P4): exact identities, so a wrong sign on e or a dropped spin
    # channel fails by far more than rounding
    near = pts[:2]                      # the 25 A probe is outside every orbital
    I = alie_at_points(mol, mf, near) * HA2EV
    good = bool(np.all(np.isfinite(I)) and np.all(I > 0))
    ok &= good
    print(f"  {'PASS' if good else 'FAIL'}  ALIE positive at the probes "
          f"({', '.join(f'{x:.2f}' for x in I)} eV)")
    for label, atom, spin in (("H2, RKS, one doubly occupied orbital",
                               "H 0 0 0; H 0 0 0.74", 0),
                              ("H atom, UKS, one alpha orbital", "H 0 0 0", 1)):
        m1, f1 = run_scf(atom, basis=basis, spin=spin)
        e = np.atleast_2d(np.asarray(f1.mo_energy))
        occ = np.atleast_2d(np.asarray(f1.mo_occ)) > 0
        want = float(-e[occ].max())
        got = alie_at_points(m1, f1, near)
        good = bool(np.allclose(got, want, rtol=0, atol=1e-10))
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  ALIE = -e_1 everywhere: {label} "
              f"(max dev {np.abs(got - want).max():.1e} Ha)")

    # colormap orientation. Caught two inverted sequential maps when added;
    # the convention is intense/warm = electron-rich on every scale.
    def _lum(c):
        return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]
    cmaps = build_maps()
    for name, cmap in sorted(cmaps.items()):
        lo, hi = np.array(cmap(0.0)[:3]), np.array(cmap(1.0)[:3])
        if name in SEQUENTIAL:          # t=0 is V_min, which must be darkest
            good = _lum(lo) < _lum(hi)
            why = "V_min darkest"
        else:                           # t=0 is the negative end, which must be warm
            good = lo[0] > lo[2] and hi[2] > hi[0]
            why = "negative end warm, positive end cool"
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  cmap {name:9s} {why}")

    print("selftest:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ----------------------------------------------------------------------- main
def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true",
                    help="validate the ESP sign and magnitude, then exit")
    ap.add_argument("--smiles", help="build 3D coordinates from SMILES (RDKit)")
    ap.add_argument("--file", help="read a structure file: .mol .sdf .cdxml .xyz "
                                   ".pdb — i.e. what a drawing program exports")
    ap.add_argument("--record", type=int, default=None,
                    help="1-based record to compute from a multi-record SDF/CDXML; "
                         "without it such a file is refused")
    ap.add_argument("--geometry-origin", default=None, metavar="LABEL",
                    help="how a --file or --atoms geometry was made, e.g. "
                         "'dft:b3lyp/def2-tzvp/water', 'mmff94s', 'etkdg_mmff94s' "
                         "or 'idealised' (any geometry_method label). It enters "
                         "analysis_key; an SDF declares it instead as an "
                         "'espmap_geometry_origin' (or 'geometry') property. "
                         "Undeclared geometries are keyed 'undeclared'")
    ap.add_argument("--seed", type=lambda s: int(s, 0), default=0xC0FFEE,
                    help="ETKDG random seed; fixed so runs are reproducible")
    ap.add_argument("--nconfs", type=int, default=20,
                    help="conformers to embed; the lowest-energy one is used")
    ap.add_argument("--optimize", action="store_true",
                    help="DFT geometry optimisation after the force field "
                         "(slow; measured worth 0.40 kcal/mol on a sigma-hole, up to 6.5 on small anions)")
    ap.add_argument("--molecule", default="bromobenzene", choices=sorted(MOLECULES))
    ap.add_argument("--atoms", help="PySCF atom string; overrides --molecule")
    ap.add_argument("--charge", type=int, default=None,
                    help="net charge; read from SMILES when --smiles is used")
    ap.add_argument("--spin", type=int, default=0, help="2S, i.e. n_alpha - n_beta")
    ap.add_argument("--scale", choices=["auto", "diverging", "sequential"],
                    default="auto",
                    help="auto picks diverging only if the surface changes sign")
    ap.add_argument("--yes", action="store_true",
                    help="proceed past the basis-size cost guard")
    ap.add_argument("--check", action="store_true",
                    help="run the structure admissibility check and exit "
                         "without computing anything")
    ap.add_argument("--force", action="store_true",
                    help="compute even if the structure check rejects")
    ap.add_argument("--require-stereo", action="store_true",
                    help="refuse a structure whose input leaves any stereo element "
                         "open (a centre, E/Z, an allene, an atropisomer or helix); "
                         "3D coordinates fix every element")
    ap.add_argument("--basis", default="def2-tzvp")
    ap.add_argument("--density-fit", action="store_true",
                    help="RI-J/K: faster SCF up to ~600 nao; moves V by <= 0.02 kcal/mol (measured)")
    ap.add_argument("--gpu", action="store_true",
                    help="SCF on GPU4PySCF (NVIDIA only); validated with --density-fit")
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--solvent", default="vacuum",
                    choices=sorted(SOLVENTS),
                    help="PCM continuum. The map is the potential of the "
                         "PCM-POLARISED SOLUTE only; the PCM surface charges' "
                         "own potential (the reaction field) is left out unless "
                         "--include-reaction-field. Measured under that "
                         "convention: gas->water moves the chlorobenzene "
                         "sigma-hole by -3.6 kcal/mol. Default vacuum, stated "
                         "in the output either way")
    ap.add_argument("--include-reaction-field", action="store_true",
                    help="with --solvent, add the potential of the PCM surface "
                         "charges to the map. A different quantity: it changes "
                         "analysis_key and the sign of the solvent effect on the "
                         "chlorobenzene sigma-hole (review 3 H7)")
    ap.add_argument("--iso", type=float, default=ISO_BADER)
    ap.add_argument("--spacing", type=float, default=0.26,
                    help="bohr; V_max converged to 0.04 kcal/mol at this value")
    ap.add_argument("--margin", type=float, default=5.0,
                    help="bohr of vacuum around the nuclei; below 4 risks clipping")
    ap.add_argument("--no-strict", action="store_true",
                    help="render even if the isosurface fails its integrity check")
    ap.add_argument("--cmap", default=None,
                    help="default: vik for a diverging scale, lajolla for a "
                         "sequential one (ions)")
    ap.add_argument("--force-cmap", action="store_true",
                    help="allow a diverging map on a sequential scale or vice versa")
    ap.add_argument("--vmin", type=float, default=None,
                    help="explicit lower bound for a SEQUENTIAL scale (ions)")
    ap.add_argument("--vmax", type=float, default=None,
                    help="explicit upper bound for a SEQUENTIAL scale (ions)")
    ap.add_argument("--lim", type=float, default=None,
                    help="range in kcal/mol: symmetric +/-LIM on a diverging "
                         "scale. Omitting it auto-scales to this molecule, "
                         "which makes molecules incomparable")
    ap.add_argument("--elev", type=float, default=20)
    ap.add_argument("--azim", type=float, default=172)
    ap.add_argument("--size", type=int, default=700)
    ap.add_argument("--ambient", type=float, default=0.62)
    ap.add_argument("--out", default="esp.png")
    ap.add_argument("--html", help="also write an interactive 3D view (needs plotly)")
    ap.add_argument("--npz", help="also save the surface (verts, faces, V, atoms) "
                                  "so it can be re-rendered without an SCF")
    ap.add_argument("--json", help="also write the settings and results")
    ap.add_argument("--sites", help='JSON {"by": "map"|"index", "groups": [{label, atoms}]}: '
                                    "reaction sites to report instead of one per heavy "
                                    "atom; atoms by atom-map number or by atom index in "
                                    "file order (heavy atoms only), '3+H' adds 3's hydrogens")
    a = ap.parse_args()

    if a.selftest:
        raise SystemExit(selftest())
    if a.html:
        try:                        # fail before the SCF, not after it
            import plotly  # noqa: F401
        except ImportError:
            raise SystemExit("--html needs plotly: pip install plotly") from None
    if a.optimize:
        try:                        # same reason: the optimiser runs AFTER the first SCF
            import berny  # noqa: F401  (pyberny; pyscf.geomopt.berny_solver needs it)
        except ImportError:
            raise SystemExit("--optimize needs pyberny: pip install 'espmap[opt]' "
                             "(or pip install pyberny)") from None
    if a.geometry_origin is not None:
        if not (a.file or a.atoms) or a.smiles:
            raise SystemExit("--geometry-origin describes supplied coordinates; it "
                             "needs --file or --atoms (SMILES geometries are "
                             "generated here and labelled accordingly)")
        if declared_geometry_method(a.geometry_origin) is None:
            raise SystemExit(f"--geometry-origin {a.geometry_origin!r} is not a "
                             f"recognised origin; use a geometry_method label, e.g. "
                             f"'dft:b3lyp/def2-tzvp/water', 'mmff94s', 'etkdg_mmff94s' "
                             f"or 'idealised'")
    if a.include_reaction_field and a.solvent == "vacuum":
        raise SystemExit("--include-reaction-field needs --solvent: in the gas "
                         "phase there is no reaction field to add")
    if a.lim == 0.0:
        raise SystemExit("--lim 0 is not a range")
    err = cmap_preflight(a.scale, a.cmap, a.force_cmap, available=set(build_maps()))
    if err:
        raise SystemExit(err)

    charge, spin, i_sub = a.charge or 0, a.spin, 6
    sh_reason = None
    pre_opt_sha1 = None
    gi = None
    if a.smiles and a.file:
        raise SystemExit("--smiles and --file are mutually exclusive")

    # --atoms and --molecule coordinates are gated too (review 4 N3): their
    # bonds are perceived from the coordinates, as for an XYZ file.
    atoms_spec = (None if (a.smiles or a.file) else
                  a.atoms or monosubstituted_benzene(*MOLECULES[a.molecule]))
    if a.smiles or a.file or atoms_spec:
        probe, probs = parse_for_check(smiles=a.smiles, path=a.file,
                                       charge=a.charge or 0, record=a.record,
                                       atoms=atoms_spec)
        if probe is None:
            print_check(BLOCK, probs, a.smiles or str(a.file or "--atoms"))
            raise SystemExit(1)
        from rdkit import Chem
        verdict, findings = check_structure(Chem.RemoveHs(probe), a.basis,
                                            spin=a.spin, maps_from=probe,
                                            smiles=a.smiles, charge=a.charge,
                                            require_stereo=a.require_stereo)
        print_check(verdict, findings, Chem.MolToSmiles(Chem.RemoveHs(probe)))
        if verdict == BLOCK and not a.force:
            raise SystemExit("structure rejected; fix the drawing, or pass "
                             "--force to compute anyway")
        if a.check:
            raise SystemExit(0 if verdict != BLOCK else 1)

    from rdkit import Chem as _C
    approved = _C.MolToSmiles(_C.RemoveHs(probe))
    gate = gate_report(verdict, findings, a.force)

    if a.file:
        try:
            spec, gi = from_file(a.file, n_confs=a.nconfs, seed=a.seed,
                                 charge=a.charge or 0, record=a.record)
        except ValueError as e:              # e.g. H_COUNT_MISMATCH: refuse, no traceback
            raise SystemExit(f"structure rejected: {e}") from None
        print(f"{a.file}  ->  {gi['formula']}, charge {gi['charge']:+d}"
              + (f"  (record {gi['record']['index']} of {gi['record']['of']})"
                 if gi["record"]["of"] > 1 else ""))
        print(f"  SMILES as read: {gi['smiles']}")
        print(f"  geometry: {gi['source']}")
        d = gi.get("declared_geometry")
        d_label = declared_geometry_method(d) if d else None
        if a.geometry_origin and d_label and \
                d_label != declared_geometry_method(a.geometry_origin):
            # the file's own declaration wins; a disagreeing flag is refused,
            # before any compute, rather than silently overriding it
            raise SystemExit(f"--geometry-origin {a.geometry_origin!r} contradicts the "
                             f"file's own declaration {d!r} ({d_label}); drop the flag, "
                             f"or fix the file's property")
        if d and not a.geometry_origin and d_label is None:
            print(f"note: the file declares its geometry as {d!r}, which is not a "
                  f"recognised origin; it is keyed as undeclared. Pass "
                  f"--geometry-origin to state it.", file=sys.stderr)
        if "DISCARDED" in gi["source"]:
            print("WARNING: this file is a FLAT 2D depiction, not a geometry. The "
                  "drawn layout has been replaced by a generated 3D conformer, so "
                  "the map does NOT reflect the coordinates in the file.",
                  file=sys.stderr)
    elif a.smiles:
        spec, gi = from_smiles(a.smiles, n_confs=a.nconfs, seed=a.seed)
        print(f"SMILES {a.smiles}  ->  {gi['formula']}, charge {gi['charge']:+d}, "
              f"{gi['forcefield']}")
        print(f"  {gi['n_converged']}/{gi['n_embedded']} conformers converged, "
              f"seed 0x{gi['seed']:x}")
    verified, computed, mapped_ok = None, None, None
    if gi is not None:
        # The gate inspected a separately parsed object. Check that the
        # COORDINATES about to be computed are the structure it approved --
        # otherwise the gate is a seal over a copy. Comparing gi["smiles"]
        # did not do this: it is the input string, or the reader's own output
        # (review 3, H2).
        kept = gi["source"].startswith("file 3D coordinates used as-is")
        if kept and a.optimize and not a.force:
            # checked before any compute: optimising a supplied pose discards
            # exactly what the pose was for (its placement in an environment)
            raise SystemExit("--optimize would move the coordinates this file "
                             "supplied (a pose or reference geometry). Drop "
                             "--optimize to compute the file's geometry, or "
                             "pass --force to optimise it anyway")
        if approved is not None:
            verified, computed, mapped_ok = check_identity(
                spec, probe, approved, kept, a.force, "before the SCF")
        if a.charge is None:
            charge = gi["charge"]
        elif a.charge != gi["charge"]:
            print(f"WARNING: --charge {a.charge} overrides the structure's formal "
                  f"charge {gi['charge']:+d}", file=sys.stderr)
        if gi["n_rot"] >= 1 and "as-is" not in gi["source"]:
            print(conformer_note(gi["n_rot"], gi["e_spread"]), file=sys.stderr)
        sym = [t.split()[0] for t in spec.split(";")]
        # F is deliberately excluded: fluorine is too electronegative to form
        # a sigma-hole, so reporting a "V_max on the F cap" invites a reader
        # to read noise as a halogen-bond donor.
        hal = [i for i, s in enumerate(sym) if s in ("Cl", "Br", "I")]
        i_sub = hal[0] if len(hal) == 1 else None
        if len(hal) > 1:
            sh_reason = "more than one Cl/Br/I"
        elif i_sub is None and "F" in sym:
            sh_reason = "fluorine does not form one"
            print("note: fluorine present; no sigma-hole is reported because F "
                  "does not form one.", file=sys.stderr)
    else:                           # --atoms / --molecule: gated above
        spec = atoms_spec
        # the coordinates are the input, so verify against them as given
        verified, computed, mapped_ok = check_identity(
            spec, probe, approved, True, a.force, "before the SCF")
        if a.atoms:
            i_sub, sh_reason = None, "input has no halogen index (--atoms)"
        else:
            i_sub = None if a.molecule == "benzene" else 6

    site_groups, sites_by = None, None
    if a.sites:                     # resolved now, so a bad file costs no SCF
        rows = [t.split() for t in spec.split(";")]
        maps = ({x.GetIdx(): x.GetAtomMapNum() for x in gi["mol"].GetAtoms()
                 if x.GetAtomMapNum()} if gi is not None else None)
        try:
            sites_by, site_groups = resolve_site_groups(
                load_site_groups(a.sites), [r[0] for r in rows],
                np.array([[float(v) for v in r[1:4]] for r in rows]), maps)
        except (OSError, ValueError) as e:
            raise SystemExit(f"--sites {a.sites}: {e}") from None

    natm, nao = estimate_cost(spec, a.basis, charge, spin)
    print(f"{natm} atoms, {nao} basis functions in {a.basis}")
    if nao > NAO_REFUSE and not a.yes:
        raise SystemExit(f"{nao} basis functions exceeds the {NAO_REFUSE} guard; "
                         f"an SCF this size is hours, not seconds. Use a smaller "
                         f"basis, or pass --yes.")
    if nao > NAO_WARN:
        print(f"WARNING: {nao} basis functions; expect many minutes of SCF "
              f"(253 nao took 165 s on one core).", file=sys.stderr)

    solv = None if a.solvent == "vacuum" else a.solvent
    mol, mf = run_scf(spec, basis=a.basis, xc=a.xc, charge=charge, spin=spin,
                      density_fit=a.density_fit, gpu=a.gpu,
                      solvent=solv)
    if a.optimize:
        from pyscf.geomopt.berny_solver import optimize as _opt
        print("optimising geometry...", flush=True)
        if a.gpu:
            # run_scf hands back a CPU object, so berny's gradients are CPU
            print("note: --optimize --gpu: each SCF runs on the GPU, but the "
                  "geometry gradients run on the CPU.", file=sys.stderr)
        pre_opt_sha1 = geometry_record(mol, charge, spin)["sha1"]
        mol = _opt(mf, maxsteps=50)
        spec = "; ".join(
            f"{mol.atom_symbol(i)} "
            f"{mol.atom_coord(i)[0] * BOHR:.6f} {mol.atom_coord(i)[1] * BOHR:.6f} "
            f"{mol.atom_coord(i)[2] * BOHR:.6f}" for i in range(mol.natm))
        if approved is not None:
            # the coordinates PySCF now receives are not the ones checked
            # (invariant 3): a proton can move, e.g. a gas-phase zwitterion.
            # reference=False: optimisation legitimately moves a file's atoms
            verified, computed, mapped_ok = check_identity(
                spec, probe, approved, False, a.force, "after --optimize")
        mol, mf = run_scf(spec, basis=a.basis, xc=a.xc, charge=charge,
                          density_fit=a.density_fit, gpu=a.gpu,
                          spin=spin, solvent=solv)
    dm = total_dm(mf)
    print(f"SCF {mf.e_tot:.6f} Ha, {mol.nao} basis functions, {a.xc}/{a.basis}, "
          f"{'PCM ' + a.solvent if solv else 'GAS PHASE'}"
          + (f", UNRESTRICTED (spin={spin}, <S^2>={mf.spin_square()[0]:.3f})"
             if spin else ""))
    electronic = electronic_report(mf, charge, solv)
    for p in electronic["problems"]:
        print(f"WARNING: {p}", file=sys.stderr)

    verts, faces, diag = extract_surface(mol, dm, a.spacing, iso=a.iso,
                                         margin=a.margin, strict=not a.no_strict)
    print(f"grid {diag['grid']}, {diag['n_verts']:,} vertices, chi={diag['chi']} "
          f"({'closed' if diag['closed'] else 'OPEN'}, genus {diag['genus']}, "
          f"{diag['components']} piece(s)), max rho on a box face {diag['face_max_rho']:.1e}")
    # F3: diagnostic only. This number moves -5%..+7% with grid spacing while
    # V_max moves 0.5%, and is IDENTICAL between def2-SVP and def2-TZVP where
    # V_max moves 3.5 kcal/mol. It cannot validate the quantity of interest.
    print(f"  [diagnostic, not a check] integral rho dV = "
          f"{diag['rho_integral']:.2f} vs {diag['nelec']} electrons")

    V = esp_at_points(mol, dm, verts) * HA2KCAL
    rf_added = bool(solv and a.include_reaction_field)
    if rf_added:
        rf = reaction_field_at_points(mf, verts) * HA2KCAL
        print(f"reaction field on the surface {rf.min():+.2f} .. {rf.max():+.2f} "
              f"kcal/mol, added to V")
        V = V + rf
    print(f"V = {V.min():+.2f} .. {V.max():+.2f} kcal/mol"
          + (f"  [{SOLVENT_ESP[rf_added]}]" if solv else ""))
    sh_rec, sh_reason = sigma_hole_record(verts, V, mol, i_sub, sh_reason)
    if sh_rec is not None:
        sh = sh_rec["V_max_kcal_mol"]
        print(f"sigma-hole V_max = {float('nan') if sh is None else sh:+.2f} "
              f"kcal/mol  ({sh_rec['n_cap_vertices']} vertices in the cap)")

    vlim, kind, notes = choose_scale(V, lim=a.lim, kind=a.scale,
                                     vmin=a.vmin, vmax=a.vmax)
    cmap_name, cmap_warning = cmap_after_scf(kind, a.cmap, a.force_cmap)
    if cmap_warning:
        notes.append(cmap_warning)
    for n in notes:
        print(("WARNING: " if n.startswith("WARNING") else "note: ")
              + n.removeprefix("WARNING: "), file=sys.stderr)
    if a.lim is None and kind == "diverging":
        print(f"WARNING: no --lim given, auto-scaling to this molecule's own "
              f"extremum ({vlim[0]:+.1f} to {vlim[1]:+.1f} kcal/mol).\n"
              f"         Panels made this way are NOT comparable across "
              f"molecules. Pass --lim for any multi-panel figure.", file=sys.stderr)
    out_of_range = 100.0 * float(np.mean((V < vlim[0]) | (V > vlim[1])))
    if out_of_range > 0:
        print(f"note: {out_of_range:.1f}% of surface vertices fall outside the "
              f"scale and are clamped to the end colours")

    cmaps = build_maps()
    if cmap_name not in cmaps:
        raise SystemExit(f"--cmap must be one of {sorted(cmaps)}")
    img = render(verts, faces, V, cmaps[cmap_name], vlim, size=a.size,
                 elev=a.elev, azim=a.azim, ambient=a.ambient,
                 diffuse=1.0 - a.ambient)
    plt.imsave(a.out, img)
    print(f"wrote {a.out}  ({kind} scale {vlim[0]:+.1f} to {vlim[1]:+.1f} "
          f"kcal/mol, cmap {cmap_name})")

    g_method = geometry_method(
        source=(gi or {}).get("source"), forcefield=(gi or {}).get("forcefield"),
        atoms=bool(a.atoms and gi is None),
        molecule=a.molecule if gi is None and not a.atoms else None,
        optimize=(a.xc, a.basis, a.solvent) if a.optimize else None,
        # the file's readable declaration, else the flag (a contradiction was
        # refused above); unreadable file text does not mask a valid flag
        declared=next((x for x in ((gi or {}).get("declared_geometry"), a.geometry_origin)
                       if x and declared_geometry_method(x)), None))
    a_fields = analysis_fields(a.xc, a.basis, a.solvent, a.iso, a.spacing,
                               a.density_fit, g_method, reaction_field=rf_added)
    r_fields = render_fields(cmap_name, kind, vlim, a.elev, a.azim, a.size, a.ambient)
    symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
    label = (gi or {}).get("formula") or a.smiles or a.file or a.molecule
    P_A, coords_A = verts * BOHR, mol.atom_coords() * BOHR
    surface = surface_report(diag, a.iso)
    # an atom order forced past verify_mapping invalidates the input too
    valid = validity(gate, False if mapped_ok is False else verified, surface,
                     electronic=electronic)
    # what the coordinates are; settings.geometry said "generated by ETKDGv3"
    # even after --optimize had replaced them (review 4 N1)
    geom_source = (gi["source"] if gi is not None else
                   "supplied with --atoms" if a.atoms else
                   f"built-in idealised geometry ({a.molecule})")
    if a.optimize:
        geom_source += (f"; then optimised with {a.xc}/{a.basis}"
                        f"{', PCM ' + a.solvent if solv else ''} (pyberny)")
    # Atom maps survive --optimize: berny does not reorder atoms, and the
    # post-optimisation verify_mapping has just confirmed the order (N1).
    atom_maps = ({x.GetIdx(): x.GetAtomMapNum() for x in gi["mol"].GetAtoms()
                  if x.GetAtomMapNum()} if gi is not None else None)
    if gi is not None and not a.optimize:
        h_source = gi["h_source"]
    else:                           # optimised, or --atoms / --molecule
        h_source = [("optimised" if a.optimize else "supplied") if s == "H" else None
                    for s in symbols]
    sites, sites_reason = sites_block(valid, P_A, faces, V, coords_A, symbols,
                                      input_reason=input_invalid_reason(gate, verified,
                                                                        mapped_ok),
                                      atom_maps=atom_maps, h_source=h_source,
                                      groups=site_groups, source=f"caller:{sites_by}")
    valid["sites"] = sites is not None
    if sites is not None:
        generated = gi is not None and "generated" in gi["source"]
        if a.optimize:
            relaxation = f"{a.xc}/{a.basis} (pyberny), all atoms, isolated molecule"
        elif generated:
            relaxation = f"ETKDGv3 + {gi['forcefield']}, all atoms, isolated molecule"
        elif gi is not None and gi["forcefield"] not in (None, "none"):
            relaxation = f"{gi['forcefield']}, isolated molecule"
        else:
            relaxation = "none"
        n_rot = gi["n_rot"] if gi is not None else None
        sites.update(site_context(
            solvent_esp=a_fields["solvent_esp"],
            heavy_atoms=("optimised" if a.optimize else
                         "generated" if generated else "supplied_frozen"),
            h_source=h_source, relaxation=relaxation, n_rot=n_rot, surface=surface,
            identity=dict(chemical_smiles=computed,
                          # the input's mapped structure; verify_mapping checked
                          # the geometry against it atom by atom
                          mapped_smiles=approved if atom_maps else None,
                          chemical_match=verified, mapped_match=mapped_ok),
            budget=budget_rows(basis=a.basis, solvent=solv, charge=charge, n_rot=n_rot,
                               geometry_method=g_method)))
    if a.html:
        esp_figure(verts, faces, V, symbols, mol.atom_coords(), cmaps[cmap_name],
                   vlim, title=f"{label} · {a.xc}/{a.basis}", kind=kind).write_html(
            a.html, include_plotlyjs=True)
        print(f"wrote {a.html}  (interactive 3D)")
    if a.npz:
        import json
        # the atom owning each triangle, per partition rule: enough to colour
        # or regroup by site without another SCF. Only when the record has a
        # site table (tri_atom_arrays).
        tri_atom = tri_atom_arrays(valid["sites"], P_A, faces, coords_A, symbols)
        maps = np.array([(atom_maps or {}).get(i, 0) for i in range(len(symbols))])
        np.savez_compressed(a.npz, verts=verts, faces=faces, V=V,
                            symbols=np.array(symbols), coords=mol.atom_coords(),
                            vlim=np.array(vlim), kind=kind, cmap=cmap_name,
                            units="coords, verts: bohr; V, vlim: kcal/mol",
                            xc=a.xc, basis=a.basis, solvent=a.solvent,
                            solvent_esp=str(a_fields["solvent_esp"]), iso=a.iso,
                            analysis_key=record_key(a_fields),
                            render_key=record_key(r_fields),
                            validity=json.dumps(valid),
                            surface_problems=json.dumps(list(diag["problems"])),
                            geometry_sha1=geometry_record(mol, charge, spin)["sha1"],
                            atom_maps=maps, **tri_atom)
        print(f"wrote {a.npz}  (surface, bohr; V in kcal/mol)")
    if a.json:
        import json
        rec = dict(
            schema=RECORD_SCHEMA,           # docs/record-schema.json; not hashed
            validity=valid,
            analysis_key=record_key(a_fields), render_key=record_key(r_fields),
            key_fields=dict(analysis=a_fields, render=r_fields),
            geometry=dict(geometry_record(mol, charge, spin),
                          optimised_from_sha1=pre_opt_sha1,
                          # the approved structure's bonds in this atom order,
                          # only when verify_mapping confirmed that order (A1)
                          bonds=record_bonds(probe) if mapped_ok else None),
            input=dict(smiles=a.smiles, file=a.file,
                       record=(gi or {}).get("record"),
                       # the structure the gate saw; "approved" only if it passed
                       gate_smiles=approved,
                       approved_smiles=approved if not (gate and gate["forced"]) else None,
                       gate=gate,
                       computed_smiles=computed, structure_verified=verified,
                       formula=(gi or {}).get("formula"), charge=charge, spin=spin),
            settings=dict(xc=a.xc, basis=a.basis, solvent=a.solvent,
                          solvent_esp=a_fields["solvent_esp"],
                          include_reaction_field=rf_added,
                          density_fit=a.density_fit, gpu=a.gpu,
                          gpu_validated=(not gpu_untested(spin, solv, a.density_fit)
                                         if a.gpu else None),
                          iso=a.iso,
                          spacing_bohr=a.spacing, margin_bohr=a.margin,
                          grids_level=int(mf.grids.level), conv_tol=float(mf.conv_tol),
                          seed=a.seed, nconfs=a.nconfs,
                          optimize=a.optimize, geometry=geom_source),
            results=dict(e_tot_Ha=float(mf.e_tot), nao=int(mol.nao),
                         dipole=dipole_about_nuclear_centre(mol, dm),
                         V_min_kcal_mol=float(V.min()), V_max_kcal_mol=float(V.max()),
                         sigma_hole=sh_rec, sigma_hole_reason=sh_reason,
                         electronic=electronic,
                         n_verts=int(diag["n_verts"]), surface=surface,
                         surface_stats=area_stats(P_A, faces, V),
                         gipf=gipf_stats(P_A, faces, V)),
            sites=sites, sites_reason=sites_reason,
            scale=dict(kind=kind, vmin=float(vlim[0]), vmax=float(vlim[1]),
                       sign_change_ratio=sign_change_ratio(V),
                       sign_change_threshold=SIGN_CHANGE_MIN,
                       cmap=cmap_name, cmap_requested=a.cmap,
                       cmap_substituted=bool(cmap_warning),
                       cmap_forced=bool(a.cmap and a.force_cmap),
                       **scale_provenance(kind, a.lim, a.vmin, a.vmax)),
            provenance=provenance(),
            caveat=("Measured error budget (README): basis def2-SVP -> "
                    "def2-TZVPP moves the chlorobenzene sigma-hole by +3.54 "
                    "kcal/mol; solvent gas -> water by -3.60. Quote numbers only "
                    "against those margins."))
        with open(a.json, "w") as f:
            json.dump(rec, f, indent=2)
        print(f"wrote {a.json}")


if __name__ == "__main__":
    main()
