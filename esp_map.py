#!/usr/bin/env python3
"""
esp_map.py v4 — electrostatic-potential maps for organic molecules,
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
      structure is now asserted against the computed one.
    Fluorine no longer reports a sigma-hole, because it does not form one.

NEW IN v3 — the structure check (--check)

    A drawing can be perfectly legible to a chemist and still be ambiguous to
    a machine, and the pipeline will happily compute a flawless ESP map of a
    molecule you did not mean. Every structure is therefore gated BEFORE any
    geometry or SCF work:

      REJECTED  valence/charge errors, R-groups and query atoms, radicals
                (this is restricted DFT), and disconnected fragments -- a
                salt has no single ESP surface and the counter-ion dominates
                it (sodium acetate forced through reaches +113 kcal/mol,
                against +61 for acetic acid's own proton).
      WARNED    undrawn E/Z geometry; two or more stereocentres with no
                wedges; elements outside the tested set.
      NOTED     one unspecified stereocentre (the two enantiomers have
                mirror-image ESP -- identical values, reflected picture),
                reachable tautomers, net charge, isotopes.

    The stereochemistry wording is measured, not assumed: threonine's two
    diastereomers differ in ESP extrema by ~0.9 kcal/mol against a
    0.4 kcal/mol conformer-selection noise floor. Stereochemistry changes the
    SHAPE of the surface far more than the numbers on it.

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
        B3LYP/def2-SVP optimised +9.59            (C-Cl 1.7512 A)

    a 0.35 kcal/mol spread, an order of magnitude below the basis-set error.
    A force-field geometry is therefore adequate for this quantity; --optimize
    runs a DFT optimisation if you want it anyway.

    ETKDG is stochastic, so --seed is fixed by default: without it the same
    SMILES gives a different conformer every run. Measured conformer effect on
    the ESP extrema, lowest- vs highest-energy MMFF conformer:

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
    3. marching cubes at rho = 0.002 e/bohr^3 (Bader surface, ~98% of rho)
    4. ESP evaluated ONLY at the isosurface vertices, not on a full grid
       -- ~10^4 points instead of ~10^6, so this costs seconds, and the
          potential is evaluated exactly where it is displayed rather than
          interpolated onto the surface from a cube file
    5. z-buffered render with deliberately weak shading, so the colormap's
       lightness carries the data instead of fighting the lighting

COLOR CONVENTION
    Every colormap is oriented warm/red = NEGATIVE V = electron-rich, the
    entrenched chemistry reading. Fixing perception should not silently flip
    the sign convention too.

WHAT THIS IS AND IS NOT
    An instrument for making legible pictures and for order-of-magnitude
    reasoning. Quantitative V_max values carry, at minimum, a basis-set
    uncertainty; gas phase throughout; geometries are idealised unless you
    supply your own. Do not quote a number from it against a margin smaller
    than the convergence tests reported in the review.
"""
import argparse
import os
import sys

import numpy as np

BOHR = 0.52917721092
HA2KCAL = 627.5095
ISO_BADER = 0.002

# basis sets at or below double zeta; penalty measured in the review
SMALL_BASIS = {"sto-3g", "3-21g", "6-31g", "6-31g*", "def2-svp", "def2-sv(p)"}

# dielectric constants for --solvent. Measured effect on the chlorobenzene
# sigma-hole: gas +9.83 -> chloroform +7.14 -> water +6.22 kcal/mol, i.e. the
# solvent axis is CO-DOMINANT with the basis-set axis (+3.55). Do not quote a
# condensed-phase number from a gas-phase run.
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
    best, ff, spread, n_emb, n_conv = _embed_and_minimise(mol, n_confs, seed,
                                                          max_iters=max_iters)
    info = dict(forcefield=ff, n_embedded=n_emb, n_converged=n_conv,
                e_spread=spread,
                charge=Chem.GetFormalCharge(mol),
                n_radical=int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())),
                formula=Chem.rdMolDescriptors.CalcMolFormula(mol),
                n_rot=int(AllChem.CalcNumRotatableBonds(mol)), seed=seed,
                source="3D generated from SMILES by ETKDGv3", n_found=1,
                smiles=smiles, mol=mol)
    return _mol_to_spec(mol, best), info


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
    """ETKDGv3 + MMFF94s over n_confs; returns (best_conf_id, forcefield, spread)."""
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
    return (int(cids[int(np.argmin(en))]), ff, float(fin.max() - fin.min()),
            len(cids), int((ok == 0).sum()))


def from_file(path, n_confs=20, seed=0xC0FFEE, charge=0):
    """Read .mol/.sdf/.cdxml/.xyz/.pdb from a drawing program and return
    (PySCF atom string, info).

    The load-bearing distinction is 2D vs 3D. A structure drawn in ChemDraft,
    ChemDraw, Ketcher or MarvinSketch and exported as MOL/SDF/CDXML carries a
    FLAT depiction: every z is zero and it is a picture, not a geometry.
    Feeding that to an SCF would produce a planar molecule and a confidently
    wrong map. RDKit records the dimensionality on line 2 of a molblock and
    exposes it as Conformer.Is3D(); we check that AND the coordinates, because
    a file can claim 3D and be flat.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    # ONE reader, shared with the gate. v3 had a second copy here, so a fix
    # applied to parse_for_check (XYZ bond perception) silently did not reach
    # the path that actually builds the geometry.
    mol, problems, n_found = read_structure(path, charge=charge)
    if mol is None:
        raise ValueError(problems[0][2] if problems else
                         f"RDKit could not read any molecule from {path}")

    has_conf = mol.GetNumConformers() > 0
    is3d = bool(has_conf and mol.GetConformer().Is3D())
    flat = True
    if has_conf:
        c = mol.GetConformer()
        flat = max(abs(c.GetAtomPosition(i).z)
                   for i in range(mol.GetNumAtoms())) < 1e-6
    n_h_before = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")

    if is3d and not flat:
        n_supplied = mol.GetNumAtoms()
        mol = Chem.AddHs(mol, addCoords=True)
        n_h_added = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H") - n_h_before
        source = "file 3D coordinates used as-is"
        ff, spread, n_emb, n_conv = "none", 0.0, 1, 1
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
            source += f"; {n_h_added} hydrogen(s) added"
        best = -1
    else:
        mol = Chem.AddHs(mol)
        best, ff, spread, n_emb, n_conv = _embed_and_minimise(mol, n_confs, seed)
        source = ("2D depiction: drawn coordinates DISCARDED, 3D generated "
                  "by ETKDGv3" if (has_conf and (flat or not is3d))
                  else "no coordinates in file, 3D generated by ETKDGv3")

    info = dict(forcefield=ff, n_embedded=n_emb, n_converged=n_conv,
                e_spread=spread, charge=Chem.GetFormalCharge(mol),
                n_radical=int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())),
                formula=Chem.rdMolDescriptors.CalcMolFormula(mol),
                n_rot=int(AllChem.CalcNumRotatableBonds(mol)), seed=seed,
                source=source, n_found=n_found,
                smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)), mol=mol)
    return _mol_to_spec(mol, best), info


# ------------------------------------------------------- admissibility gate
# Elements this pipeline has actually been exercised on. Wider is not refused,
# only flagged: def2 covers H-Rn, but closed-shell RKS on a transition metal is
# usually wrong for reasons that have nothing to do with the basis.
TESTED_ELEMENTS = set("H B C N O F Si P S Cl Br I Se".split())

BLOCK, WARN, INFO = "BLOCK", "WARN", "INFO"


def check_structure(mol, basis="def2-tzvp", spin=0):
    """Can this drawing be translated into SMILES unambiguously, and is the
    result a molecule this pipeline can compute?

    A drawn structure can be perfectly legible to a chemist and still be
    ambiguous to a machine. Each finding below names what to change IN THE
    DRAWING, because that is where the fix belongs.

    Returns (verdict, findings) where findings is a list of
    (severity, code, message, remedy).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    out = []

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
    rad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    if rad and spin == 0:
        out.append((BLOCK, "OPEN_SHELL",
                    f"{rad} unpaired electron(s): this is a radical, and no "
                    f"--spin was given, so it would be computed with a "
                    f"RESTRICTED wavefunction that cannot represent it",
                    f"pass --spin {rad} to use unrestricted DFT (UKS), or "
                    f"complete the valences / add the missing charge"))
    elif rad:
        out.append((INFO, "OPEN_SHELL_UNRESTRICTED",
                    f"{rad} unpaired electron(s); --spin {spin} selected, so "
                    f"UKS will be used"
                    + ("" if spin == rad else
                       f". NOTE: --spin {spin} does not match the {rad} radical "
                       f"electron(s) RDKit perceives"),
                    "none needed" if spin == rad else "check the multiplicity"))
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
                   f"ESP EXTREMA differ by only ~0.9 kcal/mol (against a "
                   f"0.4 kcal/mol conformer-noise floor), but the SHAPE of the "
                   f"surface and where the features sit on it do differ.")
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


def parse_for_check(smiles=None, path=None, charge=0):
    """Parse to an RDKit mol WITHOUT generating 3D coordinates.

    The admissibility gate must run before embedding: a radical or an R-group
    structure should be rejected in milliseconds, not after twenty conformers
    have been generated and a force field has complained about them.
    Returns (mol or None, list of BLOCK findings explaining a parse failure).
    """
    from rdkit import Chem
    _quiet_rdkit()
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
    mol, problems, _ = read_structure(path, charge=charge)
    return mol, problems


def read_structure(path, charge=0):
    """The single structure reader. Returns (mol or None, problems, n_found).

    charge is needed only for XYZ, which carries coordinates and no bonds:
    bond perception has to be told the net charge or it fails. v4.0 hard-coded
    0 while telling the user to \"set the charge explicitly\", naming a remedy
    the code did not implement."""
    from rdkit import Chem
    _quiet_rdkit()
    ext = str(path).lower().rsplit(".", 1)[-1]
    try:
        if ext == "cdxml":
            mols = [m for m in Chem.MolsFromCDXMLFile(str(path)) if m is not None]
        elif ext == "sdf":
            mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False)
                    if m is not None]
        elif ext == "mol":
            mols = [Chem.MolFromMolFile(str(path), removeHs=False)]
        elif ext == "xyz":
            # An XYZ file has coordinates and NO connectivity. RDKit's reader
            # returns a bondless mol, which every downstream check reads as
            # one fragment per atom: v3 rejected EVERY xyz file as
            # "MULTIPLE_FRAGMENTS", with a message about salts.
            m = Chem.MolFromXYZFile(str(path))
            if m is not None:
                from rdkit.Chem import rdDetermineBonds
                try:
                    rdDetermineBonds.DetermineBonds(m, charge=charge)
                    Chem.SanitizeMol(m)      # AddHs needs implicit valences
                except Exception as e:
                    return None, [(BLOCK, "XYZ_CONNECTIVITY",
                                   f"could not perceive bonds from coordinates "
                                   f"at charge {charge:+d}: {e}",
                                   f"pass --charge if the species is not neutral, "
                                   f"or supply a format that carries connectivity "
                                   f"(MOL/SDF/CDXML)")], 0
            mols = [m]
        elif ext == "pdb":
            mols = [Chem.MolFromPDBFile(str(path), removeHs=False)]
        else:
            return None, [(BLOCK, "UNSUPPORTED_FORMAT",
                           f"'.{ext}' is not a format this reads",
                           "export as MOL, SDF, CDXML, XYZ or PDB")], 0
    except Exception as e:
        return None, [(BLOCK, "READ_ERROR", str(e), "check the file is intact")], 0
    mols = [m for m in mols if m is not None]
    if not mols:
        return None, [(BLOCK, "UNPARSEABLE",
                       f"RDKit read no valid molecule from {path}",
                       "re-export from the drawing program; check for R-groups "
                       "or query atoms the exporter could not represent")], 0
    return mols[0], [], len(mols)


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


def run_scf(atom, basis="def2-tzvp", xc="b3lyp", charge=0, spin=0, verbose=0,
            solvent=None, density_fit=False, gpu=False):
    """Restricted for a closed shell, UNRESTRICTED when spin != 0.

    v3 built an RKS object whatever the spin. PySCF does not refuse: on
    O with spin=2 (nelec 5,3) RKS returned -74.998688 Ha against UKS
    -75.000759, i.e. a converged, plausible, physically meaningless number
    from a restricted wavefunction describing an open shell.

    density_fit is OFF by default and the error budget was measured without
    it. Its effect on the ESP is unmeasured until bench/timing.py says so.

    gpu runs the SCF ONLY on GPU4PySCF and hands back a CPU object, so
    total_dm and everything downstream are unchanged. Validated once
    (bench/gpu_l4.json): RKS + density_fit, gas phase, decanoic acid, where
    V_min/V_max matched the CPU to <0.001 kcal/mol. UKS, PCM and runs
    without density fitting are untested on the GPU.
    """
    if gpu:
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
    if basis.lower() in SMALL_BASIS:
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

    diag = dict(grid=tuple(int(x) for x in n), n_verts=len(verts),
                face_max_rho=float(face_max), **topo,
                rho_integral=float(rho.sum() * spacing ** 3),
                nelec=int(mol.nelectron), margin=margin, spacing=spacing)

    problems, notes = surface_problems(topo, face_max, iso)
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


def choose_scale(V, lim=None, kind="auto", cmap_name="vik", vmin=None, vmax=None):
    """Decide the colour scale from the data, not from habit.

    A diverging scale symmetric about zero is right only if the surface
    actually CHANGES SIGN. For an ion it does not: benzoate runs -158.7 to
    -40.7 kcal/mol, so centring on zero maps the whole molecule to one end
    and hides 118 kcal/mol of structure.

    Returns (vlim, kind, notes).
    """
    data_min, data_max = float(V.min()), float(V.max())
    spans_zero = data_min < 0.0 < data_max
    notes = []
    if kind == "auto":
        kind = "diverging" if spans_zero else "sequential"
        if not spans_zero:
            notes.append(f"surface does not change sign ({data_min:+.1f} to {data_max:+.1f} "
                         f"kcal/mol): using a SEQUENTIAL scale over the actual "
                         f"range. A diverging scale centred on zero would put the "
                         f"whole molecule at one end of the ramp.")
    if kind == "diverging":
        if not spans_zero:
            notes.append(f"WARNING: diverging scale requested but V never changes "
                         f"sign ({data_min:+.1f} to {data_max:+.1f}); every vertex falls on "
                         f"one side of the midpoint colour.")
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
        if cmap_name not in SEQUENTIAL:
            notes.append(f"note: '{cmap_name}' is a diverging map being used on a "
                         f"sequential scale; prefer lajolla/batlow/viridis.")
    return vlim, kind, notes


# ------------------------------------------------------------------ molecules
MOLECULES = {"bromobenzene": ("Br", 1.890), "chlorobenzene": ("Cl", 1.740),
             "iodobenzene": ("I", 2.095), "fluorobenzene": ("F", 1.354),
             "benzene": ("H", 1.080)}


def scale_provenance(kind, lim=None, vmin=None, vmax=None):
    """Where the colour limits came from, and whether panels can be compared.

    Comparable means the limits were fixed by the user, so another molecule
    run with the same flags maps a colour to the same kcal/mol. The first
    --json marked every sequential scale comparable; an auto-scaled one (acetate
    -170..-80 next to another anion at -130..-30) is not. --lim only fixes a
    DIVERGING scale; --vmin/--vmax only fix a SEQUENTIAL one, and only as a
    pair (choose_scale ignores the others for that kind).
    """
    if kind == "diverging":
        fixed = lim is not None
    else:
        fixed = vmin is not None and vmax is not None
    return dict(source="user" if fixed else "auto", comparable=fixed)


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


def esp_figure(verts, faces, V, symbols, coords, cmap, vlim, title=""):
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
        colorbar=dict(title="V (kcal/mol)<br>warm = negative", len=0.75),
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


def sigma_hole(verts, V, mol, i_sub=6, pr_min=0.3, pe_max=1.6, i_bond=None):
    """V_max on the outer cap of the substituent, along the X-to-neighbour axis.

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
        d = np.linalg.norm(co - X, axis=1)
        d[i_sub] = np.inf
        i_bond = int(d.argmin())
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
    ap.add_argument("--seed", type=lambda s: int(s, 0), default=0xC0FFEE,
                    help="ETKDG random seed; fixed so runs are reproducible")
    ap.add_argument("--nconfs", type=int, default=20,
                    help="conformers to embed; the lowest-energy one is used")
    ap.add_argument("--optimize", action="store_true",
                    help="DFT geometry optimisation after the force field "
                         "(slow; measured worth 0.35 kcal/mol on a sigma-hole)")
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
    ap.add_argument("--basis", default="def2-tzvp")
    ap.add_argument("--density-fit", action="store_true",
                    help="RI-J/K: much faster SCF; ESP effect not yet measured")
    ap.add_argument("--gpu", action="store_true",
                    help="SCF on GPU4PySCF (NVIDIA only); validated with --density-fit")
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--solvent", default="vacuum",
                    choices=sorted(SOLVENTS),
                    help="PCM continuum. Measured: gas->water moves the "
                         "chlorobenzene sigma-hole by -3.6 kcal/mol, as large "
                         "as the basis-set error. Default vacuum, stated in "
                         "the output either way")
    ap.add_argument("--iso", type=float, default=ISO_BADER)
    ap.add_argument("--spacing", type=float, default=0.26,
                    help="bohr; V_max converged to 0.04 kcal/mol at this value")
    ap.add_argument("--margin", type=float, default=5.0,
                    help="bohr of vacuum around the nuclei; below 4 risks clipping")
    ap.add_argument("--no-strict", action="store_true",
                    help="render even if the isosurface fails its integrity check")
    ap.add_argument("--cmap", default="vik")
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
    a = ap.parse_args()

    if a.selftest:
        raise SystemExit(selftest())
    if a.html:
        try:                        # fail before the SCF, not after it
            import plotly  # noqa: F401
        except ImportError:
            raise SystemExit("--html needs plotly: pip install plotly") from None

    charge, spin, i_sub = a.charge or 0, a.spin, 6
    gi = None
    if a.smiles and a.file:
        raise SystemExit("--smiles and --file are mutually exclusive")

    if a.smiles or a.file:
        probe, probs = parse_for_check(smiles=a.smiles, path=a.file,
                                       charge=a.charge or 0)
        if probe is None:
            print_check(BLOCK, probs, a.smiles or str(a.file))
            raise SystemExit(1)
        from rdkit import Chem
        verdict, findings = check_structure(Chem.RemoveHs(probe), a.basis,
                                            spin=a.spin)
        print_check(verdict, findings, Chem.MolToSmiles(Chem.RemoveHs(probe)))
        if verdict == BLOCK and not a.force:
            raise SystemExit("structure rejected; fix the drawing, or pass "
                             "--force to compute anyway")
        if a.check:
            raise SystemExit(0 if verdict != BLOCK else 1)
    elif a.check:
        raise SystemExit("--check needs --smiles or --file")

    approved = None
    if a.smiles or a.file:
        from rdkit import Chem as _C
        approved = _C.MolToSmiles(_C.RemoveHs(probe))

    if a.file:
        spec, gi = from_file(a.file, n_confs=a.nconfs, seed=a.seed,
                             charge=a.charge or 0)
        print(f"{a.file}  ->  {gi['formula']}, charge {gi['charge']:+d}")
        print(f"  SMILES as read: {gi['smiles']}")
        print(f"  geometry: {gi['source']}")
        if gi["n_found"] > 1:
            print(f"WARNING: file contains {gi['n_found']} molecules; using the "
                  f"first. Split the file to map the others.", file=sys.stderr)
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
    if gi is not None:
        # The gate inspected a separately parsed object. Assert that the
        # structure about to be computed is the one that was approved --
        # otherwise the gate is a seal over a copy.
        if approved is not None:
            from rdkit import Chem as _C2
            computed = _C2.CanonSmiles(gi["smiles"])
            if computed != approved:
                raise SystemExit(
                    f"the structure that passed the check is not the one about "
                    f"to be computed:\n  approved: {approved}\n  computing: "
                    f"{computed}\nrefusing to continue")
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
        if i_sub is None and sum(1 for s in sym if s == "F"):
            print("note: fluorine present; no sigma-hole is reported because F "
                  "does not form one.", file=sys.stderr)
    elif a.atoms:
        spec, i_sub = a.atoms, None
    else:
        spec = monosubstituted_benzene(*MOLECULES[a.molecule])
        i_sub = None if a.molecule == "benzene" else 6

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
        mol = _opt(mf, maxsteps=50)
        spec = "; ".join(
            f"{mol.atom_symbol(i)} "
            f"{mol.atom_coord(i)[0] * BOHR:.6f} {mol.atom_coord(i)[1] * BOHR:.6f} "
            f"{mol.atom_coord(i)[2] * BOHR:.6f}" for i in range(mol.natm))
        mol, mf = run_scf(spec, basis=a.basis, xc=a.xc, charge=charge,
                          density_fit=a.density_fit, gpu=a.gpu,
                          spin=spin, solvent=solv)
    dm = total_dm(mf)
    print(f"SCF {mf.e_tot:.6f} Ha, {mol.nao} basis functions, {a.xc}/{a.basis}, "
          f"{'PCM ' + a.solvent if solv else 'GAS PHASE'}"
          + (f", UNRESTRICTED (spin={spin}, <S^2>={mf.spin_square()[0]:.3f})"
             if spin else ""))

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
    print(f"V = {V.min():+.2f} .. {V.max():+.2f} kcal/mol")
    if i_sub is not None:
        sh, ns = sigma_hole(verts, V, mol, i_sub=i_sub)
        print(f"sigma-hole V_max = {sh:+.2f} kcal/mol  ({ns} vertices in the cap)")

    if a.lim == 0.0:
        raise SystemExit("--lim 0 is not a range")
    vlim, kind, notes = choose_scale(V, lim=a.lim, kind=a.scale, cmap_name=a.cmap,
                                     vmin=a.vmin, vmax=a.vmax)
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
    if a.cmap not in cmaps:
        raise SystemExit(f"--cmap must be one of {sorted(cmaps)}")
    img = render(verts, faces, V, cmaps[a.cmap], vlim, size=a.size,
                 elev=a.elev, azim=a.azim, ambient=a.ambient,
                 diffuse=1.0 - a.ambient)
    plt.imsave(a.out, img)
    print(f"wrote {a.out}  ({kind} scale {vlim[0]:+.1f} to {vlim[1]:+.1f} "
          f"kcal/mol, cmap {a.cmap})")

    symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
    label = (gi or {}).get("formula") or a.smiles or a.file or a.molecule
    if a.html:
        esp_figure(verts, faces, V, symbols, mol.atom_coords(), cmaps[a.cmap],
                   vlim, title=f"{label} · {a.xc}/{a.basis}").write_html(
            a.html, include_plotlyjs=True)
        print(f"wrote {a.html}  (interactive 3D)")
    if a.npz:
        np.savez_compressed(a.npz, verts=verts, faces=faces, V=V,
                            symbols=np.array(symbols), coords=mol.atom_coords(),
                            vlim=np.array(vlim), kind=kind, cmap=a.cmap)
        print(f"wrote {a.npz}  (surface, bohr; V in kcal/mol)")
    if a.json:
        import json
        import subprocess
        import pyscf
        import rdkit
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                text=True, cwd=os.path.dirname(os.path.abspath(
                                    __file__))).stdout.strip() or None
        rec = dict(
            input=dict(smiles=a.smiles, file=a.file, approved_smiles=approved,
                       formula=(gi or {}).get("formula"), charge=charge, spin=spin),
            settings=dict(xc=a.xc, basis=a.basis, solvent=a.solvent,
                          density_fit=a.density_fit, gpu=a.gpu, iso=a.iso,
                          spacing_bohr=a.spacing, seed=a.seed, nconfs=a.nconfs,
                          optimize=a.optimize, geometry=(gi or {}).get("source")),
            results=dict(e_tot_Ha=float(mf.e_tot), nao=int(mol.nao),
                         dipole_D=float(np.linalg.norm(mf.dip_moment(
                             unit="Debye", verbose=0))),
                         V_min=float(V.min()), V_max=float(V.max()),
                         n_verts=int(diag["n_verts"]), chi=int(diag["chi"])),
            scale=dict(kind=kind, vmin=float(vlim[0]), vmax=float(vlim[1]),
                       cmap=a.cmap, **scale_provenance(kind, a.lim, a.vmin, a.vmax)),
            versions=dict(pyscf=pyscf.__version__, rdkit=rdkit.__version__,
                          espmap_commit=commit),
            caveat=("Measured error budget (README): basis def2-SVP -> "
                    "def2-TZVPP moves the chlorobenzene sigma-hole by +3.55 "
                    "kcal/mol; solvent gas -> water by -3.61. Quote numbers only "
                    "against those margins."))
        with open(a.json, "w") as f:
            json.dump(rec, f, indent=2)
        print(f"wrote {a.json}")


if __name__ == "__main__":
    main()
