"""P4 calibration set (docs/p4-spec.md section 6): reactions whose reactive
site is known from the literature, with what ALIE is expected to show,
committed BEFORE ALIE is implemented. Run from the repository root:

    python tests/data/p4_calibration/make_set.py

It rewrites tests/data/p4_calibration/*.sdf and manifest.json
deterministically.

Sites are heavy atoms, by 0-based index in the file's atom order (a site
includes the atom's hydrogens, as in the per-site partition). Statistic: the
site's minimum ALIE under the nearest-nucleus rule. Expectation types:
  - lower:            min over sites `a` < min over sites `b` (ALIE points to a);
  - not_clearly_lower: NOT (min over `a` < min over `b` - tol): ALIE must not
                      be claimed to pick `a`;
  - report_only:      no expectation; the literature says conditions decide.

`evidence` says how far each source was checked: Crossref metadata and the
abstract (`abstract`), or the full text (`full_text`). An expectation whose
evidence is `abstract` must be confirmed from the full text before ALIE's
result on it is quoted.
"""
import json
import shutil
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
for p in OUT.glob("*.sdf"):
    p.unlink()

SRC = {
    "mayr2011": dict(doi="10.1002/anie.201007100", evidence="full_text",
                     ref="Mayr, Breugst & Ofial (2011) Farewell to the HSAB Treatment of "
                         "Ambident Reactivity. Angew. Chem. Int. Ed. 50, 6470-6505"),
    "loos2003": dict(doi="10.1021/ja037317u", evidence="full_text",
                     ref="Loos, Kobayashi & Mayr (2003) Ambident Reactivity of the Thiocyanate "
                         "Anion Revisited. J. Am. Chem. Soc. 125, 14126-14132"),
    "kornblum1955": dict(doi="10.1021/ja01628a064", evidence="full_text",
                         ref="Kornblum, Smiley, Blackwood & Iffland (1955) The Mechanism of the "
                             "Reaction of Silver Nitrite with Alkyl Halides. J. Am. Chem. Soc. 77, 6269"),
    "tishkov2005": dict(doi="10.1002/anie.200501274", evidence="full_text",
                        ref="Tishkov, Schmidhammer, Roth, Riedle & Mayr (2005) Ambident Reactivity "
                            "of the Nitrite Ion Revisited. Angew. Chem. Int. Ed. 44, 4623-4626"),
    "kornblum1963": dict(doi="10.1021/ja00891a024", evidence="full_text",
                         ref="Kornblum, Berrigan & le Noble (1963) Solvation as a Factor in the "
                             "Alkylation of Ambident Anions. J. Am. Chem. Soc. 85, 1141-1147"),
    "stock1963": dict(doi="10.1016/S0065-3160(08)60277-4", evidence="metadata",
                      ref="Stock & Brown (1963) A Quantitative Treatment of Directive Effects in "
                          "Aromatic Substitution. Adv. Phys. Org. Chem. 1, 35-154"),
    "sjoberg1990": dict(doi="10.1139/v90-220", evidence="full_text",
                        ref="Sjoberg, Murray, Brinck & Politzer (1990) Average local ionization "
                            "energies on the molecular surfaces of aromatic systems as guides to "
                            "chemical reactivity. Can. J. Chem. 68, 1440-1443"),
    "murray1992": dict(doi="10.1016/0166-1280(92)85015-D", evidence="full_text",
                       ref="Murray, Brinck & Politzer (1992) Applications of calculated local "
                           "surface ionization energies to chemical reactivity. "
                           "J. Mol. Struct. THEOCHEM 255, 271-281"),
}
ANION_SETTINGS = dict(xc="b3lyp", basis="def2-tzvp", solvent="water")
NEUTRAL_SETTINGS = dict(xc="b3lyp", basis="def2-tzvp", solvent="vacuum")
CASES = []


def fixture(name):
    """A DFT geometry already committed for P3 (B3LYP/def2-TZVP, PCM water)."""
    shutil.copyfile(ROOT / "bench" / "fixtures" / f"{name}.sdf", OUT / f"{name}.sdf")
    return f"{name}.sdf", f"bench/fixtures/{name}.sdf (B3LYP/def2-TZVP, PCM water, DFT-optimised)"


def mmff(name, smiles, charge_note=""):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=5000)
    m.SetProp("_Name", name)
    m.SetProp("geometry", "mmff94s")
    w = Chem.SDWriter(str(OUT / f"{name}.sdf"))
    w.write(m)
    w.close()
    return f"{name}.sdf", "RDKit ETKDG + MMFF94s, seed 0xC0FFEE" + charge_note


def case(name, smiles, file_and_source, settings, reaction, sources, expect, sites, why):
    f, prov = file_and_source
    m = next(iter(Chem.SDMolSupplier(str(OUT / f), removeHs=False)))
    # the file's atom order must be the SMILES order, so site indices mean what they say
    ref = Chem.MolFromSmiles(smiles)
    for i in range(ref.GetNumAtoms()):
        assert m.GetAtomWithIdx(i).GetSymbol() == ref.GetAtomWithIdx(i).GetSymbol(), (name, i)
    CASES.append(dict(name=name, smiles=smiles, file=f, geometry_source=prov,
                      settings=settings, reaction=reaction,
                      sources=[dict(key=k, **SRC[k]) for k in sources],
                      sites=sites, expect=expect, why=why))


# --- ambident anions -----------------------------------------------------------
# Full texts read 2026-10-23. Only thiocyanate has one kinetic site that holds
# across conditions; for the enolate and nitrite the site depends on the
# counterion, solvent and alkylating agent, so they are report_only.
case("thiocyanate", "[S-]C#N", fixture("thiocyanate"), ANION_SETTINGS,
     "kinetic control favours S: kS/kN = 1e2-1e3 in SN2 alkylations and 1e3-1e4 for "
     "activation-controlled carbocations; R-NCS forms only under thermodynamic control "
     "(isothiocyanates are more stable) or with diffusion-controlled cations (kS/kN 2-10)",
     ["loos2003", "mayr2011"], dict(type="lower", a=["S"], b=["N"]),
     dict(S=[0], C=[1], N=[2]),
     "ESP points to N; the 2026-10-22 pilot (prototype ALIE) found its minimum on S "
     "(7.10 vs 8.99 eV)")
case("enolate", "C=C[O-]", fixture("enolate"), ANION_SETTINGS,
     "free enolate: O-alkylation is kinetically favoured (lower intrinsic barrier; Houk & "
     "Paddon-Row and Lee et al., reviewed in Mayr 2011 section 3.4.2), C-alkylation is "
     "thermodynamically favoured; synthetic C-alkylation reflects metal-enolate aggregates, "
     "counterion and solvent",
     ["mayr2011"], dict(type="report_only"),
     dict(C_CH2=[0], C_CO=[1], O=[2]),
     "no single kinetic site for a free enolate in PCM water; the pilot's ALIE minimum on "
     "CH2 (5.71 vs 7.47 eV on O) is a property to report, not a prediction to grade")
case("nitrite", "N(=O)[O-]", fixture("nitrite"), ANION_SETTINGS,
     "SN2 O/N ratio depends on solvent and reagent: MeOSO2Me + KNO2/18-crown-6 30:70 "
     "(EtOH) to 92:8 (THF); MeI 54:46 (DMF, DMSO); NaNO2 + primary/secondary RBr in DMF "
     "55-65 % nitroalkane, 30-33 % nitrite ester; diffusion-controlled carbocations favour O",
     ["tishkov2005", "kornblum1955", "mayr2011"], dict(type="report_only"),
     dict(N=[0], O=[1, 2]),
     "the selectivity itself moves with conditions, so no descriptor can be graded on one "
     "site; the pilot's ALIE tie (O 7.86, N 7.89 eV) is consistent with that")
# --- electrophilic aromatic substitution (classic surface-ALIE use) ------------
# Sjoberg 1990 (HF/6-31G*, STO-3G geometries, rho = 0.002) reports surface ALIE
# minima over o/p for aniline, toluene and fluorobenzene and over meta for
# benzaldehyde and nitrobenzene; Murray 1992 adds pyridine (meta) and warns
# that a ranking can change with the basis (pyridine N-oxide: para at 6-31G*,
# ortho at 3-21G). Their ALIE is Hartree-Fock (Koopmans); espmap's is
# Kohn-Sham, so the functional row of the error budget applies.
case("toluene", "Cc1ccccc1", mmff("toluene", "Cc1ccccc1"), NEUTRAL_SETTINGS,
     "EAS on toluene: ortho/para directing, activated",
     ["stock1963", "sjoberg1990"], dict(type="lower", a=["ortho", "para"], b=["meta"]),
     dict(ortho=[2, 6], meta=[3, 5], para=[4]),
     "Sjoberg 1990 Table 1: 11.93 (ortho), 11.96 (para), 12.08 (meta) eV; a narrow gap")
case("fluorobenzene", "Fc1ccccc1", mmff("fluorobenzene", "Fc1ccccc1"), NEUTRAL_SETTINGS,
     "EAS on fluorobenzene: ortho/para directing, deactivated",
     ["stock1963", "sjoberg1990"], dict(type="lower", a=["ortho", "para"], b=["meta"]),
     dict(ortho=[2, 6], meta=[3, 5], para=[4]),
     "Sjoberg 1990 Table 1: minima over para (12.31) and ortho (12.44 eV), none over meta")
case("nitrobenzene", "O=[N+]([O-])c1ccccc1", mmff("nitrobenzene", "O=[N+]([O-])c1ccccc1"),
     NEUTRAL_SETTINGS, "EAS on nitrobenzene: meta directing, strongly deactivated",
     ["stock1963", "sjoberg1990", "murray1992"], dict(type="lower", a=["meta"], b=["ortho", "para"]),
     dict(ortho=[4, 8], meta=[5, 7], para=[6]),
     "Sjoberg 1990 Table 1 and Murray 1992 Fig. 1: minima over meta (13.03, 13.04 eV)")
case("anisole", "COc1ccccc1", mmff("anisole", "COc1ccccc1"), NEUTRAL_SETTINGS,
     "EAS on anisole: ortho/para directing, activated",
     ["stock1963"], dict(type="lower", a=["ortho", "para"], b=["meta"]),
     dict(ortho=[3, 7], meta=[4, 6], para=[5]),
     "not in Sjoberg 1990; the experimental o/p preference is the expectation")
case("chlorobenzene", "Clc1ccccc1", mmff("chlorobenzene", "Clc1ccccc1"), NEUTRAL_SETTINGS,
     "EAS on chlorobenzene: ortho/para directing, deactivated",
     ["stock1963"], dict(type="lower", a=["ortho", "para"], b=["meta"]),
     dict(ortho=[2, 6], meta=[3, 5], para=[4]),
     "not in Sjoberg 1990 (fluorobenzene is); the experimental o/p preference is the expectation")
# --- conditions decide ---------------------------------------------------------
case("phenoxide", "[O-]c1ccccc1", mmff("phenoxide", "[O-]c1ccccc1",
                                         " (an anion: MMFF is a known weak point, see CLAUDE.md)"),
     ANION_SETTINGS,
     "allyl/benzyl halides: O only in DMF, THF, MeOH, EtOH, glymes, dioxane and toluene; "
     "in water 49:41 (allyl chloride) and 65:24 (benzyl chloride) O:C, in phenol 22:78; "
     "C-alkylation is attributed to hydrogen bonding at O",
     ["kornblum1963"], dict(type="report_only"),
     dict(O=[0], ortho=[2, 6], para=[4]),
     "the solvent decides through specific hydrogen bonds, which a PCM continuum does not "
     "represent; no ground-state number should be read as the answer")

(OUT / "manifest.json").write_text(json.dumps(
    {"about": __doc__.split("\n\n")[0], "cases": CASES}, indent=1) + "\n", encoding="utf-8")
print(len(CASES), "calibration cases")
