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
    "mayr2011": dict(doi="10.1002/anie.201007100", evidence="abstract",
                     ref="Mayr, Breugst & Ofial (2011) Farewell to the HSAB Treatment of "
                         "Ambident Reactivity. Angew. Chem. Int. Ed. 50, 6470-6505"),
    "loos2003": dict(doi="10.1021/ja037317u", evidence="abstract",
                     ref="Loos, Kobayashi & Mayr (2003) Ambident Reactivity of the Thiocyanate "
                         "Anion Revisited. J. Am. Chem. Soc. 125, 14126-14132"),
    "kornblum1955": dict(doi="10.1021/ja01628a064", evidence="abstract",
                         ref="Kornblum, Smiley, Blackwood & Iffland (1955) The Mechanism of the "
                             "Reaction of Silver Nitrite with Alkyl Halides. J. Am. Chem. Soc. 77, 6269"),
    "tishkov2005": dict(doi="10.1002/anie.200501274", evidence="abstract",
                        ref="Tishkov, Schmidhammer, Roth, Riedle & Mayr (2005) Ambident Reactivity "
                            "of the Nitrite Ion Revisited. Angew. Chem. Int. Ed. 44, 4623-4626"),
    "kornblum1963": dict(doi="10.1021/ja00891a024", evidence="abstract",
                         ref="Kornblum, Berrigan & le Noble (1963) Solvation as a Factor in the "
                             "Alkylation of Ambident Anions. J. Am. Chem. Soc. 85, 1141-1147"),
    "stock1963": dict(doi="10.1016/S0065-3160(08)60277-4", evidence="metadata",
                      ref="Stock & Brown (1963) A Quantitative Treatment of Directive Effects in "
                          "Aromatic Substitution. Adv. Phys. Org. Chem. 1, 35-154"),
    "sjoberg1990": dict(doi="10.1139/v90-220", evidence="metadata",
                        ref="Sjoberg, Murray, Brinck & Politzer (1990) Average local ionization "
                            "energies on the molecular surfaces of aromatic systems as guides to "
                            "chemical reactivity. Can. J. Chem. 68, 1440-1443"),
    "murray1992": dict(doi="10.1016/0166-1280(92)85015-D", evidence="metadata",
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


# --- ambident anions (the calibration pilot's three) --------------------------
case("enolate", "C=C[O-]", fixture("enolate"), ANION_SETTINGS,
     "SN2 alkylation with alkyl halides goes to C; hard electrophiles (silyl, acyl) to O",
     ["mayr2011"], dict(type="lower", a=["C_CH2"], b=["O"]),
     dict(C_CH2=[0], C_CO=[1], O=[2]),
     "ESP points to O; the pilot found ALIE min on the CH2 carbon (5.71 vs 7.47 eV)")
case("thiocyanate", "[S-]C#N", fixture("thiocyanate"), ANION_SETTINGS,
     "SN2 alkylation gives R-SCN (S); N attack only for highly stabilised carbocations",
     ["loos2003", "mayr2011"], dict(type="lower", a=["S"], b=["N"]),
     dict(S=[0], C=[1], N=[2]),
     "ESP points to N; the pilot found ALIE min on S (7.10 vs 8.99 eV)")
case("nitrite", "N(=O)[O-]", fixture("nitrite"), ANION_SETTINGS,
     "SN2 alkylation gives mainly nitroalkane (N); SN1/Ag+ conditions give alkyl nitrite (O)",
     ["kornblum1955", "tishkov2005", "mayr2011"],
     dict(type="not_clearly_lower", a=["N"], b=["O"], tol_eV=0.05),
     dict(N=[0], O=[1, 2]),
     "the pilot found O 7.86 vs N 7.89 eV: a tie. ALIE must not be claimed to pick N")
# --- electrophilic aromatic substitution (classic surface-ALIE use) ------------
case("anisole", "COc1ccccc1", mmff("anisole", "COc1ccccc1"), NEUTRAL_SETTINGS,
     "EAS on anisole: ortho/para directing, activated",
     ["stock1963", "sjoberg1990"], dict(type="lower", a=["ortho", "para"], b=["meta"]),
     dict(ortho=[3, 7], meta=[4, 6], para=[5]),
     "ortho/para partial rate factors far exceed meta")
case("chlorobenzene", "Clc1ccccc1", mmff("chlorobenzene", "Clc1ccccc1"), NEUTRAL_SETTINGS,
     "EAS on chlorobenzene: ortho/para directing, deactivated",
     ["stock1963", "sjoberg1990"], dict(type="lower", a=["ortho", "para"], b=["meta"]),
     dict(ortho=[2, 6], meta=[3, 5], para=[4]),
     "a halogen deactivates but still directs ortho/para")
case("nitrobenzene", "O=[N+]([O-])c1ccccc1", mmff("nitrobenzene", "O=[N+]([O-])c1ccccc1"),
     NEUTRAL_SETTINGS, "EAS on nitrobenzene: meta directing, strongly deactivated",
     ["stock1963", "sjoberg1990", "murray1992"], dict(type="lower", a=["meta"], b=["ortho", "para"]),
     dict(ortho=[4, 8], meta=[5, 7], para=[6]),
     "meta is the least deactivated position; a hard case for any one-number descriptor")
# --- conditions decide ---------------------------------------------------------
case("phenoxide", "[O-]c1ccccc1", mmff("phenoxide", "[O-]c1ccccc1",
                                         " (an anion: MMFF is a known weak point, see CLAUDE.md)"),
     ANION_SETTINGS, "phenoxide alkylation: O in aprotic solvents, C(ortho/para) in protic ones",
     ["kornblum1963"], dict(type="report_only"),
     dict(O=[0], ortho=[2, 6], para=[4]),
     "the solvent decides, so no single ground-state number should be read as the answer")

(OUT / "manifest.json").write_text(json.dumps(
    {"about": __doc__.split("\n\n")[0], "cases": CASES}, indent=1) + "\n", encoding="utf-8")
print(len(CASES), "calibration cases")
