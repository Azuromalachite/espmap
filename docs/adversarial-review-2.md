# Adversarial review 2 — `esp_map.py` v3 → v4

Reviewer 3, 2026-10-22. Second pass, on the code added since review 1: SMILES and structure-file input, the admissibility gate, and scale selection — v3 grew from 814 to 1078 lines, most of it written fast.

Environment: Python 3.13.16, numpy 2.5.3, pyscf 2.14.0, rdkit 2026.09.1, scikit-image 0.26.0.
Standards applied: the project's own `computational-chemistry-properties` and `molecular-surface-visualisation` skills.
Every finding reproduced by execution. Revised artifact: `esp_map.py` (v4).

---

## Lead finding

**`--spin` was accepted, passed to a *restricted* solver, and returned a converged, plausible, physically meaningless number.** v3 built `dft.RKS` whatever the spin. PySCF does not refuse:

```
O atom, def2-SVP, spin=2  ->  nelec = (5, 3), which RKS cannot represent
  RKS (what v3 did)   -74.998688 Ha
  UKS (correct)       -75.000759 Ha        difference +1.3 kcal/mol
```

The energy is not merely 1.3 kcal/mol off — it is a restricted wavefunction describing an open shell, so every property derived from it, the whole ESP surface included, is unfounded. Nothing warned.

This was compounded by the gate's own advice. The radical rejection told users to *"rerun with `--spin` and a UKS build"* — pointing at a flag that produced garbage and a build that did not exist. The gate and the solver disagreed about what the pipeline could do.

**Correction.** `run_scf` selects `dft.UKS` when `spin != 0`; a `total_dm()` helper returns α+β for either solver, since both the density grid and the ESP electronic term want the total density; `<S²>` is printed. The gate now rejects a radical only when no `--spin` is given, and the remedy names the correct flag.
**Validation.** Methyl radical, `--spin 1`: `<S²> = 0.754` against the exact doublet value 0.75, V = −11.82 to +12.54 kcal/mol. Closed-shell path byte-identical to v3.

---

## Findings

### G2 · The σ-hole axis was hard-coded to a benzene ring `[HIGH]` · CONFIRMED EXECUTION FAILURE

`sigma_hole()` took its axis as `X − mean(atom_coords()[:6])` — "the first six atoms are the ring". True only for the built-in halobenzenes. It was called for **any** structure with exactly one halogen, including everything arriving by SMILES or file.

| molecule | axis error | v3 reported | correct axis |
| --- | --- | --- | --- |
| `Clc1ccccc1` | 4.3° | +6.45 | +6.45 |
| `CCCl` | 30.6° | −1.08 | −0.65 |
| `ClCCCCCCCC` | 25.4° | −1.70 | −1.70 |
| `CC(=O)OCCCl` | **40.5°** | **−0.05** | **+1.84** |

The last row is a sign flip on the reported quantity — the tool said the σ-hole was slightly negative where it is positive.

**Correction.** The axis is the bond vector from the substituent's nearest heavy neighbour to the substituent, which is what a σ-hole is defined along; `i_bond` may be passed explicitly.
**Validation.** `CC(=O)OCCCl` now reports **+1.72** (same sign and magnitude as the hand-computed +1.84; the residual is grid and geometry). Chlorobenzene unchanged at +6.70.

### G3 · `.xyz` was advertised and completely broken `[HIGH]` · CONFIRMED EXECUTION FAILURE

An XYZ file carries coordinates and **no connectivity**. RDKit returns a bondless mol, which the fragment check reads as one fragment per atom. Every XYZ file was therefore rejected:

```
water.xyz -> bonds 0, fragments 3, SMILES 'O.[HH].[HH]'
[BLOCK] MULTIPLE_FRAGMENTS — 3 disconnected fragments: O, [HH], [HH]
        fix: draw one molecule. A salt or solvate has no single ESP surface...
```

100 % failure rate on a listed format, with a remedy about salts that points the user nowhere.

**Correction.** `rdDetermineBonds.DetermineBonds` followed by `SanitizeMol` (needed before `AddHs` can compute implicit valences).
**Validation.** The same file now reads as `O`, geometry "file 3D coordinates used as-is", V = −44.79 to +48.40 kcal/mol — the correct water range.

### G4 · Two file readers; the fix reached only one `[HIGH]` · CONFIRMED EXECUTION FAILURE

`parse_for_check` and `from_file` each had their own reader. The G3 fix applied to the gate's copy and not to the one that builds geometry, so for a moment the gate **accepted** an XYZ file that the geometry step then crashed on — a worse failure than rejecting it.

This is the *one name, two definitions* class, and it is the reason G5 matters too.

**Correction.** A single `read_structure(path)` returning `(mol, problems, n_found)`; `parse_for_check` and `from_file` both call it.
**Validation.** One `MolsFromCDXMLFile` occurrence in the file (was two). All seven input paths exercised: `--smiles`, `.cdxml`, `.mol` 3D, `.mol` 2D, `.mol` 3D-without-H, `.xyz`, built-in.

### G5 · Two embedding implementations `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT

`_embed_and_minimise` was written for `from_file` while `from_smiles` kept its own copy of the same ETKDG + MMFF logic — two `EmbedMultipleConfs` call sites, free to drift in seed handling, force-field fallback or conformer selection.

**Correction.** One implementation, called by both.
**Validation.** One `EmbedMultipleConfs` occurrence (was two); `--smiles CCO` twice gives byte-identical PNGs, and `from_smiles` returns identical coordinates across calls.

### G6 · `--lim` on a sequential scale re-imposed symmetry `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

The sequential branch exists because an ion's surface never changes sign. But:

```
anion, data -162.0 .. -44.0
  --lim 200   ->  sequential scale  -200.0 .. +200.0
```

More than half the ramp is unreachable by any vertex, and the result is the symmetric range the branch was written to prevent. Directly contrary to `molecular-surface-visualisation` §1.

**Correction.** `--lim` is documented and enforced as the symmetric *diverging* range; on a sequential scale it is ignored with a note. `--vmin`/`--vmax` added for explicit sequential bounds, with `vmin < vmax` asserted.
**Validation.** Five combinations exercised; `--lim 200` on the anion now yields −161.9 to −45.0 with the note, `--vmin -170 --vmax -40` yields exactly that.

### G7 · The gate sealed a copy `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT

`parse_for_check` and `from_smiles`/`from_file` parse the input **independently**. The gate approved one object; the SCF computed another. They agreed in every case tested, but nothing asserted it — and for a multi-molecule file, two readers selecting differently is exactly the failure that would not announce itself.

**Correction.** The canonical SMILES of the computed structure is asserted against the approved one before any SCF.
**Validation.** Fault-injected by monkey-patching the geometry step to return a different molecule: refuses with `approved: Clc1ccccc1 / computing: CCO`. Healthy inputs pass.

**My first version of this correction was wrong and blocked valid input.** It compared the gate's *canonical* SMILES against the raw user string stored by `from_smiles`, so `C[CH2]` was rejected as "not the same structure as `[CH2]C]`". Both sides now go through `CanonSmiles`. The independence that matters is that the two *mols* were parsed separately, not that the comparison differs.

### G8 · The cost guard crashed on every open-shell molecule `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

`estimate_cost` built `gto.M` without `spin`, and PySCF refuses an electron count inconsistent with the spin:

```
RuntimeError: Electron number 17 and spin 0 are not consistent
```

So the guard that exists to prevent an expensive mistake was itself the first thing to fail once UKS made radicals reachable.

**Correction.** `estimate_cost(spec, basis, charge, spin)`.
**Validation.** Ethyl radical with `--spin 1` reports 53 basis functions and proceeds.

### G9 · Fluorine was reported as having a σ-hole `[LOW]` · METHODOLOGICAL CONCERN

`F` was in the halogen list, so fluorobenzene printed a "σ-hole V_max". Fluorine is too electronegative to form one; the printed number is noise on the F cap, and naming it invites a reader to treat it as a halogen-bond donor.

**Correction.** F excluded from σ-hole reporting, with an explicit note saying why.
**Validation.** `Fc1ccccc1` prints the note and no σ-hole; `Clc1ccccc1` unchanged.

---

## Compliance with the project's own skills

Checked v3 against the two standards it is meant to embody.

**`computational-chemistry-properties`** — complied on §4 (isosurface integrity: both detectors present and firing), §5 (geometry provenance stated), §6 (parse → check → embed, in that order), §2 (`--selftest` asserts the physics). **Violated §6's premise**: the skill qualifies the radical rejection with *"if the method is restricted/closed-shell"*, and v3's method was restricted while its CLI accepted `--spin` — the gate and the solver held different beliefs about the pipeline. G1 reconciles them.

**`molecular-surface-visualisation`** — complied on §1 (scale type chosen from the data), §7 (orientation asserted in code — 11 colormaps pass), §4 (specular off). **Violated §1 through the back door**: G6, where a user-supplied `--lim` silently reimposed the symmetric range on an ion.

Both violations are the same shape: a correct rule implemented in one place and bypassed by a parameter that reached the same decision by another route.

---

## What it got right

- **The admissibility gate is well-designed and its ordering is correct.** Parse → check → embed was the right call and it is what made G1 and G3 cheap to find: failures surface at the gate with a named code rather than as a traceback from inside a force field.
- **Each finding naming a remedy in the *source structure*** rather than in the code is the detail that makes the gate usable by a chemist rather than a programmer.
- **Detecting 2D-vs-3D from both the dimensionality flag and the coordinates**, rather than trusting either. A file can claim 3D and be flat, and this is the one check standing between a drawn structure and a planar SCF.
- **Freezing heavy atoms when adding hydrogens to a 3D file**, verified to reproduce the all-atom result exactly (+6.61 both ways). The obvious shortcut — a full force-field relax — would have silently altered a user's supplied geometry.
- **Measuring the stereochemistry claim before writing the warning.** The intuitive wording ("diastereomers will differ") was wrong; 0.9 kcal/mol against a 0.4 noise floor is what got written down.
- **The colormap-orientation assertions in `--selftest`.** They caught three inverted maps in review 1 and are the reason G6 was the only visualisation defect this round.

---

### G10 · A residual in my own G3 fix, found while verifying `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

`rdDetermineBonds.DetermineBonds(m, charge=0)` hard-coded neutrality, while the failure message told the user to *"set the charge explicitly"* — and `--charge` never reached the reader. A charged XYZ was unreadable and the remedy named a capability that did not exist.

This is **the same defect as G1**, committed in the fix for G3: an error message pointing at a flag the code ignores. Reproduced on an acetate XYZ — rejected identically with and without `--charge -1`.

**Correction.** `charge` threaded through `read_structure` → `parse_for_check` / `from_file`, and the message now names `--charge`, which works.
**Validation.** Acetate XYZ with `--charge -1` reads as `CC(=O)[O-]`, is accepted, and correctly routes to a sequential scale (−179.8 to −…). Without it, the message names a remedy that exists. Neutral XYZ unchanged.

---

## Verdict

**Blocks:** G1 (meaningless open-shell wavefunction), G2 (sign-flipped σ-hole on any non-benzene molecule), G3 and G4 (a listed input format broken, and the near-miss where the gate accepted what the geometry step could not read). All fixed and fault-injected.

**Can wait:** G5 and G7 are structural hygiene with no observed wrong output. G9 is a wording fix.

**Not addressed:** functional and solvent dependence remain untested, as in review 1. `--optimize` is still exercised only on one molecule. The gate cannot tell you the drawing is the molecule you *meant* to draw — a cleanly drawn wrong tautomer passes with a note.

**The structural lesson, and it is different from last time.** Review 1 was about checks that measured something adjacent to the conclusion. This round is about **the same rule implemented twice**: two file readers, two embedders, a gate parsing one object while the SCF computed another, and a solver choice that the CLI and the gate disagreed about. Every one of those let a correct fix land in one copy and not the other — G4 is the pure case, where the XYZ repair reached the gate and not the geometry path, turning a clean rejection into a crash. When a rule matters, it should exist once and be called twice.
