# Adversarial review — `esp_map.py` and the claims in `esp-colormaps.md`

Reviewer 3, 2026-10-22. Reviewing my own code and my own document from yesterday.
Environment: Python 3.13.16, numpy 2.5.3, pyscf 2.14.0, scikit-image 0.26.0, matplotlib 3.11.2, colorspacious, cmcrameri.
Revised artifact: `esp_map.py` (v2). Every finding below was reproduced by execution, not inferred.

---

## RETRACTION — read this first

**`esp-colormaps.md` §2.2 claims:**

> Chlorobenzene's σ-hole is **+6.4** kcal/mol. The apparent neutral of the RGB map is **+9.6**. So the entire σ-hole of chlorobenzene sits *below the perceived zero* and renders as neutral-to-electron-rich — the electrophilic site that makes halogen bonding work reads as its opposite.

**That claim is not supported and must be withdrawn.** It was built on a def2-SVP number with a 3.2 kcal/mol margin and no uncertainty attached. Basis-set convergence moves the σ-hole by **+3.55 kcal/mol** — larger than the margin — and moves it in the direction that destroys the claim:

| basis | nao | σ-hole V<sub>max</sub> | SCF time |
| --- | --- | --- | --- |
| def2-SVP | 127 | **+6.35** | 20 s |
| def2-TZVP | 253 | **+9.83** | 165 s |
| def2-TZVPP | 298 | **+9.90** | 224 s |

At convergence the σ-hole is **+9.90** and the map's apparent neutral is **+9.6**. The σ-hole does not fall below the perceived zero; it lands on it, to within 0.3 kcal/mol.

**What survives, stated at its real strength:** the RGB map renders chlorobenzene's σ-hole in very nearly its *least saturated* colour — it paints the molecule's electrophilic site the colour the map uses for "nothing here". That is still a real indictment, and it is now robust (0.3 kcal/mol coincidence against a 0.2 kcal/mol combined non-basis uncertainty). But it is a **coincidence between one molecule and one colormap**, not a general result, and I presented it as the decisive general case. The generalisable finding is §2.1 (chroma at zero), which needed no chemistry at all.

I also gave this claim the most rhetorical weight in the summary — "the decisive example". That emphasis was wrong twice over: the claim was unsupported, and the finding that *is* solid was filed beneath it.

---

## Findings

### F1 · The σ-hole number was quoted without an uncertainty `[BLOCKER]` · CONFIRMED EXECUTION FAILURE

The retraction above is the instance; this is the class. v1 quoted V<sub>max</sub> to 0.1 kcal/mol from a double-zeta basis, gas phase, idealised geometry, and compared it against a 3.2 kcal/mol margin. No convergence test of any kind was run before the number was published.

Error budget, now measured:

| source | effect on chlorobenzene σ-hole V<sub>max</sub> |
| --- | --- |
| **basis set** (def2-SVP → def2-TZVPP) | **+3.55** |
| grid spacing (0.26 → 0.16 bohr) | +0.04 |
| σ-hole cone definition (9 combinations) | 0.00 |
| C–Cl bond length, ±0.01 Å (dV/dR ≈ −14 kcal·mol⁻¹·Å⁻¹) | ±0.15 |
| functional | **not tested** |
| solvent | **not tested** (gas phase only) |

**The basis error is molecule-specific, so it cannot be calibrated away.** Re-running bromobenzene at the new default: σ-hole **+16.05** (def2-TZVP) against **+16.81** (def2-SVP) — a shift of **−0.76**, opposite in sign and a fifth the magnitude of chlorobenzene's **+3.55**. A reviewer who had checked only bromobenzene would have concluded the double-zeta basis was fine. The global range also moves (−11.7 to +24.5 at TZVP, against −12.9 to +21.1 at SVP), so the published figure's ±21.2 scale is itself a def2-SVP artifact — harmless for a colormap comparison, since all six panels share it, but it is not the molecule's converged range.

**Correction.** Default basis raised to def2-TZVP; `SMALL_BASIS` warning names the measured +3.5 kcal/mol penalty on stderr; module docstring states that V<sub>max</sub> carries at minimum a basis-set uncertainty and must not be quoted against a smaller margin.
**Validation.** `run_scf('He 0 0 0', basis='def2-svp')` warns, `def2-tzvp` is silent (T9, confirmed).

---

### F2 · "Apparent neutral at +9.6 kcal/mol" is not a property of the colormap `[HIGH]` · CONFIRMED EXECUTION FAILURE

The chroma minimum sits at a *fractional position along the ramp*. Converting it to kcal/mol multiplies by the chosen range, so the headline number is a property of the colormap **times my choice of ±21.2**:

| colormap | t at C′<sub>min</sub> | V at LIM=10 | V at LIM=21.2 | V at LIM=50 |
| --- | --- | --- | --- | --- |
| Rainbow | 0.637 | +2.7 | +5.8 | +13.7 |
| Red–green–blue | 0.727 | +4.5 | **+9.6** | +22.7 |
| RdBu | 0.488 | −0.2 | −0.5 | −1.2 |
| vik | 0.504 | +0.1 | +0.2 | +0.4 |
| PuOr | 0.492 | −0.2 | −0.3 | −0.8 |
| PINDR div | 0.524 | +0.5 | +1.0 | +2.4 |

The document quotes +9.6 kcal/mol with no range stated. Against a molecule with a wider range the same colormap "displaces the neutral" by twice as much.

**Correction.** The range-independent statement is `t = 0.727` — the perceptual neutral sits 73 % along the ramp instead of at the midpoint. Any physical-units version must carry the range.
**Validation.** Recompute at three ranges and confirm `t` is invariant while V scales linearly (confirmed, table above).

**A second, worse problem with the same number.** "Chroma minimum" is a *proxy* for "where a reader perceives neutral", and I never validated the proxy — no reader was asked, and no psychophysical source was cited. The defensible form of §2.1 does not need it: at the true midpoint `t = 0.5`, chroma is **33.7** (rainbow) and **39.4** (RGB) against **1.1–3.3** for the four proper maps. That is a direct measurement of the colormap at the data's zero, with no inference about perception. It should have been the headline.

---

### F3 · The validation check presented as passing cannot fail usefully `[HIGH]` · CONFIRMED EXECUTION FAILURE

`esp-colormaps.md` lists "Density integral vs electron count — 76.93 / 76, ±1 %" in its validation table. That check is insensitive to the quantity it was offered as evidence for.

Same wavefunction, varying only grid spacing (chlorobenzene, def2-TZVP, 58 electrons):

| spacing / bohr | ∫ρ dV | σ-hole V<sub>max</sub> |
| --- | --- | --- |
| 0.34 | 55.06 (−5.1 %) | +9.39 |
| 0.30 | 57.12 (−1.5 %) | +9.72 |
| 0.26 | 59.47 (+2.5 %) | +9.83 |
| 0.20 | **61.97 (+6.8 %)** | +9.84 |
| 0.16 | 57.06 (−1.6 %) | +9.87 |

The integral swings over a 12-point range **non-monotonically** while the quantity of interest moves 0.5 %. It is sampling the nuclear cusp, where the grid is far too coarse, and tells you nothing about the valence density at ρ = 0.002 where the surface lives. Worse: it is **59.46 at def2-SVP and 59.47 at def2-TZVP** — identical to two decimals across the basis change that moves V<sub>max</sub> by 3.5 kcal/mol. The "±1 %" in the document was an artifact of which spacing I happened to pick.

**Correction.** Relabelled `[diagnostic, not a check]` in the output and in the docstring; it gates nothing.
**Validation.** The real check is now `--selftest` (F7), which exercises the ESP itself.

---

### F4 · No guard that the isosurface is closed and interior to the grid box `[HIGH]` · CONFIRMED EXECUTION FAILURE

`margin` was hard-coded at 5.0 bohr, unexposed and unchecked. If the ρ = 0.002 surface reaches a box face, marching cubes returns an **open** mesh — the box face becomes a flat sheet painted with ESP values, indistinguishable from real structure — and every number in the run stays plausible:

| margin / bohr | verts | ρ on box faces | χ | σ-hole V<sub>max</sub> | global V<sub>max</sub> |
| --- | --- | --- | --- | --- | --- |
| 5.0 | 9,984 | 9.7e−05 | **2** | **+6.35** | +20.80 |
| 3.5 | 9,182 | 4.7e−03 | 0 | **+4.69** (−26 %) | +20.81 |
| 2.5 | 5,529 | 2.6e−02 | 1 | nan | +20.77 |
| 2.0 | 3,670 | 8.0e−02 | −4 | nan | +20.76 |

At margin 3.5 the σ-hole is wrong by 26 % and nothing reports it. At tighter margins the cone yields `nan` while `global V_max` keeps printing a healthy-looking +20.8.

**Correction.** `extract_surface()` now computes two independent detectors — maximum density on the six box faces, and the Euler characteristic (χ = 2 for a single closed genus-0 surface) — and raises by default. `--margin` exposed; `--no-strict` is the documented escape hatch.
**Validation.** T3: `--margin 3.5` raises with **both** detectors firing. T4: `--margin 5.0` succeeds with χ=2 and reproduces +6.35. T7: `--no-strict` warns loudly and reports the damaged +4.69 so the user sees it.

---

### F5 · `int1e_grids` allocation unguarded, quadratic in basis size `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT

v1 used a fixed `chunk=1500`, allocating `(chunk, nao, nao)` float64:

| nao | v1 allocation | v2 chunk | v2 allocation |
| --- | --- | --- | --- |
| 141 (bromobenzene/SVP) | 0.24 GB | 2000 | 0.32 GB |
| 298 (chlorobenzene/TZVPP) | 1.07 GB | 703 | 0.50 GB |
| 900 (a 40-atom drug/TZVP) | **9.72 GB** | 77 | 0.50 GB |
| 2000 | **48.0 GB** | 16 | 0.51 GB |

**Correction.** `chunk = clamp(mem_budget_gb·1e9 / (8·nao²), 16, 2000)`.
**Validation.** T8, confirmed. Honest caveat: the floor of 16 means the budget is exceeded above nao ≈ 1976 (0.51 GB at nao=2000). The floor is deliberate — a chunk of 0 would hang — and the overshoot is 2 %.

---

### F6 · The CLI default contradicted the document's own rule `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT

`esp-colormaps.md` §3 item 3 says to fix the range across a figure. v1's `--lim` defaulted to **per-molecule auto-scaling**, silently, which is precisely what that rule forbids. Compounding it, `lim = a.lim or max(...)` makes `--lim 0` falsy, so an explicit zero fell through to auto-scaling with no message.

**Correction.** Auto-scaling now prints a stderr warning naming the consequence; `--lim 0` raises; out-of-range vertices are reported as a clamped percentage when `--lim` is given.
**Validation.** T5 (`--lim 0` rejected), T6 (omitted `--lim` warns and still renders), confirmed.

---

### F7 · No runtime check of the thing most likely to be wrong `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT

The v1 sign error — `int1e_grids` returns the positive Coulomb kernel, so the electronic term must be subtracted — was caught by hand, by noticing the far field failed to decay. Nothing in the shipped code would have caught it. A sign-inverted ESP map is perfectly plausible-looking: smooth, structured, publishable, and wrong everywhere.

**Correction.** `--selftest` runs water and asserts three probes against published windows, plus the sign asymmetry (lone pair negative *and* donor H positive).
**Validation.** T1 passes on healthy code. T2 re-introduces the v1 `+=` and all four assertions fail, with the far field reading **+265.87** instead of ~0 — the same signature that exposed it originally. Fault injection confirms the guard catches the real historical bug, not a strawman.

---

### F8 · Global normal flip keyed on numerical noise `[LOW]` · CONFIRMED STATIC CODE DEFECT, measured immaterial

`if np.mean(vn_s[:,2]) > 0: vn_s = -vn_s`. On a closed surface the mean normal is ~0, so the test statistic is 3–4 orders of magnitude below the per-vertex scale (|mean|/RMS = 1.8e−03 at elev 20/azim 172, 6.9e−04 at 12/−62) and its *sign* flips between viewpoints.

**Measured before recommending:** it changed nothing. Global-flip and per-vertex orientation agree on **100.0 %** of visible vertices at both viewpoints, and the lit fraction is identical (34.1 % and 27.6 %). The published figures are unaffected. It is worth removing because the open-mesh case from F4 is where the mean normal stops being ~0 and the flip starts doing something real.

**Correction.** Per-vertex orientation toward the camera.
**Validation.** Agreement measured at 100 % on closed meshes (no visual regression); the guard in F4 now prevents the open-mesh case from arising silently.

---

### F9 · A cause asserted without diagnosis `[LOW]` · METHODOLOGICAL CONCERN

`esp-colormaps.md` states: *"Iodobenzene currently fails the σ-hole extraction because def2-SVP uses an ECP for iodine and the core density changes where the isosurface falls near the halogen."* I never tested that. The observed failure was `V[cone].max()` on an empty selection; the ECP explanation is a plausible story I wrote down as if it were a finding.

**Correction.** Reworded to state the symptom and that the cause is not diagnosed. `sigma_hole()` now returns `(nan, 0)` with an explicit vertex count rather than raising on an empty cone, so the next person sees *which* step failed.

---

## Concerns that tested negative — report these too

Three things I flagged as likely defects and could not substantiate. Recording them so they are not re-raised:

- **σ-hole cone thresholds.** V<sub>max</sub> was **identical (+9.84)** across `pr_min ∈ {0, 0.3, 0.6} × pe_max ∈ {1.0, 1.6, 2.2}` — nine combinations, 76 to 576 vertices. The maximum is interior to even the tightest cone, so the hand-tuned thresholds do not bias it. Not a defect. (The function still reports its vertex count, because a *different* molecule could put the max on the boundary.)
- **Grid spacing.** V<sub>max</sub> converges monotonically to within 0.04 kcal/mol by 0.26 bohr. The default was already adequate. Not a defect.
- **Idealised geometry.** dV/dR ≈ −14 kcal·mol⁻¹·Å⁻¹, and the built-in C–Cl 1.740 Å is within 0.001 Å of the microwave value, so geometry contributes ≈ ±0.15 kcal/mol. Materially irrelevant next to F1 — and the agreement is luck, not method, which the docstring now says.

---

## On the modelling skills

Asked to involve them "if necessary". Most of `mechanistic-model-construction` does **not** apply and I will not pretend otherwise: there is no dynamics, no population, nothing fitted, so dimensionless groups, timescale elimination, the CME→ODE ladder and closure-naming have no referent in a single-point electronic-structure calculation.

Three of its rules apply with full force, and all three were broken:

- *"Hand off to the identifiability skill before any parameter is quoted."* I quoted V<sub>max</sub> = +6.4 against a 3.2 margin with no uncertainty. That is F1, and it is the blocker.
- *"Check that the model you tested is the model you are making a claim about."* The claim was about chlorobenzene's σ-hole; the model tested was a double-zeta gas-phase calculation whose error on that specific quantity exceeds the claim's margin.
- *"A caveat that is not where the reader meets the model is not a caveat."* v1 did disclose the small basis and the idealised geometry — in a "Limits" section at the *bottom* of the document, three sections after the claim that depended on them. The disclosure was real and useless.

`mechanistic-fit-identifiability` has no purchase here: nothing is calibrated against data, so there is no parameter-vs-target arithmetic to run.

---

## What it got right

Naming the hard decisions, not the easy ones:

- **Evaluating the ESP at surface vertices instead of on a volumetric grid.** This is the correct architecture and it is not the conventional one. It is ~100× cheaper and removes the cube-file interpolation error entirely, because the potential is computed exactly where it is displayed. It survives review untouched.
- **Validating the physics before rendering anything.** The water probes caught the sign error. The instinct to check a chemically known answer first is why there was a correct figure at all — the failure was not checking *convergence* as well as *correctness*.
- **Rasterising once and recolouring.** Guarantees the six panels differ only in colormap. Without it, a reviewer could not distinguish a colormap effect from a geometry or lighting difference, and the central figure would prove nothing.
- **Suppressing specular highlights and weakening diffuse shading.** Correct and underappreciated: a white highlight at an arbitrary location on an ESP surface reads as a potential extremum.
- **Orienting every candidate colormap warm = negative.** Fixing perception without silently flipping the sign convention is the kind of restraint that makes a comparison legitimate.
- **Reporting the greyscale result that went against the thesis** (RdBu ΔE′ = 3, worse than the crude RGB map's 26) instead of dropping it.

---

## Verdict

**Blocks:** F1 is a blocker and has been acted on — the claim is retracted above, and the defaults, warnings and docstring now make the same mistake harder. F2's physical-units framing and F3's fake validation must be corrected in `esp-colormaps.md`; both are now patched. F4 is a blocker for anyone *else* running the tool, and is fixed and fault-injected.

**Can wait:** F5 matters only above ~300 basis functions. F8 is cosmetic on closed meshes and measured immaterial. F9 is a one-line honesty fix, already applied.

**Not addressed, and should be before any number from this is quoted:** functional dependence (B3LYP only) and solvent (gas phase only). Both were untested in v1 and remain untested in v2; the docstring now says so rather than implying a completeness that was never there.

**The structural lesson.** Every substantive finding here is the same defect wearing different clothes: *a check that measures something adjacent to what the conclusion needs.* The density integral checked core sampling and was offered as evidence for valence accuracy. The water probes checked correctness and were offered as evidence for convergence. The chroma minimum measured colorimetry and was offered as evidence about perception. In each case the check passed, honestly, and licensed nothing.
