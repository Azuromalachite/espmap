# Adversarial review: the espmap generalisation roadmap

Reviewer 3, 2026-10-22. Subject: the 10-item roadmap (comparison mode → Colab → annotations → JSON/CSV → ensembles → protonation → browser UI → fast backend → more properties → more DFT). Checked against `esp_map.py` (1,209 lines), `README.md`, `CLAUDE.md`, and `docs/adversarial-review-1.md` and `-2.md`.

Labels: `[V]` verified in the repo · `[M]` measured or computed by me · `[I]` inferred · `[U]` not verifiable here.

## The one thing to do first

**Time one real target molecule before committing to any ordering.** The roadmap is built around C18 pheromones, BDSF and decanoic acid, but every timing in the repo is from a 6-ring aromatic with 141–298 basis functions. The pheromone acetate has about **900**. Nothing in the roadmap prices that, and the answer reorders the list (F1).

## Findings

### F1 · The compute cost of the target molecules is never priced · HIGH

**The data.**
- `[V]` Measured SCF times, one core, B3LYP:
  - bromobenzene, def2-SVP (141 nao): 20 s (README)
  - chlorobenzene, def2-TZVP (253 nao): 165 s (review 1)
  - chlorobenzene, def2-TZVPP (298 nao): 224 s (review 1)
- `[V]` `run_scf` uses no density fitting. `grids.level = 3`. The default basis is def2-TZVP.
- `[M]` Basis-function counts for def2-TZVP (H = 6, C/O = 31, Cl = 37). The formula reproduces chlorobenzene's 253 exactly.

**Projected SCF time per conformer.** `[M]` The scaling exponent fitted from the 141 → 253 timings is about 3.6, from two points on aromatics only, so the times below are order-of-magnitude `[I]`.

| molecule | formula | nao | ≈ SCF, one core |
| --- | --- | --- | --- |
| chlorobenzene (measured) | C₆H₅Cl | 253 | 165 s |
| decanoic acid | C₁₀H₂₀O₂ | 492 | ~30 min |
| BDSF | C₁₂H₂₂O₂ | 566 | ~50 min |
| Z3,Z13-18:OH | C₁₈H₃₄O | 793 | ~3 h |
| Z3,Z13-18:OAc | C₂₀H₃₆O₂ | 898 | ~4–5 h |

**Why it matters.** `[I]` Item 5 (`--ensemble 10`) on the three pheromone variants comes to roughly 100 h of SCF. Items 2 and 7 (Colab and a browser UI) promise "click Run, see the ESP" for exactly these molecules. Item 8 (a fast mode) is ranked 8th but is a **prerequisite** for items 2, 5 and 7 on the stated targets.

**Correction.** Move a labelled preview mode to just after item 1. The cheapest version is already in the repo: def2-SVP plus `mf.density_fit()`. `SMALL_BASIS` already prints a "do not quote a number" warning for this case. See F8 on why xTB is not a drop-in replacement.

**Validation.** Time decanoic acid and the acetate at SVP and TZVP, each with and without density fitting. Report the σ-hole or V_min shift caused by density fitting. `[U]` This could not be run here: Python is not installed on this machine `[M]`.

### F2 · "Determine an appropriate common scale" collides with invariant 7 · HIGH

**The conflict.**
- `[V]` Invariant 7 says diverging and symmetric scales are used **only when the field changes sign**. An ion's surface never changes sign, and centring on zero collapses its range into one colour.
- `[I]` Item 6 proposes exactly the series that breaks this: acid ↔ carboxylate. The repo ships `examples/acetate.xyz` `[V]`.
  - A neutral member needs a diverging scale. The anion needs a sequential one.
  - A shared symmetric scale wide enough for the anion washes out the neutral.
  - A shared sequential scale discards the sign of the neutral molecule's field, which is the chemistry the reader needs.

The roadmap leaves the rule for this case to be decided later. That is the "threshold deferred" defect class.

**Correction.** Specify the rule for mixed-sign series before writing `compare`. One defensible rule:
- one shared diverging scale for every member whose field changes sign;
- the all-one-sign members in a separate, labelled sub-panel with their own shared sequential scale;
- never one bar spanning both.

**Validation.** Add a test where acetic acid plus acetate must produce two colour bars, with the orientation assertions from `--selftest` applied to each.

### F3 · The ensemble case rests on glycerol and inverts the repo's own finding · MEDIUM

**Where the evidence comes from.**
- `[V]` The "conformer ~1 kcal/mol" row in the error budget was measured on **glycerol** (README l.136): three carbons, across 8 kcal/mol of MMFF energy.
- `[I]` Carrying that over to an 18-carbon chain with about 15 rotatable bonds is extrapolation. CLAUDE.md itself says "validating on one molecule licenses nothing about another."

**Two specific problems.**
1. **Sampling is too thin, and "lowest energy" is likely misleading.** `[V]` The default is `--nconfs 20`, and the run picks the lowest-energy MMFF conformer. `[I]` Twenty conformers badly undersample a C18 chain. The gas-phase lowest-MMFF conformer of a long chain also tends to be folded, so it is not representative of the extended chain.
2. **The proposed headline reports the robust quantity.** The roadmap's example is "V_max = 42.1 ± 0.7 across 10 conformers." The repo's measured lesson is the opposite: **extrema are robust, shape is not**. The ± on extrema will mostly be small. The values that should vary are the shape descriptors: dipole, positive and negative surface area, and the variance of V over the surface.

**Correction.**
- Report shape descriptors across the ensemble, with the extrema as a secondary output.
- Raise `nconfs` for flexible molecules, scaled by rotatable-bond count.
- Warn when the selected conformer is compact (radius of gyration well below the extended chain).

**Validation.** On one pheromone, compare the spread of V_min across conformers with the spread of dipole and surface areas. If the extrema spread exceeds 1 kcal/mol, the glycerol number does not transfer and the error budget needs a new row.

### F4 · `--svg` cannot mean what the roadmap implies · MEDIUM

**The renderer is raster.** `[V]` `render` is the project's own z-buffer rasteriser, and the image is written with `plt.imsave` (l.1203). A molecule rendered this way can only appear in an SVG as an embedded bitmap.

**What the flags can do.** `[I]` SVG or vector output is real only for the colour bar and the annotations. `--transparent` needs alpha from the z-buffer's empty pixels, which is cheap. `--dpi` here means `--size`, which already exists `[V]`.

**Correction.** Rename the item to "transparent PNG + separate vector colour bar and labels."

### F5 · The JSON report is ranked 4th but item 1 depends on it · MEDIUM

**The dependency.** Item 1 itself promises "write a machine-readable results table." `[V]` At present the only numeric output is two `print` lines (l.1175, 1178).

**Cost.** `[I]` Most fields are already in memory or close to it:
- the mesh area comes straight from the marching-cubes faces;
- `mf.dip_moment()` is a single call;
- formula and charge are known from the gate.

**Correction.** Build the report first, as the data layer that `compare` writes.

### F6 · Alignment is ill-defined for the flexible series it targets · LOW–MEDIUM

**The problem.** `[I]` "Scaffold alignment" works for rigid homologues. For E/Z isomers of a C18 chain, the chain *is* the scaffold, and its 3D path differs between conformers, so a whole-molecule fit is arbitrary.

**Correction.** Align on the maximum common substructure of the **head group** (OH, OAc or CHO, plus the two or three nearest carbons). Fix the camera relative to that.

### F7 · Batch mode must keep "gate before compute" for the whole set · LOW

**The rule.** `[V]` Invariants 1–3 require one reader and one embedder, and the structure the gate approved must be the structure that is computed.

**What the roadmap misses.** `[I]` In a batch, gate **all N** molecules before starting any SCF. Otherwise a malformed row 7 surfaces after hours of compute on rows 1–6. The roadmap lists "validate every structure" first but does not say this ordering is required.

**Correction.** Make `compare` run the gate over all inputs, print one combined report, and start computing only if every row passes or `--force` is given.

### F8 · An xTB "fast mode" does not produce the same quantity · LOW

**The difference.** `[I]` GFN-xTB is a minimal-basis tight-binding method. An ESP taken from its density, or from its atomic charges, is a different approximation, not a cheaper version of the same calculation. It would need its own error row before even a "preview" label is honest. It would also add a dependency outside PySCF.

**Correction.** Make def2-SVP plus density fitting the preview mode (F1). Treat xTB as a separate, later backend.

## What it got right

These are the hard calls, checked:

- **"More DFT methods" ranked last.** `[V]` This matches the measurements: the functional spread is 0.69 kcal/mol, against +3.55 for basis and −3.61 for solvent.
- **"Comparable" is the right framing.** `[V]` The code already warns that auto-scaled panels "are NOT comparable" (l.1188–1191). The roadmap turns that warning into a feature.
- **Protonation stays advisory.** No automatic pKa prediction is proposed, and `[M]` the repo has no ionisable-group detection today. The proposal is honest about its own uncertainty.
- **The gate stays visible in any UI.** This correctly identifies the gate as the repo's strongest feature, and E/Z checking really is implemented (`UNSPECIFIED_DOUBLE_BOND`, l.383).
- **3D surfaces are not averaged across conformers.** This is the right instinct and avoids a conceptually incoherent object.
- **The factual claims about the repo all check out:**
  - configurable basis and functional, PCM with hexane among the solvents (l.152–154);
  - the gate's checks;
  - charge-aware scaling;
  - PNG plus stdout as the only outputs;
  - the dependency list matches `pyproject.toml`.

## Summary

| # | finding | severity | change to roadmap |
| --- | --- | --- | --- |
| F1 | 900-nao targets, cost unpriced | HIGH | measure first; preview mode moves up to #2 |
| F2 | shared scale vs invariant 7 on acid/anion series | HIGH | specify the mixed-sign rule before `compare` |
| F3 | ensemble evidence from glycerol; ± on the robust quantity | MEDIUM | report shape descriptors; scale `nconfs` |
| F4 | SVG of a raster render | MEDIUM | vector colour bar and labels only |
| F5 | JSON ranked below its dependant | MEDIUM | build it as part of #1 |
| F6 | scaffold alignment of flexible chains | LOW–MED | align on the head-group MCS |
| F7 | batch gating order | LOW | gate all rows, then compute |
| F8 | xTB ≠ cheap DFT | LOW | SVP + density fitting as preview |

**Revised order:**
1. timing measurement
2. report layer
3. comparison mode, with the F2 rule and the F7 gating
4. preview mode
5. output polish
6. ensembles (shape descriptors)
7. protonation advisory
8. Colab
9. browser UI
10. other properties
11. more functionals

## Sources

All sources are in-repo:
- `esp_map.py`: l.148–154 (solvents), 197–215 (embedding), 356–397 (stereo and tautomer gate), 572–600 (`run_scf`), 1175–1203 (outputs)
- `README.md`: l.39–47 (timings), l.136 (glycerol conformer row)
- `docs/adversarial-review-1.md`: l.17–21 (basis timings)
- `CLAUDE.md`: invariants 1–3 and 7, "Deliberately absent"

No external literature was retrieved. The xTB point (F8) is general knowledge, marked `[I]`.
