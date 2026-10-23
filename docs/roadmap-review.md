# Adversarial review: the espmap generalisation roadmap

Reviewer 3. Subject: the 10-item roadmap (comparison mode → Colab → annotations → JSON/CSV → ensembles → protonation → browser UI → fast backend → more properties → more DFT). Checked against `esp_map.py`, `README.md`, `CLAUDE.md`, and `docs/adversarial-review-1.md` and `-2.md`.

Labels: `[V]` verified in the repo · `[M]` measured or computed by me · `[I]` inferred · `[U]` not verifiable here.

**Revision 2 (2026-10-22).** The measurements this review asked for have been run (`bench/`). This revision replaces estimates with them, retracts four of my claims, adds F9 (fixed and verified), and adds a reactivity-descriptors direction.

**Status at v0.2.0 (2026-10-23).** Item 1 of the revised order is done for
what the record needed first. `--json` writes a versioned record,
`espmap.record/1`, whose key tree is `docs/record-schema.json`. It carries:
- area-weighted surface statistics (`results.surface_stats`) and GIPF
  statistics (`results.gipf`);
- the dipole and the computed geometry;
- a `validity` envelope;
- per-site statistics under two partition rules (`docs/per-site-spec.md`,
  stages P0 to P3).

Six fixture records are in `fixtures/records/`. F5's radius of gyration,
stereo descriptors and conformer ID are still missing. The rest of the
revised order stands; the reactivity descriptors (item 6) are the spec's
P4. The findings below are as of revision 2.

## Retractions — read these first

1. **F1's cost table was 10–20× too pessimistic.** I extrapolated one-core aromatic timings to the C18 targets, with a scaling exponent fitted on two points, and projected 4–5 h for the acetate. On the author's 16-thread laptop it took **875 s** `[M]`, and decanoic acid took 211 s, not ~30 min. The recommendation "measure first" stands. The conclusion I drew from the estimate does not: a fast mode is **not** a prerequisite for C18 work.
2. **F1 and F8 called def2-SVP + density fitting the "cheapest" preview.** At 898 basis functions, density fitting was **slower** than exact integrals on this machine (1,175 s vs 875 s) `[M]`. It is still faster up to 566 basis functions. Preview means SVP; density fitting is a separate, size-dependent choice.
3. **F3 said the lowest-energy MMFF conformer of a long chain "tends to be folded".** For the C18 acetate it is one of the most **extended**: 87 of 100 embedded conformers are more compact `[M]`. A seconds-long RDKit check would have caught this before I wrote it.
4. **An earlier draft of this revision said "the glycerol figure transfers to the C18 extrema".** That was measured only on the 8 conformers the χ = 2 check let through. The 2 it refused were the folded ones, and once recovered (F9) they widened the V_min spread from 1.1 to **5.2 kcal/mol** `[M]`. The check was removing exactly the cases that contradict the claim.

## The one thing to do first

**Extend the report layer (`--json`) with surface and shape descriptors, then build comparison mode on it.** Speed is no longer the blocker (F1). The measured conformer spread (F3) shows the numbers worth comparing between molecules are shape descriptors, and nothing records them yet.

## Status of the roadmap items

| item | state on branch `bench/timing` |
| --- | --- |
| fast / preview path | **done**: `--basis def2-svp`, optional `--density-fit`, optional `--gpu` (validated, F8) |
| Colab | **done**: `notebooks/espmap_colab.ipynb`, run on an L4 |
| JSON report | **partial**: `--json` has identity, settings, energy, dipole, V range, versions, caveat; no areas, Rg, stereo or conformer ID |
| 3D output | **done**: `--html` interactive view; `--npz` cached surface |
| comparison mode | not built; the notebook's series cell already enforces F2's rule |
| ensembles, protonation, alignment, browser UI, other properties | not built |

## Findings

### F1 · The compute cost of the target molecules was never priced · HIGH → resolved by measurement

**Measured** `[M]` on a Ryzen 7 260 (16 threads) in WSL2 with 7.4 GB RAM, PySCF 2.14.0, B3LYP (`bench/results.jsonl`, `bench/env.json`):

| molecule | nao (TZVP) | SVP exact / RI | TZVP exact | TZVP RI |
| --- | --- | --- | --- | --- |
| chlorobenzene | 253 | — | 39 s (150 s on 1 core; repo: 165 s) | 21 s |
| decanoic acid | 492 | 66 / 48 s | 211 s | 97 s |
| BDSF | 566 | 76 / 59 s | 286 s | 129 s |
| Z3,Z13-18:OH | 793 | — / 103 s | — | 813 s (machine shared) |
| Z3,Z13-18:OAc | 898 | 182 / 146 s | **875 s** | **1,175 s** idle; 1,613 s shared |

- **GPU** `[M]`: decanoic acid, TZVP with density fitting on a Colab L4, took 22.5 s, 4.3× faster than the 16-thread CPU (`bench/gpu_l4.json`).
- **Why density fitting slows down at 898 nao** `[I]`: peak memory was 2.9 GB, below PySCF's default `max_memory` of 4 GB. The integrals were probably sent out of core by that cap rather than by the machine's RAM. **Validation:** re-run with `mf.max_memory` raised to 6 GB.
- **Scratch disk** `[M]`: the first C18 density-fitting runs failed with `OSError 28`, because WSL's `/tmp` is a 3.7 GB tmpfs. `PYSCF_TMPDIR` must point at real disk.

**The basis step is large on the targets** `[M]`. Going from def2-SVP to def2-TZVP moves V_max by **+2.82** (decanoic acid), **+2.71** (BDSF), **+2.65** (C18 acetate) and **+1.53** (C18 alcohol), and V_min by +1.18, +1.01, +0.29 and +2.37. The alcohol, in the same file, was first left out; it is the molecule that shows the shift is molecule-specific (review 3 H8). Differences between related pheromones may be smaller than this. **Validation:** def2-TZVPP on one C18, which is cheap on a GPU, to show whether TZVP is converged.

### F2 · "Determine an appropriate common scale" collides with invariant 7 · HIGH

**The conflict.**
- `[V]` Invariant 7 allows a diverging, symmetric scale **only when the field changes sign**. An ion's surface never changes sign, and centring on zero collapses its range into one colour.
- `[I]` Item 6 proposes exactly the series that breaks this: acid ↔ carboxylate. Measured `[M]`: acetic acid spans −39.9 … +57.3 kcal/mol; acetate spans −171.8 … −79.9.

**Correction.** Specify the rule for mixed-sign series before writing `compare`:
- one shared diverging scale for the members whose field changes sign;
- the one-sign members in a separate, labelled panel with a shared sequential scale;
- never one colour bar spanning both.

**Status.** `[V]` The notebook's series cell now detects the mix and refuses it, naming the members to split off. Tested on acetic acid + acetate. `compare` must implement the two-panel version.

### F3 · The ensemble case rested on glycerol · MEDIUM → partly confirmed by measurement

**Measured** `[M]` (`bench/conformers.jsonl`, `bench/f9_recheck.jsonl`). C18 acetate: 100 conformers embedded; the lowest-energy one plus 9 spanning the range of radius of gyration were sent to B3LYP/def2-SVP with density fitting. The old surface check refused the 2 most compact; they were recomputed after the F9 fix.

| descriptor | 8 extended conformers | all 10 (adds the 2 folded) |
| --- | --- | --- |
| V_min | 1.10 kcal/mol | **5.20** |
| V_max | 1.18 kcal/mol | **2.48** |
| dipole | 0.77 D (43 % of mean) | 0.77 D |
| negative surface area | 68.8–78.5 Å² | **68.8–117.5 Å²** (60 % of mean) |
| total surface area | 413–424 Å² | 401–424 Å² |

- **Among extended conformers, the glycerol figure holds:** extrema move ~1.1 kcal/mol, while shape descriptors move far more.
- **Folding breaks it.** The folded conformer with genus 3 has V_min −35.70 against −39.8 to −40.9 for the rest, and twice the negative area. Intramolecular contact changes the local environment, so extrema are not conformer-robust for a chain that can fold. The error budget's ~1 kcal/mol conformer row must not be quoted for flexible molecules.
- **The folded conformers are 4.8 and 6.4 kcal/mol above the minimum in MMFF, gas phase** `[M]`. Whether they matter depends on the environment (solvent, receptor). That is an argument for reporting the spread, not for discarding them.
- **Retracted:** "lowest-MMFF conformer tends to be folded", and the earlier "glycerol figure transfers" (see the retractions above).

**Correction.**
- Report shape descriptors and extrema across the ensemble, with the spread, and flag folded members explicitly (genus > 0 is a cheap folding marker).
- Scale `nconfs` with the rotatable-bond count.
- Sample by shape (radius of gyration), not only by energy, or folded shapes are never seen.

### F4 · `--svg` cannot mean what the roadmap implies · MEDIUM

**The renderer is raster.** `[V]` `render` is the project's own z-buffer rasteriser, written out with `plt.imsave`. A molecule rendered this way can appear in an SVG only as an embedded bitmap.

**What the flags can do.** `[I]` Only the colour bar and labels can be real vector graphics. `--transparent` needs alpha from the z-buffer's empty pixels, which is cheap. `--dpi` here means `--size`, which already exists `[V]`.

**Note.** `[V]` The `--html` view now exists. It is interactive and shows a real colour bar, which covers most of what a vector figure was wanted for on screen.

**Correction.** Rename the item to "transparent PNG + separate vector colour bar and labels."

### F5 · The JSON report is ranked 4th but item 1 depends on it · MEDIUM

**Status.** `[V]` `--json` exists. It does not yet carry surface areas, radius of gyration, stereo descriptors or a conformer ID, which F3 shows are the fields worth comparing.

**Correction.** Extend `--json` first, as the data layer `compare` writes. `bench/conformers.py` already computes the areas and V statistics to copy.

### F6 · Alignment is ill-defined for the flexible series it targets · LOW–MEDIUM

**The problem.** `[I]` "Scaffold alignment" works for rigid homologues. For E/Z isomers of a C18 chain, the chain *is* the scaffold, and its 3D path differs between conformers, so a whole-molecule fit is arbitrary.

**Also seen.** `[V]` With the fixed default camera, decanoic acid rendered end-on down its chain in the Colab run.

**Correction.**
- Align on the maximum common substructure of the **head group** (OH, OAc or CHO, plus the two or three nearest carbons).
- For single molecules, default the camera to the principal axis with the head group on a fixed side.

### F7 · Batch mode must keep "gate before compute" for the whole set · LOW

**The rule.** `[V]` Invariants 1–3 require one reader and one embedder, and the structure the gate approved must be the structure that is computed.

**What the roadmap misses.** `[I]` In a batch, gate **all N** molecules before starting any SCF. Otherwise a malformed row 7 surfaces after the compute on rows 1–6.

**Correction.** Make `compare` run the gate over all inputs, print one combined report, and compute only if every row passes or `--force` is given.

### F8 · An xTB "fast mode" does not produce the same quantity · LOW → preview validated

**The difference.** `[I]` GFN-xTB is a minimal-basis tight-binding method. Its ESP is a different approximation, not a cheaper version of the same calculation, and it would add a dependency outside PySCF.

**Measured** `[M]`:
- **Density fitting vs exact:** 7 pairs, 253–898 nao, SVP and TZVP. |ΔV| ≤ 0.018 kcal/mol, with V_max shifted upward in every pair (+0.006 on the C18 acetate at TZVP).
- **Total energy shift** from density fitting: 6×10⁻⁵ to 3×10⁻⁴ Ha. So a validation gate must test V, not energy.
- **GPU vs CPU**, both density-fitted: |ΔV| < 0.001 kcal/mol (`bench/gpu_l4.json`).

**Correction.**
- The preview is def2-SVP; that is where the time is saved.
- Density fitting enters the error budget as its own row, at ≤ 0.02 kcal/mol, and is worth enabling up to ~600 nao on this hardware.
- Treat xTB as a separate, later backend with its own error row.

### F9 · The surface check refuses valid folded surfaces · MEDIUM (new)

**The problem.** `[V]` Invariant 8 required Euler characteristic χ = 2. `[M]` The two most compact C18 acetate conformers were refused at χ = −2 and −4. Their box-face density was below the isovalue, so neither surface was clipped. χ = 2 conflates an **open** mesh (clipped) with a **closed mesh of higher genus**, such as a folded chain touching itself at the isovalue.

**Why it matters.** An ensemble mode built on this check silently drops the most compact shapes, which differ most in shape descriptors (F3).

**Correction.** Test topology directly:
- an edge used by one triangle means open: refuse;
- an edge shared by more than two triangles means non-manifold: refuse;
- a closed surface of genus > 0 passes, with a note;
- more than one closed piece passes, with a warning.

**Validation — done** `[M]`. Recomputed with the same geometry (`bench/recheck_refused_conformers.py`):

| conformer | χ | pieces | open edges | non-manifold edges | genus | max ρ on box face |
| --- | --- | --- | --- | --- | --- | --- |
| 40 | −2 | 1 | 0 | 0 | 2 | 7.5×10⁻⁶ |
| 72 | −4 | 1 | 0 | 0 | 3 | 7.3×10⁻⁶ |

Both are single closed surfaces far inside the box (isovalue 2×10⁻³): valid, and wrongly refused before. The fix is on `bench/timing` (`mesh_topology`, `surface_problems`), with 6 fast tests covering a torus, two separate pieces, an open mesh, box-face density on a closed mesh, a sphere-like mesh and a non-manifold edge. The existing slow test that a clipped surface is refused is unchanged.

## Direction: reactivity descriptors (after comparison mode)

ESP shows where the electrostatic potential is negative or positive, which is not the same as where a molecule reacts. For teaching nucleophilic and electrophilic sites, and as the "free substrate" baseline when asking why an enzyme picks one site over another, espmap's "property evaluated on the surface vertices" design extends directly. In order of cost:

1. **Average local ionization energy (ALIE)** on the same surface. Where electrons are most easily removed: electrophilic attack, C–H oxidation sites. Computed from the existing orbitals at almost no extra cost.
2. **HOMO / LUMO amplitude** on the surface. Where electrons are donated and accepted. Same SCF.
3. **Fukui functions (f⁺, f⁻) and the dual descriptor.** Where the density actually changes on adding or removing an electron. Two extra SCFs at fixed geometry, and the anion and cation are **open-shell**, so this is the first real use of the UKS path.

**Calibration set: cases where ESP alone points to the wrong site.** Run each, predict from the map, then compare with the known outcome:
- enolate: hard electrophiles go to O, alkyl halides mostly to C;
- thiocyanate (S vs N);
- nitrite (N vs O).

A descriptor earns its place by getting these right where ESP does not.

**Limits, stated in the tool:**
- **Enzyme selectivity cannot be computed from the free substrate.** Geometry, active-site electrostatics, protonation and dynamics need the protein (QM/MM), which is out of scope. ESP, ALIE, frontier orbitals and Fukui functions are intrinsic electronic *descriptors*, not reaction preferences. espmap provides a free-substrate electronic baseline; enzyme selectivity emerges from how that baseline combines with substrate geometry, active-site electrostatics, protonation, dynamics and reaction barriers.
- **Charged species** (protonated carbonyls, carboxylates) need implicit solvent and their own scale (F2). Gas-phase cations exaggerate electrophilicity.
- **Cofactors** (PLP Schiff bases and similar) are large and often zwitterionic: feasible on a GPU, but advanced examples.

## What it got right

These are the hard calls, checked:

- **"More DFT methods" ranked last.** `[V]` This matches the measurements: the functional spread is 0.69 kcal/mol, against +3.54 for basis and −3.60 for solvent. On the targets the basis step alone is +1.53 to +2.82 `[M]`.
- **"Comparable" is the right framing.** `[V]` The code already warned that auto-scaled panels "are NOT comparable". The roadmap turns that warning into a feature.
- **Protonation stays advisory.** No automatic pKa prediction is proposed.
- **The gate stays visible in any UI.** The notebook kept it (`[V]`, the structure check runs before any SCF).
- **3D surfaces are not averaged across conformers.**
- **Local chemistry vs global shape.** The roadmap's split holds for extended conformers; folding moves the extrema too (F3).

## Summary

| # | finding | severity | state |
| --- | --- | --- | --- |
| F1 | cost unpriced | HIGH | **measured**; estimates retracted; C18 TZVP ≈ 15 min (CPU) |
| F2 | shared scale vs invariant 7 | HIGH | rule enforced in notebook; `compare` to implement |
| F3 | ensemble evidence from glycerol | MEDIUM | **measured**: holds for extended conformers (1.1); folding widens V_min to 5.2; two claims retracted |
| F4 | SVG of a raster render | MEDIUM | open; `--html` covers on-screen use |
| F5 | JSON ranked below its dependant | MEDIUM | partial; add shape descriptors |
| F6 | alignment of flexible chains | LOW–MED | open; plus default orientation |
| F7 | batch gating order | LOW | open (with `compare`) |
| F8 | xTB ≠ cheap DFT | LOW | **measured**: SVP preview; density fitting ≤ 0.02, size-dependent speed |
| F9 | χ = 2 refuses folded surfaces | MEDIUM | **fixed and verified** on both refused conformers |

**Revised order:**
1. extend `--json` with shape descriptors (F5): the record, surface and site
   statistics are done in v0.2.0; Rg, stereo and conformer ID remain
2. comparison mode with the F2 rule, F7 gating and head-group alignment (F6)
3. ~~fix the surface check (F9)~~ done
4. ensembles reporting shape descriptors (F3)
5. protonation advisory
6. reactivity descriptors: ALIE → HOMO/LUMO → Fukui
7. output polish (F4)
8. browser UI
9. other properties
10. more functionals

## Sources

All in-repo:
- `bench/results.jsonl`, `bench/results_clean_rerun.jsonl`, `bench/env.json` (timings, ESP, environment)
- `bench/conformers.jsonl` (conformer spread, the two refused surfaces)
- `bench/gpu_l4.json` (GPU validation)
- `esp_map.py`, `README.md`, `CLAUDE.md` (invariants 1–3, 7, 8), `docs/adversarial-review-1.md`

No external literature was retrieved. The xTB, ALIE and Fukui statements are general knowledge, marked `[I]` where they make a claim.
