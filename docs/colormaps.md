# ESP maps without the rainbow

> **CORRECTED 2026-10-22.** An adversarial review retracted this document's §2.2 "what that costs, concretely" claim about chlorobenzene, downgraded its density-integral validation row, and restated the apparent-neutral figures in range-independent form. Corrections are marked inline; the reasoning is in `esp-adversarial-review.md`. The code is now v2.

**Question:** how hard is it to render 3D electrostatic-potential maps for organic molecules without the conventional red–green–blue scheme?

**Answer:** the 3D part took 23 seconds of compute and about 200 lines. The colormap is a one-line change. The convention is not hard to fix — it is just never fixed, because every package ships it as the default.

Everything below was computed and measured here, 2026-10-21, not quoted. Code: `esp_map.py`. Figures: `esp-colormap-comparison.png`, `esp-colormap-diagnostics.png`.

---

## 1 · How hard, concretely

Pipeline for bromobenzene, B3LYP/def2-SVP, on one CPU:

| Step | Tool | Time |
| --- | --- | --- |
| SCF single point (141 basis functions) | PySCF | 20.0 s |
| Electron density on a 81×70×39 grid (221k points) | PySCF `eval_rho` | 0.7 s |
| Marching cubes at ρ = 0.002 e·bohr⁻³ → 10,400 verts / 20,796 tris | scikit-image | <0.1 s |
| ESP at the 10,400 **surface vertices** | PySCF `int1e_grids` | 2.5 s |
| Render 700², supersampled | ~120 lines of numpy | 1.2 s |
| **Total** | | **23.3 s** |

Install is one line: `pip install pyscf scikit-image matplotlib cmcrameri colorspacious`. No Gaussian licence, no GaussView, no VMD, no cube files on disk.

**The implementation choice that matters most, borrowed rather than invented.** Do not compute the ESP on a volumetric grid and then interpolate it onto the surface, which is what the cube-file workflow does. Compute the density on a grid, extract the isosurface, then evaluate the ESP **only at the surface vertices** — 10⁴ points instead of 10⁶. Quantitative molecular-surface methods evaluate properties on polygonal surfaces in this way (Bulat et al. 2010; Lu & Chen 2012; see [prior-art.md](prior-art.md)); the measurement below is this project's. That is the difference between 2.5 seconds and most of an hour, and it removes the interpolation error entirely, because the potential is evaluated exactly where it is displayed.

### Validation

Physics checked before any picture was made, because an ESP map with a sign error looks perfectly plausible:

| Check | Computed | Expected |
| --- | --- | --- |
| Water, lone-pair lobe at 1.5 Å | **−49.3** kcal/mol | −70 to −30 (selftest window) |
| Water, 1.0 Å beyond H | **+59.1** | +35 to +80 (selftest window) |
| Water, far field at 25 Å | **−0.2** | → 0 |
| Bromobenzene σ-hole V<sub>max</sub> | **+16.8** | +12 to +16 |
| Bromobenzene Br equatorial belt V<sub>min</sub> | **−12.9** | negative belt |
| σ-hole trend, Cl → Br | **+6.4 → +16.8** | increases down the group |

**Retracted from this table:** "Density integral vs electron count, 76.93 / 76, ±1 %". That check cannot fail usefully — it varies −5 % to +6.8 % *non-monotonically* with grid spacing while V<sub>max</sub> moves 0.5 %, and is identical to two decimals (59.46 vs 59.47) across a basis change that moves V<sub>max</sub> by 3.5 kcal/mol. It samples the nuclear cusp, not the valence density at ρ = 0.002 where the surface lives. The ±1 % was an artifact of the spacing I happened to pick. v2 reports it as a diagnostic and gates nothing on it; the real check is `--selftest`.

**Missing from this table, and it was the one that mattered:** basis-set convergence. Nothing here tests it, and §2.2 below depended on it.

Two errors caught on the way, both of which would have produced a confident, wrong figure: `int1e_grids` returns the *positive* Coulomb kernel, so the electronic term must be subtracted (the far-field value failing to decay is what exposed it); and my first validation probe sat on the hydrogen side of water, not the lone-pair side.

---

## 2 · Why the convention is wrong — measured, not asserted

The test molecule is bromobenzene because its σ-hole is the ideal stress case: a **small positive cap, 1.5 % of the surface area**, sitting inside a negative belt *on the same atom*. It is the feature that explains halogen bonding, and it is exactly what a bad colormap destroys.

All six colormaps below are oriented warm/red = negative = electron-rich, so the entrenched chemistry reading is preserved. Fixing perception should not silently flip the sign convention too.

Perceptual units are CAM02-UCS ΔE′ (≈2 = just noticeable):

| Colormap | J′ at V=0 | **C′ at V=0** | 0↔neg | ·CVD | greyscale sign |
| --- | --- | --- | --- | --- | --- |
| Rainbow (GaussView-style) | 90 | **34** | 82 | 59 | 15 |
| Red–green–blue (3-stop) | 68 | **39** | 67 | **13** | 26 |
| RdBu (ColorBrewer) | 98 | 2 | 80 | 73 | **3** |
| vik (Crameri) | 93 | 1 | 77 | 73 | **7** |
| PuOr | 98 | 2 | 66 | 64 | 23 |
| PINDR div-5→div-1 | 88 | 3 | 49 | 40 | 16 |

### 2.1 The convention puts maximum saturation on the least meaningful value

A diverging scale must be *least* saturated at its midpoint, so the zero crossing recedes and the extremes advance. The conventional maps do the opposite: **C′ = 34 and 39 at V = 0**, against 1–3 for every properly constructed diverging map, and within a hair of their own maximum chroma (43–45). Green at zero is the loudest colour on the surface. The eye is drawn hardest to the place where nothing is happening.

This is the cleanest single statement of the problem, and it has nothing to do with colour blindness.

### 2.2 The perceptual zero is in the wrong place — in normal vision

Locating each map's chroma minimum, which is where a reader perceives "neutral":

| Colormap | apparent neutral | error |
| --- | --- | --- |
| Rainbow | **+5.8** kcal/mol | +5.8 |
| Red–green–blue (3-stop) | **+9.6** kcal/mol | +9.6 |
| RdBu | −0.5 | ~0 |
| vik | +0.2 | ~0 |
| PuOr | −0.3 | ~0 |
| PINDR div | +1.0 | ~0 |

The three-stop RGB map's green midpoint is *not* its least-saturated colour; the neutral point falls about ten kcal/mol into the electrophilic region. Deuteranopia barely moves it (+8.6), so this is a bias in the colormap, not an accessibility problem.

> **RETRACTED.** This paragraph previously read: *"Chlorobenzene's σ-hole is +6.4 kcal/mol. The apparent neutral of the RGB map is +9.6. So the entire σ-hole of chlorobenzene sits below the perceived zero and renders as neutral-to-electron-rich… A convention that resolves a feature in one halide and inverts it in the next is not a convention worth keeping."*
>
> **It is not supported.** The +6.4 came from def2-SVP. Converged: **+9.83** (def2-TZVP), **+9.90** (def2-TZVPP) — a **+3.54** kcal/mol basis-set error against a claimed margin of 3.2. The σ-hole does not fall below the apparent neutral; it lands on it, within 0.3 kcal/mol.
>
> **What survives:** the RGB map renders chlorobenzene's σ-hole in very nearly its *least saturated* colour — it paints the molecule's electrophilic site the colour the map uses for "nothing here". That is still damning and now robust, but it is a coincidence between **one molecule and one colormap**, not a general result. §2.1 is the general finding and needed no chemistry.

**Two caveats on the numbers above.** The displacement is a fraction of the ramp — `t = 0.727` for the RGB map, `t = 0.637` for the rainbow — so the kcal/mol figures hold only at this molecule's range of ±21.2 and scale with it (the RGB map's apparent neutral is +4.5 at a ±10 range, +22.7 at ±50). And "apparent neutral" uses chroma minimum as a *proxy* for perceived neutrality, which is not validated here. The range-independent, proxy-free version of the finding is §2.1.

### 2.3 Red–green loses the neutral/electron-rich distinction entirely

Separation between V = 0 and the negative extreme: **67 → 13** under simulated deuteranopia. Red and green both collapse toward olive, so the *location of the zero crossing* disappears. The sign information in the three-stop map is carried almost entirely on the one axis ~8 % of men cannot resolve.

Rainbow loses less (82 → 59) only because it passes through cyan and yellow on the way.

### 2.4 The honest trade-off: nothing symmetric survives greyscale

This one went against me. The best perceptual maps are the **worst** in print: RdBu separates its two ends by ΔE′ = **3** in greyscale and vik by **7** — effectively invisible — while the crude RGB map manages 26.

That is not a flaw to fix, it is a theorem. A diverging map symmetric in lightness *must* render ±V as the same grey; symmetry is what makes the centre recede, and it is what destroys greyscale. If the figure will be printed in one colour, you cannot encode sign in the colormap at all. Add a zero isoline, contour the surface, or accept asymmetry (PuOr, at 23, is mildly asymmetric and survives best).

---

## 3 · What to use

**Default: `vik` reversed, or `PuOr`.** Both have C′ ≈ 1–2 at zero, a light centre, monotonic lightness toward each end, and no red–green dependence. `vik` if the figure is for screen; `PuOr` if it may be printed or if CVD readers are a stated concern.

**Non-negotiable regardless of hue:**

1. **Symmetric range about zero — for a surface that actually changes sign.** `vmin = −vmax`, otherwise the neutral colour is not at neutral potential and every statement about sign is unreliable. **This rule was stated unconditionally and that was wrong.** An ion's surface has no sign change at all: benzoate runs **−158.7 to −40.7** kcal/mol, entirely negative. Centring a diverging map on zero then maps the whole molecule to one end of the ramp and collapses 117 kcal/mol of real structure into one flat colour (see `esp-from-smiles.png`, third panel). A charged species needs a **sequential** scale over its actual range, oriented so the most electron-rich region is the most intense. v3 detects the no-sign-change case and switches automatically.
2. **State the range numerically** on the colourbar, in kcal/mol or au. Auto-scaling per molecule makes a weak σ-hole and a strong one look identical.
3. **Fix the range across a figure** when comparing molecules. Bromobenzene and chlorobenzene both top out near +21 kcal/mol — but on the *hydrogens*, not the halogen, so per-molecule auto-scaling silently rescales the feature you are actually discussing.
4. **Weaken the lighting.** On a lit 3D surface the shading multiplies lightness, and a colormap that uses lightness to carry sign is fighting the renderer. `ambient = 0.62, diffuse = 0.38, specular = 0` keeps geometry readable without swamping the data. Kill the specular highlight — it is a white patch at an arbitrary location that reads as a potential extremum.

**PINDR's own `div` ramp works** (C′ = 3 at zero, apparent neutral +1.0) and is the right choice if these figures sit alongside other project material. Two caveats, both real: it is a 5-stop ramp interpolated in sRGB rather than a perceptual space, so it is slightly less uniform than `vik`; and its sign separation is lower (62 vs 48 normal, but 40 vs 73 for the 0↔negative distinction), because its endpoints are less extreme. For a standalone chemistry figure `vik` is the better instrument. For a figure that has to look like it belongs to the thesis, `div` is defensible and I would use it.

---

## 4 · Usage

```bash
python esp_map.py --selftest                       # run this first
python esp_map.py --smiles "Clc1ccccc1" --lim 25 --out esp.png
python esp_map.py --file drawn.cdxml --lim 25 --out esp.png     # drawing-program export
python esp_map.py --smiles "c1ccccc1C(=O)[O-]" --cmap lajolla --out benzoate.png
python esp_map.py --molecule bromobenzene --cmap vik --lim 25 --out esp.png
python esp_map.py --atoms "O 0 0 0.117; H 0 0.757 -0.469; H 0 -0.757 -0.469" \
                  --cmap PuOr --lim 50 --out water.png
```

**`--check` gates every structure before any geometry or SCF work.** A drawing can be perfectly legible to a chemist and still be ambiguous to a machine, and the pipeline will otherwise compute a flawless map of a molecule you did not mean. **Rejected:** valence and charge errors, R-groups and query atoms, radicals unless `--spin` is given (then run unrestricted, UKS), and disconnected fragments — a salt has no single ESP surface and the counter-ion swamps it, sodium acetate forced through reaches **+113** kcal/mol against **+61** for acetic acid's own proton. **Warned:** undrawn E/Z geometry, two or more stereocentres with no wedges, untested elements. **Noted:** a single unspecified stereocentre, reachable tautomers, net charge, isotopes. Each finding names what to change *in the drawing*.

The stereochemistry wording is measured rather than assumed: threonine's two diastereomers differ in mean ESP extrema by **1.6–1.9 kcal/mol** at def2-SVP, against a spread of up to **1.25** between MMFF seeds of one diastereomer (`bench/error_budget.jsonl`), and enantiomers are mirror-identical by symmetry. An earlier "0.9 against a 0.4 noise floor" had no data behind it and did not reproduce (review 3 H8). Stereochemistry moves the numbers about as much as conformer choice, and the *shape* of the surface more.

`--file` reads `.mol .sdf .cdxml .xyz .pdb` — what ChemDraft, ChemDraw, Ketcher and MarvinSketch export — so you can draw a structure and map it. The load-bearing distinction is **2D vs 3D**: a drawn structure exported as MOL/SDF/CDXML is a flat picture with every *z* zero, and feeding it to an SCF would produce a planar molecule and a confidently wrong map. A 2D file therefore has its coordinates discarded, a conformer generated, and says so loudly; a file carrying real 3D coordinates is used as-is, with any missing hydrogens added and relaxed against frozen heavy atoms (verified to reproduce the all-atom result exactly, +6.61 either way).

`--smiles` builds 3D coordinates with RDKit (ETKDGv3 + MMFF94s, lowest of 20 conformers, fixed seed). Measured against two reference geometries for the chlorobenzene σ-hole: idealised bond lengths **+9.83**, SMILES+MMFF **+9.94**, B3LYP/def2-SVP optimised **+9.54** — a 0.40 kcal/mol spread, an order of magnitude below the basis-set error, so a force-field geometry is adequate for this quantity on this molecule (`bench/error_budget.jsonl`). It is not for small anions: MMFF moves their V_min/V_max by up to 6.5 kcal/mol (`bench/anion_geometry.jsonl`), and `--optimize` is the remedy there. Net charge and radical count are read from the SMILES; an open-shell species is refused unless `--spin` is given, and then run unrestricted (UKS) rather than silently as closed-shell. Conformer effect on the ESP *extrema* is small — ~1 kcal/mol across glycerol's full 8 kcal/mol conformer range, since the extrema are set by local chemistry — but the *shape* of the surface is of course conformer-specific, so the map is one conformer and says so. That ~1 kcal/mol holds only while the molecule does not fold: two folded C18 acetate conformers moved V_min by up to 5.2 kcal/mol (`bench/conformers.jsonl`, `bench/f9_recheck.jsonl`).

`--selftest` asserts the ESP sign and magnitude against known water values, and the orientation of every colormap; the v1 sign error was caught by hand and would otherwise have shipped silently. `--lim` sets the symmetric range in kcal/mol — omitting it auto-scales to the molecule and now says so loudly. `--margin` controls the vacuum around the nuclei; below ~4 bohr the isosurface clips and v2 refuses to proceed (at 3.5 bohr the σ-hole was wrong by 26 % with no other symptom). `--ambient` raises toward 1.0 to let the colormap rather than the light carry the data. `--iso` defaults to the Bader 0.002 e·bohr⁻³ surface.

Built-in geometries cover the halobenzenes and benzene; anything else goes in through `--atoms` as a PySCF atom string. Iodobenzene returns an empty σ-hole cap and therefore `nan`. **The cause is not diagnosed** — an earlier version of this document asserted it was the iodine ECP shifting the core density; that was a plausible story, not a finding, and no test was run.

### Limits of what was done here

- Geometries are constructed from standard bond lengths, **not optimised**. Measured sensitivity: dV/dR ≈ −14 kcal·mol⁻¹·Å⁻¹, so ±0.01 Å costs ≈ ±0.15 kcal/mol. The built-in C–Cl 1.740 Å happens to sit within 0.001 Å of the microwave value — luck, not method.
- **Basis set and solvent are co-dominant.** This originally read "basis-sensitive at the few-kcal/mol level — enough to matter *if* you wanted to publish it". It did not matter conditionally; it falsified the claim outright: +6.35 (def2-SVP) → +9.90 (def2-TZVPP). v2 defaults to def2-TZVP and warns on anything double-zeta.
- **Solvent is now measured and it is not a small axis.** Gas phase → PCM water moves the chlorobenzene σ-hole from **+9.83 to +6.22** kcal/mol (chloroform: +7.14) — as large as the basis-set error, and in the opposite direction. Gas phase remains the default and is stated in every run; `--solvent water` is available. For the benzoate anion the surface span widens from 117 to 148 kcal/mol.
- **Functional dependence is measured and it is small.** B3LYP +9.83, PBE0 +9.99, M06-2X +9.51, HF +10.20 — a 0.69 kcal/mol spread at def2-TZVP. ωB97X-D is unavailable in this PySCF build, so range-separated functionals with dispersion remain untested.
- The CVD figures are *simulations* via `colorspacious`, which models deuteranomaly at 100 % severity. They are a design check, not evidence about any particular reader.
