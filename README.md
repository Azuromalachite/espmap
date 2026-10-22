# espmap

Electrostatic-potential maps for organic molecules, from a SMILES string or a
structure drawn in any editor — with a colour scale chosen from the data rather
than from habit.

```bash
pip install -e .
python esp_map.py --selftest
python esp_map.py --smiles "Clc1ccccc1" --lim 25 --out esp.png
```

No install needed to try it: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Azuromalachite/espmap/blob/master/notebooks/espmap_colab.ipynb)
(type a SMILES, then *Runtime → Run all*; a rotatable 3D surface, an image and a zip of results come back).

![six colormaps on one molecule](docs/esp-colormap-comparison.png)

*The same bromobenzene surface, identical geometry, viewpoint and lighting in
every panel — only the colormap differs. Rows: normal vision, simulated
deuteranopia, greyscale.*

> **Scope.** This is research code, written for one project and published
> because the pipeline is reusable. It works and it is tested, but it is **not
> a maintained package**: there is no support commitment and no stability
> guarantee. Fork it rather than depend on it.

---

## What it does

1. DFT single point (PySCF)
2. electron density on a regular grid
3. marching cubes at ρ = 0.002 e·bohr⁻³ — the Bader surface, ~98 % of the density
4. **the ESP evaluated only at the isosurface vertices**, not on a volumetric grid
5. a z-buffered render with deliberately weak shading

Step 4 is the one worth borrowing. The conventional route writes two cube
files — density and potential — then interpolates the second onto an
isosurface of the first. Evaluating the potential at the ~10⁴ surface vertices
instead of ~10⁶ grid points is about 100× cheaper *and* removes the
interpolation error entirely, because the potential is computed exactly where
it is displayed. Bromobenzene, end to end, on one core:

| step | time |
| --- | --- |
| SCF (141 basis functions) | 20.0 s |
| density on an 81×70×39 grid | 0.7 s |
| marching cubes → 10,400 vertices | <0.1 s |
| **ESP at those vertices** | **2.5 s** |
| render | 1.2 s |

No Gaussian, no GaussView, no cube files on disk.

## Speed: preview, density fitting, GPU

Measured on a 16-thread laptop CPU (Ryzen 7 260, 7.4 GB in WSL2), B3LYP, SCF
time only (`bench/`):

| molecule | basis functions (TZVP) | def2-SVP | def2-TZVP |
| --- | --- | --- | --- |
| chlorobenzene | 253 | — | 39 s |
| decanoic acid | 492 | 66 s | 211 s |
| C18 pheromone acetate | 898 | 182 s | 875 s |

- **`--basis def2-svp`** is the preview: right for pictures, about +2.7 kcal/mol
  low on V<sub>max</sub> for these molecules, so do not quote numbers from it.
- **`--density-fit`** changed V<sub>min</sub>/V<sub>max</sub> by ≤ 0.02 kcal/mol in every pair
  measured, and halves the SCF time up to ~600 basis functions. At 898 it was
  *slower* (1,175 s vs 875 s), most likely because its integrals went to disk;
  time it on your machine before relying on it.
- **`--gpu`** runs the SCF on an NVIDIA GPU via
  [GPU4PySCF](https://github.com/pyscf/gpu4pyscf). On a Colab L4 it matched the
  CPU to < 0.001 kcal/mol and was 4.3× faster than the laptop above
  (`bench/gpu_l4.json`). Validated for closed-shell, gas-phase, density-fitted
  runs only.

## Outputs

```bash
python esp_map.py --smiles "CCO" --out esp.png \
    --html esp_3d.html --npz surface.npz --json results.json
```

- `--html`: an interactive 3D view (rotate, hover to read V, see-through
  toggle). Same colour map and scale rules as the image; works offline.
- `--npz`: the surface (vertices, faces, V, atoms) so it can be re-rendered
  without another SCF.
- `--json`: settings, energy, dipole, V range, package versions, git commit
  and the error-budget caveat, for anything you will compare later.

## Input

```bash
python esp_map.py --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"   # caffeine
python esp_map.py --file drawn.cdxml                        # ChemDraft / ChemDraw
python esp_map.py --file structure.mol --check              # gate only, no compute
```

`--file` reads `.mol .sdf .cdxml .xyz .pdb`. A structure **drawn** in an editor
and exported as MOL/SDF/CDXML is a picture, so its coordinates are discarded,
a conformer is generated, and the run says so loudly. A file is treated as a
drawing when its header says 2D, or when it says 3D but every *z* is zero and
every bond has the same length across two or more element pairs. A planar
*geometry* — benzene optimised in the *xy* plane, water from an `.xyz` — is
kept; `.xyz` and `.pdb` are never drawings. A flat file with a 3D header and
one element pair (nitrite, nitrate, CO₂) cannot be told apart from a drawing
and is refused with `DIMENSION_AMBIGUOUS`. A file carrying real 3D
coordinates is used as-is, with any missing hydrogens added and relaxed
against frozen heavy atoms.

A file with more than one record (an SDF of docking poses) is refused; pass
`--record N` to compute record *N*. The JSON records which one under
`input.record`.

## The structure gate

Every structure is checked **before** any geometry or SCF work. A drawing can
be perfectly legible to a chemist and still be ambiguous to a machine, and the
pipeline would otherwise compute a flawless map of a molecule you did not mean.

| | |
| --- | --- |
| **Rejected** | valence and charge errors; R-groups and query atoms; radicals without `--spin`; disconnected fragments |
| **Warned** | undrawn E/Z geometry; two or more unwedged stereocentres; untested elements |
| **Noted** | one unspecified stereocentre; reachable tautomers; net charge; isotopes |

Each finding names what to change *in the structure*, not in the code. A salt
is rejected rather than averaged: forcing sodium acetate through reaches
**+113 kcal/mol** against **+61** for acetic acid's own proton — the
counter-ion *is* the map.

## The colour scale

**The general argument is not new.** The case against rainbow colour maps in
science is twenty years old and thoroughly made — see
[Crameri, Shephard & Heron, *Nature Communications* **11**, 5444 (2020)](https://doi.org/10.1038/s41467-020-19160-7)
and [Moreland's work on diverging maps](https://www.kennethmoreland.com/color-maps/).
What this repository adds is the measurement for *this* field, where the
defaults still ignore it.

Perceptual units are CAM02-UCS ΔE′ (≈2 is just noticeable), measured on a real
ESP surface:

| colormap | **chroma at V = 0** | lightness at V = 0 | 0↔negative, deuteranopia |
| --- | --- | --- | --- |
| rainbow (GaussView-style) | **33.7** | 90 | 59 |
| red–green–blue (3-stop) | **39.4** | 68 | **13** |
| RdBu | 2.4 | 98 | 73 |
| vik (Crameri) | 1.1 | 93 | 73 |
| PuOr | 2.3 | 98 | 64 |

A diverging scale must be *least* saturated where the data is least
meaningful. The conventional maps are near their own maximum chroma at zero,
so the eye is pulled hardest to the value that means nothing is happening. And
because the midpoint colour is not the least-chromatic colour, the perceived
neutral is **displaced to t = 0.727** along the ramp rather than 0.5.

Three rules the tool enforces:

- **Symmetric about zero — only if the field changes sign.** An ion's surface
  never does. Benzoate runs −158.7 to −40.7 kcal/mol; a diverging map centred
  on zero collapses 117 kcal/mol of structure into one flat colour. `espmap`
  detects this and switches to a sequential scale.
- **No specular highlight.** A white patch at a position set by the camera
  reads as a potential extremum. Shading is held to a narrow 0.62–1.0 span so
  the colormap's own lightness carries the data.
- **Orientation asserted in code.** All eleven colormaps are checked by
  `--selftest`; three shipped inverted before that check existed.

Honest consequence, measured: **no symmetric diverging map survives
greyscale** — RdBu separates its ends by ΔE′ = 3 and vik by 7, because
symmetry in lightness is exactly what makes ±V print identically. If the
figure will print in one colour, add a zero isoline rather than changing hue.

## Error budget

A property computed to three decimals is worthless without one. Measured on
the chlorobenzene σ-hole V<sub>max</sub>, B3LYP unless stated:

| axis | effect | note |
| --- | --- | --- |
| **basis set** (def2-SVP → def2-TZVPP) | **+3.55** | +6.35 → +9.90; dominant |
| **solvent** (gas → PCM water) | **−3.61** | +9.83 → +6.22; co-dominant |
| solvent (gas → chloroform) | −2.69 | +7.14 |
| functional (B3LYP / PBE0 / M06-2X / HF) | 0.69 | spread, +9.51 to +10.20 |
| conformer (lowest vs highest MMFF) | ~1 | glycerol; C18 acetate: 1.1 over 8 extended conformers, **5.2** on V<sub>min</sub> once folded ones are included |
| density fitting (`--density-fit`) | ≤ 0.02 | 7 pairs, always upward on V<sub>max</sub> |
| GPU vs CPU (`--gpu`) | < 0.001 | decanoic acid, density-fitted |
| stereochemistry (diastereomers) | ~0.9 | against a 0.4 noise floor |
| geometry source (MMFF / idealised / DFT) | 0.35 | +9.59 to +9.94 |
| grid spacing (0.26 → 0.16 bohr) | 0.04 | converged |
| σ-hole region definition | 0.00 | max is interior to the cone |

**The basis-set error is molecule-specific and cannot be calibrated away**: the
same def2-SVP → def2-TZVP change moves bromobenzene by **−0.76**, opposite in
sign and a fifth the magnitude, and moves V<sub>max</sub> of decanoic acid, BDSF and a C18
pheromone acetate by +2.65 to +2.82. Checking one molecule licenses nothing
about another.

**Conformation moves shape more than extrema, until the molecule folds.**
Across eight extended conformers of the C18 acetate, V<sub>min</sub>/V<sub>max</sub> varied by
~1.1 kcal/mol while the dipole varied 43 % and the negative surface area 13 %.
The two folded conformers moved V<sub>min</sub> by up to 5.2 kcal/mol and the negative
area by 60 % (`bench/conformers.jsonl`, `bench/f9_recheck.jsonl`). For a
flexible molecule, sample folded shapes before quoting an extremum.

Two consequences for anyone quoting a number from this:

- the default basis is **def2-TZVP**, and anything double-zeta prints a warning
  naming the measured penalty;
- the default is **gas phase**, stated in every run, and `--solvent water`
  changes the answer by as much as the basis does.

## Limits

- `ωB97X-D` is unavailable in this PySCF build, so range-separated functionals
  with dispersion are **untested**.
- Geometries from `--smiles` are force-field (MMFF94s), adequate for local
  electronic properties at 0.35 kcal/mol but **not** for conformational
  energies or barriers. `--optimize` runs a DFT optimisation (needs
  `pip install 'espmap[opt]'`, which adds pyberny) and is exercised on one
  molecule only.
- The gate checks that a structure is *unambiguous*. It cannot tell you it is
  the molecule you meant to draw: a cleanly drawn wrong tautomer passes with a
  note.
- Fluorine deliberately reports no σ-hole, because it does not form one.
- Only closed-shell and simple open-shell (`--spin`) species; no excited
  states, no periodic systems, no transition metals validated.
- `--gpu` is validated only for closed-shell, gas-phase, density-fitted runs;
  solvent and open-shell runs on the GPU are untested.
- A folded molecule can give a closed surface with handles (genus > 0). That
  is accepted with a note; an *open* (clipped) surface is still refused.

## Tests

```bash
pip install -e ".[test]"
pytest                      # 100 fast tests, ~5 s, no SCF
pytest --runslow -m slow    # 13 tests with real DFT, ~3 min
```

Every test is a defect that actually shipped at some point, lifted from the two
adversarial reviews in `docs/`. The reviews are included deliberately: they
record what was measured, what was wrong before, and what remains untested.
Among the findings they document — a sign-flipped Coulomb kernel caught only
by far-field decay, a validation check that could not fail, a clipped
isosurface that changed a value by 26 % silently, and a conclusion whose margin
was smaller than its own basis-set error.

## Credit

Built on [PySCF](https://pyscf.org), [RDKit](https://rdkit.org),
[scikit-image](https://scikit-image.org) and
[cmcrameri](https://github.com/callumrollo/cmcrameri). The `vik`, `broc`,
`cork`, `lajolla` and `batlow` maps are Fabio Crameri's *Scientific colour
maps*; please cite them if you use the figures.

MIT licensed. If it is useful in published work, a link back is appreciated but
not required.
