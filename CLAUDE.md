# espmap — working notes for Claude Code

Read this before changing anything. It records what is measured, what is
deliberately absent, and which defects have already shipped once so they are
not reintroduced.

## What this is

A single-module tool: SMILES or a structure file → DFT → electrostatic
potential on the ρ = 0.002 e·bohr⁻³ isosurface → a rendered image with a
colour scale chosen from the data.

`esp_map.py` (~1200 lines) is the whole program. There is no package layout
and that is intentional — it stays readable end to end, and `pyproject.toml`
exposes it as a `py-modules` entry plus an `espmap` console script.

**Status: research code, unmaintained.** Correctness matters; API stability
does not. Do not add a compatibility shim for a signature you changed.

## Run it

```bash
pip install -e ".[test]"
python esp_map.py --selftest           # physics + all 11 colormap orientations
pytest                                 # 28 fast tests, ~1.5 s, no SCF
pytest --runslow -m slow               # 11 tests with real DFT, ~4 min
python esp_map.py --smiles "Clc1ccccc1" --lim 25 --out esp.png
```

`--selftest` is the fastest signal that something fundamental broke. Run it
before and after any change to `esp_at_points`, `run_scf`, `total_dm` or
`build_maps`.

## Architecture, in the order data flows

| stage | functions | notes |
| --- | --- | --- |
| read | `read_structure`, `parse_for_check`, `from_smiles`, `from_file` | **one** reader, **one** embedder — see Invariants |
| gate | `check_structure`, `print_check` | runs **before** geometry or SCF |
| electronic | `run_scf`, `total_dm`, `density_on_grid`, `esp_at_points` | RKS, or UKS when `spin != 0` |
| surface | `extract_surface`, `euler_characteristic` | integrity-checked, raises by default |
| property | `sigma_hole` | axis from the bond, not from a ring |
| scale | `choose_scale`, `SEQUENTIAL` | diverging only if the field changes sign |
| render | `look_at`, `render`, `build_maps` | own z-buffer rasteriser, no VTK |

### The one idea worth preserving

`esp_at_points` is called on the **isosurface vertices**, never on a
volumetric grid. ~10⁴ points instead of ~10⁶: about 100× cheaper, and it
removes the cube-file interpolation error because the potential is computed
exactly where it is drawn. If you ever find yourself writing a `.cube` file,
stop and re-read this paragraph.

## Invariants — do not regress these

1. **One reader, one embedder.** `read_structure` is the only place a file is
   parsed; `_embed_and_minimise` is the only place a conformer is generated.
   Both existed twice before, and a fix landed in one copy and not the other
   (the XYZ repair reached the gate but not the geometry path, so the gate
   accepted a file the next step crashed on).
2. **Gate before compute.** Parse → check → embed → SCF. Checking after
   embedding wastes conformer generation on structures that were never
   computable and buries the real error under force-field complaints.
3. **The approved structure is the computed structure.** `main` asserts the
   canonical SMILES of what the SCF will run against what the gate passed.
   Both sides go through `CanonSmiles`; comparing a canonical form against a
   raw input string produced a false positive that blocked valid input.
4. **`int1e_grids` returns the positive Coulomb kernel**, so the electronic
   term is **subtracted**. A sign flip gives a smooth, structured, publishable,
   entirely wrong map. The far-field decay is the detector — a flipped kernel
   leaves an undecayed monopole (+265.87 instead of ~0 kcal/mol at 25 Å).
5. **`total_dm(mf)` is the only consumer of `make_rdm1`.** UKS returns a
   `(2, nao, nao)` stack; the density grid and the ESP both want α+β.
6. **Colormap orientation is asserted, not reasoned about.** Warm = negative =
   electron-rich on diverging scales; V_min darkest on sequential ones. Three
   maps shipped inverted, one of them introduced by a "fix" based on
   reasoning. `--selftest` checks all eleven.
7. **Diverging and symmetric only when the field changes sign.** An ion's
   surface never does, and centring on zero collapses its entire range into
   one colour. `--lim` is the diverging range; `--vmin`/`--vmax` are the
   sequential ones.
8. **The isosurface must be closed and interior to its box** (χ = 2, and no
   density above the isovalue on a box face). A clipped surface changed a
   measured V_max by 26 % with every other number in the run still plausible.

## Measured error budget — chlorobenzene σ-hole V_max

Do not re-derive these; do not quote a number against a smaller margin.

| axis | effect (kcal/mol) |
| --- | --- |
| basis set, def2-SVP → def2-TZVPP | **+3.55** |
| solvent, gas → PCM water | **−3.61** |
| solvent, gas → chloroform | −2.69 |
| functional, B3LYP / PBE0 / M06-2X / HF | 0.69 spread |
| conformer, lowest vs highest MMFF | ~1 |
| stereochemistry, diastereomers | ~0.9 (0.4 noise floor) |
| geometry source, MMFF / idealised / DFT | 0.35 |
| grid spacing, 0.26 → 0.16 bohr | 0.04 |
| σ-hole region definition | 0.00 |

**The basis error is molecule-specific.** The same change moves bromobenzene
by **−0.76** — opposite sign, a fifth the size. Validating on one molecule
licenses nothing about another.

Two counterintuitive results, both measured, both load-bearing in the warning
text: **extrema are robust, shape is not.** Conformation and stereochemistry
move V_max by ~1 kcal/mol because extrema are set by local chemistry, while
the spatial arrangement of the features changes a great deal.

## Deliberately absent

- **ωB97X-D** — not in this PySCF build, so range-separated + dispersion is
  untested. Do not claim functional convergence beyond the four tested.
- **Excited states, periodic systems, transition metals.** `TESTED_ELEMENTS`
  is the validated set; anything outside warns rather than refuses.
- **`--optimize`** exists but is exercised on one molecule.
- **Fluorine reports no σ-hole**, because it does not form one. This is a
  deliberate exclusion in `main`, not an oversight.
- **Intent checking.** The gate verifies a structure is *unambiguous*, never
  that it is the molecule the user meant. A cleanly drawn wrong tautomer
  passes with a note.

## Conventions

- Potential in **kcal/mol** at the boundary; Hartree internally (`HA2KCAL`).
- Coordinates in **bohr** internally, Ångström at the CLI and in atom strings.
- Warm/red = negative = electron-rich. This is the chemistry reading and it is
  the opposite of the physics one; say which is in force in any caption.
- Every warning names what to change **in the structure**, not in the code.

## If you change the science

Any change to basis, functional, grid, surface definition or solvent handling
invalidates the table above. Re-measure the affected axis and update both this
file and `README.md` — the numbers appear in both and in the module docstring.

`docs/adversarial-review-1.md` and `docs/adversarial-review-2.md` record every
defect found so far, with the reproduction for each. If a test fails in a way
that looks familiar, it is probably in there.

## Good first tasks

- Close ωB97X-D by swapping in a PySCF build that has it, or by using a
  different range-separated functional and saying which.
- Exercise `--optimize` on three more molecules and report the spread.
- A conformer-ensemble mode: Boltzmann-weight the surface rather than picking
  the lowest-energy conformer. The measurement says extrema barely move, so
  the value is in the *shape*, which is exactly what a single conformer gets
  wrong.
- Reduce the slow suite's runtime; 4 minutes is near the limit of what anyone
  runs locally before pushing.
