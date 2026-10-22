# espmap — working notes for Claude Code

Read this before changing anything. It records what is measured, what is
deliberately absent, and which defects have already shipped once so they are
not reintroduced.

## What this is

A single-module tool: SMILES or a structure file → DFT → electrostatic
potential on the ρ = 0.002 e·bohr⁻³ isosurface → a rendered image with a
colour scale chosen from the data.

`esp_map.py` (~1700 lines) is the whole program. There is no package layout
and that is intentional — it stays readable end to end, and `pyproject.toml`
exposes it as a `py-modules` entry plus an `espmap` console script.

**Status: research code, unmaintained.** Correctness matters; API stability
does not. Do not add a compatibility shim for a signature you changed.

## Run it

```bash
pip install -e ".[test]"
python esp_map.py --selftest           # physics + all 11 colormap orientations
pytest                                 # 100 fast tests, ~5 s, no SCF
pytest --runslow -m slow               # 13 tests with real DFT, ~3 min
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
| surface | `extract_surface`, `mesh_topology`, `surface_problems` | integrity-checked, raises by default |
| property | `sigma_hole` | axis from the bond, not from a ring |
| scale | `choose_scale`, `resolve_cmap`, `cmap_preflight`, `cmap_after_scf` | diverging only if the field changes sign; the map follows the scale |
| render | `look_at`, `render`, `build_maps` | own z-buffer rasteriser, no VTK |

### The implementation choice worth preserving

`esp_at_points` is called on the **isosurface vertices**, never on a
volumetric grid. ~10⁴ points instead of ~10⁶: about 100× cheaper, and it
removes the cube-file interpolation error because the potential is computed
exactly where it is drawn. If you ever find yourself writing a `.cube` file,
stop and re-read this paragraph. It is an established approach to
quantitative molecular surfaces (Bulat et al. 2010; Lu & Chen 2012), not an
espmap invention: describe it that way, and keep `docs/prior-art.md` current
when a feature builds on published methods.

## Invariants — do not regress these

1. **One reader, one embedder.** `read_structure` is the only place a file is
   parsed; `_embed_and_minimise` is the only place a conformer is generated.
   Both existed twice before, and a fix landed in one copy and not the other
   (the XYZ repair reached the gate but not the geometry path, so the gate
   accepted a file the next step crashed on).
2. **Gate before compute.** Parse → check → embed → SCF. Checking after
   embedding wastes conformer generation on structures that were never
   computable and buries the real error under force-field complaints.
3. **The approved structure is the computed structure.** `verify_computed`
   reads `spec`, the coordinates PySCF receives, and nothing else:
   connectivity (hydrogens included) perceived from distances must match the
   approval, and stereo assigned from 3D must be one completion of what the
   approval specified. The check it replaced compared the input string with
   itself (review 3, H2) and passed an injected acetaldehyde geometry under an
   ethanol approval. Bond orders are deliberately not compared: the SCF never
   sees them, and comparing them refused valid resonance forms (`[N-]=C=S`,
   DMSO). A mismatch refuses; with `--force` it warns and the JSON records
   `structure_verified: false`.
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
   sequential ones. "Changes sign" is `sign_change_ratio(V) > SIGN_CHANGE_MIN`
   (0.05), never `V.min() < 0 < V.max()`: that knife-edge let 0.02 kcal/mol
   decide the scale (review 3 H11). Anything that asks whether V changes sign,
   the notebook included, calls the helper.
8. **The isosurface must be closed and interior to its box**: no edge used by
   only one triangle, and no density above the isovalue on a box face. A
   clipped surface changed a measured V_max by 26 % with every other number in
   the run still plausible. Do **not** go back to requiring χ = 2: a folded
   chain can give a closed surface of higher genus (two compact C18 acetate
   conformers had χ = −2 and −4 with nothing clipped). `mesh_topology` tells the
   two apart; genus > 0 passes with a note.

## Measured error budget — chlorobenzene σ-hole V_max

Do not re-derive these; do not quote a number against a smaller margin.

| axis | effect (kcal/mol) |
| --- | --- |
| basis set, def2-SVP → def2-TZVPP | **+3.54** |
| solvent, gas → PCM water (polarised solute, the default) | **−3.60** |
| solvent, gas → chloroform (polarised solute) | −2.69 |
| solvent, gas → water / chloroform, + reaction field (`--include-reaction-field`) | +0.05 / −0.34 |
| functional, B3LYP / PBE0 / M06-2X / HF | 0.69 spread |
| conformer, lowest vs highest MMFF | ~1 (extended; 5.2 on V_min with folded C18 conformers, see below) |
| stereochemistry, threonine diastereomers (def2-SVP) | 1.6–1.9 (seed spread up to 1.25) |
| geometry source, MMFF / idealised / DFT | 0.40 (chlorobenzene); **up to 6.5** on small anions |
| grid spacing, 0.26 → 0.16 bohr | 0.04 |
| density fitting (`--density-fit`), 7 pairs up to 898 nao | ≤ 0.02, always upward on V_max |
| GPU4PySCF (`--gpu`) vs CPU, decanoic acid, both RI | < 0.001 |
| σ-hole region definition | 0.00 |

**The solvent row is a convention, not a fact about the molecule.** The
default map is the PCM-polarised solute alone; adding the PCM surface
charges' potential cancels almost all of the shift (+0.05 in water, def2-TZVP,
`bench/solvent_convention.jsonl`). Never quote −3.60 without naming the
convention, and never compare records whose `solvent_esp` differs (it is in
`analysis_key`).

**The basis error is molecule-specific.** The same change moves bromobenzene
by **−0.76** — opposite sign, a fifth the size. Validating on one molecule
licenses nothing about another.

Two counterintuitive results, both measured, both load-bearing in the warning
text: **extrema are robust, shape is not.** Conformation and stereochemistry
move V_max by ~1 kcal/mol because extrema are set by local chemistry, while
the spatial arrangement of the features changes a great deal.

**This holds only while the molecule does not fold onto itself.** Across 8
extended conformers of a C18 acetate the extrema moved ~1.1 kcal/mol; adding
the 2 folded ones (genus 2 and 3 surfaces) widened V_min to 5.2 and V_max to
2.5 kcal/mol, and the negative surface area from 69–79 to 69–118 Å²
(`bench/conformers.jsonl`, `bench/f9_recheck.jsonl`). Intramolecular contact
is local chemistry too. Do not quote the ~1 kcal/mol conformer row for a
flexible molecule without sampling folded shapes.

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

`docs/adversarial-review-1.md`, `-2.md` and `-3.md` record every
defect found so far, with the reproduction for each. If a test fails in a way
that looks familiar, it is probably in there.

## Around the module

- **`--gpu`** lives inside `run_scf`: `to_gpu()` before the kernel and
  `to_cpu()` after it, so `total_dm` never sees a CuPy array. Keep the switch
  there; a notebook that sets up its own GPU SCF crashed on exactly this.
- **`bench/`** holds the measurements behind the speed and error-budget
  numbers: `timing.py` (one process per row, resumable), `conformers.py`,
  `env.json` (the environment they were taken in), `gpu_l4.json`. Quote from
  these files, not from memory; add a row rather than overwrite one.
- **`notebooks/`**: edit `build_notebook.py`, then regenerate
  `espmap_colab.ipynb`. Never hand-edit the `.ipynb`. Every code cell is plain
  Python so the notebook can be executed locally before anyone opens Colab.
- **`docs/roadmap-review.md`** is the current priority order, with measured
  status per roadmap item.

## Good first tasks

See `docs/roadmap-review.md` (revised order) for the main line of work. Smaller
ones:

- Close ωB97X-D by swapping in a PySCF build that has it, or by using a
  different range-separated functional and saying which.
- Exercise `--optimize` on three more molecules and report the spread.
- Raise PySCF's `max_memory` for `--density-fit` runs and re-time the C18
  acetate: at 898 nao density fitting was slower than exact (1,175 s vs
  875 s) with peak memory under the default 4 GB cap.
- Reduce the slow suite's runtime; 4 minutes is near the limit of what anyone
  runs locally before pushing.
