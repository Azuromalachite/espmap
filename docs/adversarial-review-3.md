# Adversarial review 3 — espmap v0.1.0

Reviewer 3, 2026-10-22. Third pass, and the first on the whole repository as released: `esp_map.py`, `tests/`, `bench/`, `notebooks/`, `README.md`, `CLAUDE.md` and `docs/`.

**Reviewed:** release v0.1.0, commit `4bddf7beb3cb4a5bf960d14cc88cd24eddab118e`, from a read-only export of that tree. Every number below comes from that tree. The working checkout was not touched.
**Environment:** WSL2 Ubuntu, Python 3.12.15, numpy 2.5.4, scipy 1.18.1, pyscf 2.14.0, rdkit 2026.09.1, scikit-image 0.26.0, matplotlib 3.11.2, plotly 7.1.0, cmcrameri present, 16 threads, no GPU. These are the versions in `bench/env.json`, so I did not test any other library version. The GPU path was checked by reading the code and was not run.
**Standards applied:** the project's own invariants (CLAUDE.md), its `computational-chemistry-properties` and `molecular-surface-visualisation` rules, and the defect classes from reviews 1 and 2.
**Reproduction scripts:** all in `review3/work/`: `e_fast.py` (no SCF), `e_cli.sh` (end-to-end CLI runs, def2-SVP), `e_scf.py` (DF, PCM, dipole), `e_conf.py` (conformer provenance), `e_bench.py` (recomputes the README from `bench/*.jsonl`), and `e_fix.py` / `e_fix3.py` (fault injection of the proposed corrections). Total compute was about 6 minutes.

Evidence labels: **[X]** verified by execution · **[S]** verified in source · **[I]** inferred.

---

## Baseline: what passes as shipped

| check | result |
| --- | --- |
| `pytest -q` | **45 passed**, 11 skipped, 1.50 s [X]. README says 37 and CLAUDE.md says 28 (H12) |
| `esp_map.py --selftest` | PASS: water −49.34 / +59.14 / −0.22 kcal/mol; all 11 colormaps oriented [X] |
| `--selftest` with the kernel sign flipped *inside* `esp_at_points` (monkeypatched) | **returns 1**: +4189.76 / +3634.62 / **+265.87** [X]. The guard works on the real function, not only on a copy. |
| `build_notebook.py` → `.ipynb` | **byte-identical** to the shipped notebook (sha1 `859c214d8c`). All 8 code cells parse. No stored outputs, every `execution_count` is null [X] |
| conformer provenance (`bench/conformers.jsonl`) | all 10 conformer ids, Rg values and MMFF ΔE reproduce to 6 dp. The lowest-energy conformer is still id 4 [X] |

The physics core is sound. The defects are at the edges: the inputs it accepts, the outputs it writes, and claims that have drifted away from the data.

---

## Lead finding

### H1 · A `.pdb` file computes a different molecule, and every check passes `[BLOCKER]` · CONFIRMED EXECUTION FAILURE

`read_structure` reads PDB with `Chem.MolFromPDBFile(path, removeHs=False)` (esp_map.py:536–537). A PDB file without CONECT bond orders is the normal case for a ligand from a QM program, a docking run or a crystal structure. Proximity bonding then makes every bond single. Each aromatic carbon is left with three neighbours, and RDKit fills the gap with implicit hydrogens. I wrote chlorobenzene with all 12 atoms and explicit H, in the `HETATM` layout (`review-3/fx/chlorobenzene.pdb`):

```
$ esp_map.py --file chlorobenzene.pdb --check
--- structure check — ClC1CCCCC1 ---
  nothing to flag: unambiguous, closed-shell, single molecule
verdict: ACCEPTED

$ esp_map.py --file chlorobenzene.pdb --basis def2-svp --json pdb.json
chlorobenzene.pdb  ->  C6H11Cl, charge +0
  geometry: file 3D coordinates used as-is; 6 hydrogen(s) added
SCF -694.986512 Ha ... V = -19.23 .. +31.89 kcal/mol
sigma-hole V_max = +1.07 kcal/mol
```

The pipeline computed **chlorocyclohexane with a flat ring**. It took the five real hydrogens from the file and placed six more by `AddHs(addCoords=True)`. It reported a σ-hole of +1.07, where chlorobenzene's is +6.35 at this basis. Every check passed:

- **The gate** saw a valid closed-shell molecule, because the reader handed it one.
- **The approved-versus-computed assertion** (invariant 3) compared the reader's output with itself (H2).
- **The JSON** records `approved_smiles: "ClC1CCCCC1"`, `formula: "C6H11Cl"`, with nothing to show that it is not the file's molecule.

This is review 2's G3 again (`.xyz` lost its connectivity) for the one format G3 did not cover. The program has one reader, so this is not a two-copies defect (G4). The reader has a hole for one format, and the gate cannot see it because the gate reads through that same reader.

**Correction.**
1. Send `.pdb` through the same route as `.xyz`: take its coordinates, run `rdDetermineBonds.DetermineBonds(m, charge=charge)`, then `SanitizeMol`. Fault-injected [X]: the same file comes back as `Clc1ccccc1`.
2. Add an independent guard in `from_file` that costs one line. If the file already contains explicit hydrogens and `AddHs` adds more, raise `BLOCK H_COUNT_MISMATCH`. The reader's bond orders then disagree with the file. Fault-injected [X]: it trips on the PDB (5 H in the file, 11 after `AddHs`). It does not trip on `examples/real_3d.mol` (5 → 5) or `examples/real_3d_noH.mol` (0 → 5, a legitimate H-less file).

**Validation.** Add a fast test that writes this PDB fixture and asserts the result is `Clc1ccccc1`, or a BLOCK. Never `ClC1CCCCC1`.

---

## Findings

### H2 · The approved-versus-computed assertion cannot fail on the SMILES path, and never looks at the coordinates `[HIGH]` · CONFIRMED STATIC CODE DEFECT, fault-injected

Review 2's G7 fix made invariant 3 an assertion (esp_map.py:1279–1286). Here is what each side is built from:

- `approved = MolToSmiles(RemoveHs(probe))`, where `probe = MolFromSmiles(a.smiles)` (esp_map.py:1253).
- `computed = CanonSmiles(gi["smiles"])`, and `from_smiles` sets `gi["smiles"] = smiles`, **the raw input string** (esp_map.py:189).

Both sides are the same deterministic function of `a.smiles`. The assertion compares the user's input with itself. On the `--file` path both sides come out of `read_structure` on the same file, which is why H1 passes. On neither path does the check read `spec`, the coordinates PySCF actually receives. Review 2 tested it by swapping the molecule returned by the geometry step. That changes `gi`. A defect that changes only the coordinates is invisible to it. Fixture [X]:

| injected defect | shipped check | SMILES perceived from the coordinates sent to PySCF | proposed check |
| --- | --- | --- | --- |
| none (CCO) | pass | `CCO` | pass |
| none (acetate) | pass | `CC(=O)[O-]` | pass |
| coordinates of CC=O, approval of CCO | **pass** | `CC=O` | **fail** |
| coordinates of the E isomer, approval of `C/C=C\C` | **pass** | `C/C=C/C` | **fail** |

This is review 1's *seal over a copy* in a new form: the seal now covers the input string.

**Correction.** Perceive the computed structure from `spec`. Build an XYZ block, call `DetermineBonds(charge)`, then `AssignStereochemistryFrom3D`. Compare the result with `approved`: in full where the approved structure specifies stereo, stereo-stripped where it does not. Write the perceived, stereo-specified SMILES to the JSON (see H6). If bond perception fails, record `structure_verified: false` rather than skipping the check silently.
**Validation.** The four rows above, as a fast test. Both injected rows must refuse, and both healthy rows must pass.

### H3 · On an ion, the default colormap paints a negative surface blue, which the stated convention reads as positive `[HIGH]` · CONFIRMED EXECUTION FAILURE

`choose_scale` correctly picks a sequential scale when the field does not change sign. But `main` still renders with `a.cmap`, which defaults to `vik`, a diverging map (esp_map.py:1199, 1385). `choose_scale` only adds a note (esp_map.py:941–943). On acetate [X]:

```
$ esp_map.py --smiles "CC(=O)[O-]" --basis def2-svp --html ... --json ...
note: note: 'vik' is a diverging map being used on a sequential scale; prefer lajolla/batlow/viridis.
V = -171.79 .. -79.86 kcal/mol
wrote acetate.png  (sequential scale -171.8 to -79.9 kcal/mol, cmap vik)
```

| V (kcal/mol) | t | vik colour |
| --- | --- | --- |
| −170 | 0.00 | dark red |
| −125 | 0.50 | **white, the map's neutral** |
| −80 | 1.00 | **dark blue** |

The rendered PNG (`review-3/out/acetate.png`) shows the methyl end dark blue. The README, the notebook header and the HTML colourbar title (`"V (kcal/mol)<br>warm = negative"`, esp_map.py:1008) all tell the reader that blue means positive. The surface is −80 kcal/mol everywhere in that region. The PNG has no colourbar to correct the impression.

The sequential branch was added so that an ion would not be shown with a misleading zero (invariant 7). The default flags bring the misleading zero back, at −125 kcal/mol instead of at 0. The notebook's display cell switches to `lajolla` (build_notebook.py:200). The `esp.png` and `esp_3d.html` it zips for download come from the CLI with the default `vik`, because the command in cell 4 never passes `--cmap`.

**Correction.** Make `--cmap` default to `None`. Resolve it after `choose_scale`: `vik` when the scale is diverging, `lajolla` when it is sequential. If the user explicitly asks for a diverging map on a sequential scale, raise rather than add a note, or require `--force-cmap`. Make the HTML colourbar title depend on the scale kind: "darker = more negative" for sequential.
**Validation.** `--smiles "CC(=O)[O-]" --json` must record `scale.cmap == "lajolla"`. A fast test should assert that `main`'s resolved cmap is in `SEQUENTIAL` whenever `kind == "sequential"`, unless the user explicitly overrode it.

### H4 · A clipped surface rendered with `--no-strict` produces a JSON indistinguishable from a valid one `[HIGH]` · CONFIRMED EXECUTION FAILURE

`--json` is offered "for anything you will compare later" (README:88–89). Water, def2-SVP, `--margin 2.0 --no-strict` [X]:

| | `--margin 2.0 --no-strict` | default margin |
| --- | --- | --- |
| stderr | `isosurface integrity failed: … CLIPPED … OPEN (266 edges …)` | — |
| JSON `results.V_min` | **−34.51** | −44.63 (**23 % error**) |
| JSON integrity fields | `"chi": -2` | `"chi": 2` |

The JSON records nothing else: no strict flag, no margin, no `face_max_rho`, no `open_edges`, no integrity verdict. And `chi = −2` is exactly the value of `bench/f9_recheck.jsonl` conformer 40, a **valid** closed genus-2 surface that the code now accepts (esp_map.py:675–704). The only integrity field kept in the JSON is the one the code itself stopped using as a criterion. It cannot separate a damaged surface from a valid one.

The same pattern covers the gate. `--smiles "CC(=O)[O-].[Na+]" --force` is REJECTED with MULTIPLE_FRAGMENTS, then computed. The JSON writes the rejected structure under the key **`approved_smiles`**, with no gate verdict and no list of findings [X]. A downstream reader sees an approved structure.

**Correction.** Add to `results`: `surface = {closed, open_edges, nonmanifold_edges, components, genus, face_max_rho, margin_bohr, strict, problems: [...]}`. Add to `input`: `gate = {verdict, forced: bool, findings: [[severity, code], …]}`. Rename the key to `gate_smiles`, or write `approved_smiles` only when the verdict is not BLOCK.
**Validation.** The clipped water run must produce `surface.problems` that is non-empty and `strict: false`. The Na-acetate run must produce `gate.verdict == "BLOCK"` and `forced: true`. A healthy run must produce empty problems.

### H5 · A planar molecule with real 3D coordinates is discarded as a "2D depiction" `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

`from_file` declares a file flat when `max|z| < 1e-6` (esp_map.py:250–254). It discards the coordinates whenever `flat or not is3d` (esp_map.py:257, 273–278). A planar molecule in the xy plane is flat in that sense. Many QM outputs and builders write planar molecules exactly that way. [X]:

| file | header | shipped |
| --- | --- | --- |
| MMFF chlorobenzene rotated into z = 0, written as MOL | `RDKit 3D` | **DISCARDED**, ETKDG regenerated, *"2D depiction"* |
| water in the xy plane, `.xyz` | — | **DISCARDED**, CLI prints *"this file is a FLAT 2D depiction, not a geometry"* |
| the same chlorobenzene tilted 10° | `3D` | used as-is |
| the same coordinates, `.xyz` with 6 dp (max\|z\| = 4.5e−6) | — | used as-is |

The last two rows show that the threshold is a knife-edge on the number of decimals in the file. An XYZ or PDB file is never a drawing, yet it can be classified as one. The README promises *"a file carrying real 3D coordinates is used as-is"*. The JSON then records `geometry: "2D depiction: drawn coordinates DISCARDED"` for a user-supplied DFT geometry. The numeric damage is small for rigid planar molecules: the geometry axis is 0.35 kcal/mol. The provenance is wrong, and the user's geometry is silently replaced.

**Correction.** `.xyz` and `.pdb` are never depictions. For MOL/SDF/CDXML, believe the header when it says 2D. When the header says 3D and z ≡ 0, call the file a depiction only if the molecule has more than one bond type and every bond length is equal within 3 %. Drawing programs draw every bond at the same length; a real geometry has C–H ≈ 1.09, C–C ≈ 1.39 and C–Cl ≈ 1.74 Å. Fault-injected [X] on 9 files: planar 3D MOL → kept; xy-plane water XYZ → kept; `examples/drawn_2d.mol`, `examples/drawn.cdxml`, an RDKit 2D depiction and **a 2D depiction relabelled 3D** → still discarded; `real_3d.mol`, `real_3d_noH.mol` and `water.xyz` → kept. My first version of this rule misclassified water (two identical O–H bonds). The bond-type condition is what fixed it.
**Validation.** Those 9 files as a parametrised fast test.

### H6 · The JSON leaves out the reported number, mislabels another, and over-states comparability `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

Read as a data contract (esp_map.py:1412–1431):

- **The σ-hole V_max is printed and never written.** It is the quantity every error-budget row is about. The record holds the method and leaves out the one result a paper would quote.
- **`dipole_D` of an ion depends on where the molecule happens to sit.** Acetate as embedded gives 6.423 D. The same molecule translated by 5 Å gives **30.348 D**, with an identical energy [X]. PySCF takes the origin at (0,0,0), and the ETKDG coordinates sit wherever ETKDG put them. For a charged species the number has no meaning unless the origin is stated.
- **`V_min`/`V_max` have no unit in the key**, while `e_tot_Ha` and `spacing_bohr` do. The NPZ holds bohr coordinates and kcal/mol potentials with no unit, basis or functional metadata at all.
- **The stereoisomer actually computed is not recorded.** Threonine with both centres unspecified gives `C[C@H](N)[C@H](O)…` for seeds 0xC0FFEE, 2, 3, 5 and 7, and its enantiomer for seeds 1, 4 and 6 [X]. The JSON says `CC(N)C(O)C(=O)O` in every case. Enantiomers have mirror-image ESP, so the extrema are unaffected. The *map* is a different picture depending on the seed, and the record cannot say which.
- **`comparable: true` depends only on the colour flags** (`scale_provenance`, esp_map.py:953–967). Two runs at `--lim 25` are "comparable" across def2-SVP and def2-TZVP, gas phase and PCM water, or different `--iso` values. Those are exactly the axes the README says move the numbers by 3.5 kcal/mol.

**Correction.**
- Add `sigma_hole = {V_max_kcal_mol, n_cap_vertices, i_sub, i_bond}`, or `null` with the reason (no halogen / F / more than one halogen).
- Rename the keys to `V_min_kcal_mol` and `V_max_kcal_mol`.
- For a charged species, compute the dipole about the centre of nuclear charge and write `dipole_origin`. Or write `null` with a note.
- Write the coordinate-perceived stereo SMILES from H2.
- Rename `comparable` to `colour_scale_fixed`. Add `compare_key = sha1(xc, basis, solvent, iso, spacing, cmap, vlim)`, so a consumer can only compare runs whose keys match.
- Store `units` and the method in the NPZ.

**Validation.** Translating acetate by 5 Å must leave `dipole_D` unchanged to 1e−6. Chlorobenzene must produce a σ-hole value in the JSON.

### H7 · The solvent error-budget row measures a convention, and the convention is not stated `[MEDIUM]` · METHODOLOGICAL CONCERN, measured

With `--solvent`, `run_scf` polarises the solute density with PCM (esp_map.py:608–612). `esp_at_points` then evaluates only the solute's nuclear and electronic potential. The potential of the PCM surface charges (the reaction field) is left out. That is a legitimate convention: the ESP of the polarised solute. But it is the convention the "−3.61, co-dominant with basis" row depends on, and it appears nowhere in the README, the JSON or the CLI help. Chlorobenzene, def2-SVP, PCM water, reaction field summed over the PCM point charges (`e_scf.py` S8) [X]:

| | σ-hole V_max (kcal/mol) |
| --- | --- |
| gas | +6.36 |
| PCM water, polarised solute only (what espmap reports) | **+3.25** (Δ −3.10) |
| PCM water, solute + reaction field at that vertex | **+6.93** (Δ **+0.58**) |

The reaction-field potential on the ρ = 0.002 surface spans **−18.3 to +13.0 kcal/mol**, larger than the effect the row reports. Under the other convention, the sign of the "solvent effect" reverses. A caveat on my measurement: the charges are discretised points on the cavity, and the vertices sit close to some of them, so the reaction-field values are order-of-magnitude. That still settles the question of which convention the row depends on.

**Correction.** Name the convention in the solvent row, in the `--solvent` help text and in the JSON: `settings.solvent_esp = "polarised solute only (no reaction field)"`. Optionally add `--include-reaction-field`, using `mf.with_solvent` charges.
**Validation.** Re-measure the solvent row under both conventions at def2-TZVP and print both in the README table.

### H8 · Claims recomputed from `bench/` do not all match, and seven of the eleven budget rows have no data in the repository `[MEDIUM]` · CONFIRMED by recomputation (`e_bench.py`)

| claim (where) | data | verdict |
| --- | --- | --- |
| def2-SVP → TZVP moves V_max of "decanoic acid, BDSF and a C18 pheromone acetate by +2.65 to +2.82" (README:186–187); "about +2.7 low for these molecules" (README:65); "+2.7 [M]" (roadmap-review:179) | decanoic +2.81/+2.83, BDSF +2.69/+2.71, Z3Z13-18OAc +2.64/+2.65, **Z3Z13-18OH +1.53** (V_min **+2.37**), all from `results.jsonl` | **The fourth molecule in the same file is left out.** The measured range is +1.53 to +2.83. The omitted row is the one that agrees with the README's own message: the shift is molecule-specific. |
| `--density-fit` "halves the SCF time up to ~600 basis functions" (README:68) | TZVP: 0.53, 0.46, 0.45 · **SVP: 0.72, 0.78, 0.80** | true at TZVP only |
| DF ≤ 0.02 kcal/mol, 7 pairs, V_max always upward (README:177) | 7 pairs, max \|ΔV\| 0.0176, V_max upward 7/7 | **holds.** I extended it to the σ-hole: Δ = **+0.0003** (SVP) and **+0.0002** (TZVP) [X] |
| GPU < 0.001 kcal/mol (README:178) | gpu_l4.json ΔV −0.00055/+0.00068 | holds, against a DF CPU run |
| Conformer: 1.1 over 8 extended; 5.2 on V_min once folded; dipole 43 %; negative area 13 % / 60 % | 1.098 / 1.175; 5.204; 43 % of mean; 13 % of the 8-conformer mean, 60 % of the 10-conformer mean | holds. Different denominators between 13 % and 60 %; say so. |
| basis +3.55, solvent −3.61 / −2.69, functional 0.69, stereo 0.9 / 0.4, geometry 0.35, grid 0.04, region 0.00 | **no file in `bench/`** | prose only (reviews 1–2, docstring). I reproduced +6.35 (SVP) and +9.83 (TZVP) myself [X]; the other rows cannot be checked. |

**Correction.** Report +1.53 to +2.83 across the four molecules, and the V_min shifts too. Change "halves" to "halves at def2-TZVP, saves 20–28 % at def2-SVP". Commit a `bench/error_budget.py` and its `.jsonl` for the seven unbacked rows. Until then, mark those rows `[prose]` in the README table.
**Validation.** A test that recomputes every number in the README table from `bench/*.jsonl`. The branch after this release adds a claim-consistency test, but it checks that the documents agree with each other, not that they agree with the data.

### H9 · The `--gpu` contract is written down and not enforced `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT

The `run_scf` docstring, README:74–75 and README:218 say `--gpu` is validated only for closed-shell, gas-phase, density-fitted runs. `run_scf` (esp_map.py:595–621) and `main` accept every combination silently. The notebook makes the untested combinations the easy path:

- The notebook's metadata requests an L4 accelerator, and `GPU = "auto"` turns the GPU on whenever there is one (build_notebook.py:89–90).
- A student who picks `SOLVENT = water`, sets `SPIN > 0` or unticks `DENSITY_FIT` therefore lands on the untested path, with no message.
- `--optimize --gpu` passes the CPU object back from `to_cpu()` into `berny_solver.optimize`, so the gradients run on the CPU. Nothing says so.
- Both validation runs in `gpu_l4.json` were recorded at commits `2313965` and `183ccd7`, not at the release. Nothing in the release tree records which commit it validated, so it cannot be checked from the tree.

**Correction.** In `run_scf`, warn when `gpu` is used with `spin`, a solvent or `not density_fit`, naming the untested combination. Pass `gpu_validated: bool` through to the JSON. The notebook should print the same warning, or force density fitting when `USE_GPU` is on.
**Validation.** A fast test with gpu4pyscf stubbed: `run_scf(..., gpu=True, solvent="water")` must emit the warning.

### H10 · The notebook: the install marker ignores BRANCH, the cache key ignores the colour range, an all-ion series is not refused, and the limits text is stale `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT (one item by fixture)

- **The marker file hides a branch change.** `ESPMAP_READY` is checked before cloning and is not keyed to `BRANCH` (build_notebook.py:97–105). After a first run, changing `BRANCH` in cell 3 neither clones nor reinstalls, and the old code keeps running without any message.
- **The cache key omits the colour range.** `settings` (build_notebook.py:142–144) leaves out `LIM`, but the command passes `--lim` (build_notebook.py:161). Change `COLOUR_RANGE` and re-run, and cell 4 reuses the folder. The zipped `esp.png` and `esp_3d.html` keep the first run's range, while the display cell shows the new one. The cache is also judged valid by `surface.npz` alone (build_notebook.py:150). `main` writes the NPZ before the JSON (esp_map.py:1399–1434), so an interrupted run leaves a folder the cache will trust [I].
- **An all-ion series is drawn on a symmetric scale.** `if one_sign and len(one_sign) < len(surfs)` (build_notebook.py:280) refuses mixed series only. A series made entirely of ions passes, and is drawn on a symmetric ±60 diverging scale: the exact failure the markdown cell says it refuses. Logic lifted verbatim, acetate + benzoate → `refused=False` [X].
- **The limits text is stale.** It says *"ESP extrema move little with conformation"* (build_notebook.py:357), which contradicts README:190–195 (5.2 kcal/mol once folded) and the CLI's own conformer note. It says `preview` is *"about ±3.5 kcal/mol in absolute values"* (build_notebook.py:57), turning one signed, molecule-specific measurement into a symmetric bound the data does not support (H8: +1.53 to +3.55, and −0.76 for bromobenzene).
- **The runtime claim has no basis.** *"Decanoic acid ~1 min"* on CPU (build_notebook.py:33) is the 16-thread laptop figure (48 s, SVP + DF). Colab's free CPU runtime has 2 vCPUs [I].

**Correction.**
- Write `BRANCH` and the resolved commit into the marker, and reinstall on mismatch.
- Add `LIM` and the espmap commit to the cache key, and make `settings.json` (written last) the cache sentinel.
- Refuse when `one_sign` is non-empty at all, or render an all-ion series on a shared *sequential* scale.
- Copy the conformer and basis wording from the README.

**Validation.**
- With an all-ion series, the cell must raise.
- Changing `COLOUR_RANGE` must create a new job folder.
- Changing `BRANCH` must print a reinstall.

### H11 · `choose_scale` is a knife-edge at zero `[LOW]` · CONDITIONAL EDGE CASE

`spans_zero = data_min < 0.0 < data_max` (esp_map.py:906). A cation surface from −0.01 to +120 gets a diverging ±120 scale, which uses **50 %** of the ramp. From +0.01 to +120 it gets a sequential scale over 100 % [X, verbatim function]. A 0.02 kcal/mol change, far below the 3.5 kcal/mol basis error, decides the scale kind. The scale kind also decides the JSON's `kind` and `comparable` fields. A large ion with a slightly negative patch is precisely the case invariant 7 was written for.

**Correction.** Call the field sign-changing only if the minority sign reaches more than a set fraction of the larger extremum, for example `min(|Vmin|,|Vmax|) / max(...) > 0.05`, and print the ratio. Record the ratio in the JSON.
**Validation.** Parametrised test at ratios 0.01, 0.04, 0.06 and 0.5.

### H12 · Tests that guard something adjacent, and no test of `main` `[LOW]` · CONFIRMED STATIC CODE DEFECT

- `euler_characteristic` occurs once in `esp_map.py`, at its own definition (esp_map.py:669). Nothing calls it. `test_euler_characteristic_*` (test_fast.py:242–251) test dead code. CLAUDE.md still lists it as the surface stage, and invariant 8 still says "χ = 2", which the code no longer requires (genus > 0 is accepted).
- `test_a_flipped_coulomb_kernel_is_caught` (test_slow.py:36–52) writes the broken formula out by hand and shows that it fails. It never injects the defect into `esp_at_points` or `selftest`. I did, and the guard does fire (baseline table). The test should do the same.
- `test_isosurface_is_closed_with_adequate_margin` asserts `chi == 2`, not `closed`.
- `main` has no test at all. That leaves untested: invariant 2 (gate before compute), invariant 3 (H2), the JSON contract (H4, H6), the cmap resolution (H3), and `--density-fit`, which is not exercised by any test, fast or slow.
- Test counts: 28 in CLAUDE.md, 37 in the README, **45 collected**.

**Correction.**
- Delete `euler_characteristic` and its tests, or route `mesh_topology` through it.
- Turn the kernel test into a monkeypatch of `esp_at_points` and assert `selftest() == 1`.
- Add one fast end-to-end test of `main` on H₂ / STO-3G (seconds), with `--json`, asserting the keys from H4 and H6.
- Add one slow DF-versus-exact test on the water probes.

**Validation.** Each new test must be red on v0.1.0 and green after its fix.

### H13 · Documentation contradicts the code it ships with `[LOW]` · CONFIRMED STATIC CODE DEFECT

| statement | location | reality |
| --- | --- | --- |
| DF "effect on the ESP is unmeasured" / "ESP effect not yet measured" | esp_map.py:586–587, 1182 | measured: ≤ 0.02 (README), σ-hole +0.0003 (H8) |
| "extrema are robust, shape is not"; conformer row "~1" | CLAUDE.md error-budget table and the paragraph under it | 5.2 kcal/mol once folded (README, docstring, roadmap-review retraction 4) |
| "radicals (this is restricted DFT)" rejected; "open-shell species is refused" | esp_map.py:37–38; docs/colormaps.md §4 | UKS with `--spin` since v4 |
| "gas phase throughout" | esp_map.py:136 | `--solvent` exists |
| water reference windows "−40 to −55", "+40 to +60" | docs/colormaps.md validation table | selftest asserts [−70, −30] and [+35, +80], and test_slow.py calls them "the published range" with no source cited |
| "esp_map.py v4" | esp_map.py:3 | the package is 0.1.0, and the JSON carries neither number |

**Correction.** Fix the statements above. The claim-consistency test on the next branch should also scan the module docstring, the CLI help strings and `docs/colormaps.md`.

### H14 · Reproducibility: what a third party could not reconstruct from the outputs `[LOW]` · CONFIRMED STATIC CODE DEFECT

- **The JSON has no geometry.** It holds no coordinates and no geometry hash. Coordinates exist only in the NPZ, and only if `--npz` was passed. ETKDG output depends on the RDKit version, so a seed is not a substitute for coordinates.
- **The JSON omits versions and settings that move results.** It records pyscf and rdkit only, not numpy, scipy, scikit-image, Python, BLAS or the thread count. It also leaves out `grids.level`, `conv_tol` and `margin`.
- **`espmap_commit` is `git rev-parse HEAD` in the module's directory** (esp_map.py:1409–1411). It is `null` for any non-git install (verified: `null` for this export [X]). It reports a clean HEAD for a modified working tree, because no dirty flag is recorded. For a non-editable install inside some other repository it would report that repository's HEAD [I]. The package version, 0.1.0, is not written at all.
- **Bench rows carry no per-row code identity.** `env.json` was captured once, at `553e72c` on `bench/timing`, and `capture_env.py` itself says `timing.py` "re-reads itself for every row", meaning it could change between rows. The GPU rows are values transcribed from a printed table (stated honestly in the file).
- **The notebook pins PySCF but not RDKit** (build_notebook.py:101). The GPU check's `geom_sha` will therefore differ on another RDKit version, and the cell prints `DIFFERENT` and carries on comparing energies of different geometries.

**On the positive side**, verified: the seeded embedding reproduces `bench/conformers.jsonl` exactly (every conformer id, Rg and MMFF ΔE to 6 dp), and the def2-SVP/TZVP σ-holes of review 1 (+6.35 / +9.83) reproduce to 0.01.

**Correction.** Add the following to the JSON:
- `geometry = {xyz_angstrom, sha1}`;
- the full version set and `threads`;
- `espmap_version` from `importlib.metadata`;
- `git describe --dirty`, with `null` written as "unknown", never as a clean commit.

Add `commit` to every bench row. In the notebook, pin `rdkit==` next to `pyscf==`, and have the GPU check cell stop when the geometry hashes differ.
**Validation.** Re-running a JSON's recorded geometry through `--atoms` must reproduce `V_min`/`V_max` to 1e−6.

---

## Prior findings: were they actually fixed?

| | status in v0.1.0 |
| --- | --- |
| F1–F9 (review 1) | fixed and still holding. `--selftest` catches the sign flip injected into the real function. Clipping is refused by default (but see H4 for `--no-strict`). |
| G1, G2, G3, G4, G5, G6, G8, G9, G10 (review 2) | fixed and still holding |
| **G7, the gate sealed a copy** | **not effectively fixed.** The assertion exists, but on the SMILES path it compares the input string with itself, and on the file path it compares the reader with itself (H2). H1 is the consequence. |
| roadmap F9, folded surfaces refused | fixed by `mesh_topology`. The JSON did not follow: it still records only `chi` (H4). |

---

## What it got right

- **The guard that matters most is real.** I injected the historical kernel-sign flip into `esp_at_points` itself, not into a copy. `--selftest` returned 1 and showed the original +265.87 far-field signature. Many projects have a self-test that cannot fail. This one can, and it fails on the bug it was written for.
- **ESP at the surface vertices.** Re-verified: no cube files anywhere in the tree, and the interpolation error is absent by construction.
- **Seeded, reproducible geometry.** Ten conformers out of a 100-conformer C18 ensemble, regenerated bit-for-bit from a seed. Few benchmark directories can be regenerated this way, and it is why this review could verify the conformer claims instead of trusting them.
- **The benchmark keeps its failures.** The two `OSError` rows from a full `/tmp` stay in `results.jsonl`. The contaminated 1,613 s DF row is kept beside a clean re-run, and the README quotes the clean one. Conformers 40 and 72 were recorded as refused, then re-examined (F9) instead of dropped, and the retraction in roadmap-review says that the check had been removing exactly the cases that contradicted the claim.
- **`mesh_topology` separates the two things χ conflated.** An open mesh (one edge, one triangle) and a closed higher-genus surface now give different answers. The torus fixtures in the tests show the distinction was understood rather than patched.
- **`gpu_l4.json` states its own limits:** `not_covered`, and a note on the precision of transcribed values. H9 is about the code not enforcing what that file says, not about the file.
- **Density fitting was measured before it was recommended.** I extended the measurement to the σ-hole: Δ = 0.0003 kcal/mol. The README's claim was right, and stronger than it said.

---

## Summary table

| # | finding | severity | evidence |
| --- | --- | --- | --- |
| H1 | `.pdb` without CONECT → chlorobenzene computed as flat-ring chlorocyclohexane; gate, assertion and JSON all pass | **BLOCKER** | CONFIRMED EXECUTION FAILURE |
| H2 | approved-versus-computed assertion is tautological on SMILES, same-reader on files, blind to coordinates | HIGH | CONFIRMED STATIC, fault-injected |
| H3 | default `vik` on an ion's sequential scale: −80 kcal/mol painted blue ("positive") | HIGH | CONFIRMED EXECUTION FAILURE |
| H4 | `--no-strict` clipped run (V_min 23 % off) and `--force`d BLOCK both produce clean-looking JSON | HIGH | CONFIRMED EXECUTION FAILURE |
| H5 | planar 3D geometry discarded as a "2D depiction"; knife-edge on file precision | MEDIUM | CONFIRMED EXECUTION FAILURE |
| H6 | JSON: no σ-hole; ion dipole origin-dependent (6.4 → 30.3 D); no units; stereoisomer unrecorded; `comparable` ignores method | MEDIUM | CONFIRMED EXECUTION FAILURE |
| H7 | solvent row depends on an unstated no-reaction-field convention; the other convention reverses its sign | MEDIUM | METHODOLOGICAL CONCERN, measured |
| H8 | basis range omits the 4th molecule (+1.53); "halves" false at SVP; 7/11 budget rows have no data | MEDIUM | CONFIRMED by recomputation |
| H9 | GPU validated envelope not enforced; notebook defaults into it | MEDIUM | CONFIRMED STATIC |
| H10 | notebook: install marker ignores BRANCH; cache ignores range; all-ion series not refused; stale limits | MEDIUM | CONFIRMED STATIC + fixture |
| H11 | `choose_scale` knife-edge at 0 | LOW | CONDITIONAL EDGE CASE |
| H12 | dead code under test; kernel test not injected; no test of `main`, the JSON or DF | LOW | CONFIRMED STATIC |
| H13 | docstrings and docs contradict the code and the README | LOW | CONFIRMED STATIC |
| H14 | JSON/bench lack geometry, full versions, dirty flag, per-row commit | LOW | CONFIRMED STATIC |

## Verdict

**Blocks the release:**
- **H1.** An advertised input format computes the wrong molecule and every safeguard approves it.
- **H2**, because it is the reason H1 was invisible. Fixing H2 properly (perceive from `spec`) would not have caught H1 alone, since both sides would agree on the wrong molecule. That is why H1 needs the reader fix *and* the H-count guard.
- **H3 and H4** block anyone who uses the image or the JSON as the record: one paints negative as positive under default flags, the other hides a damaged surface.

**Should be fixed before the next number is quoted:** H6 (σ-hole missing from the record), H7 (state the solvent convention), H8 (report the fourth molecule).

**Can wait:** H5, H9–H14.

**The structural lesson, and it differs again from the last two rounds.** Review 1 found checks that measured something next to the conclusion. Review 2 found rules implemented twice. This round, the checks inside the pipeline are mostly sound; **what fails is everything at the boundary.** A format the reader half-understands; an assertion that compares the input with itself; a colormap default that ignores the scale decision made one line earlier; a JSON that keeps the topology field the code stopped trusting. Each boundary takes a correct internal result and passes on something less. The guard for that is not another internal check. It is to derive what goes *out* (the record, the image, the approval) from the same object that was computed, and to test `main` end to end.
