# Adversarial review 4 — espmap at master `6c66bbf`

Reviewer 3 stance, 2026-10-22. Fourth pass. It covers the ~20 PRs merged since v0.1.2: the fixes for review 3's H1–H14, and the per-site work (P1 atom identity, P2 partition and statistics, caller-defined sites). P1, P2 and `--sites` had not been reviewed independently before this pass. I wrote none of the code.

**Reviewed:** master at `6c66bbf` (merge of PR #8), from a separate worktree on branch `review/4`. No repository file was modified. The only additions are this file and `docs/review-4/`.

**Environment:** WSL2 Ubuntu, `~/.venvs/espmap` (Python 3.12, numpy 2.5, pyscf 2.14.0, rdkit 2026.09.1, pyberny). These are the versions in `bench/env.json`, and I tested no other library version. All runs used `OMP_NUM_THREADS=4`. The memory budget was shared, so every SCF here is STO-3G and the largest molecule is glycine. Numbers marked *sto-3g* are evidence about logic, not about chemistry. The GPU path was read, not run (no GPU).

**Reproductions:** in `docs/review-4/`.
- `e_fast.py`: no SCF.
- `e_scf.sh` + `e_json.py`: CLI runs and record inspection.
- `e_glycine.sh`, `e_glycine_noopt.sh` + `e_opt.py`: the lead finding.
- `e_misc.py`: dipole, GPU warnings, kernel-sign injection.
- `e_nb.py`: notebook regeneration and a cross-cell name check.
- `e_bench.py`: per-site rows recomputed from `bench/site_partition.jsonl`.
- `e_mapstereo.py`, `e_pdbnoh.py`.

Selected records are in `docs/review-4/out/` and fixtures in `docs/review-4/fx/`. The scripts source an `env.sh` that holds absolute scratch paths, so adapt it before re-running.

Evidence labels: **[X]** verified by execution · **[S]** verified in source · **[I]** inferred.

---

## Summary verdict

The H1–H14 fixes are mostly real. I reproduced each original failure against current code and it no longer occurs:
- the PDB is read as chlorobenzene;
- an E-for-Z geometry is refused;
- acetate is drawn with `lajolla`;
- a clipped or forced record says `validity.surface/input: false` and withholds sites;
- the dipole is translation-invariant to 1e-6 D;
- the kernel-sign injection still trips `--selftest`.

The partition itself is sound. Areas are conserved, the site extrema equal the global extrema, and symmetric sites agree within mesh noise.

**The lead finding is the H2 defect class coming back through a door H2 did not cover.** `verify_computed` checks `spec` *before* `--optimize` replaces it. For gas-phase glycine zwitterion at sto-3g, the optimiser moves a proton and the program computes neutral glycine. The record still says:
- `structure_verified: true`;
- `computed_smiles: "[NH3+]CC(=O)[O-]"`;
- `validity: {input: true, surface: true, sites: true}`;
- the same `analysis_key` as the unoptimised zwitterion.

The two records differ by 25–29 kcal/mol in both extrema [X]. `--optimize` is the remedy the module docstring itself prescribes for the small anions that make up the P3 calibration set.

The other findings sit at the record boundary again:
- `analysis_key` omits the geometry method;
- an ungated `--atoms` record is `validity.input: true` (the optional-presence idiom);
- the atom-level partition has no budget, and on ethanol it varies by up to 100 % of an atom's area;
- the JSON has drifted from `docs/per-site-spec.md` §6, which P3 is meant to freeze.

Counts: **1 HIGH, 5 MEDIUM, 7 LOW** new findings. **H1–H14:** 8 FIXED (H1, H3, H4, H5, H6, H7, H9, H11) and 6 PARTIAL (H2, H8, H10, H12, H13, H14). None regressed.

---

## Baseline

| check | result |
| --- | --- |
| `pytest -q` | **241 passed**, 19 skipped (18 slow + 1 git-unreadable worktree), 9.8 s [X]. CLAUDE.md and README still say "100 fast tests" and "13 slow"; 18 slow are collected (N10) |
| `esp_map.py --selftest` | PASS: water −49.34 / +59.14 / −0.22 kcal/mol; 11 colormaps oriented [X] |
| `--selftest` with the electronic term *added* inside `esp_at_points` (monkeypatched) | **returns 1**; far field **+265.87** [X]. Invariant 4's detector still works on the real function |
| `build_notebook.py` → `.ipynb` | byte-identical to the committed git blob. The checkout differs only by CRLF (`core.autocrlf=true`). 8 code cells, all `execution_count` null, no outputs; every free name is bound by an earlier cell [X] |
| light slow tests (`-k "symmetric_nitrite or generated_polar or geometry_reproduces or json_record_carries or anion_is_drawn"`) | 6 passed, 22.9 s [X] |

---

## Disposition of review 3's H1–H14

| # | review 3 finding | verdict | evidence |
| --- | --- | --- | --- |
| H1 | `.pdb` read as chlorocyclohexane | **FIXED** | `review-3/fx/chlorobenzene.pdb` now reads as `Clc1ccccc1`, and `from_file` gives C6H5Cl [X]. The same file with its H removed is refused with `PDB_CONNECTIVITY` and a remedy, not guessed [X]. |
| H2 | approved-vs-computed check tautological | **PARTIAL** | All four review-3 rows behave: CCO and acetate pass; acetaldehyde coordinates under an ethanol approval and E-for-Z are refused [X]. **But the check runs on the pre-`--optimize` `spec`. With `--optimize`, PySCF receives coordinates nothing verified (N1).** |
| H3 | `vik` on an ion's sequential scale | **FIXED** | acetate, default flags: `sequential`, `lajolla` (r7) [X]; `--cmap vik` is substituted and recorded (slow test passes) [X] |
| H4 | clipped / forced records look valid | **FIXED** | water at `--margin 2.0 --no-strict`: `validity.surface: false`, `sites: null`, problems listed [X]. Forced Na-acetate: `gate.verdict BLOCK`, `forced: true`, `approved_smiles: null`, `validity.input: false`, `sites: null` [X]. Residual: the NPZ carries no validity (N8). |
| H5 | planar 3D input discarded | **FIXED** | all six review-3 fixtures classified correctly: planar MOL and xy-plane water kept; `rd2d`, `rd2d_flagged3d` drawings [X] |
| H6 | σ-hole, units, dipole origin, stereo, `comparable` | **FIXED, with a new defect class** | Acetate dipole: **2.432360 D at both 0 and +5 Å**. PySCF's origin-0 value is 5.680 → 29.613 D [X]. σ-hole and unit-suffixed keys present [X]. `comparable` was replaced by `analysis_key`/`render_key`, but **`analysis_key` omits how the geometry was made (N2).** |
| H7 | solvent convention unstated | **FIXED** | `solvent_esp` is in settings and `analysis_key`; `--include-reaction-field` exists; the reaction-field reader is checked against PySCF's `e_solvent` (slow test) [S] |
| H8 | claims vs bench | **PARTIAL** | `test_claims.py` recomputes every README/CLAUDE.md budget row from `bench/` and passes [X]. The echoes were not updated (N10): <ul><li>the notebook says "about 2.7 kcal/mol low on the molecules benchmarked" (data +1.53 to +2.82);</li><li>the notebook says density fitting "roughly halves the SCF time" under its def2-SVP default (data: 20–28 %);</li><li>`docs/colormaps.md:144` still quotes the retracted +9.59 / 0.35;</li><li>`roadmap-review.md:179` says "+2.7 [M]".</li></ul> |
| H9 | `--gpu` contract not enforced | **FIXED** | `run_scf(gpu=True, solvent="water")` warns about PCM and about exact integrals before failing on the missing gpu4pyscf [X]. `gpu_validated` is in the JSON [S]. The notebook prints the same warnings [S]. |
| H10 | notebook marker, cache key, ion series, stale text | **PARTIAL** | Fixed: <ul><li>the marker is keyed to BRANCH;</li><li>the cache key holds `lim` and the commit;</li><li>`settings.json` is the sentinel and is written last;</li><li>an all-ion series uses a shared sequential scale via `sign_change_ratio` [S].</li></ul> Not fixed: the stale "about 2.7" and "roughly halves" above. Edge [S]: if every series member fails, `series_scale([])` raises `min() iterable argument is empty`, which is re-raised as "ions and sign-changing molecules cannot share…". The message is misleading but the cell does stop. |
| H11 | `choose_scale` knife-edge | **FIXED** | ratios 0.010 / 0.040 → sequential; 0.060 / 0.500 → diverging; −0.01..+120 → sequential [X]. Note: the ratio is still an extremum ratio, so a single vertex decides it. The knife-edge moved to 5 % rather than going away. That is acceptable, and it was what review 3 proposed. |
| H12 | tests guard something adjacent | **PARTIAL** | Fixed: `main` is now tested end to end (fast with a stop-before-SCF hook, slow on sto-3g) [S]. Not fixed: <ul><li>`euler_characteristic` (`esp_map.py:1053`) is still dead and still has two tests (`test_fast.py:683–692`);</li><li>`test_a_flipped_coulomb_kernel_is_caught` still hand-writes the broken formula (`test_slow.py:36–53`) instead of injecting into `esp_at_points`. I did the injection, and it trips;</li><li>`test_isosurface_is_closed_with_adequate_margin` still asserts `chi == 2` (`test_slow.py:59`), the criterion invariant 8 retired;</li><li>no DF-vs-exact test exists.</li></ul> |
| H13 | docs contradict code | **PARTIAL** | DF "unmeasured", "gas phase throughout" and the v4 header are gone [S]. `docs/colormaps.md:144` still says "B3LYP/def2-SVP optimised +9.59 — a 0.35 kcal/mol spread". These are exactly the figures review 3 H8 retracted. `test_claims` does not scan `colormaps.md`, although review 3 asked it to. |
| H14 | reproducibility | **PARTIAL** | Fixed: <ul><li>JSON `geometry` round-trips through `--atoms` to 1e-6 (slow test passes) [X];</li><li>versions, threads and `git describe --dirty` are recorded [S].</li></ul> Not fixed: per-row bench identity. **All 59 bench rows written since the fix say `"code"/"commit": "unknown"`** (anion_geometry 15, error_budget 15, isosurface_enclosed 4, site_partition 25), because WSL git cannot read a Windows-created worktree. The 6 `solvent_convention.jsonl` rows behind the −3.60 headline have no code field at all [X]. |

---

## Lead finding

### N1 · `--optimize` computes a structure nothing verified, and the record says it was verified `[HIGH]` · CONFIRMED EXECUTION FAILURE

`main` verifies `spec` at `esp_map.py:2185` (`verify_computed`) and at `esp_map.py:2201` (`verify_mapping`). Then, at `esp_map.py:2260–2274`, `--optimize` replaces `spec` with the berny-optimised geometry and runs the SCF again. Neither check runs again. Invariant 3 says "`verify_computed` reads `spec`, the coordinates PySCF receives". After `--optimize`, the coordinates PySCF receives are not the ones that were read.

Glycine zwitterion, B3LYP/STO-3G, gas phase (`e_glycine.sh`, 58 s, 475 MB peak) [X]:

| | `--smiles "[NH3+]CC(=O)[O-]"` | same with `--optimize` |
| --- | --- | --- |
| gate | ACCEPTED | ACCEPTED |
| `input.structure_verified` | true | **true** |
| `input.computed_smiles` | `[NH3+]CC(=O)[O-]` | **`[NH3+]CC(=O)[O-]`** |
| `validity` | all true | **all true** |
| `analysis_key` | `62f86d4d6b…` | **`62f86d4d6b…`** |
| V_min / V_max (kcal/mol) | −73.77 / +78.45 | **−49.01 / +49.05** |
| structure perceived from the record's own `geometry` | zwitterion | **`NCC(=O)O`, neutral glycine** |
| `verify_computed(record geometry, approval)` | — | **False**: "the connectivity perceived from the coordinates is not the approved structure's" |

The gas-phase zwitterion is not a minimum, and the proton transfers. The program then reports a map of a different tautomer under the approval and the verification of the original one, with a site block attached.

There is a second, quieter loss on the same path. For a mapped input, `main` drops the atom maps whenever `a.optimize` is set (`esp_map.py:2347–2354`).

`[CH3:1][OH:2] --optimize` (r3) [X] produces:
- `sites.atoms[*].atom_map` all `null`;
- `identity.mapped_smiles: null`;
- `identity.mapped_match: true`;
- `geometry_provenance.hydrogens: {"optimised": 4}`.

Berny does not reorder atoms, so the maps are still valid; they are discarded. In addition:
- `settings.geometry` still reads "3D generated from SMILES by ETKDGv3";
- `--sites "3+H"` groups were resolved from pre-optimisation hydrogen parents (`esp_map.py:2234–2244`).

In the glycine case one of those parents changed.

**Why it matters now.** `esp_map.py:81` tells the user to "Use --optimize" for small anions, because the geometry axis is worth up to 6.5 kcal/mol there. Nitrite, thiocyanate and the enolate are the P3 fixture set.

**Correction.**
1. After the optimisation, re-run `verify_computed(spec_opt, RemoveHs(probe))` and `verify_mapping(spec_opt, probe, reference=False)`. A mismatch refuses. With `--force` it warns, and the record gets `structure_verified: false`, `validity.input: false` and `sites: null`.
2. Record both identities: `geometry.method = "dft_optimised"` together with the pre-optimisation sha1.
3. Keep `atom_maps` and resolve `--sites` against the *final* geometry. Set `h_source` to `"optimised"` only for atoms that actually moved.
4. For a kept pose (`file 3D coordinates used as-is`), refuse `--optimize`, or freeze the supplied atoms. Optimising a pose discards exactly what the pose was for.

**Validation.**
- The glycine command above must exit non-zero. With `--force`, its JSON must have `structure_verified: false` and `validity.input: false`.
- `[CH3:1][OH:2] --optimize --json` must keep `atom_map` [1, 2].
- A fast test can stand in for berny: monkeypatch `pyscf.geomopt.berny_solver.optimize` to return a molecule built from the neutral glycine coordinates in `docs/review-4/out/r8.json`.

---

## Findings

### N2 · `analysis_key` omits how the geometry was made, the largest measured axis on anions `[MEDIUM]` · METHODOLOGICAL CONCERN, confirmed by execution

`analysis_fields` (`esp_map.py:1695–1709`) hashes `xc`, `basis`, `solvent`, `solvent_esp`, `spacing`, `density_fit` and the surface definition. It does not hash the source of the coordinates: ETKDG+MMFF94s, file frozen with MMFF-relaxed hydrogens, `--atoms`, or DFT-optimised. The spec (§5) defines the key as the answer to "are these numbers methodologically comparable?".

The project's own error budget makes the geometry source the largest axis on small anions: **up to 6.5 kcal/mol**, which is above the +3.54 basis row. Basis *is* in the key. The glycine pair in N1 shares `analysis_key 62f86d4d6b` across a 24.8 / 29.4 kcal/mol difference in V_min / V_max [X]. Part of that is the tautomer change, but an anion MMFF-vs-DFT pair would also share a key at a 6.5 kcal/mol difference [I, from `bench/anion_geometry.jsonl`].

`geometry.sha1` does not cover this: two *different* molecules always have different hashes, and comparing different molecules under one method is what the key is for.

**Correction.**
- Add `geometry_method` to `analysis_fields`, normalised to one of `etkdg_mmff94s` / `file_frozen` / `file_frozen+mmff94s_h` / `atoms_string` / `dft_opt:<xc>/<basis>[/<solvent>]`.
- Alternatively, keep the key as it is and rename it `electronic_key`, so it does not claim methodological comparability.

**Validation.** r8 and r9 (`docs/review-4/out/`) must get different keys. Two `--smiles` runs of different molecules must still share one.

### N3 · An ungated `--atoms` record is `validity.input: true` `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

`validity` (`esp_map.py:1712–1716`) computes `input = not (gate and gate["forced"]) and verified is not False`. This is the optional-presence guard. On `--atoms` and `--molecule` there is no gate (`gate = None`) and no verification (`verified = None`), so the absence of both reads as a pass.

Two waters 8 Å apart, `--atoms` (r5) [X], give:
- `validity: {input: true, surface: true, sites: true}`;
- `gate: null`, `structure_verified: null`;
- `surface.components: 2`, `problems: []`;
- a site table `O0`, `O3`.

The gate would have refused this input as `MULTIPLE_FRAGMENTS`. The "2 separate closed surfaces" warning reaches stderr only; in the record, `components > 1` is not a problem (`esp_map.py:1107–1110`).

This matters for the contract, because spec §7 test 9 recommends "an explicit `--atoms` string" for the symmetric P3 fixtures. Those fixtures will come out `validity.input: true` with `chemical_smiles: null`, and Arrowfield's documented rule (`if not record["validity"]["sites"]: refuse`) will accept them.

**Correction.** Make `validity.input` tri-state, `null` when no gate ran, and add `validity.input_checked: bool`. Better: perceive the structure from `spec` (`DetermineBonds` with `--charge`) and run `check_structure` on it, so `--atoms` is gated too. Either way, set `validity.sites = false` when `components > 1` and nothing explains it.

**Validation.** On r5, `validity.input` must not be `true`, and the record must name the reason.

### N4 · The atom-level partition is unbudgeted and dominated by the rule `[MEDIUM]` · METHODOLOGICAL CONCERN, measured

The per-site budget (README, CLAUDE.md; `test_claims._site_rows`) measures **groups**: a heavy atom with its hydrogens. In that setting the two rules mostly exchange area *inside* a group, and the exchange cancels. The record also publishes an `atoms` level with both rules.

Mapped ethanol, B3LYP/STO-3G (r1) [X]:

| atom | area, nearest / radius-scaled (Å²) | V̄, nearest / radius-scaled (kcal/mol) |
| --- | --- | --- |
| C0, C1 | **0.00 / 11.32**, **0.00 / 8.00** | null / −3.27, null / +3.93 |
| O2 | 9.58 / 12.78 (+33 %) | −30.85 / −24.39 (**+6.46**) |
| H (×6) | 11.3–11.7 / 7.6–8.0 (−32 %) | e.g. H8 +16.75 / +23.71 (**+6.97**) |
| **groups** C0, C1, O2 | Δ −0.52, +0.65, −0.13 | Δ **+0.32, +0.04, −0.80** |

Under the primary rule, both sp³ carbons are "buried". Under the other rule they own 8–11 Å². The atom-level deltas are about 8× the group-level ones. The record shows both values, so nothing is hidden. But no budget row covers this level, and the only guidance in the spec (§1, "why groups are needed") is marked `[I]`.

**Correction.**
- Add an atom-level row to `bench/site_partition.py`: the same rows, before the hydrogens are folded into their parents.
- State in the spec and the README that atom values are diagnostic and that `reaction_sites` is the quotable level.
- Optionally mark `atoms` in the schema as `"level": "diagnostic"`.

**Validation.** A committed atom-level row, and a `test_claims` entry that recomputes it.

### N5 · The record has drifted from the contract P3 is meant to freeze `[MEDIUM]` · CONFIRMED STATIC CODE DEFECT (from r1.json)

`docs/per-site-spec.md` §6 against what `main` writes:

| spec | record | effect on a consumer |
| --- | --- | --- |
| `sites.budget: [{axis, …, source}]` (§5, §6), with rows selected by the run's settings | **absent** | Arrowfield's "display a `prose` row as unbacked" has nothing to read |
| `sites.partition = {primary, compared, unit, radii}` | `{primary, unit, radius_scaled}` (`esp_map.py:1905–1907`) | renamed key; a reader of `compared`/`radii` gets `KeyError` or `None` |
| `geometry_provenance = {heavy_atoms, hydrogens: {supplied, generated, generated_relaxed}, relaxation, caveat}` | `relaxation` absent; `hydrogens` holds only the keys that occur (`{"generated_relaxed": 6}`) | absent reads as zero only if the consumer guesses it does |
| `identity.mapped_smiles = "<computed, maps kept>"` | the *input* (`approved`), `esp_map.py:2376` | names a different thing than the spec |
| NPZ `tri_atom`, `tri_atom_compared` | `tri_atom_nearest_nucleus`, `tri_atom_radius_scaled` | renamed |
| record shape `espmap.record/1` (identity · geometry · calculation · validity · …) | no record-level schema or version; still `input/settings/results` | nothing for a consumer to branch on |
| `validity.sites` false "if the site block failed its own checks" | the block has no checks; `sites = block is not None` (`esp_map.py:2358`) | conservation and nullability are asserted only in tests |

The spec header still says P1, P2 and `--sites` are "in review", although all three are merged.

**Correction.** Pick one source of truth before P3. Either update §6 to the code as it is, or change the code, then add `schema: "espmap.record/1"` at the top level. Add a fast test that builds the key tree of a `main` record (H₂O/STO-3G, or the stop-before-SCF hook with a fake surface) and compares it to a JSON schema committed next to the spec. The comparison must run in both directions, so that a key the spec declares and the code never writes fails too.

**Validation.** The schema test fails on today's record (`budget`, `compared`, `radii`, `relaxation`).

### N6 · `--sites` naming changes meaning with the input, silently `[MEDIUM]` · CONDITIONAL EDGE CASE, confirmed by execution

`resolve_site_groups` (`esp_map.py:1835–1837`) reads names as atom-map numbers when the input is mapped, and as **any** RDKit atom index otherwise. Three consequences [X, `e_fast.py` N3/N5]:

1. **The same file means different atoms.** On unmapped ethanol, `["1+H"]`, written for map 1 (C at index 0), resolves without error to **C index 1 and its hydrogens** (`[1, 6, 7]`). This happens if an adapter writes map-based groups and the maps are lost on the way. Map numbers are 1-based and indices 0-based, so the wrong atom is always a neighbour, and the result is always plausible.
2. **A hydrogen can be named.** `[5]` resolves to a lone H, and `"5+H"` to the same single atom. Spec §1 says heavy-atom indices only.
3. **"index" is not the heavy-atom index.** For an SDF that interleaves its hydrogens (`fx/ethanol_interleaved.sdf`, order C H H H C H H O H), O is index **7** and is labelled `O7`. Spec §1's claim that "the 0-based heavy-atom index in input order … RDKit index i is PySCF atom i for every heavy atom" holds only when the file lists heavy atoms first.

**Correction.**
- Require `{"by": "map" | "index", "groups": [...]}` in the `--sites` file.
- Refuse `by: map` on an unmapped input, and refuse an index that names a hydrogen.
- Write `by` into `reaction_sites[*].source`.
- Correct the spec wording to "RDKit atom index in file order, hydrogens included".

**Validation.** Parametrised fast tests for all three cases above.

### N7 · A molecule with no heavy atom crashes after the SCF `[LOW]` · CONFIRMED EXECUTION FAILURE

`site_parents` (`esp_map.py:1797–1803`) takes `argmin` over an empty heavy-atom list. `sites_block` runs on every run, with or without `--json` (`esp_map.py:2355`).

`--atoms "H 0 0 0; H 0 0 0.74"` and `--smiles "[H][H]"` (sto-3g) both write the PNG, then die with `ValueError: attempt to get argmin of an empty sequence`, with no JSON written [X]. H₂ passes the gate.

**Correction.** Return `None` parents when there are no heavy atoms. Then either write `sites: null` with a reason, or give each H its own default group.

**Validation.** The H₂ run writes a JSON.

### N8 · The NPZ carries the partition but not the validity `[LOW]` · CONFIRMED EXECUTION FAILURE

`--margin 2.0 --no-strict` on water (r4) writes a JSON with `validity.surface: false` and `sites: null` [X]. The NPZ from the same run carries `tri_atom_*` arrays and no `validity`, `problems`, gate, geometry sha1 or atom maps (`esp_map.py:2383–2398`).

Spec §6 offers the NPZ so Arrowfield can "redo a partition … without another SCF". That route bypasses "no site block on an invalid record".

**Correction.** Write `validity` (as JSON text), `surface_problems`, `geometry_sha1` and `atom_maps` into the NPZ. Omit `tri_atom_*` when `validity.surface` is false.

**Validation.** On r4's NPZ, `validity` is present and false.

### N9 · Basis names are not normalised, so one calculation gets three keys `[LOW]` · CONFIRMED EXECUTION FAILURE

PySCF treats `def2-tzvp`, `def2tzvp` and `DEF2_TZVP` as one basis (43 nao on water for each) [X]. `analysis_fields` lowercases but keeps the punctuation, which gives keys `055431fc19`, `be7935e2ac` and `af290272b0` [X].

The same gap disables the double-zeta warning: `"def2svp" in SMALL_BASIS` is False (`esp_map.py:165, 960`), so `--basis def2svp` runs without the +3.5 kcal/mol warning [X].

**Correction.** Normalise once, by stripping `-`, `_` and spaces and lowercasing (or use PySCF's own basis-name formatter), and use the result for both the key and `SMALL_BASIS`.

**Validation.** All three spellings give one key, and `def2svp` warns.

### N10 · Stale documentation and missing invariants `[LOW]` · CONFIRMED STATIC CODE DEFECT

- **Test counts.** "100 fast tests … 13 tests with real DFT" in `CLAUDE.md:25–26` and `README.md:305–306`. Collected: 241 fast, 18 slow. The README still says the tests come from "the two adversarial reviews"; there are three.
- **The notebook.**
  - "about 2.7 kcal/mol low on the molecules benchmarked" (`build_notebook.py:57`); the data range is +1.53 to +2.82.
  - "roughly halves the SCF time" (`:63`) under its def2-SVP default, where the data gives 20–28 %.
- **`docs/colormaps.md:144`.** Retracted +9.59 / 0.35 (H13).
- **`docs/roadmap-review.md:179`.** "+2.7 [M]".
- **`docs/prior-art.md:25–26`.** Global stats and per-atom surfaces listed as "planned (PR #24)"; both are merged.
- **`tests/test_claims.py` docstring.** Says the chlorobenzene budget is "pinned once below" with no data file. It is now recomputed from `bench/`.
- **Module docstring (`esp_map.py:87–88`).** Paracetamol and glycerol conformer numbers are quoted as measured with no `bench/` row. The README marks glycerol as "prose"; the docstring does not.
- **`CLAUDE.md`.**
  - The architecture table has no row for the sites stage (`partition_surface`, `site_report`, `verify_mapping`, `resolve_site_groups`).
  - Invariant 3 does not mention the `mapped` check, which spec "When this lands" said P1 would add.
  - "One geometry in, one record out" is not an invariant, although the spec said P0 would make it one.

**Correction.**
- Fix the items above.
- Extend `test_claims` to scan `notebooks/build_notebook.py`, `docs/colormaps.md` and `docs/roadmap-review.md` for the retracted figures. It already reads the notebook for "5.2".
- Derive the test counts in the docs from `pytest --collect-only`, or delete them.

### N11 · `V_area_sd` is cited as GIPF statistics but is not GIPF's σ² `[LOW]` · METHODOLOGICAL CONCERN

`docs/prior-art.md:25` and spec "Prior art" trace the global statistics to GIPF (Murray et al. 1994). GIPF's variance is σ²_tot = σ₊² + σ₋², each taken about its own sign's mean, together with the balance ν. `area_stats` reports one standard deviation about the overall mean (`esp_map.py:1776`). A consumer comparing it with Multiwfn's "overall variance" gets a different quantity under a borrowed name [S].

**Correction.** Either compute σ₊², σ₋², σ²_tot and ν, or reword the provenance to "area-weighted SD about the mean (not GIPF σ²_tot)".

### N12 · The per-site bench validates the implementation against its own prototype `[LOW]` · METHODOLOGICAL CONCERN

`bench/site_partition.py` says that `site_report` "reproduces" the bench prototype to 8.4e-10. Both were written in the same session with the same algorithm (centroid → argmin, nearest-heavy-atom parent, the same weighted quantile). This is the *checker shares the checked's assumption* pattern: it proves the refactor, not the partition.

The independent checks are the symmetry rows, and they pass (see "What holds"). Say so in the docstring, and do not cite the 8.4e-10 as validation.

### N13 · CLAUDE.md quotes the per-site budget without its protocol `[LOW]` · METHODOLOGICAL CONCERN

CLAUDE.md quotes the per-site budget without its protocol: "≤ 3.0 on the small anions … ≤ 3.95 on the C18 acetate … conformers move the median C18 site by 5.1". From the rows [X]:
- the anion rows are B3LYP/def2-TZVP **in PCM water** on DFT geometries;
- the C18 rows are def2-SVP+RI **in the gas phase**;
- the 5.1 includes the two folded conformers. Over the 8 extended conformers alone the median is **2.69** (max 5.16).

The README states some of this; CLAUDE.md, the file a future session reads first, states none of it.

**Correction.** Add the basis, the solvent and "folded included; 2.69 extended only" to the CLAUDE.md paragraph.

---

## What holds

- **The partition conserves area, and its extrema are the global extrema.** Conservation is ≤ 1e-12 in every one of the 25 bench rows. On the nitrite record the site areas sum to the total area to 1e-9, and min(site V_min) equals the global V_min (slow test) [X].
- **Symmetric sites agree within mesh noise, on geometries built symmetric** [X, `e_bench.py`]:
  - benzene, six C sites: means 2.378–2.566 kcal/mol, areas 18.50–18.93 Å² at 0.26 bohr, tightening to 2.487–2.539 at 0.16;
  - nitrite: O1/O2 means −147.455 / −147.442 at 0.26 bohr.

  This is the independent detector the spec asked for, and it is the right one.
- **The atom-map gate does not invent stereochemistry.** Mapped isopropanol, mapped 2,3-dimethylbutane and mapped citric acid give no `UNSPECIFIED_TETRAHEDRAL` finding [X]. The decision to strip maps for chemical identity and to check the correspondence separately in `verify_mapping` is the hard correct call. It is also why a correctly mapped isopropanol is not refused.
- **`verify_mapping` separates two failure modes.** Neighbour-by-index catches a reordered geometry. The reference-to-file comparison catches a consistent swap of graph-equivalent groups, which no graph test can see (`test_fast.py` swap tests) [S].
- **The H-count guard and bond perception treat PDB/XYZ honestly.** A hydrogen-free PDB is refused with a remedy instead of being completed into a plausible molecule [X].
- **The far-field detector still guards the most consequential bug.** Injecting the kernel-sign flip into the real `esp_at_points` gives `selftest() == 1` at +265.87 [X].
- **The dipole fix is right for the right reason.** About the centre of nuclear charge, acetate is 2.432360 D wherever it sits [X].
- **`test_claims.py` recomputes the numbers instead of string-matching them.** Every chlorobenzene budget row and every per-site row is computed from `bench/` and then required in the prose. That is the H8 correction done properly for the files it covers.
- **The notebook is a real artifact.** It regenerates byte-for-byte from its builder, holds no stored state, and every cell binds before use [X].

---

## Summary table

| # | finding | severity | evidence |
| --- | --- | --- | --- |
| N1 | `--optimize` computes an unverified structure (glycine zwitterion → neutral) recorded as verified and valid; maps dropped | **HIGH** | CONFIRMED EXECUTION FAILURE |
| N2 | `analysis_key` omits the geometry method (6.5 kcal/mol axis); r8/r9 share a key 25–29 kcal/mol apart | MEDIUM | METHODOLOGICAL, executed |
| N3 | ungated `--atoms` → `validity.input: true`; two waters give a valid site table | MEDIUM | CONFIRMED EXECUTION FAILURE |
| N4 | atom-level partition unbudgeted; ethanol C area 0 ↔ 11 Å², atom ΔV̄ up to 7 vs ≤ 0.8 for groups | MEDIUM | METHODOLOGICAL, measured (sto-3g) |
| N5 | record vs spec §6: `budget`, `partition.compared/radii`, `relaxation`, NPZ names, no record schema | MEDIUM | CONFIRMED STATIC |
| N6 | `--sites` meaning flips between maps and indices; H indices accepted; "heavy-atom index" false for interleaved files | MEDIUM | CONDITIONAL EDGE CASE, executed |
| N7 | H₂ / heavy-atom-free input crashes after the SCF | LOW | CONFIRMED EXECUTION FAILURE |
| N8 | NPZ carries the partition, not the validity | LOW | CONFIRMED EXECUTION FAILURE |
| N9 | basis spellings → three keys; `def2svp` skips the DZ warning | LOW | CONFIRMED EXECUTION FAILURE |
| N10 | stale counts, notebook/colormaps/roadmap figures; CLAUDE.md lacks the P1/P2 stage and invariants | LOW | CONFIRMED STATIC |
| N11 | `V_area_sd` attributed to GIPF but is not σ²_tot | LOW | METHODOLOGICAL |
| N12 | P2 validated against its own prototype; symmetry rows are the real check | LOW | METHODOLOGICAL |
| N13 | per-site budget quoted without basis/solvent; 5.1 is 2.69 without folded conformers | LOW | METHODOLOGICAL, recomputed |

## Verdict

**Blocks P3 (real fixtures, the Arrowfield contract):**
- **N1.** The documented remedy for the calibration anions is the one path whose output is never verified, and its record says otherwise.
- **N3.** The spec's recommended way to build symmetric fixtures yields `validity.input: true` with no gate.
- **N5.** P3 *is* the contract freeze, and the code and the spec disagree on six keys.

**Before any number is quoted from a record:** N2 (key scope) and N4 (atom level).

**Can wait:** N6–N13 and the PARTIAL residue of H8, H10, H12, H13 and H14.

**The structural lesson.** Review 3 found a check that compared the input with itself. Its fix compares the coordinates with the approval, but at one moment in the run. A later step that rewrites the coordinates escapes it. The general form is that **a verification is only as good as the last mutation after it**. Wherever `spec` can change (`--optimize` today; perhaps rotamer sampling or a protonation step tomorrow), the identity checks must run on the final `spec`, as close as possible to `gto.M`. The record must also derive `structure_verified` from that final check, not from an earlier one.

---

## Heavy reproductions requested

None of the findings above needs a heavy run to be established. These would size N1, N2 and N4, and close residue I could not run under the 1 GB / 4-thread limit:

1. **N1 on the P3 anions.**
   - Command: `esp_map.py --smiles "C=C[O-]" --optimize --basis def2-tzvp --solvent water --json enolate_opt.json`, and the same for `N(=O)[O-]` and `[S-]C#N`.
   - Expected cost: ~1–2 GB, 5–20 min each.
   - Decides: whether any calibration anion changes connectivity under optimisation, and the size of the MMFF→DFT shift that N2 says shares one `analysis_key`.
2. **N4 at production settings.**
   - Command: re-run `bench/site_partition.py` rows `grid|*|0.26`, plus three of the C18 conformers, with an atom-level (unfolded) statistics block added.
   - Expected cost: def2-TZVP anions ~1 GB each; C18 at def2-SVP+RI ~3 GB, ~5 min/conformer.
   - Decides: whether atom-level deltas stay ~8× the group-level ones beyond STO-3G ethanol.
3. **Full slow suite.**
   - Command: `pytest --runslow -m slow`, 18 tests.
   - Expected cost: ~3–4 min, peak memory unknown, def2-SVP chlorobenzene and bromobenzene.
   - Decides: the 12 slow tests I did not run (6 of 18 ran here).
4. **H14 residue.** Re-run one row of each of `error_budget.py`, `solvent_convention.py` and `site_partition.py` from a git checkout that WSL can read.
   - Decides: whether the recorded values reproduce at a known commit, so the 59 `"unknown"` rows can be tied to code.
