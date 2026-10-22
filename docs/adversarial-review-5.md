# Adversarial review 5 — espmap at master `1ef8e28`, and P3 at `8783322`

Reviewer 3 stance, 2026-10-23. This is the fifth pass and the last before `espmap.record/1` freezes. It covers the review-4 fixes N1–N13, the record contract (`docs/record-schema.json`), the P3 fixture generator and tests (draft PR #44, `feat/p3-fixtures`), and the six P3 fixture records. The records landed at `23dede9` while this review was under way. I first ran the generator's settings at STO-3G, then checked the real records once they existed (see "The P3 records at `23dede9`"). I wrote none of the code.

**Reviewed:**
- master `1ef8e28` in the worktree on branch `review/5`;
- P3 `8783322`, then `23dede9`, in a read-only checkout. `23dede9` adds only `fixtures/records/*`. Every P3 run used a copy in WSL home, never the checkout.

No repository file was modified. The only additions are this file and `docs/review-5/`. Nothing is committed.

**Environment:** WSL2 Ubuntu, `~/.venvs/espmap` with Python 3.12, numpy 2.5.4, pyscf 2.14.0, rdkit 2026.09.1 and pyberny. These are the versions in `bench/env.json`, and I tested no other library version. Every process ran alone at `OMP_NUM_THREADS=4`, peaking at 0.49 GB or less. The SCFs are STO-3G, except the anion HOMO and nitrite symmetry probes, which use def2-SVP. Numbers marked *sto-3g* are evidence about logic, not chemistry.

**Reproductions** are in `docs/review-5/`. Their `env.sh` holds absolute scratch paths, so adapt it before re-running.
- `e_fast.py`: no-SCF probes. Key collisions, budget selection, basis recognition, `coordinate_kind`, `--charge`/`--spin` through the gate.
- `e_records.py`: 13 CLI records, each checked against the schema and for agreement between its fields.
- `e_revert.py`: each review-4 fix reverted as a minimal mutation, with its tests run.
- `e_claims.py`: what `test_claims` actually pins.
- `e_scf.py`: anion HOMO signs at def2-SVP, and ν with one positive triangle.
- `e_forced.py`: a forced mapping failure, and the forced glycine `--optimize`.
- `run_p3.sh`: the P3 generator's settings at STO-3G, plus `tests/test_fixtures.py` on the results.
- `run_nitrite.sh`: nitrite O1/O2 symmetry against frame, spacing and basis.
- `e_p3records.py` (via `run_p3rec.sh`): the six real fixture records at `23dede9`, read-only, plus `tests/test_fixtures.py` on them and the installed-version check.

Outputs are in `docs/review-5/out/`.

Evidence labels: **[X]** verified by execution · **[S]** verified in source · **[I]** inferred.

---

## Summary verdict

The review-4 fixes are real and, with two exceptions in the claims test, guarded.
- I reverted each fix as a minimal mutation. A test fails for 9 of the 13 fixes, and for N12, N4 and N6 too once the selection is right (`e_revert.txt`).
- The headline N1 is closed end to end. The gas-phase glycine zwitterion with `--optimize --force` now records:
  - `structure_verified: false`;
  - `validity.input: false`;
  - `sites: null`;
  - `geometry_method: dft_opt:b3lyp/sto3g`;
  - `optimised_from_sha1` [X].
- Without `--force` it is refused [X, fast test].
- Thread count does not move a record: 4 vs 1 thread changes no number by more than 7.4e-10, and the meshes are identical [X].

**The lead finding is review 4's lead lesson again, one level up.** N2 put "how the geometry was made" into `analysis_key`, but only for the routes where espmap makes the geometry itself. Every `--file` geometry is `file_frozen`, whatever produced it.

The P3 fixtures are all `--file`, so the six records about to be frozen will carry `geometry_method: "file_frozen"` whether the coordinates are B3LYP/def2-TZVP/PCM (the four small molecules) or MMFF94s (the C18 pair).
- A DFT nitrite and an MMFF nitrite get **the same `analysis_key`** [X]. The MMFF geometry has N–O bonds of 1.424 Å against DFT's 1.255 Å, and the project's own budget puts this axis at up to 6.5 kcal/mol.
- `budget_rows` selects the "MMFF → DFT, 6.5 kcal/mol" row by `charge < 0` alone. It attaches that row to the DFT-geometry nitrite, thiocyanate and enolate fixtures, and to every `--optimize` run [X].
- The fixture test that checks the budget recomputes it with the same function, so it cannot fail.

The SDF already states the method (`geometry = "B3LYP/def2-TZVP, PCM water, …"`), and espmap ignores it.

The rest sits at the record boundary, and most of it is about what the record says the calculation was:
- **The gate never sees `--charge` or `--spin`.** `--smiles c1ccccc1 --charge 2` computes the benzene dication. The record says `computed_smiles: "c1ccccc1"`, `chemical_match: true` and validity all true [X].
- **An unbound gas-phase anion is warned on stderr only.** At def2-SVP, gas-phase acetate, nitrite and the enolate all have HOMO > 0, and their records are fully valid with site tables [X].
- **The frozen schema is wrong for a branch the code writes.** `results.surface.genus` is typed `int`, but every open surface writes `null` [X].
- **Default site labels are 0-based indices even on mapped input.** The P3 nitrite's "O1" is atom map 2, while the P3 C18's "C3-H" is map 3: two label conventions in one fixture set [X]. This is the N6 class at the label level.
- **`provenance.espmap_version` says 0.1.0 on all six fixture records**, while the code is 0.1.2-82 (R18). It comes from the stale editable install's metadata, not from the module that ran [X].

The real P3 records (`23dede9`) match the schema. The author's symmetry claims check out (benzene 0.0753, nitrite 0.0073), and 57/57 fixture tests pass [X]. The record-level findings above are all present in them.

Counts: **1 HIGH, 6 MEDIUM, 11 LOW** new findings. **N1–N13:** 9 FIXED (N1, N3, N4, N6, N7, N8, N9, N11, N12) and 4 PARTIAL (N2, N5, N10, N13). None regressed.

---

## Baseline

| check | result |
| --- | --- |
| `pytest -q` (fast + claims) | **314 passed**, 27 skipped (slow, plus P3 tests absent on master), 11.0 s [X] |
| `esp_map.py --selftest` | PASS, 11 colormaps oriented, 3.9 s, 0.18 GB [X] |
| determinism, acetate def2-SVP PCM, 4 / 4 / 1 threads | max numeric difference 3.7e-11 (4 vs 4) and 7.4e-10 (4 vs 1); `faces` identical; `analysis_key` identical [X] (`e_det.txt`) |
| P3 small fixtures at STO-3G (`mk.argv(..., "sto-3g")`) + `tests/test_fixtures.py -k` the four small molecules | 4 records, 5–12 s, ≤ 0.48 GB each; **31 passed**, 4 skipped [X] |

I did not re-run the full slow suite (26/26 on `1ef8e28` is already measured). The slow tests I touched, all STO-3G, ran inside `e_revert.py`.

---

## Disposition of review 4's N1–N13

"Guarded" means that reverting the fix as a minimal mutation of today's code makes at least one test fail (`e_revert.py`, `out/e_revert.txt`). I ran each one; none of these verdicts comes from reading the test.

| # | finding | verdict | guarded? | evidence and residue |
| --- | --- | --- | --- | --- |
| N1 | `--optimize` computed an unverified structure | **FIXED** | yes. Removing the post-optimisation `check_identity` fails `test_optimisation_that_moves_a_proton_is_refused`. Dropping maps under `--optimize` fails the slow `test_optimize_keeps_atom_maps…` | End to end, forced glycine gives `structure_verified: false`, `validity.input: false`, `computed_smiles: null`, sites withheld [X]. The re-check runs on `--file` and `--atoms` too, because `approved` is now always set (`esp_map.py:2533`) [S]. Residue: when verification fails, the record holds `computed_smiles: null` and does not say *what* was computed (it was neutral glycine). R8 is the same class on the mapping side. |
| N2 | `analysis_key` omitted the geometry method | **PARTIAL** | yes. Removing `geometry_method` from `analysis_fields` fails a key test | Closed for SMILES, `--atoms`, `--molecule` and `--optimize`. **Open for every `--file`: DFT and MMFF geometries share `file_frozen` and one key (R1).** One calculation also gets several keys (R6). |
| N3 | ungated `--atoms` was `validity.input: true` | **FIXED** | yes. Reverting `validity` to the optional-presence form fails `test_validity_input`. Letting the `--atoms` gate keep only the largest fragment fails `test_two_molecules_in_atoms…` | DIMENSION_AMBIGUOUS is reachable for CO₂, N₃⁻, BF₃, NO₃⁻ and O₃ written flat in a 3D-header MOL [X]. It is not reachable through `--atoms`, which never calls `coordinate_kind` [X], so the remedy text ("non-zero z, or save as .xyz") is right for every input that can trip it. New in this class: the gate checks the parsed charge, not `--charge` (R2). |
| N4 | atom level unbudgeted | **FIXED** | yes (README 11.9 → 10.9 fails `test_every_per_site_row…`) | "Buried" is defined the same way in the test and the prose: zero nearest-nucleus area, heavy atoms only, 22 heavy atoms in C₂₀H₃₆O₂ [S]. |
| N5 | record drifted from spec | **PARTIAL** | yes, for the shape: renaming `partition.compared` fails the fast schema test, and disabling `site_checks` fails `test_a_site_block_that_fails…` | The schema is checked against **one** valid record. The open-surface branch violates it (R4). The verdict vocabulary in spec §6 is one the code never writes (R10). |
| N6 | `--sites` naming flipped | **FIXED** | yes (whole fast file: allowing an H index, or `by: map` on unmapped input, each fails a test). A `-k "site"` selection misses them, which is worth knowing when running them | Residue: the *default* labels are index-based on mapped input (R5). |
| N7 | H₂ crashed after the SCF | **FIXED** | yes | H₂ writes its record with `sites: null` and a reason [X]. Residue: its NPZ still carries `tri_atom_*` (R9). |
| N8 | NPZ lacked validity | **FIXED** | yes (writing `tri_atom_*` unconditionally fails the slow NPZ test) | Residue: the condition is `input and surface`, not `validity.sites` (R9). |
| N9 | basis spellings, three keys | **FIXED** | yes | Residue: `SMALL_BASIS` misses cc-pVDZ, 6-31G(d), 6-31+G*, def2-SVPD and pcseg-1, so none of them warns; 6-31G(d) and 6-31g* are one basis with two keys (R13). |
| N10 | stale docs | **PARTIAL** | yes for the listed figures (re-inserting "about 2.7 kcal/mol" fails `test_stale_figures…`) | New drift has appeared (R16). |
| N11 | `V_area_sd` mis-cited as GIPF | **FIXED** | yes (variance about the overall mean fails the GIPF test) | ν flips from null to 0.0 on one triangle (R14). The weighting is stated in the record (`gipf.weighting`) and in `prior-art.md` [S]. |
| N12 | bench validated against its own prototype | **FIXED** | yes (changing `2.378` in the bench docstring fails `test_symmetry_rows…`) | P3 turns this check into a fixture tolerance that is not calibrated (R7). |
| N13 | per-site budget quoted without protocol | **PARTIAL** | **no for 2 of 9 numbers.** CLAUDE.md "2.69 → 2.96" and "5.1 → 5.7" both pass the whole claims file | The prose is right today. The guard is a substring test that `−2.69` (chloroform) and `5.16` satisfy (R12). |

---

## Lead finding

### R1 · A file geometry's method is invisible to the key and the budget, and the P3 fixtures are all file geometries `[HIGH]` · CONFIRMED EXECUTION FAILURE

`geometry_method` (`esp_map.py:1796–1822`) returns `file_frozen` for every file whose coordinates were kept, whatever produced them (`:1814–1816`). `budget_rows` (`:2192–2200`) selects the anion-geometry row on `charge < 0` alone and has no geometry argument.

`e_fast.py` [X]:

| input | how the coordinates were made | `geometry_method` | `analysis_key` |
| --- | --- | --- | --- |
| P3 `fixtures/geometries/nitrite.sdf` | B3LYP/def2-TZVP, PCM water (the SDF's own `geometry` property) | `file_frozen` | `85d0b5cbef` |
| `fx/nitrite_mmff.xyz` | ETKDG + MMFF94s, N–O 1.424 Å (DFT: 1.255) | `file_frozen` | **`85d0b5cbef`** |

**On the real records at `23dede9`** [X, `e_p3records.txt`]:
- all six record `geometry_method: file_frozen`;
- the four small molecules share `analysis_key 85d0b5cbef…`. That is **the very key** my MMFF nitrite XYZ gets at the same settings (table above), so a consumer's force-field anion run is declared comparable with the DFT-geometry fixture;
- nitrite, thiocyanate and the enolate carry the `anion geometry` row.

The STO-3G dry run (`run_p3.sh`) gives the same: all four small fixtures record `geometry_method: file_frozen`. Nitrite, thiocyanate and the enolate each carry the budget row `anion geometry`: "geometry MMFF → DFT … 6.5 kcal/mol", although their geometry *is* the DFT one. `--optimize` runs carry the same row (`e_fast.txt`). The C18 pair is also `file_frozen` (real records, key `49feb9f5e5…`) [X]. Their coordinates are byte-for-byte the MMFF conformers that `bench/site_partition.py` computed through `from_smiles` (`etkdg_mmff94s`), so the fixture and the budget row it cites have *different* keys (R6).

`tests/test_fixtures.py:111` (`test_fixture_budget_is_the_one_its_settings_select`) asserts `rec["sites"]["budget"] == em.budget_rows(...)` with the same inputs. It recomputes the producer with the producer, so it cannot fail. This is the *checker shares the checked's assumption* pattern.

**Why HIGH.** Review 4 rated the in-program version of this MEDIUM (N2). The fixtures are the exemplars `espmap.record/1` is frozen with, and Arrowfield is written against them. A consumer that trusts `analysis_key` to mean "methodologically comparable" will compare a DFT-geometry anion with an MMFF-geometry one, on the axis the project measured as its largest for anions. The fixture is also labelled with a 6.5 kcal/mol uncertainty it does not have. Once frozen, the vocabulary of `geometry_method` cannot grow a distinction without a schema version.

**Correction.**
1. Give files an origin. Read a declared property: `espmap_geometry_method`, or the `geometry` tag that `bench/fixture_geometries.py` already writes. Normalise it to the same vocabulary (`dft:b3lyp/def2tzvp/water`, `mmff94s`). If the file declares nothing, write `file_frozen:undeclared` rather than plain `file_frozen`. The key then refuses to equate two undeclared files with two declared ones.
2. Make `budget_rows` take `geometry_method`:
   - attach `anion geometry` only when the method is a force field or undeclared;
   - attach the `basis` row (def2-SVP → def2-TZVPP) only to double-zeta runs. Today a def2-TZVP run carries it too [X: benzene dry run].
3. Replace the tautological fixture test with expected row ids written out per fixture. Nitrite, for example: `{"basis"?, "functional", "grid", "partition small anions", "grid mean"}`, with no `anion geometry` row.
4. Enumerate the `geometry_method` vocabulary in `docs/record-schema.json` (see R10).

**Validation.**
- The two nitrite files above get different `geometry_method` values and different keys.
- The P3 nitrite record carries no `anion geometry` row.
- `--optimize` drops it too.
- The rewritten fixture test fails if `budget_rows` is reverted.

---

## Findings

### R2 · The gate never sees `--charge` or `--spin`, and the record says the approved species was computed `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

`check_structure` runs on the parsed probe (`esp_map.py:2425`), whose charge is the drawing's formal charge. For SMILES and files, `--charge` overrides it later (`:2479–2484`) with a stderr WARNING only. `--spin` with no perceived radical produces no finding at all (`check_structure`, `rad == 0` branch).

`e_records.py` [X]:

| run | gate | record |
| --- | --- | --- |
| `--smiles c1ccccc1 --charge 2` | ACCEPTED, nothing flagged | `computed_smiles: "c1ccccc1"`, `identity.chemical_match: true`, **validity all true**, site table, `geometry.charge: 2`, V **+213 … +295** |
| `--smiles c1ccccc1 --spin 2` | ACCEPTED, "closed-shell" | triplet UKS, validity all true, V −10.5 … +14.4 |
| `--smiles CC(=O)[O-] --charge 0 --spin 1` | ACCEPTED (findings: TAUTOMERS, CHARGED −1) | reaches compute (`e_fast.txt`) |

Invariant 3 says "the approved structure is the computed structure". The charge is half of what the SCF receives (nuclei, charge, spin), and it is the half the gate did not check. The record contradicts itself: `chemical_smiles` carries formal charge 0 while `geometry.charge` is +2. Spec §6's rule that a consumer "never infers validity from stderr" is broken by design here.

`--atoms` and XYZ are not affected, because `_perceive_bonds` uses `--charge`.

**Correction.**
- When `--charge` differs from the structure's formal charge, either:
  - refuse with a BLOCK `CHARGE_OVERRIDE` naming both values (remedy: "draw the charge on the atom"); or
  - with `--force`, record `gate.findings += [BLOCK, CHARGE_OVERRIDE]`, `forced: true`, `validity.input: false`.
- Do the same for `--spin ≠ radical count`: today that is INFO, or nothing when the count is zero. It should be at least WARN, and written to the record.
- Note that `structure_verified` for `--atoms` is true by construction, because the approval was perceived from the same coordinates. Say so in spec §6 so that "verified" is not read as an independent check on that route.

**Validation.**
- `--smiles c1ccccc1 --charge 2` exits non-zero.
- With `--force`, the record has `validity.input: false` and a `CHARGE_OVERRIDE` finding.
- `--smiles O=O --spin 2` (a deliberate triplet) still runs, with the finding recorded.

### R3 · An unbound gas-phase anion is a stderr warning, not a record field `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

`unbound_anion_warning` (`esp_map.py:1113–1129`) is printed at `:2561–2563` and never written. The default solvent is vacuum. `e_scf.py homo` at **def2-SVP** [X]:

| anion | gas HOMO (eV) | PCM water HOMO (eV) |
| --- | --- | --- |
| acetate | **+0.17** | −5.37 |
| nitrite | **+0.40** | −5.60 |
| enolate | **+1.38** | −4.06 |
| thiocyanate | −0.12 | −5.43 |
| F⁻ / OH⁻ | +4.25 / +4.48 | −3.80 / −2.99 |

In every case where the warning fires, the record still says validity all true and carries a site table. STO-3G gas nitrite (`nitrite_xyz.json`, +7.15 eV) and F⁻ (`fluoride_gas.json`, +18.3 eV) are the same [X]. The code's own warning says such a map "tracks the basis set, not the molecule". The P3 anions are computed in PCM water and are not affected. A consumer's own gas-phase anion run is affected, and it looks exactly like a P3 record.

Whether this persists at def2-TZVP is a heavy reproduction (below). The defect is that the condition is not recorded, not its frequency.

**Correction.**
- Write `results.homo_eV` (and `lumo_eV`).
- Add the condition to a record-level `problems` list. Alternatively, add `validity.electronic: bool`, false for an anion with HOMO > 0 in the gas phase. If the field is added, `sites` follows it, as it does for `surface`.
- This is a record change, so it belongs before the freeze.

**Validation.** `--smiles CC(=O)[O-] --basis def2-svp --json` gives `validity.electronic: false`, or a problem string naming the HOMO. With `--solvent water` it gives true.

### R4 · The frozen schema types `results.surface.genus` as non-null; every open surface writes null `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

`mesh_topology` sets `genus = None` when the mesh is not closed (`esp_map.py:1165`). `docs/record-schema.json:79` declares `"genus": "int"`.

`--smiles O --margin 2.0 --no-strict` (`water_open.json`) [X] gives:

```
schema problems: ['record.results.surface.genus: null, expected int']
```

The slow schema test runs one valid record (mapped CH₃Cl), so it never exercises the branch. I ran 12 other branches through `tree_problems`:
- UKS;
- forced BLOCK;
- `--file` XYZ;
- PCM with the reaction field;
- caller sites;
- an element without a Bondi radius;
- `--optimize`;
- H₂;
- an unbound anion;
- `--charge` and `--spin` overrides.

All 12 match the schema [X] (`e_records.txt`). The open surface is the only violation I found.

**Correction.**
- Change it to `"genus": "int?"` and state why: null means "not closed", not "unknown".
- Add the open-surface record to the schema test. It is cheap: water at STO-3G takes 2 s. So are the forced and UKS records.

**Validation.** `tree_problems` is empty for `water_open.json`, and the parametrised schema test covers at least the open, forced and UKS branches.

### R5 · Default site labels are 0-based indices on mapped input; P3 mixes two label conventions `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

Default reaction-site labels are `f"{symbol}{index}"` (`esp_map.py:2069`) whether or not the input carries atom maps.

In the committed P3 nitrite record (`23dede9`) and in my dry run [X], atom index 1 has `atom_map: 2`, and its site is labelled `"O1"` (`e_p3records.txt`, `p3dry_nitrite_sto3g.json`). The enolate's `"O2"` is map 3, and benzene's `"C0"` is map 1. The C18 fixtures, by contrast, label by map: `"C3-H"` is map 3.

A reader of the fixture set meets "O1" (index 1, map 2) next to "C3-H" (map 3). This is exactly the 0-based/1-based ambiguity that review 4 N6 closed for `--sites` *input*: a lost map always lands on a plausible neighbour. It is now frozen into the output labels instead. `tests/test_fixtures.py:145` hard-codes `m["O1"]`, `m["O2"]`.

**Correction (decide before the freeze).** Either:
- label default sites by map when the input is mapped (`O2`, `O3` for nitrite) and by index otherwise; or
- make the label unambiguous in both cases (`O#1` for an index, `O:2` for a map, SMILES-style).

Record which scheme was used in `reaction_sites[*].source` (`default:map` / `default:index`), as `caller:*` already does.

**Validation.** On the mapped nitrite, the label's number equals `atom_map` for every default site. The fixture test reads labels from `atoms[*].atom_map`, not from literals.

### R6 · One calculation gets several `analysis_key`s `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

The key hashes *labels*, not what the labels resolve to [X, `e_fast.txt`].

| one calculation | spellings or routes | keys |
| --- | --- | --- |
| B3LYP (VWN-RPA) | `b3lyp`, `b3lypg`, `HYB_GGA_XC_B3LYP`: one functional and the same energy to 1e-8 Ha | **3** (`3f69c228…`, `fd65a6e2…`, `cb05c012…`) |
| idealised chlorobenzene | `--molecule chlorobenzene` vs the identical atom string through `--atoms` | 2 (`idealised_builtin` / `atoms_string`) |
| any SMILES record | re-run through its own `geometry` with `--atoms`, the H14 reproduction route that `test_json_geometry_reproduces_the_run_through_atoms` exercises | 2 (`etkdg_mmff94s` / `atoms_string`). That test compares V and sha1, never the key |
| the C18 conformers | P3 `--file` vs `bench/site_partition.py` `from_smiles`, the same coordinates | 2 (`file_frozen` / `etkdg_mmff94s`) [I] |

The xc case also runs the other way [I]: `b3lyp` resolves through PySCF's `B3LYP_WITH_VWN5` configuration, which is unset here. A site that sets it computes a different functional (`b3lyp5`, 0.037 Ha away on water) under the same key. `provenance.pyscf` is recorded, but it is not in the key.

The second and third rows are what `geometry_method` was asked to do, which is to name the *route*. The table shows that a route is not a method.

**Correction.**
- Hash the resolved functional, e.g. `dft.libxc.parse_xc(xc)` coefficients plus `hybrid_coeff`, or a canonical name from them, rather than `xc.lower()`.
- For geometry, name the method, not the channel:
  - `--molecule` → `idealised`;
  - `--atoms` → `undeclared` (as in R1);
  - a run reproduced from a record → the original record's method, carried by an optional `--geometry-method` flag that the reproduction route sets.

**Validation.** `b3lyp`, `b3lypg` and `HYB_GGA_XC_B3LYP` give one key and `b3lyp5` another. The H14 round-trip test also asserts key equality.

### R7 · P3's symmetry tolerances are not calibrated against mesh noise `[LOW]` · CONDITIONAL EDGE CASE, measured

*Downgraded from MEDIUM once the real records landed.* At the fixture settings, the committed nitrite gives |O1 − O2| = **0.0073** and benzene a C-site range of **0.0753** [X]. Both are far inside their tolerances, so the test passes with margin today. What follows is about why the tolerances are what they are, not about a failure.

`tests/test_fixtures.py:145` asserts `|O1 − O2| < 0.1` kcal/mol, citing one bench value (0.013: def2-TZVP, PCM water, 0.26 bohr, the planar frame). `run_nitrite.sh` computes the same quantity at the settings I am allowed to run. These are site means, nearest-nucleus rule, PCM water [X]:

| basis | frame | 0.26 bohr | 0.16 bohr |
| --- | --- | --- | --- |
| STO-3G | planar (bench) | 0.070 | 0.097 |
| STO-3G | tilted (P3 input) | 0.084 | 0.012 |
| def2-SVP | planar | 0.083 | **0.119** |
| def2-SVP | tilted | 0.058 | 0.026 |

The geometry is symmetric to 2e-4 Å (`e_fast.txt`), so the differences are mesh noise. They depend on how the molecule sits in the grid, and they are not monotone in spacing. Eight measurements span 0.012–0.119, and one exceeds the tolerance. The bench 0.013 is one draw from this distribution, and the tilt itself changes the draw: the tilted and planar V_min differ by 0.25 kcal/mol at def2-SVP and 0.26 bohr. The tolerance is also tighter than the project's own "grid mean" budget row, 0.15.

At def2-TZVP the two draws I can see (0.013 bench and 0.007 fixture) are much smaller than my STO-3G and def2-SVP spread. The real noise at the fixture setting may therefore be small. But two draws are not a distribution, and a re-run in another frame or at another spacing could land anywhere in the range above. Benzene's 0.5 against draws of 0.19 and 0.075 has more room, but it rests on the same few draws.

**Correction.** Calibrate both tolerances from a noise distribution: the same geometry under about 10 random rigid rotations at the fixture's own basis and spacing. Use its maximum, or its 99th percentile, plus a margin. Record the distribution in `bench/` and cite it in the test. This is the "a resampling scheme must reproduce the observed statistic" check, applied to the one independent validation of the partition.

**Validation.** A committed `bench` row with the rotation spread, and the test tolerance read from it.

### R8 · A forced atom-order failure loses its reason `[LOW]` · CONFIRMED EXECUTION FAILURE (fault-injected)

The run was CCO with `verify_mapping` fault-injected to fail, plus `--force` (`e_forced_mapping.txt`) [X]:
- `validity.input: false` (correct);
- `structure_verified: true`, `gate.forced: false`;
- `sites_reason: "input not valid (forced past the gate, or structure not verified)"`. Neither cause in that message is what happened.

With `sites: null`, `identity.mapped_match` disappears, and no other field records the mapping result. The record is internally inconsistent: valid by every field it shows, invalid by the summary.

The mapping failure is hard to reach from the CLI today (berny does not reorder atoms, and `from_file` keeps the order), so this is LOW. It is still the field P3's consumers branch on.

**Correction.**
- Add `input.mapping_verified: bool?` next to `structure_verified`.
- Make `sites_reason` name the failing check (`gate forced` / `structure not verified` / `atom order not verified`).
- When `verify_computed` fails, also write the perceived SMILES, so the forced glycine record says it computed `NCC(=O)O`.

**Validation.** The injected run's record has `mapping_verified: false`, and its reason names the atom order.

### R9 · The NPZ partition is keyed on `input and surface`, not on `validity.sites`; absence has three encodings `[LOW]` · CONFIRMED EXECUTION FAILURE

`esp_map.py:2689` writes `tri_atom_*` when the input and surface are valid. `sites_block` also withholds a block for "no heavy atom" and for a failed `site_checks`. H₂ therefore writes an NPZ with both `tri_atom_*` arrays next to a JSON whose `validity.sites` is false (`e_records.txt`, `h2`) [X]. This is the N8 principle, partition only where the record has sites, applied to two of its three conditions.

Absence is also encoded three ways for one consumer:
- `null` in JSON (`partition.compared`, `atom_map`);
- `-1` in `tri_atom_radius_scaled` when a rule is unavailable ([AlH₄]⁻: every entry −1) [X];
- `0` in the NPZ `atom_maps`.

**Correction.**
- Key the NPZ arrays on `valid["sites"]`.
- Pick one absence convention before the freeze, and state it in the schema's `npz` comment. Omitting the array is the cleanest, since the schema already lists which arrays a valid run has.

**Validation.** The H₂ NPZ has no `tri_atom_*`. The [AlH₄]⁻ NPZ has no `tri_atom_radius_scaled`, or the schema documents −1.

### R10 · The record's vocabularies are typed from memory: spec §6 and a fast test use `"ACCEPT"`, which the code never writes `[LOW]` · CONFIRMED STATIC CODE DEFECT

`gate_report` writes `verdict` as the severity constant: `"INFO"`, `"WARN"` or `"BLOCK"` (`esp_map.py:1524`; every record in `e_records.txt` says `INFO`). `docs/per-site-spec.md:421` shows `"verdict": "ACCEPT"`, and so do `tests/test_fast.py:181,183`. A consumer written from the spec that branches on `verdict == "ACCEPT"` never matches.

The schema types every vocabulary as `"str"`:
- `verdict`;
- `reaction_sites[*].source`;
- `geometry_method`;
- `geometry_provenance.heavy_atoms`;
- `h_source`;
- `scale.kind`;
- `scale.source`.

So `record_tree` cannot catch this kind of drift.

**Correction.**
- Enumerate each vocabulary in `record-schema.json`, e.g. `"enum:INFO|WARN|BLOCK"`, and make `record_tree` check it.
- Fix the spec example and the test fixtures to values the code writes.
- Consider whether `"INFO"` is the verdict a consumer should read for "accepted"; renaming it after the freeze costs a version.

**Validation.** `record_tree` fails on a record whose verdict is `"ACCEPT"`.

### R11 · `coordinate_kind` discards real planar geometries with two near-equal bond types `[LOW]` · CONFIRMED EXECUTION FAILURE

The drawing test (`esp_map.py:310–315`) calls a z = 0 3D file a drawing when every bond is within 3 % of the mean and there are at least two element pairs.

C₆F₆ at DFT-like lengths (C–C 1.394, C–F 1.335 Å; maximum deviation 2.2 %) and C₂F₄ (1.327 / 1.319) are both read as `drawing`. Their coordinates are **discarded** and an ETKDG conformer is computed instead (`e_fast.txt`) [X]. The record shows it (`settings.geometry` says DISCARDED, and `geometry_method` is `etkdg_mmff94s`), so this is LOW. A supplied perfluoroaromatic pose is replaced silently in a pipeline.

The single-pair set CO₂, N₃⁻, BF₃, NO₃⁻ and O₃ gives DIMENSION_AMBIGUOUS [X]. That answers the brief's N3 question: the BLOCK is reachable for any D∞h or D3h single-pair molecule written flat. The remedy text is right for files, and `--atoms` is never affected.

**Correction.** For MOL/SDF with a 3D header, treat a file as a drawing only when bond lengths are uniform *and* equal to a drawing program's standard length (≈ 1.5 Å, or 0.825 Å in ChemDraw's units), or *and* the 2D header flag is set. Otherwise refuse as ambiguous rather than discard.

**Validation.** `fx/C6F6.mol` is read as `geometry`, and the review-3 H5 drawing fixtures are still read as `drawing`.

### R12 · `test_claims` substring checks do not pin two of the CLAUDE.md per-site numbers `[LOW]` · CONFIRMED EXECUTION FAILURE

`test_claude_md_quotes_the_per_site_budget` (`tests/test_claims.py:442`) checks `value in CLAUDE.md` over the whole file. `e_claims.py` replaces each value in the per-site paragraph and re-runs the check [X]:
- "5.1" (conformer C18) still passes, because `5.16` contains it;
- "2.69" (extended conformers) still passes, because the chloroform row's `−2.69` contains it.

The mutation run confirms it: CLAUDE.md "2.69 → 2.96" and "5.1 → 5.7" both pass all 66 claims tests (`e_revert.txt`). The same `in TEXT[name]` idiom is used at `:200` for the README and CLAUDE.md budget rows. I did not audit those collisions.

**Correction.** Match tokens, not substrings (`re.search(rf"(?<![\d.]){re.escape(v)}(?!\d)", para)`), and search the paragraph the number belongs to, not the whole file.

**Validation.** Both mutations above fail the test.

### R13 · Double-zeta recognition and basis aliases are incomplete (N9 residue) `[LOW]` · CONFIRMED EXECUTION FAILURE

`SMALL_BASIS` (`esp_map.py:173`) lists six names. cc-pVDZ, aug-cc-pVDZ, 6-31G(d), 6-31+G*, 6-31G**, def2-SVPD and pcseg-1 are all double zeta, and none warns [X]. `6-31G(d)` and `6-31g*` resolve to the same basis (18 nao on water) and get two keys [X].

**Correction.**
- Classify by zeta from the basis itself, e.g. count contracted valence s functions on carbon, or keep an explicit table with a test over PySCF's basis aliases.
- Canonicalise `(d)` → `*` and `(d,p)` → `**`.

**Validation.** Each name above warns, and the 6-31G(d)/6-31g* pair shares a key.

### R14 · ν jumps from null to 0.0 on one triangle `[LOW]` · CONFIRMED EXECUTION FAILURE

On a synthetic all-negative sphere, one triangle (5.8e-4 of the area) set to +0.01 kcal/mol changes `gipf` [X] (`e_nu.txt`):
- from `nu: null, sigma2_pos: null`;
- to `nu: 0.0, sigma2_pos: 0.0, V_pos_mean: 0.01`.

`sign_change_ratio` is then 0.0001, so the scale stays sequential, but the GIPF block reports a balance. This is the knife-edge class that review 3 H11 removed from the scale.

**Correction.** Below a minimum area fraction for either sign (say the same 5 %, or one stated in the record), report that sign's statistics and ν as null with a reason. Always record `area_pos_A2`/`area_neg_A2` next to ν. They are already in `surface_stats`; a pointer in the `gipf` block is enough.

**Validation.** The one-triangle case gives `nu: null`.

### R15 · Prior-art gaps `[LOW]` · METHODOLOGICAL CONCERN [S]

`docs/prior-art.md` and `references.bib` do not cite:
- **Bondi 1964 and Mantina et al. 2009.** The radius-scaled partition rule depends on these radii (`BONDI`, `esp_map.py:1857`; spec §2 names both), and the record writes `"radii": "Bondi; Mantina 2009 for B"`.
- **xyz2mol (Kim & Kim 2015).** RDKit's `rdDetermineBonds` implements it, and every XYZ, PDB and `--atoms` input, and every `verify_computed` call, goes through it. Citing RDKit alone does not name the method.
- **The positive-HOMO criterion for an unbound anion** (`unbound_anion_warning`). It is established practice in DFT anion work, and it is presented without a citation [I: I did not run a literature search].

Nothing on the "specific to espmap" list is contradicted by the code. The radius-scaled rule is a multiplicatively weighted Voronoi assignment. Before "two partition rules computed together" is claimed, someone should check whether Multiwfn's atomic surface decomposition already offers it [I].

**Correction.** Add the three references and a row for weighted-Voronoi surface partitioning.

### R16 · Documentation drift since review 4 `[LOW]` · CONFIRMED STATIC CODE DEFECT

- `CLAUDE.md:13`: "~1700 lines"; it is 2761.
- `CLAUDE.md:213`: "`-1.md`, `-2.md` and `-3.md` record every defect found so far"; review 4 exists.
- `CLAUDE.md:243`: "4 minutes" for the slow suite; it measured 129 s on `1ef8e28`.
- `CLAUDE.md:192` and `README.md:289`: `--optimize` "exercised on one molecule". It has now been run on nitrite, SCN⁻ and the enolate at def2-TZVP/water (brief §A), plus methanol and glycine in tests.
- Spec §6 (`per-site-spec.md:432`) lists when `validity.sites` is false, but omits "no heavy atom" (N7).
- "Extended" means genus 0 in `test_claims._site_rows`. The 8/2 split holds in the data [X: 8 genus-0 rows and 2 at genus 2 and 3], but the README never states the definition.

### R17 · P3 fixture tests do not pin provenance or values `[LOW]` · CONDITIONAL EDGE CASE

The committed records are clean on this point. All six say `code: v0.1.2-82-g878332291704` (that is `8783322`, not dirty), with one pyscf, one rdkit, one numpy and 4 threads [X]. The concern is what stops the next regeneration from breaking that:
- `make_records.py:155` skips any fixture whose JSON exists. A partial re-run therefore mixes records from different code revisions, and no test asserts that the six share `provenance.code`, `pyscf` and `rdkit`.
- `code` is `"unknown"` whenever the run cannot read git. My dry run from a copy got `"unknown"` [X], and nothing would refuse such a record.
- No fixture test pins a computed number, so a regeneration passes whatever it produces. The determinism I measured (≤ 7.4e-10 across thread counts) and the author's reported Windows-vs-WSL agreement (≤ 1e-10) both say a value pin would be safe. I have not seen the second run's files; the claim is consistent with my measurement but unverified.

**Correction.**
- Assert in `test_fixtures.py` that `code` is not `"unknown"`, does not end in `-dirty`, and is equal across the six records.
- Pin V_min/V_max and one site mean per fixture to about 1e-6.
- Make the generator refuse to mix revisions: compare `provenance.code` of the existing JSONs before skipping.

### R18 · `provenance.espmap_version` is the installed metadata, not the code that ran: every fixture says 0.1.0 `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

This is the Publishing session's candidate finding, verified.

`code_identity` (`esp_map.py:1717–1745`) reads `importlib.metadata.version("espmap")`. All six P3 records say `espmap_version: "0.1.0"`, while `pyproject.toml` (master and P3) says 0.1.2 and `code` says `v0.1.2-82-g8783322…` [X]. The WSL venv holds `espmap-0.1.0.dist-info`, an editable install of `/mnt/c/Users/user/Documents/espmap` made at 0.1.0 (`direct_url.json`) [X].

**It does not even describe the same tree.** In my run, `esp_map` was imported from a copy in `~/p3copy`, while the version came from the editable install of the Windows checkout [X]. The field names whatever `espmap` distribution happens to be installed in the interpreter, regardless of which file ran. The docstring says so ("for a checkout `code` is the identity to trust"), but the record does not, and `espmap_version` is the field a consumer reads first. 0.1.0 predates reviews 3 and 4: a reader triaging by version would treat the fixtures as pre-fix output.

**Before the freeze? Yes.** It is a field-semantics question, cheap now and a version bump later. The fix is about five lines.

**Correction.**
1. Single-source the version in the module: `__version__ = "0.1.2"` in `esp_map.py`, with `pyproject.toml` `dynamic = ["version"]` and `[tool.setuptools.dynamic] version = {attr = "esp_map.__version__"}`. Record `espmap_version = esp_map.__version__`, read from the module that actually ran.
2. If the installed metadata is wanted at all, keep it under an honest name, `installed_distribution_version`, nullable.
3. Add a fast test: `em.__version__` equals the version `tomllib` reads from `pyproject.toml`.
4. Add to `test_fixtures.py`: `espmap_version` is consistent with `code` (the `git describe` tag prefix equals `v` + version).

**Validation.**
- With the stale 0.1.0 editable install still in the venv, a record made from this tree says `espmap_version: "0.1.2"`.
- Bumping `pyproject.toml` without `__version__` (impossible once single-sourced), or the reverse, fails the fast test.
- The six fixture records need regenerating, or a documented one-line patch of `provenance.espmap_version`. Regenerating is cleaner, and the records are deterministic.

---

## The P3 records at `23dede9`

I read the six committed records (`e_p3records.py`, `out/e_p3records.txt`), with no SCF.

**What I verified** [X]:

| author's claim | verified |
| --- | --- |
| 57/57 fixture tests pass | **57 passed** |
| benzene C sites agree to 0.075 | 0.0753 (nearest nucleus; radius-scaled also 0.0753) |
| nitrite O1/O2 to 0.007 | 0.0073 (radius-scaled 0.0027) |
| made from a clean clone at `8783322` | `code: v0.1.2-82-g878332291704` on all six, not dirty |
| Windows-worktree re-run agrees to ≤ 1e-10 with identical keys | not checkable (files not supplied); consistent with my 7.4e-10 thread test |

**Every record is schema-clean.** `tree_problems` is empty, the NPZ file lists equal `schema["npz"]`, keys recompute, and the validity is all true. The C18 conformer is folded to genus 3, the extended one is genus 0, and both have 15 buried heavy atoms (the README says 14–15). The tests cover all of this, so it confirms them rather than adding to them.

**What the records freeze that the tests do not see:**
- **R1 on the real data.** All six are `file_frozen`. The four small molecules share key `85d0b5cbef`, which an MMFF nitrite XYZ also gets. Nitrite, SCN⁻ and the enolate carry the 6.5 kcal/mol `anion geometry` row on their DFT geometries.
- **R5 on the real data.** The nitrite sites `O1`/`O2` are maps 2/3, while the C18 sites `C1-H`…`C18-H` are maps 1–18.
- **The C18 sites overlap, and the record does not say so.** Atom index 2 (carbonyl O, map 32) is in both `carbonyl O` and `acetyl`. The site areas sum to **1.041** (extended) and **1.027** (folded) of the surface. A consumer that sums `reaction_sites` areas, or computes "fraction of surface", is off by 3–4 %, and no field signals it (Before the freeze, item 4).
- **R18.** All six say `espmap_version: "0.1.0"`.

The small-anion records have no unbound-HOMO issue (R3), because they are in PCM water; the C18 ester is neutral.

---

## What holds

- **N1 is closed where it shipped.** The re-verification runs after every replacement of `spec`. It runs for SMILES, `--file` and `--atoms` alike, and forced or refused, the record agrees with itself on the glycine case [X]. "A verification is only as good as the last mutation after it" is now in the code, not just in review 4.
- **The schema test is two-directional and real.** Renaming one key fails it. Twelve record branches I generated beyond its one fixture pass it. The one failure (R4) is a nullability annotation, not a shape drift [X].
- **`site_checks` is a real gate.** Disabling it fails a test. Conservation holds to 1e-9 on every record I made, including caller sites and UKS [X].
- **The partition and the record are deterministic across thread counts**, to 7.4e-10 [X]. That is the precondition for freezing fixtures, and it was not a given with PCM and BLAS threading.
- **The `--sites` naming fix is complete for input.** `by` is mandatory, hydrogens are refused, and `by: map` on unmapped input is refused, each guarded by a test [X].
- **The hard decision in P3 is right.** It computes each fixture in the condition its geometry was made in, rather than one condition for all six. The alternative would compute a PCM-optimised anion in the gas phase (unbound, R3) or the C18 in a protocol no budget row covers. The cost is that the two halves of the fixture set are not comparable, and `analysis_key` already says so. It will say so correctly once R1 is fixed.
- **`chain_maps` fails loudly** on a branched chain or a chain that is not C18 (`assert`, `make_records.py`) [S]. The tilt is a rigid rotation, checked by an inter-atomic-distance test at 1e-3 Å, and it leaves `verify_mapping`'s reference check comparing the file with itself [S, X: dry run `mapped_match: true`].
- **The kernel-sign detector still guards invariant 4.** `--selftest` passes. Review 4's injection result stands, and nothing since has touched `esp_at_points`.

---

## Summary table

| # | finding | severity | evidence |
| --- | --- | --- | --- |
| R1 | file geometries all `file_frozen`: DFT and MMFF nitrite share a key; the 6.5 kcal/mol anion-geometry row is attached to DFT and `--optimize` geometries; the P3 budget test is tautological | **HIGH** | CONFIRMED EXECUTION FAILURE |
| R2 | gate ignores `--charge`/`--spin`; a benzene dication is recorded as verified `c1ccccc1`, validity all true | MEDIUM | CONFIRMED EXECUTION FAILURE |
| R3 | unbound gas-phase anion (acetate, nitrite, enolate HOMO > 0 at def2-SVP) not in the record | MEDIUM | CONFIRMED EXECUTION FAILURE |
| R4 | schema types `surface.genus` int; an open surface writes null | MEDIUM | CONFIRMED EXECUTION FAILURE |
| R5 | default site labels are 0-based indices on mapped input; P3 mixes index and map labels | MEDIUM | CONFIRMED EXECUTION FAILURE |
| R6 | one calculation, several keys: xc aliases, `--molecule`/`--atoms`, the H14 reproduction route | MEDIUM | CONFIRMED EXECUTION FAILURE |
| R7 | P3 symmetry tolerances uncalibrated: noise 0.012–0.119 at STO-3G/def2-SVP; real fixture 0.007 | LOW | CONDITIONAL EDGE CASE, measured |
| R8 | forced atom-order failure: reason lost, `sites_reason` misattributes | LOW | CONFIRMED (fault-injected) |
| R9 | NPZ partition on H₂; −1 / 0 / null absence encodings | LOW | CONFIRMED EXECUTION FAILURE |
| R10 | `"ACCEPT"` in spec and tests, never written; vocabularies not in the schema | LOW | CONFIRMED STATIC |
| R11 | C₆F₆/C₂F₄ real geometries discarded as drawings; DIMENSION_AMBIGUOUS set | LOW | CONFIRMED EXECUTION FAILURE |
| R12 | `test_claims` substring check: 5.1 and 2.69 unpinned | LOW | CONFIRMED EXECUTION FAILURE |
| R13 | cc-pVDZ, 6-31G(d) etc. not recognised as DZ; alias keys | LOW | CONFIRMED EXECUTION FAILURE |
| R14 | ν null → 0.0 on one triangle | LOW | CONFIRMED EXECUTION FAILURE |
| R15 | Bondi, Mantina, xyz2mol uncited | LOW | METHODOLOGICAL |
| R16 | CLAUDE.md, README and spec drift | LOW | CONFIRMED STATIC |
| R17 | P3 tests pin no provenance or values; generator can mix revisions | LOW | CONDITIONAL EDGE CASE |
| R18 | `espmap_version` from installed metadata: all six fixtures say 0.1.0, code is 0.1.2-82 | MEDIUM | CONFIRMED EXECUTION FAILURE |

## Verdict

**Blocks the freeze** (each changes a key, a vocabulary or a field of `espmap.record/1`): R1, R4, R5, R10 and R18, plus the record half of R2 and R3. R1, R5 and R18 are present in the six committed fixture records, which will need regenerating once they are fixed; the records are deterministic, so that is cheap.

**Before the P3 PR merges:** R17 (pin provenance and values) and the site-overlap field. R7 can follow with a calibration row.

**Can wait:** R6 (key normalisation can be tightened inside version 1, because a key is recomputable from `key_fields`), and R7, R8, R9 and R11–R16.

**The structural lesson.** Review 4's lesson was that a verification is only as good as the last mutation after it. Review 5's is its mirror: **a field is only as good as the inputs that can reach it.** `geometry_method` was fixed on the paths espmap controls and left blind on the path its own fixtures take. The gate approves a species that `--charge` then changes. The schema was validated on the one branch the test ran. In each case the fix was tested on the input that motivated it. The next class to hunt is any field whose value is set by one route and read as if every route set it.

---

## Before the freeze

Record changes to make before `espmap.record/1` is frozen, ranked by the cost of making them later.

1. **Geometry origin (R1).** Read a declared method from the file (SDF property), write `file_frozen:undeclared` otherwise, and select the anion-geometry and basis budget rows by `geometry_method` and zeta. Without this, the P3 records freeze a key that equates DFT and MMFF anion geometries, and a budget that claims 6.5 kcal/mol of uncertainty for the DFT ones.
2. **Electronic validity in the record (R3, R2).** Add `results.homo_eV`/`lumo_eV` and a record-level `problems` list, or `validity.electronic`. Record a `CHARGE_OVERRIDE`/`SPIN_OVERRIDE` finding in `gate.findings` with `validity.input: false`. These are new fields: free now, a version bump later.
3. **Schema correctness and vocabularies (R4, R10).**
   - `results.surface.genus: "int?"`.
   - Enumerate `gate.verdict`, `reaction_sites[*].source`, `geometry_method`, `heavy_atoms`, `h_source` and `scale.kind`/`source`, and have `record_tree` check them.
   - Fix the `"ACCEPT"` example.
   - Extend the schema test to the open, forced and UKS branches (STO-3G, seconds each).
4. **Site label and site-set conventions (R5).**
   - Map-based default labels on mapped input (or an unambiguous syntax), with `source: default:map|default:index`.
   - While deciding, also fix the order of `reaction_sites` (caller order for caller sites, ascending index for defaults; today's behaviour, but unwritten).
   - Add a field saying whether the sites tile the surface (`reaction_sites_partition_surface: bool`, or a per-site `overlaps: [labels]`). The committed C18 records' `acetyl` and `carbonyl O` groups overlap on map 32, and their site areas sum to 1.041 and 1.027 of the surface [X]; no field says so. My own first probe for such a field matched `sha1_covers`, which is the substring trap again.
5. **Version identity (R18).** Single-source `__version__` in `esp_map.py`, record it as `espmap_version`, and test it against `pyproject.toml`. Then regenerate the six fixtures (they say 0.1.0 today).
6. **Absence encodings and mapping verification (R8, R9).** Add `input.mapping_verified: bool?` and a specific `sites_reason`. Key the NPZ partition on `validity.sites`, and use one absence convention (omit the array, not −1 or 0).
7. **Functional identity in the key (R6).** Hash the resolved functional, not the label. A key is recomputable from `key_fields`, so this can follow the freeze if it must, but `key_fields.analysis.xc` should already carry the canonical name.

---

## Heavy reproductions requested

None of the findings above needs a heavy run to be established. These would size R3 and R7 and close the P3 contract.

1. **R7, calibrate the nitrite and benzene symmetry tolerances.**
   - Command: tilted nitrite and benzene at the fixture settings, `--basis def2-tzvp --solvent water --spacing 0.26`, under 10 random rigid rotations each (`--atoms`, the `tilted()` helper with random axes). Record `|O1−O2|` and the C-site range.
   - Expected cost: nitrite ~0.6 GB and ~30 s per run; benzene ~1 GB and ~2 min per run; ~25 min in total.
   - Decides: whether `< 0.1` and `< 0.5` are inside the noise, and what tolerance the P3 tests should cite.
2. **R3 at production basis.**
   - Command: `run_scf` gas phase, B3LYP/def2-TZVP, on acetate, nitrite, the enolate and SCN⁻, printing the HOMO; also def2-TZVPP for acetate.
   - Expected cost: ≤ 1 GB, ~2–5 min total.
   - Decides: whether the default (gas, def2-TZVP) run of a common anion produces an unbound-HOMO record. That decides whether `validity.electronic` is a corner case or the usual case.
3. **R1 on the C18 pair.**
   - Command: the P3 `Z3Z13-18OAc_conf4` input through `--file`, and the same coordinates through `bench/site_partition.py`'s `from_smiles` route, def2-SVP + RI, gas.
   - Expected cost: ~3 GB, ~5 min each, so wait for the heavy slot.
   - Decides: that identical V comes with two `analysis_key`s, i.e. that the fixture and the budget row it cites are, by the key, incomparable.
