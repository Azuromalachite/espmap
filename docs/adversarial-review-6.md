# Adversarial review 6: espmap master `dc9d9ef` and PR #73 (`feat/ensemble-selection`, `d51ae9e`)

Reviewer 3 stance, 2026-10-23. This is the first review of v0.3 code.

**What was reviewed**
- Master `dc9d9ef`, in a detached worktree (`espmap-review6`).
- PR #73 `d51ae9e`, in a second detached worktree (`espmap-pr73`).

**What was changed:** nothing in the repository. The only additions are this file and `docs/review-6/`. Nothing is committed.

**Environment**
- WSL2 Ubuntu, `~/.venvs/espmap`, the library versions in `bench/env.json`.
- Every process ran at 4 threads and used under 0.6 GB.
- No SCF was larger than def2-SVP on 4 heavy atoms. The ensemble probes use RDKit only.

**Reproductions** are in `docs/review-6/`. `env.sh` holds absolute scratch paths, so adapt it before re-running. Prefixes:
- `e_*`: PR #73;
- `d_*`: the review-5 disposition;
- `a_*`: R11 and `--require-stereo`;
- `c_*`: ALIE, claims and prior art.

Outputs are in `docs/review-6/out/`, and the fixture SDFs in `docs/review-6/fx/`.

Evidence labels: **[X]** verified by running it · **[S]** verified by reading the source · **[I]** inferred.

## Disposition (added 2026-10-23, after #73 merged at 292a954)

The review below is unedited. This section records what happened to its findings.

**Fixed in #73.** The counterexamples were committed as strict xfails before each fix (c92a448, then 9dd0029), and the fixes followed (2887cd2, 64a5494).

| finding | outcome |
| --- | --- |
| S1 (HIGH) | Fixed, with one correction to the review (below). |
| S2 | Fixed. |
| S4 | Rg ranking fixed. The re-spread half is deferred to §8, with its strict xfail kept. |
| S5 | Fixed. |
| S6 | Fixed. This includes master's `--record N` numbering. |

- **S1.** `configuration_signs` reads, per conformer, three elements RDKit doesn't perceive: the helicity of a trans double bond in a ring of 8 or fewer atoms, a stereogenic cumulene's axis, and an atropisomeric biaryl axis.
  - Supplied members must match record 1's.
  - A generated ensemble from a 3D file keeps only conformers in the file's configuration (`n_other_configuration`).
  - `--require-stereo` leaves small-ring trans bonds open (`small_ring_trans_bonds`). Hostile cases E18–E21.
  - Arrowfield's review of the gate change found that E/Z had been read as ring-trans. E/Z is a CIP label, so `Cl/C1=C/CCCCCC1`, a *cis* ring, was refused. Cis/trans is now read relative to the ring atoms (64a5494).
- **S2.** Every supplied record is gated, and must share record 1's formal charge, radical count and heavy-atom maps. Hostile cases E22–E24.
- **S4.** Selection ranks by RDKit's mass-weighted Rg, as `bench/conformers.py` did, and reproduces the measured C18 members exactly. The re-spread half is deferred to the §8 measurements: the C18 lowest-energy conformer collides with an Rg pick, and the measured rows used the bench's top-up.
- **S5.** A molecule reaching `SYMMETRY_MATCH_CAP` (10 000) is refused, not truncated.
- **S6.** Unreadable records keep their file number (`RECORD_UNREADABLE`). Hostile case E25.

**One finding corrected.** The S1 table's tetra-ortho biaryl pair, `Cc1cccc(C)c1-c1c(C)cccc1[C@H](C)O`, is not a counterexample. Its 2,6-dimethylphenyl end is symmetric about the axis, so the two "atropisomers" are one compound, and accepting them is correct. The trans-cyclooctene and allene rows stand.

**Deferred**
- **S3** (contact threshold): waits on heavy reproduction 1 and the spec's §8 measurements. The constants stay provisional.
- **S7, S8** (R11: z-noise drawings, a 2D header with z ≠ 0): Arrowfield, each in its own PR, counterexamples first.
- **S9** (ALIE |ε| on unbound anions): its own PR.
- **S10–S26** (LOW), and the "Before `espmap.ensemble/1` is frozen" list: open. The list feeds step 4 (the runner and the record).
- **Heavy reproductions 1–3:** not run. Each needs the shared WSL slot.

---

## Summary verdict

**Review 5 is closed and guarded.**
- I reverted each fix with a minimal mutation (29 in all). Every mutation fails a test that passes on the unmutated tree (`out/d_revert.txt`).
- 14 of R1–R18 are FIXED and 4 are PARTIAL (R8, R10, R13, R16). None regressed.
- Five fixes left a residue of their own class (S19–S24), all LOW.

**PR #73 should not merge as is.** The lead finding is the one the ensemble spec was written to prevent: **a generated ensemble can average over diastereomers** (S1).
- `--require-stereo` and `verify_computed` share RDKit's stereo model. The model cannot see planar chirality (trans-cyclooctene) or a stereogenic axis (an allene, a biaryl) once the input is 3D.
- In those cases ETKDG picks a configuration per embedding, and every member is "verified" [X]:
  - `O[C@H]1CC/C=C/CCC1` with K = 10 gives 5 members whose ring helicities are `[+,−,+,−,−]`;
  - a 3D allene SDF gives `[+,−,−,+]`.
- The supplied route accepts the same diastereomer pairs, and a tetra-ortho biaryl pair too.

**The supplied route gates record 1 only.**
- One "ensemble" can hold three redox states of the same nuclei (ethyl cation, radical and anion; TEMPO⁺, TEMPO•, TEMPO⁻), and conflicting atom maps [X].
- Which state the ensemble takes depends on the record order: radical first is refused, cation first is accepted.

**The two shape rules don't reproduce the protocol they cite.**
- **Contact rule.** It calls 36 of 100 C18 conformers folded, 34 of them by contacts only 5–6 bonds apart. It also calls bench conformer 57 folded, although its measured genus is 0. Under it, the quoted "extended" spread of 2.69 / 5.16 becomes 2.51 / 4.44 [X].
- **Selection.** `_rg_spread` picks 3 of 10 C18 members that differ from the bench's at the same seed, because its Rg is a different quantity [X].

**Counts:** 1 HIGH, 8 MEDIUM, 17 LOW new findings.

---

## Baseline

| check | result |
| --- | --- |
| master fast suite | **1077 passed**, 62 skipped, 17 xfailed, 13.9 s, 0.29 GB [X] (`out/fast_baseline.txt`) |
| PR #73 fast suite | **1113 passed**, 62 skipped, 10 xfailed, 19.6 s [X] |
| master + small STO-3G slow subset | 1090 passed, 34 s [X] (`out/d_revert.txt`) |
| committed R11 set (105) | 97 ok, 8 safe but not preferred, **0 unsafe** [X] (`out/a_sets.txt`) |
| R11 held-out (58) | 57 ok, 1 not preferred (`draw_indole_1.25` refused) [X] |
| committed stereo set (106) | 339 tests passed [X] |
| stereo held-out (51, Publishing session's set), with the flag | 51/51 [X] |
| PR determinism, octanol K = 8, two processes | identical member sha1 `3cf96a45b5e9` [X] (`out/e_select_*.txt`) |
| bench C18 conformers re-embedded by the PR's embedder | all 10 Rg values match `bench/conformers.jsonl` to 1e-4 Å [X] |

I did not run the slow suite. Its 760 / 8 xfail at v0.2.1 is consistent with the 768 tests v0.2.1 collects [X, `out/c_counts.txt`].

---

## Disposition of review 5's R1–R18

"Guarded" means that, in a `git archive` copy of `dc9d9ef`, a minimal mutation reverting the fix fails a test that passes on the unmutated copy (`d_revert.py`, `out/d_revert.txt`). I checked that the module was imported from the copy, not from the stale editable install.

| # | finding | verdict | guarded? | evidence and residue |
| --- | --- | --- | --- | --- |
| R1 | file geometry origin; budget rows | **FIXED** | yes (3 mutations: anion row on charge alone, basis row always, declared origin ignored) | `geometry_method` is still `"str"` in `record-schema.json:39` [S]. A declared origin is trusted by design. |
| R2 | gate ignored `--charge`/`--spin` | **FIXED** | yes (no CHARGE_OVERRIDE, no SPIN_OVERRIDE, `main` passing `charge=None`) | Forced `--smiles O --charge 2` gives `validity.input` false and sites withheld [X]. `sites_reason` doesn't name CHARGE_OVERRIDE, though it is in `gate.findings`. PR #73 reopens this class for records 2…N (S2). |
| R3 | unbound anion not recorded | **FIXED** | yes | F⁻ gas: `validity.electronic` false, `UNBOUND_ANION_GAS` [X]. The HOMO > 0 criterion is a knife-edge (SCN⁻ −0.12 eV) [I]. See S9 for ALIE. |
| R4 | genus null vs int | **FIXED** | yes | |
| R5 | index labels on mapped input | **FIXED** | yes | |
| R6 | one calculation, several keys | **FIXED** (key) | yes; this includes `test_b3lyp_with_vwn5_is_a_different_key` | All 20 xc spellings collapse correctly. With `B3LYP_WITH_VWN5` set, the key becomes `b3lyp5` (ΔE 3.73e-2 Ha) [X] (`out/d_keys.txt`). Residue: S20. |
| R7 | uncalibrated symmetry tolerances | **FIXED** | yes | 1.5 × the max over 10 rotations [S] |
| R8 | forced mapping failure lost its reason | **PARTIAL** | yes (2 mutations) | `mapping_verified` is deferred to the next schema version. Residue: S21. |
| R9 | NPZ keyed on input and surface | **FIXED** (`tri_atom`) | yes | "`tri_atom` present iff `validity.sites`" holds on 8 records: valid, H₂, AlH₄⁻, forced BLOCK, forced charge, F⁻ gas, forced verify failure, ethanol file. No −1 anywhere [X] (`out/d_records.txt`). Residue: S22. |
| R10 | vocabularies typed from memory | **PARTIAL** | yes | `ACCEPT` is gone. `geometry_method` and `budget[*].source` are still `"str"` [S]. |
| R11 | C₆F₆/C₂F₄ discarded | **FIXED** | yes (94 failures when reverted) | New edge cases: S6, S7, S8, S15–S17. |
| R12 | substring claims | **FIXED** | yes (2.69→2.96 and 5.1→5.7 each fail) | |
| R13 | double-zeta recognition | **PARTIAL** | yes | Residue: S23. |
| R14 | ν null → 0 on one triangle | **FIXED** (ν) | yes | The same knife-edge has moved to σ²_tot: S19. |
| R15 | prior-art gaps | **FIXED** [S] | n/a | ALIE row is stale: S26. |
| R16 | doc drift | **PARTIAL** | n/a | `CLAUDE.md:13` says "~3000 lines"; the file has 3575 [X]. |
| R17 | P3 pins nothing | **FIXED** | yes (V_min +1e-4 fails) | |
| R18 | version from installed metadata | **FIXED** | yes | Under the stale editable 0.2.0 install, the record says 0.2.1, from the module that ran [X] (`out/d_version.txt`). Residue: S24. |

---

## Lead finding

### S1 · Ensembles mix diastereomers through stereo elements RDKit does not perceive, on both routes `[HIGH]` · CONFIRMED EXECUTION FAILURE

**Where**
- PR `esp_map.py`:
  - `plan_ensemble` gate: `:1758`, `_gate(..., not supplied)`;
  - generated member check: `:1776–1784`;
  - `supplied_members` stereo check: `:1692–1708`.
- Master: `coordinates_fix_stereo` `:661`, `stereo_left_open` `:804`.

**How it happens**
- Spec §1 says one stereoisomer per ensemble. The generated route relies on `--require-stereo` to guarantee it, and the supplied route on stereo assigned from 3D.
- Both rely on `FindPotentialStereo` / `AssignStereochemistryFrom3D`, which perceive tetrahedral centres and ordinary E/Z only.
- `coordinates_fix_stereo` declares every element fixed once the input is 3D.
- The generated route then **discards the 3D coordinates** (`RemoveAllConformers`) and re-embeds from the graph. The graph does not carry an axis or a planar configuration.

**Reproduction [X]** (`e_supplied.py`, `out/e_supplied.txt`). "Sign" is the sign of the torsion that defines the element.

| input | route | outcome |
| --- | --- | --- |
| `O[C@H]1CC/C=C/CCC1` (trans-cyclooctenol, centre fixed), K = 10 | generated from SMILES | ACCEPTED, 5 members, helicity `[+1, −1, +1, −1, −1]`: two diastereomers |
| `C[C@H](O)C=C=CC`, one (R, Ra) conformer as a 3D SDF, K = 10 | generated from a 3D file | ACCEPTED, 4 members, allene axis `[+1, −1, −1, +1]` |
| (R, Rp) + (R, Sp) trans-cyclooctenol | supplied SDF | ACCEPTED, identity `O[C@H]1CC/C=C/CCC1` |
| (R, Ra) + (R, Sa) allene | supplied | ACCEPTED |
| (R, Ra) + (R, Sa) tetra-ortho biaryl `Cc1cccc(C)c1-c1c(C)cccc1[C@H](C)O` | supplied | ACCEPTED |

Over 30 embeddings, ETKDG splits the trans-cyclooctenol helicity 16/14, the allene 18/12 and the biaryl 17/13.

The R11/stereo half (`out/a_stereo.txt`):
- `--require-stereo` finds nothing open on the trans-cyclooctenol, an o-tBu anilide C–N axis, an in-ring (7-ring) biaryl axis, or 1,12-dimethylbenzo[c]phenanthrene;
- the heavy-atom RMSD between the two trans-cyclooctenol helicities is 0.48 Å, under `DEDUP_RMSD` = 0.5, so dedup can also merge them.

**What holds.** Tetrahedral centres and E/Z survive the 3D-file route. Threonine from a 3D SDF or an XYZ gives all 9 members as (2S,3R) [X].

**Why HIGH.**
- This is E13's worst outcome ("an average over diastereomers"), reached by two routes the hostile set does not try.
- trans-Cyclooctenes are a mainstream compound class (bioorthogonal TCO reagents), whose axial and equatorial diastereomers are distinct compounds.
- Nothing in the plan says the members differ. The runner will average them.

**Correction**
1. After embedding, fingerprint each member's configuration over the elements RDKit misses:
   - the signed C–C=C–C torsion of every trans double bond in a ring of ≤ 8 atoms;
   - the axis torsion of every cumulene chain (`cumulene_chains`);
   - every biaryl or aryl–X axis (`atropisomer_axes`, widened as in S-e of `out/a_stereo.txt`).
2. For generated ensembles, require all members to share the first member's fingerprint. For a 3D-file input, require the file's fingerprint, or refuse generated ensembles from 3D input and point the user to `supplied`.
3. For supplied ensembles, require every record to share record 1's fingerprint (`MEMBER_STEREO`).
4. Extend `stereo_left_open` with trans-cycloalkenes (ring ≤ 8), so that `--require-stereo` refuses `O[C@H]1CC/C=C/CCC1` from SMILES. No SMILES can fix the helicity.

**Validation**
- Each row of the table is refused, or yields members with a single sign.
- Threonine from 3D still passes.
- Add these cases to `tests/data/ensemble` as strict xfails before the fix.

---

## Findings on PR #73

### S2 · Supplied records 2…N are never gated: three redox states of one skeleton form one ensemble `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** PR `supplied_members` `:1645–1716`. `check_structure` runs on record 1 only (`plan_ensemble` `:1758`).

**How it happens**
- `verify_computed` takes bond orders and charges from the approval (record 1), so a record with the same nuclei and a different electron count matches.
- `plan_ensemble` returns no per-member charge or spin, so the runner will compute every member as record 1's species.

**Reproduction [X]** (`e_supplied2.py`, `e_select.py tempo`; `out/e_supplied2.txt`, `out/e_select_a.txt`)

| supplied SDF | outcome |
| --- | --- |
| ethyl cation, radical, anion on one skeleton | **ACCEPTED**, 3 members, identity `[CH2+]C` |
| the same records with the radical first | REFUSED `OPEN_SHELL` (only record 1 is gated) |
| TEMPO⁺ (oxoammonium), TEMPO•, TEMPO⁻ | **ACCEPTED**, identity `…[N+]1=O` |
| maps 1,2,3 then the same molecule mapped 3,2,1 | ACCEPTED |
| maps, then duplicate maps (map 5 ×3) | ACCEPTED; the same record as record 1 is refused `DUPLICATE_ATOM_MAP` |
| maps, then no maps | ACCEPTED |

This is review 5 R2 one level up: invariant 2 ("every input is gated") applied to the first record only. The verdict depends on the record order.

**Correction**
- Run `check_structure` on every record.
- Require every record's formal charge, radical count and atom maps (heavy atoms) to equal record 1's. Refuse with `MEMBER_CHARGE`, `MEMBER_SPIN` or `MEMBER_MAPS`.
- Carry the approved charge and spin in the plan.

**Validation:** all four ACCEPTED rows are refused, and E7 still passes.

### S3 · The contact rule calls a third of C18 conformers folded by local contacts and disagrees with measured genus `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** PR `intramolecular_contacts` / `CONTACT_MIN_BONDS` `:1447–1489`.

**Reproduction [X]** (`e_select.py c18`, `e_contacts.py`, `e_extended.py`; outputs `out/e_select_b.txt`, `out/e_contacts.txt`, `out/e_extended.txt`). The PR's embedder reproduces the bench's 100 conformers (seed 0xC0FFEE, Rg matching to 1e-4).

**Rule vs measured genus on the 10 bench conformers**

| conf | genus | PR folded | contacts (bonds apart) |
| --- | --- | --- | --- |
| 40 | 2 | yes | 20 |
| 72 | 3 | yes | 13, 5 |
| **57** | **0** | **yes** | 5 (O3–C8, across the Z C=C), 5 (C11–C16, saturated) |
| the other 7 | 0 | no | |

**Over all 100 conformers**
- 36 are called folded, and **34 of them only by contacts 5–6 bonds apart**.
- Contact pairs by distance in bonds: 5 bonds, 24 saturated (the g+g− motif E14 says is not folding) and 29 across a cis C=C; 6 bonds, 3; 13 bonds, 1; 20 bonds, 1.
- The only long-range contacts are exactly bench conformers 40 and 72.

**Consequence for the quoted numbers.** The "extended" population changes from the 8 genus-0 conformers to 7, and the site spread changes with it:

| population | median site spread | max |
| --- | --- | --- |
| genus 0, as quoted | 2.69 | 5.16 |
| PR rule | **2.51** | **4.44** |

**Why it matters**
- E14 passes on its one fixture by 0.08 Å, while the population it stands for fails.
- `CONTACT_MIN_BONDS` → 7 and `CONTACT_SLACK` → 0.3 both pass every test (`out/e_mutate.txt`, M7 and M8), so the constants are unpinned.
- An idealised (non-MMFF) g+g− chain is called folded at 3.88 Å against a 3.90 Å limit (`out/e_select_a.txt`), so the verdict depends on force-field relaxation.

**Correction**
- Treat the rule as unmeasured (spec §8) and say so in the docstring.
- On this data, ≥ 7 bonds apart reproduces genus > 0 on all 10 measured conformers and gives a folded rate of 2/100. Measure before choosing.
- Pin it with a test over the 10 bench conformer ids: re-embed, then assert `folded == (genus > 0)`. The test takes seconds and needs no SCF.

**Validation:** that test fails today on conf 57 and passes after the change.

### S4 · `_rg_spread` is not the rule the measured conformer rows used `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** PR `:1529–1547` (its docstring cites `bench/conformers.py select()`) and `select_ensemble` `:1594`.

**How it differs from the bench**
- The bench ranks by `Descriptors3D.RadiusOfGyration`: mass-weighted, all atoms, over all 100 conformers.
- The PR ranks by heavy-atom, unweighted Rg, over converged, deduplicated conformers.

**Reproduction [X]** (`out/e_select_b.txt`), C18 at the same seed:
- PR: `[4, 40, 57, 28, 27, 54, 22, 67, 26, 72]`;
- bench, and the bench's `select()` re-run here: `[4, 28, 34, 40, 53, 54, 57, 67, 72, 82]`.

Three members differ. Conformer 82, the bench's most extended, is dropped. E4 therefore cannot reproduce the 5.1 / 2.69 / 5.2 rows it is meant to check.

**A second defect in the same function.** When the lowest-energy conformer sits at an Rg extreme, the top-up adds the next most compact conformer instead of re-spreading. Octanol, K = 3, gives Rg `[3.26, 2.51, 2.53]` (`out/e_select_a.txt`). The spec says the lowest-energy conformer is usually among the most extended, so this case is the common one.

**Correction**
- Either use the bench's Rg, or re-measure and drop the citation.
- Re-spread over the remaining range when a pick collides.

**Validation:** a dry E4 selection equals the bench ids, or the docstring and a `bench/` row say why not.

### S5 · The symmetry-match cap truncates silently and keeps symmetry copies `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** PR `select_ensemble` `:1570–1572` (`max_matches = 1000`).

**Reproduction [X]** (`out/e_select_b.txt`). `(tBuOCH2)3CH` has 1296 heavy-atom automorphisms.
- One conformer, relabelled by an automorphism past the cap, gives **RMSD 0.820 Å** with the cap and 0.000 Å without it, against `DEDUP_RMSD` 0.5.
- On 30 conformers, `n_unique` is 11 with the cap and 10 with 2000.
- The protocol records `n_symmetry_matches: 1000`, with no sign that the list was truncated.

This is E6's worst outcome ("symmetry copies inflating n"). Four CF₃ groups or three tBu groups are enough to reach it.

**Correction**
- If `len(matches) == max_matches`, refuse, or fall back to a canonical-rank alignment.
- Record `symmetry_matches_truncated`.

**Validation:** this molecule gives 10, or a refusal.

### S6 · `read_structure` drops unparseable records and renumbers the rest, on master and in supplied ensembles `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** master `esp_map.py:1215` (`mols = [m for m in mols if m is not None]`) and `:1233`. PR `supplied_members` `:283–298`.

**Reproduction [X]**

| input | outcome |
| --- | --- |
| `--record 2` on records `pose_1`, *unparseable*, `pose_3` | reads `pose_3` |
| `--record 3` on the same file | `RECORD_OUT_OF_RANGE` |
| a two-record file whose second record doesn't parse | read as one record, no `MULTIPLE_RECORDS` (invariant 9 bypassed) |
| supplied ensemble of 3 records, record 2 unparseable | **ACCEPTED, 2 members**; `MEMBER_UNREADABLE` never fires |

Sources: `out/a_r11b.txt` for the first three rows, `out/e_supplied.txt` §8 for the last. Spec §4 says excluded members are "never silently dropped".

**Correction:** keep `None` placeholders, number records as the file does, and name each unreadable record.

**Validation:** the supplied case is refused naming record 2, and `--record 2` names it.

---

## Findings on master (R11, `--require-stereo`, ALIE)

### S7 · A flat 3D-header drawing takes its chirality from z noise, and passes `--require-stereo` `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** `esp_map.py:1243`. `_coords = drawing` is set, but the chiral tags assigned from 3D are kept. `stereo_left_open` and the embedder treat them as specified.

**Reproduction [X]** (`out/a_stereo2.txt`, `out/a_stereo3.txt`). Unwedged alanine and threonine depictions with a 3D header and z noise, 40 seeds per amplitude; "pass" means the flag accepted the drawing:

| noise | alanine passes | threonine passes |
| --- | --- | --- |
| ±0.02 Å | 11/40 | 3/40 |
| ±0.05 Å | 26/40 | 19/40 |

- Two seeds of the same threonine drawing compute two diastereomers, both `verified: True`.
- Without the flag, `UNSPECIFIED_TETRAHEDRAL` disappears.

This is where the known missing chirality-volume floor matters more than stated: it affects every flat drawing with a 3D header, not only near-planar geometries.

**Correction:** for a `drawing`, re-perceive stereo from wedges only. Add a volume floor for 3D files.

**Validation:** every noisy threonine blocks under the flag and warns without it.

### S8 · A 2D header is ignored as soon as any z ≠ 0 `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** `esp_map.py:345–351`.

**Reproduction [X]** (`out/a_r11.txt`, `out/a_r11b.txt`), three files with a 2D header:

| file | read as | cost |
| --- | --- | --- |
| threonine drawing, wedge atoms at z ±0.8, heavy bonds 1.500 Å | `geometry` ("file 3D coordinates used as-is; 9 hydrogen(s) added") | 3 |
| chlorobenzene drawing, one atom at z = 0.3 | `geometry` | 3 |
| real 2-butanol | `geometry` (the manifest prefers `drawing`) | 2 |

**Correction:** a 2D header means `drawing`, or `ambiguous` if the file isn't flat.

**Validation:** all three cases give drawing or ambiguous.

### S9 · ALIE uses |ε|, so an unbound gas-phase anion yields a plausible positive Ī `[MEDIUM]` · CONFIRMED EXECUTION FAILURE

**Where:** `esp_map.py:1601`; `docs/p4-spec.md:186,234`.

**Reproduction [X]** (`c_alie_unbound.py`, `out/c_alie_unbound.txt`), def2-SVP, gas phase:

| anion | HOMO | I_min with \|ε\| | I_min with −ε |
| --- | --- | --- | --- |
| enolate | +1.38 eV | 2.19 eV | −0.04 eV, at another vertex |
| nitrite | +0.40 eV | 2.52 eV | 1.95 eV |

- In PCM water the two conventions agree (difference 0.000).
- Because of |ε|, positivity test 3 and the selftest's "ALIE positive" check cannot fail.
- The |ε| convention follows Sjoberg 1990, so this is a domain-of-validity defect, not a transcription error.

**Correction:** use −ε_i, or refuse any positive occupied ε. Withhold ALIE when `validity.electronic` is false.

**Validation:** the gas enolate gives a null ALIE with a problem code.

---

## Low

**S10 · A supplied record can swap diastereotopic groups' indices** [X] (`out/e_supplied2.txt` §9).
- Record 2 is record 1 with the two diastereotopic isopropyl methyls' indices swapped. Atom 4 goes from pro-S to pro-R (CIP with a ¹³C tag), and the ensemble is ACCEPTED.
- `verify_mapping` checks element and neighbours, which a consistent swap preserves.
- **Correction:** for graph-equivalent atoms on a prochiral centre, compare the pro-R/pro-S label by index across records.

**S11 · Isotope policy depends on record order and is inconsistent** [X] (`out/e_supplied2.txt`, `out/e_supplied.txt`).
- CD₃OH + CH₃OH is accepted in either order, and the identity is whichever comes first.
- (R)- and (S)-ethanol-1-d are refused as `MEMBER_STEREO`, although their ESP is identical and `--require-stereo` exempts isotope-only stereo.
- **Correction:** pick one rule (strip isotopes from identity, or require equal labels) and state it.

**S12 · H₂ crashes `plan_ensemble`** [X] (`out/e_supplied.txt` §13).
- `smiles="[H][H]"` gives `LinAlgError` in `shape_descriptors` (`:1500`): there are no heavy atoms.
- **Correction:** return null descriptors, or refuse with a code.

**S13 · `supplied_members` is quadratic** [X] (`out/e_scale.txt`).
- It re-parses the whole file once per record. Timings: 20 records 0.9 s, 100 records 6.7 s, 200 records 28.2 s. 1000 poses would take about 12 min [I].
- **Correction:** read the file once.

**S14 · Several PR guards are untested** [X] (`out/e_mutate.txt`).
- Each of these mutations passes both ensemble test files:
  - M1: remove the `computed != ref_smiles` backstop;
  - M3: skip the generated per-member `verify_computed`/`verify_mapping`;
  - M7, M8: change the contact constants.
- M2, M4, M5, M6, M9 and M10 are caught.
- `test_contacts_follow_the_stated_rule` reads `em.CONTACT_MIN_BONDS`, so it checks the code against itself.
- **Correction:** pin each guard with a case that reaches it.

**S15 · Text after the last `$$$$` misaligns `_molfile_header_dims`** [X] (`out/a_r11.txt`).
- A C₂F₄ SDF with one trailing comment line is read as a drawing, and its real geometry is discarded.
- **Correction:** align header dims to the records RDKit actually returns.

**S16 · Flatness is the maximum deviation from the least-squares plane** [X] (`out/a_stereo3.txt`).
- A drawing with |z| ≤ 0.080 Å has a plane deviation of 0.101 Å, so it is read as `geometry`. This happened in 3 of 80 drawings at ±0.09 Å.
- **Correction:** use the smaller of the plane and z = 0 deviations, or the RMS.

**S17 · Two real species are refused or discarded at the thresholds** [X] (`out/a_r11.txt`).
- Dicyanamide (N–C 1.31 Å, −11 % from reference) is refused as `ambiguous`.
- Heavy-atom-only disiloxane (Si–O–Si 144°, outside `SP3_REAL_ANGLE` 90–118°) is discarded as a drawing.
- **Correction:** widen the angle window for O between heavy atoms, or refuse rather than discard outside it.

**S18 · The ALIE selftest misses several real errors** [X] (`out/c_alie_mut.txt`).
- These mutations pass `--selftest`: dropping β, giving β the α energies, skipping the last chunk, and a chunk slice off by one.
- The weight mutations pass everything: a constant weight cancels in the ratio, so p4-spec §10.1's claim that the identity "catches a wrong occupation weight" is false.
- The slow test uses atol 1e-6, while the spec says 1e-8.
- **Correction:** add an OH-radical selftest with ≥ 3 probe points against a dense reference, and fix the §10.1 wording.

**S19 · σ²_tot drops a sign at the GIPF floor (R14 residue)** [X] (`out/d_gipf.txt`).
- At 4.9 % positive area σ²_tot is 8.64; at 5.0 % it is 17.65.
- **Correction:** σ²_tot is null whenever one sign is null, unless that sign's area is 0.

**S20 · `geometry_method` label ignores `B3LYP_WITH_VWN5` (R6 residue)** [X] (`out/d_keys.txt`).
- With VWN5 set, `--optimize` geometries are labelled `dft:b3lyp` while the key says `b3lyp5`.
- **Correction:** build the label from `xc_name`.

**S21 · The fallback `computed_smiles` doesn't parse (R8 residue)** [X] (`out/d_records.txt`, `out/d_smiles.txt`).
- A forced ethanol run with swapped coordinates records `[CH][OH2]C`, which RDKit rejects.
- `sites_reason` also names a mapping failure that is only a consequence of the connectivity failure.
- **Correction:** record `computed_smiles_kind`, and name the mapping cause only when `verify_computed` passed.

**S22 · NPZ schema (R9 residue)** [X].
- `schema["npz"]` lists conditional arrays as unconditional, so the H₂, F⁻, forced and AlH₄⁻ NPZ files all differ from it.
- NPZ `atom_maps` still uses 0 for "no map".
- **Correction:** mark the conditional arrays (e.g. `tri_atom_radius_scaled?`), and document or replace the 0 convention.

**S23 · Bases outside the list are not "small" (R13 residue)** [X].
- STO-6G, pcseg-0, 6-31++G, 6-31G(2d) and 6-31G(2df,p) give no warning and no budget row.
- **Correction:** classify by zeta, as R13 asked.

**S24 · Master still says 0.2.1 (R18 class)** [X] (`out/d_verrule.txt`).
- `git describe` gives `v0.2.1-24`, and the project's own `version_agrees_with_code` returns False.
- When `code` is "unknown", a master record cannot be told from v0.2.1.
- **Correction:** bump to `0.2.2.dev0` after each tag.

**S25 · Claims against data** [X] (`out/c_claims*.txt`, `out/c_r11_counts.txt`).
- `README.md:182–184` says "100 … seven"; the R11 set has 105 cases and 8 refusals. The v0.2.1 notes and the roadmap are right.
- The bromobenzene −0.76 basis shift (CLAUDE.md, README) has no `bench/` source.
- The README's "6.5 … nitrite" is thiocyanate's value; nitrite's maximum is 5.64.
- `README.md:252`: the span is 118, not 117.
- The density-fitting energy range in the roadmap is 6.2e-5 to 4.6e-4 Ha, not "6×10⁻⁵ to 3×10⁻⁴".
- README timing "141 basis functions 20.0 s" doesn't state that it is def2-SVP.
- `CLAUDE.md:13` line count (see R16).
- Every other bench-derived number checked reproduces, including the E14 margin of 0.584 Å.

**S26 · ALIE prior-art row** [S, Crossref X] (`out/c_crossref.txt`).
- `prior-art.md:34` still says "not implemented".
- It omits Politzer, Murray & Bulat 2010, *J. Mol. Model.* 16:1731, doi:10.1007/s00894-010-0709-5.
- It doesn't cite Multiwfn for surface ALIE, and doesn't state the |ε| convention.
- The existing bib entries match Crossref.

Each LOW item's validation is the inverse of its reproduction.

---

## What holds

- **Every review-5 fix is real and guarded.** Each of the 29 reverts fails a test.
- **Determinism.** Embedding and selection are deterministic across processes. The PR's embedder reproduces the bench conformers exactly, so the §8 measurements can be re-run on known members.
- **Ordinary stereo on the supplied route** [X]. These are all refused for the right reason:
  - a meso/chiral tartaric acid pair, E/Z pent-2-ene and one inverted centre (`MEMBER_STEREO`);
  - imidazole tautomers whose H counts match but sit on other atoms (`MEMBER_IDENTITY`);
  - a depiction with z noise under a 3D header (`MEMBER_IDENTITY`);
  - a reordered record (`MEMBER_ATOM_ORDER`);
  - a record missing hydrogens (`MEMBER_MISSING_H`).
- **The InChIKey is only a label, as it should be.** Identity is checked by SMILES, so E8 tautomers are refused.
- **Refusals happen before any embedding** (`test_refused_before_embedding`).
- **Mirror merging is correct.** It applies only to achiral molecules, meso included. `_is_achiral` ignores atom maps.
- **Never pads.** Benzene gives 1 member, the cyclohexane chair 1, butane 2.
- **ALIE is correct for bound systems.** It agrees with a dense orbital loop to ≤ 7e-15 Ha (H₂, water, OH•, NH₂•), is bit-identical across chunk sizes, and Σn|ψ|² = ρ.
- **The R11 and stereo sets are intact.** 0 unsafe on 105 + 58 cases, and 51/51 on the stereo held-out set.

---

## Summary table

| # | finding | severity | evidence |
| --- | --- | --- | --- |
| S1 | diastereomer mix via unperceived stereo (trans-cyclooctene, allene, biaryl), generated and supplied | **HIGH** | [X] |
| S2 | supplied records 2…N ungated: redox states, maps; order-dependent | MEDIUM | [X] |
| S3 | contact rule: 36/100 folded, conf 57 genus 0 called folded, extended 2.69→2.51 | MEDIUM | [X] |
| S4 | `_rg_spread` ≠ bench rule: 3/10 members differ; compact top-up | MEDIUM | [X] |
| S5 | symmetry cap truncates silently, duplicate kept | MEDIUM | [X] |
| S6 | unparseable records dropped and renumbered (master `--record`, supplied) | MEDIUM | [X] |
| S7 | z-noise drawings: stereo from noise passes `--require-stereo` | MEDIUM | [X] |
| S8 | 2D header ignored when z ≠ 0 | MEDIUM | [X] |
| S9 | ALIE \|ε\| on unbound anions; positivity test tautological | MEDIUM | [X] |
| S10–S26 | see Low | LOW | [X]/[S] |

## Verdict

**PR #73: changes required before merge.** Fix S1 (configuration fingerprint), S2 (gate every record) and S6 (no silent drops). Mark S3 and S4 as unmeasured, or fix them, since they decide which geometries become members. S5 and S12–S14 are small fixes.

**Master:** S7–S9 are worth fixing before the ensemble runner (step 4), because they feed it. S7 and S8 choose the configuration of a generated member, and S9 matters once ALIE enters records. The rest can wait.

**The structural lesson.** Reviews 4 and 5 were about the last mutation and the inputs a field can be reached from. This review's is about the **model a check shares with what it checks**. `--require-stereo`, `verify_computed` and the supplied stereo check all ask RDKit what the stereo is, so where RDKit is blind all three agree, and they agree wrongly. The same holds for the shape rules: each was validated on one fixture (E14) instead of the population it classifies.

---

## Before `espmap.ensemble/1` is frozen

Ranked by the cost of changing later.

1. **A configuration fingerprint per member (S1).** Record it in `members[*]`, and add `identity.configuration` beyond the SMILES. Without it, a record cannot show that its members are one stereoisomer.
2. **Approved charge, spin and atom maps at the top level, required equal across members (S2).** The spec's `identity` block holds only SMILES and InChIKey.
3. **`folded` must name its rule and parameters (S3).** Use `{min_bonds, slack}`, plus `genus_excess` agreement, and keep both populations. The current constants change the "extended" statistics.
4. **`protocol.selection` must say how Rg is defined, whether the symmetry list was truncated, and the dedup population (S4, S5).**
5. **`excluded` must include unreadable supplied records by their file record number (S6).**
6. **Isotope policy in `identity` (S11).**

---

## Heavy reproductions requested

1. **Do local contacts ever coincide with excess genus?**
   - Command: `bench/recheck_refused_conformers.py`-style single points, def2-SVP + RI, gas, 0.26 bohr, on 8 of the 34 C18 conformers that the PR calls folded by 5–6-bond contacts only (ids from `out/e_contacts.txt`), recording genus and site means.
   - Cost: about 3 GB and 4–5 min each, about 40 min in all.
   - Decides: the `CONTACT_MIN_BONDS` threshold (S3), and whether those conformers belong in the "extended" population.
2. **What a diastereomer mix does to the numbers.**
   - Command: trans-cyclooctenol (R, Rp) vs (R, Sp), one MMFF conformer each, def2-SVP, PCM water, site means.
   - Cost: about 1 GB and 3 min each.
   - Decides: how large the averaging error is in kcal/mol, and so how much S1 moves a site, and whether `MEMBER_STEREO` is protecting anything measurable.
3. **Bromobenzene −0.76 (S25).**
   - Command: `--molecule bromobenzene`, σ-hole at def2-SVP vs def2-TZVP, gas.
   - Cost: about 1.2 GB and 3 min.
   - Decides: whether the "molecule-specific" claim reproduces; it gives the claim a `bench/` row.
