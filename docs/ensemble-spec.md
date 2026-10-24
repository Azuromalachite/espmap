# Ensembles and comparison records: specification

Status (v0.3.0 candidate, after review 7):
- **`espmap.ensemble/1` is implemented** (`--ensemble`, `run_ensemble`,
  `assemble_ensemble`) and its schema, `docs/ensemble-schema.json`, is
  frozen with v0.3.0.
- **`espmap.compare/1`** (`espmap compare`, `docs/compare-schema.json`) is a
  release candidate for v0.3.0.
- **Every generic threshold is provisional unless it names a measurement**:
  the convergence tolerance (`CONV_TOL_KCAL`, 1 kcal/mol) and the comparison
  rule (`espmap-0.3-provisional`, k = 2, s_min = 0.9) are not measured; the
  folded rule and the deduplication RMSD are (§8.2).

History: revision 2 of this document was a draft reviewed by Arrowfield as
the record's main consumer (stereo, identity, the resolved formula and nine
smaller points) before any code. This is the theme of v0.3 (the user's
decision of 2026-10-23): conformer ensembles as first-class records, and
paired site comparisons between records, each with its uncertainty stated.
The hostile cases in §7 were written **before** any code, as R11 and
`--require-stereo` were, and were committed before the first implementation
PR. Reviews 6 and 7 added their own counterexamples the same way.

## Boundary

espmap measures, it does not interpret. A comparison record may say that a
site contrast is *resolved* or *unresolved* under the uncertainty it states.
It never says which site reacts, which compound binds, or what a difference
means biologically. Those claims belong to the consumers (a pheromone atlas,
a teaching tool, a medicinal-chemistry notebook), and no feature here exists
because one consumer needs it.

## Why

A single geometry hides the largest measured axis for a flexible molecule.
On the C18 acetate (`bench/conformers.jsonl`, `bench/f9_recheck.jsonl`,
`bench/site_partition.jsonl`):
- the median site mean moves 5.1 kcal/mol across 10 conformers with the two
  folded ones included, and 2.69 over the 8 extended ones (CLAUDE.md, site
  budget);
- V_min moves 5.2 kcal/mol once folded conformers are included;
- the lowest-MMFF conformer, which espmap computes today, is one of the most
  extended: 87 of 100 embedded conformers are more compact.

So a comparison between two flexible molecules from one conformer each can
report a difference smaller than either molecule's own conformational
spread. That's the failure this spec designs out.

## 1. What an ensemble is

An **ensemble** is a set of conformers of **one** molecule, every one of
them computed under the **same** `analysis_key`, plus the protocol that chose
them. Two sources:

- **Generated** (`--ensemble K` with `--smiles` or a single-record file).
  `_embed_and_minimise` stays the only embedder (invariant 1). It already
  embeds `n_confs` conformers and keeps the lowest; it will return all of
  them, and a new `select_ensemble` chooses K.
- **Supplied** (`--ensemble supplied` with a multi-record SDF). Today a
  multi-record file is a BLOCK unless `--record N` (invariant 9). In ensemble
  mode it's accepted only when every record passes `verify_computed` against
  **one** approved structure and `verify_mapping` gives the same atom order,
  so a file mixing molecules, tautomers or protonation states is refused,
  not averaged.

"One geometry in, one record out" still holds per member: each member is a
complete `espmap.record` with its own `validity` and `geometry.sha1`. The
ensemble record refers to them by sha1; it doesn't replace them.

**One stereoisomer per ensemble.** Averaging diastereomers is the quietest
way an ensemble can be wrong:
- a **generated** ensemble requires fully specified input: the
  `--require-stereo` rule (`stereo_left_open`, master 3746cad) must find
  nothing open, or the run is refused before embedding. Otherwise ETKDG picks
  a configuration per embedding, and the "ensemble" mixes them;
- a **supplied** ensemble requires every member to be the **same**
  completion. Stereo is assigned from each member's 3D coordinates, and all
  members must agree. `verify_computed` alone accepts *some* completion per
  record, so it isn't enough here.

**Identity at the top.** The ensemble record carries the approved canonical
isomeric SMILES and its InChIKey. Every member's computed structure (stereo
from 3D included) must match the **SMILES**. The InChIKey is a label for
archives, not the check: standard InChI merges mobile-H tautomers, so
2-hydroxypyridine and 2-pyridone share one key, and an InChIKey check would
accept the tautomer mix that E8 must refuse (found by the hostile set's own
fixture test).

## 2. Selection (generated ensembles)

1. **Embed and minimise** `n_embed` conformers (default 100; the protocol
   records it) with ETKDGv3 + MMFF94s, as today.
2. **Deduplicate** by symmetry-aware heavy-atom RMSD (`rdMolAlign.GetBestRMS`)
   below `dedup_rmsd` (default 0.35 Å, measured in §8.2), keeping the lower
   MMFF energy. Methyl and hydroxyl rotors don't make a new conformer.
3. **Choose K for DFT** by a declared rule. The default is the one already
   validated in `bench/conformers.py`: the lowest-MMFF conformer, plus K−1
   spread evenly across the radius-of-gyration range, so the selection spans
   shape rather than repeating near-duplicate minima. The alternative is
   Butina clustering on RMSD, taking the lowest-energy member per cluster.
   The rule's name and parameters go into the record.
4. **Seed and order are fixed**, so the same input and protocol give the
   same members in the same order on any machine, up to RDKit and PySCF
   versions, which are recorded.

**K requested, K obtained.** Deduplication can leave fewer distinct
conformers than asked for (benzene leaves one). The record keeps both
`k_requested` and `k_obtained`, and never pads with duplicates.

**Energy windows and weights.** MMFF energies choose and order the members,
but **never weight** them. Their spread on the C18 acetate is 10.5 kcal/mol
over 100 conformers, so Boltzmann weights from MMFF would put nearly all the
weight on one member, from a force field whose relative conformer energies
are not validated here. Ensemble statistics are **unweighted** by default.

`--weight dft` (Boltzmann over the members' SCF energies, at the stated
temperature) requires **DFT-optimised** members (`--optimize` per member,
all under one `analysis_key`). Single-point energies on MMFF geometries
carry force-field strain, so they aren't conformer energies, and weighting
with them is refused. The record reports the effective number of members
(Kish: (Σw)² / Σw²), so a near-single-member ensemble is visible.

## 3. Shape descriptors, and extended vs folded

Per member, recorded:
- radius of gyration (heavy atoms) and maximum heavy-atom extent;
- asphericity, from the gyration tensor;
- surface genus and component count (`mesh_topology`, as today);
- **intramolecular contacts:** pairs of heavy atoms at least 6 bonds apart
  and closer than the sum of their Bondi radii plus 0.5 Å (measured, §8.2;
  the first draft said 5 bonds).

**Folded** is a property with a stated rule, not a label someone assigns:
- a member is **folded** if it has an intramolecular contact. The contact
  rule (≥ 6 bonds apart, Bondi sum + 0.5 Å) excludes the 1,5 pairs of an
  ordinary gauche kink, which are 4 bonds apart, and the 5-bond pairs of a
  g+g− kink or a cis-alkene neighbourhood, none of which closed a surface
  loop in the conformers measured (§8.2); E14 tests that a single g+g− kink
  isn't called folded;
- surface genus is a **diagnostic**, not part of the rule, because a ring's
  own genus can't be predicted before computing it: a macrocycle may or may
  not open a hole at ρ = 0.002. Each member records `genus_excess`, its genus
  minus that of the most extended member (largest Rg). §8 checks the contact
  rule against it on known members, and the record lists any member where
  the two disagree;
- the record reports the folded count, and the statistics both over all
  members and over the extended members alone, because the measured
  conformer row is quoted both ways (5.1 vs 2.69).

## 4. The ensemble record (`espmap.ensemble/1`)

The schema is `docs/ensemble-schema.json`, reviewed by a consumer (Arrowfield,
against the Atlas's `spa.espmap_ensemble/0`) before the runner was written, and
**frozen with v0.3.0**: a change is `espmap.ensemble/2`, never an edit in place.
Its status fields say what is measured and what is not: `folded_rule` and the
selection's dedup threshold are measured (§8.2); the convergence tolerance is
provisional.
Members are labelled `conf_00`, `conf_01`, ... in selection order, zero-based,
so files already written under the Atlas's index stay valid. The Atlas's
`protocol` string is its own directory label and is not this record's
`protocol` object.

```
schema: espmap.ensemble/1
subject: opaque caller label from --subject, copied verbatim, never read
ensemble_key: hash of (schema, identity, analysis_key, protocol, ordered
              member ids, excluded conformers) (review 7 R7-2)
analysis_key: the members' shared key (all members must share it)
protocol: {source: generated|supplied, n_embed, seed, dedup_rmsd,
           selection: {rule, k, ...}, weighting: none|dft, temperature_K?}
identity: {smiles (canonical, isomeric), inchikey}
k_requested, k_obtained
members: [{member_id, embed_conf_id, geometry_sha1, record_file, record_sha256,
           npz_file?, npz_sha256?, mmff_dE?, scf_energy_Ha, optimized,
           rg_A, extent_A, asphericity, genus, genus_excess, contacts,
           folded, validity}]
excluded: [{index, reason, folded, rg_A}]   # never silently dropped
convergence: {...}                          # §5
results: global V_min / V_max / surface stats per member, and over members
sites: per site and per rule: per-member means, then over members:
       mean, sd, min, max, range, n; and the same over extended members
validity: {members, convergence (true | false | null: not judgeable), unbiased}
```

- **Members are verifiable.** `member_id` is the sha1 of the member's
  geometry in the **verified-mapping atom order**, so a permuted copy of the
  same conformer gets the same id. It is stable **within the recorded RDKit
  and PySCF versions**: ETKDG output changes between RDKit releases. Across
  versions, `embed_conf_id` with the protocol seed is the identity. Each
  member's own `espmap.record` file (and NPZ, when written) sits beside the
  ensemble record with its sha256 listed, so an archive can check every
  member without trusting the ensemble record.
- **Statistics conventions:** `sd` uses ddof = 1; `range` is max − min over
  the members in the stated population.
- **Sites are keyed by their definition,** atom maps included, never by a
  label alone (review 5 R5). Each statistic states its population: all
  members, or extended only.
- **`subject`** is an opaque string the caller supplies (an archive's own
  compound ID, say). espmap copies it and never parses it, so no consumer's
  vocabulary enters the record. It is not part of `ensemble_key`: two
  callers computing the same thing get the same key.

**Exclusions are reported with their shape.** A member whose surface fails
integrity, or whose input or electronic validity is false, is excluded and
listed. `validity.unbiased` is false **whenever any excluded member is folded**
(named so that true reads as fine, like every other validity flag; `bias` in
the first draft).
That rule is conservative and testable: "disproportionately" wouldn't be. That is exactly the failure espmap
shipped once: the χ = 2 check refused precisely the two folded C18
conformers, and the ensemble looked converged at 1.1 kcal/mol until they
were recovered (F9, `bench/f9_recheck.jsonl`).

## 5. Convergence diagnostics

For each site statistic, over the members in selection order:
- the running mean and range as members are added;
- a bootstrap interval of the site mean over members;
- `converged`: true when adding the last third of the members moved no site
  mean by more than `conv_tol` (to be measured, §6), and false otherwise.
  The record says which site failed.

A non-converged ensemble still writes its record. Comparisons that use it
say so (§6).

## 6. Comparison records (`espmap compare`)

```
espmap compare A.json B.json --sites correspondence.json --out cmp.json
```

`A` and `B` are single-geometry records or ensemble records.

**Refusals before anything is compared:**
- the `analysis_key` fields that set the numbers (method, basis, solvent and
  convention, surface, grid, density fitting, geometry protocol) differ
  between A and B. A comparison across levels of theory is refused, never
  computed with a footnote;
- a site in the correspondence file doesn't resolve in A or in B;
- either record's `validity` is false where the compared quantity needs it.

**Mixed sources and versions:**
- an ensemble compared with a single record of a flexible molecule is
  **allowed**. The single side contributes sd = 0 and carries
  `CONFORMERS_NOT_SAMPLED`. Generated and supplied sources may be compared
  only if their `geometry_method` (part of the `analysis_key`) matches;
- two records with equal `analysis_key` but different espmap versions are
  **allowed**. Both versions are recorded and the comparison carries
  `VERSIONS_DIFFER`. The key is meant to say when numbers are comparable;
  the version flag lets a reader check that claim.

**Site correspondence** is explicit and chemical:
`{"by": "map"|"smarts", "pairs": [{"name", "A": <site>, "B": <site>}]}`.
- **By map:** a pair asserts homology, so mapped atoms must have the same
  element in A and B, or the pair is refused.
- **By SMARTS:** the match must be unique in each molecule. Matches that are
  symmetry-equivalent (equal canonical rank, `breakTies=False`), such as
  acetate's two oxygens, count as one site, the union of the equivalent
  atoms. Any other multiple match is refused.
- **By atom index** across two molecules: not offered, because indices carry
  no chemistry (review 5 R5).

**Per site pair, under each partition rule:**
- `delta_mean`: B − A difference of site means (ensemble means for
  ensembles);
- `conformational`: each side's member SD and range (zero, and labelled
  single-geometry, for a single record);
- `partition`: the difference between the two rules' `delta_mean`;
- `sign_consistency`: for two ensembles, the fraction of member pairs (A_i,
  B_j) whose difference has the sign of `delta_mean`;
- `resolved`, by a stated formula, so any reader re-deciding gets the same
  answer:

      U = k · sqrt(sd_A² + sd_B²) + |partition|
      resolved  ⇔  |delta_mean| > U  and  sign_consistency ≥ s_min

  - `sd_A`, `sd_B` are each side's member standard deviation (ddof = 1),
    not a standard error. The members are a deliberately spread sample of
    shapes, not independent draws, so the question is whether the two
    distributions separate, not whether two means are estimated precisely;
  - the partition sensitivity is systematic, not random, so it is added
    linearly rather than in quadrature;
  - `k` (default 2) and `s_min` (default 0.9) are provisional until §8
    measures the false-resolution rate;
  - for a single-geometry side, sd = 0 and the comparison carries
    `CONFORMERS_NOT_SAMPLED`;
- `uncertainty_included`: the components that went into the decision. In
  v0.3 that's conformational and partition. Method, basis, solvent and
  geometry source are **named as not included** (they arrive as a formal
  object in v0.4) and never silently omitted.
- `rule`: the resolved/unresolved rule, its thresholds and their `bench/`
  source, so a reader can re-decide with a different rule.

The comparison record carries `subject_A` and `subject_B` copied from the
inputs, and the sha256 of each input record.

**Additions from the consumer review** (Arrowfield, against the Atlas's
prospective analysis plan, before the code):
- a site may include its own hydrogens, `"n+H"` by map or `include_h` by
  SMARTS, recombined exactly from the members' per-atom rows (triangles are
  owned per atom); single-atom values stay diagnostic;
- per-member rows per side (`conformer`, `folded`, both rules), so any other
  decision rule, a bootstrap or an extended-only reading can be recomputed
  from the record alone;
- the partition used by the rule is named (`delta_between_rules`); two other
  readings are recorded as descriptive numbers (`between_side_means`,
  `max_member`), never used by the rule;
- a `global_extrema` block for V_min and V_max, shaped like a site;
- per ensemble side, the excluded count, their folded flags,
  `validity.unbiased` (flag `SELECTION_BIASED` when false) and convergence;
- `resolved_under` names the rule, so `resolved` is never read as a
  consumer's own verdict; a single-geometry side says
  `measured: not_measured`.

A comparison with a single-geometry side (a single record, or an ensemble
that kept one member) is **allowed but never resolved** (the user's
decision, 2026-10-24, an invariant of the code): delta, partition and the
descriptive numbers are reported, that side's conformational uncertainty is
`not_measured`, `resolved` is false with a `not_resolvable` reason, and the
comparison carries `CONFORMERS_NOT_SAMPLED`. Refusing such a comparison
belongs to a teaching layer, not to the measurement.

Two more gates on `resolved` (review 7), each an invariant of the code, each
stated in `not_resolvable` when it applies:
- **A site must have a value on every member under both rules** (R7-4). A
  folded member can bury a site. Dropping it would shrink n silently, so
  such a site is `evaluable: false`, `not_evaluable` names the side,
  conformer and rule, and its numbers are reported unresolved.
- **An ensemble side's convergence must be true** (R7-5). `false`, or `null`
  (not judgeable, as with K = 2), reports the numbers unresolved: "ensemble
  convergence not established".

Before anything is compared, both files are validated (R7-8). An unknown
schema, a foreign file or a missing field is `COMPARE_SCHEMA_UNSUPPORTED`.
A member whose sha256, analysis_key or geometry differs from what its parent
lists is `COMPARE_MEMBER_CHANGED` or `COMPARE_MEMBER_MISMATCH`.
`docs/compare-schema.json` is the record's JSON Schema.

## 7. Hostile cases, written first

These 29 cases become `tests/data/ensemble/` with a manifest (allowed outcomes, the
preferred outcome, the cost of each) before any code, as for R11.
Expectations marked *measure* get their numbers from §8 first.

| # | input | what must happen | worst wrong outcome |
| --- | --- | --- | --- |
| E1 | benzene, `--ensemble 10` | one member after dedup, ensemble of 1, stated | ten copies averaged as if independent |
| E2 | cyclohexane | chair; a twist-boat only if within protocol and distinct | duplicate chairs counted twice |
| E3 | 18-crown-6 or a cyclic peptide | genus 1 from the ring is not "folded" | every member called folded |
| E4 | C18 acetate, K = 10 | folded members identified (contact or excess genus); stats both ways | folded members averaged in without a flag |
| E5 | E4 with the folded members' surfaces failing integrity | excluded, listed, `validity.unbiased` false | a converged-looking ensemble of only extended members |
| E6 | a symmetric molecule whose conformers differ by a symmetry operation | deduplicated by symmetry-aware RMSD | symmetry copies inflating n |
| E7 | supplied SDF: 20 docking poses of one ligand | accepted as an ensemble, protocol "supplied" | refused, or computed as pose 1 only |
| E8 | supplied SDF mixing a ligand and its tautomer | refused before compute | a mixed "ensemble" |
| E9 | supplied SDF, same molecule, atom order permuted in one record | refused, or reordered by verified mapping | site rows from the wrong atoms |
| E10 | members under two different bases | refused | a mixed-level ensemble |
| E11 | `--weight dft` with one member 5 kcal/mol below the rest | weights computed, Kish n ≈ 1 reported | a weighted mean presented as an ensemble of K |
| E12 | MMFF energies offered as weights | refused (not a weighting option) | a near-single-member average from force-field energies |
| E13 | threonine SMILES with no stereo, `--ensemble` | refused before embedding (stereo left open) | an average over diastereomers |
| E14 | a C18 member with a single g+g− kink | not folded | an ordinary kink counted as folding |
| E15 | supplied SDF, one connectivity, one centre inverted in one member | refused (members are not one stereoisomer) | a mixed-stereoisomer ensemble |
| E16 | `--weight dft` on single-point energies of MMFF geometries | refused (members not optimised) | weights carrying force-field strain |
| E17 | a supplied SDF whose identity doesn't match the ensemble's InChIKey | refused | a record whose identity and members disagree |
| C1 | compare A with itself | delta 0, unresolved | resolved |
| C2 | compare A and B at different bases | refused | a cross-level difference |
| C3 | a site pair whose SMARTS matches twice in B | refused | an arbitrary match |
| C4 | a difference smaller than the partition sensitivity | unresolved | resolved |
| C5 | a mean difference whose sign flips across 30 % of member pairs | unresolved (sign) | resolved on means alone |
| C6 | two single-geometry C18 records | allowed; `CONFORMERS_NOT_SAMPLED`; resolved only beyond partition | resolved with no conformer caveat |
| C7 | a non-converged ensemble vs a converged one | compared, with the non-convergence on the record | silent |
| C8 | a correspondence by atom index | refused (by map or SMARTS only) | index-based pairing |
| C9 | a map pair whose mapped atoms differ in element | refused | a non-homologous "pair" |
| C10 | a SMARTS site matching acetate's two equivalent oxygens | one site, their union | refused, or one arbitrary oxygen |
| C11 | equal `analysis_key`, different espmap versions | compared, `VERSIONS_DIFFER`, both versions recorded | silent |
| C12 | an ensemble vs a single record of the same flexible molecule | compared; single side `CONFORMERS_NOT_SAMPLED` | refused, or the single side treated as sampled |

## 8. Measured before thresholds are fixed

None of these defaults is quotable until measured into `bench/`:
- `dedup_rmsd`: the RMSD below which two MMFF minima give the same site
  means within the grid row (0.15 kcal/mol), on the C18 acetate and a
  smaller flexible molecule;
- `conv_tol` and K: site-mean convergence against K on the C18 acetate
  (K = 5, 10, 20), extended and folded separately;
- the `resolved` rule's thresholds: false-resolution rate on pairs that are
  the same molecule (C1-like) across seeds;
- contact distance (Bondi + 0.5 Å): its agreement with excess surface genus.
  The E14 fixture (an otherwise all-anti octadecyl acetate with one g+g−
  kink) has its closest pair ≥ 5 bonds apart 0.58 Å beyond the Bondi sum, so
  it clears the rule by only 0.08 Å. Real conformers with a kink near a cis
  alkene may not. Measure the false-folded rate on extended members before
  fixing the margin;
  on the existing folded and extended C18 members.

Each run is heavy (K SCFs per molecule). One job at a time in WSL, logged
like `bench/timing.py`: resumable, one process per member.

### 8.1 Measured so far without an SCF (`bench/ensemble_rdkit.jsonl`)

RDKit only, seed 0xC0FFEE, 100 embedded, espmap's own embedder and
selection. These settle the structural half of three questions; none of
the constants is changed here (thresholds change in their own PRs).

**Dedup threshold.** Truth for chains: the backbone torsion string (anti,
g+, g−), canonical under reversal and reflection. A false merge is a
torsion class lost.

| threshold (Å) | n-butane | n-pentane | n-hexane | n-heptane | 1-hexanol |
| --- | --- | --- | --- | --- | --- |
| 0.25 | 0 | 0 | 0 | 0 | 0 |
| 0.35 | 0 | 0 | 0 | 0 | 0 |
| **0.5 (current)** | 0 | **1** | **1** | **2** | **2** |
| 0.75 | 1 | 2 | 7 | 12 | 13 |

(false merges; classes present: 2, 4, 9, 17, 19.) 0.5 Å loses real
conformers from pentane on; 0.35 Å loses none here. Below 0.35 Å a few
same-class pairs stay separate (hexane, heptane: 1; hexanol: up to 10,
largely hydroxyl-rotor variants a 3-state backbone class does not
distinguish). Cyclohexane has no backbone torsions, so the class truth
does not apply to it. **Not yet settled:** whether a merged pair's site
means differ by more than the grid row (0.15 kcal/mol), the spec's actual
criterion: `bench/ensemble_dft.py dedup_esp` runs the 5 pairs.

**Contact rule vs measured genus** (the 10 bench C18 conformers; genus from
`conformers.jsonl` χ and `f9_recheck.jsonl`), and the folded rate over all
100:

| min bonds | slack 0.3 Å | slack 0.5 Å |
| --- | --- | --- |
| **5 (current)** | 10/10 agree, 17 % folded | 9/10 (conformer 57 genus 0, called folded), 36 % |
| 6 | 10/10, 2 % | 10/10, 5 % |
| 7 | 10/10, 2 % | 10/10, 2 % |

At 5 bonds, 53 of the 58 contacts in 100 conformers are 5 bonds apart (the
g+g− and cis-alkene neighbourhoods review 6 S3 found). The rules ≥ 6 and
≥ 7 bonds differ only on conformers 3, 7 and 70 (6-bond contacts):
`bench/ensemble_dft.py genus` computes their genus, with four 5-bond-only
ones as controls, and decides between them.

**Re-spread after a collision** (review 6 S4), largest Rg gap left among the
selected members, Å:

| case | bench rule | re-spread | Rg range |
| --- | --- | --- | --- |
| C18 acetate, K = 10 | 1.16 | 1.19 | 2.53 |
| 1-octanol, K = 3 | 0.74 | 0.64 | 0.76 |
| 1-octanol, K = 5 | 0.36 | 0.25 | 0.76 |
| 1-hexanol, K = 6 | 0.20 | 0.17 | 0.61 |

Re-spreading helps at small K and does nothing for the C18 selection the
measured conformer rows used, so it stays deferred: a change would move
those rows' members for no gain on them.

### 8.2 Measured with an SCF (`bench/ensemble_dft.jsonl`)

**Contact rule vs surface genus.** Seven more C18 conformers were computed
(B3LYP/def2-SVP + RI, gas, 0.26 bohr, as the bench conformer rows):

| conformer | contacts | genus | ≥ 5 bonds | ≥ 6 bonds | ≥ 7 bonds |
| --- | --- | --- | --- | --- | --- |
| 70 | 6 bonds apart | **2** | folded | folded | **missed** |
| 3, 7 | 6 bonds apart | 0 | folded | folded | extended |
| 2, 5, 8, 11 | 5 bonds apart only | 0 | **folded** | extended | extended |

Over the 17 conformers with a measured genus, ≥ 5 bonds (current) is wrong on
7. ≥ 7 bonds is wrong on 1, but that one is a folded member counted as
extended, the error that contaminates the extended statistics. ≥ 6 bonds is
wrong on 2, both conservative (a genus-0 member left out of the extended
population), and misses no folded member. **The measurements favour ≥ 6
bonds with the 0.5 Å slack**; the change is its own PR.

**Adopted** (the threshold PR after #81): `CONTACT_MIN_BONDS` 6 and
`DEDUP_RMSD` 0.35 Å. Records state both as measured, citing this section.

**Dedup threshold vs site means.** Every pair that 0.5 Å merges and 0.35 Å
keeps differs in its largest site mean by more than the grid row
(0.15 kcal/mol), under both rules:

| pair | torsion classes | largest site-mean difference (kcal/mol) |
| --- | --- | --- |
| n-pentane 22 / 25 | (+,−) / (+,+) | 0.71 |
| n-hexane 5 / 76 | (+,−,a) / (+,a,a) | 0.48 |
| n-heptane 41 / 90 | same class | 0.59 |
| n-heptane 59 / 20 | (+,−,a,+) / (+,+,a,−) | 0.86 |
| n-heptane 85 / 62 | (+,−,a,a) / (+,a,a,−) | 0.78 |

So 0.5 Å merges conformers the ESP can tell apart. **The measurements favour
0.35 Å or lower**; the change is its own PR, and E2 and E6's expectations
(at most 2 members) still hold at 0.35 Å (`ensemble_rdkit.jsonl`).

**Review 6 S1's size.** (R)-trans-cyclooct-4-enol, one MMFF conformer per
ring helicity, B3LYP/def2-SVP, PCM water: V_min −45.63 vs −40.74, V_max
55.21 vs 60.83 kcal/mol, and per-site means differing by up to 12.2 kcal/mol
(C5). An ensemble mixing the two helicities would average numbers 5–12
kcal/mol apart; the configuration check (`configuration_signs`) is not a
formality.

**Review 6 S25.** Bromobenzene σ-hole V_max, gas, idealised geometry:
def2-SVP 16.81, def2-TZVP 16.05, def2-TZVPP 16.08 kcal/mol. SVP → TZVPP is
**−0.74**, against the −0.76 quoted in `CLAUDE.md` and the README: the claim
reproduces within 0.02 and now has a bench row.

## 9. Order of work

1. This spec (docs PR).
2. The hostile set (§7): data and expectations, strict xfail; no code.
3. `select_ensemble` and the supplied-ensemble identity checks, with
   `_embed_and_minimise` returning all conformers; no SCF changes.
4. The ensemble runner and `espmap.ensemble/1` (members as ordinary
   records; resumable).
5. `espmap compare` and `espmap.compare/1`.
6. The §8 measurements, then the thresholds, in their own PRs, with the
   README and `CLAUDE.md` rows.

Review 6 runs at step 3, the first code PR, as the roadmap schedules it for
the start of new record work.

## Not in this spec

- Uncertainty components beyond conformational and partition (v0.4).
- Site fingerprints and patch topology (v0.5).
- Batch execution across many molecules, HPC and caching (v0.6). Step 4's
  runner is resumable per member, not a batch system.
- Free energies, populations from MD, or solvent-explicit sampling.
- ALIE in ensemble records. ALIE's record integration is deferred; when it
  returns, it reuses §4's per-site shape.
