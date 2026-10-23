# Ensembles and comparison records: specification

Status: **revision 1, draft for review.** Nothing here is implemented. This
is the theme of v0.3 (the user's decision of 2026-10-23): conformer ensembles
as first-class records, and paired site comparisons between records, each
with its uncertainty stated. The hostile cases in §7 are written **before**
any code, as R11 and `--require-stereo` were, and become a committed set
before the first implementation PR.

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

## 2. Selection (generated ensembles)

1. **Embed and minimise** `n_embed` conformers (default 100; the protocol
   records it) with ETKDGv3 + MMFF94s, as today.
2. **Deduplicate** by symmetry-aware heavy-atom RMSD (`rdMolAlign.GetBestRMS`)
   below `dedup_rmsd` (default 0.5 Å, to be measured, §6), keeping the lower
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

**Energy windows and weights.** MMFF energies choose and order the members,
but **never weight** them. Their spread on the C18 acetate is 10.5 kcal/mol
over 100 conformers, so Boltzmann weights from MMFF would put nearly all the
weight on one member, from a force field whose relative conformer energies
are not validated here. Ensemble statistics are **unweighted** by default.
`--weight dft` (Boltzmann over the members' own SCF energies, at the stated
temperature) is allowed only when every member's energy comes from the same
`analysis_key`. The record then says so, and reports the effective number of
members (Kish: (Σw)² / Σw²), so a near-single-member ensemble is visible.

## 3. Shape descriptors, and extended vs folded

Per member, recorded:
- radius of gyration (heavy atoms) and maximum heavy-atom extent;
- asphericity, from the gyration tensor;
- surface genus and component count (`mesh_topology`, as today);
- **intramolecular contacts:** pairs of heavy atoms at least 5 bonds apart
  and closer than the sum of their Bondi radii plus 0.5 Å.

**Folded** is a property with a stated rule, not a label someone assigns:
- a member is **folded** if it has an intramolecular contact, or if its
  surface genus exceeds the genus of the molecule's own rings. A macrocycle
  that is merely large enough to open a hole isn't "folded" (§7, case E3);
- the record reports the folded count, and the statistics both over all
  members and over the extended members alone, because the measured
  conformer row is quoted both ways (5.1 vs 2.69).

## 4. The ensemble record (`espmap.ensemble/1`)

```
schema: espmap.ensemble/1
subject: opaque caller label from --subject, copied verbatim, never read
ensemble_key: hash of (member analysis_key, protocol)
analysis_key: the members' shared key (all members must share it)
protocol: {source: generated|supplied, n_embed, seed, dedup_rmsd,
           selection: {rule, k, ...}, weighting: none|dft, temperature_K?}
members: [{member_id, embed_conf_id, geometry_sha1, record_file, record_sha256,
           mmff_dE?, scf_energy_Ha, rg_A, extent_A, asphericity, genus,
           contacts, folded, validity}]
excluded: [{index, reason, folded, rg_A}]   # never silently dropped
convergence: {...}                          # §5
results: global V_min / V_max / surface stats per member, and over members
sites: per site and per rule: per-member means, then over members:
       mean, sd, min, max, range, n; and the same over extended members
validity: {members, convergence, bias}
```

- **Members are verifiable.** `member_id` is stable: the geometry sha1,
  which the same input and protocol reproduce. Each member's own
  `espmap.record` file is written beside the ensemble record, and its sha256
  is listed, so an archive can check every member without trusting the
  ensemble record.
- **Sites are keyed by their definition,** atom maps included, never by a
  label alone (review 5 R5). Each statistic states its population: all
  members, or extended only.
- **`subject`** is an opaque string the caller supplies (an archive's own
  compound ID, say). espmap copies it and never parses it, so no consumer's
  vocabulary enters the record. It is not part of `ensemble_key`: two
  callers computing the same thing get the same key.

**Exclusions are reported with their shape.** A member whose surface fails
integrity, or whose input or electronic validity is false, is excluded and
listed. If the exclusions fall disproportionately on folded members, the
ensemble's `validity.bias` is false. That is exactly the failure espmap
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

**Site correspondence** is explicit and chemical:
`{"by": "map"|"smarts", "pairs": [{"name", "A": <site>, "B": <site>}]}`.
SMARTS matching must be unique in each molecule, or the pair is refused.
Correspondence by atom index across two molecules is not offered, because
indices carry no chemistry (review 5 R5).

**Per site pair, under each partition rule:**
- `delta_mean`: B − A difference of site means (ensemble means for
  ensembles);
- `conformational`: each side's member SD and range (zero, and labelled
  single-geometry, for a single record);
- `partition`: the difference between the two rules' `delta_mean`;
- `sign_consistency`: for two ensembles, the fraction of member pairs (A_i,
  B_j) whose difference has the sign of `delta_mean`;
- `resolved`: true only if |delta_mean| exceeds the stated uncertainty
  **and** `sign_consistency` ≥ 0.9 (both thresholds to be measured, §8);
- `uncertainty_included`: the components that went into the decision. In
  v0.3 that's conformational and partition. Method, basis, solvent and
  geometry source are **named as not included** (they arrive as a formal
  object in v0.4) and never silently omitted.
- `rule`: the resolved/unresolved rule, its thresholds and their `bench/`
  source, so a reader can re-decide with a different rule.

The comparison record carries `subject_A` and `subject_B` copied from the
inputs, and the sha256 of each input record.

A comparison of two single-geometry records of flexible molecules is
allowed. Its conformational component is "not measured", and it can't be
`resolved` beyond the partition uncertainty alone. The record says so in a
stable code (`CONFORMERS_NOT_SAMPLED`).

## 7. Hostile cases, written first

These become `tests/data/ensemble/` with a manifest (allowed outcomes, the
preferred outcome, the cost of each) before any code, as for R11.
Expectations marked *measure* get their numbers from §8 first.

| # | input | what must happen | worst wrong outcome |
| --- | --- | --- | --- |
| E1 | benzene, `--ensemble 10` | one member after dedup, ensemble of 1, stated | ten copies averaged as if independent |
| E2 | cyclohexane | chair; a twist-boat only if within protocol and distinct | duplicate chairs counted twice |
| E3 | 18-crown-6 or a cyclic peptide | genus 1 from the ring is not "folded" | every member called folded |
| E4 | C18 acetate, K = 10 | folded members identified (contact or excess genus); stats both ways | folded members averaged in without a flag |
| E5 | E4 with the folded members' surfaces failing integrity | excluded, listed, `validity.bias` false | a converged-looking ensemble of only extended members |
| E6 | a symmetric molecule whose conformers differ by a symmetry operation | deduplicated by symmetry-aware RMSD | symmetry copies inflating n |
| E7 | supplied SDF: 20 docking poses of one ligand | accepted as an ensemble, protocol "supplied" | refused, or computed as pose 1 only |
| E8 | supplied SDF mixing a ligand and its tautomer | refused before compute | a mixed "ensemble" |
| E9 | supplied SDF, same molecule, atom order permuted in one record | refused, or reordered by verified mapping | site rows from the wrong atoms |
| E10 | members under two different bases | refused | a mixed-level ensemble |
| E11 | `--weight dft` with one member 5 kcal/mol below the rest | weights computed, Kish n ≈ 1 reported | a weighted mean presented as an ensemble of K |
| E12 | MMFF energies offered as weights | refused (not a weighting option) | a near-single-member average from force-field energies |
| C1 | compare A with itself | delta 0, unresolved | resolved |
| C2 | compare A and B at different bases | refused | a cross-level difference |
| C3 | a site pair whose SMARTS matches twice in B | refused | an arbitrary match |
| C4 | a difference smaller than the partition sensitivity | unresolved | resolved |
| C5 | a mean difference whose sign flips across 30 % of member pairs | unresolved (sign) | resolved on means alone |
| C6 | two single-geometry C18 records | allowed; `CONFORMERS_NOT_SAMPLED`; resolved only beyond partition | resolved with no conformer caveat |
| C7 | a non-converged ensemble vs a converged one | compared, with the non-convergence on the record | silent |
| C8 | a correspondence by atom index | refused (by map or SMARTS only) | index-based pairing |

## 8. Measured before thresholds are fixed

None of these defaults is quotable until measured into `bench/`:
- `dedup_rmsd`: the RMSD below which two MMFF minima give the same site
  means within the grid row (0.15 kcal/mol), on the C18 acetate and a
  smaller flexible molecule;
- `conv_tol` and K: site-mean convergence against K on the C18 acetate
  (K = 5, 10, 20), extended and folded separately;
- the `resolved` rule's thresholds: false-resolution rate on pairs that are
  the same molecule (C1-like) across seeds;
- contact distance (Bondi + 0.5 Å): its agreement with excess surface genus
  on the existing folded and extended C18 members.

Each run is heavy (K SCFs per molecule). One job at a time in WSL, logged
like `bench/timing.py`: resumable, one process per member.

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
