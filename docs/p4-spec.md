# P4: reactivity descriptors on the surface — specification

Status: **revision 1. Deferred (2026-10-23).** The user made v0.3 about
ensembles and comparison records (`docs/ensemble-spec.md`). `alie_at_points`
itself is a standalone PR (#70); the record/2 integration below waits, and
returns framed as a measured surface property, not a reactivity predictor.
The calibration set (§6) stays committed. Originally: nothing here is
implemented. This is
stage P4 of `docs/per-site-spec.md` §8 ("ALIE, then frontier-orbital
amplitudes, then Fukui functions"), and the theme of v0.3. It fixes **what**
is computed and **what must be measured before any number is quoted**, before
the code is written, as R11 did for the drawing rule.

This revision specifies **ALIE** in full. Frontier amplitudes and Fukui
functions are outlined in §8 and get their own revision once ALIE has landed.

## Why

ESP shows where the potential is negative or positive. It does not show where
electrons are most easily given up. ALIE does, and on the ρ = 0.002 surface
its minima are an established guide to the sites of electrophilic aromatic
substitution: ortho/para over meta for NH₂, CH₃ and F, and meta for CHO and NO₂
(HF/6-31G*, `sjoberg1990`; `murray1992`). That published, checkable record is
the reason ALIE comes first.

**The ambident anions are a weaker test than first planned.** The 2026-10-22
pilot (§2) found ESP and ALIE pointing to different sites on all three:

| anion | ESP points to | pilot ALIE points to | what the literature says, read in full |
| --- | --- | --- | --- |
| thiocyanate | N | **S** | S: kS/kN = 10²–10³ in SN2 alkylation (`loos2003`) |
| enolate | O | C | conditions decide: O is kinetically favoured for a free enolate, C is thermodynamic (`mayr2011` §3.4.2) |
| nitrite | O | a 0.03 eV tie | conditions decide: SN2 O/N from 30:70 (EtOH) to 92:8 (THF) (`tishkov2005`); 55–65 % nitro in DMF (`kornblum1955`) |

So only thiocyanate grades a descriptor. For the enolate and nitrite, the site
moves with the counterion, solvent and alkylating agent, so no ground-state
descriptor can be "right" about them in general. Mayr's review argues that
these selectivities are set by intrinsic barriers and product stabilities
(Marcus theory), not by hard/soft matching (`mayr2011`). They are reported,
not graded (§6).

## 1. Definition

The **average local ionisation energy** at a point **r** (Sjoberg et al.
1990, `sjoberg1990`; on molecular surfaces, Bulat et al. 2010, `bulat2010`):

    Ī(r) = Σ_i n_i |ψ_i(r)|² |ε_i|  /  Σ_i n_i |ψ_i(r)|²

over **occupied** orbitals i, with occupation n_i and orbital energy ε_i.

- **Closed shell (RKS):** n_i = 2 for each occupied spatial orbital.
- **Open shell (UKS):** sum over α and β orbitals separately, n_i = 1 each.
  The denominator is then ρ_α + ρ_β, consistent with `total_dm` (invariant 5).
  The pilot's prototype handled RKS only.
- **Orbital energies are Kohn–Sham.** Ī therefore depends on the functional
  far more than ESP does: KS orbital energies are not ionisation energies.
  Absolute Ī values carry that dependence. Rankings between sites are the
  claim P4 makes, and §5 measures how far they survive a change of functional.
- **Units:** eV at the boundary, Hartree internally, as for V.

**Evaluated at the isosurface vertices,** never on a volumetric grid, like
`esp_at_points` (the implementation choice in `CLAUDE.md`). It needs the AO
values at the vertices (`numint.eval_ao`) and the MO coefficients: one pass,
chunked by memory like `esp_at_points`, and no new SCF.

**Established, not new.** Ī on a ρ-isosurface and its surface minima are
standard in the quantitative-surface literature, and Multiwfn computes them
(`lu2012mt`, `lu2012multiwfn`). `docs/prior-art.md` gains a row in the same
PR as the code.

## 2. Evidence already in hand

The pilot is uncommitted scratch work: a prototype ALIE (RKS only), three
anions, B3LYP/def2-TZVP, five conditions each. Its rankings held across:
- geometry (MMFF vs DFT);
- solvent (gas vs PCM water);
- basis (def2-TZVP vs def2-TZVPD);
- the two partition rules.

It also found:
- **Gas-phase anions are unbound in this basis** (positive HOMO), which makes
  their absolute ALIE a basis artefact. espmap already flags this
  (`UNBOUND_ANION_GAS`, review 5 R3), and ALIE inherits the consequence (§4).
- **Site minima and site means can disagree** (nitrite's HOMO maximum; the
  enolate's V_min). So every per-site descriptor reports the mean and the
  extreme together.
- **Not measured:** the functional axis. It is the largest expected axis for
  ALIE and the first row of §5.

The pilot's scripts and raw numbers move into `bench/` with the calibration
PR (§6), regenerated with the real implementation, not the prototype.

## 3. What a run reports (`espmap.record/2`)

New keys make this a new record schema. `espmap.record/2` also carries
`input.mapping_verified` (deferred from review 5 R8, spec §6). The key list
below is for review by Arrowfield before any code is written.

- **`settings.descriptors`:** the list computed, e.g. `["esp", "alie"]`.
  `--alie` turns ALIE on; ESP is always computed.
- **`results.alie`:**
  `{ "I_min_eV", "I_max_eV", "I_mean_eV", "I_p05_eV", "I_min_vertex",
  "I_min_atom" }`. Statistics are area-weighted over triangles, as
  `area_stats` does for V. p05 is the low tail, because low Ī is the reactive
  end. Vertex and atom name where the minimum is.
- **Per atom and per reaction site,** under **both** partition rules, as for V:
  `alie: { "mean_eV", "p05_eV", "min_eV" }`, and `delta` between the rules.
  "Lead with means" holds here too (the pilot's enolate).
- **NPZ:** an `ALIE` array per vertex, in eV, beside `V`.
- **`analysis_key`** is unchanged in its fields. Ī depends on exactly what V
  depends on: method, basis, solvent and convention, surface and geometry.
  Whether ALIE was computed is not a property of the numbers V has.
- **`render_key`** gains the property rendered (§7).

## 4. Validity

ALIE follows the record's `validity` exactly as V's per-site table does:
- **`validity.electronic` false** (an unbound gas-phase anion): `results.alie`
  is still written, labelled by the same problem code, but **no per-site ALIE
  table**. A positive HOMO makes the outer orbital's contribution a basis
  artefact. The remedy text is the one R3 already gives (`--solvent`).
- **`validity.input` or `validity.surface` false:** as for V's sites, no
  per-site table, and `sites_reason` says why.
- **Open shell:** allowed. The record states RKS or UKS, and ⟨S²⟩ is recorded
  so that spin contamination is visible.

## 5. Error budget: measured before any ALIE number is quoted

None of these numbers exists yet. Each row is a `bench/` measurement, with a
`tests/test_claims.py` check, in the error-budget PR (§9). The README and
`CLAUDE.md` quote **no ALIE value** until its row is in.

| axis | measured on | statistic |
| --- | --- | --- |
| **functional** (B3LYP / PBE0 / M06-2X / HF) | the calibration set | spread of site min and mean; **whether any site ranking flips** |
| basis (def2-SVP → def2-TZVPP; + diffuse for anions) | the calibration set | as above |
| solvent (gas → PCM water) | the anions, and one neutral | as above |
| geometry source (MMFF vs DFT) | the anions | as above |
| grid spacing (0.26 → 0.16 bohr) | chlorobenzene and one anion | site mean and min |
| partition rule | every site | `delta`, as for V |
| density fitting | one molecule over 500 nao | global and site |

The **functional row decides what P4 may claim.** If site rankings survive
all four functionals, rankings are quotable and absolute values carry the
spread. If a ranking flips, that case is reported as functional-dependent, in
the record's caveat and in Arrowfield's lesson.

## 6. Calibration set, before the code

Counterexamples first (the `counterexamples-first` skill, and R11's lesson).
Before ALIE is implemented, a calibration PR commits:
- **`tests/data/p4_calibration/manifest.json`**, with per molecule:
  - the reaction and its conditions;
  - the site that reacts;
  - the **literature source**, with a DOI verified on Crossref;
  - what ALIE is **expected** to do: `"points_to": site`, `"tie"`, or
    `"does_not_point_to": site`.
- the molecules' geometries, with their provenance.
- a test that runs once ALIE exists. Until then it is `xfail(strict=True)`,
  so ALIE's first PR must turn each expectation green or justify it in the
  open.

Contents, as committed in the calibration PR (#68), with every source checked
on Crossref and read in full except Stock & Brown 1963:

| molecule | known chemistry | expectation for ALIE |
| --- | --- | --- |
| thiocyanate | SN2: S, kS/kN 10²–10³ (`loos2003`) | S below N |
| toluene, fluorobenzene | EAS o/p; ALIE minima over o/p in `sjoberg1990` | o/p below meta |
| nitrobenzene | EAS meta; ALIE minima over meta in `sjoberg1990`, `murray1992` | meta below o/p |
| anisole, chlorobenzene | EAS o/p (`stock1963`) | o/p below meta |
| enolate | conditions decide (`mayr2011` §3.4.2) | report only |
| nitrite | conditions decide (`tishkov2005`, `kornblum1955`) | report only |
| phenoxide | O in most solvents; C with water, phenol, fluorinated alcohols, through hydrogen bonding (`kornblum1963`) | report only: PCM cannot represent the cause |

The published ALIE work is Hartree–Fock (Koopmans). espmap's ALIE uses
Kohn–Sham orbitals, so an aromatic case that fails under B3LYP is a finding
for the functional row of §5, not a reason to edit the case. `murray1992`
also shows a ranking changing with the basis (pyridine N-oxide: para at
6-31G*, ortho at 3-21G), which is why §5 measures the basis axis on rankings. **Arrowfield is asked to keep a held-out calibration set**,
built from the same kind of literature but not shared until ALIE is pushed.

## 7. Rendering

- `--property {esp,alie}` chooses the field drawn, default `esp`.
- **ALIE is sequential, never diverging.** It is always positive, so
  invariant 7 applies by construction.
- **Orientation, asserted in `--selftest`:** low Ī, the reactive end, is
  drawn **warm/dark**, matching ESP's "warm = where an electrophile goes".
  The caption states this. Invariant 6 is extended to the ALIE maps:
  orientation asserted, not reasoned about.
- `render_key` includes the property, so an ESP picture and an ALIE picture
  are never "comparable".

## 8. Later revisions: frontier amplitudes and Fukui functions (outline)

- **Frontier amplitudes:** |ψ_HOMO|² (and LUMO) at the vertices.
  - A degenerate HOMO (thiocyanate's π pair) is summed over the degenerate
    set, with the tolerance recorded.
  - HOMO−1 is reported when it lies within a stated window. The pilot's
    enolate HOMO−1 is O-centred and 1 eV lower.
- **Fukui f⁻ / f⁺:** vertical, ρ(N) − ρ(N∓1) at the same geometry, the
  radical in UKS with ⟨S²⟩ recorded. This is the first descriptor needing a
  second SCF, so its cost goes in `bench/timing.py`.
  - Condensed values are optional. If present, they name their scheme
    (Hirshfeld) and carry a second scheme beside them. The pilot's nitrite
    ties under Hirshfeld and favours N under Mulliken.
- Both reuse §3's shape: a site gains a field per descriptor.

## 9. Order of work

1. **This spec** (docs PR). Arrowfield reviews the record/2 key list in §3.
2. **Review 6 brief**, written now and run when step 4 opens, so that the
   whole repository is reviewed as v0.3 starts (`docs/roadmap-review.md`).
3. **Calibration set** (§6): data, expectations and citations; strict
   xfail; no code.
4. **ALIE implementation:** `alie_at_points` beside `esp_at_points`, RKS and
   UKS; `--alie`; record/2 keys; NPZ; `--property alie` rendering.
   - `--selftest` gains the one-orbital identity (below) and the ALIE
     colormap orientation.
   - Prior-art row.
5. **Error budget** (§5): `bench/alie_budget.py` and `.jsonl`, claims tests,
   and the README / `CLAUDE.md` rows. ALIE values are quotable only after this.
6. **Calibration result:** the xfail list shrinks or the expectation is
   challenged in the open, never edited to pass.

## 10. Acceptance tests

1. **One-orbital identity:** for any system with a single occupied orbital
   (H₂ or He in HF), Ī(r) = |ε₁| at every vertex, to 1e-10, in any basis.
   This is exact, and catches a wrong occupation weight or a missing |·|.
2. **RKS vs UKS:** a closed-shell molecule computed as UKS with spin 0 gives
   the same Ī as RKS, to 1e-8.
3. **Positivity:** Ī > 0 at every vertex of a bound system.
4. **Area-weighted statistics** reuse `area_stats` (the same weighting as V),
   with conservation tests as for V.
5. **Validity:** an unbound gas-phase anion writes `results.alie` and no
   per-site ALIE; the NPZ follows `validity.sites`, as for the partition
   arrays (R9).
6. **Schema:** `docs/record-schema.json` for record/2 is checked in both
   directions against a real `--alie` record, as for record/1.
7. **Orientation:** the ALIE colormap draws the minimum darkest/warm
   (`--selftest`).
8. **Calibration:** §6's strict expectations.

## Not in this spec

- Transition states, reaction barriers and product stability. The pilot's
  nitrite result suggests these are what decide its alkylation site. They
  are out of scope for ground-state descriptors.
- Enzyme or condensed-phase selectivity beyond PCM (`roadmap-review.md`,
  "Limits").
- Ensembles: one geometry per run, as for P0–P3.
