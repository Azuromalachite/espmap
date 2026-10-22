# Per-site output: specification

Status: **draft, revision 3. P0 is done; P2 (#24), P1 (#27) and caller-defined sites (#28) are in review.** This is the
contract that espmap's `--json` will grow into and that Arrowfield's v0.0
fixtures are written against. It extends item 1 of `docs/roadmap-review.md`
("extend `--json` with shape descriptors") from whole-molecule numbers to
per-atom ones.

Labels as in the roadmap review: `[V]` verified in the repo · `[M]` measured ·
`[I]` inferred · `[U]` not verifiable yet, a test must settle it.

**Revision 3 (2026-10-22).** Brought up to date with v0.1.1 and reviewed for
Arrowfield's needs:
- Done since revision 2: supplied-hydrogen freezing (#6), H1 (#7), H2 (#10).
  §0, §1, §4 and §8 now describe the code as it is.
- **Two comparability keys**, not one: `analysis_key` for the numbers and
  `render_key` for the picture (§5). A colormap must never decide whether two
  calculations can be compared.
- **A `validity` envelope** at the top of the record (§6), so a consumer reads
  one field instead of inferring validity from stderr or χ.
- **One geometry in, one record out** (Boundary): a file with more than one
  record is refused unless one is chosen.
- **H3 added** as a prerequisite for anything a learner sees (§0).
- **Surface definition stated** in the record (§6): ρ = 0.002, not SASA.
- **Robust site statistics** (area-weighted percentiles) beside the extrema
  (§3).
- **Mapped-hydrogen refusal is a v1 limitation**, and every hydrogen carries
  enough provenance to be promoted to a persistent identity later (§1, §4).
- The two RDKit unknowns about atom maps are now measured (§1).

Revision 2 rebuilt revision 1 on `docs/adversarial-review-3.md`: identity from
the coordinates (H2), no site block on a damaged surface (H4), budget rows with
sources (H8), geometry in the record (H14), atom and reaction-site levels,
both partition rules reported.

## Boundary

espmap decides what a defensible site descriptor is. Arrowfield decides how to
teach with it. espmap never reads AF3 or GNINA output, never knows a protein is
involved, and never sees a docking score.

**One geometry in, one record out.** A run computes one molecule in one
geometry and writes one record. A file holding several records (a GNINA
multi-pose SDF) is a gate BLOCK (`MULTIPLE_RECORDS`) unless `--record N`
names one, and the record then says which (done, #17). Before that, such a
file warned and used the first record: fine for a person at a terminal, but in
a pipeline 20 poses would silently have become an analysis of pose 1.
Arrowfield's adapter splits multi-pose files into one file per pose.

**Input for a pose is SDF (or MOL) with explicit bond orders.** PDB and XYZ
are accepted, but their bond orders are perceived from geometry, not given
(H1, fixed in #7), so the mapped structure Arrowfield intended is not
guaranteed. Arrowfield's adapters assign bond orders from the mapped ligand
definition before writing the SDF.

## Prior art

Nothing in this spec is new theory, and it claims none; `docs/prior-art.md`
and `docs/references.bib` hold the full record, with DOIs. In that file's
terms:

- **Established.** Quantitative analysis of a property on a polygonal
  molecular surface (areas, extrema, area-weighted statistics) is Bulat et al.
  `bulat2010` and Lu & Chen `lu2012mt`. Global surface-ESP statistics
  (positive and negative areas, mean, variance) go back to the GIPF approach
  of Murray et al. `murray1994gipf`. The ρ = 0.002 e·bohr⁻³ surface is Bader's
  `bader1987`, `bader1990`. The σ-hole is `clark2007`, `murray2009`,
  `politzer2013`.
- **Adapted.** Local, per-atom molecular surfaces and their ESP statistics are
  Multiwfn's quantitative molecular-surface analysis `lu2012multiwfn`,
  `lu2012mt`, `lu2024multiwfn`. Our nearest-nucleus and radius-scaled
  triangle partitions, and the reaction sites built from them, are versions
  of that idea. They were chosen to be simple to state and test, not because
  either is believed to be the correct partition (§2).
- **Possibly specific, unsearched.** Reporting both partitions with their
  difference in every record, carrying measured budget rows with their data
  files, and withholding site values behind a validity envelope are
  engineering choices. No prior-art search has been thorough enough to
  support calling them new, and this spec does not.

One correction the prior-art pass made, worth repeating where surfaces are
discussed: the ~98 % of the electron density often quoted for the "Bader
surface" is Bader's figure for the 0.001 au envelope. For the 0.002 surface
used here, `bench/isosurface_enclosed.jsonl` measures 98.5–99.1 %.

## 0. Prerequisites

| finding | what this needs from it | status |
| --- | --- | --- |
| supplied-H freezing | a pose's supplied hydrogens are not moved | **done**, #6 |
| H1 | PDB read through bond perception; `H_COUNT_MISMATCH` refuses a reader that guessed wrong | **done**, #7 (v0.1.1) |
| H2 | identity checked against `spec`, the coordinates PySCF receives (`verify_computed`); §1 extends it to atom maps | **done**, #10 |
| H3 | a one-signed field (every anion in the calibration set) gets a sequential colormap by default, never `vik`. Not needed by the site block; needed by Arrowfield, since a learner shown −80 kcal/mol in blue has been taught the wrong sign | **done**, #11, #12, #14 |
| H4 | `results.surface = {closed, open_edges, …, strict, problems}` and `input.gate = {verdict, forced, findings}` | **done**, #11 |
| H5 | planar 3D input kept; a single-element-pair molecule (nitrite) with z ≡ 0 is refused as ambiguous rather than guessed. Must land before any calibration fixture | **done**, #17 |
| H6 | unit-suffixed keys, `analysis_key` / `render_key` (§5), dipole with origin, σ-hole, `surface.definition`, and the `validity` envelope (§6) | **done**, #18 |
| H7 | `settings.solvent_esp` naming the PCM convention, part of `analysis_key` | **done**, #19 (both conventions measured; `--include-reaction-field`) |
| H14 | `geometry = {xyz_angstrom, sha1}` and full version provenance. A site's `V_min_at_A` means nothing without its coordinates | **done**, #20 |

## 1. Site identity and atom maps

Two levels, and only the first is a partition:

| level | what it is | who defines it |
| --- | --- | --- |
| `atoms` | every atom, hydrogens included, with its own share of the surface | the partition rule (§2) |
| `reaction_sites` | a named group of atoms; its values are aggregated from its atoms' triangles | the caller, or the default below |

**Why groups are needed.** Over an sp³ carbon almost all of the ρ = 0.002
surface lies outside the hydrogens `[I]`, so a chain carbon such as C3 gets
nearly zero area of its own. A C–H oxidation or desaturation question is about
C3 + H3a + H3b, not about the carbon nucleus. But "carbon plus its hydrogens"
is not the right unit for everything: a carbonyl O, an enolate C, or an
O–C=O group are sites too. So espmap does not fix one chemistry-shaped
definition.

- **Default groups**, used when the caller gives none: one group per heavy
  atom with its attached hydrogens, labelled by element and index. This keeps
  the CLI output meaningful on its own.
- **Caller groups** (`--sites groups.json`): a list of `{label, atoms}`, where
  `atoms` are atom-map numbers if the input is mapped, otherwise heavy-atom
  indices; hydrogens are included by naming their parent with `"+H"`
  (`{"label": "C3-H", "atoms": ["3+H"]}`). espmap checks only that every named
  atom exists. It does not know what an enolate is. Arrowfield's lessons
  define the groups.
- Groups may overlap (the carbonyl O inside a `CO2` group and on its own). Each
  group's values are computed from its own triangles, so overlap
  double-counts nothing.

Rules for both levels:

- A heavy atom is keyed by `index`, the 0-based heavy-atom index in input order.
  `[V]` `Chem.AddHs` appends hydrogens after the heavy atoms and
  `_mol_to_spec` writes atoms in RDKit index order, so RDKit index i is PySCF
  atom i for every heavy atom.
- `atom_map` is the input's atom-map number, or `null` if the input has none.
  It is metadata: it never changes what is computed.
- **Maps are all or nothing on heavy atoms.** Some heavy atoms mapped and
  others not is a gate BLOCK (`PARTIAL_ATOM_MAP`, remedy "map every heavy atom
  or none"). A site table with holes is worse than no maps.
- A repeated map number is a gate BLOCK (`DUPLICATE_ATOM_MAP`).
- **A map number on a hydrogen is a gate BLOCK in v1** (`MAPPED_HYDROGEN`).
  This is a limitation of v1, not a principle. The two hydrogens on one carbon
  are not always interchangeable: pro-R versus pro-S abstraction in P450
  oxidation or stereospecific desaturation is a question Arrowfield will want
  to ask. So every hydrogen already carries `index`, `parent`, coordinates and
  `h_source` (§4, §6), which is what a later schema needs to promote a
  hydrogen to a persistent identity. Until then hydrogens belong to their
  parent's site.
- A buried site (area 0 under the partition) is **still reported**, with
  `area_A2: 0` and `null` potentials. Dropping it would make "this carbon is not
  on the surface" indistinguishable from "this carbon was not computed".

### Identity checks, from the coordinates

`verify_computed` (#10) already checks the chemical identity against `spec`
`[V]`: connectivity, hydrogens included, perceived from distances and matched
to the approval's skeleton; stereo assigned from 3D and accepted if it
completes what the approval specified. It deliberately does not compare bond
orders, because the SCF never sees them and comparing them refused valid
resonance forms (`[N-]=C=S`, DMSO). Atom maps add one check:

| check | compares | must hold |
| --- | --- | --- |
| `chemical` | `verify_computed`, maps stripped from the approval first | as today |
| `mapped` | the **identity correspondence**: heavy atom i in `spec` against heavy atom i in the approval, by element and by the set of heavy neighbours' indices; each hydrogen against its parent's index | i → i is itself a match of the skeleton |

`verify_computed` searches for *some* match, so on a molecule with
graph-equivalent atoms (the two methyls of an isopropyl) it would accept a
geometry in which they were swapped. In a pose they are not equivalent: map 5
must point at the coordinates the adapter put there. `mapped` therefore checks
the correspondence the maps rely on, not just that one exists. It is the same
signature as "map 3 is the same graph node after the round trip": element and
mapped neighbours. Bond orders are left out for the same reason as in
`chemical`. Both results go into the JSON. A failure refuses, as invariant 3
does.

Measured `[M]` (RDKit 2026.09.1, the `espmap` venv):
- Atom maps change canonical SMILES: `[CH3:1][CH2:2][OH:3]` and `CCO` are
  different strings, and renumbering the maps even changes the atom order of
  the output. `chemical` must strip maps on both sides.
- RDKit writes maps to the V2000 atom block's mapping column and reads them
  back from that column alone, with no `molAtomMapNumber` property needed.
  Arrowfield's adapters can write maps there. Hydrogens added by `AddHs` get
  map 0.
- `verify_computed` passes a mapped approval against its own geometry as it
  stands; its `computed_smiles` then carries the maps.

## 2. Surface partition

**Partition triangles, not vertices.** Each triangle goes to one atom, chosen
by its centroid. Atom areas then sum exactly to the total area, and no area is
counted twice across a boundary. **This is an invariant**, tested numerically
under both rules (§7, test 1). Reaction sites are unions of their atoms'
triangles.

Two rules, both computed in **every** run. Partitioning is post-processing on
a surface that already exists, so the second rule costs nothing:

| rule | assigns a triangle to | role |
| --- | --- | --- |
| `nearest_nucleus` | argminᵢ dᵢ | primary: transparent and easy to explain |
| `radius_scaled` | argminᵢ dᵢ / rᵢ | comparison |

rᵢ is the Bondi van der Waals radius, with Mantina et al. (2009) for B, which
Bondi does not list. Every element in `TESTED_ELEMENTS` must have a radius. An
element outside it gets `radius_scaled: null` with a note, not a guessed
radius.

**Both rules' values and their difference are reported for every atom and
every site** (§5). That makes the partition an observable uncertainty axis,
measured on every molecule, rather than a hidden implementation choice or a
one-off benchmark. Neither rule is claimed to be the chemically correct one.

**Measure before refining.** No exact weighted-Voronoi clipping of triangles
in v1. If the two rules agree closely on the calibration set and the C18
acetate, that settles it. If they differ materially, that difference is the
finding, and a finer partition is investigated then.

On a folded surface (genus > 0) a triangle can be nearest to an atom on the
other side of a contact. That is correct, and it is how the contact shows up
in the site table, but it makes those sites not comparable with an extended
conformer's (§5, `topology`).

Before anything is quoted against it, add an error-budget row for the
partition, with its data committed under `bench/`. This is the per-site
analogue of the "σ-hole region definition" row, and unlike that row it must
have a data file (H8).

## 3. Area weighting and statistics

Use exactly what `bench/conformers.py:surface_stats` already does `[V]`:
triangle area × mean V over the triangle's three vertices. Every per-site mean
and every area is area-weighted. **Never average over vertices.**
Marching-cubes vertex density follows grid orientation and curvature, not
chemistry.

Per atom and per site, under each rule:

| field | what | robust? |
| --- | --- | --- |
| `V_area_mean_kcal_mol` | area-weighted mean | yes |
| `V_p05_kcal_mol`, `V_p95_kcal_mol` | area-weighted 5th and 95th percentiles of the per-triangle V | yes: the "reactive patch" value |
| `V_min_kcal_mol`, `V_max_kcal_mol` | extrema over the site's triangles' vertices, with position in Å | **no** |

Extrema are kept because V is computed exactly at the vertices and the global
extrema are the numbers the error budget is about. On a small site, though,
an extremum can be one vertex on a boundary, and it moves with mesh
orientation and small geometric changes. **Arrowfield leads with the mean and
the percentiles, and shows the extremum as a raw value.** A vertex on a
boundary belongs to triangles of two sites and counts for both sites'
extrema.

`[U]` Whether 5 % is the right patch size is not measured: on a site of a few
Å² it is a handful of triangles. Measure the percentiles' grid sensitivity
alongside the areas (§7, test 10) and change the cut-off then if needed.

**No site block on an invalid record.** If `validity.surface` or
`validity.input` is false (§6), the record holds `"sites": null` and the
reason. A clipped surface moved one extremum by 23 % with every other number
plausible (H4), and per-site values near the box face would be worse than the
global ones.

`[U]` How far per-site areas move with grid spacing is not measured. The
existing 0.04 kcal/mol spacing row is for V_max only, and it has no data file
(H8). Measure 0.26 → 0.16 bohr on site areas, commit the rows, and only then
quote an area difference.

## 4. Geometry and hydrogen provenance

`from_file` freezes every atom the file supplied, hydrogens included, and
relaxes only the hydrogens `AddHs` created (#6) `[V]`. Before #6 the supplied
hydrogens moved too: a test pose's O–H moved 0.79 Å. If MMFF94s has no
parameters, the added hydrogens keep their `addCoords` positions and the record
says so. For PDB and XYZ, `H_COUNT_MISMATCH` (#7) refuses a file whose
perceived bond orders would need more hydrogens than it supplied.

Each record carries a summary:

```json
"geometry_provenance": {
  "heavy_atoms": "supplied_frozen",
  "hydrogens": {"supplied": 2, "generated": 34, "generated_relaxed": 34},
  "relaxation": "MMFF94s, generated hydrogens only, isolated molecule",
  "caveat": "heavy-atom pose from the input; generated hydrogens placed without any protein environment"
}
```

and every hydrogen in `atoms` carries `h_source`: `"supplied"`,
`"generated"` or `"generated_relaxed"`. `heavy_atoms` is `"generated"` on the
SMILES path.

**Hydrogen sensitivity is flagged per site.** A site containing a polar
hydrogen (on O, N or S) that was generated, not supplied, gets
`h_generated_polar: true`. `[I]` For a terminal –CH₂–OH, the O–H rotamer in a
pocket is set by a protein hydrogen bond, and isolated relaxation will pick a
different one. The ESP around O is sensitive to it. espmap records the fact, so
Arrowfield can say "the heavy-atom pose came from the protein model; this O–H
was placed without the protein". Sampling O–H rotamers is a later, separate
option.

## 5. Uncertainty and comparability

espmap cannot put an error bar on each site from a single run. It can report
what was measured, and it must point to what was not.

- **Partition sensitivity**, per atom and per site: every quantity under
  `nearest_nucleus` and under `radius_scaled`, and `delta` (radius-scaled
  minus nearest-nucleus). Measured in this run. `delta` is `null` for a
  quantity that is `null` under either rule. Reporting both values and not
  just the delta lets a reader see a site whose area halves between rules
  while its mean barely moves.
- `budget` at run level: see §6. A measured shift is evidence about how
  sensitive the method is, recorded with the molecule it was measured on. The
  record does not turn chlorobenzene's +3.54 basis shift into a ±3.54
  interval on every site.
- `conformer_scope`: `"single conformer"` together with the rotatable-bond
  count, plus the measured warning when a folded geometry is possible (C18
  acetate, V_min spread 1.1 → 5.2 kcal/mol once folded members are included)
  `[M]`, `bench/conformers.jsonl`, `bench/f9_recheck.jsonl`.
- `topology`: genus and number of pieces, copied from `results.surface`.
  Genus > 0 means the molecule touches itself, and sites near the contact are
  not comparable with an extended conformer's.

**Two keys, not one.** Review 3 (H6) proposed a single `compare_key` over the
method *and* the colour settings. That would make two identical calculations
"not comparable" because one was drawn with another colormap.

| key | hashes | answers |
| --- | --- | --- |
| `analysis_key` | `xc`, `basis`, `solvent`, `solvent_esp` (H7), `iso`, `spacing`, `density_fit`, surface definition; names normalised first (lower-case, `vacuum` → null) | are these numbers methodologically comparable? |
| `render_key` | `cmap`, `vmin`/`vmax` (or `lim`), scale kind, camera, image size, lighting | do these pictures share one visual scale? |

**Schema versions are not in `analysis_key`.** Revision 3 put the partition
and site schema versions there. That would make every record written before
P2 stop matching every record written after it, even for global V_min and
V_max, which a partition cannot change. Instead the sites block carries its
own `schema` (`espmap.sites/1`, §6).

Arrowfield compares global values between records whose `analysis_key`
matches, site values only when `sites.schema` matches too, and shows images
side by side only when `render_key` also matches. `scale.comparable` depended
on the colour flags alone; H6 (#18) renames it `colour_scale_fixed` and moves
the comparison question to these keys.

## 6. JSON schema (v1)

The record as a whole moves towards one versioned shape, `espmap.record/1`:

```
identity · geometry · calculation · validity · surface · properties · sites · uncertainty · provenance
```

The existing `input` / `settings` / `results` keys are mapped into it when H4,
H6 and H14 land. This spec fixes the parts the site block reads.

**Validity first.** A consumer reads one field and never infers validity from
stderr, χ, a filename or console text:

```json
"validity": {"input": true, "surface": true, "sites": true},
"gate": {"verdict": "ACCEPT", "forced": false, "findings": []},
"surface": {
  "definition": {"property": "electron_density", "isovalue_e_bohr3": 0.002},
  "closed": true, "open_edges": 0, "nonmanifold_edges": 0,
  "components": 1, "genus": 0, "strict": true, "problems": []
}
```

- `validity.input` is false if the gate verdict was BLOCK (forced) or
  `structure_verified` is false.
- `validity.surface` is false if `surface.problems` is non-empty.
- `validity.sites` is false if either is, or the site block failed its own
  checks. `"sites"` is then `null` with a reason.
- `surface.definition` states that every area here is on the ρ = 0.002 e·bohr⁻³
  electron-density isosurface. **It is not a solvent-accessible surface.** A
  large site area means much of that atom's density boundary is exposed; it
  does not mean the atom is accessible to a protein or a solvent molecule,
  which is a question for Arrowfield's geometry layer.

Arrowfield's rule is then `if not record["validity"]["sites"]: refuse`, with
no chemistry-specific defensive logic of its own.

The site block:

```json
"sites": {
  "schema": "espmap.sites/1",
  "identity": {
    "chemical_smiles": "<computed, maps stripped>",
    "mapped_smiles": "<computed, maps kept>",
    "chemical_match": true,
    "mapped_match": true
  },
  "partition": {"primary": "nearest_nucleus", "compared": "radius_scaled",
                "unit": "triangle", "radii": "Bondi; Mantina 2009 for B"},
  "geometry_provenance": { "...": "see §4" },
  "solvent_esp": "<copied from settings.solvent_esp (H7)>",
  "atoms": [
    {
      "index": 0,
      "atom_map": 1,
      "element": "O",
      "parent": null,
      "h_source": null,
      "nearest_nucleus": {
        "area_A2": 0.0, "area_pos_A2": 0.0, "area_neg_A2": 0.0,
        "V_area_mean_kcal_mol": 0.0,
        "V_p05_kcal_mol": 0.0, "V_p95_kcal_mol": 0.0,
        "V_min_kcal_mol": 0.0, "V_min_at_A": [0.0, 0.0, 0.0],
        "V_max_kcal_mol": 0.0, "V_max_at_A": [0.0, 0.0, 0.0]
      },
      "radius_scaled": { "...": "same fields" },
      "delta": {"area_A2": 0.0, "V_area_mean_kcal_mol": 0.0,
                "V_p05_kcal_mol": 0.0, "V_p95_kcal_mol": 0.0,
                "V_min_kcal_mol": 0.0, "V_max_kcal_mol": 0.0}
    }
  ],
  "reaction_sites": [
    {
      "label": "O1",
      "source": "default",
      "atoms": [0],
      "h_generated_polar": false,
      "nearest_nucleus": { "...": "same fields as an atom" },
      "radius_scaled": { "...": "same fields" },
      "delta": { "...": "same fields" }
    }
  ],
  "conformer_scope": "single conformer; 0 rotatable bonds",
  "topology": {"genus": 0, "components": 1},
  "budget": [
    {"axis": "basis def2-SVP -> def2-TZVP", "V_max_kcal_mol": [1.53, 2.82],
     "measured_on": "decanoic acid, BDSF, Z3Z13-18:OAc, Z3Z13-18:OH",
     "source": "bench/results.jsonl"},
    {"axis": "solvent gas -> PCM water, polarised solute", "V_max_kcal_mol": -3.60,
     "measured_on": "chlorobenzene", "source": "bench/solvent_convention.jsonl"},
    {"axis": "partition rule, site means", "kcal_mol": 3.95,
     "measured_on": "Z3Z13-18:OAc, 10 conformers", "source": "bench/site_partition.jsonl"}
  ]
}
```

The zeros are **placeholders for type and shape only**. Arrowfield's v0.0
fixtures use this shape with values marked as fabricated until real runs
replace them.

`parent` is the heavy-atom index for a hydrogen, `null` otherwise; `h_source`
is set for hydrogens only. Coordinates come from `geometry.xyz_angstrom` by
`index`. `source` on a site is `"default"` or `"caller"`.

**Budget rows.** Each row names its source: a `bench/` file, or `"prose"` for
a figure with no data behind it. Since H8 (#25) every row of the σ-hole
budget has a `bench/` file and `test_claims.py` recomputes each one from it;
the per-site rows (partition, conformer, grid) do the same from
`bench/site_partition.jsonl`. Arrowfield must display a `"prose"` row as
unbacked. Rows are selected by the run's settings: the basis row appears for
a def2-SVP run, the solvent row for a gas-phase run, and so on, and every
solvent row names its convention (`solvent_esp`). The rows above illustrate
the shape with the current measured values.

**NPZ.** `--npz` gains `tri_atom` and `tri_atom_compared`: the atom index of
every triangle under each rule. Sites are unions of atoms, so this is enough
to colour a render by site, to regroup with other sites, or for Arrowfield to
redo a partition, all without another SCF. Units and method go in the NPZ as
H6 requires.

## 7. Acceptance tests

Fast (no SCF, on synthetic meshes like the existing `mesh_topology` tests, or
on RDKit alone):

1. **Conservation.** Atom areas sum to the total area (rel. 1e-12) under both
   rules, and the default sites partition the surface the same way.
2. Renumbering atom maps leaves `chemical_smiles` unchanged, changes
   `mapped_smiles`, and permutes `atoms` and `reaction_sites` to match. After
   an SDF write and read, each map number has the same neighbourhood
   signature as before.
3. Duplicate, partial and hydrogen maps are each a gate BLOCK. A multi-record
   SDF is a BLOCK without `--record`, and accepted with it.
4. **Equivalent-atom swap detector.** Swap the coordinates of two
   graph-equivalent heavy atoms (the methyls of an isopropyl) in `spec` after
   the gate: `chemical_match` stays true, `mapped_match` must be false and the
   run refused. This is the case `verify_computed` alone cannot see.
5. Heavy atoms and supplied hydrogens move < 1e-6 Å through hydrogen
   addition; generated hydrogens are the only atoms that move (#6 has the
   O–H case); `h_source` is right for every hydrogen.
6. A buried atom appears with area 0 and null potentials.
6b. A caller group naming an atom that does not exist is refused before the
   SCF. Overlapping groups each get their own values.
7. `validity.sites` false and `"sites": null` with a reason on a surface with
   `problems`, on a forced BLOCK, and on `structure_verified: false`.
8. The minimum of the site `V_min` values equals the global `V_min`, and the
   same for `V_max` (exact; they are the same vertices).
8b. Two records that differ only in `--cmap` have equal `analysis_key` and
   different `render_key`; two that differ only in `--basis` the reverse.

Slow (real DFT):

9. **Symmetry detector.** Sites related by symmetry agree: the two O of
   nitrite, the two O of acetate, the six C of benzene. Use geometries that are
   symmetric by construction (an explicit `--atoms` string, or
   `monosubstituted_benzene("H", …)`), not force-field output, so any
   asymmetry comes from the partition or the mesh and not from the geometry.
   Set the tolerance from the observed mesh noise and record it. This plays
   the role that far-field decay plays for the kernel sign (invariant 4): a
   partition bug gives plausible numbers, but it cannot give equal numbers to
   equivalent atoms by accident.
10. Grid-spacing convergence of site areas and percentiles, 0.26 → 0.16 bohr,
    committed as rows under `bench/`.
11. Partition sensitivity on the C18 acetate, extended and folded conformers,
    committed under `bench/`.

## 8. Implementation order

Each step lands with its tests, as a separate change.

**Done:** supplied-hydrogen freezing (#6), H1 (#7, v0.1.1), H2 (#10).

**P0, boundaries.** Make one record impossible to misread:
1. H3: the colormap follows the scale kind.
2. H4: `gate`, `surface` and the `validity` envelope.
3. H5: planar 3D input kept (before any calibration fixture).
4. H6 and H7: unit-suffixed keys, `analysis_key` / `render_key`,
   `solvent_esp`.
5. H14: geometry and full provenance in the record.
6. One record per run: `MULTIPLE_RECORDS` and `--record N`.

**P1, identity.** Atom-map gate rules, the `mapped` check, SDF round-trip
tests (§1; tests 2–4).

**P2, sites.** Triangle partition under both rules (§2–3; tests 1, 6–8), then
reaction-site grouping (test 6b), then the slow tests 9–11, with every
measurement committed under `bench/`.

**P3, the contract.** Freeze `espmap.sites/1` and produce real fixture records
for benzene, nitrite, thiocyanate, an enolate model and a C18 pheromone.
**Arrowfield v0.0 starts here.**

**P4, more descriptors.** ALIE, then frontier-orbital amplitudes, then Fukui
functions, each reusing the same geometry, surface, atom maps, partition and
reaction sites, so a site gains a field per descriptor rather than a new
integration model.

## When this lands

Adding a partition is new science under "If you change the science" in
`CLAUDE.md`. The partition row and the grid row for site areas go into the
error-budget tables in `CLAUDE.md` and `README.md`, with their `bench/` files,
in the same change. Invariant 3 in `CLAUDE.md` (rewritten for
`verify_computed` in #10) gains the `mapped` check in P1. "One geometry in,
one record out" becomes a new invariant in P0.

## Not in this spec

- How ALIE, orbital amplitudes and Fukui functions are computed. They reuse
  this partition and this `atoms` / `reaction_sites` shape (P4), which is why
  the partition has to be settled first.
- Anything about AF3, GNINA or licensing. Whether AlphaFold Server output may
  feed docking matters to Arrowfield's adapters, not here; espmap only sees an
  SDF. Those claims came from an external summary and are not verified in
  this repo.
- Ensembles. One run means one geometry. Spreads across conformers or poses
  are computed by whoever runs the ensemble, from these records.
