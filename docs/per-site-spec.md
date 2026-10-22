# Per-site output: specification

Status: **draft, revision 2, nothing implemented.** This is the contract that
espmap's `--json` will grow into and that Arrowfield's v0.0 fixtures are
written against. It extends item 1 of `docs/roadmap-review.md` ("extend
`--json` with shape descriptors") from whole-molecule numbers to per-atom ones.

Labels as in the roadmap review: `[V]` verified in the repo · `[M]` measured ·
`[I]` inferred · `[U]` not verifiable yet, a test must settle it.

**Revision 2 (2026-10-22).** Rewritten against `docs/adversarial-review-3.md`,
which landed while revision 1 was being drafted. Revision 1 built on four
things that review shows are broken:
- The identity split rested on the invariant-3 assertion, which compares the
  input with itself (H2), so the "mapped identity" check would have been
  tautological too. Both checks now perceive from the coordinates sent to
  PySCF (§1).
- Per-site numbers would have been written for a clipped surface with nothing
  in the record to say so (H4). They are now withheld (§3).
- The `budget` block quoted error-budget rows that have no data in `bench/`
  (H8). It now carries a source for each row (§6).
- Site positions are meaningless without the geometry, which the JSON does not
  hold (H14). The geometry is now a prerequisite (§0).

Other changes: two levels, atoms and reaction sites, instead of one
carbon-centric site definition (§1); values under both partition rules, not
only their difference (§5); an implementation order (§8); partial and
hydrogen atom maps are refused; buried sites,
boundary vertices and radii for every tested element are defined; a
hydrogen-sensitivity flag, an NPZ field, two consistency tests and symmetric
test geometries are added; the example values no longer reuse a measured
number.

## Boundary

espmap decides what a defensible site descriptor is. Arrowfield decides how to
teach with it. espmap never reads AF3 or GNINA output, never knows a protein is
involved, and never sees a docking score. It receives one molecule with 3D
coordinates and returns numbers per site.

**Input for a pose is SDF (or MOL) with explicit bond orders.** Not PDB: a PDB
ligand without CONECT records is read with every bond single and computes a
different molecule while every check passes (review 3, H1). Not XYZ: bond
orders are perceived, not given, so the mapped structure Arrowfield intended
is not guaranteed. Arrowfield's adapters assign bond orders from the mapped
ligand definition before writing the SDF.

## 0. Prerequisites from review 3

The site block is built on the record and the checks below. None of them
exists in a usable form today, so they land first, in this order:

| finding | what the site block needs from it |
| --- | --- |
| H1 | the PDB reader fix and the `H_COUNT_MISMATCH` guard, so a pose with missing bond orders is refused, not silently recomputed as another molecule |
| H2 | structure perceived from `spec` (the coordinates PySCF receives), with stereo from 3D. §1 extends this to atom maps |
| H4 | `results.surface = {closed, open_edges, …, strict, problems}` and `input.gate = {verdict, forced, findings}` |
| H5 | planar 3D input kept. Nitrite and thiocyanate, two of the three calibration anions, are planar or linear, and a pose written in the xy plane would be discarded and regenerated |
| H6 | unit-suffixed keys and `compare_key`. The site block follows the same naming |
| H7 | `settings.solvent_esp` naming the PCM convention. Every site potential under `--solvent` is the polarised solute only |
| H14 | `geometry = {xyz_angstrom, sha1}`. A site's `V_min_at_A` means nothing without the coordinates it refers to |

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
- **A map number on a hydrogen is a gate BLOCK** (`MAPPED_HYDROGEN`).
  Hydrogens belong to their parent's site and carry no identity of their own.
  This also avoids depending on whether `RemoveHs` keeps mapped hydrogens
  `[U]`.
- A repeated map number is a gate BLOCK (`DUPLICATE_ATOM_MAP`).
- A buried site (area 0 under the partition) is **still reported**, with
  `area_A2: 0` and `null` potentials. Dropping it would make "this carbon is not
  on the surface" indistinguishable from "this carbon was not computed".

### Identity checks, perceived from the coordinates

H2's correction perceives the computed structure from `spec`: XYZ block,
`DetermineBonds(charge)`, `AssignStereochemistryFrom3D`. Because of the index
correspondence above, atom i of that perceived molecule is input atom i, so
maps can be transferred onto it by index. Two checks follow, both against the
**perceived** molecule, never the input:

| check | compares | must hold |
| --- | --- | --- |
| `chemical` | perceived canonical isomeric SMILES, maps stripped, against the gate-approved SMILES, maps stripped | equal (stereo-stripped where the approved structure leaves stereo unspecified, as H2 specifies) |
| `mapped` | for every heavy atom: (map, element, sorted maps of heavy neighbours, bond orders to them) in the perceived molecule against the input | identical for every atom |

`chemical` alone would pass a pose where two atoms were swapped between the
reader and the coordinates. `mapped` alone would pass nothing extra on an
unmapped input. `mapped` is a neighbourhood signature, so it also covers
map-to-chemistry consistency: map 3 after the round trip must be the same
graph node (same element, same mapped neighbours, same bond orders), not just
some unique atom that carries the number 3. Both are written to the JSON with their result. A failure is a
refusal, the same as invariant 3.

`[I]` Atom maps appear in RDKit's SMILES output, so a mapped SMILES and its
unmapped form would compare as different strings. That is why `chemical`
strips maps on both sides. Confirm by test.

`[U]` Whether `SDMolSupplier` reads atom-map numbers from a V2000 atom block's
mapping column, or only from a `molAtomMapNumber` property, decides how
Arrowfield must write its SDFs. Test this before writing the adapter.

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

On a folded surface (genus > 0) a triangle can be nearest to an atom on the
other side of a contact. That is correct, and it is how the contact shows up
in the site table, but it makes those sites not comparable with an extended
conformer's (§5, `topology`).

Before anything is quoted against it, add an error-budget row for the
partition, with its data committed under `bench/`, measured on the
calibration set and the C18 acetate. This is the per-site analogue of the
"σ-hole region definition" row, and unlike that row it must have a data file
(H8).

## 3. Area weighting and extrema

Use exactly what `bench/conformers.py:surface_stats` already does `[V]`:
triangle area × mean V over the triangle's three vertices. Every per-site mean
and every area is area-weighted. **Never average over vertices.**
Marching-cubes vertex density follows grid orientation and curvature, not
chemistry.

Extrema per site (`V_min_kcal_mol`, `V_max_kcal_mol`) are taken over the
vertices of the triangles assigned to that site, because V is computed exactly
at the vertices. They are not area-weighted: an extremum is a point property.
A vertex on a boundary belongs to triangles of two sites and counts for both.
Each extremum also reports its position in Å.

**No site block on a damaged surface.** If `results.surface.problems` is
non-empty (a `--no-strict` run), or the gate verdict was BLOCK and `--force`
was used, the record holds `"sites": null` and the reason. A clipped surface
moved one extremum by 23 % with every other number plausible (H4), and per-site
values near the box face would be worse than the global ones.

`[U]` How far per-site areas move with grid spacing is not measured. The
existing 0.04 kcal/mol spacing row is for V_max only, and it has no data file
(H8). Measure 0.26 → 0.16 bohr on site areas, commit the rows, and only then
quote an area difference.

## 4. Geometry and hydrogen provenance

**Today's behaviour** `[V]` (`from_file`, `esp_map.py:257-271`): when a 3D
file is missing some hydrogens, they are added with `addCoords=True`, every
heavy atom is fixed, and MMFF94s minimises the rest. **The rest includes the
hydrogens the file did supply.** A GNINA pose with its polar hydrogens kept
would therefore have its O–H hydrogen re-optimised along with the added C–H
hydrogens. That is the exact case where protein context mattered.

**Required change** (PR #6): freeze every supplied atom, hydrogens included,
and relax only the hydrogens that were added. If MMFF94s has no parameters, the added
hydrogens keep their `addCoords` positions and the record says so. H1's
`H_COUNT_MISMATCH` guard runs first: if the file has explicit hydrogens and
`AddHs` wants to add more to an atom that already has them, the bond orders
are wrong and nothing is relaxed.

Each record carries:

```json
"geometry_provenance": {
  "heavy_atoms": "supplied_frozen",
  "hydrogens": {"supplied": 2, "generated": 34, "generated_relaxed": 34},
  "relaxation": "MMFF94s, generated hydrogens only, isolated molecule",
  "caveat": "heavy-atom pose from the input; generated hydrogens placed without any protein environment"
}
```

`heavy_atoms` is `"generated"` on the SMILES path. `generated_relaxed` is a
subset of `generated`.

**Hydrogen sensitivity is flagged per site.** A site containing a polar hydrogen
(on O, N or S) that was generated, not supplied, gets
`h_generated_polar: true`. `[I]` For a
terminal –CH₂–OH, the O–H rotamer in a pocket is set by a protein hydrogen
bond, and isolated relaxation will pick a different one. The ESP around O is
sensitive to it. espmap records the fact. Arrowfield decides what to say about
it. Sampling O–H rotamers is a later, separate option.

## 5. Uncertainty fields

espmap cannot put an error bar on each site from a single run. It can report
what was measured, and it must point to what was not.

- **Partition sensitivity**, per atom and per site: every quantity under
  `nearest_nucleus` and under `radius_scaled`, and `delta` (radius-scaled
  minus nearest-nucleus). Measured in this run. `delta` is `null` for a
  quantity that is `null` under either rule. Reporting both values and not
  just the delta lets a reader see a site whose area halves between rules
  while its mean barely moves.
- `budget` at run level: see §6. Rows come with their source; the record does
  not attach ±3.55 to every site.
- `conformer_scope`: `"single conformer"` together with the rotatable-bond
  count, plus the measured warning when a folded geometry is possible (C18
  acetate, V_min spread 1.1 → 5.2 kcal/mol once folded members are included)
  `[M]`, `bench/conformers.jsonl`, `bench/f9_recheck.jsonl`.
- `topology`: genus and number of pieces, copied from `results.surface`.
  Genus > 0 means the molecule touches itself, and sites near the contact are
  not comparable with an extended conformer's.

Arrowfield, not espmap, compares two runs' site values against these numbers,
and only between runs whose `compare_key` (H6) matches.

## 6. JSON schema (v1)

Added to the existing `--json` record as a new top-level key, so the existing
fields change only as review 3 already requires.

```json
"sites": {
  "schema": "espmap.sites/1",
  "identity": {
    "chemical_smiles": "<perceived, maps stripped>",
    "mapped_smiles": "<perceived, maps transferred by index>",
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
      "nearest_nucleus": {
        "area_A2": 0.0, "area_pos_A2": 0.0, "area_neg_A2": 0.0,
        "V_area_mean_kcal_mol": 0.0,
        "V_min_kcal_mol": 0.0, "V_min_at_A": [0.0, 0.0, 0.0],
        "V_max_kcal_mol": 0.0, "V_max_at_A": [0.0, 0.0, 0.0]
      },
      "radius_scaled": { "...": "same fields" },
      "delta": {"area_A2": 0.0, "V_area_mean_kcal_mol": 0.0,
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
    {"axis": "basis def2-SVP -> def2-TZVP", "V_max_kcal_mol": [1.53, 2.83],
     "measured_on": "decanoic acid, BDSF, Z3Z13-18:OAc, Z3Z13-18:OH",
     "source": "bench/results.jsonl"},
    {"axis": "solvent gas -> PCM water", "V_max_kcal_mol": -3.61,
     "measured_on": "chlorobenzene", "source": "prose"}
  ]
}
```

The zeros are **placeholders for type and shape only**. Arrowfield's v0.0
fixtures use this shape with values marked as fabricated until real runs
replace them.

**Budget rows.** Each row names its source: a `bench/` file, or `"prose"` for
the seven rows review 3 found with no data (H8). Arrowfield must display a
`"prose"` row as unbacked. Rows are selected by the run's settings: the basis
row appears for a def2-SVP run, the solvent row for a gas-phase run, and so
on. The numbers come from one table in the module, and the claim-consistency
test checks that table against `CLAUDE.md` and `README.md`. The two rows above
are illustrations of the shape; the basis range is the one H8 recomputed.

`parent` is the heavy-atom index for a hydrogen, `null` otherwise. `source`
is `"default"` or `"caller"`.

**NPZ.** `--npz` gains `tri_atom` and `tri_atom_compared`: the atom index of
every triangle under each rule. Sites are unions of atoms, so this is enough
to colour a render by site, to regroup with other sites, or for Arrowfield to
redo a partition, all without another SCF. Units and method go in
the NPZ as H6 requires.

## 7. Acceptance tests

Fast (no SCF, on synthetic meshes like the existing `mesh_topology` tests, or
on RDKit alone):

1. **Conservation.** Atom areas sum to the total area (rel. 1e-12) under both
   rules, and the default sites partition the surface the same way.
2. Renumbering atom maps leaves `chemical_smiles` unchanged, changes
   `mapped_smiles`, and permutes `atoms` and `reaction_sites` to match. After
   an SDF write and read, each map number has the same neighbourhood
   signature as before.
3. Duplicate, partial and hydrogen maps are each a gate BLOCK.
4. **Index swap detector.** Swap two non-equivalent heavy atoms' coordinates in
   `spec` after the gate: `mapped_match` must be false and the run refused.
   This is the H2 lesson applied to maps. A check that never reads `spec`
   cannot catch it.
5. Heavy atoms and supplied hydrogens move < 1e-6 Å through hydrogen
   addition. Generated hydrogens are the only atoms that move.
6. A buried atom appears with area 0 and null potentials.
6b. A caller group naming an atom that does not exist is refused before the
   SCF. Overlapping groups each get their own values.
7. `"sites": null` with a reason on a surface with `problems`, and on a forced
   BLOCK.
8. The minimum of the site `V_min` values equals the global `V_min`, and the
   same for `V_max` (exact; they are the same vertices).

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
10. Grid-spacing convergence of site areas, 0.26 → 0.16 bohr, committed as
    rows under `bench/`.

## 8. Implementation order

Each step lands with its tests, as a separate change:

1. **Supplied-hydrogen freezing:** a regression test first (it must fail on
   today's code), then the fix: remember which atoms existed before
   `AddHs()`, freeze all of them, relax only the new ones. It depends on
   nothing else, and it is a correctness bug today, with or without
   Arrowfield. **In review: PR #6** (the test fails on master with the
   supplied O–H moving 0.79 Å, and passes with the fix).
2. **Review 3's H1 and H2:** the PDB reader, `H_COUNT_MISMATCH`, and identity
   perceived from `spec`. Atom-map round-trip checks extend H2's perceived
   molecule, so they cannot come first.
3. **Atom-map gate rules and the `mapped` check** (§1), tests 2–4.
4. **Record fields from review 3** the site block reads: H4 (surface,
   gate), H14 (geometry), H6 (units, `compare_key`), H7 (`solvent_esp`).
5. **Triangle partition under both rules** (§2–3), tests 1, 6–8.
6. **Reaction-site grouping** above the atom level, default and caller
   groups, test 6b.
7. **Slow tests:** symmetry on constructed geometries, grid convergence of
   areas, and the partition budget row, committed under `bench/`.
8. **JSON schema and fixture outputs** from real runs on the calibration set.
9. Only then: Arrowfield v0.0 scaffolded against those fixtures.

## When this lands

Adding a partition is new science under "If you change the science" in
`CLAUDE.md`. The partition row and the grid row for site areas go into the
error-budget tables in `CLAUDE.md` and `README.md`, with their `bench/` files,
in the same change. Invariant 3 in `CLAUDE.md` is reworded to the two
perceived checks of §1.

## Not in this spec

- ALIE, frontier-orbital amplitudes and Fukui functions. They will reuse the
  same partition and the same `atoms` / `reaction_sites` shape, with one extra field each,
  which is why the partition has to be settled first.
- Anything about AF3, GNINA or licensing. Whether AlphaFold Server output may
  feed docking matters to Arrowfield's adapters, not here; espmap only sees an
  SDF. Those claims came from an external summary and are not verified in
  this repo.
- Ensembles. One run means one geometry. Spreads across conformers or poses
  are computed by whoever runs the ensemble, from these records.
