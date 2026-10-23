# Prior art and intellectual provenance

espmap's development record (three adversarial reviews, `bench/` data, the
roadmap) documents where its *implementation* decisions came from. This file
does the same for the *chemistry and methods* it builds on, so that no reader
mistakes an established technique for an espmap invention. Every reference is
in [`references.bib`](references.bib) with a DOI, checked against Crossref or
the publisher on 2026-10-22; the few details that could not be checked are
marked `UNVERIFIED` there.

The rule this file enforces: **espmap claims no new theory.** Its
contributions are engineering, validation and measurement choices, listed
in the last section, and none of them has had a prior-art search thorough
enough to support a novelty claim.

## Established methods espmap uses

| concept | in espmap | prior work |
| --- | --- | --- |
| Electrostatic potential as a guide to reactivity and interactions | the whole program | Politzer & Murray 2002 `politzer2002` |
| ρ = 0.002 e·bohr⁻³ as the molecular surface | `ISO_BADER`, `extract_surface` | Bader et al. 1987 `bader1987`; Bader 1990 `bader1990` |
| Quantitative analysis of a property *on a polygonal molecular surface*: areas, extrema, statistics | `extract_surface` + `esp_at_points` | Bulat et al. 2010 `bulat2010`; Lu & Chen 2012 `lu2012mt` |
| Evaluating the potential at surface vertices rather than interpolating from a volumetric grid | `esp_at_points` on the mesh | Bulat et al. 2010 and Lu & Chen 2012 discuss the cost of evaluating properties at surface vertices, artificial extrema from uneven discretisation, and grid-spacing effects |
| The σ-hole and V_s,max along the R–X axis | `sigma_hole` | Clark et al. 2007 `clark2007` (the term); Murray, Lane & Politzer 2009 `murray2009`; Politzer, Murray & Clark 2013 `politzer2013` |
| Local / per-atom molecular surfaces and their ESP statistics | `partition_surface`, `site_report` (per-site records, #24) | Multiwfn's quantitative molecular surface analysis, Lu & Chen 2012 `lu2012multiwfn`, `lu2012mt`; Lu 2024 `lu2024multiwfn` |
| Global surface-ESP statistics: σ₊², σ₋² (each about its own sign's mean), σ²_tot, ν | `gipf_stats` → `results.gipf`, area-weighted over triangles where GIPF averaged over points. `area_stats`' `V_area_sd` is an SD about the overall mean, not σ²_tot | Murray et al. 1994, GIPF `murray1994gipf`; Bulat et al. 2010 |
| Positive/negative areas, area-weighted mean and percentiles | `area_stats` → `results.surface_stats` and every site (#24) | Bulat et al. 2010 `bulat2010`; Lu & Chen 2012 `lu2012mt` |
| Van der Waals radii for the radius-scaled partition | `BONDI`, `partition_surface(..., "radius_scaled")`; the record's `"radii": "Bondi; Mantina 2009 for B"` | Bondi 1964 `bondi1964`; Mantina et al. 2009 `mantina2009` (B, which Bondi does not list) |
| Assigning space to the nearest of weighted sites | the radius-scaled rule is a multiplicatively weighted Voronoi assignment (distance / r_Bondi) | Aurenhammer 1991 `aurenhammer1991` |
| Bonds and bond orders perceived from coordinates (xyz2mol) | RDKit's `rdDetermineBonds`, for XYZ, PDB, `--atoms` and every `verify_computed` | Kim & Kim 2015 `kim2015xyz2mol` |
| Expected bond lengths from covalent radii | RDKit's radii table, scaled for bond order, to ask whether a flat MOL/SDF file can be a real geometry (`coordinate_kind`, review 5 R11) | Cordero et al. 2008 `cordero2008` |
| Atropisomer axes classed by rotation barrier: below about 20 kcal/mol they interconvert within a minute and are not stereo for practical purposes | `--require-stereo` treats tri- and tetra-ortho biaryl axes (peri fusions weighted 2) as stable, an approximation from the graph alone (`atropisomer_axes`) | LaPlante et al. 2011 `laplante2011` |
| A positive HOMO marks an anion the basis set cannot bind | `unbound_anion_warning` → `results.electronic.problems` (`UNBOUND_ANION_GAS`) | Galbraith & Schaefer 1996 `galbraith1996`, with the comment by Rösch & Trickey 1997 `rosch1997`; Jensen 2010 `jensen2010` |
| Average local ionisation energy | P4, specified in `docs/p4-spec.md` (not implemented) | Sjoberg et al. 1990 `sjoberg1990`; Murray, Brinck & Politzer 1992 `murray1992`; Bulat et al. 2010 |
| Marching cubes | scikit-image `measure.marching_cubes` | Lorensen & Cline 1987 `lorensen1987`; Lewiner et al. 2003 `lewiner2003` |
| Perceptual colormaps; diverging maps neutral at zero | `build_maps`, the colormap policy | Crameri et al. 2020 `crameri2020`, `crameri_scm`; Moreland 2009 `moreland2009` |

The ~98 % often quoted for the density surface is Bader's figure for the
**0.001** au envelope (H atoms in hydrocarbons, `bader1990`). For the 0.002
surface espmap uses, `bench/isosurface_enclosed.jsonl` measures 98.5–99.1 %
of the electrons inside it for water, chlorobenzene, acetate and decanoic acid
(B3LYP/def2-TZVP, molecular quadrature grid).

## Experimental chemistry the P4 calibration rests on

`tests/data/p4_calibration` grades ALIE against reactive sites known from
experiment, so those sites are cited like any method. Each was read in full
except Stock & Brown, and the manifest records how far each was read.

| case | what the source establishes | reference |
| --- | --- | --- |
| thiocyanate | kinetic S-alkylation, kS/kN 10²–10³ (SN2) and 10³–10⁴ (activation-controlled cations) | Loos, Kobayashi & Mayr 2003 `loos2003` |
| nitrite | the O/N ratio of SN2 alkylation depends on solvent and reagent | Tishkov et al. 2005 `tishkov2005`; Kornblum et al. 1955 `kornblum1955` |
| enolate, and ambident anions in general | O kinetically, C thermodynamically, for a free enolate; ambident selectivity is set by intrinsic barriers and product stability rather than hard/soft matching | Mayr, Breugst & Ofial 2011 `mayr2011` |
| phenoxide | O-alkylation in most solvents, C-alkylation in water, phenol and fluorinated alcohols, attributed to hydrogen bonding | Kornblum, Berrigan & le Noble 1963 `kornblum1963` |
| substituted benzenes | ortho/para vs meta directing effects (partial rate factors) | Stock & Brown 1963 `stock1963` |

## Methods and software espmap runs

| what | reference |
| --- | --- |
| PySCF (SCF, PCM, integrals) | `sun2018pyscf`, `sun2020pyscf` |
| GPU4PySCF (`--gpu`) | `li2025gpu4pyscf`, `wu2025gpu4pyscf` |
| B3LYP | `becke1993`, `stephens1994` |
| def2 basis sets | `weigend2005` |
| PCM (C-PCM, as PySCF implements it) | `tomasi2005`, `barone1998` |
| Gaussian-smeared PCM surface charges (used by `reaction_field_at_points`) | `lange2010`, `york1999` (cited in PySCF's `solvent/pcm.py`) |
| RDKit (parsing, gate, embedding) | `rdkit` |
| ETKDG conformer generation | `riniker2015`, `wang2020` |
| MMFF94 / MMFF94s | `halgren1996`, `halgren1999` |
| scikit-image | `vanderwalt2014` |
| pyberny (`--optimize`) | `pyberny` |

## What espmap adapts, and what is specific to it

**Adapted, not new.** Evaluating the ESP at the vertices of a marching-cubes
mesh is an implementation of the polygonal-surface approach above, inside a
small PySCF pipeline. espmap's part is the measurement: the ~100× cost
saving against a volumetric grid, and the grid-spacing convergence row of
its error budget.

**Specific to espmap as far as we know, but unsearched.** None of these has
had a prior-art search, so none is claimed as novel:
- a measured error budget whose every row is recomputed from `bench/` data
  by a test (`tests/test_claims.py`);
- a structure gate run before compute, with a check that the computed
  coordinates are the approved structure (`verify_computed`);
- surface-integrity checks that separate a clipped surface from a folded,
  higher-genus one (`mesh_topology`, `surface_problems`);
- the colormap policy: diverging only when the field changes sign, with the
  orientation of every map asserted by test;
- separate `analysis_key` and `render_key` comparability keys, and a
  `validity` envelope in the JSON record;
- for the per-site work (PR #24, after the methods above): triangle-based
  atom ownership under two partition rules computed together, partition
  sensitivity reported per run, and caller-defined reaction sites. The
  radius-scaled rule itself is a weighted Voronoi assignment
  (`aurenhammer1991`), not new. Whether Multiwfn's atomic surface
  decomposition already offers both rules side by side has not been
  checked (review 5 R15); do that before claiming the pairing.

Before any of these appears in a paper as a contribution, it needs a
targeted literature search, recorded here.

## Development note

AI-assisted coding tools (Claude, via Claude Code) were used for
implementation, benchmarking scripts and the adversarial reviews; commits
made that way carry a `Co-Authored-By` line. The scientific specification,
acceptance criteria and interpretation are the author's responsibility.
Every change is tested against the invariants in `CLAUDE.md` and the data in
`bench/`. AI tools are not sources for chemistry: every established concept
above is cited to the literature.
