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
| Average local ionisation energy | roadmap (v0.3) | Sjoberg et al. 1990 `sjoberg1990`; Bulat et al. 2010 |
| Marching cubes | scikit-image `measure.marching_cubes` | Lorensen & Cline 1987 `lorensen1987`; Lewiner et al. 2003 `lewiner2003` |
| Perceptual colormaps; diverging maps neutral at zero | `build_maps`, the colormap policy | Crameri et al. 2020 `crameri2020`, `crameri_scm`; Moreland 2009 `moreland2009` |

The ~98 % often quoted for the density surface is Bader's figure for the
**0.001** au envelope (H atoms in hydrocarbons, `bader1990`). For the 0.002
surface espmap uses, `bench/isosurface_enclosed.jsonl` measures 98.5–99.1 %
of the electrons inside it for water, chlorobenzene, acetate and decanoic acid
(B3LYP/def2-TZVP, molecular quadrature grid).

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
  sensitivity reported per run, and caller-defined reaction sites.

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
