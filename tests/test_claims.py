"""Measured numbers must say the same thing everywhere they are quoted.

Within one day a conclusion changed (folding widens the C18 V_min spread to
5.2 kcal/mol) and the README, CLAUDE.md and roadmap review were updated while
the runtime warning and module docstring kept the old ~1 kcal/mol claim. These
tests make that drift fail CI.

Where a number comes from bench/, it is COMPUTED from the data here, so the
data is the single source and the prose is checked against it. That now
includes every row of the chlorobenzene budget (review 3 H8); the few
figures with no data file are marked as prose where they are quoted.
No SCF, no network: reads text and JSON only.
"""
import json
from pathlib import Path

import pytest

import esp_map as em

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "bench"
TEXT = {name: (ROOT / name).read_text(encoding="utf-8")
        for name in ("README.md", "CLAUDE.md", "docs/roadmap-review.md")}
DOC = em.__doc__
NOTE = em.conformer_note(15, 10.5)


def _jsonl(name):
    return [json.loads(l) for l in (BENCH / name).read_text().splitlines() if l.strip()]


def _timing_rows():
    rows = {r["key"]: r for r in _jsonl("results.jsonl") if "error" not in r}
    for r in _jsonl("results_clean_rerun.jsonl"):      # the clean re-run supersedes
        rows[r["key"]] = r
    return rows


# --------------------------------------------------- numbers from bench/ data
def _conformer_vmin():
    ok = [r["V_min"] for r in _jsonl("conformers.jsonl") if "V_min" in r]
    folded = [r["V_min"] for r in _jsonl("f9_recheck.jsonl")]
    return ok, ok + folded


def test_conformer_spreads_match_the_data():
    extended, everything = _conformer_vmin()
    assert len(extended) == 8 and len(everything) == 10
    assert f"{max(extended) - min(extended):.1f}" == "1.1"
    assert f"{max(everything) - min(everything):.1f}" == "5.2"


@pytest.mark.parametrize("where, text", [
    ("README.md", TEXT["README.md"]),
    ("CLAUDE.md", TEXT["CLAUDE.md"]),
    ("module docstring", DOC),
    ("runtime conformer note", NOTE),
    # review 3 H10: the notebook said extrema "move little with conformation"
    ("notebook", (ROOT / "notebooks" / "build_notebook.py").read_text(encoding="utf-8")),
], ids=lambda v: v if len(str(v)) < 40 else "")
def test_folded_conformer_spread_is_quoted_everywhere(where, text):
    assert "5.2" in text, f"{where} does not quote the folded-conformer V_min spread"


def test_density_fitting_bound_matches_the_data():
    rows = _timing_rows()
    shifts = []
    for k, ex in rows.items():
        if "|exact|16t" in k and k.replace("|exact|", "|df|") in rows:
            df = rows[k.replace("|exact|", "|df|")]
            shifts += [abs(df["V_min"] - ex["V_min"]), abs(df["V_max"] - ex["V_max"])]
    assert len(shifts) == 2 * 7, "expected 7 density-fitting/exact pairs"
    assert max(shifts) <= 0.02
    for name in ("README.md", "CLAUDE.md"):
        assert "≤ 0.02" in TEXT[name], f"{name} lost the density-fitting bound"


def _basis_steps(field):
    """SVP -> TZVP shift for EVERY molecule with both rows: review 3 H8 found
    the quoted range left out the one molecule (the C18 alcohol) that most
    supports the README's own point."""
    rows = _timing_rows()
    mols = sorted({k.split("|")[0] for k in rows
                   if f"{k.split('|')[0]}|def2-svp|df|16t" in rows
                   and f"{k.split('|')[0]}|def2-tzvp|df|16t" in rows})
    return {m: rows[f"{m}|def2-tzvp|df|16t"][field] - rows[f"{m}|def2-svp|df|16t"][field]
            for m in mols}


def test_basis_step_on_targets_matches_the_data():
    vmax, vmin = _basis_steps("V_max"), _basis_steps("V_min")
    assert len(vmax) == 4, sorted(vmax)
    lo, hi = f"{min(vmax.values()):+.2f}", f"{max(vmax.values()):+.2f}"
    assert (lo, hi) == ("+1.53", "+2.82")
    nlo, nhi = f"{min(vmin.values()):+.2f}", f"{max(vmin.values()):+.2f}"
    assert (nlo, nhi) == ("+0.29", "+2.37")
    for name in ("README.md", "docs/roadmap-review.md"):
        for x in (lo, hi, nhi):
            assert x in TEXT[name], f"{name}: SVP->TZVP range lacks {x}"
    assert f"{lo} to {hi}" in TEXT["README.md"]


def test_density_fitting_speed_claim_matches_the_data():
    """'Halves the SCF time' held at def2-TZVP only (review 3 H8)."""
    rows = _timing_rows()
    ratio = {k: rows[k.replace("|exact|", "|df|")]["t_scf"] / ex["t_scf"]
             for k, ex in rows.items()
             if "|exact|16t" in k and k.replace("|exact|", "|df|") in rows}
    svp = [r for k, r in ratio.items() if "|def2-svp|" in k]
    tz = [r for k, r in ratio.items() if "|def2-tzvp|" in k and rows[k]["nao"] <= 600]
    assert max(tz) <= 0.55
    assert f"{round(100 * (1 - max(svp)))}–{round(100 * (1 - min(svp)))} %" == "20–28 %"
    assert "20–28 %" in TEXT["README.md"]


def test_c18_timing_matches_the_data():
    t = _timing_rows()["Z3Z13-18OAc|def2-tzvp|exact|16t"]["t_scf"]
    assert f"{t:.0f}" == "875"
    for name in ("README.md", "CLAUDE.md"):
        assert "875" in TEXT[name]


def test_gpu_agreement_matches_the_record():
    g = json.loads((BENCH / "gpu_l4.json").read_text())
    assert max(abs(g["diff"]["V_min"]), abs(g["diff"]["V_max"])) < 0.001
    for name in ("README.md", "CLAUDE.md"):
        assert "< 0.001" in TEXT[name]


# --------------------------- the chlorobenzene budget, recomputed from bench/
# Review 3 H8: seven of the eleven rows had no data in the repository, and two
# of them (stereo, geometry) did not reproduce when measured. Every row is now
# computed here from its file and must appear in the README's budget table.
def _budget():
    return {r["key"]: r for r in _jsonl("error_budget.jsonl")}


def _solvent(basis, solvent):
    return next(r for r in _jsonl("solvent_convention.jsonl")
                if r["basis"] == basis and r["solvent"] == solvent)


def _budget_table():
    t = TEXT["README.md"]
    return t[t.index("## Error budget"):t.index("Two rows changed when")]


def _rows():
    b = _budget()
    sh = lambda k: b[k]["sigma_hole"]                      # noqa: E731
    funcs = [sh("basis|def2-tzvp")] + [sh(f"functional|{x}") for x in ("pbe0", "m062x", "hf")]
    geoms = [sh("basis|def2-tzvp"), sh("geometry|mmff"), sh("geometry|dft")]
    cones = {c["sigma_hole"] for c in b["basis|def2-tzvp"]["cones"]}
    gas, water = _solvent("def2-tzvp", None), _solvent("def2-tzvp", "water")
    chcl3 = _solvent("def2-tzvp", "chloroform")
    st = {}
    for k, r in b.items():
        if k.startswith("stereo|"):
            st.setdefault(k.split("|")[1], []).append(r)
    mean = lambda rs, f: sum(r[f] for r in rs) / len(rs)   # noqa: E731
    (a, b2) = st.values()
    stereo = [abs(mean(a, f) - mean(b2, f)) for f in ("V_min", "V_max")]
    seed = max(max(r[f] for r in rs) - min(r[f] for r in rs)
               for rs in (a, b2) for f in ("V_min", "V_max"))
    an = {r["key"]: r for r in _jsonl("anion_geometry.jsonl") if "geometry" in r}
    anion = max(abs(an[f"{s}|mmff|water"][f] - an[f"{s}|dft|water"][f])
                for s in ("nitrite", "thiocyanate", "enolate") for f in ("V_min", "V_max"))
    return {
        "basis": f"{sh('basis|def2-tzvpp') - sh('basis|def2-svp'):+.2f}",
        "solvent water": f"{water['solute_only']['sigma_hole'] - gas['solute_only']['sigma_hole']:+.2f}".replace("-", "−"),
        "solvent chloroform": f"{chcl3['solute_only']['sigma_hole'] - gas['solute_only']['sigma_hole']:+.2f}".replace("-", "−"),
        "reaction field": f"{water['with_reaction_field']['sigma_hole'] - gas['solute_only']['sigma_hole']:+.2f}",
        "functional": f"{max(funcs) - min(funcs):.2f}",
        "grid": f"{abs(sh('grid|0.16') - sh('basis|def2-tzvp')):.2f}",
        "region": f"{max(cones) - min(cones):.2f}",
        "geometry": f"{max(geoms) - min(geoms):.2f}",
        "stereo": f"{min(stereo):.1f}–{max(stereo):.1f}",
        "stereo seeds": f"{seed:.2f}",
        "anion geometry": f"{anion:.1f}",
    }


@pytest.mark.parametrize("row, expect", [
    ("basis", "+3.54"), ("solvent water", "−3.60"), ("solvent chloroform", "−2.69"),
    ("reaction field", "+0.05"), ("functional", "0.69"), ("grid", "0.04"),
    ("region", "0.00"), ("geometry", "0.40"), ("stereo", "1.6–1.9"),
    ("stereo seeds", "1.25"), ("anion geometry", "6.5"),
])
def test_every_budget_row_is_recomputed_from_bench(row, expect):
    got = _rows()[row]
    assert got == expect, f"{row}: data gives {got}, the docs say {expect}"
    assert expect in _budget_table(), f"README budget table lacks {row} = {expect}"


@pytest.mark.parametrize("name", ["README.md", "CLAUDE.md"])
def test_chlorobenzene_budget_is_quoted_identically(name):
    rows = _rows()
    for key in ("basis", "solvent water", "geometry", "stereo"):
        assert rows[key] in TEXT[name], f"{name} lacks {key} = {rows[key]}"


def test_retracted_budget_figures_are_gone():
    """Figures that did not reproduce must not survive where a number is QUOTED:
    the two budget tables, the module docstring and the runtime gate text. The
    README paragraph that explains the change may name them."""
    c = TEXT["CLAUDE.md"]
    quoted = {"README budget table": _budget_table(),
              "CLAUDE.md budget table": c[c.index("## Measured error budget"):
                                          c.index("**The solvent row is a convention")],
              "module docstring": DOC,
              "gate text": __import__("inspect").getsource(em.check_structure)}
    for name, text in quoted.items():
        assert "0.4 noise floor" not in text and "0.9 kcal/mol" not in text, name
        assert "+3.55" not in text and "3.61" not in text and "| 0.35 |" not in text, name


def test_small_basis_warning_quotes_the_basis_penalty():
    import inspect
    assert "+3.5 kcal/mol" in inspect.getsource(em.run_scf)


# --------------------------------------------- retracted claims stay retracted
@pytest.mark.parametrize("where, text, phrase", [
    ("runtime conformer note", NOTE, "barely move"),
    ("module docstring", DOC, "nearly conformation-independent"),
    ("roadmap review", TEXT["docs/roadmap-review.md"], "gives the intrinsic preference"),
], ids=lambda v: v if len(str(v)) < 40 else "")
def test_retracted_wording_does_not_return(where, text, phrase):
    assert phrase not in text, f"{where} reintroduced a retracted claim: {phrase!r}"


# ------------------------------------------------- prior art and citation
PRIOR = (ROOT / "docs" / "prior-art.md").read_text(encoding="utf-8")
BIB = (ROOT / "docs" / "references.bib").read_text(encoding="utf-8")


def _bib_keys():
    import re
    return re.findall(r"^@\w+\{([^,\s]+),", BIB, flags=re.M)


def test_every_reference_cited_in_prior_art_is_in_the_bibliography():
    import re
    keys = set(_bib_keys())
    cited = set(re.findall(r"`([a-z_]+\d{4}[a-z0-9]*|rdkit|pyberny|crameri_scm)`", PRIOR))
    assert cited, "prior-art.md cites nothing"
    assert not cited - keys, f"cited but not in references.bib: {sorted(cited - keys)}"
    assert not keys - cited, f"in references.bib but never cited: {sorted(keys - cited)}"


def test_every_bibliography_entry_has_a_doi():
    import re
    entries = re.split(r"^@", BIB, flags=re.M)[1:]
    missing = [e.split(",")[0] for e in entries if "doi" not in e]
    assert not missing, f"no DOI: {missing}"


def test_citation_cff_version_matches_the_package():
    import re
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    v_cff = re.search(r"^version: (\S+)", cff, re.M).group(1)
    v_pkg = re.search(r'^version = "([^"]+)"', pyproject, re.M).group(1)
    assert v_cff == v_pkg


def test_enclosed_density_figure_matches_the_data():
    """The README said "~98 %", which is Bader's figure for the 0.001 au
    envelope (H in hydrocarbons). The 0.002 figure espmap quotes is its own."""
    f = [r["inside_0.002"] for r in _jsonl("isosurface_enclosed.jsonl")]
    rng = f"{100 * min(f):.1f}–{100 * max(f):.1f} %"
    assert rng == "98.5–99.1 %"
    assert rng in TEXT["README.md"] and rng in PRIOR
    assert "~98 %" not in TEXT["README.md"].split("## Credit")[0]
    assert "~98%" not in DOC


@pytest.mark.parametrize("phrase", ["worth borrowing", "trick worth keeping",
                                    "idea worth preserving"])
def test_borrowed_method_is_not_presented_as_an_invention(phrase):
    for name in ("README.md", "CLAUDE.md"):
        assert phrase not in TEXT[name], name
    assert phrase not in (ROOT / "docs" / "colormaps.md").read_text(encoding="utf-8")


# ------------------------------------------------- per-site budget (P2 rows)
def _site_table():
    t = TEXT["README.md"]
    return t[t.index("**Per-site statistics**"):t.index("**The solvent rows measure")]


def _site_rows():
    """The per-site rows, recomputed from bench/site_partition.jsonl with the
    definitions the README states: site means, nearest-nucleus rule, 0.26 bohr."""
    import statistics
    rows = _jsonl("site_partition.jsonl")
    nn = lambda r: r["sites"]["nearest_nucleus"]          # noqa: E731

    def delta(r, q):
        rs = r["sites"]["radius_scaled"]
        return [abs(rs[s][q] - nn(r)[s][q]) for s in nn(r) if q in nn(r)[s] and q in rs[s]]
    small = [r for r in rows if r["part"] == "grid" and r["spacing_bohr"] == 0.26
             and r["molecule"] != "chlorobenzene"]
    c18 = [r for r in rows if r["part"] == "c18"]
    sites = list(nn(c18[0]))
    per_site = [max(abs(r["sites"]["radius_scaled"][s]["V_area_mean"] - nn(r)[s]["V_area_mean"])
                    for r in c18) for s in sites]
    spread = [max(nn(r)[s]["V_area_mean"] for r in c18) - min(nn(r)[s]["V_area_mean"] for r in c18)
              for s in sites]
    grid = {}
    for r in rows:
        if r["part"] == "grid":
            grid.setdefault(r["molecule"], {})[r["spacing_bohr"]] = nn(r)
    g = lambda q: max(abs(v[0.16][s][q] - v[0.26][s][q])  # noqa: E731
                      for v in grid.values() for s in v[0.26] if q in v[0.26][s])
    return {
        "partition small anions": f"{max(max(delta(r, 'V_area_mean')) for r in small):.1f}",
        "partition C18": f"{max(per_site):.2f}",
        "partition C18 median": f"{statistics.median(per_site):.2f}",
        "partition p05": f"{max(max(delta(r, 'V_p05')) for r in small + c18):.1f}",
        "conformer C18": f"{statistics.median(spread):.1f}",
        "grid mean": f"{g('V_area_mean'):.2f}",
        "grid p05": f"{g('V_p05'):.2f}",
        "grid V_min": f"{g('V_min'):.2f}",
    }


@pytest.mark.parametrize("row, expect", [
    ("partition small anions", "3.0"), ("partition C18", "3.95"),
    ("partition C18 median", "0.89"), ("partition p05", "8.9"), ("conformer C18", "5.1"),
    ("grid mean", "0.15"), ("grid p05", "0.72"), ("grid V_min", "0.54"),
])
def test_every_per_site_row_is_recomputed_from_bench(row, expect):
    got = _site_rows()[row]
    assert got == expect, f"{row}: data gives {got}, the docs say {expect}"
    assert expect in _site_table(), f"README per-site table lacks {row} = {expect}"


def test_site_areas_conserve_the_surface_in_every_row():
    for r in _jsonl("site_partition.jsonl"):
        assert max(abs(v) for v in r["sites"]["conservation"].values()) < 1e-12, r["key"]


def test_claude_md_quotes_the_per_site_budget():
    rows = _site_rows()
    for key in ("partition small anions", "partition C18", "conformer C18", "grid mean"):
        assert rows[key] in TEXT["CLAUDE.md"], f"CLAUDE.md lacks {key} = {rows[key]}"


# ------------------------------------------- stale figures stay gone (review 4 N10)
STALE = [
    ("+9.59", "the retracted DFT-geometry sigma-hole (now +9.54, error_budget.jsonl)"),
    ("0.35 kcal/mol spread", "the retracted geometry row (now 0.40)"),
    ("about 2.7 kcal/mol", "the three-molecule basis claim (now +1.53 to +2.82)"),
    ("+2.7 `[M]`", "the three-molecule basis claim (now +1.53 to +2.82)"),
    ("roughly halves", "density fitting halves at def2-TZVP only"),
    ("100 fast tests", "a test count that went stale; the docs give none now"),
    ("13 tests with real DFT", "a test count that went stale; the docs give none now"),
]
SCANNED = {
    "README.md": TEXT["README.md"],
    "CLAUDE.md": TEXT["CLAUDE.md"],
    "docs/roadmap-review.md": TEXT["docs/roadmap-review.md"],
    "docs/colormaps.md": (ROOT / "docs" / "colormaps.md").read_text(encoding="utf-8"),
    "notebooks/build_notebook.py": (ROOT / "notebooks" / "build_notebook.py").read_text(
        encoding="utf-8"),
    "module docstring": DOC,
}


@pytest.mark.parametrize("phrase, why", STALE, ids=[p for p, _ in STALE])
def test_stale_figures_do_not_come_back(phrase, why):
    """Review 4 N10: retracted or superseded figures survived in the notebook,
    colormaps.md and the roadmap review after the tables were fixed."""
    hits = [name for name, text in SCANNED.items() if phrase in text]
    assert not hits, f"{phrase!r} ({why}) is still quoted in {hits}"


def test_prior_art_does_not_call_merged_work_planned():
    assert "planned (v0.2 per-site records" not in PRIOR
