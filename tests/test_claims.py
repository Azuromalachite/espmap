"""Measured numbers must say the same thing everywhere they are quoted.

Within one day a conclusion changed (folding widens the C18 V_min spread to
5.2 kcal/mol) and the README, CLAUDE.md and roadmap review were updated while
the runtime warning and module docstring kept the old ~1 kcal/mol claim. These
tests make that drift fail CI.

Where a number comes from bench/, it is COMPUTED from the data here, so the
data is the single source and the prose is checked against it. Numbers with
no data file in the repo (the chlorobenzene budget) are pinned once below.
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


def test_basis_step_on_targets_matches_the_data():
    rows = _timing_rows()
    d = [rows[f"{m}|def2-tzvp|df|16t"]["V_max"] - rows[f"{m}|def2-svp|df|16t"]["V_max"]
         for m in ("decanoic_acid", "BDSF", "Z3Z13-18OAc")]
    lo, hi = f"{min(d):.2f}", f"{max(d):.2f}"
    assert (lo, hi) == ("2.65", "2.82")
    for name in ("README.md", "docs/roadmap-review.md"):
        assert lo in TEXT[name] and hi in TEXT[name], f"{name}: SVP->TZVP V_max range"


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


# ------------------------------------------- numbers with no data file here
# chlorobenzene sigma-hole budget, from review 1; pinned once, quoted in three
BASIS_CHLOROBENZENE, SOLVENT_WATER = "+3.55", "−3.61"


@pytest.mark.parametrize("name", ["README.md", "CLAUDE.md"])
def test_chlorobenzene_budget_is_quoted_identically(name):
    assert BASIS_CHLOROBENZENE in TEXT[name] and SOLVENT_WATER in TEXT[name]


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
