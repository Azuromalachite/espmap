"""Repository hygiene the other repos learned the hard way (2026-11-07 sweep).

- Evidence a review cites is tracked by git, not merely present on a disk
  (PINDR 66dd7e0: a *.npz ignore rule silently dropped two cited arrays).
- Text files are LF in git and in every checkout, so `git describe --dirty`
  answers the same under Windows git and WSL git (what-tf 20993b8). espmap's
  code identity feeds compare's dirty-build refusal and bench/acceptance.py.
- A run of hours shows its progress while it runs (the v0.3.1 release run
  printed nothing for hours; memory: long runs show progress).
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# a cited review file that was never kept in the repository, and why: a review's
# outputs were often run artefacts, and the review text names them as they were
NOT_KEPT = {
    "docs/review-3/out/acetate.png": "review 3's render of its acetate command; *.png is ignored outside docs/ "
                                     "and the run's out/ was never committed (the review text says so)",
}


def _tracked():
    out = subprocess.run(["git", "ls-files", "docs"], cwd=ROOT, capture_output=True, text=True, check=True)
    return set(out.stdout.split())


def test_every_review_file_a_doc_cites_is_tracked():
    tracked = _tracked()
    cited = set()
    for md in sorted((ROOT / "docs").glob("*.md")):
        cited |= {f"docs/{m}" for m in re.findall(r"review-\d+/[A-Za-z0-9_./-]+\.[A-Za-z0-9]+",
                                                   md.read_text(encoding="utf-8"))}
    missing = sorted(p for p in cited - tracked if p not in NOT_KEPT)
    assert not missing, f"cited in docs/*.md but not tracked by git: {missing}"
    assert not (set(NOT_KEPT) & tracked), "a file listed as not kept is tracked: remove it from NOT_KEPT"


def test_text_files_are_lf_in_git_and_in_every_checkout():
    attrs = ROOT / ".gitattributes"
    assert attrs.exists(), "no .gitattributes: core.autocrlf decides line endings per machine"
    rules = attrs.read_text(encoding="utf-8").split("\n")
    assert "* text=auto eol=lf" in rules
    eol = subprocess.run(["git", "ls-files", "--eol"], cwd=ROOT, capture_output=True, text=True, check=True)
    assert not [line for line in eol.stdout.splitlines() if line.startswith("i/crlf")], "CRLF in the index"


def test_the_suite_run_streams_its_progress(tmp_path):
    """bench/acceptance.run_suite writes each line of the run to its log as it
    arrives, with a timestamped heartbeat, and still parses pytest's summary."""
    sys.path.insert(0, str(ROOT / "bench"))
    import acceptance
    script = ("import sys, time\n"
              "print('tests/test_a.py::test_one PASSED', flush=True)\n"
              "time.sleep(1.5)\n"
              "print('tests/test_a.py::test_two PASSED', flush=True)\n"
              "print('==== 2 passed, 1 skipped in 1.50s ====', flush=True)\n")
    log = tmp_path / "suite.log"
    r = acceptance.run_suite("abc1234", cmd=[sys.executable, "-c", script], log=log, heartbeat_s=0.5)
    text = log.read_text(encoding="utf-8")
    assert "test_one PASSED" in text and "test_two PASSED" in text
    assert re.search(r"heartbeat \d{4}-\d\d-\d\dT\d\d:\d\d:\d\d", text)
    assert (r["passed"], r["skipped"], r["returncode"]) == (2, 1, 0)
