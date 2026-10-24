"""Review 7 R7-14: a release is tagged only on evidence about exactly that
code. bench/check_release.py reads the release's acceptance file and the
paths that differ between the accepted commit and the tag, and lists every
problem. Counterexamples first: each case is one way a tag could otherwise
go out on evidence about other code."""
import importlib.util
from pathlib import Path

import pytest

_P = Path(__file__).resolve().parents[1] / "bench" / "check_release.py"


def _mod():
    if not _P.exists():
        pytest.fail("bench/check_release.py does not exist yet")
    spec = importlib.util.spec_from_file_location("check_release", _P)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _acc(**kw):
    a = dict(schema="espmap.acceptance/1", target_release="0.3.1",
             acceptance_code="v0.3.0-12-g0123456789ab", declared_code_version="0.3.1",
             member_compute_codes=["v0.3.0-12-g0123456789ab"], self_contained=True, all_pass=True,
             suite=dict(code="v0.3.0-12-g0123456789ab", passed=1500, failed=0, errors=0,
                        command="pytest --runslow"),
             cases={"E1": dict(verdict="PASS")})
    a.update(kw)
    return a


def _problems(acc=None, version="0.3.1", changed=()):
    return _mod().problems(acc if acc is not None else _acc(), version, list(changed))


def test_clean_evidence_passes():
    assert _problems(changed=["README.md", "bench/acceptance_v0.3.1.json"]) == []


@pytest.mark.parametrize("why, acc, version, changed", [
    ("code changed after acceptance", None, "0.3.1", ["esp_map.py"]),
    ("a test changed after acceptance", None, "0.3.1", ["tests/test_fast.py"]),
    ("a schema changed after acceptance", None, "0.3.1", ["docs/record-schema.json"]),
    ("a fixture changed after acceptance", None, "0.3.1", ["fixtures/records/benzene.json"]),
    ("another release's acceptance", _acc(target_release="0.3.0"), "0.3.1", []),
    ("the code declared another version", _acc(declared_code_version="0.3.1.dev0"), "0.3.1", []),
    ("a case failed", _acc(all_pass=False), "0.3.1", []),
    ("members from other code", _acc(self_contained=False), "0.3.1", []),
    ("dirty code accepted", _acc(acceptance_code="v0.3.0-12-g0123456789ab-dirty"), "0.3.1", []),
    ("unknown code accepted", _acc(acceptance_code="unknown"), "0.3.1", []),
    ("no suite result", _acc(suite=None), "0.3.1", []),
    ("the suite failed", _acc(suite=dict(code="v0.3.0-12-g0123456789ab", passed=10, failed=1, errors=0,
                                         command="pytest --runslow")), "0.3.1", []),
    ("the suite ran without the slow tests", _acc(suite=dict(code="v0.3.0-12-g0123456789ab", passed=10,
                                                             failed=0, errors=0, command="pytest -q")),
     "0.3.1", []),
    ("the suite ran on other code", _acc(suite=dict(code="v0.3.0-9-gfedcba987654", passed=10, failed=0,
                                                    errors=0, command="pytest --runslow")), "0.3.1", []),
    ("not an acceptance file", {"hello": 1}, "0.3.1", []),
])
def test_a_tag_on_other_evidence_is_refused(why, acc, version, changed):
    assert _problems(acc, version, changed), why


def test_the_accepted_commit_is_read_from_its_describe():
    m = _mod()
    assert m.commit_of("v0.3.0-12-g0123456789ab") == "0123456789ab"
    assert m.commit_of("v0.3.0") == "v0.3.0"
    assert m.commit_of("unknown") is None and m.commit_of("v0.3.0-1-gabc1234-dirty") is None
