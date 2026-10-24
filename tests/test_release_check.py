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
                        command="pytest --runslow", returncode=0),
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


# --- the user's review of the v0.3.1 branch: a failed suite process --------------
_SUITE_OK = dict(code="v0.3.0-12-g0123456789ab", passed=100, failed=0, errors=0,
                 command="pytest --runslow")


def test_a_suite_whose_process_failed_does_not_license_a_tag():
    """pytest exited 1 (a collection error, a crash, an interrupted run) but
    the parsed counts look clean: the counts are not the verdict."""
    assert _problems(_acc(suite=dict(_SUITE_OK, returncode=1)))


def test_a_suite_without_a_return_code_does_not_license_a_tag():
    assert _problems(_acc(suite=dict(_SUITE_OK)))


def test_acceptance_with_a_failed_suite_exits_nonzero(tmp_path, monkeypatch):
    """bench/acceptance.py --suite must fail the command when the suite
    fails, whatever the acceptance cases did."""
    import sys
    import types
    spec = importlib.util.spec_from_file_location("acceptance", _P.parent / "acceptance.py")
    acc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(acc)
    stub = types.SimpleNamespace(_BY_ID={"E1": {}}, _run=lambda case, d: Path(d).mkdir(parents=True))
    monkeypatch.setitem(sys.modules, "test_ensemble_set", stub)
    monkeypatch.setattr(acc.em, "provenance", lambda: dict(
        code="v0.3.0-12-g0123456789ab", espmap_version="0.3.1", python="3", numpy="2", scipy=None,
        scikit_image=None, pyscf=None, rdkit=None, threads=4, blas=None))
    monkeypatch.setattr(acc, "run_suite", lambda code: dict(_SUITE_OK, code=code, returncode=1,
                                                            summary="interrupted"))
    rc = acc.main(["--release", "0.3.1", "--work", str(tmp_path / "w"), "--only", "E1",
                   "--out", str(tmp_path / "a.json"), "--suite"])
    assert rc != 0
