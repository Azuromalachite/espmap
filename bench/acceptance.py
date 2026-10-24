"""Release acceptance: every slow ensemble case of the hostile set, computed
by the code being accepted, written as one self-contained JSON.

Review 7 R7-11: bench/acceptance_v030_ensembles.json had "head": "" and
espmap_version 0.2.1, and its cases were judged against member records
computed at two earlier commits, which nothing in the file said. This file
records which code accepted, which code computed every member, and the
environment, and says whether they are one revision.

    python bench/acceptance.py --release 0.3.1 --work DIR --suite [--only E3 E4]

--suite runs the full suite (--runslow) on the same checkout first and
records it; bench/check_release.py, which CI runs on every v* tag, requires
it (review 7 R7-14).

Each case runs in DIR/<case> through the hostile set's own runner
(tests/test_ensemble_set.py _run). A case directory made by another
protocol, or before run manifests, is refused by the runner (R7-1); use a
fresh DIR for a fresh acceptance. Refuses unknown or dirty code: an
acceptance of uncommitted code accepts nothing. Writes
bench/acceptance_v<release>.json (one file per release; never overwritten:
pass --out to write elsewhere).
"""
import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "tests"))
import esp_map as em  # noqa: E402

SLOW = ("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E10", "E11", "E19")


def member_codes(case_dir):
    """{code: [conformer]} over the member records in a case directory."""
    out = {}
    for js in sorted((case_dir / "conformers").glob("conf_*.json")):
        if js.name.endswith(".part.json"):
            continue
        code = json.loads(js.read_text(encoding="utf-8")).get("provenance", {}).get("code")
        out.setdefault(code, []).append(js.stem)
    return out


def run_suite(code):
    """The full suite, slow tests included, on this checkout: counts as
    pytest reports them (review 7 R7-14: a tag needs it on exactly its code)."""
    import re
    import subprocess
    cmd = [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "--runslow", "-q"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=HERE.parent, capture_output=True, text=True)
    last = next((line for line in reversed(r.stdout.splitlines()) if " in " in line), "")
    n = {k: int(m[1]) if (m := re.search(rf"(\d+) {k}", last)) else 0
         for k in ("passed", "failed", "errors", "skipped", "xfailed", "xpassed")}
    n["errors"] += int(m[1]) if (m := re.search(r"(\d+) error\b", last)) else 0
    return dict(code=code, command="pytest -p no:cacheprovider --runslow -q", returncode=r.returncode,
                wall_s=round(time.time() - t0, 1), summary=last.strip("= "), **n)


def main(argv=None):
    import test_ensemble_set as ts
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--release", required=True, help="the release this accepts, e.g. 0.3.0")
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--suite", action="store_true",
                    help="also run the full suite (--runslow) first and record it; "
                         "bench/check_release.py requires it")
    a = ap.parse_args(argv)
    prov = em.provenance()
    code = prov["code"]
    if code == "unknown" or code.endswith("-dirty"):
        raise SystemExit(f"code identity {code!r}: commit the code being accepted first")
    out = a.out or HERE / f"acceptance_v{a.release}.json"
    if out.exists() and a.out is None:
        raise SystemExit(f"{out} exists; acceptance files are not overwritten (pass --out)")
    suite = None
    if a.suite:
        print("full suite (--runslow) ...", flush=True)
        suite = run_suite(code)
        print(suite["summary"], flush=True)
    cases = {}
    for cid in a.only or SLOW:
        d = a.work / cid
        t0 = time.time()
        try:
            ts._run(ts._BY_ID[cid], d)
            verdict = "PASS"
        except Exception as e:                                   # noqa: BLE001
            verdict = f"FAIL {type(e).__name__}: {str(e)[:300]}"
        row = dict(verdict=verdict, wall_s=round(time.time() - t0, 1),
                   member_compute_codes=member_codes(d) if (d / "conformers").exists() else {})
        ens = d / "ensemble.json"
        if ens.exists():
            r = json.loads(ens.read_text(encoding="utf-8"))
            row.update(ensemble_key=r["ensemble_key"], analysis_key=r["analysis_key"],
                       k=[r["k_requested"], r["k_obtained"]],
                       excluded={e["conformer"]: e.get("code") for e in r["excluded"]},
                       folded=sum(m["folded"] for m in r["members"]),
                       relabelled=sum(bool(m.get("relabelled")) for m in r["members"]),
                       converged=r["convergence"]["converged"],
                       not_converged=len(r["convergence"]["not_converged"]),
                       not_judgeable=len(r["convergence"].get("not_judgeable", [])),
                       kish_n=r.get("kish_n"), validity=r["validity"])
        cases[cid] = row
        print(cid, verdict, f"{row['wall_s']} s", flush=True)
    computed_by = sorted({c for row in cases.values() for c in row["member_compute_codes"]},
                         key=str)
    rec = dict(
        schema="espmap.acceptance/1", target_release=a.release,
        acceptance_code=code, declared_code_version=prov["espmap_version"],
        member_compute_codes=computed_by,
        self_contained=computed_by in ([code], []),
        environment={k: prov[k] for k in ("python", "numpy", "scipy", "scikit_image",
                                          "pyscf", "rdkit", "threads", "blas")},
        all_pass=all(r["verdict"] == "PASS" for r in cases.values()),
        suite=suite,
        cases=cases)
    out.write_text(json.dumps(rec, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {out}: all_pass {rec['all_pass']}, self_contained {rec['self_contained']}")
    return 0 if rec["all_pass"] and rec["self_contained"] else 1


if __name__ == "__main__":
    sys.exit(main())
