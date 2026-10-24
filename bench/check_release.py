"""Is this tag's evidence about this tag's code? (review 7 R7-14)

    python bench/check_release.py v0.3.1          # in a checkout of the tag

Reads bench/acceptance_v<version>.json and compares the commit that
accepted it with the tag. Only files that cannot change a number or a test
may differ between the two: prose (*.md), CITATION.cff, and the acceptance
file itself, which is written after the run. Lists every problem and exits
1 if there is any. CI runs it on every pushed v* tag (.github/workflows).
"""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# what may differ between the accepted commit and the tag
MAY_DIFFER = (re.compile(r"[^/]*\.md$|.*/[^/]*\.md$"), re.compile(r"CITATION\.cff$"),
              re.compile(r"bench/acceptance_v[^/]+\.json$"))


def commit_of(code):
    """The commit a git-describe string names: the hash after -g, the tag
    itself for an exact tag, None for unknown or dirty code."""
    if not code or code == "unknown" or code.endswith("-dirty"):
        return None
    m = re.fullmatch(r"(v\d+\.\d+\.\d+)(?:-\d+-g([0-9a-f]{7,40}))?", code)
    if m is None:
        return None
    return m[2] or m[1]


def problems(acc, version, changed):
    """Every reason this acceptance does not license tagging `version`, given
    the paths that differ between the accepted commit and the tag."""
    if not isinstance(acc, dict) or acc.get("schema") != "espmap.acceptance/1":
        return ["not an espmap.acceptance/1 file"]
    out = []
    if acc.get("target_release") != version:
        out.append(f"accepts {acc.get('target_release')!r}, not {version!r}")
    if acc.get("declared_code_version") != version:
        out.append(f"the accepted code declared {acc.get('declared_code_version')!r}, not {version!r}")
    code = acc.get("acceptance_code")
    if commit_of(code) is None:
        out.append(f"accepted code {code!r} is unknown or dirty")
    if acc.get("all_pass") is not True:
        out.append("not every acceptance case passed")
    if acc.get("self_contained") is not True:
        out.append("members were computed by other code than the accepted code")
    s = acc.get("suite")
    if not isinstance(s, dict):
        out.append("no full-suite result")
    else:
        if s.get("code") != code:
            out.append(f"the suite ran on {s.get('code')!r}, not the accepted {code!r}")
        if "--runslow" not in (s.get("command") or ""):
            out.append("the suite ran without the slow tests (--runslow)")
        if s.get("failed") != 0 or s.get("errors") != 0 or not s.get("passed"):
            out.append(f"the suite did not pass: {s.get('passed')} passed, {s.get('failed')} failed, "
                       f"{s.get('errors')} errors")
    moved = [p for p in changed if not any(rx.fullmatch(p) for rx in MAY_DIFFER)]
    if moved:
        out.append("changed after acceptance: " + ", ".join(sorted(moved)))
    return out


def main(argv=None):
    tag = (argv or sys.argv[1:])[0]
    version = tag.lstrip("v")
    path = ROOT / "bench" / f"acceptance_v{version}.json"
    if not path.exists():
        print(f"{tag}: no {path.relative_to(ROOT)}")
        return 1
    acc = json.loads(path.read_text(encoding="utf-8"))
    base = commit_of(acc.get("acceptance_code"))
    changed = []
    if base:
        r = subprocess.run(["git", "-C", str(ROOT), "diff", "--name-only", base, tag],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"{tag}: cannot diff the accepted commit {base}: {r.stderr.strip()}")
            return 1
        changed = r.stdout.split()
    bad = problems(acc, version, changed)
    for p in bad:
        print(f"{tag}: {p}")
    print(f"{tag}: {'NOT releasable' if bad else 'evidence is about this code'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
