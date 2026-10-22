"""Review 5: revert each review-4 fix (as a minimal mutation of the fixed
code) in a scratch copy, and run the tests that claim to guard it.
    python e_revert.py <repo> <scratch_copy>
A fix is GUARDED if at least one selected test fails with the mutation in.
"""
import shutil
import subprocess
import sys
from pathlib import Path

REPO, MUT = Path(sys.argv[1]), Path(sys.argv[2])
ONLY = sys.argv[3:]

F, S = "-m not slow", "--runslow"
MUTATIONS = [
    # id, file, old, new, pytest args
    ("N1 re-verify after --optimize", "esp_map.py",
     'verified, computed, mapped_ok = check_identity(\n                spec, probe, approved, False, a.force, "after --optimize")',
     "pass", ["tests/test_fast.py", "-k", "optimis or optimiz"]),
    ("N1 keep atom maps under --optimize", "esp_map.py",
     "if x.GetAtomMapNum()} if gi is not None else None)\n    if gi is not None and not a.optimize:",
     "if x.GetAtomMapNum()} if gi is not None and not a.optimize else None)\n    if gi is not None and not a.optimize:",
     [S, "tests/test_slow.py", "-k", "optimize_keeps_atom_maps"]),
    ("N2 geometry_method in analysis_key", "esp_map.py",
     "                geometry_method=geometry_method,\n", "",
     ["tests/test_fast.py", "-k", "key or geometry_method"]),
    ("N3 validity.input needs a gate", "esp_map.py",
     'input=gate is not None and not gate["forced"] and verified is not False',
     'input=not (gate and gate["forced"]) and verified is not False',
     ["tests/test_fast.py", "-k", "validity"]),
    ("N3 --atoms gated (gate ignores fragments of --atoms)", "esp_map.py",
     '        m = _perceive_bonds(m, charge, "ATOMS")\n        return (None, m) if isinstance(m, list) else (m, [])',
     '        m = _perceive_bonds(m, charge, "ATOMS")\n        if not isinstance(m, list):\n            from rdkit import Chem as _K\n            m = max(_K.GetMolFrags(m, asMols=True), key=lambda f: f.GetNumAtoms())\n        return (None, m) if isinstance(m, list) else (m, [])',
     ["tests/test_fast.py", "-k", "atoms"]),
    ("N5 record key renamed (partition.compared)", "esp_map.py",
     'partition=dict(primary="nearest_nucleus",\n                               compared=',
     'partition=dict(primary="nearest_nucleus",\n                               radius_scaled=',
     ["tests/test_fast.py", "-k", "schema or partition"]),
    ("N5 site_checks disabled", "esp_map.py",
     "    problems = site_checks(rep, P, faces, V)\n", "    problems = []\n",
     ["tests/test_fast.py", "-k", "site"]),
    ("N6 hydrogen named as a site atom", "esp_map.py",
     '            if symbols[i] == "H":\n                raise', '            if False:\n                raise',
     ["tests/test_fast.py", F]),
    ("N6 by:map on unmapped input", "esp_map.py",
     "        if not atom_maps:\n            raise ValueError('\"by\": \"map\"",
     "        if False:\n            raise ValueError('\"by\": \"map\"",
     ["tests/test_fast.py", F]),
    ("N7 no-heavy-atom early return", "esp_map.py",
     '    if all(s == "H" for s in symbols):\n        return None, ("no heavy atom',
     '    if False:\n        return None, ("no heavy atom',
     ["tests/test_fast.py", "-k", "heavy or hydrogen"]),
    ("N8 tri_atom only for valid records", "esp_map.py",
     '        if valid["input"] and valid["surface"]:\n            for r in PARTITION_RULES:',
     '        if True:\n            for r in PARTITION_RULES:',
     [S, "tests/test_slow.py", "-k", "npz_of_an_invalid"]),
    ("N9 basis_key normalisation", "esp_map.py",
     '    return "".join(c for c in name.lower() if c not in "-_ ")', "    return name.lower()",
     ["tests/test_fast.py", "-k", "basis or spelled or budget"]),
    ("N10 stale figure back in the notebook", "notebooks/build_notebook.py",
     None, "# about 2.7 kcal/mol low on the molecules benchmarked\n",
     ["tests/test_claims.py", "-k", "stale"]),
    ("N11 GIPF variance about the overall mean", "esp_map.py",
     "        mu = float((a[m] * vf[m]).sum() / a[m].sum())",
     "        mu = float((a * vf).sum() / a.sum())",
     ["tests/test_fast.py", "-k", "gipf"]),
    ("N13 CLAUDE.md per-site protocol (2.69 -> 2.96)", "CLAUDE.md",
     "the median is 2.69", "the median is 2.96",
     ["tests/test_claims.py", "-k", "claude_md or per_site"]),
    ("N4 atom-level number (11.9 -> 10.9 in README)", "README.md",
     "up to 11.9", "up to 10.9", ["tests/test_claims.py", "-k", "per_site"]),
    ("N12 symmetry rows (bench docstring 2.378 -> 2.278)", "bench/site_partition.py",
     "span 2.378", "span 2.278", ["tests/test_claims.py"]),
    ("N13b full claims file (2.69 -> 2.96 in CLAUDE.md)", "CLAUDE.md",
     "the median is 2.69", "the median is 2.96", ["tests/test_claims.py"]),
    ("N13c CLAUDE.md conformer median 5.1 -> 5.7", "CLAUDE.md",
     "Conformers move the median site by 5.1", "Conformers move the median site by 5.7",
     ["tests/test_claims.py"]),
    # probes for checks the brief asked about; no review-4 fix behind these
    ("probe: verify_computed always True", "esp_map.py",
     "    verified, computed, why = verify_computed(spec, _rdkit_chem().RemoveHs(probe))",
     "    verified, computed, why = True, approved, ''", ["tests/test_fast.py", F]),
    ("probe: verify_mapping neighbour check off", "esp_map.py",
     "        if got != want:\n            return False, (f\"atom order differs: atom {i}",
     "        if False:\n            return False, (f\"atom order differs: atom {i}",
     ["tests/test_fast.py", F]),
    ("probe: --charge override ignored by gate (already the case)", "esp_map.py",
     None, None, []),
]


def run(args):
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", *args],
                       cwd=MUT, capture_output=True, text=True)
    last = [ln for ln in p.stdout.strip().splitlines() if ln.strip()][-1:] or ["(no output)"]
    return p.returncode, last[0]


if MUT.exists():
    shutil.rmtree(MUT)
shutil.copytree(REPO, MUT, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.png"))
print("baseline:", run(["tests/test_fast.py", "tests/test_claims.py"]))
for mid, fname, old, new, args in MUTATIONS:
    if ONLY and not any(o in mid for o in ONLY):
        continue
    if not args:
        continue
    path = MUT / fname
    text = path.read_text(encoding="utf-8")
    if old is None:
        mutated = new + text
    else:
        n = text.count(old)
        if n != 1:
            print(f"{mid}: mutation text found {n} times, SKIPPED")
            continue
        mutated = text.replace(old, new)
    path.write_text(mutated, encoding="utf-8")
    try:
        code, last = run(args)
    finally:
        path.write_text(text, encoding="utf-8")
    verdict = "GUARDED (a test fails)" if code == 1 else (
        "NOT GUARDED (selected tests pass)" if code == 0 else f"pytest exit {code}")
    print(f"{mid:55} -> {verdict:34} | {' '.join(args)} | {last}")
