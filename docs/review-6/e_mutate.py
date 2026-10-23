"""Review 6, PR #73: minimal mutations of the ensemble guards; does a test fail?
usage: e_mutate.py PR_DIR WORK_DIR. Copies the PR tree's esp_map.py + tests per
mutation into WORK_DIR and runs the two ensemble test files there."""
import shutil
import subprocess
import sys
from pathlib import Path

PR, WORK = Path(sys.argv[1]), Path(sys.argv[2])
SRC = (PR / "esp_map.py").read_text()
MUT = {
    "M1 drop the computed != ref_smiles backstop": ("        if computed != ref_smiles:", "        if False:"),
    "M2 supplied: skip verify_mapping": ("        mapped, why_m = verify_mapping(spec, first)", "        mapped, why_m = True, ''"),
    "M3 generated: skip per-member identity": ("            if not (ok and mapped):", "            if False:"),
    "M4 generated: no require_stereo": ("_gate(probe, smiles, basis, spin, charge, not supplied)", "_gate(probe, smiles, basis, spin, charge, False)"),
    "M5 no mirror probe": ("for probe in (heavy, mirror) if achiral else (heavy,):", "for probe in (heavy,):"),
    "M6 symmetry cap 1": ("maxMatches=max_matches)", "maxMatches=1)"),
    "M7 CONTACT_MIN_BONDS 7": ("CONTACT_MIN_BONDS = 5 ", "CONTACT_MIN_BONDS = 7 "),
    "M8 CONTACT_SLACK 0.3": ("CONTACT_SLACK = 0.5 ", "CONTACT_SLACK = 0.3 "),
    "M9 DEDUP_RMSD 0.25": ("DEDUP_RMSD = 0.5 ", "DEDUP_RMSD = 0.25 "),
    "M10 supplied: no MEMBER_NOT_3D": ("        if not coordinates_fix_stereo(m):\n            out.append(_ensemble_block(\n                \"MEMBER_NOT_3D\"",
                                       "        if False:\n            out.append(_ensemble_block(\n                \"MEMBER_NOT_3D\""),
}
for name, (old, new) in MUT.items():
    assert SRC.count(old) == 1, (name, SRC.count(old))
    d = WORK / name.split()[0]
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(PR / "tests", d / "tests")
    (d / "esp_map.py").write_text(SRC.replace(old, new))
    for f in ("pyproject.toml", "conftest.py"):
        if (PR / f).exists():
            shutil.copy(PR / f, d / f)
    for extra in ("docs", "bench", "fixtures", "examples", "README.md", "CLAUDE.md", "CITATION.cff"):
        s = PR / extra
        if s.is_dir():
            shutil.copytree(s, d / extra)
        elif s.exists():
            shutil.copy(s, d / extra)
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-x",
                        "tests/test_ensemble_selection.py", "tests/test_ensemble_set.py"],
                       cwd=d, capture_output=True, text=True)
    tail = [l for l in r.stdout.splitlines() if l.strip()][-1]
    failed = [l.split(" ")[1] for l in r.stdout.splitlines() if l.startswith("FAILED")]
    print(f"{name}: {'CAUGHT' if r.returncode else 'NOT CAUGHT'}  ({tail}) {failed[:2]}", flush=True)
    shutil.rmtree(d)
