"""Review 6, disposition of review 5 R1-R18: each fix reverted as a minimal
mutation of master dc9d9ef, in a scratch copy (never the worktree), then the
fast suite (plus named small STO-3G slow tests) run against it.

usage: d_revert.py MASTER_TREE WORK_DIR [ids...]
Prints, per mutation, the tests that fail beyond the unmutated baseline.
"""
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

MASTER, WORK = Path(sys.argv[1]), Path(sys.argv[2])
ONLY = set(sys.argv[3:])
PY = str(Path.home() / ".venvs/espmap/bin/python")
SLOW_K = ("other_record_branches or hydrogen_molecule or npz_of_an or forced_charge "
          "or deliberate_triplet or gas_phase_unbound_anion_record or forced_ungated "
          "or record_matches_the_committed or double_zeta")

# (id, file, old, new, run_slow)
M = [
    ("R1a-anion-row-on-charge", "esp_map.py",
     '"anion geometry": charge < 0 and not geometry_method.startswith("dft:"),',
     '"anion geometry": charge < 0,', False),
    ("R1b-basis-row-always", "esp_map.py",
     'keep = {"basis": basis_key(basis) in SMALL_BASIS,', 'keep = {"basis": True,', False),
    ("R1c-declared-origin-ignored", "esp_map.py",
     '        base = label or "undeclared"\n', '        base = "undeclared"\n', False),
    ("R2a-no-charge-override", "esp_map.py",
     "    if charge is not None and charge != formal:", "    if False:", True),
    ("R2b-no-spin-override", "esp_map.py",
     "    elif spin and spin != rad:", "    elif False:", True),
    ("R2c-main-drops-charge", "esp_map.py",
     "smiles=a.smiles, charge=a.charge,", "smiles=a.smiles, charge=None,", True),
    ("R3a-electronic-always-valid", "esp_map.py",
     'electronic=not electronic["problems"],', "electronic=True,", True),
    ("R3b-sites-ignore-electronic", "esp_map.py",
     '    if not valid.get("electronic", True):', "    if False:", True),
    ("R4-genus-int", "docs/record-schema.json", '"genus": "int?"', '"genus": "int"', True),
    ("R5-labels-by-index", "esp_map.py",
     "        by_map = bool(atom_maps) and all(", "        by_map = False and all(", False),
    ("R6a-xc-key-is-label", "esp_map.py",
     "    name = xc_name(xc)\n    if name:\n        return name\n",
     "    return xc.lower()\n    name = xc_name(xc)\n    if name:\n        return name\n", False),
    ("R6b-atoms-route-label", "esp_map.py",
     '    if atoms:\n        return label or "undeclared"',
     '    if atoms:\n        return "atoms_string:" + (label or "undeclared")', False),
    ("R7-noise-bench-ignored", "tests/test_fixtures.py",
     "    return 1.5 * max(spreads)", "    return 0.1 if molecule == 'nitrite' else 0.5", False),
    ("R8a-reason-omits-mapping", "esp_map.py",
     "    if mapped_ok is False:\n        # the record",
     "    if False:\n        # the record", False),
    ("R8b-no-perceived-smiles", "esp_map.py",
     "return False, _perceived_smiles(g, Chem.GetFormalCharge(approved)), (",
     "return False, None, (", False),
    ("R9a-tri-atom-always", "esp_map.py",
     "    if not has_sites:\n        return {}", "    if False:\n        return {}", True),
    ("R9b-minus-one-fill", "esp_map.py",
     '        if t is not None:\n            out[f"tri_atom_{r}"] = t',
     '        out[f"tri_atom_{r}"] = t if t is not None else -np.ones(len(faces), int)', True),
    ("R10a-verdict-str", "docs/record-schema.json",
     '"verdict": "enum:INFO|WARN|BLOCK"', '"verdict": "str"', True),
    ("R10b-tree-ignores-enum", "tests/record_tree.py",
     '        if kind.startswith("enum:"):', '        if False:', True),
    ("R11-flat-geometry-discarded", "esp_map.py",
     '    if not flat_length_problems(mol, lengths):\n        return "geometry"',
     '    if not flat_length_problems(mol, lengths):\n        return "drawing"', False),
    ("R12a-claude-2.69", "CLAUDE.md",
     "extended ones alone the median is 2.69", "extended ones alone the median is 2.96", False),
    ("R12b-claude-5.1", "CLAUDE.md",
     "move the median site by 5.1 **", "move the median site by 5.7 **", False),
    ("R13a-no-pople-fold", "esp_map.py",
     '    return k.replace("(d,p)", "**").replace("(d)", "*")', "    return k", False),
    ("R13b-ccpvdz-dropped", "esp_map.py",
     '"def2-svp", "def2-svpd", "cc-pvdz",', '"def2-svp", "def2-svpd",', False),
    ("R14-gipf-no-floor", "esp_map.py",
     "GIPF_MIN_AREA_FRACTION = 0.05", "GIPF_MIN_AREA_FRACTION = 0.0", False),
    ("R17-fixture-value", "fixtures/records/nitrite.json",
     None, None, False),                       # special-cased below
    ("R18a-version-skew", "esp_map.py",
     '__version__ = "0.2.1"', '__version__ = "0.2.0"', False),
    ("R18b-version-from-metadata", "esp_map.py",
     "    return dict(espmap_version=__version__, code=code,",
     "    return dict(espmap_version=installed, code=code,", False),
    ("overlap-flag", "esp_map.py",
     "sites_tile_surface=bool((counts == 1).all()))", "sites_tile_surface=True)", False),
]


def run(tree, slow):
    args = [PY, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-rf", "--tb=no"]
    if slow:
        args += ["--runslow", "-k", f"not slow or {SLOW_K}"]
    t = time.time()
    p = subprocess.run(args, cwd=tree, capture_output=True, text=True)
    failed = sorted(set(re.findall(r"^FAILED (\S+)", p.stdout, re.M)))
    tail = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else p.stderr[-300:]
    return failed, tail, time.time() - t


tree = WORK / "tree"
if tree.exists():
    shutil.rmtree(tree)
tree.mkdir(parents=True)
# a clean copy of master's tracked files (git archive), no .git
# (the worktree's .git points at a Windows path, so the archive is made on the
#  host first: git -C espmap-review6 archive -o WORK/master.tar dc9d9ef)
subprocess.run(["tar", "-xf", str(WORK / "master.tar"), "-C", str(tree)], check=True)
p = subprocess.run([PY, "-c", "import esp_map;print(esp_map.__file__)"], cwd=tree,
                   capture_output=True, text=True)
print("module imported from the copy:", p.stdout.strip())
base = {}
for slow in (False, True):
    f, tail, dt = run(tree, slow)
    base[slow] = set(f)
    print(f"BASELINE slow={slow}: {tail} ({dt:.0f} s); failing: {f}")
sys.stdout.flush()

for mid, fname, old, new, slow in M:
    if ONLY and mid not in ONLY:
        continue
    path = tree / fname
    orig = path.read_text(encoding="utf-8")
    if mid == "R17-fixture-value":
        import json
        rec = json.loads(orig)
        rec["results"]["V_min_kcal_mol"] += 1e-4
        text = json.dumps(rec, indent=1)
    else:
        n = orig.count(old)
        if n != 1:
            print(f"{mid}: old string occurs {n} times, SKIPPED")
            continue
        text = orig.replace(old, new)
    path.write_text(text, encoding="utf-8")
    try:
        f, tail, dt = run(tree, slow)
    finally:
        path.write_text(orig, encoding="utf-8")
    new_f = sorted(set(f) - base[slow])
    verdict = "GUARDED" if new_f else "NOT GUARDED"
    print(f"{mid} [{'fast+slow-subset' if slow else 'fast'}] {verdict} ({dt:.0f} s): {tail}")
    for x in new_f:
        print(f"    {x}")
    sys.stdout.flush()
