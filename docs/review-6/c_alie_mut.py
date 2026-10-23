"""Review 6: mutate alie_at_points in a scratch copy; does --selftest or any
ALIE test fail?  usage: c_alie_mut.py MASTER SCRATCH"""
import shutil, subprocess, sys
from pathlib import Path
REPO, MUT = Path(sys.argv[1]), Path(sys.argv[2])
PY = sys.executable
A = "            w = (ao @ Cs[:, occ]) ** 2 * ns[occ]    # n_i |psi_i|^2 at each point\n"
MUTS = [
 ("none (control)", None, None),
 ("W1 no occupation weight", A, A.replace(" * ns[occ]", "")),
 ("W2 weight n^2", A, A.replace("* ns[occ]", "* ns[occ] ** 2")),
 ("W3 UKS weighted 2 like RKS", A, A.replace("* ns[occ]", "* (2.0 if len(C) == 2 else ns[occ])")),
 ("S1 signed e (no abs)", "num += w @ np.abs(es[occ])", "num += w @ es[occ]"),
 ("S2 -e (Koopmans) instead of |e|", "num += w @ np.abs(es[occ])", "num += w @ -es[occ]"),
 ("B1 beta channel dropped", "for Cs, es, ns in zip(C, e, n):", "for Cs, es, ns in zip(C[:1], e[:1], n[:1]):"),
 ("B2 beta uses alpha energies", "for Cs, es, ns in zip(C, e, n):", "for Cs, es, ns in zip(C, (e[0],) * len(C), n):"),
 ("B3 denominator alpha only", "            den += w.sum(axis=1)\n", "            den += w.sum(axis=1) if Cs is C[0] else 0\n"),
 ("K1 last chunk skipped", "for s in range(0, len(pts), chunk):\n        ao = numint.eval_ao(mol, pts[s:s + chunk], deriv=0)\n        num",
  "for s in range(0, len(pts) - 1, chunk):\n        ao = numint.eval_ao(mol, pts[s:s + chunk], deriv=0)\n        num"),
 ("K2 chunk overlap (s+chunk-1 slice)", "out[s:s + chunk] = np.divide", "out[s:s + chunk - 1] = np.divide(num[:-1], den[:-1], out=np.full(len(num) - 1, np.nan), where=den[:-1] > 0) if len(num) > 1 else None\n        _ = np.divide"),
 ("V1 virtuals included", "            occ = ns > 0\n", "            occ = ns >= 0\n"),
]
TESTS = [["-m", "pytest", "-q", "-p", "no:cacheprovider", "tests/test_fast.py", "-k", "alie"],
         ["-m", "pytest", "-q", "-p", "no:cacheprovider", "--runslow", "tests/test_slow.py", "-k", "alie"]]
for label, old, new in MUTS:
    if MUT.exists(): shutil.rmtree(MUT)
    shutil.copytree(REPO, MUT, ignore=shutil.ignore_patterns(".git", "docs", "bench", "notebooks", "__pycache__", "fixtures"))
    f = MUT / "esp_map.py"; t = f.read_text()
    if old is not None:
        assert t.count(old) == 1, (label, t.count(old)); f.write_text(t.replace(old, new))
    st = subprocess.run([PY, "esp_map.py", "--selftest"], cwd=MUT, capture_output=True, text=True)
    alie_lines = [l.strip() for l in st.stdout.splitlines() if "ALIE" in l]
    res = [st.stdout.strip().splitlines()[-1] if st.stdout.strip() else "no output: " + st.stderr[-200:]]
    for args in TESTS:
        r = subprocess.run([PY, *args], cwd=MUT, capture_output=True, text=True)
        res.append(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-200:])
    print(f"== {label}\n   selftest: {res[0]} | fast -k alie: {res[1]} | slow -k alie: {res[2]}")
    for l in alie_lines: print("     " + l)
