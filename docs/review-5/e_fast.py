"""Review 5, no-SCF probes. Run from the repo root with the repo on sys.path:
    python e_fast.py <repo> <p3_checkout> <fxdir>
"""
import io
import json
import sys
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

REPO, P3, FX = (Path(x) for x in sys.argv[1:4])
sys.path.insert(0, str(REPO))
import numpy as np  # noqa: E402
import esp_map as em  # noqa: E402

FX.mkdir(exist_ok=True)


def hdr(t):
    print(f"\n=== {t}")


A = dict(xc="b3lyp", basis="def2-tzvp", solvent="water", iso=0.002, spacing=0.26,
         density_fit=False, geometry_method="etkdg_mmff94s")


def key(**kw):
    return em.record_key(em.analysis_fields(**{**A, **kw}))[:10]


# ---------------------------------------------------------------- R: xc aliases
hdr("xc spellings -> analysis_key, and what PySCF computes for each")
from pyscf import dft, gto  # noqa: E402
import pyscf  # noqa: E402
print("pyscf", pyscf.__version__, "B3LYP_WITH_VWN5 =",
      getattr(pyscf.__config__, "B3LYP_WITH_VWN5", "(unset)"))
mol = gto.M(atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587", basis="sto-3g", verbose=0)
for xc in ["b3lyp", "B3LYP", "b3lypg", "B3LYP-G", "b3lyp5", "b3lyp/vwn5",
           "HYB_GGA_XC_B3LYP", ".2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWN"]:
    try:
        parsed = dft.libxc.parse_xc(xc)
        hyb = dft.libxc.hybrid_coeff(xc)
        fns = tuple((int(i), round(float(c), 4)) for i, c in parsed[1])
    except Exception as e:
        fns, hyb = f"ERR {e}", None
    print(f"  {xc!r:52} key {key(xc=xc)}  hyb={hyb}  terms={fns}")
print("  energies (sto-3g water), to show which spellings are one functional:")
for xc in ["b3lyp", "b3lypg", "b3lyp5", "HYB_GGA_XC_B3LYP"]:
    mf = dft.RKS(mol)
    mf.xc = xc
    print(f"    {xc:20} E = {mf.kernel():.8f}")

# ---------------------------------------------------- geometry_method channels
hdr("geometry_method: one calculation, two keys / two methods, one key")
g_mol = em.geometry_method(molecule="chlorobenzene")
g_atm = em.geometry_method(atoms=True)
print(f"  --molecule chlorobenzene -> {g_mol}  key {key(geometry_method=g_mol)}")
print(f"  same coordinates via --atoms -> {g_atm}  key {key(geometry_method=g_atm)}")
print(f"  identical atom strings: "
      f"{em.monosubstituted_benzene(*em.MOLECULES['chlorobenzene'])[:40]}...")
# a SMILES record re-fed through its own geometry (the H14 reproduction route)
g_sm = em.geometry_method("3D generated from SMILES by ETKDGv3", "MMFF94s")
print(f"  --smiles CCO -> {g_sm} key {key(geometry_method=g_sm)}; "
      f"its geometry via --atoms -> key {key(geometry_method=g_atm)}")

# DFT nitrite file vs an MMFF nitrite file: both file_frozen
dft_file = P3 / "fixtures" / "geometries" / "nitrite.sdf"
spec_d, gi_d = em.from_file(str(dft_file))
spec_m, gi_m = em.from_smiles("N(=O)[O-]")
mmff_xyz = FX / "nitrite_mmff.xyz"
rows = [r.split() for r in spec_m.split(";")]
mmff_xyz.write_text(f"{len(rows)}\nnitrite MMFF94s\n" + "\n".join(" ".join(r) for r in rows) + "\n")
spec_x, gi_x = em.from_file(str(mmff_xyz), charge=-1)
for lab, gi in [("P3 nitrite.sdf (B3LYP/def2-TZVP/PCM geometry)", gi_d),
                ("nitrite_mmff.xyz (MMFF94s geometry)", gi_x)]:
    gm = em.geometry_method(gi["source"], gi["forcefield"])
    print(f"  {lab:48} -> {gm:14} key {key(geometry_method=gm)}")
from rdkit import Chem  # noqa: E402
sdf = Chem.SDMolSupplier(str(dft_file), removeHs=False)[0]
print(f"  the SDF itself declares: geometry = {sdf.GetProp('geometry')!r} (ignored by espmap)")


def bond_lengths(spec):
    X = np.array([[float(v) for v in r.split()[1:4]] for r in spec.split(";")])
    return np.linalg.norm(X[1] - X[0]), np.linalg.norm(X[2] - X[0])


print("  N-O lengths DFT file %.4f %.4f; MMFF %.4f %.4f" % (*bond_lengths(spec_d), *bond_lengths(spec_x)))

hdr("budget_rows for the P3 nitrite fixture (DFT geometry, PCM water, def2-TZVP)")
ids = [r["id"] for r in em.budget_rows(basis="def2-tzvp", solvent="water", charge=-1, n_rot=0)]
print("  ", ids)
print("   'anion geometry' (MMFF -> DFT, 6.5) attached to a DFT geometry:", "anion geometry" in ids)
ids_opt = [r["id"] for r in em.budget_rows(basis="def2-tzvp", solvent="water", charge=-1, n_rot=0)]
print("   same row after --optimize (budget_rows has no geometry argument):", "anion geometry" in ids_opt)
print("   cation, MMFF geometry (charge +1): anion-geometry row present:",
      "anion geometry" in [r["id"] for r in em.budget_rows(basis="def2-tzvp", solvent=None, charge=1, n_rot=0)])

# ------------------------------------------------------------ SMALL_BASIS gaps
hdr("double-zeta recognition and basis aliases")
for b in ["cc-pVDZ", "aug-cc-pVDZ", "6-31G(d)", "6-31+G*", "6-31G**", "def2-SVPD", "pcseg-1",
          "6-31g*"]:
    try:
        nao = gto.M(atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587", basis=b, verbose=0).nao
    except Exception as e:
        nao = f"ERR {type(e).__name__}"
    print(f"  {b:12} basis_key {em.basis_key(b):10} in SMALL_BASIS {em.basis_key(b) in em.SMALL_BASIS!s:5} "
          f"nao(water) {nao}  key {key(basis=b)}")

# ------------------------------------------------ DIMENSION_AMBIGUOUS and C6F6
hdr("planar 3D-header MOL files with z = 0: what the reader decides")


def molblock(atoms, bonds, charges=()):
    lines = ["", "     RDKit          3D", "",
             f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000"]
    for s, (x, y) in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{0.0:10.4f} {s:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for i, j, o in bonds:
        lines.append(f"{i:3d}{j:3d}{o:3d}  0")
    for i, q in charges:
        lines.append(f"M  CHG  1 {i:3d} {q:3d}")
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def ring(n, r, start=0):
    return [(r * np.cos(np.radians(start + 360 * k / n)), r * np.sin(np.radians(start + 360 * k / n)))
            for k in range(n)]


cases = {
    "CO2": molblock([("O", (-1.16, 0)), ("C", (0, 0)), ("O", (1.16, 0))], [(1, 2, 2), (2, 3, 2)]),
    "azide": molblock([("N", (-1.18, 0)), ("N", (0, 0)), ("N", (1.18, 0))],
                      [(1, 2, 2), (2, 3, 2)], [(1, -1), (2, 1), (3, -1)]),
    "BF3": molblock([("B", (0, 0))] + [("F", p) for p in ring(3, 1.31)],
                    [(1, 2, 1), (1, 3, 1), (1, 4, 1)]),
    "nitrate": molblock([("N", (0, 0))] + [("O", p) for p in ring(3, 1.25)],
                        [(1, 2, 2), (1, 3, 1), (1, 4, 1)], [(1, 1), (3, -1), (4, -1)]),
    "ozone": molblock([("O", (0, 0)), ("O", (1.272 * np.sin(np.radians(58.9)), 1.272 * np.cos(np.radians(58.9)))),
                       ("O", (-1.272 * np.sin(np.radians(58.9)), 1.272 * np.cos(np.radians(58.9))))],
                      [(1, 2, 2), (1, 3, 1)], [(1, 1), (3, -1)]),
}
cc, cf = 1.394, 1.335       # C6F6, B3LYP-class bond lengths
c6 = ring(6, cc)
f6 = ring(6, cc + cf)
cases["C6F6 (real geometry)"] = molblock(
    [("C", p) for p in c6] + [("F", p) for p in f6],
    [(i + 1, (i + 1) % 6 + 1, 2 if i % 2 == 0 else 1) for i in range(6)] + [(i + 1, i + 7, 1) for i in range(6)])
cc2, cf2 = 1.327, 1.319      # C2F4
cases["C2F4 (real geometry)"] = molblock(
    [("C", (-cc2 / 2, 0)), ("C", (cc2 / 2, 0))]
    + [("F", (-cc2 / 2 - cf2 * np.cos(np.radians(57)), s * cf2 * np.sin(np.radians(57)))) for s in (1, -1)]
    + [("F", (cc2 / 2 + cf2 * np.cos(np.radians(57)), s * cf2 * np.sin(np.radians(57)))) for s in (1, -1)],
    [(1, 2, 2), (1, 3, 1), (1, 4, 1), (2, 5, 1), (2, 6, 1)])
for name, block in cases.items():
    p = FX / (name.split()[0] + ".mol")
    p.write_text(block)
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        m, probs, _ = em.read_structure(str(p))
    if m is None:
        print(f"  {name:22} BLOCK {probs[0][1]}")
    else:
        print(f"  {name:22} read as '{m.GetProp('_coords')}'"
              + ("  -> coordinates DISCARDED, ETKDG conformer computed instead"
                 if m.GetProp('_coords') == 'drawing' else ""))
print("  remedy text names 'give any atom a non-zero z (or save as .xyz)';"
      " --atoms never reaches coordinate_kind:")


class Reached(Exception):
    pass


def run_main(argv, stop=True):
    em_run = em.estimate_cost
    if stop:
        def s(*a, **k):
            raise Reached
        em.estimate_cost = s
    sys.argv = ["esp_map.py", *argv]
    out, err = io.StringIO(), io.StringIO()
    res = "returned"
    try:
        with redirect_stdout(out), redirect_stderr(err):
            em.main()
    except Reached:
        res = "REACHED COMPUTE"
    except SystemExit as e:
        res = f"SystemExit({str(e.code)[:90]!r})"
    finally:
        em.estimate_cost = em_run
    return res, out.getvalue(), err.getvalue()


co2 = "O -1.16 0 0; C 0 0 0; O 1.16 0 0"
print("  --atoms CO2 (z = 0):", run_main(["--atoms", co2])[0])

# ------------------------------------------------ the gate never sees --charge
hdr("--charge / --spin overrides after the gate approved a different species")
for argv in (["--smiles", "c1ccccc1", "--charge", "2"],
             ["--smiles", "CC(=O)[O-]", "--charge", "0", "--spin", "1"],
             ["--smiles", "c1ccccc1", "--spin", "2"],
             ["--smiles", "c1ccccc1", "--charge", "2", "--check"]):
    res, out, err = run_main(argv)
    verdict = [ln for ln in out.splitlines() if ln.startswith("verdict")]
    flagged = [ln.strip() for ln in out.splitlines() if ln.strip().startswith("[")]
    print(f"  {' '.join(argv):48} gate {verdict} findings {flagged} -> {res}")
    if "WARNING" in err:
        print("      stderr:", [ln for ln in err.splitlines() if "WARNING" in ln][:2])

# ------------------------------------------------ H3+ and the no-heavy-atom path
hdr("H3+ via --atoms")
h3 = "H 0 0 0; H 0.87 0 0; H 0.435 0.7535 0"
print("  --charge 1:", run_main(["--atoms", h3, "--charge", "1"])[0])

# ------------------------------------------------ nu knife-edge
hdr("gipf nu: no positive triangle vs one tiny positive triangle")
from skimage import measure  # noqa: E402
g = np.linspace(-2, 2, 41)
X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
v, f, _, _ = measure.marching_cubes(X ** 2 + Y ** 2 + Z ** 2, level=1.0, spacing=(0.1,) * 3)
f = f.astype(int)
V = -50.0 - 10 * v[:, 2]
gs0 = em.gipf_stats(v, f, V)
i = int(np.argmax(v[:, 0]))       # one vertex made barely positive
V2 = V.copy()
tri = f[np.any(f == i, axis=1)][0]
V2[tri] = 0.01
gs1 = em.gipf_stats(v, f, V2)
a = em.triangle_areas(v, f)
print(f"  all negative: nu={gs0['nu']}  sigma2_pos={gs0['sigma2_pos_kcal2_mol2']}")
print(f"  one triangle at +0.01 ({a[np.all(np.isin(f, tri), axis=1)].sum() / a.sum():.2e} of the area): "
      f"nu={gs1['nu']}  sigma2_pos={gs1['sigma2_pos_kcal2_mol2']}  V_pos_mean={gs1['V_pos_mean_kcal_mol']}")

# --------------------------------------------- site overlap: is it flagged?
hdr("overlapping caller sites (P3 C18 'carbonyl O' and 'acetyl' share map 32)")
sites = json.loads((P3 / "fixtures" / "geometries" / "Z3Z13-18OAc.sites.json").read_text())
seen = {}
for gr in sites["groups"]:
    for t in gr["atoms"]:
        seen.setdefault(str(t).split("+")[0], []).append(gr["label"])
print("  atoms in more than one group:", {k: v for k, v in seen.items() if len(v) > 1})
rep_keys = set(em.site_report.__code__.co_varnames)
print("  does any record field say the sites overlap / tile the surface? ",
      any(k in json.dumps(json.loads((REPO / "docs" / "record-schema.json").read_text()))
          for k in ("overlap", "tile", "covers")))
