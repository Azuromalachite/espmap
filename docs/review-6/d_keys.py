"""Review 6, R6/R13 residue: one key per method across spellings, and a
setting that changes what a name means (PySCF B3LYP_WITH_VWN5).
usage: d_keys.py CODE_DIR [vwn5]   (vwn5: set the PySCF config before import)"""
import sys
sys.path.insert(0, sys.argv[1])
VWN5 = len(sys.argv) > 2
if VWN5:
    import pyscf.__config__ as cfg
    cfg.B3LYP_WITH_VWN5 = True
import esp_map as em
from pyscf import gto, dft

print("module:", em.__file__, "B3LYP_WITH_VWN5 set:", VWN5)
XC = ["b3lyp", "B3LYP", "b3lypg", "HYB_GGA_XC_B3LYP", "b3lyp5", "B3LYP5",
      "M06-2X", "m062x", "pbe0", "PBE0", "hf", "HF", "cam-b3lyp", "CAMB3LYP", "wb97x",
      "0.5*HF + 0.5*B88, LYP", "0.5*HF+0.5*B88,LYP", "b3p86", "pbe,pbe", "PBE"]
for xc in XC:
    try:
        print(f"xc {xc!r:28} key {em.xc_key(xc)!r:24} label {em._dft_label(xc, 'def2-svp', None)!r}")
    except Exception as e:
        print(f"xc {xc!r:28} ERROR {type(e).__name__}: {e}")

# what 'b3lyp' actually computes in this process (water, sto-3g)
mol = gto.M(atom="O 0 0 0; H 0.96 0 0; H -0.24 0.93 0", basis="sto-3g", verbose=0)
E = {}
for xc in ("b3lyp", "b3lypg", "b3lyp5"):
    mf = dft.RKS(mol); mf.xc = xc; mf.verbose = 0
    E[xc] = mf.kernel()
print("E(b3lyp) - E(b3lypg) = %.3e Ha;  E(b3lyp) - E(b3lyp5) = %.3e Ha"
      % (E["b3lyp"] - E["b3lypg"], E["b3lyp"] - E["b3lyp5"]))
print("geometry_method(optimize=('b3lyp','def2-svp','vacuum')) =",
      em.geometry_method(optimize=("b3lyp", "def2-svp", "vacuum")))
print("declared 'dft:b3lyp/def2-svp' ->", em.declared_geometry_method("dft:b3lyp/def2-svp"))
print("declared 'B3LYP/def2-TZVP, PCM water' ->",
      em.declared_geometry_method("B3LYP/def2-TZVP, PCM water"))

# bases: key, small?, and whether PySCF resolves the name (nao on CH4)
B = ["def2-svp", "def2svp", "DEF2_SVP", "def2-SVP", "def2-tzvp", "def2-TZVP",
     "6-31G(d)", "6-31g*", "6-31G(d,p)", "6-31g**", "6-31+G(d)", "6-31+g*",
     "6-31+G(d,p)", "6-31++G(d,p)", "6-31G(2df,p)", "6-31G(2d)", "6-311G(d)", "6-311g*",
     "cc-pVDZ", "ccpvdz", "aug-cc-pVDZ", "jun-cc-pvdz", "cc-pvdz-f12", "pcseg-0",
     "pcseg-1", "aug-pcseg-1", "def2-svpd", "def2-SV(P)", "sto-3g", "STO-6G",
     "3-21g", "3-21+g", "6-31++g", "6-31g", "ano-rcc-vdzp", "dzp", "def2-tzvpp"]
ch4 = "C 0 0 0; H 0.63 0.63 0.63; H -0.63 -0.63 0.63; H -0.63 0.63 -0.63; H 0.63 -0.63 -0.63"
for b in B:
    try:
        nao = gto.M(atom=ch4, basis=b, verbose=0).nao
    except Exception as e:
        nao = f"unresolved ({type(e).__name__})"
    print(f"basis {b!r:16} key {em.basis_key(b)!r:14} small {em.basis_key(b) in em.SMALL_BASIS!s:5} nao(CH4) {nao}")
