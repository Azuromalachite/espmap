"""Generate notebooks/espmap_colab.ipynb (ColabFold-style layout).

Edit this file, not the .ipynb:  python notebooks/build_notebook.py notebooks/espmap_colab.ipynb
"""
import json
import sys

BRANCH = "master"
cells = []


def md(text):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")})


def code(text):
    cells.append({"cell_type": "code", "execution_count": None,
                  "metadata": {"cellView": "form"}, "outputs": [],
                  "source": text.strip("\n")})


md(f"""
# espmap · electrostatic potential maps

<a href="https://colab.research.google.com/github/Azuromalachite/espmap/blob/{BRANCH}/notebooks/espmap_colab.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

**Quick start:** type a SMILES in the first cell, then **Runtime → Run all**. You get a rotatable 3D surface, an image, and a zip of everything.

SMILES → structure check → DFT (PySCF) → electrostatic potential on the ρ = 0.002 e·bohr⁻³ surface → a colour scale chosen from the data.

**Reading the colours:** warm/red = **negative** = electron-rich; cool/blue = **positive** = electron-poor (the chemistry convention).

**Runtime:** CPU is fine for small molecules. Decanoic acid took about a minute at `preview` on a 16-thread laptop; Colab's free CPU runtime has 2 vCPUs, so expect several times that. For C12–C20 molecules or the publication setting, use a GPU: *Runtime → Change runtime type → L4*. On an L4 the SCF ran 4.3× faster than a 16-thread laptop and matched the CPU to < 0.001 kcal/mol.

Research code from [github.com/Azuromalachite/espmap](https://github.com/Azuromalachite/espmap) (branch `{BRANCH}`). Limits, error budget and citations are at the bottom.
""")

code("""
#@title 1 · Input molecule, then hit `Runtime` → `Run all`
SMILES = "CCCCCCCCCC(=O)O"  #@param {type:"string"}
#@markdown A name for the results folder (letters, digits, `_`). A short hash of the molecule and settings is appended.
NAME = "decanoic_acid"  #@param {type:"string"}
#@markdown Net charge: `auto` takes it from the SMILES (e.g. `CC(=O)[O-]` is −1).
CHARGE = "auto"  #@param ["auto", "-2", "-1", "0", "1", "2"]
#@markdown Unpaired electrons (2S). 0 for ordinary closed-shell molecules.
SPIN = 0  #@param {type:"integer"}

import re
SMILES = "".join(SMILES.split())          # stray spaces/newlines from copy-paste
NAME = re.sub(r"\\W+", "", NAME) or "molecule"
assert SMILES, "Type a SMILES string in the SMILES field first."
print(f"molecule: {SMILES}   name: {NAME}   charge: {CHARGE}   spin: {SPIN}")
""")

code("""
#@title 2 · Advanced settings (the defaults are fine)
#@markdown **Quality.** `preview` = def2-SVP: right for pictures and teaching, but it read V_max about 2.7 kcal/mol low on the molecules benchmarked, and the shift is molecule-specific (bromobenzene −0.76, the other sign). `publication` = def2-TZVP: quote numbers only from this.
QUALITY = "preview (def2-SVP)"  #@param ["preview (def2-SVP)", "publication (def2-TZVP)"]
SOLVENT = "vacuum"  #@param ["vacuum", "water", "hexane", "chloroform", "toluene", "methanol", "acetonitrile", "dmso"]
FUNCTIONAL = "b3lyp"  #@param ["b3lyp", "pbe0", "m06-2x"]
#@markdown **Colour range.** `auto` fits this molecule only. For panels you want to compare, set the same number (kcal/mol, symmetric ±) for every molecule.
COLOUR_RANGE = "auto"  #@param {type:"string"}
#@markdown **Speed.** Density fitting changed V by ≤ 0.02 kcal/mol in our benchmark and roughly halves the SCF time. GPU `auto` uses one if the runtime has it.
DENSITY_FIT = True  #@param {type:"boolean"}
GPU = "auto"  #@param ["auto", "on", "off"]
#@markdown Conformers tried by the force field (the lowest-energy one is computed) and the random seed.
NCONFS = 20  #@param {type:"integer"}
SEED = "0xC0FFEE"  #@param {type:"string"}

BASIS = "def2-tzvp" if QUALITY.startswith("publication") else "def2-svp"
LIM = None if COLOUR_RANGE.strip().lower() in ("", "auto") else float(COLOUR_RANGE)
print(f"{FUNCTIONAL}/{BASIS}, {SOLVENT}, density fitting {'on' if DENSITY_FIT else 'off'}, "
      f"GPU {GPU}, colour range {'auto' if LIM is None else f'±{LIM:g}'}")
""")

code(f"""
#@title 3 · Install (once per runtime: ~1 min, +3 min for the GPU part)
import os, subprocess, sys, time
BRANCH = "{BRANCH}"  #@param {{type:"string"}}
WORK = "/content"
REPO = os.path.join(WORK, "espmap")

def _sh(cmd):
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"install step failed: {{cmd}}\\n{{p.stdout[-1500:]}}\\n{{p.stderr[-1500:]}}")
    return p.stdout

HAS_GPU = subprocess.run("nvidia-smi", shell=True, capture_output=True).returncode == 0
USE_GPU = {{"on": True, "off": False}}.get(GPU, HAS_GPU)
if USE_GPU and not HAS_GPU:
    raise RuntimeError("GPU = on, but this runtime has no GPU. "
                       "Runtime → Change runtime type → L4 (or set GPU to auto/off).")

def write_install_marker(path, branch, commit):
    import os
    with open(path, "w") as f:
        f.write(branch + " " + commit + os.linesep)


def install_marker_current(path, branch):
    # the marker names the branch it was installed from: a changed BRANCH
    # must reinstall, not keep running the old code (review 3 H10)
    import os
    if not os.path.isfile(path):
        return False
    parts = open(path).read().split()
    return len(parts) == 2 and parts[0] == branch


# the marker is written only after the install SUCCEEDS, so re-running is instant
t0 = time.time()
MARKER = os.path.join(WORK, "ESPMAP_READY")
if not install_marker_current(MARKER, BRANCH):
    if os.path.isfile(MARKER):
        print(f"BRANCH is now {{BRANCH}}, not the installed one: reinstalling")
    if os.path.isdir(REPO):
        import shutil
        shutil.rmtree(REPO)
    _sh(f"git clone -q --depth 1 -b {{BRANCH}} https://github.com/Azuromalachite/espmap.git {{REPO}}")
    # PySCF and RDKit pinned to bench/env.json: the error budget was measured
    # with these, and ETKDG geometries (so geometry hashes) depend on RDKit
    _sh(f"pip install -q -e '{{REPO}}[viz]' 'pyscf==2.14.0' 'rdkit==2026.9.1'")
    out = _sh(f"cd {{REPO}} && {{sys.executable}} esp_map.py --selftest")
    if "selftest: PASS" not in out:
        raise RuntimeError("self-test FAILED; do not trust results from this install:\\n" + out)
    write_install_marker(MARKER, BRANCH, _sh(f"git -C {{REPO}} rev-parse --short HEAD").strip())
    print("espmap installed, self-test PASS")
ESPMAP_COMMIT = open(MARKER).read().split()[1]       # part of every result's cache key
if USE_GPU and not os.path.isfile(os.path.join(WORK, "ESPMAP_GPU_READY")):
    # Colab ships cupy-cuda13x; gpu4pyscf-cuda12x brings cupy-cuda12x. Two CuPys conflict.
    _sh("pip uninstall -y -q cupy-cuda13x")
    _sh("pip install -q gpu4pyscf-cuda12x cutensor-cu12")
    open(os.path.join(WORK, "ESPMAP_GPU_READY"), "w").close()
    print("GPU4PySCF installed")
sys.path.insert(0, REPO)
os.chdir(REPO)
import esp_map as em, pyscf, rdkit
gpu_name = _sh("nvidia-smi --query-gpu=name --format=csv,noheader").strip() if HAS_GPU else "none"
print(f"ready in {{time.time() - t0:.0f}} s · espmap {{BRANCH}} @ {{ESPMAP_COMMIT}} · pyscf {{pyscf.__version__}} · rdkit {{rdkit.__version__}} · GPU: {{gpu_name}}"
      f"{{' (in use)' if USE_GPU else ''}}")
if not HAS_GPU:
    print("note: CPU runtime. Fine for small molecules; for C12+ or 'publication' quality, "
          "switch to a GPU runtime.")
""")

code("""
#@title 4 · Run (structure check → DFT → surface)
import hashlib, json, shutil
from rdkit import Chem

mol0 = Chem.MolFromSmiles(SMILES)
if mol0 is None:
    raise ValueError(f"RDKit cannot read the SMILES {SMILES!r}. Check brackets, ring "
                     "numbers and charges, e.g. acetate is CC(=O)[O-].")
canon = Chem.MolToSmiles(mol0)
heavy = mol0.GetNumHeavyAtoms()
untested = sorted({a.GetSymbol() for a in mol0.GetAtoms()} - em.TESTED_ELEMENTS)
if untested:
    print(f"WARNING: {', '.join(untested)} outside the validated element set; results untested.")
if USE_GPU:
    # GPU = auto puts every GPU runtime on --gpu; say when that leaves the validated path
    for why in em.gpu_untested(SPIN, SOLVENT, DENSITY_FIT):
        print(f"WARNING: GPU with {why} is untested (validated: closed shell, gas phase, "
              f"density fitting). Set GPU = off to stay on the CPU path.")
if heavy > 14 and not USE_GPU:
    print(f"note: {heavy} heavy atoms on a CPU runtime; expect several minutes "
          f"({'tens of minutes' if BASIS == 'def2-tzvp' else 'a few minutes'}).")

# everything that changes the files in the folder: the colour range changes
# esp.png and esp_3d.html, and the code version can change every number
settings = dict(smiles=canon, charge=CHARGE, spin=SPIN, basis=BASIS, xc=FUNCTIONAL,
                solvent=SOLVENT, density_fit=DENSITY_FIT, gpu=USE_GPU,
                nconfs=NCONFS, seed=SEED, lim=LIM, espmap=ESPMAP_COMMIT)
tag = hashlib.sha1(json.dumps(settings, sort_keys=True).encode()).hexdigest()[:5]
JOB = os.path.join(WORK, f"{NAME}_{tag}")          # same molecule + settings → same folder
os.makedirs(JOB, exist_ok=True)
npz = os.path.join(JOB, "surface.npz")

# settings.json is written LAST, after espmap succeeded: an interrupted run
# leaves surface.npz without it and is computed again, never reused
if os.path.isfile(os.path.join(JOB, "settings.json")):
    print(f"reusing {JOB} (same molecule and settings were already computed)")
else:
    cmd = [sys.executable, "esp_map.py", "--smiles", SMILES, "--basis", BASIS,
           "--xc", FUNCTIONAL, "--solvent", SOLVENT, "--spin", str(SPIN),
           "--nconfs", str(NCONFS), "--seed", SEED, "--yes",
           "--out", os.path.join(JOB, "esp.png"), "--html", os.path.join(JOB, "esp_3d.html"),
           "--npz", npz, "--json", os.path.join(JOB, "results.json")]
    cmd += [] if CHARGE == "auto" else ["--charge", CHARGE]
    cmd += ["--density-fit"] if DENSITY_FIT else []
    cmd += ["--gpu"] if USE_GPU else []
    cmd += [] if LIM is None else ["--lim", str(LIM)]
    t0 = time.time()
    with open(os.path.join(JOB, "log.txt"), "w") as log:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in p.stdout:
            log.write(line)
            if not line.startswith("  [diagnostic"):
                print(line, end="")
        p.wait()
    if p.returncode:
        shutil.rmtree(JOB)   # leave no half-finished folder for the cache to trust
        raise RuntimeError("espmap stopped (see the messages above). The structure check "
                           "explains rejections; fix the SMILES rather than forcing it.")
    json.dump(settings, open(os.path.join(JOB, "settings.json"), "w"), indent=2)
    print(f"\\ndone in {time.time() - t0:.0f} s → {JOB}")
""")

code("""
#@title 5 · Display { run: "auto" }
#@markdown Changing these re-draws instantly from the saved surface; nothing is recomputed.
VIEW = "3D interactive"  #@param ["3D interactive", "image", "both"]
#@markdown Colour range for display: `auto`, or ± a number in kcal/mol (use the same number across molecules you compare).
DISPLAY_RANGE = "auto"  #@param {type:"string"}
DISPLAY_CMAP = "auto"  #@param ["auto", "vik", "RdBu", "PuOr", "broc", "cork", "lajolla", "batlow", "viridis"]
ELEVATION = 20  #@param {type:"slider", min:-90, max:90, step:5}
AZIMUTH = 172  #@param {type:"slider", min:0, max:360, step:4}
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display

z = np.load(os.path.join(JOB, "surface.npz"))
V = z["V"]
lim = None if DISPLAY_RANGE.strip().lower() in ("", "auto") else float(DISPLAY_RANGE)
vlim, kind, notes = em.choose_scale(V, lim=lim)
for n in notes:
    print("note:", n)
print(f"V = {V.min():+.1f} … {V.max():+.1f} kcal/mol · {kind} scale {vlim[0]:+.1f} … {vlim[1]:+.1f}"
      + ("" if lim is not None or kind != "diverging" else
         " · auto range: NOT comparable with other molecules"))
cmap_name = em.resolve_cmap(kind, None if DISPLAY_CMAP == "auto" else DISPLAY_CMAP)
cmap = em.build_maps()[cmap_name]
if VIEW in ("3D interactive", "both"):
    fig = em.esp_figure(z["verts"], z["faces"], V, [str(s) for s in z["symbols"]],
                        z["coords"], cmap, vlim, title=f"{NAME} · {FUNCTIONAL}/{BASIS} · {SOLVENT}", kind=kind)
    fig.update_layout(height=560)
    fig.show()
if VIEW in ("image", "both"):
    img = em.render(z["verts"], z["faces"], V, cmap, vlim, size=560,
                    elev=ELEVATION, azim=AZIMUTH)
    plt.figure(figsize=(5.6, 5.6)); plt.imshow(img); plt.axis("off")
    plt.title(f"{NAME}   ({vlim[0]:+.0f} … {vlim[1]:+.0f} kcal/mol, {em.scale_caption(kind)})", fontsize=9)
    plt.show()
""")

code("""
#@title 6 · Package and download
#@markdown Zips the image, interactive HTML, surface data, settings, log and citations.
SAVE_TO_DRIVE = False  #@param {type:"boolean"}
CITATIONS = r'''@article{pyscf2020, author = {Sun, Qiming and others}, title = {Recent developments in the {PySCF} program package}, journal = {J. Chem. Phys.}, volume = {153}, pages = {024109}, year = {2020}}
@article{crameri2020, author = {Crameri, Fabio and Shephard, Grace E. and Heron, Philip J.}, title = {The misuse of colour in science communication}, journal = {Nat. Commun.}, volume = {11}, pages = {5444}, year = {2020}}
% RDKit: https://www.rdkit.org   GPU4PySCF (if the GPU was used): https://github.com/pyscf/gpu4pyscf
% espmap: https://github.com/Azuromalachite/espmap
'''
open(os.path.join(JOB, "citation.bib"), "w").write(CITATIONS)
open(os.path.join(JOB, "README.txt"), "w").write(
    "esp.png        static image (warm/red = negative = electron-rich)\\n"
    "esp_3d.html    interactive 3D view; open in any browser, works offline\\n"
    "surface.npz    vertices and faces (bohr), V (kcal/mol), atoms: re-render without DFT\\n"
    "results.json   settings, energies, V range, versions, error-budget caveat\\n"
    "settings.json  the notebook form values; log.txt the full run log\\n\\n"
    "Measured error budget (chlorobenzene sigma-hole): basis def2-SVP -> def2-TZVPP "
    "+3.55 kcal/mol; solvent gas -> water -3.61. Quote numbers only against these margins.\\n")
zip_path = shutil.make_archive(JOB, "zip", JOB)
print(f"packed {zip_path}")
try:
    from google.colab import files
    files.download(zip_path)
    print("download started; if nothing appears, disable your ad blocker and run this cell again")
except ImportError:
    print("not running in Colab; the zip is at", zip_path)
if SAVE_TO_DRIVE:
    from google.colab import drive
    drive.mount("/content/drive")
    dest = "/content/drive/MyDrive/espmap"
    os.makedirs(dest, exist_ok=True)
    shutil.copy(zip_path, dest)
    print("saved to", dest)
""")

md("""
## Optional · Compare a series on one colour scale

Panels are only comparable when they share **one** colour range. Ions need care: an ion's surface never changes sign, so it cannot share a symmetric scale with a neutral molecule honestly. A series of neutral molecules shares a symmetric diverging scale; a series made only of ions shares one *sequential* scale over all their values; a series that mixes the two is **refused** rather than drawn misleadingly.
""")

code("""
#@title Compare a series (skipped unless RUN_SERIES is ticked)
RUN_SERIES = False  #@param {type:"boolean"}
SERIES = "ethanol:CCO, acetaldehyde:CC=O, acetic acid:CC(=O)O, ethyl acetate:CCOC(C)=O"  #@param {type:"string"}
SHARED_RANGE = 60  #@param {type:"number"}


def series_scale(Vs, shared_range):
    # one colour scale for a whole series, by the CLI's own sign rule:
    # sign-changing members share a symmetric diverging range, one-signed
    # members (ions) share a sequential range over all their values, and a
    # mix of the two gets no single honest scale (review 3 H10)
    one_sign = [em.sign_change_ratio(V) <= em.SIGN_CHANGE_MIN for V in Vs]
    if all(one_sign):
        return "sequential", (float(min(V.min() for V in Vs)), float(max(V.max() for V in Vs)))
    if any(one_sign):
        raise ValueError("ions and sign-changing molecules cannot share one colour scale; "
                         "compare them as two separate series")
    return "diverging", (-shared_range, shared_range)


if RUN_SERIES:
    items = [s.split(":", 1) for s in SERIES.split(",")]
    surfs = []
    for name, smi in items:
        name, smi = name.strip(), "".join(smi.split())
        out = os.path.join(WORK, "series", re.sub(r"\\W+", "_", name))
        os.makedirs(out, exist_ok=True)
        cmd = [sys.executable, "esp_map.py", "--smiles", smi, "--basis", BASIS, "--xc", FUNCTIONAL,
               "--solvent", SOLVENT, "--lim", str(SHARED_RANGE), "--yes",
               "--out", f"{out}/esp.png", "--html", f"{out}/esp_3d.html", "--npz", f"{out}/surface.npz"]
        cmd += ["--density-fit"] if DENSITY_FIT else []
        cmd += ["--gpu"] if USE_GPU else []
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode:
            print(f"{name}: FAILED\\n{p.stdout[-800:]}{p.stderr[-800:]}")
            continue
        zz = np.load(f"{out}/surface.npz")
        surfs.append((name, zz))
        print(f"{name:15s} V = {zz['V'].min():+.1f} … {zz['V'].max():+.1f} kcal/mol")
    one_sign = [n for n, zz in surfs if em.sign_change_ratio(zz["V"]) <= em.SIGN_CHANGE_MIN]
    try:
        kind, vl = series_scale([zz["V"] for _, zz in surfs], SHARED_RANGE)
    except ValueError as e:
        raise RuntimeError(f"{e}. One-signed (ions?): {', '.join(one_sign)}.") from None
    over = [n for n, zz in surfs if np.abs(zz["V"]).max() > SHARED_RANGE]
    if kind == "diverging" and over:
        print(f"note: {', '.join(over)} exceed ±{SHARED_RANGE}; those colours saturate. Raise SHARED_RANGE.")
    cm = em.build_maps()[em.resolve_cmap(kind)]
    fig, axes = plt.subplots(1, len(surfs), figsize=(4 * len(surfs), 4.6), squeeze=False)
    for ax, (name, zz) in zip(axes[0], surfs):
        ax.imshow(em.render(zz["verts"], zz["faces"], zz["V"], cm, vl, size=420))
        ax.set_title(name); ax.axis("off")
    sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(*vl))
    fig.colorbar(sm, ax=axes[0].tolist(), shrink=0.7,
                 label=f"V (kcal/mol), {em.scale_caption(kind)}")
    fig.suptitle(f"{FUNCTIONAL}/{BASIS}, one {kind} scale {vl[0]:+.0f} … {vl[1]:+.0f} kcal/mol", y=0.99)
    plt.show()
    print("interactive views:", os.path.join(WORK, "series", "*", "esp_3d.html"))
else:
    print("skipped (tick RUN_SERIES to compare a series)")
""")

md("""
## Optional · Check the GPU against the CPU reference

Re-runs one row of the CPU benchmark (`bench/timing.py`: decanoic acid, B3LYP/def2-TZVP, density-fitted) through espmap's `--gpu` path and compares. Pass: energy within 1×10⁻⁵ Ha, V_min/V_max within 0.05 kcal/mol, same geometry hash. Recorded result: `bench/gpu_l4.json` (L4: PASS, 4.3×).
""")

code("""
#@title GPU vs CPU check (skipped unless RUN_GPU_CHECK is ticked; needs a GPU runtime)
RUN_GPU_CHECK = False  #@param {type:"boolean"}
if RUN_GPU_CHECK:
    if not HAS_GPU:
        raise RuntimeError("No GPU on this runtime: Runtime → Change runtime type → L4.")
    if not os.path.isfile(os.path.join(WORK, "ESPMAP_GPU_READY")):
        raise RuntimeError("Set GPU to auto or on in cell 2, then re-run cell 3 to install GPU4PySCF.")
    REF = dict(geom_sha="26dcaf3519", nao=492, e_tot=-543.8103149269617,
               V_min=-37.277894106462526, V_max=60.079848397342694, t_scf=97.4)
    spec, _ = em.from_smiles("CCCCCCCCCC(=O)O")
    sha = hashlib.sha1(spec.encode()).hexdigest()[:10]
    if sha != REF["geom_sha"]:
        # a different geometry makes every comparison below meaningless
        raise RuntimeError(f"geometry hash {sha} != reference {REF['geom_sha']}: this "
                           f"RDKit embeds a different decanoic acid, so GPU and CPU "
                           f"would be compared on different molecules. Check the "
                           f"rdkit pin in cell 3.")
    t = time.perf_counter()
    m, mf = em.run_scf(spec, basis="def2-tzvp", density_fit=True, gpu=True)
    t_scf = time.perf_counter() - t
    dm = em.total_dm(mf)
    verts, faces, diag = em.extract_surface(m, dm, 0.26)
    Vg = em.esp_at_points(m, dm, verts) * em.HA2KCAL
    got = dict(geom_sha=sha, nao=m.nao, e_tot=float(mf.e_tot), V_min=float(Vg.min()),
               V_max=float(Vg.max()), t_scf=t_scf)
    tol = dict(e_tot=1e-5, V_min=0.05, V_max=0.05)
    print(f"{'':10s} {'CPU ref':>18s} {'GPU':>18s} {'diff':>12s}")
    for k in REF:
        r, g = REF[k], got[k]
        if isinstance(r, str) or k == "nao":
            print(f"{k:10s} {r!s:>18s} {g!s:>18s} {'same' if r == g else 'DIFFERENT':>12s}")
        else:
            flag = "" if k not in tol else ("  ok" if abs(g - r) <= tol[k] else "  FAIL")
            print(f"{k:10s} {r:18.6f} {g:18.6f} {g - r:+12.6f}{flag}")
    print(f"\\nSCF speed-up over the 16-thread laptop CPU: {REF['t_scf'] / t_scf:.1f}x on {gpu_name}")
else:
    print("skipped (tick RUN_GPU_CHECK on a GPU runtime)")
""")

md(f"""
## Limits, troubleshooting, citations

**What the numbers are worth** (measured on chlorobenzene's σ-hole; see the repo README):

| change | effect on V_max (kcal/mol) |
| --- | --- |
| basis def2-SVP → def2-TZVPP | +3.55 |
| solvent gas → water (PCM) | −3.61 |
| functional B3LYP / PBE0 / M06-2X / HF | 0.69 spread |
| density fitting (this notebook's default) | ≤ 0.02 |

The basis error is molecule-specific (bromobenzene: −0.76). Validating on one molecule licenses nothing about another.

**Limits.**
- **One conformer:** the lowest-energy force-field one. Across extended conformers of a C18 acetate the extrema moved ~1 kcal/mol, but folded conformers moved V_min by 5.2 kcal/mol, and the *shape* of the surface changes either way. For a flexible molecule, sample folded shapes before quoting a number.
- **Protonation is your choice:** a carboxylic acid at pH 7 is mostly carboxylate (`[O-]`), and its ESP is completely different. Draw the state you mean.
- **Out of scope:** transition metals, excited states, periodic systems; ωB97X-D is untested.
- **Fluorine shows no σ-hole**, because it does not form one.
- **GPU validation so far:** closed-shell, gas phase, density-fitted only. Solvent and open-shell runs on the GPU are untested.

**Troubleshooting.**
- *"structure rejected"*: read the check; it names what to fix in the SMILES (charges, E/Z, valence).
- *Install failed*: Runtime → Disconnect and delete runtime, then Run all again.
- *Download did nothing*: disable your ad blocker and re-run cell 6, or tick SAVE_TO_DRIVE.
- *Slow*: switch to a GPU runtime, keep `preview` quality, keep density fitting on.
- Bugs: [github.com/Azuromalachite/espmap/issues](https://github.com/Azuromalachite/espmap/issues)

**Cite** PySCF (Sun et al., *J. Chem. Phys.* 153, 024109, 2020), RDKit, GPU4PySCF if used, and, for the colour maps, Crameri et al. (*Nat. Commun.* 11, 5444, 2020). The zip includes `citation.bib`.

*Notebook layout follows [ColabFold](https://github.com/sokrypton/ColabFold)'s conventions. Branch `{BRANCH}`.*
""")

nb = {"nbformat": 4, "nbformat_minor": 0,
      "metadata": {"colab": {"provenance": [], "gpuType": "L4"},
                   "kernelspec": {"name": "python3", "display_name": "Python 3"},
                   "language_info": {"name": "python"}, "accelerator": "GPU"},
      "cells": []}
for c in cells:
    src = c["source"].split("\n")
    c["source"] = [l + "\n" for l in src[:-1]] + [src[-1]]
    nb["cells"].append(c)
json.dump(nb, open(sys.argv[1], "w", encoding="utf-8"), indent=1, ensure_ascii=False)
open(sys.argv[1], "a", encoding="utf-8").write("\n")
print(f"wrote {len(cells)} cells")
