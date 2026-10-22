"""Per-site partition measurements for per-site-spec P2 (docs/per-site-spec.md).

The site block will quote per-atom and per-site ESP statistics. Before any of
them is quoted, three things have to be measured, and this script records them:

  grid       How far do site areas, means, percentiles and V_min move when the
             isosurface grid is refined 0.26 -> 0.20 -> 0.16 bohr? The
             existing 0.04 kcal/mol spacing row is for V_max only.
  symmetry   Do symmetry-equivalent sites agree (benzene's six C, nitrite's two
             O) on geometries that are symmetric BY CONSTRUCTION? A partition
             bug gives plausible numbers but cannot give equal numbers to
             equivalent atoms by accident (the role far-field decay plays for
             the kernel sign).
  atoms      The grid cases at 0.26 bohr again, with every atom's own share
             (hydrogens NOT folded into their parents): the record's `atoms`
             level, which review 4 N4 found unbudgeted. Each row's `sites`
             equals the matching grid row's, so both levels are one surface.
             `c18atoms` does the same for the C18 conformers.
  c18        How large is the partition sensitivity (nearest-nucleus vs
             radius-scaled) on a flexible molecule, against the conformer
             spread of the same site? Same conformers as bench/conformers.py
             (seed, selection) plus the two folded ones from f9_recheck.jsonl.

The partition here is the PROTOTYPE of P2 (triangle centroid -> nearest
nucleus or nearest by d / r_Bondi; hydrogens folded into their parent heavy
atom; area-weighted statistics), kept so the recorded rows stay reproducible
from this file. Calculations use esp_map's own run_scf / extract_surface /
esp_at_points.

What validates the partition is the symmetry rows, which do not depend on
how the partition was written: sites equivalent by construction agree within
mesh noise. Benzene's six C site means span 2.378–2.566 kcal/mol at 0.26 bohr,
tightening to 2.487–2.539 at 0.16. Nitrite's O1 and O2 means differ by 0.013
kcal/mol at 0.26 bohr (the grid rows, built symmetric). tests/test_claims.py
recomputes all three from site_partition.jsonl.

P2's esp_map.site_report also reproduces this prototype: recomputing the
"grid|nitrite|0.26" row with esp_map agreed to 8.4e-10 over every site and
both rules. That is not a validation (review 4 N12). Both were written in one
session with the same algorithm, so it shows only that the refactor kept the
prototype's behaviour, assumptions included.

Provenance of the recorded rows: every row says commit "unknown", because
they ran in a worktree created from Windows that WSL's git cannot read. They
were computed at espmap d19b031 with esp_map.py unmodified (2026-10-22, the
c18 part resumed across a WSL restart). New rows use em.code_identity().

    python bench/site_partition.py grid symmetry   # ~15 min
    python bench/site_partition.py c18             # ~1 h (12 SCFs, def2-SVP + RI)
    python bench/site_partition.py atoms           # ~5 min (4 SCFs)
    python bench/site_partition.py c18atoms        # ~1 h, as c18

Appends to bench/site_partition.jsonl; rows already present are skipped.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import esp_map as em  # noqa: E402

OUT = HERE / "site_partition.jsonl"
BOHR = em.BOHR
BONDI = {"H": 1.20, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
         "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
         "Se": 1.90}                          # B from Mantina 2009
RULES = ("nearest_nucleus", "radius_scaled")
SPACINGS = (0.26, 0.20, 0.16)

# B3LYP/def2-TZVP geometries optimised in the gas phase (Cartesian BFGS on
# analytic gradients, gtol 3e-4 Ha/A) during the 2026-10-22 calibration pilot.
DFT_GEOMETRIES = {
    "enolate": ("C=C[O-]", -1,
                "C -0.535848 -0.101164 -0.459954; C 0.513443 0.067412 0.419871; "
                "O 1.526819 -0.653583 0.642148; H -1.329381 0.637814 -0.509415; "
                "H -0.592469 -0.966650 -1.115194; H 0.418227 1.016035 1.024274"),
    "thiocyanate": ("[S-]C#N", -1,
                    "S 1.497604 -0.135408 0.0; C -0.162617 -0.130344 0.0; "
                    "N -1.335037 -0.126766 0.0"),
}


def c2v_nitrite(r=1.2577, angle_deg=116.4):
    """Nitrite symmetric by construction, at the pilot's DFT bond length and
    angle (the optimised geometry is C2v only to the optimiser's tolerance)."""
    h = np.radians(angle_deg) / 2
    return (f"N 0.0 0.0 0.0; O {-r * np.sin(h):.6f} {-r * np.cos(h):.6f} 0.0; "
            f"O {r * np.sin(h):.6f} {-r * np.cos(h):.6f} 0.0")


# ----------------------------------------------------------------- partition
def parse_spec(spec):
    rows = [t.split() for t in spec.split(";")]
    return [r[0] for r in rows], np.array([[float(v) for v in r[1:4]] for r in rows])


def parents(symbols, X):
    heavy = [i for i, s in enumerate(symbols) if s != "H"]
    return np.array([i if s != "H" else
                     heavy[int(np.argmin(np.linalg.norm(X[heavy] - X[i], axis=1)))]
                     for i, s in enumerate(symbols)])


def owners(verts_A, faces, X, symbols, rule):
    """Atom index owning each triangle, by its centroid."""
    cen = verts_A[faces].mean(axis=1)
    d = np.linalg.norm(cen[:, None, :] - X[None, :, :], axis=2)
    if rule == "radius_scaled":
        d = d / np.array([BONDI[s] for s in symbols])[None, :]
    return np.argmin(d, axis=1)


def weighted_quantile(v, w, q):
    o = np.argsort(v)
    c = np.cumsum(w[o]) / w.sum()
    return float(v[o][min(np.searchsorted(c, q), len(v) - 1)])


def stats(verts_A, faces, V, mask):
    """Area-weighted statistics over the triangles in mask (spec section 3)."""
    tri = verts_A[faces[mask]]
    a = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    if a.sum() == 0:
        return dict(area_A2=0.0)
    vf = V[faces[mask]].mean(axis=1)
    vv = V[np.unique(faces[mask])]
    return dict(area_A2=float(a.sum()),
                area_neg_A2=float(a[vf < 0].sum()),
                V_area_mean=float((a * vf).sum() / a.sum()),
                V_p05=weighted_quantile(vf, a, 0.05),
                V_p95=weighted_quantile(vf, a, 0.95),
                V_min=float(vv.min()), V_max=float(vv.max()))


def sites(spec, verts_bohr, faces, V):
    """{rule: {site_label: stats}} with H folded into its parent heavy atom."""
    symbols, X = parse_spec(spec)
    P = verts_bohr * BOHR
    par = parents(symbols, X)
    out = {}
    for rule in RULES:
        owner = par[owners(P, faces, X, symbols, rule)]
        out[rule] = {f"{symbols[h]}{h}": stats(P, faces, V, owner == h)
                     for h in sorted(set(par.tolist()))}
    tri = P[faces]
    total = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum()
    out["conservation"] = {r: float(sum(s.get("area_A2", 0.0) for s in out[r].values()) / total - 1.0)
                           for r in RULES}
    return out


def atom_sites(spec, verts_bohr, faces, V):
    """{rule: {atom_label: stats}} for every atom, hydrogens NOT folded: the
    record's `atoms` level, which review 4 N4 found unbudgeted."""
    symbols, X = parse_spec(spec)
    P = verts_bohr * BOHR
    return {rule: {f"{s}{i}": stats(P, faces, V, owners(P, faces, X, symbols, rule) == i)
                   for i, s in enumerate(symbols)} for rule in RULES}


# --------------------------------------------------------------------- rows
def provenance():
    import pyscf
    import rdkit
    return dict(commit=em.code_identity()["code"], pyscf=pyscf.__version__,
                rdkit=rdkit.__version__)


def done_keys():
    if not OUT.exists():
        return set()
    return {json.loads(l)["key"] for l in OUT.read_text().splitlines() if l.strip()}


def write(row):
    with OUT.open("a") as f:
        f.write(json.dumps(row) + "\n")


def surfaces_at(spec, charge, basis, solvent, spacings, density_fit=False):
    """One SCF, then one surface per spacing: the density does not depend on it."""
    mol, mf = em.run_scf(spec, basis=basis, charge=charge, solvent=solvent,
                         density_fit=density_fit)
    dm = em.total_dm(mf)
    for h in spacings:
        t = time.perf_counter()
        verts, faces, diag = em.extract_surface(mol, dm, h)
        V = em.esp_at_points(mol, dm, verts) * em.HA2KCAL
        yield h, verts, faces, V, diag, time.perf_counter() - t


def run_grid_and_symmetry(parts, prov):
    done = done_keys()
    cases = []
    if "grid" in parts:
        for name, (smi, q, spec) in DFT_GEOMETRIES.items():
            cases.append(("grid", name, smi, q, spec, "def2-tzvp", "water", SPACINGS))
        cases.append(("grid", "nitrite", "N(=O)[O-]", -1, c2v_nitrite(), "def2-tzvp", "water", SPACINGS))
        cases.append(("grid", "chlorobenzene", "Clc1ccccc1", 0,
                      em.monosubstituted_benzene(*em.MOLECULES["chlorobenzene"]),
                      "def2-tzvp", None, SPACINGS))
    if "symmetry" in parts:
        cases.append(("symmetry", "benzene", "c1ccccc1", 0,
                      em.monosubstituted_benzene(*em.MOLECULES["benzene"]),
                      "def2-tzvp", None, SPACINGS))
    if "atoms" in parts:
        # the grid cases at 0.26 bohr again, now also per atom; their `sites`
        # must equal the grid rows', which ties the two levels to one surface
        for name, (smi, q, spec) in DFT_GEOMETRIES.items():
            cases.append(("atoms", name, smi, q, spec, "def2-tzvp", "water", (0.26,)))
        cases.append(("atoms", "nitrite", "N(=O)[O-]", -1, c2v_nitrite(), "def2-tzvp", "water", (0.26,)))
        cases.append(("atoms", "chlorobenzene", "Clc1ccccc1", 0,
                      em.monosubstituted_benzene(*em.MOLECULES["chlorobenzene"]),
                      "def2-tzvp", None, (0.26,)))
    for part, name, smi, q, spec, basis, solvent, spacings in cases:
        todo = [h for h in spacings if f"{part}|{name}|{h}" not in done]
        if not todo:
            continue
        for h, verts, faces, V, diag, dt in surfaces_at(spec, q, basis, solvent, todo):
            row = dict(key=f"{part}|{name}|{h}", part=part, molecule=name, smiles=smi,
                       charge=q, xc="b3lyp", basis=basis, solvent=solvent or "gas",
                       spacing_bohr=h, n_verts=int(diag["n_verts"]), genus=diag["genus"],
                       V_min=float(V.min()), V_max=float(V.max()),
                       sites=sites(spec, verts, faces, V), atoms=spec,
                       t_surface=dt, **prov)
            if part == "atoms":
                row["atom_sites"] = atom_sites(spec, verts, faces, V)
            write(row)
            print(f"{part:8s} {name:13s} h={h:.2f}  n_verts={row['n_verts']:6d}  "
                  f"conservation={max(abs(v) for v in row['sites']['conservation'].values()):.1e}",
                  flush=True)


def run_c18(prov, atoms=False):
    sys.argv = sys.argv[:1]
    from conformers import MOLECULES, select
    from rdkit.Chem import AllChem, Descriptors3D
    smi = MOLECULES["Z3Z13-18OAc"]
    _, info = em.from_smiles(smi, n_confs=100)        # same embedding and seed
    mol = info["mol"]
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    cids = [c.GetId() for c in mol.GetConformers()]
    energy = np.array([AllChem.MMFFGetMoleculeForceField(mol, props, confId=c).CalcEnergy()
                       for c in cids])
    rg = np.array([Descriptors3D.RadiusOfGyration(mol, confId=c) for c in cids])
    picks = [cids[i] for i in select(rg, energy, 10)]
    picks += [c for c in (40, 72) if c not in picks]   # the folded pair (f9_recheck)
    done = done_keys()
    for cid in picks:
        key = f"{'c18atoms' if atoms else 'c18'}|Z3Z13-18OAc|{cid}"
        if key in done:
            continue
        spec = em._mol_to_spec(mol, cid)
        for h, verts, faces, V, diag, dt in surfaces_at(spec, 0, "def2-svp", None, (0.26,),
                                                        density_fit=True):
            i = cids.index(cid)
            row = dict(key=key, part="c18atoms" if atoms else "c18", molecule="Z3Z13-18OAc",
                       smiles=smi, conf_id=cid,
                       charge=0, xc="b3lyp", basis="def2-svp", density_fit=True, solvent="gas",
                       spacing_bohr=h, mmff_dE=float(energy[i] - energy.min()), rg=float(rg[i]),
                       n_verts=int(diag["n_verts"]), genus=diag["genus"],
                       V_min=float(V.min()), V_max=float(V.max()),
                       sites=sites(spec, verts, faces, V), atoms=spec, t_surface=dt, **prov)
            if atoms:
                row["atom_sites"] = atom_sites(spec, verts, faces, V)
            write(row)
            print(f"c18 conf {cid:3d}  genus={row['genus']}  rg={row['rg']:.2f}  "
                  f"V={row['V_min']:+.1f}..{row['V_max']:+.1f}", flush=True)


if __name__ == "__main__":
    parts = set(sys.argv[1:]) or {"grid", "symmetry"}
    prov = provenance()
    run_grid_and_symmetry(parts, prov)
    if "c18" in parts:
        run_c18(prov)
    if "c18atoms" in parts:
        run_c18(prov, atoms=True)
