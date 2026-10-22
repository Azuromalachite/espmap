"""Conformer experiment for roadmap-review F3.

The error-budget row "conformer ~1 kcal/mol" was measured on glycerol
(3 carbons). This asks whether it transfers to a C18 chain, and whether
the roadmap's proposed ensemble output (± on V_max) reports the quantity
that actually varies. Three questions:

  1. Do ESP extrema stay within ~1 kcal/mol across conformers?
  2. Do shape descriptors (dipole, +/- surface area, sigma of V) vary more?
  3. Is the lowest-MMFF conformer -- what espmap computes today -- compact
     (folded) relative to the rest of the ensemble?

Conformers come from esp_map.from_smiles, i.e. the ONE embedder
(invariant 1); this script only selects among them. DFT is def2-SVP with
density fitting: this is a SPREAD measurement, so the basis offset cancels
to first order, but the absolute values are not quotable.

    python bench/conformers.py            # Z3Z13-18OAc, 100 embedded, 10 DFT
    python bench/conformers.py --dry      # RDKit only: energies and Rg
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from timing import MOLECULES  # noqa: E402

BOHR2_TO_A2 = 0.529177210903 ** 2


def select(rg, energy, k):
    """Lowest-energy conformer, plus k-1 spread evenly across the Rg range.

    Spanning Rg rather than taking the k lowest energies is deliberate: the
    question is how much SHAPE moves the ESP, and the k lowest MMFF minima
    of a long chain are often near-duplicates.
    """
    best = int(np.argmin(energy))
    order = np.argsort(rg)
    picks = [best]
    for q in np.linspace(0, len(order) - 1, k - 1):
        i = int(order[int(round(q))])
        if i not in picks:
            picks.append(i)
    for i in order:                         # top up if Rg ties collapsed picks
        if len(picks) >= k:
            break
        if int(i) not in picks:
            picks.append(int(i))
    return picks


def surface_stats(verts, faces, V):
    tri = verts[faces]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0],
                                         tri[:, 2] - tri[:, 0]), axis=1)
    vf = V[faces].mean(axis=1)
    w = area / area.sum()
    mean = float((w * vf).sum())
    return dict(area_A2=float(area.sum() * BOHR2_TO_A2),
                area_pos_A2=float(area[vf > 0].sum() * BOHR2_TO_A2),
                area_neg_A2=float(area[vf < 0].sum() * BOHR2_TO_A2),
                V_mean=mean,
                V_sigma=float(np.sqrt((w * (vf - mean) ** 2).sum())))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecule", default="Z3Z13-18OAc", choices=MOLECULES)
    ap.add_argument("--nconfs", type=int, default=100)
    ap.add_argument("--k", type=int, default=10, help="conformers sent to DFT")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--dry", action="store_true", help="RDKit stage only")
    ap.add_argument("--out", default=str(HERE / "conformers.jsonl"))
    a = ap.parse_args()

    import esp_map as em
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors3D

    smi = MOLECULES[a.molecule]
    probe, _ = em.parse_for_check(smiles=smi, path=None, charge=0)
    verdict, findings = em.check_structure(Chem.RemoveHs(probe), a.basis, spin=0)
    if verdict == em.BLOCK:                 # gate before compute (invariant 2)
        raise SystemExit(f"gate BLOCK: {findings}")

    t = time.perf_counter()
    _, info = em.from_smiles(smi, n_confs=a.nconfs)
    mol = info["mol"]
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    cids = [c.GetId() for c in mol.GetConformers()]
    energy = np.array([AllChem.MMFFGetMoleculeForceField(
        mol, props, confId=c).CalcEnergy() for c in cids])
    rg = np.array([Descriptors3D.RadiusOfGyration(mol, confId=c) for c in cids])
    heavy = [x.GetIdx() for x in mol.GetAtoms() if x.GetAtomicNum() > 1]
    extent = np.array([
        max(np.linalg.norm(p - q) for p in pos for q in pos)
        for pos in (mol.GetConformer(c).GetPositions()[heavy] for c in cids)])
    best = int(np.argmin(energy))
    rg_rank = int((rg < rg[best]).sum())
    ens = dict(molecule=a.molecule, smiles=smi, n_rot=info["n_rot"],
               n_embedded=len(cids), seed=info["seed"],
               t_embed=time.perf_counter() - t,
               mmff_dE_range=float(energy.max() - energy.min()),
               rg_min=float(rg.min()), rg_max=float(rg.max()),
               extent_min=float(extent.min()), extent_max=float(extent.max()),
               lowest_E_rg=float(rg[best]), lowest_E_extent=float(extent[best]),
               lowest_E_rg_rank=f"{rg_rank} of {len(cids)} are more compact")
    print("ENSEMBLE " + json.dumps(ens), flush=True)
    if a.dry:
        return

    out = Path(a.out)
    done = set()
    if out.exists():
        done = {json.loads(l).get("conf_id") for l in out.read_text().splitlines()
                if l.strip() and json.loads(l).get("molecule") == a.molecule}
    with out.open("a") as f:
        if not done:
            f.write(json.dumps({"ensemble": ens}) + "\n")
        for i in select(rg, energy, a.k):
            if cids[i] in done:
                continue
            spec = em._mol_to_spec(mol, cids[i])
            row = dict(molecule=a.molecule, conf_id=cids[i], basis=a.basis,
                       density_fit=True, is_lowest_E=(i == best),
                       mmff_dE=float(energy[i] - energy[best]),
                       rg=float(rg[i]), extent=float(extent[i]))
            t = time.perf_counter()
            try:
                m, mf = em.run_scf(spec, basis=a.basis, density_fit=True)
                dm = em.total_dm(mf)
                d = mf.dip_moment(unit="Debye", verbose=0)
                verts, faces, diag = em.extract_surface(m, dm, 0.26)
                V = em.esp_at_points(m, dm, verts) * em.HA2KCAL
                row.update(dipole_D=float(np.linalg.norm(d)), chi=diag["chi"],
                           V_min=float(V.min()), V_max=float(V.max()),
                           **surface_stats(verts, faces, V))
            except Exception as e:          # record, don't hide
                row["error"] = repr(e)
            row["t_wall"] = time.perf_counter() - t
            f.write(json.dumps(row) + "\n")
            f.flush()
            print(f"conf {cids[i]:3d}  dE={row['mmff_dE']:5.2f}  rg={row['rg']:.2f}  "
                  + (f"ERROR {row['error'][:120]}" if "error" in row else
                     f"V={row['V_min']:+.2f}..{row['V_max']:+.2f}  "
                     f"mu={row['dipole_D']:.2f} D  {row['t_wall']:.0f}s"),
                  flush=True)


if __name__ == "__main__":
    main()
