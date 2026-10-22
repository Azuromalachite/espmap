"""Input geometries for the per-site-spec P3 fixtures, as atom-mapped SDF.

P3 freezes espmap.sites/1 and produces real records for benzene, nitrite,
thiocyanate, an enolate model and a C18 pheromone. This writes their INPUT
geometries, so the records can be regenerated and so P1's map round trip has
real files to test against.

  small molecules   B3LYP/def2-TZVP optimised in PCM water, the condition the
                    records are computed in (Cartesian BFGS on analytic
                    gradients; espmap's --optimize needs the optional pyberny).
                    MMFF is not used for them: the 2026-10-22 calibration pilot
                    found MMFF nitrite N-O 1.424 A against 1.258 from DFT.
  C18 pheromone     Z3Z13-18OAc, two MMFF94s conformers from the same embedding
                    as bench/conformers.py (seed 0xC0FFEE, 100 conformers): the
                    lowest-energy one (conf 4, extended) and the most compact
                    folded one (conf 72, genus 3). A DFT optimisation at TZVP is
                    days, and the geometry row was measured on neutral molecules.

Every heavy atom carries an atom-map number, 1..n in SMILES order; hydrogens
carry none (spec section 1). RDKit writes the maps to the V2000 mapping column.
The P3 records do not use these files as they stand. fixtures/make_records.py
renumbers the C18 maps along the chain (map n = Cn, not map 1 = acetyl CH3). It
also rotates nitrite out of z = 0, which the gate refuses as
DIMENSION_AMBIGUOUS. The coordinates here stay the source either way.

    python bench/fixture_geometries.py      # writes bench/fixtures/*.sdf + geometries.json
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import esp_map as em  # noqa: E402

FIX = HERE / "fixtures"
SMALL = {"benzene": "c1ccccc1", "nitrite": "N(=O)[O-]", "thiocyanate": "[S-]C#N",
         "enolate": "C=C[O-]"}


def build(sym, X):
    return "; ".join(f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(sym, X))


def optimise(sym, X0, charge, basis="def2-tzvp", solvent="water"):
    from scipy.optimize import minimize

    def f(x):
        mol, mf = em.run_scf(build(sym, x.reshape(-1, 3)), basis=basis, charge=charge,
                             solvent=solvent)
        return mf.e_tot, (mf.nuc_grad_method().kernel() / em.BOHR).ravel()
    r = minimize(f, X0.ravel(), jac=True, method="BFGS", options=dict(gtol=3e-4, maxiter=80))
    return r.x.reshape(-1, 3), r


def mapped_mol(info, X):
    """The from_smiles mol (atom order = spec order) with maps and coordinates X."""
    from rdkit import Chem
    mol = Chem.Mol(info["mol"])
    mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(X):
        conf.SetAtomPosition(i, p.tolist())
    conf.Set3D(True)
    mol.AddConformer(conf, assignId=True)
    n = 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1:
            n += 1
            a.SetAtomMapNum(n)
    return mol


def write_sdf(path, mol, props):
    from rdkit import Chem
    for k, v in props.items():
        mol.SetProp(k, str(v))
    w = Chem.SDWriter(str(path))
    w.write(mol)
    w.close()


def bonds(sym, X):
    heavy = [i for i, s in enumerate(sym) if s != "H"]
    return {f"{sym[a]}{a}-{sym[b]}{b}": round(float(np.linalg.norm(X[a] - X[b])), 4)
            for a in heavy for b in heavy if a < b and np.linalg.norm(X[a] - X[b]) < 1.9}


def main():
    FIX.mkdir(exist_ok=True)
    js = FIX / "geometries.json"
    record = json.loads(js.read_text()) if js.exists() else {}
    for name, smi in SMALL.items():
        if (FIX / f"{name}.sdf").exists() and name in record:   # resumable
            continue
        spec, info = em.from_smiles(smi)
        sym = [t.split()[0] for t in spec.split(";")]
        X0 = np.array([[float(v) for v in t.split()[1:4]] for t in spec.split(";")])
        X, r = optimise(sym, X0, info["charge"])
        from rdkit import Chem
        # the optimised geometry must still be the molecule asked for (invariant 3)
        ok, computed, why = em.verify_computed(build(sym, X), Chem.RemoveHs(info["mol"]))
        props = dict(espmap_fixture=name, smiles=smi, geometry="B3LYP/def2-TZVP, PCM water, "
                     "Cartesian BFGS gtol 3e-4 Ha/A", converged=r.success, iterations=r.nit)
        write_sdf(FIX / f"{name}.sdf", mapped_mol(info, X), props)
        record[name] = dict(props, bonds_mmff=bonds(sym, X0), bonds_dft=bonds(sym, X),
                            structure_verified=bool(ok), computed_smiles=computed, note=why)
        print(f"{name:12s} converged={r.success} nit={r.nit} verified={ok}  {bonds(sym, X)}",
              flush=True)
        js.write_text(json.dumps(record, indent=1) + "\n")    # after each molecule

    smi = "CC(=O)OCC/C=C\\CCCCCCCC/C=C\\CCCC"
    _, info = em.from_smiles(smi, n_confs=100)
    for cid, label in ((4, "lowest_energy_extended"), (72, "folded_genus3")):
        X = info["mol"].GetConformer(cid).GetPositions()
        props = dict(espmap_fixture=f"Z3Z13-18OAc_{label}", smiles=smi,
                     geometry=f"MMFF94s conformer {cid} of 100 (from_smiles seed 0xC0FFEE), "
                              f"as bench/conformers.py")
        write_sdf(FIX / f"Z3Z13-18OAc_conf{cid}.sdf", mapped_mol(info, X), props)
        record[f"Z3Z13-18OAc_conf{cid}"] = props
        print(f"Z3Z13-18OAc conf {cid} ({label}) written", flush=True)
    (FIX / "geometries.json").write_text(json.dumps(record, indent=1) + "\n")


if __name__ == "__main__":
    main()
