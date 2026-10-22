"""How much of the electron density lies inside the rho = 0.002 surface?

The README said "~98 % of the density" for 0.002 e/bohr^3. The ~98 % figure
in the literature is Bader's (Atoms in Molecules, 1990) for the 0.001 au
envelope, quoted there for H atoms in hydrocarbons; no source was found for
0.002. This measures both, for the molecules espmap is used on.

The fraction is integrated on the SCF's own molecular quadrature grid (Becke
partitioning, PySCF grids level 3), not on the uniform marching-cubes grid,
whose integral of rho is off by -5..+7 % with spacing (see extract_surface).
Here the integral of rho over the whole grid matches the electron count to
within 2e-4 (all four molecules).

    python bench/isosurface_enclosed.py      # writes isosurface_enclosed.jsonl
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import esp_map as em  # noqa: E402

OUT = HERE / "isosurface_enclosed.jsonl"
MOLECULES = {"water": "O", "chlorobenzene": "Clc1ccccc1", "acetate": "CC(=O)[O-]",
             "decanoic acid": "CCCCCCCCCC(=O)O"}


def main():
    from pyscf.dft import numint
    rows = []
    for name, smi in MOLECULES.items():
        spec, info = em.from_smiles(smi)
        mol, mf = em.run_scf(spec, basis="def2-tzvp", charge=info["charge"])
        dm = em.total_dm(mf)
        g = mf.grids
        rho = numint.eval_rho(mol, numint.eval_ao(mol, g.coords), dm, xctype="LDA")
        total = float(g.weights @ rho)
        row = dict(molecule=name, smiles=smi, basis="def2-tzvp", xc="b3lyp",
                   nelectron=int(mol.nelectron), integral_rho=total,
                   code=em.code_identity()["code"])
        for iso in (0.001, 0.002):
            row[f"inside_{iso}"] = float(g.weights[rho >= iso] @ rho[rho >= iso]) / total
        rows.append(row)
        print(name, {k: round(v, 4) for k, v in row.items() if k.startswith("inside")},
              flush=True)
    OUT.write_text("".join(json.dumps(r) + "\n" for r in rows))


if __name__ == "__main__":
    main()
