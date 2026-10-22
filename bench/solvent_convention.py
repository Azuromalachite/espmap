"""Solvent error-budget row under both ESP conventions (review 3 H7).

The "gas -> PCM water -3.61 kcal/mol" row was measured as the potential of
the PCM-polarised solute alone. The PCM surface charges have a potential of
their own (the reaction field); adding it is the other legitimate convention,
and review 3 measured it reversing the sign of the shift at def2-SVP. This
re-measures the row at both bases, for water and chloroform, under both
conventions, from one SCF per (basis, solvent).

Same molecule and settings as the row: idealised chlorobenzene (C-Cl
1.740 A), B3LYP, iso 0.002, spacing 0.26 bohr, margin 5.0 bohr.

    python bench/solvent_convention.py          # appends to solvent_convention.jsonl
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import esp_map as em  # noqa: E402

OUT = HERE / "solvent_convention.jsonl"


def done():
    if not OUT.exists():
        return set()
    return {(r["basis"], r["solvent"]) for r in map(json.loads, OUT.read_text().splitlines())}


def main():
    spec = em.monosubstituted_benzene("Cl", 1.740)
    have = done()
    for basis in ("def2-svp", "def2-tzvp"):
        for solvent in (None, "water", "chloroform"):
            if (basis, solvent) in have:
                continue
            t0 = time.time()
            mol, mf = em.run_scf(spec, basis=basis, solvent=solvent)
            dm = em.total_dm(mf)
            v, _, diag = em.extract_surface(mol, dm, 0.26, iso=em.ISO_BADER, margin=5.0)
            V = em.esp_at_points(mol, dm, v) * em.HA2KCAL
            row = dict(molecule="chlorobenzene (idealised, C-Cl 1.740 A)", xc="b3lyp",
                       basis=basis, solvent=solvent, spacing_bohr=0.26,
                       e_tot_Ha=float(mf.e_tot), n_verts=int(len(v)),
                       solute_only=dict(sigma_hole=em.sigma_hole(v, V, mol, i_sub=6)[0],
                                        V_min=float(V.min()), V_max=float(V.max())))
            if solvent:
                rf = em.reaction_field_at_points(mf, v) * em.HA2KCAL
                W = V + rf
                row["with_reaction_field"] = dict(
                    sigma_hole=em.sigma_hole(v, W, mol, i_sub=6)[0],
                    V_min=float(W.min()), V_max=float(W.max()))
                row["reaction_field_on_surface"] = [float(rf.min()), float(rf.max())]
            row["t_wall"] = round(time.time() - t0, 1)
            with OUT.open("a") as f:
                f.write(json.dumps(row) + "\n")
            print(basis, solvent, row["solute_only"]["sigma_hole"],
                  row.get("with_reaction_field", {}).get("sigma_hole"), flush=True)


if __name__ == "__main__":
    main()
