"""The P3 fixture records (docs/per-site-spec.md section 8): real espmap
records, under the frozen contract (docs/record-schema.json,
espmap.record/1), that Arrowfield v0.0 is written against.

  small molecules  benzene, nitrite, thiocyanate, enolate. B3LYP/def2-TZVP in
                   PCM water, the condition their geometries were optimised in
                   (bench/fixtures/*.sdf, bench/fixture_geometries.py). Default
                   sites: each heavy atom with its hydrogens.
  C18 pheromone    Z3Z13-18:OAc, the lowest-energy extended conformer (4) and
                   the most compact folded one (72, genus 3). B3LYP/def2-SVP +
                   RI in the gas phase, the protocol of the C18 budget rows
                   (bench/site_partition.jsonl). Caller sites by atom map.

Atom maps on the C18 fixtures follow the CHAIN, not the SMILES: map n is
carbon Cn (C1 bears the ester oxygen), the acetyl carbons are 21 (CH3) and
22 (C=O), the ester O 31 and the carbonyl O 32. bench/fixtures numbers maps
in SMILES order (map 1 = acetyl CH3), which no chemist reads as C1; the
coordinates here are those files' coordinates, only the maps change.

Nitrite's DFT geometry lies exactly in z = 0 with two equal N-O bonds, which
v0.2.0's structure gate refused as DIMENSION_AMBIGUOUS (H5: a drawing looks
the same; review 5 R11 now keeps it). Its fixture input is that geometry
rigidly rotated out of the plane, the remedy the gate named: every internal
coordinate, and so the ESP, is unchanged; only the frame differs from
bench/fixtures/nitrite.sdf.

Every record is computed from a committed file with every atom supplied, so
nothing is embedded or relaxed: geometry_provenance.relaxation is "none".

    python fixtures/make_records.py                      # the fixtures (~10 min, 4 threads)
    python fixtures/make_records.py --basis sto-3g --out DIR   # dry run, minutes

One process per record. A record whose JSON exists is skipped only if the same
code revision made it; a real run refuses an unknown or dirty code identity.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GEOM = HERE / "geometries"
BENCH_FIX = ROOT / "bench" / "fixtures"

C18_SITES = {"by": "map", "groups": (
    [{"label": f"C{n}-H", "atoms": [f"{n}+H"]} for n in range(1, 19)]
    + [{"label": "head O", "atoms": [31]}, {"label": "carbonyl O", "atoms": [32]},
       {"label": "acetyl", "atoms": ["21+H", 22, 32]}])}

SMALL = dict(basis="def2-tzvp", solvent="water", density_fit=False, sites=None)
C18 = dict(basis="def2-svp", solvent=None, density_fit=True, sites="Z3Z13-18OAc.sites.json")
FIXTURES = {
    "benzene": dict(SMALL, file=BENCH_FIX / "benzene.sdf"),
    "nitrite": dict(SMALL, file=GEOM / "nitrite.sdf", source=BENCH_FIX / "nitrite.sdf",
                    transform="tilt"),
    "thiocyanate": dict(SMALL, file=BENCH_FIX / "thiocyanate.sdf"),
    "enolate": dict(SMALL, file=BENCH_FIX / "enolate.sdf"),
    "Z3Z13-18OAc_extended": dict(C18, file=GEOM / "Z3Z13-18OAc_conf4.sdf",
                                 source=BENCH_FIX / "Z3Z13-18OAc_conf4.sdf", transform="chain"),
    "Z3Z13-18OAc_folded": dict(C18, file=GEOM / "Z3Z13-18OAc_conf72.sdf",
                               source=BENCH_FIX / "Z3Z13-18OAc_conf72.sdf", transform="chain"),
}
TILT = (0.6, 0.0, 0.8)          # rotation axis (unit), angle below: any non-trivial one


def tilted(X, angle_deg=35.0):
    """X rotated rigidly about TILT: Rodrigues' formula."""
    import numpy as np
    k, t = np.array(TILT), np.radians(angle_deg)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
    return X @ R.T


def chain_maps(mol):
    """Atom index -> map for an n-alkenyl acetate: Cn along the chain from the
    ester O, 21/22 for the acetyl carbons, 31/32 for the ester/carbonyl O."""
    def heavy_nbrs(a):
        return [n for n in a.GetNeighbors() if n.GetAtomicNum() > 1]
    o_ester = next(a for a in mol.GetAtoms()
                   if a.GetSymbol() == "O" and len(heavy_nbrs(a)) == 2)
    acyl = next(n for n in heavy_nbrs(o_ester)
                if any(x.GetSymbol() == "O" and x.GetIdx() != o_ester.GetIdx()
                       for x in heavy_nbrs(n)))
    o_carbonyl = next(x for x in heavy_nbrs(acyl)
                      if x.GetSymbol() == "O" and x.GetIdx() != o_ester.GetIdx())
    methyl = next(x for x in heavy_nbrs(acyl) if x.GetSymbol() == "C")
    maps = {o_ester.GetIdx(): 31, o_carbonyl.GetIdx(): 32, acyl.GetIdx(): 22,
            methyl.GetIdx(): 21}
    prev, cur, n = o_ester, next(x for x in heavy_nbrs(o_ester) if x.GetIdx() != acyl.GetIdx()), 1
    while cur is not None:
        maps[cur.GetIdx()] = n
        nxt = [x for x in heavy_nbrs(cur) if x.GetIdx() != prev.GetIdx()]
        assert len(nxt) <= 1, "the chain branches: not an n-alkenyl acetate"
        prev, cur, n = cur, (nxt[0] if nxt else None), n + 1
    assert sorted(m for m in maps.values() if m < 21) == list(range(1, 19)), maps
    return maps


def write_input(src, dst, transform):
    """The fixture's input file from its bench/fixtures source: maps renumbered
    along the chain ("chain"), or rotated out of the z = 0 plane ("tilt")."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    mol = Chem.MolFromMolFile(str(src), removeHs=False)
    for k, v in Chem.SDMolSupplier(str(src), removeHs=False)[0].GetPropsAsDict().items():
        mol.SetProp(k, str(v))
    if transform == "chain":
        maps = chain_maps(mol)
        for a in mol.GetAtoms():
            a.SetAtomMapNum(maps.get(a.GetIdx(), 0))
        mol.SetProp("atom_maps", "chain numbering: map n = Cn; acetyl CH3 21, C=O 22; "
                                 "ester O 31, carbonyl O 32 (fixtures/make_records.py)")
        # an ETKDG conformer minimised with MMFF94s: the same method label a
        # --smiles run gives it, so the two share an analysis_key (review 5 R1/R6)
        mol.SetProp("espmap_geometry_origin", "etkdg_mmff94s")
    elif transform == "tilt":
        conf = mol.GetConformer()
        for i, p in enumerate(tilted(conf.GetPositions())):
            conf.SetAtomPosition(i, Point3D(*p))
        mol.SetProp("frame", f"{src.name} rotated rigidly by 35 deg about {TILT} so no "
                             "plane is z = 0 (DIMENSION_AMBIGUOUS); internal coordinates "
                             "unchanged (fixtures/make_records.py)")
    w = Chem.SDWriter(str(dst))
    w.write(mol)
    w.close()


def argv(name, fx, out, basis=None):
    a = [sys.executable, str(ROOT / "esp_map.py"), "--file", str(fx["file"]),
         "--basis", basis or fx["basis"], "--yes", "--size", "500",
         "--out", str(out / f"{name}.png"), "--json", str(out / f"{name}.json"),
         "--npz", str(out / f"{name}.npz")]
    if fx["solvent"]:
        a += ["--solvent", fx["solvent"]]
    if fx["density_fit"]:
        a.append("--density-fit")
    if fx["sites"]:
        a += ["--sites", str(GEOM / fx["sites"])]
    return a


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--basis", help="override every fixture's basis (dry runs only)")
    p.add_argument("--out", default=str(HERE / "records"))
    p.add_argument("names", nargs="*", help="fixtures to make (default: all)")
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    GEOM.mkdir(exist_ok=True)
    (GEOM / C18["sites"]).write_text(json.dumps(C18_SITES, indent=1) + "\n")
    for name, fx in FIXTURES.items():
        if "source" in fx:              # deterministic, so rewriting is a no-op unless changed
            write_input(fx["source"], fx["file"], fx["transform"])
    # Review 5 R17: one set of fixtures, one code revision. Skipping a record
    # that exists used to mix revisions on a partial re-run, and a run that
    # cannot name its code (outside git) or runs edited code is refused.
    sys.path.insert(0, str(ROOT))
    import esp_map as em
    code = em.code_identity()["code"]
    if a.basis is None and (code == "unknown" or code.endswith("-dirty")):
        sys.exit(f"code identity is {code!r}: run from a clean git checkout, so every "
                 f"record names the code that made it")
    failed = []
    for name in a.names or FIXTURES:
        js = out / f"{name}.json"
        if js.exists():
            old = json.loads(js.read_text())["provenance"]["code"]
            if old != code:
                sys.exit(f"{js.name} was made by {old}, this is {code}: delete the old "
                         f"records so the set comes from one revision")
            print(f"{name}: exists, made by this revision, skipped", flush=True)
            continue
        print(f"{name}: computing", flush=True)
        if subprocess.run(argv(name, FIXTURES[name], out, a.basis)).returncode:
            failed.append(name)
    if failed:
        sys.exit(f"failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
