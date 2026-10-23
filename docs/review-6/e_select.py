"""Review 6, PR #73: select_ensemble, shape_descriptors and the contact rule
against the measured C18 data. No SCF. usage: e_select.py PR_DIR [part ...]"""
import hashlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, sys.argv[1])
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D, rdMolTransforms as T

import esp_map as em

PARTS = set(sys.argv[2:]) or {"c18", "matches", "k2", "det", "kink", "tempo"}
BENCH = Path(sys.argv[1]) / "bench"
em._quiet_rdkit()


def c18():
    smi = "CC(=O)OCC/C=C\\CCCCCCCC/C=C\\CCCC"
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    t = time.perf_counter()
    best, ff, spread, n_emb, n_conv, conv = em._embed_and_minimise(mol, em.N_EMBED, 0xC0FFEE)
    print(f"embedded {n_emb}, converged {n_conv}, {time.perf_counter() - t:.1f} s")
    # the bench's genus per conf_id (conformers.jsonl chi, f9_recheck genus)
    rows = [json.loads(l) for l in (BENCH / "conformers.jsonl").read_text().splitlines() if l.strip()]
    bench = {r["conf_id"]: r for r in rows if "conf_id" in r}
    for l in (BENCH / "f9_recheck.jsonl").read_text().splitlines():
        if l.strip():
            r = json.loads(l)
            bench[r["conf_id"]].update(genus=r["genus"], chi=r["chi"])
    print("conf  bench_rg  rg_mass_allatoms  genus  folded  n_contacts  min_bonds  closest(d - limit)")
    topo = Chem.GetDistanceMatrix(mol)
    for cid in sorted(bench):
        b = bench[cid]
        rg_all = Descriptors3D.RadiusOfGyration(mol, confId=cid)
        genus = b.get("genus", (2 - b["chi"]) // 2 if "chi" in b else None)
        s = em.shape_descriptors(mol, cid)
        con = em.intramolecular_contacts(mol, cid)
        mb = min((int(topo[i, j]) for i, j, *_ in con), default=None)
        worst = min((d - lim for *_, d, lim in con), default=None)
        print(f"{cid:4d}  {b['rg']:.4f}    {rg_all:.4f}            {genus}      {s['folded']!s:5}  "
              f"{len(con):3d}         {mb}          {worst if worst is None else round(worst, 2)}")
    # the whole embedded population: folded rate, and how local the contacts are
    allc = []
    for cid, e in conv:
        con = em.intramolecular_contacts(mol, cid)
        allc.append((cid, e, con))
    folded = [x for x in allc if x[2]]
    only_local = [x for x in folded if max(int(topo[i, j]) for i, j, *_ in x[2]) <= 6]
    print(f"population: {len(folded)}/{len(allc)} converged conformers 'folded'; "
          f"{len(only_local)} of them have every contact only 5-6 bonds apart")
    # what the PR selects vs what the bench selected
    t = time.perf_counter()
    sel = em.select_ensemble(mol, conv, 10)
    print(f"select_ensemble K=10: {time.perf_counter() - t:.1f} s; n_unique={sel['n_unique']}, "
          f"symmetry matches={sel['n_symmetry_matches']}")
    pr_ids = [c for c, _ in sel["members"]]
    print(f"PR members (conf ids, selection order): {pr_ids}")
    print(f"bench members (conformers.jsonl):        {sorted(bench)}")
    print(f"shared: {sorted(set(pr_ids) & set(bench))}")
    sys.path.insert(0, str(BENCH))
    import conformers as bc
    cids = [c.GetId() for c in mol.GetConformers()]
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    en = np.array([AllChem.MMFFGetMoleculeForceField(mol, props, confId=c).CalcEnergy() for c in cids])
    rg = np.array([Descriptors3D.RadiusOfGyration(mol, confId=c) for c in cids])
    print(f"bench select() re-run here: {sorted(cids[i] for i in bc.select(rg, en, 10))}")
    print("PR members: folded flags", [em.shape_descriptors(mol, c)["folded"] for c in pr_ids])


def matches():
    """The symmetry-match cap. (tBuOCH2)3CH: 3! * 6^3 = 1296 heavy-atom automorphisms."""
    smi = "CC(C)(C)OCC(COC(C)(C)C)COC(C)(C)C"
    heavy = Chem.MolFromSmiles(smi)
    n_all = len(heavy.GetSubstructMatches(heavy, uniquify=False, useChirality=False, maxMatches=100000))
    print(f"{smi}: {n_all} heavy-atom automorphisms; cap {1000}")
    mol = Chem.AddHs(heavy)
    *_, conv = em._embed_and_minimise(mol, 30, 0xC0FFEE)
    # the same conformer under a symmetry relabelling: rotate one tBu's methyls
    h = Chem.RemoveHs(mol)
    from rdkit.Chem import rdMolAlign
    cid = conv[0][0]
    perm = list(h.GetSubstructMatches(h, uniquify=False, useChirality=False, maxMatches=100000))
    missing = [m for m in perm if m not in
               h.GetSubstructMatches(h, uniquify=False, useChirality=False, maxMatches=1000)]
    print(f"automorphisms outside the first 1000: {len(missing)}")
    m = missing[-1]
    probe = Chem.Mol(h)
    c = probe.GetConformer(cid)
    P = h.GetConformer(cid).GetPositions()
    for q, t in enumerate(m):
        c.SetAtomPosition(t, P[q].tolist())
    maps_cap = [list(zip(mm, range(h.GetNumAtoms())))
                for mm in h.GetSubstructMatches(h, uniquify=False, useChirality=False, maxMatches=1000)]
    maps_all = [list(zip(mm, range(h.GetNumAtoms()))) for mm in perm]
    r_cap = rdMolAlign.GetBestRMS(Chem.Mol(probe), h, cid, cid, map=maps_cap)
    r_all = rdMolAlign.GetBestRMS(Chem.Mol(probe), h, cid, cid, map=maps_all)
    print(f"identical conformer, relabelled by an automorphism past the cap: "
          f"RMSD {r_cap:.3f} A with the cap, {r_all:.3f} A with every automorphism "
          f"(DEDUP_RMSD {em.DEDUP_RMSD})")
    for cap in (1000, 2000):
        t = time.perf_counter()
        s = em.select_ensemble(mol, conv, 30, max_matches=cap)
        print(f"  30 conformers, max_matches={cap}: n_unique {s['n_unique']}, "
              f"{time.perf_counter() - t:.1f} s")


def k2():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCCCCCO"))
    *_, conv = em._embed_and_minimise(mol, em.N_EMBED, 0xC0FFEE)
    every = em.select_ensemble(mol, conv, em.N_EMBED)
    rg = {c: em.shape_descriptors(mol, c)["rg"] for c, _ in every["members"]}
    print(f"octanol: {every['n_unique']} unique; Rg {min(rg.values()):.2f}-{max(rg.values()):.2f}")
    for k in (2, 3):
        s = em.select_ensemble(mol, conv, k)
        print(f"  K={k}: Rg of members {[round(rg[c], 2) for c, _ in s['members']]}")


def det():
    p = em.plan_ensemble(smiles="CCCCCCCCO", ensemble=8)
    h = hashlib.sha1("|".join(m["spec"] for m in p["members"]).encode()).hexdigest()[:12]
    print(f"octanol K=8 members sha1 {h}  ids {[m['embed_conf_id'] for m in p['members']]}")


def kink():
    """Idealised all-anti C18 chain with one local double-gauche motif, MMFF-free:
    which motifs does the contact rule call folded?"""
    smi = "CCCCCCCCCCCCCCCCCC"
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    conf = mol.GetConformer()
    chain = list(range(18))
    for motif, angles in (("all anti", {}), ("g+g-", {7: 60, 8: -60}), ("g+g+", {7: 60, 8: 60}),
                          ("g+g+g+", {7: 60, 8: 60, 9: 60}), ("g+ t g+", {7: 60, 9: 60}),
                          ("g+ t g-", {7: 60, 9: -60}), ("g+g+g+g+", {7: 60, 8: 60, 9: 60, 10: 60})):
        for k in range(1, 16):
            T.SetDihedralDeg(conf, chain[k - 1], chain[k], chain[k + 1], chain[k + 2],
                             float(angles.get(k, 180)))
        con = em.intramolecular_contacts(mol)
        topo = Chem.GetDistanceMatrix(mol)
        print(f"  {motif:10s} folded={bool(con)!s:5}  contacts "
              f"{[(int(topo[i, j]), round(d, 2)) for i, j, d, _ in con][:4]}  "
              f"Rg {em.shape_descriptors(mol)['rg']:.2f}")


def tempo():
    """Three redox states of one nitroxide, same nuclei: radical, oxoammonium, aminoxide."""
    m = Chem.AddHs(Chem.MolFromSmiles("CC1(C)CCCC(C)(C)N1[O]"))
    AllChem.EmbedMolecule(m, randomSeed=3)
    AllChem.MMFFOptimizeMolecule(m)
    out = []
    for smi in ("CC1(C)CCCC(C)(C)[N+]1=O", "CC1(C)CCCC(C)(C)N1[O]", "CC1(C)CCCC(C)(C)N1[O-]"):
        t = Chem.AddHs(Chem.MolFromSmiles(smi))
        assert [a.GetAtomicNum() for a in t.GetAtoms()] == [a.GetAtomicNum() for a in m.GetAtoms()]
        c = Chem.Conformer(m.GetConformer())
        c.Set3D(True)
        t.AddConformer(c, assignId=True)
        out.append(t)
    p = Path(sys.argv[1]).parent / "espmap-review6" / "docs" / "review-6" / "fx" / "ens" / "s3e_tempo_redox.sdf"
    w = Chem.SDWriter(str(p))
    for t in out:
        w.write(t)
    w.close()
    try:
        r = em.plan_ensemble(path=p, ensemble="supplied")
        print(f"TEMPO+ / TEMPO. / TEMPO-: ACCEPTED, {len(r['members'])} members, identity {r['identity']['smiles']}")
    except em.EnsembleRefused as e:
        print("TEMPO: REFUSED", [c for _, c, *_ in e.findings])


for name in ("tempo", "kink", "k2", "det", "matches", "c18"):
    if name in PARTS:
        print(f"===== {name}")
        globals()[name]()
        sys.stdout.flush()
