"""docs/ensemble-spec.md section 8, the parts that need no SCF.

Three provisional constants of the ensemble code, measured with RDKit only
(the embedder and selection are espmap's own, invariant 1):

  dedup     DEDUP_RMSD (0.5 A). Ground truth for chains: the conformer's
            backbone torsion string (a, g+, g-), canonical under reversal
            and, for an achiral molecule, under reflection. For each
            threshold, count the distinct conformers select_ensemble keeps
            and compare with the distinct torsion classes among the
            converged conformers: a FALSE MERGE is a class lost, a FALSE
            SPLIT a class kept twice.
  respread  _rg_spread (bench/conformers.py's rule) against re-spreading
            after a collision with the lowest-energy conformer (review 6
            S4): which members each picks, and the largest Rg gap left.
  contacts  the folded rule (CONTACT_MIN_BONDS, CONTACT_SLACK) against the
            measured surface genus of the 10 bench C18 conformers
            (conformers.jsonl chi; f9_recheck.jsonl genus), and the folded
            rate over all 100 (review 6 S3).

    python bench/ensemble_rdkit.py          # appends to bench/ensemble_rdkit.jsonl

Rows already present (same kind and key) are skipped; add rows, never
overwrite them.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import esp_map as em  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Descriptors3D, rdMolTransforms  # noqa: E402

OUT = HERE / "ensemble_rdkit.jsonl"
SEED = 0xC0FFEE
THRESHOLDS = (0.05, 0.1, 0.15, 0.25, 0.35, 0.5, 0.75)
C18 = "CC(=O)OCC/C=C\\CCCCCCCC/C=C\\CCCC"


def done():
    if not OUT.exists():
        return set()
    return {(r["kind"], r["key"]) for r in map(json.loads, OUT.read_text().splitlines())}


def emit(row):
    with OUT.open("a") as f:
        f.write(json.dumps(row) + "\n")
    print(row["kind"], row["key"], {k: v for k, v in row.items()
                                    if k in ("n_classes", "n_unique", "false_merges",
                                             "false_splits", "agree", "folded_rate")})


# ------------------------------------------------------------------- dedup
def chain(mol):
    """The longest heavy-atom path: the backbone whose torsions classify."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    dm = Chem.GetDistanceMatrix(mol)
    i, j = max(((i, j) for i in heavy for j in heavy), key=lambda p: dm[p])
    return list(Chem.GetShortestPath(mol, i, j))


def torsion_class(mol, cid, path, achiral):
    conf = mol.GetConformer(cid)
    s = []
    for k in range(len(path) - 3):
        d = rdMolTransforms.GetDihedralDeg(conf, *path[k:k + 4])
        s.append("a" if abs(d) > 120 else ("+" if d > 0 else "-"))
    s = tuple(s)
    forms = {s, s[::-1]}
    if achiral:
        flip = {"+": "-", "-": "+", "a": "a"}
        forms |= {tuple(flip[x] for x in f) for f in list(forms)}
    return min(forms)


def dedup_rows(name, smiles, n_embed=100):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    *_, conv = em._embed_and_minimise(mol, n_embed, SEED)
    path = chain(mol)
    achiral = em._is_achiral(mol)
    cls = {c: torsion_class(mol, c, path, achiral) for c, _ in conv}
    truth = set(cls.values())
    for t in THRESHOLDS:
        key = f"{name} {t}"
        if ("dedup", key) in done():
            continue
        sel = em.select_ensemble(mol, conv, len(conv), dedup_rmsd=t)
        kept = [cls[c] for c, _ in sel["members"]]
        emit(dict(kind="dedup", key=key, molecule=name, smiles=smiles, threshold_A=t,
                  n_embed=n_embed, n_converged=len(conv), n_classes=len(truth),
                  n_unique=sel["n_unique"], false_merges=len(truth - set(kept)),
                  false_splits=len(kept) - len(set(kept)), backbone_torsions=len(path) - 3,
                  mirror_merged=achiral, seed=SEED))


# ---------------------------------------------------------------- respread
def respread(rg, energy, k):
    """Review 6 S4's alternative: when a pick collides, re-spread the k-1
    over the conformers other than the lowest-energy one."""
    best = int(np.argmin(energy))
    rest = [int(i) for i in np.argsort(rg, kind="stable") if int(i) != best]
    picks = [best]
    for q in np.linspace(0, len(rest) - 1, k - 1):
        i = rest[int(round(q))]
        if i not in picks:
            picks.append(i)
    return picks


def max_gap(rg, picks):
    v = np.sort(np.asarray(rg)[picks])
    return float(np.max(np.diff(v))) if len(v) > 1 else None


def respread_rows():
    for name, smiles, k in (("C18 acetate", C18, 10), ("1-octanol", "CCCCCCCCO", 3),
                            ("1-octanol", "CCCCCCCCO", 5), ("1-hexanol", "CCCCCCO", 6)):
        key = f"{name} K={k}"
        if ("respread", key) in done():
            continue
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        *_, conv = em._embed_and_minimise(mol, 100, SEED)
        sel = em.select_ensemble(mol, conv, len(conv))
        kept = sel["members"]
        rg = [Descriptors3D.RadiusOfGyration(mol, confId=c) for c, _ in kept]
        e = [x for _, x in kept]
        a, b = em._rg_spread(rg, e, k), respread(rg, e, k)
        emit(dict(kind="respread", key=key, molecule=name, k=k, n_unique=len(kept),
                  bench_rule=[kept[i][0] for i in a], respread=[kept[i][0] for i in b],
                  collided=len(set(a) ^ set(b)) > 0,
                  bench_max_gap_A=max_gap(rg, a), respread_max_gap_A=max_gap(rg, b),
                  rg_range_A=float(max(rg) - min(rg)), seed=SEED))


# ---------------------------------------------------------------- contacts
def measured_genus():
    g = {}
    for r in map(json.loads, (HERE / "conformers.jsonl").read_text().splitlines()):
        if "conf_id" in r and r.get("chi") is not None:
            g[r["conf_id"]] = (2 - r["chi"]) // 2
    for r in map(json.loads, (HERE / "f9_recheck.jsonl").read_text().splitlines()):
        g[r["conf_id"]] = r["genus"]
    return g


def contact_rows():
    mol = Chem.AddHs(Chem.MolFromSmiles(C18))
    *_, conv = em._embed_and_minimise(mol, 100, SEED)
    genus = measured_genus()
    topo = Chem.GetDistanceMatrix(mol)
    real_min, real_slack = em.CONTACT_MIN_BONDS, em.CONTACT_SLACK
    for min_bonds in (5, 6, 7, 8):
        for slack in (0.3, 0.5):
            key = f"C18 min_bonds={min_bonds} slack={slack}"
            if ("contacts", key) in done():
                continue
            em.CONTACT_MIN_BONDS, em.CONTACT_SLACK = min_bonds, slack
            try:
                folded = {c: bool(em.intramolecular_contacts(mol, c)) for c, _ in conv}
                bonds = [int(topo[i, j]) for c, _ in conv for i, j, *_ in em.intramolecular_contacts(mol, c)]
            finally:
                em.CONTACT_MIN_BONDS, em.CONTACT_SLACK = real_min, real_slack
            per = {str(c): dict(genus=genus[c], folded=folded[c]) for c in sorted(genus)}
            emit(dict(kind="contacts", key=key, molecule="Z3Z13-18OAc", min_bonds=min_bonds,
                      slack_A=slack, measured=per,
                      agree=sum(v["folded"] == (v["genus"] > 0) for v in per.values()),
                      n_measured=len(per), folded_rate=sum(folded.values()) / len(folded),
                      n_conformers=len(folded),
                      contact_bond_separation=dict(zip(*np.unique(bonds, return_counts=True)))
                      if bonds else {}, seed=SEED))


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


if __name__ == "__main__":
    _emit = emit
    emit = lambda row: _emit(_jsonable(row))  # noqa: E731
    for name, smiles in (("n-butane", "CCCC"), ("n-pentane", "CCCCC"), ("n-hexane", "CCCCCC"),
                         ("n-heptane", "CCCCCCC"), ("1-hexanol", "CCCCCCO"),
                         ("cyclohexane", "C1CCCCC1")):
        dedup_rows(name, smiles)
    respread_rows()
    contact_rows()
