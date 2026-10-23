"""The ensemble runner and espmap.ensemble/1 (docs/ensemble-spec.md section 9,
step 4), without an SCF: member runs are replaced by a stub that writes what
a member run writes (record and NPZ, under the *.part names it is given),
so resume, failure, exclusion and the statistics are tested in seconds. The
real end-to-end cases are the hostile set's slow ones."""
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

import esp_map as em

_SCHEMA = json.loads((Path(__file__).resolve().parents[1] / "docs" / "ensemble-schema.json")
                     .read_text(encoding="utf-8"))


def _arg(argv, flag):
    return argv[argv.index(flag) + 1]


class Stub:
    """A member run. V_min is -10 - (member number), so statistics are known;
    fail/invalid/key name the conformers to break."""

    def __init__(self, fail=(), invalid=(), key_for=None, kill_after_json=()):
        self.calls, self.fail, self.invalid = [], set(fail), set(invalid)
        self.key_for, self.kill = key_for or {}, set(kill_after_json)

    def __call__(self, argv, cwd, env, timeout):
        sdf = Path(_arg(argv, "--file"))
        label = sdf.stem
        self.calls.append((label, env["PYSCF_TMPDIR"]))
        if label in self.fail:
            return 1, 2.0, "Traceback ...\nRuntimeError: SCF not converged"
        n = int(label.split("_")[1])
        m = Chem.MolFromMolFile(str(sdf), removeHs=False)
        sha1 = em.spec_sha1(em._mol_to_spec(m))
        rec = dict(schema="espmap.record/1", analysis_key=self.key_for.get(label, "a" * 40),
                   key_fields=dict(analysis=dict(xc="b3lyp", basis="def2svp")),
                   validity=dict(input=True, surface=label not in self.invalid,
                                 electronic=True, sites=label not in self.invalid),
                   geometry=dict(sha1=sha1), settings=dict(optimize=False),
                   provenance=dict(espmap_version=em.__version__, code="0123abc",
                                   rdkit="x", pyscf="y"),
                   results=dict(e_tot_Ha=-100.0 - 0.001 * n, V_min_kcal_mol=-10.0 - n,
                                V_max_kcal_mol=20.0 + n, surface=dict(genus=0)),
                   sites=dict(atoms=[dict(index=0, atom_map=0, element="C",
                                          nearest_nucleus=dict(V_area_mean_kcal_mol=1.0 + n),
                                          radius_scaled=dict(V_area_mean_kcal_mol=2.0 + n))],
                              reaction_sites=[]))
        Path(_arg(argv, "--json")).write_text(json.dumps(rec), encoding="utf-8")
        if label in self.kill:
            return -9, 1.0, "killed"
        np.savez_compressed(_arg(argv, "--npz"), geometry_sha1=sha1)
        return 0, 3.0, ""


@pytest.fixture(scope="module")
def plan():
    return em.plan_ensemble(smiles="CCCCCCO", ensemble=4)


def test_members_are_written_run_and_indexed(plan, tmp_path):
    stub = Stub()
    rec = em.run_ensemble(plan, tmp_path, ["--basis", "def2-svp"], runner=stub, log=lambda s: None)
    labels = [m["conformer"] for m in rec["members"]]
    assert labels == ["conf_00", "conf_01", "conf_02", "conf_03"]          # zero-based
    assert len({tmp for _, tmp in stub.calls}) == 4                         # own PYSCF_TMPDIR each
    for m in rec["members"]:
        js = tmp_path / m["record"]
        assert m["record_sha256"] == __import__("hashlib").sha256(js.read_bytes()).hexdigest()
        assert m["geometry_sha1"] == json.loads(js.read_text())["geometry"]["sha1"]
        assert (tmp_path / m["npz"]).exists()
    # the SDFs carry the geometry origin, so it reaches analysis_key
    sdf = next(iter(Chem.SDMolSupplier(str(tmp_path / "conformers" / "conf_00.sdf"), removeHs=False)))
    assert sdf.GetProp("espmap_geometry_origin") == "etkdg_mmff94s"
    assert json.loads((tmp_path / "ensemble.json").read_text()) == rec
    log = [json.loads(x) for x in (tmp_path / "run_log.jsonl").read_text().splitlines()]
    assert [x["conformer"] for x in log] == labels and all(x["rc"] == 0 for x in log)


def test_a_finished_member_is_not_recomputed(plan, tmp_path):
    em.run_ensemble(plan, tmp_path, runner=Stub(), log=lambda s: None)
    again = Stub()
    em.run_ensemble(plan, tmp_path, runner=again, log=lambda s: None)
    assert again.calls == []


def test_a_run_killed_after_its_json_is_failed_and_rerun(plan, tmp_path):
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(kill_after_json={"conf_01"}),
                          log=lambda s: None)
    assert not (tmp_path / "conformers" / "conf_01.json").exists()      # never renamed in
    assert [e["conformer"] for e in rec["excluded"]] == ["conf_01"]
    again = Stub()
    rec = em.run_ensemble(plan, tmp_path, runner=again, log=lambda s: None)
    assert [c for c, _ in again.calls] == ["conf_01"] and rec["excluded"] == []


def test_a_failed_member_is_excluded_with_its_shape_and_reason(plan, tmp_path):
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(fail={"conf_02"}), log=lambda s: None)
    (e,) = rec["excluded"]
    assert e["conformer"] == "conf_02" and e["code"] == "MEMBER_RUN_FAILED"
    assert "SCF not converged" in e["reason"]
    assert e["rg_angstrom"] == plan["members"][2]["rg"] and e["folded"] == plan["members"][2]["folded"]
    assert rec["k_obtained"] == 3 and rec["validity"]["members"] is False
    assert rec["validity"]["unbiased"] is (not plan["members"][2]["folded"])


def test_an_excluded_folded_member_makes_the_ensemble_biased(plan, tmp_path):
    folded = [i for i, m in enumerate(plan["members"]) if m["folded"]]
    if not folded:
        pytest.skip("no folded member in this plan")
    label = em.member_label(folded[0])
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(invalid={label}), log=lambda s: None)
    assert [x["code"] for x in rec["excluded"]] == ["MEMBER_INVALID"]
    assert rec["validity"]["unbiased"] is False


def test_members_under_two_protocols_are_refused(plan, tmp_path):
    with pytest.raises(em.EnsembleRefused) as e:
        em.run_ensemble(plan, tmp_path, runner=Stub(key_for={"conf_03": "b" * 40}),
                        log=lambda s: None)
    assert e.value.findings[0][1] == "ENSEMBLE_MIXED_ANALYSIS"


def test_statistics_over_all_and_extended_members(plan, tmp_path):
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(), log=lambda s: None)
    vmin = rec["results"]["V_min_kcal_mol"]["all"]
    assert vmin["mean"] == pytest.approx(-11.5) and vmin["range"] == pytest.approx(3.0)
    assert vmin["sd"] == pytest.approx(np.std([-10, -11, -12, -13], ddof=1))
    ext = [i for i, m in enumerate(plan["members"]) if not m["folded"]]
    assert rec["results"]["V_min_kcal_mol"]["extended"]["n"] == len(ext)
    site = rec["sites"]["atom 0 map 0 C"]["nearest_nucleus"]
    assert site["all"]["mean"] == pytest.approx(2.5) and len(site["all"]["bootstrap95_mean"]) == 2
    c = rec["convergence"]["per_site"]["atom 0 map 0 C | nearest_nucleus"]
    assert c["running_mean"] == pytest.approx([1.0, 1.5, 2.0, 2.5])
    assert c["shift_last_third"] == pytest.approx(0.5) and c["converged"] is True


def test_dft_weights_report_the_effective_number(tmp_path):
    plan = em.plan_ensemble(smiles="CCCCCCO", ensemble=4, weight="dft", optimize=True)
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(), log=lambda s: None)
    w = [m["weight"] for m in rec["members"]]
    assert sum(w) == pytest.approx(1) and 1 <= rec["kish_n"] <= 4
    assert rec["kish_n"] == pytest.approx(1 / sum(x * x for x in w))


def test_the_record_matches_the_schema_both_ways(plan, tmp_path):
    """Every required key is written and every written key is declared, at
    the top level and per member (no jsonschema in the venv)."""
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(), log=lambda s: None)
    for obj, sch in ((rec, _SCHEMA), (rec["members"][0], _SCHEMA["properties"]["members"]["items"]),
                     (rec["identity"], _SCHEMA["properties"]["identity"]),
                     (rec["protocol"], _SCHEMA["properties"]["protocol"]),
                     (rec["protocol"]["selection"],
                      _SCHEMA["properties"]["protocol"]["properties"]["selection"]),
                     (rec["validity"], _SCHEMA["properties"]["validity"])):
        assert set(sch["required"]) <= set(obj), set(sch["required"]) - set(obj)
        assert set(obj) <= set(sch["properties"]), set(obj) - set(sch["properties"])
    assert rec["schema"] == "espmap.ensemble/1"
    assert rec["espmap"]["commit"] == "0123abc"


def test_a_supplied_ensemble_keeps_the_file_numbers(tmp_path):
    plan = em.plan_ensemble(path=Path(__file__).parent / "data" / "ensemble" /
                            "E7_poses_one_ligand.sdf", ensemble="supplied")
    plan = dict(plan, members=plan["members"][:3])
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(), log=lambda s: None, subject="spa.cmpd.x")
    assert [m["source_record"] for m in rec["members"]] == [1, 2, 3]
    assert rec["protocol"]["selection"]["dedup_atoms"] is None and rec["subject"] == "spa.cmpd.x"
