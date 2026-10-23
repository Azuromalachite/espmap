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

    def __init__(self, fail=(), invalid=(), key_for=None, kill_after_json=(), site=None,
                 optimize=False):
        self.optimize = optimize                       # write an optimised geometry, as --optimize does
        self.calls, self.fail, self.invalid = [], set(fail), set(invalid)
        self.key_for, self.kill = key_for or {}, set(kill_after_json)
        self.site = site or (lambda n: 1.0 + n)       # member number -> nearest-nucleus mean

    def __call__(self, argv, cwd, env, timeout):
        sdf = Path(_arg(argv, "--file"))
        label = sdf.stem
        self.calls.append((label, env["PYSCF_TMPDIR"]))
        if label in self.fail:
            return 1, 2.0, "Traceback ...\nRuntimeError: SCF not converged"
        n = int(label.split("_")[1])
        m = Chem.MolFromMolFile(str(sdf), removeHs=False)
        sha1 = em.spec_sha1(em._mol_to_spec(m))
        # --optimize: the record's geometry is the optimised one, and it names
        # the input it was optimised from
        geom = (dict(sha1=em._sha1(sha1 + "optimised"), optimised_from_sha1=sha1) if self.optimize
                else dict(sha1=sha1))
        rec = dict(schema="espmap.record/1", analysis_key=self.key_for.get(label, "a" * 40),
                   key_fields=dict(analysis=dict(xc="b3lyp", basis="def2svp")),
                   validity=dict(input=True, surface=label not in self.invalid,
                                 electronic=True, sites=label not in self.invalid),
                   geometry=geom, settings=dict(optimize=self.optimize),
                   provenance=dict(espmap_version=em.__version__, code="0123abc",
                                   rdkit="x", pyscf="y"),
                   results=dict(e_tot_Ha=-100.0 - 0.001 * n, V_min_kcal_mol=-10.0 - n,
                                V_max_kcal_mol=20.0 + n, surface=dict(genus=0)),
                   sites=dict(atoms=[dict(index=0, atom_map=0, element="C",
                                          nearest_nucleus=dict(V_area_mean_kcal_mol=self.site(n)),
                                          radius_scaled=dict(V_area_mean_kcal_mol=self.site(n) + 1.0))],
                              reaction_sites=[]))
        Path(_arg(argv, "--json")).write_text(json.dumps(rec), encoding="utf-8")
        if label in self.kill:
            return -9, 1.0, "killed"
        np.savez_compressed(_arg(argv, "--npz"), geometry_sha1=geom["sha1"])
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


def test_a_deterministic_failure_is_not_retried_on_a_rerun(plan, tmp_path):
    """Rerunning until a member passes would be selection on outcome."""
    em.run_ensemble(plan, tmp_path, runner=Stub(fail={"conf_02"}), log=lambda s: None)
    again = Stub()
    rec = em.run_ensemble(plan, tmp_path, runner=again, log=lambda s: None)
    assert again.calls == []
    (e,) = rec["excluded"]
    assert e["conformer"] == "conf_02" and e["attempts"] == 1
    assert e["earlier_failures"][0].startswith("deterministic, rc 1")


def test_an_infrastructure_failure_is_retried_and_its_attempts_kept(plan, tmp_path):
    em.run_ensemble(plan, tmp_path, runner=Stub(kill_after_json=set()), log=lambda s: None)
    # a kill with no record written: the stub's -9 path writes the JSON first,
    # so model the OOM kill directly
    class OOM(Stub):
        def __call__(self, argv, cwd, env, timeout):
            if Path(_arg(argv, "--file")).stem == "conf_01":
                return -9, 50.0, "Killed"
            return super().__call__(argv, cwd, env, timeout)
    for f in (tmp_path / "conformers").glob("conf_01.*json"):
        f.unlink()
    rec = em.run_ensemble(plan, tmp_path, runner=OOM(), log=lambda s: None)
    assert [e["conformer"] for e in rec["excluded"]] == ["conf_01"]
    again = Stub()
    rec = em.run_ensemble(plan, tmp_path, runner=again, log=lambda s: None)
    assert [c for c, _ in again.calls] == ["conf_01"] and rec["excluded"] == []
    m = next(m for m in rec["members"] if m["conformer"] == "conf_01")
    assert m["attempts"] == 3 and len(m["earlier_failures"]) == 1
    assert m["earlier_failures"][0].startswith("infrastructure, rc -9")


@pytest.mark.parametrize("rc, err, kind", [
    (0, "", "ok"),
    (-9, "Killed", "infrastructure"),
    (137, "", "infrastructure"),
    (1, "MemoryError", "infrastructure"),
    (1, "structure rejected", "deterministic"),
    (1, "surface not closed", "deterministic"),
])
def test_failure_kinds(rc, err, kind):
    assert em.failure_kind(rc, err) == kind


# --- the --ensemble command line ---------------------------------------------------
def test_member_runs_get_the_users_settings_and_not_the_ensemble_flags():
    argv = ["--smiles", "CCO", "--ensemble", "4", "--basis", "def2-svp", "--density-fit",
            "--ensemble-dir=out", "--subject", "x", "--seed", "7", "--solvent", "water",
            "--weight", "dft", "--optimize", "--check", "--sites", "s.json", "--json", "r.json"]
    assert em.member_argv(argv) == ["--basis", "def2-svp", "--density-fit", "--solvent",
                                    "water", "--optimize", "--sites", "s.json"]


def _cli(*args, cwd):
    import subprocess
    import sys
    r = subprocess.run([sys.executable, em.__file__, *args], cwd=cwd, capture_output=True,
                       text=True, timeout=300)
    return r.returncode, r.stdout + r.stderr


@pytest.mark.parametrize("args, says", [
    (["--smiles", "CC(O)C(N)C(=O)O", "--ensemble", "4"], "STEREO_REQUIRED"),
    (["--smiles", "CCCCCCO", "--ensemble", "4", "--weight", "mmff"], "WEIGHT_MMFF"),
    (["--smiles", "CCO", "--ensemble", "2", "--force"], "--force does not apply"),
    (["--smiles", "CCO", "--ensemble", "two"], "give a number of conformers"),
    (["--atoms", "O 0 0 0; H 0 0 1; H 0 1 0", "--ensemble", "2"], "needs --smiles or --file"),
])
def test_the_command_line_refuses_before_computing(args, says, tmp_path):
    rc, out = _cli(*args, cwd=tmp_path)
    assert rc != 0 and says in out
    assert not (tmp_path / "ensemble").exists()


def test_check_plans_without_computing(tmp_path):
    rc, out = _cli("--smiles", "CCCCCCO", "--ensemble", "4", "--check", cwd=tmp_path)
    assert rc == 0 and "4 member(s)" in out and not (tmp_path / "ensemble").exists()


def test_a_supplied_file_refusal_names_the_record(tmp_path):
    rc, out = _cli("--file", str(Path(__file__).parent / "data" / "ensemble" /
                                 "E8_ligand_and_tautomer.sdf"), "--ensemble", "supplied",
                   cwd=tmp_path)
    assert rc != 0 and "MEMBER_IDENTITY" in out and "record 4" in out


def test_the_record_states_the_measured_thresholds(plan, tmp_path):
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(), log=lambda s: None)
    fr, sel = rec["folded_rule"], rec["protocol"]["selection"]
    assert (fr["min_bonds"], fr["slack_A"], fr["status"]) == (6, 0.5, "measured")
    assert fr["source"] == sel["dedup_source"] == "docs/ensemble-spec.md 8.2"
    assert (sel["dedup_rmsd_A"], sel["dedup_status"]) == (0.35, "measured")


def test_the_record_carries_relabelling_and_the_isotope_rule(tmp_path):
    """Review 6 S10 and S11 in espmap.ensemble/1: which members were relabelled
    to member 1's diastereotopic labels, and that isotope labels must match."""
    plan = em.plan_ensemble(path=Path(__file__).parent / "data" / "ensemble" /
                            "E7_poses_one_ligand.sdf", ensemble="supplied")
    plan = dict(plan, members=plan["members"][:4])
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(), log=lambda s: None)
    assert [bool(m["relabelled"]) for m in rec["members"]] == [False, True, False, True]
    assert rec["identity"]["isotopes"] == "must_match" and rec["identity"]["isotope_labels"] == []


# --- convergence: not judgeable is not failed (the v0.3.0 acceptance run) -------------
def _conv_rec(plan, n, tmp_path, site=None):
    p = dict(plan, members=plan["members"][:n])
    return em.run_ensemble(p, tmp_path, runner=Stub(site=site), log=lambda s: None)


@pytest.mark.parametrize("n", [1, 2])
def test_too_few_members_is_not_judgeable_not_failed(plan, tmp_path, n):
    rec = _conv_rec(plan, n, tmp_path)
    assert rec["convergence"]["converged"] is None
    assert rec["validity"]["convergence"] is None
    assert rec["validity"]["members"] is True


def test_three_or_more_members_that_converge_are_true(plan, tmp_path):
    rec = _conv_rec(plan, 4, tmp_path, site=lambda n: 1.0 + 0.1 * n)
    assert rec["convergence"]["converged"] is True and rec["validity"]["convergence"] is True


def test_three_or_more_members_that_do_not_converge_are_false(plan, tmp_path):
    rec = _conv_rec(plan, 4, tmp_path, site=lambda n: 9.0 if n == 3 else 0.0)
    assert rec["convergence"]["converged"] is False and rec["validity"]["convergence"] is False


def test_not_judgeable_convergence_is_null_in_validity_and_allowed_by_the_schema(plan, tmp_path):
    """The invariant: convergence.converged is null => validity.convergence is null."""
    types = _SCHEMA["properties"]["validity"]["properties"]["convergence"]["type"]
    assert "null" in (types if isinstance(types, list) else [types])
    for n in (1, 2, 4):
        rec = _conv_rec(plan, n, tmp_path / str(n))
        if rec["convergence"]["converged"] is None:
            assert rec["validity"]["convergence"] is None


@pytest.mark.parametrize("origin, moves", [("etkdg_mmff94s", False), ("etkdg_uff", False),
                                           ("mmff94s", False), (None, True),
                                           ("dft:b3lyp/def2tzvp", True), ("undeclared", True)])
def test_optimize_may_move_only_a_generated_geometry(origin, moves):
    """The guard protects a pose or a reference geometry from --optimize. A
    member file the ensemble wrote from its own embedder (etkdg_mmff94s) is
    neither, and --weight dft needs those members optimised."""
    gi = {"source": "file 3D coordinates used as-is", "declared_geometry": origin}
    assert em.optimize_would_move_supplied(gi) is moves



def test_optimised_members_are_complete_and_resumed(tmp_path):
    """Found by the v0.3.0 acceptance run: with --optimize every member ran
    (rc 0) and was then counted failed, because its record's geometry is the
    optimised one. The input is named by optimised_from_sha1."""
    plan = em.plan_ensemble(smiles="CCCCCCO", ensemble=3, weight="dft", optimize=True)
    rec = em.run_ensemble(plan, tmp_path, runner=Stub(optimize=True), log=lambda s: None)
    assert rec["excluded"] == [] and rec["k_obtained"] == 3
    assert all(m["optimized"] for m in rec["members"])
    again = Stub(optimize=True)
    em.run_ensemble(plan, tmp_path, runner=again, log=lambda s: None)
    assert again.calls == []                           # resumed, not recomputed



class BuriedOne(Stub):
    """Members with a second atom that has no surface under nearest-nucleus
    in every member (as most heavy atoms of a chain do): never judgeable."""

    def __call__(self, argv, cwd, env, timeout):
        rc = super().__call__(argv, cwd, env, timeout)
        js = Path(_arg(argv, "--json"))
        if js.exists():
            rec = json.loads(js.read_text())
            rec["sites"]["atoms"].append(dict(index=1, atom_map=0, element="C",
                                              nearest_nucleus=dict(V_area_mean_kcal_mol=None),
                                              radius_scaled=dict(V_area_mean_kcal_mol=3.0)))
            js.write_text(json.dumps(rec))
        return rc


def test_an_unjudgeable_site_does_not_hide_a_verdict(plan, tmp_path):
    """Found by the v0.3.0 acceptance run: E4 had 7 sites not converged and
    reported null, because atoms buried under nearest-nucleus are never
    judgeable. Overall: false if any judgeable site fails, else true if any
    is judgeable, else null; the unjudgeable sites stay listed."""
    p = dict(plan, members=plan["members"][:4])
    ok = em.run_ensemble(p, tmp_path / "ok", runner=BuriedOne(site=lambda n: 1.0 + 0.1 * n),
                         log=lambda s: None)
    assert ok["convergence"]["converged"] is True and ok["validity"]["convergence"] is True
    assert ok["convergence"]["not_judgeable"] == ["atom 1 map 0 C | nearest_nucleus"]
    bad = em.run_ensemble(p, tmp_path / "bad", runner=BuriedOne(site=lambda n: 9.0 if n == 3 else 0.0),
                          log=lambda s: None)
    assert bad["convergence"]["converged"] is False and bad["validity"]["convergence"] is False


# --- review 7 (the user's, 2026-10-24): counterexamples before any fix -----------------
class ProtocolStub(Stub):
    """A member run whose analysis_key follows the protocol flags it is given,
    as a real member record's does."""

    def __call__(self, argv, cwd, env, timeout):
        rc = super().__call__(argv, cwd, env, timeout)
        js = Path(_arg(argv, "--json"))
        if js.exists():
            rec = json.loads(js.read_text())
            proto = [a for a in argv if a not in ("--file", "--json", "--npz", "--out")
                     and not a.endswith((".sdf", ".json", ".npz", ".png"))]
            rec["analysis_key"] = em._sha1(" ".join(proto))
            js.write_text(json.dumps(rec))
        return rc


def _rerun_is_safe(plan, tmp_path, first, second):
    """A rerun under another protocol in the same directory must refuse, or
    recompute every member; never silently reuse what the first computed."""
    em.run_ensemble(plan, tmp_path, first, runner=ProtocolStub(), log=lambda s: None)
    again = ProtocolStub()
    try:
        rec = em.run_ensemble(plan, tmp_path, second, runner=again, log=lambda s: None)
    except em.EnsembleRefused:
        return
    assert len(again.calls) == len(plan["members"]), "members reused under another protocol"
    assert rec["analysis_key"] == em._sha1(" ".join(second))


def test_r7_e1_a_rerun_under_another_basis_does_not_reuse_members(plan, tmp_path):
    _rerun_is_safe(plan, tmp_path, ["--basis", "def2-svp"], ["--basis", "def2-tzvp"])


def test_r7_e2_a_rerun_in_solvent_does_not_reuse_gas_members(plan, tmp_path):
    _rerun_is_safe(plan, tmp_path, ["--basis", "def2-svp"], ["--basis", "def2-svp", "--solvent", "water"])


def test_r7_e3_a_failure_under_one_protocol_does_not_suppress_another(plan, tmp_path):
    em.run_ensemble(plan, tmp_path, ["--basis", "def2-svp"], runner=Stub(fail={"conf_01"}),
                    log=lambda s: None)
    again = ProtocolStub()
    try:
        rec = em.run_ensemble(plan, tmp_path, ["--basis", "def2-tzvp"], runner=again,
                              log=lambda s: None)
    except em.EnsembleRefused:
        return
    assert "conf_01" in [c for c, _ in again.calls]
    assert all(e["conformer"] != "conf_01" for e in rec["excluded"])


def test_r7_e4_two_molecules_have_different_ensemble_keys(tmp_path):
    keys = []
    for smi in ("CO", "CN"):                           # one selection outcome each: 1 member
        p = em.plan_ensemble(smiles=smi, ensemble=4)
        keys.append(em.run_ensemble(p, tmp_path / smi, runner=Stub(), log=lambda s: None)["ensemble_key"])
    assert keys[0] != keys[1]


def test_r7_e5_two_supplied_sets_of_one_molecule_have_different_keys(tmp_path):
    path = Path(__file__).parent / "data" / "ensemble" / "E7_poses_one_ligand.sdf"
    p = em.plan_ensemble(path=path, ensemble="supplied")
    keys = [em.run_ensemble(dict(p, members=p["members"][s]), tmp_path / str(i), runner=Stub(),
                            log=lambda s: None)["ensemble_key"]
            for i, s in enumerate((slice(0, 3), slice(3, 6)))]
    assert keys[0] != keys[1]


class OptimiseTo(Stub):
    """--optimize whose final geometry is the plan member `target[label]`'s:
    the record carries those coordinates, as a real optimised record does."""

    def __init__(self, plan, target, **kw):
        super().__init__(optimize=True, **kw)
        self.plan, self.target = plan, target

    def __call__(self, argv, cwd, env, timeout):
        rc = super().__call__(argv, cwd, env, timeout)
        label = Path(_arg(argv, "--file")).stem
        js = Path(_arg(argv, "--json"))
        if js.exists():
            rec = json.loads(js.read_text())
            src = self.plan["members"][self.target.get(label, int(label.split("_")[1]))]
            rows = [r.split() for r in src["spec"].split(";")]
            rec["geometry"].update(symbols=[r[0] for r in rows],
                                   xyz_angstrom=[[float(x) for x in r[1:]] for r in rows],
                                   sha1=em.spec_sha1(src["spec"]))
            js.write_text(json.dumps(rec))
            np.savez_compressed(_arg(argv, "--npz").replace(".part.npz", ".part.npz"),
                                geometry_sha1=rec["geometry"]["sha1"])
        return rc


@pytest.mark.xfail(strict=True, reason="R7-3: shape comes from the pre-optimisation plan")
def test_r7_e6_the_record_describes_the_optimised_geometry(tmp_path):
    p = em.plan_ensemble(smiles="CC(=O)OCC/C=C\CCCCCCCC/C=C\CCCC", ensemble=10, weight="dft", optimize=True)
    folded = [i for i, m in enumerate(p["members"]) if m["folded"]]
    extended = [i for i, m in enumerate(p["members"]) if not m["folded"]]
    if not folded or not extended:
        pytest.skip("this plan has no folded/extended pair")
    src, dst = extended[0], folded[0]
    rec = em.run_ensemble(p, tmp_path, runner=OptimiseTo(p, {em.member_label(src): dst}),
                          log=lambda s: None)
    m = next(x for x in rec["members"] + rec["excluded"] if x["conformer"] == em.member_label(src))
    assert m["folded"] is True                         # the final geometry is the folded one
    assert m["rg_angstrom"] == pytest.approx(p["members"][dst]["rg"])


@pytest.mark.xfail(strict=True, reason="R7-3: two starts at one minimum are both weighted")
def test_r7_e7_two_starts_at_one_minimum_are_not_double_counted(tmp_path):
    p = em.plan_ensemble(smiles="CCCCCCO", ensemble=4, weight="dft", optimize=True)
    rec = em.run_ensemble(p, tmp_path, runner=OptimiseTo(p, {"conf_01": 0}), log=lambda s: None)
    finals = [m["geometry_sha1"] for m in rec["members"]]
    assert len(finals) == len(set(finals)), "one minimum counted twice"
    dup = [e for e in rec["excluded"] if e["conformer"] == "conf_01"]
    assert dup and "minimum" in dup[0]["reason"]
    assert sum(m["weight"] for m in rec["members"]) == pytest.approx(1.0)


def test_r7_1_a_directory_without_a_manifest_is_not_resumed(plan, tmp_path):
    """Member runs from before the manifest existed: what computed them is
    unknown, so nothing is resumed."""
    em.run_ensemble(plan, tmp_path, ["--basis", "def2-svp"], runner=Stub(), log=lambda s: None)
    (tmp_path / em.RUN_MANIFEST).unlink()
    with pytest.raises(em.EnsembleRefused) as e:
        em.run_ensemble(plan, tmp_path, ["--basis", "def2-svp"], runner=Stub(), log=lambda s: None)
    assert e.value.findings[0][1] == "ENSEMBLE_RESUME_PROTOCOL_MISMATCH"


def test_r7_1_the_same_protocol_resumes(plan, tmp_path):
    em.run_ensemble(plan, tmp_path, ["--basis", "def2-svp"], runner=Stub(), log=lambda s: None)
    again = Stub()
    em.run_ensemble(plan, tmp_path, ["--basis", "def2-svp"], runner=again, log=lambda s: None)
    assert again.calls == []
