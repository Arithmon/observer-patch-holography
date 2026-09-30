"""Adversarial controls for the same-metric stationary field history."""
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_whitney_maxwell_dynamics as verifier
needs_lake = pytest.mark.skipif(verifier.lake_path() is None, reason="Lean lake executable is unavailable")


def test_independent_stationary_volume_replay():
    result = verifier.verify(verifier.load())
    assert result["field_variations_per_history"] == 68
    assert result["events_per_history"] == 805
    assert result["exact_global_stiffness_bound"] == 24


@needs_lake
def test_fresh_producer_passes_independent_quadrature_and_event_replay():
    import whitney_maxwell_dynamics as producer
    fresh = producer.build()
    result = verifier.verify(fresh)
    assert result["instrumented_slices_per_history"] == 3
    assert fresh["stability_certificate"] == verifier.load()["stability_certificate"]


@needs_lake
def test_exact_pipeline():
    current = verifier.generate_stability_certificate()
    q = lambda a: {"a": str(a), "b": "0"}
    matrix = [[q(x) for x in row] for row in ((2, 1, 0), (1, 3, 1), (0, 1, 2))]
    problem, *_ = verifier.exact_source_problem()
    unseen = deepcopy(problem)
    unseen["sourceSha256"] = "unseen"
    unseen["targets"].append({"name": "unseen_spd_3", "matrix": matrix,
                              "certificate": None})
    result, witness = verifier.run_lean_certificate_problem(unseen)
    assert result["sourceSha256"] == "unseen" and len(witness) == 64
    duplicate = verifier.canonical(unseen).replace(b'"sourceSha256":"unseen"',
        b'"sourceSha256":"forged","sourceSha256":"unseen"')
    with pytest.raises(ValueError, match="Lean certificate producer failed"):
        verifier.run_lean_certificate_problem(duplicate)
    target = next(row for row in problem["targets"] if row["name"] == "stability_24")
    forged = deepcopy(current["certificates"]["stability_24"])
    forged["lower"][0][1] = {"a": "1", "b": "0"}
    target["certificate"] = forged
    with pytest.raises(ValueError, match="Lean certificate producer failed"):
        verifier.run_lean_certificate_problem(problem)
    for field in ("transportedGradient", "transportedCurl"):
        mutant = deepcopy(unseen)
        mutant["assembly"][field][0][0]["a"] = "2"
        with pytest.raises(ValueError, match="Lean certificate producer failed"):
            verifier.run_lean_certificate_problem(mutant)
    mutant = deepcopy(unseen)
    mutant["assembly"]["edgeReindex"][0]["sign"] *= -1
    with pytest.raises(ValueError, match="Lean certificate producer failed"):
        verifier.run_lean_certificate_problem(mutant)


def replace(path, value):
    def mutate(packet):
        row = packet
        for key in path[:-1]:
            row = row[key]
        row[path[-1]] = value
    return mutate


def event_mutation(op, mutate):
    def apply(packet):
        event = next(e for e in packet["executions"][0]["events"] if e["op"] == op)
        mutate(event)
    return apply


@pytest.mark.parametrize("mutation", [
    replace(["scope"], "PHYSICAL_CONTINUUM_CONFIRMED"),
    lambda p: p["pins"].pop("Lean/Screen/WhitneyMaxwellDynamics.lean"),
    replace(["pins", "code/electromagnetism/verify_cone_whitney_bridge.py"], "0"*64),
    replace(["numeric_policy", "atol"], 1.0),
    replace(["stability_certificate", "kernel_witness_sha256"], "0"*64),
    replace(["stability_certificate", "source_sha256"], "0"*64),
    lambda p: p["stability_certificate"]["source_manifest"].popitem(),
    replace(["dynamics", "source"], "dynamical charged matter from the OPH source"),
    replace(["dynamics", "energy"], "raw velocity energy in arbitrary gauge"),
    replace(["executions", 0, "instrumented_slices"], 65),
    replace(["executions", 0, "writable_slots"], 42),
    replace(["executions", 0, "gauge"], 0),
    replace(["executions", 0, "decoded", 2, 13], "100"),
    replace(["executions", 0, "decode_event_ids", 0], 0),
    replace(["executions", 0, "projection", "projected_E0", 0], "0"),
    lambda p: p["executions"][0]["metrics"]["ampere_hat"].__delitem__(slice(0, 12)),
    event_mutation("inputs", lambda e: e["writes"].__setitem__("h", "1")),
    event_mutation("gauss_project", lambda e: e["reads"].__delitem__("rho/0/0")),
    event_mutation("advance", lambda e: e["reads"]["d/1/13"].__setitem__("writer", 0)),
    event_mutation("advance", lambda e: e["reads"]["d/1/13"].__setitem__("value", "100")),
    event_mutation("advance", lambda e: e["parents"].append(999)),
    event_mutation("advance", lambda e: e["writes"].__setitem__("x/13", "0")),
    event_mutation("probe", lambda e: e["writes"].__setitem__("x/0", "100")),
    event_mutation("feedback", lambda e: e["reads"].__delitem__("clock")),
    event_mutation("decode", lambda e: e["writes"].__setitem__("d/0/13", "100")),
    event_mutation("public", lambda e: e["writes"].__setitem__("public/E/0/0", "100")),
    replace(["continuation", "scope"], "all time steps authenticated by serial instrument"),
    replace(["continuation", "A", 20, 0], 100.0),
    replace(["continuation", "rho_load", 20, 0], 1.0),
    replace(["continuation", "J_load", 20, 0], 1.0),
    replace(["continuation", "energy", 20], 100.0),
    replace(["continuation", "source_work", 0], 0.0),
    replace(["scalar_controls", "histories", 1, "A", 16], "0"),
    lambda p: p["nonclaims"].pop(),
])
def test_false_green_mutations_fail(mutation):
    original = verifier.load()
    packet = deepcopy(original)
    mutation(packet)
    assert json.dumps(packet, sort_keys=True) != json.dumps(original, sort_keys=True)
    with pytest.raises(ValueError):
        verifier.verify(packet)


def test_invisible_nonbinary_solver_write_rejected():
    packet = verifier.load()
    ex = packet["executions"][0]
    event = next(e for e in ex["events"] if e["op"] == "gauss_project")
    old = Q(event["writes"]["projection_z/0"])
    forged = old+Q(1, 10**1000)
    assert float(old) == float(forged) and forged != old
    event["writes"]["projection_z/0"] = str(forged)
    ex["projection"]["z"][0] = str(forged)
    with pytest.raises(ValueError, match="exact float64 output encoding"):
        verifier.verify(packet)


def test_fresh_transitive_custody_after_prior_success(monkeypatch):
    packet = verifier.load()
    verifier.verify(packet)
    original = Path.read_bytes
    upstream = verifier.ROOT / "Lean/Screen/SerialMaxwellReadout.lean"
    def tampered(path):
        content = original(path)
        return content+b"\n-- altered transitive provider\n" if path == upstream else content
    monkeypatch.setattr(Path, "read_bytes", tampered)
    # Neither the new receipt nor the unchanged immediate parent's hash is
    # enough to authenticate all of the immediate parent's current providers.
    with pytest.raises(ValueError, match="source pin"):
        verifier.verify(packet)


def test_parent_replay_allows_platform_roundoff_after_full_verification(monkeypatch):
    original = verifier.geometry.verify
    calls = []
    def rounded(parent):
        result = original(parent)
        calls.append(True)
        for key in ("volume_action", "gauss_max_abs", "ampere_hat_max_abs"):
            result[key] = float(np.nextafter(result[key], np.inf))
        return result
    monkeypatch.setattr(verifier.geometry, "verify", rounded)
    assert verifier.verify(verifier.load())["field_variations_per_history"] == 68
    assert calls == [True]


@pytest.mark.parametrize("mutation", [
    lambda r: r.__setitem__("tetrahedra", 20.0),
    lambda r: r.__setitem__("gauge_histories", True),
    lambda r: r.__setitem__("field_variations_per_history", 67),
    lambda r: r.__setitem__("volume_action", r["volume_action"] + 1e-5),
    lambda r: r.__setitem__("gauss_max_abs", True),
    lambda r: r.__setitem__("ampere_hat_max_abs", float("nan")),
    lambda r: r.__setitem__("unverified", 1),
    lambda r: r.pop("volume_action"),
])
def test_parent_summary_keeps_exact_census_and_numeric_error_guards(mutation):
    fresh = verifier.load()["parent_replay"]
    recorded = deepcopy(fresh)
    mutation(recorded)
    with pytest.raises(ValueError):
        verifier.verify_parent_summary(recorded, fresh)


def test_exact_variational_map_and_modified_hamiltonian():
    q, next_q, p, h, lam = sp.symbols("q next_q p h lam", positive=True)
    discrete_l = (next_q-q)**2/(2*h)-h*lam*(q*q+q*next_q+next_q*next_q)/6
    solution = sp.solve(sp.Eq(p, -sp.diff(discrete_l, q)), next_q)[0]
    next_p = sp.diff(discrete_l, next_q).subs(next_q, solution)
    f = 1+h*h*lam/6
    a, b, cc = (1-h*h*lam/3)/f, h/f, -h*lam*(1-h*h*lam/12)/f
    assert sp.simplify(solution-a*q-b*p) == 0
    assert sp.simplify(next_p-cc*q-a*p) == 0
    assert sp.simplify(a*a-b*cc) == 1
    # Consistent canonical initial data give second-order time convergence.
    errors = []
    exact = np.array([[np.cos(1), np.sin(1)], [-np.sin(1), np.cos(1)]])
    for count in (16, 32, 64):
        step = 1/count
        z = step*step
        factor = 1+z/6
        matrix = np.array([[1-z/3, step], [-step*(1-z/12), 1-z/3]])/factor
        errors.append(np.linalg.norm(np.linalg.matrix_power(matrix, count)-exact))
    assert 3.99 < errors[0]/errors[1] < 4.01
    assert 3.99 < errors[1]/errors[2] < 4.01


def test_stable_initial_value_problem_can_have_singular_endpoint_problem():
    h, lam = Q(1, 2), Q(12)  # z=3, strictly inside the IVP stability window.
    f, middle = 1+h*h*lam/6, 2-2*h*h*lam/3
    assert f > 0 and middle == 0 and h*h*lam < 12
    # With zero endpoints, every midpoint solves J=0; none solves J=1.
    for q1 in (Q(-3), Q(0), Q(7, 2)):
        assert f*0-middle*q1+f*0 == h*h*0
        assert f*0-middle*q1+f*0 != h*h*1


def test_explicit_utf8_file_reads(monkeypatch):
    original = Path.read_text
    def windows_default(path, *args, **kwargs):
        if not args and "encoding" not in kwargs:
            kwargs["encoding"] = "cp1252"
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", windows_default)
    assert verifier.verify(verifier.load())["field_variations_per_history"] == 68

@pytest.mark.parametrize("mutation", [
    lambda p: p.update(certificates={}),
    lambda p: p["certificates"]["stability_24"]["diagonal"][0].update(a="-999", b="0"),
])
def test_lean_free_certificate_payload_mutation_rejected(monkeypatch, mutation):
    monkeypatch.setattr(verifier, "lake_path", lambda: None)
    packet = deepcopy(verifier.load()["stability_certificate"])
    mutation(packet)
    with pytest.raises(ValueError, match="payload binding"):
        verifier.certify_stability(packet)


def test_exact_forms_match_independent_quadrature(monkeypatch):
    monkeypatch.setattr(verifier, "lake_path", lambda: None)
    targets = verifier.certify_stability(verifier.load()["stability_certificate"])
    vertices, boundary_edges, boundary_faces = verifier.geometry.source_mesh()
    vertices = [(0, 0, 0)] + vertices
    edges = [(0, u+1) for u in range(12)] + [(u+1, v+1) for u, v in boundary_edges]
    faces = [tuple(u+1 for u in face) for face in boundary_faces] + [(0, u+1, v+1) for u, v in boundary_edges]
    tets = [(0, *(u+1 for u in face)) for face in boundary_faces]
    *_, mass, face_mass = verifier.geometry.quadrature(vertices, edges, faces, tets)
    curl = verifier.geometry.coboundary(edges, faces)
    stiffness = curl.T @ face_mass @ curl
    verifier.check_exact_quadrature(targets, vertices, edges, faces, tets, mass, face_mass, stiffness)
    with pytest.raises(ValueError, match="exact-to-quadrature mass"):
        verifier.check_exact_quadrature(targets, vertices, edges, faces, tets, 2*mass, face_mass, stiffness)
    with pytest.raises(ValueError, match="exact-to-quadrature"):
        verifier.check_exact_quadrature(targets, vertices, edges, faces, tets[1:], mass, face_mass, stiffness)


@needs_lake
def test_payload_rewrite_fails_kernel_even_with_rehashed_receipt(tmp_path):
    import ast
    import re
    witness = (verifier.LEAN_ROOT / "Screen/WhitneyGeneratedCertificate.lean").read_text()
    match = re.search(r"^def edge_massData : .+ := (.+)$", witness, re.M)
    payload = ast.literal_eval(match[1].replace("#[", "["))
    payload[2][0][0] = (-999, 1, 0, 1)
    forged = repr(payload).replace("[", "#[")
    witness = witness[:match.start(1)] + forged + witness[match.end(1):]
    path = tmp_path / "forged.lean"
    path.write_text(witness)
    command, environment = verifier.lean_environment()
    result = verifier.native_run([*command, str(path)], environment)
    assert result.returncode != 0 and "decide" in result.stdout

@pytest.mark.parametrize("relative", ["Lean/Screen/WhitneyCurlRow3.lean", "Lean/Screen/LocalFaceMaxwellAction.lean"])
def test_transitive_consumer_source_mutation_rejected_without_lake(monkeypatch, relative):
    read = Path.read_bytes
    changed = verifier.ROOT / relative
    monkeypatch.setattr(verifier, "lake_path", lambda: None)
    monkeypatch.setattr(Path, "read_bytes", lambda path: read(path) + (b"\n-- mutation\n" if path == changed else b""))
    with pytest.raises(ValueError, match="stale certificate source manifest"):
        verifier.certify_stability(verifier.load()["stability_certificate"])

@pytest.mark.parametrize("raw", ['{"x":1,"x":2}', '{"x":1,"x":1}', '{"x":NaN}', '{"x":Infinity}'])
def test_ambiguous_json_rejected(tmp_path, raw):
    path = tmp_path / "bad.json"
    path.write_text(raw, encoding="utf-8")
    with pytest.raises(ValueError):
        verifier.load(path)
