"""Independent source/metric/causal replay and fail-closed receipt controls."""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path

import pytest


HERE = Path(__file__).resolve().parent


def module(name):
    spec = importlib.util.spec_from_file_location(name, HERE/(name+".py"))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


producer = module("source_net_causet")
verifier = module("verify_source_net_causet")


@pytest.fixture(scope="module")
def small():
    return producer.build_level(4)


def test_entire_frozen_packet_replayed_independently():
    packet = verifier.load()
    result = verifier.verify(packet)
    assert result["accepted"] is True
    assert result["native_physical_spacetime_selected"] is False
    assert [r["source_records"] for r in result["levels"]] == [27, 125, 512, 2197]
    assert [r["positive_inner_cone"] for r in result["levels"]] == [False, False, False, True]
    assert result["levels"][-1]["reads_per_full_trace"] == 1176764


@pytest.mark.parametrize("n", [True, False, 3, 8, 4.0, "4", None, -1])
def test_regulator_domain_rejects(n):
    with pytest.raises(ValueError):
        producer.fibonacci(n)
    with pytest.raises(ValueError):
        verifier.graph_and_source(n)


@pytest.mark.parametrize("a,b", [([0, 1], [0, 0, 0]), ([0, 0, 0], [False, 0, 0]),
                                  ([0, 0, 0.0], [0, 0, 0]), (None, [0, 0, 0]),
                                  ("123", [0, 0, 0])])
def test_source_control_domain(a, b):
    with pytest.raises(ValueError):
        producer.source_control(a, b)


@pytest.mark.parametrize("z", [[1, 0, 0, 0, 0, 0], [0]*5, [False]*6, [0.0]*6, None])
def test_source_section_domain(z):
    with pytest.raises(ValueError):
        producer.source_currents(z)


@pytest.mark.parametrize("x,expected", [((0, 0), 0), ((-1, 1), 1), ((2, -1), 1),
    ((-2, 1), -1), ((1, -1), -1), ((F(-1618034, 1000000), 1), -1),
    ((F(-1618033, 1000000), 1), 1)])
def test_independent_algebraic_signs(x, expected):
    assert producer.sign(x) == expected
    assert verifier.sgn((x[0]+F(x[1], 2), F(x[1], 2))) == expected


@pytest.mark.parametrize("n", [4, 5])
def test_all_edges_and_actual_source_records_match_independent_incidence(n):
    q, p, values, sites, _, graph = producer.geometry(n)
    qq, pp, _, _, _, other, records, _ = verifier.graph_and_source(n)
    assert (q, p) == (qq, pp)
    assert graph == other
    assert all(i in row for i, row in enumerate(graph))
    for b, z, currents in records:
        assert tuple(z) == producer.source_control([values[i][0] for i in b], b)
        assert tuple(currents) == producer.source_currents(z)
    # At least one pair distinguishes actual golden readback from the
    # auxiliary uniform-cell labels: the latter do not define these edges.
    uniform = [[j for j, t in enumerate(sites)
                if sum(((s[k]*p) % q-(t[k]*p) % q)**2 for k in range(3)) <= q]
               for s in sites]
    assert graph != uniform


def test_h_and_quadrature_H_are_separate_and_q13_inner_cone_nonvacuous():
    rows = verifier.load()["levels"]
    final = rows[-1]
    # These are analytic, exact Q(phi) inequalities, not floating geometry.
    for row in rows:
        h2 = tuple(F(v) for v in row["whole_cube_fill_h_squared_over_L2_Qphi"])
        q = row["q"]
        assert row["positive_inner_cone"] == (producer.sign(producer.sub((1, 0), producer.scale(4*q, h2))) > 0)
        assert F(row["quadrature_displacement_H_squared_over_L2"]) == F(12, q*q)
        # Replacing h by the assignment bound H would wrongly erase the
        # useful small-q inner cone (4 q H²/L² = 48/q > 1 here).
        assert 4*q*F(row["quadrature_displacement_H_squared_over_L2"]) > 1
    assert F(final["certified_inner_speed_lower"]) > F(3, 10)


def test_schedule_audit_and_intervention_are_distinct(small):
    forward, reverse = small["forward_execution"], small["within_layer_reversed_execution"]
    assert forward["layer_value_sha256"] == reverse["layer_value_sha256"]
    assert forward["audit_trace_sha256"] != reverse["audit_trace_sha256"]
    assert small["center_plus_one_intervention"]["layer_value_sha256"] != forward["layer_value_sha256"]
    assert small["intervention_response"][0]["support_count"] == 1
    assert small["intervention_response"][1]["support_count"] > 1
    assert small["intervention_response"][-1]["positive_integer_delta_sum"] > 1


def test_exhaustive_small_event_order_width_and_height(small):
    """Explicit all-pairs event order gives layer antichains and wait chains."""
    _, _, _, sites, _, graph = producer.geometry(4)
    hops = [verifier.breadth(graph, i) for i in range(len(sites))]
    layers = small["layer_steps"]
    events = [(j, i) for j in range(layers+1) for i in range(len(sites))]
    relation = {(u, v) for u in events for v in events
                if u[0] < v[0] and hops[u[1]][v[1]] <= v[0]-u[0]}
    assert all(((j, i), (j+1, i)) in relation for i in range(len(sites)) for j in range(layers))
    assert all(a[0] != b[0] for a, b in relation)
    assert small["exact_width"] == len(sites)
    assert small["exact_height_in_events"] == layers+1
    # Actual non-neighbors cannot signal in one layer merely because the
    # serialized hash audit has placed their writes in a definite order.
    assert any(((0, i), (1, j)) not in relation for i in range(len(sites)) for j in range(len(sites)))


MUTATIONS = [
    ("exact_width", 1), ("exact_height_in_events", 99), ("positive_inner_cone", True),
    ("layer_time_equals_radius", False), ("word_length_bound", 0),
    ("source_records_sha256", "0"*64), ("neighbors_including_wait_sha256", "0"*64),
    ("certified_inner_speed_lower", "1"), ("h_over_a_upper", "0"),
    ("h_over_a_upper", "2/3"), ("quadrature_displacement_H_squared_over_L2", "1/3"),
    ("q", True), ("record_count", 27.0), ("radius_squared_over_L2", "1/4"),
]


@pytest.mark.parametrize("key,value", MUTATIONS)
def test_geometric_and_scope_forgery_rejected(small, key, value):
    row = deepcopy(small)
    row[key] = value
    with pytest.raises(ValueError):
        verifier.verify_level(row)


@pytest.mark.parametrize("which", ["late_value", "audit_schedule", "erased_intervention", "wrong_writer_hash",
    "source_word", "missing_wait", "wrong_interval", "wrong_weight", "count_auxiliary_reads",
    "zero_cone_missing", "intervention_support", "intervention_amplitude", "hidden_new_field"])
def test_rehashed_scientific_forgery_rejected(small, which):
    row = deepcopy(small)
    if which == "late_value":
        row["forward_execution"]["layer_value_sums"][-1] += 1
    elif which == "audit_schedule":
        row["within_layer_reversed_execution"] = deepcopy(row["forward_execution"])
    elif which == "erased_intervention":
        row["center_plus_one_intervention"] = deepcopy(row["forward_execution"])
    elif which == "wrong_writer_hash":
        # Digesting a syntactically valid forged write does not authenticate
        # the writer against reconstructed register versions.
        row["forward_execution"]["audit_trace_sha256"] = verifier.hashed([[2, 0], [], [0, 3, [0, 0], 1]])
    elif which == "source_word":
        row["source_record_examples"][-1][-1][0] *= -1
    elif which == "missing_wait":
        graph = producer.geometry(4)[-1]
        graph = [[j for j in r if j != i] for i, r in enumerate(graph)]
        row["neighbors_including_wait_sha256"] = verifier.hashed(graph)
        row["forward_execution"] = verifier.execute(graph, row["layer_steps"])
    elif which == "wrong_interval":
        row["center_intervals"][-1]["inclusive_event_count"] += 1
    elif which == "wrong_weight":
        row["center_intervals"][-1]["count_over_L4_times_sqrt_q"] = "1"
    elif which == "count_auxiliary_reads":
        row["center_intervals"][0]["inclusive_event_count"] += row["forward_execution"]["authenticated_read_count"]
    elif which == "zero_cone_missing":
        row["reachability_probes"][-1]["missing_ids_sha256"] = verifier.hashed([])
        row["reachability_probes"][-1]["exact_finite_cone_missing_count"] = -1
    elif which == "intervention_support":
        row["intervention_response"][1]["support_ids_sha256"] = verifier.hashed([row["intervention_source_id"]])
    elif which == "intervention_amplitude":
        row["intervention_response"][-1]["positive_integer_delta_sum"] = 1
    else:
        row["physical_causal_order_selected"] = True
    with pytest.raises(ValueError):
        verifier.verify_level(row)


@pytest.mark.parametrize("data", [b'{"a":1,"a":2}', b'{"a":1.0}', b'{"a":NaN}',
    b'{"a":1e9999}', b'{"a":Infinity}', b'{"a":-Infinity}', b'\xff', b' '*2000001])
def test_strict_loader(tmp_path, data):
    path = tmp_path/"bad.json"
    path.write_bytes(data)
    with pytest.raises(ValueError):
        verifier.load(path)


@pytest.mark.parametrize("which", ["pin", "source_path", "scope", "clock", "extra", "level_count", "type"])
def test_receipt_preflight_rejects_before_expensive_replay(which, monkeypatch):
    packet = verifier.load()
    if which == "pin":
        packet["source_pins"][verifier.PIN_PATHS[0]] = "0"*64
    elif which == "source_path":
        value = packet["source_pins"].pop(verifier.PIN_PATHS[0])
        packet["source_pins"]["../closure/source_net.py"] = value
    elif which == "scope":
        packet["scope"]["native_repair_selected"] = True
    elif which == "clock":
        packet["scope"]["model_clock_supplied"] = 1
    elif which == "extra":
        packet["physical_fit"] = "PASS"
    elif which == "level_count":
        packet["levels"].pop()
    else:
        packet = []
    def forbidden_replay(row):
        pytest.fail("preflight should reject before finite replay")
    monkeypatch.setattr(verifier, "verify_level", forbidden_replay)
    with pytest.raises(ValueError):
        verifier.verify(packet)


def test_no_untrusted_producer_import_or_float_receipt_tokens():
    source = (HERE/"verify_source_net_causet.py").read_text(encoding="utf-8")
    assert "import source_net_causet" not in source
    assert "importlib" not in source
    packet = verifier.load()
    assert packet["scope"]["finite_energy_signal_or_quantum_field_established"] is False
    assert verifier.load() == json.loads(producer.OUTPUT.read_bytes())
