"""Independent finite-loop Kraus execution and local stopped-clock bounds."""
from fractions import Fraction as F
from copy import deepcopy
import hashlib

import numpy as np
import pytest

from . import local_preparation as p


def _native_apply(rho, seam, law, moving=None):
    """Direct primitive grouping, independent of production superoperators."""
    states = p.native.states(binary=True)
    positions = {state: i for i, state in enumerate(states)}
    reference = len(rho)//66
    source = rho.reshape(66, reference, 66, reference)
    answer = np.zeros_like(source)
    for coin in (0, 1):
        groups = {}
        for i, state in enumerate(states):
            out, (_, current) = p.native.native(state, seam, coin)
            if moving is not None and bool(current) != moving:
                continue
            label = current if law == "C" else state[seam[0]]-state[seam[1]]
            groups.setdefault(label, []).append((i, positions[out]))
        for mapping in groups.values():
            inputs, outputs = np.array(mapping).T
            answer[outputs[:, None], :, outputs[None, :], :] += source[inputs[:, None], :, inputs[None, :], :]/2
    return answer.reshape(rho.shape)


def _tick(rho, phase, law, *, deadline, rounds, herald=False):
    """Literal chronological wait: selector timeout, native repair, tick test."""
    seams = p.quantum_mixing.seams()
    if phase == "warmup":
        rejection = F(2, 32)**rounds
        weights = {edge: (1-rejection)/30 for edge in seams}
        cap = p.CORE
    else:
        cap = p.CAPS[0]
        weights = {}
        for observer in p.CAPS:
            available = [edge for edge in seams if not set(edge).intersection(observer)]
            # Source-independent binary rejection arithmetic.
            denominator = 1 << (len(available)-1).bit_length()
            rejection = F(denominator-len(available), denominator)**rounds
            for edge in available:
                weights[edge] = weights.get(edge, F())+(1-rejection)/(2*len(available))
        rejection = 1-sum(weights.values())
    stopped = [np.zeros_like(rho), np.zeros_like(rho)]
    timeout, live = np.zeros_like(rho), rho
    for _ in range(deadline):
        timeout += float(rejection)*live
        remote = np.zeros_like(rho)
        for edge, weight in weights.items():
            if set(edge).intersection(cap):
                if herald:
                    for bit in (0, 1):
                        stopped[bit] += float(weight)*_native_apply(live, edge, law, bool(bit))
                else:
                    stopped[1] += float(weight)*_native_apply(live, edge, law)
            else:
                remote += float(weight)*_native_apply(live, edge, law)
        live = remote
    return stopped[0], stopped[1], timeout+live


def _literal(rho, law, warmup_ticks=1, completion_ticks=1, deadline=2, rounds=1):
    output = {flag: np.zeros_like(rho) for flag in p.FLAGS}
    live = rho
    for _ in range(warmup_ticks):
        _, live, timeout = _tick(live, "warmup", law, deadline=deadline, rounds=rounds)
        output["warmup_timeout"] += timeout
    output["herald_failure"], live, timeout = _tick(live, "warmup", law, deadline=deadline,
                                                   rounds=rounds, herald=True)
    output["warmup_timeout"] += timeout
    for _ in range(completion_ticks):
        _, live, timeout = _tick(live, "completion", law, deadline=deadline, rounds=rounds)
        output["completion_timeout"] += timeout
    output["success"] = live
    return output


def _entangled_input():
    basis = p.native.states(binary=True)
    i = next(i for i, x in enumerate(basis) if x[0] == x[2] == 1)
    j = next(i for i, x in enumerate(basis) if x[1] == x[2] == 1)
    vector = np.zeros(132, complex)
    vector[2*i], vector[2*j+1] = .6, .8j
    return np.outer(vector, vector.conj())


def test_source_reference_and_all_matrix_unit_premises():
    audit = p.exact_audit()
    assert audit["core_containing_cap_algebras"] == 1024
    assert audit["completion_tick_fraction"] == "5/24"
    assert audit["mixed_attempt_gap_lower"] == "1/600"
    assert all(row["signed_difference_blocks_preserved"] for row in audit["laws"].values())
    assert p.Counter(p.reference()) == {F(2, 225): 45, F(1, 36): 20, F(2, 45): 1}
    assert p.reference() != p.common_caps.reference()


def test_exact_budget_includes_outward_deadlines_selector_loss_and_normalization():
    b = p.budget()
    assert b["warmup_ticks"] == 496
    assert b["completion_ticks"] == 19200
    assert b["completion_incident_fraction"] == F(5, 24)
    assert b["private_per_tick_attempt_deadline_counter_required"]
    h = b["per_tick_attempt_deadline"]
    for phase, ticks in (("warmup", b["warmup_ticks"]+1), ("completion", b["completion_ticks"])):
        raw = ticks*(h*b[f"{phase}_selector_timeout_per_slot"]+(1-b[f"{phase}_incident_fraction"])**h)
        rounded = b[f"{phase}_cutoff_probability_upper"]
        assert raw <= rounded < raw+F(1, 2**80)
    assert b["complete_flagged_reference_halftrace_upper"] < F(6, 10**8)
    assert b["conditional_reference_uniform_halftrace_upper"] < F(8, 10**7)
    lo, hi = b["success_probability_interval"]
    assert lo < F(5, 33) < hi
    assert b["operator_sandwich_radius"] == b["complete_flagged_reference_halftrace_upper"]/(F(5, 33)*min(p.reference()))
    assert b["total_native_attempts_upper"] == b["work_native_attempts_upper"]+b["auxiliary_native_attempts_upper"]
    assert b["total_native_attempts_upper"] == 659460352
    assert b["expected_work_native_attempts_upper"] < 15621
    assert b["expected_auxiliary_native_attempts_upper"] < 90000
    assert b["expected_total_native_attempts_upper"] == b["expected_work_native_attempts_upper"]+b["expected_auxiliary_native_attempts_upper"]
    assert b["controller_register_bits_subtotal"] == sum(b["controller_register_bits"].values()) == 51


@pytest.mark.parametrize("law", ("C", "E"))
def test_literal_complete_reference_history_includes_coherent_failure_channels(law):
    rho = _entangled_input()
    actual = p.execute(rho, law=law, warmup_ticks=1, completion_ticks=1, deadline=2, selector_rounds=1)
    expected = _literal(rho, law)
    for flag in p.FLAGS:
        np.testing.assert_allclose(actual[flag], expected[flag], rtol=0, atol=2e-15)
        assert np.linalg.eigvalsh(actual[flag]).min() > -2e-15
        assert np.trace(actual[flag]).real > 0
    reference = lambda a: np.trace(a.reshape(66, 2, 66, 2), axis1=0, axis2=2)
    np.testing.assert_allclose(reference(sum(actual.values())), reference(rho), rtol=0, atol=2e-14)
    # Timeouts retain a coherent quantum state; replacing them by scalar mass
    # or a fresh preparation would fail this independently enumerated equality.
    failure = actual["warmup_timeout"]
    assert np.linalg.norm(failure-np.diag(np.diag(failure))) > .05


def test_c_and_e_complete_failure_maps_are_not_interchangeable():
    rho = _entangled_input()
    c = p.execute(rho, law="C", warmup_ticks=1, completion_ticks=1, deadline=2, selector_rounds=1)
    e = p.execute(rho, law="E", warmup_ticks=1, completion_ticks=1, deadline=2, selector_rounds=1)
    assert np.linalg.norm(c["warmup_timeout"]-e["warmup_timeout"]) > 1e-4


def test_mixed_selector_acceptance_bias_and_finite_deadline_are_real():
    selected = p.selectors(1)
    ideal, _, _ = p._schedule()
    acceptance = 1-selected["completion_timeout_probability"]
    weights = selected["completion_accepted_seam_weights"]
    assert any(weights[e] != acceptance*ideal[e] for e in ideal)
    rho = _entangled_input()
    one = p.execute(rho, warmup_ticks=0, completion_ticks=0, deadline=1, selector_rounds=1)
    assert np.trace(one["warmup_timeout"]).real > .7
    assert abs(sum(np.trace(x).real for x in one.values())-1) < 1e-14


def test_population_execution_all_inputs_and_arbitrary_reference_certificate():
    b = p.budget()
    c, e = p.population_maps(), p.population_maps(law="E")
    for flag in p.FLAGS:
        np.testing.assert_array_equal(c[flag], e[flag])
    np.testing.assert_allclose(sum(c.values()).sum(axis=0), 1, rtol=0, atol=3e-12)
    rho = _entangled_input()
    ref = np.trace(rho.reshape(66, 2, 66, 2), axis1=0, axis2=2)
    out = p.execute(rho, warmup_ticks=b["warmup_ticks"], completion_ticks=b["completion_ticks"],
                    deadline=b["per_tick_attempt_deadline"], selector_rounds=b["selector_rounds"])
    q = float(F(5, 33))
    raw = np.array(list(map(float, p.incident.reference(p.CORE)[0])))
    ideal = {"success": np.kron(q*np.diag(list(map(float, p.reference()))), ref),
             "herald_failure": np.kron(np.diag(1/66-q*raw), ref),
             "warmup_timeout": np.zeros_like(rho), "completion_timeout": np.zeros_like(rho)}
    distance = sum(np.linalg.norm(out[flag]-ideal[flag], "nuc") for flag in p.FLAGS)/2
    assert distance < float(b["complete_flagged_reference_halftrace_upper"])


def test_stopped_completion_similarity_not_self_adjoint_or_global_clock():
    weights, fraction, gamma = p._schedule()
    a, r = np.zeros((66, 66)), np.zeros((66, 66))
    basis = p.native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    for seam, weight in weights.items():
        target = a if set(seam).intersection(p.CAPS[0]) else r
        for j, x in enumerate(basis):
            for coin in (0, 1):
                y, _ = p.native.native(x, seam, coin)
                target[positions[y], j] += float(weight)/2
    m = np.eye(66)-r
    stopped = np.linalg.solve(m.T, a.T).T
    assert np.linalg.norm(stopped-stopped.T) > .01
    counts, scale = p.cut.completion_counts(p.CORE)
    indices = np.arange(66)*67
    projection = counts[indices, :][:, indices].toarray()/scale
    vals, vecs = np.linalg.eigh(m)
    inverse_root = (vecs/np.sqrt(vals))@vecs.T
    b = inverse_root@a@inverse_root
    assert np.linalg.eigvalsh(b).min() > -2e-14
    assert np.linalg.eigvalsh(b-projection).max() <= 1-float(gamma)+2e-14
    for n in (1, 600):
        bound = float(fraction)**-.5*(1-float(gamma))**n
        assert np.linalg.norm(np.linalg.matrix_power(stopped, n)-projection, 2) < bound
    assert all(not set(edge).intersection(p.CORE) for edge in weights)


@pytest.mark.parametrize("kwargs", ({"precision_bits": True}, {"precision_bits": 11},
    {"deadline": True}, {"deadline": 0}, {"deadline": 1}, {"selector_rounds": 0},
    {"selector_rounds": 33}, {"selector_rounds": 1}))
def test_invalid_or_uncertifiable_budgets_fail_closed(kwargs):
    with pytest.raises(ValueError):
        p.budget(**kwargs)


@pytest.mark.parametrize("kind", ("indefinite", "distributed_negative_mass", "aggregate_hermitian_defect",
    "combined_negative_and_trace_error", "inflated", "nonfinite", "masked", "boolean", "object"))
def test_hostile_density_input_is_rejected(kind):
    rho = np.zeros((66, 66), complex)
    rho[0, 0] = 1
    if kind == "indefinite":
        rho[0, 0], rho[1, 1] = 2, -1
    elif kind == "distributed_negative_mass":
        rho -= np.eye(66)*1e-13
        rho[0, 0] += 66e-13
    elif kind == "aggregate_hermitian_defect":
        rho = np.eye(66, dtype=complex)/66+np.triu(np.ones((66, 66)), 1)*2e-14
    elif kind == "combined_negative_and_trace_error":
        rho[0, 0], rho[1, 1] = 1+1.2e-13, -6e-14
    elif kind == "inflated":
        rho *= 2
    elif kind == "nonfinite":
        rho[0, 1] = np.inf
    elif kind == "masked":
        rho = np.ma.array(rho)
    elif kind == "boolean":
        rho = rho.real.astype(bool)
    elif kind == "object":
        rho = rho.astype(object)
    with pytest.raises(ValueError):
        p.execute(rho)


def test_zero_input_is_preserved_without_normalization():
    assert all(not np.any(x) for x in p.execute(np.zeros((66, 66), complex)).values())


def test_numerical_density_screen_does_not_symmetrize_or_mutate_the_input():
    rho = np.eye(66, dtype=complex)/66
    rho[0, 1] = 1e-15
    saved = rho.copy()
    actual = p.execute(rho, warmup_ticks=0, completion_ticks=0, deadline=1, selector_rounds=1)
    expected = _literal(rho, "C", warmup_ticks=0, completion_ticks=0, deadline=1, rounds=1)
    np.testing.assert_array_equal(rho, saved)
    for flag in p.FLAGS:
        np.testing.assert_allclose(actual[flag], expected[flag], rtol=0, atol=2e-17)
        assert abs(actual[flag][0, 1]-expected[flag][0, 1]) < 1e-29
        assert actual[flag][1, 0] == expected[flag][1, 0] == 0
    assert np.linalg.norm(actual["warmup_timeout"]-actual["warmup_timeout"].conj().T) > 1e-16


def test_nonfinite_computed_eigenvalues_fail_closed(monkeypatch):
    monkeypatch.setattr(np.linalg, "eigvalsh", lambda a: np.full(len(a), np.nan))
    with pytest.raises(ValueError):
        p.execute(np.eye(66)/66)


def test_nonfinite_power_and_channel_intermediates_fail_closed(monkeypatch):
    with pytest.raises(ArithmeticError):
        p._power_sum(np.array([[np.finfo(float).max]]), 2)

    def forged(law, indices, **kwargs):
        return {flag: np.full((len(indices), len(indices)), np.inf) for flag in p.FLAGS}

    monkeypatch.setattr(p, "_block_maps", forged)
    with pytest.raises(ArithmeticError):
        p.execute(np.eye(66)/66)


@pytest.fixture(scope="module")
def receipt():
    result = p.run()
    p.compare_receipt(result, result)
    return result


def test_receipt_binds_regional_basis_and_complete_palm_anchor_dependencies(receipt):
    # cap_schedule and basis_supports call regional.source, which checks the
    # Palm anchor via event and relational native instruments. These are
    # mathematical inputs even though no global-Palm preparation is run.
    required = {
        "model/native_modular_reference/common_caps.py",
        "model/native_modular_reference/regional.py",
        "model/native_modular_reference/preparation.py",
        "model/event_circulation/construction.py",
        "model/native_circulation/construction.py",
        "model/native_relational_readout/construction.py",
        "model/native_quantum_history/construction.py",
        "model/observer/quantum_memory.py",
        "model/primitive_observer/transfer.py",
        "model/primitive_observer/runtime.py",
        "model/native_clock/clock.py",
        "model/observer/native.py",
    }
    assert required <= set(receipt["source_sha256"])
    for path in required:
        assert receipt["source_sha256"][path] == hashlib.sha256((p.ROOT/path).read_bytes()).hexdigest()


@pytest.mark.parametrize("field", ("maximum_complete_halftrace_diagnostic",
    "maximum_total_probability_residual_diagnostic", "warmup_timeout", "completion_timeout"))
@pytest.mark.parametrize("actual_side", (False, True))
def test_tiny_negative_receipt_quantities_never_enter_replay_tolerance(receipt, field, actual_side):
    forged = deepcopy(receipt)
    row = forged["population_execution"]
    if field in p.FLAGS:
        row["basis_zero_output_flag_probabilities_diagnostic"][field] = -1e-12
    else:
        row[field] = -1e-12
    # This mutation lies inside the shared numerical tolerance; only the
    # local physical-domain screen excludes it, on either comparison side.
    p.common_caps.full.corner.compare_receipt(forged, receipt)
    with pytest.raises(ValueError):
        p.compare_receipt(receipt, forged) if actual_side else p.compare_receipt(forged, receipt)


def test_receipt_mass_must_match_its_recorded_maximum_residual(receipt):
    forged = deepcopy(receipt)
    forged["population_execution"]["basis_zero_output_flag_probabilities_diagnostic"]["herald_failure"] += 1e-10
    p.common_caps.full.corner.compare_receipt(forged, receipt)
    with pytest.raises(ValueError):
        p.compare_receipt(forged, receipt)


@pytest.mark.parametrize("attack", ("cancelling_negative_flag", "nan_flag", "infinite_flag"))
def test_producer_checks_each_original_population_flag(monkeypatch, attack):
    maps = p.population_maps()
    if attack == "cancelling_negative_flag":
        maps["warmup_timeout"][0, 0] = -1e-14
        maps["success"][0, 0] += 1e-14
    else:
        maps["warmup_timeout"][0, 0] = np.nan if attack == "nan_flag" else np.inf
    monkeypatch.setattr(p, "population_maps", lambda: maps)
    with pytest.raises(ArithmeticError):
        p.run()
