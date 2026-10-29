"""Independent endpoint oracle and optimal-control checks for moving memory."""
from collections import Counter
from fractions import Fraction as F

import pytest

from . import mobile as m
from . import screening
from ..observer.source_discrimination import circuit


def _oracle(center, law, seam, coin, raw):
    columns = []
    for column in m.supports(center):
        outputs = []
        for pair in column:
            x = [int(i in pair) for i in range(12)]
            a, b = seam
            parity = (x[a]+x[b]) % 2
            low = (x[a]+x[b])//2
            left, right = low+coin*parity, low+(1-coin)*parity
            current = x[a]-left
            code = {0: 0, 1: 1, -1: 2}[current]+(4*parity if law == "E" else 0)
            if code == raw:
                x[a], x[b] = left, right
                outputs.append(circuit.encode(tuple(x)))
        columns.append(tuple(outputs))
    return tuple(columns)


def test_complete_raw_native_oracle_and_record_only_target():
    for center in m.legal_cores():
        for law in ("C", "E"):
            for seam in m.SEAMS:
                total = F()
                for coin in (0, 1):
                    for raw in range(1 << circuit.RECORD_BITS[law]):
                        images = _oracle(center, law, seam, coin, raw)
                        assert m.branch_images(center, law, seam, coin, raw) == images
                        probability = F(len(images[0]), 2*len(m.supports(center)[0]))
                        row = m.transition(center, law, seam, coin, raw)
                        assert row["weight"] == probability
                        total += probability
                        if probability:
                            current = {0: 0, 1: 1, 2: -1}[raw % 4]
                            a, b = seam
                            moved = tuple(sorted(b if p == a else a if p == b else p for p in center)) if current else center
                            assert row["native_core"] == moved
                            expected = seam if not screening.safe_exterior(moved) else moved
                            assert row["target_core"] == expected
                            assert row["reset_performed"] == (expected != moved)
                assert total == 1


def test_actual_rotation_equivariance_and_no_extra_frame_register():
    assert m.geometry_certificate()["equivariance_checks"] == 3600
    assert len(m.legal_cores()) == 60
    for center in m.legal_cores():
        base = m.prototype(center)
        charts = [p for p in m.rotations() if tuple(sorted(p[i] for i in base)) == center]
        assert len(charts) == 2
        for p in charts:
            expected = tuple(tuple(sorted(tuple(sorted(p[i] for i in pair)) for pair in column))
                             for column in screening.sharp_code(base))
            assert m.supports(center) == expected


def test_arbitrary_complex_qubit_core_occupation_lemma():
    # A family of different admissible codes, including unequal exterior
    # weights and complex row frames. This checks the universal proof's
    # square-block premise beyond the two implemented uniform prototypes.
    import numpy as np
    safe = screening.safe_exterior(m.CORE)
    frames = (np.eye(2), np.array([[0, 1], [1, 0]]),
              np.array([[1, 1j], [1j, 1]])/np.sqrt(2),
              np.array([[1, -1j], [-1j, 1]])/np.sqrt(2))
    weights = (F(1, 6), F(1, 6), F(1, 3), F(1, 3))
    V = np.zeros((66, 2), complex)
    pairs = tuple(screening.CORES)
    for j, weight, frame in zip(safe, weights, frames):
        for endpoint, row in zip(m.CORE, frame):
            V[pairs.index(tuple(sorted((endpoint, j))))] = np.sqrt(float(weight))*row
    assert np.allclose(V.conj().T@V, np.eye(2))
    for endpoint in m.CORE:
        effect = np.diag([int(endpoint in pair) for pair in pairs])
        assert np.allclose(V.conj().T@effect@V, np.eye(2)/2)
    for j in safe:
        effect = np.diag([int(j in pair) for pair in pairs])
        assert np.allclose(V.conj().T@effect@V, float(weights[safe.index(j)])*np.eye(2))


def test_sharp_q1_alone_does_not_imply_universal_walk_or_exact_correction():
    # Unequal endpoint amplitudes preserve normalization, orthogonality and
    # sharp q=1 while breaking the exact all-seam scalar-effect premise.
    import numpy as np
    pairs = tuple(screening.CORES)
    epsilon = F(1, 10)
    V = np.zeros((66, 2), complex)
    for logical, column in enumerate(m.supports(m.CORE)):
        for pair in column:
            scale = 1+epsilon if m.CORE[0] in pair else 1-epsilon
            V[pairs.index(pair), logical] = np.sqrt(float(scale)/len(column))
    assert np.allclose(V.conj().T@V, np.eye(2))
    a = np.diag([int(m.CORE[0] in pair) for pair in pairs])
    assert np.allclose(V.conj().T@a@V, float((1+epsilon)/2)*np.eye(2))
    outside = screening.safe_exterior(m.CORE)[0]
    exterior_effect = np.diag([int(outside in pair) for pair in pairs])
    gram = V.conj().T@exterior_effect@V
    assert not np.allclose(gram, np.trace(gram)*np.eye(2)/2)
    # Moving out of core endpoint0 on seam(0,5), coin0, occurs with a
    # different probability. It cannot be assigned the exact1/120 walk.
    movement = np.diag([int(0 in pair and 5 not in pair) for pair in pairs])
    moving_gram = V.conj().T@movement@V/60
    assert np.allclose(moving_gram, float((1+epsilon)/120)*np.eye(2))
    assert not np.allclose(moving_gram, np.eye(2)/120)


def test_first_hit_survival_and_mean_from_full66_walk():
    legal = m.legal_cores()
    for start in (m.CORE, (0, 2)):
        measure = {start: F(1)}
        for n in range(12):
            assert sum(measure.values()) == m.survival(n, initial=start)
            next_measure = Counter()
            for center, weight in measure.items():
                for target, p in m.native_kernel(center).items():
                    if target in legal:
                        next_measure[target] += weight*p
            measure = next_measure
        h = m.POTENTIAL[m.core_class(start)]
        assert h == 1+sum(p*m.POTENTIAL[m.core_class(target)]
                          for target, p in m.native_kernel(start).items())
    assert m.survival(1) == 1
    assert m.survival(2) == F(1799, 1800)


def test_full60_reset_walk_stationarity_and_directed_currents():
    kernels = {c: m.continuing_kernel(c) for c in m.legal_cores()}
    pi = {c: F(1, 50) if c in m.SEAMS else F(1, 75) for c in kernels}
    for target in kernels:
        assert sum(pi[c]*kernels[c].get(target, F()) for c in kernels) == pi[target]
    one_way = [(c, d) for c in kernels for d, p in kernels[c].items()
               if p and not kernels[d].get(c, F())]
    assert len(one_way) == 60
    assert sum(pi[c]*sum(p for d, p in m.native_kernel(c).items()
                        if m.core_class(d) == "antipodal") for c in kernels) == F(1, 150)


def test_finite_reset_expectation_and_all_policy_potential_bound():
    for start in (m.CORE, (0, 2)):
        measure, count = {start: F(1)}, F()
        for n in range(35):
            assert count == m.expected_resets(n, initial=start)
            assert count >= m.reset_lower_bound(n, initial=start)
            next_measure = Counter()
            for center, weight in measure.items():
                count += weight*sum(p for target, p in m.native_kernel(center).items()
                                    if m.core_class(target) == "antipodal")
                for target, p in m.continuing_kernel(center).items():
                    next_measure[target] += weight*p
            measure = next_measure
    # Every possible reset target, including an earlier intervention at a
    # legal core, obeys the bound used for arbitrary randomized policies.
    for before in screening.CORES:
        for after in m.legal_cores():
            assert m.POTENTIAL[m.core_class(after)]-m.POTENTIAL[m.core_class(before)] <= 150


def test_finite_horizon_dynamic_program_all_reset_choices():
    # Bellman recursion permits a reset to either legal class after EVERY
    # observed native step, and separately tests an initial preemptive reset.
    # Compare exact finite-horizon optimal values, not a floating asymptote.
    values = (F(), F())
    native = (m.KILLED[0]+(F(),), m.KILLED[1]+(F(1, 60),))
    assert m.expected_resets(0) == m.expected_resets(0, initial=(0, 2)) == 0
    for horizon in range(1, 174):
        reset = 1+min(values)
        continuation = (min(values[0], reset), min(values[1], reset), reset)
        values = tuple(sum(p*v for p, v in zip(row, continuation)) for row in native)
        assert values[0] >= m.reset_lower_bound(horizon)
        assert values[1] >= m.reset_lower_bound(horizon, initial=(0, 2))
        assert values == (m.expected_resets(horizon), m.expected_resets(horizon, initial=(0, 2)))
        assert values[1]-values[0] == (1-F(11, 12)**horizon)/5
        assert values[0] <= 1+min(values) and values[1] <= 1+min(values)
    certificate = m.finite_horizon_certificate()
    assert certificate["fixed_horizon_global_minimum"]
    assert certificate["verified_identities"] == ["P*a=a", "P*b=(11/12)*b", "a+(1/12)*b=g"]


def test_replay_retains_all_flags_and_reconstructs_boundary_reset():
    witness = next(row for seam in m.SEAMS for row in m.branches((0, 2), "C", seam)
                   if row["positive"] and row["reset_performed"])
    event = (witness["law"], witness["seam"], witness["coin"], witness["raw"])
    replay = m.replay((event,), initial=(0, 2))
    assert replay["terminal_core"] == witness["seam"]
    assert replay["reset_count"] == 1
    assert replay["archive"] == (event,)
    assert replay["probability"] == witness["weight"]/30
    abort = ("C", m.SEAMS[0], 0, 3)
    result = m.replay((abort,))
    assert result["halted"] and result["probability"] == 0
    assert result["archive"] == (abort,)
    with pytest.raises(ValueError, match="halted"):
        m.replay((abort, event))
    zero = next(row for row in m.branches(m.CORE, "E", m.SEAMS[0])
                if row["admitted"] and not row["positive"])
    assert zero["identity_zero_branch"] and zero["target_core"] == m.CORE


@pytest.mark.parametrize("center,law,seam,coin,raw", [
    ([0, 1], "C", (0, 1), 0, 0), ((True, 1), "C", (0, 1), 0, 0),
    ((0, 3), "C", (0, 1), 0, 0), ((1, 0), "C", (0, 1), 0, 0),
    ((0, 1), "X", (0, 1), 0, 0), ((0, 1), "C", [0, 1], 0, 0),
    ((0, 1), "C", (0, 2), 0, 0), ((0, 1), "C", (0, 1), True, 0),
    ((0, 1), "C", (0, 1), 0, True), ((0, 1), "C", (0, 1), 0, 4),
])
def test_invalid_controller_inputs_fail_closed(center, law, seam, coin, raw):
    with pytest.raises(ValueError):
        m.transition(center, law, seam, coin, raw)


@pytest.mark.parametrize("n", [-1, True, 1.0, "1"])
def test_attempt_counts_are_exact_nonnegative_integers(n):
    for fn in (m.survival, m.expected_resets, m.reset_lower_bound):
        with pytest.raises(ValueError):
            fn(n)


def test_complete_certificate_boundaries():
    receipt = m.exact_certificate()
    assert sum(row["raw"] for row in receipt["raw_census"].values()) == 43200
    assert receipt["stationary_reset_probability_per_attempt"] == "1/150"
    assert receipt["full60_detailed_balance"] is False
    assert not any(receipt["physical_scope"].values())
