"""Native herald, full-flag error, and noncommuting modular controls."""
from collections import defaultdict
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import preparation as p


def test_global_herald_is_the_actual_native_cp_map_for_both_laws():
    degree, sigma = p.anchor()
    states = p.relational.quantum.STATES
    index = p.relational.quantum.INDEX
    seams = p.relational.quantum.mixing.seams()
    independently_prepared = [F(0)] * 66
    for i, state in enumerate(states):
        moves = 0
        for seam in seams:
            for coin in (0, 1):
                output, (_, current) = p.native.native(state, seam, coin)
                if current:
                    moves += 1
                    independently_prepared[index[output]] += F(1, 66*60)
        assert moves == degree[i]
    assert sum(independently_prepared) == F(5, 33)
    assert tuple(value / F(5, 33) for value in independently_prepared) == tuple(sigma)
    assert p.moving_columns("C") == p.moving_columns("E")
    assert min(sigma) == F(1, 75)
    # The no-current maps are different, though both complete the instrument.
    assert p.marked("C").no_current_counts != p.marked("E").no_current_counts
    for law in ("C", "E"):
        source = p.marked(law)
        for i, d in enumerate(degree):
            assert F(source.no_current_counts[i][i], 60) + F(d, 60) == 1


def test_exact_diagonal_trial_counts_native_events_and_all_failures():
    initial = tuple(int(i in (0, 2)) for i in range(12))
    law = {initial: F(1)}
    for _ in range(3):
        after = defaultdict(F)
        for state, mass in law.items():
            for seam in p.relational.quantum.mixing.seams():
                for coin in (0, 1):
                    output, _ = p.native.native(state, seam, coin)
                    after[output] += mass / 60
        law = after
    success, failure = defaultdict(F), defaultdict(F)
    for state, mass in law.items():
        for seam in p.relational.quantum.mixing.seams():
            for coin in (0, 1):
                output, (_, current) = p.native.native(state, seam, coin)
                (success if current else failure)[output] += mass / 60
    row = p.diagonal_trial(initial, 3)
    assert row["success"] == tuple(success[state] for state in p.relational.quantum.STATES)
    assert row["failure"] == tuple(failure[state] for state in p.relational.quantum.STATES)
    assert sum(row["success"]) + sum(row["failure"]) == 1
    assert row["flagged_error"] <= row["warmup_error"]
    assert row["conditional_error"] <= row["warmup_error"] / p.Q


def test_flagged_contractivity_pays_for_conditioning_and_probability_change():
    # A two-outcome successful work state saturates the factor 1/q exactly.
    q, eps = F(5, 33), F(1, 1000)
    ideal = (q/2, q/2, 1-q)
    actual = (q/2+eps, q/2-eps, 1-q)
    flagged = sum(abs(x-y) for x, y in zip(actual, ideal, strict=True))/2
    conditional = (abs(actual[0]/q-F(1, 2))+abs(actual[1]/q-F(1, 2)))/2
    assert flagged == eps and conditional == flagged/q > flagged
    # Keeping only the unnormalized success branch loses half the mass error.
    changed = (q/2+eps, q/2, 1-q-eps)
    full_error = sum(abs(x-y) for x, y in zip(changed, ideal, strict=True))/2
    success_only = sum(abs(x-y) for x, y in zip(changed[:2], ideal[:2], strict=True))/2
    assert full_error == eps == 2*success_only


def test_uniform_finite_budget_certifies_faithfulness_and_bounded_modular_channel():
    n = p.first_sufficient(F(1, 1000))
    assert n == 369
    row, previous = p.budget(n), p.budget(n-1)
    assert row["successful_work_reference_faithful"]
    assert row["minimum_work_eigenvalue_lower"] > F(13332, 10**6)
    assert row["modular_channel_half_diamond_error_upper"] < F(1, 1000)
    assert previous["modular_channel_half_diamond_error_upper"] >= F(1, 1000)
    assert row["total_native_attempts_per_trial"] == 370
    pmin, pmax = row["herald_success_interval"]
    assert F(2, 15) <= pmin <= p.Q <= pmax <= F(1, 6)
    assert row["herald_failure_interval"] == (1-pmax, 1-pmin)
    assert row["expected_attempts_with_full_warmup_repeated_each_trial_upper"] == F(370)/pmin
    assert p.budget(369, F(8))["modular_channel_half_diamond_error_upper"] == row["modular_channel_half_diamond_error_upper"]/2
    assert p.first_sufficient(F(1, 1000), coverage=F(1, 2)) > n
    assert not p.budget(0)["successful_work_reference_faithful"]
    assert p.budget(0)["modular_channel_half_diamond_error_upper"] is None


def test_noncommuting_herald_log_bound_and_entangled_probe_channel():
    degree, sigma = p.anchor()
    columns = p.moving_columns("C")
    # Choose an actual native coherence that survives the herald between
    # unequal Palm eigenvalues. This is a noncommuting input/control, not
    # an assertion that a diagonal-only replay checks the quantum theorem.
    i, j = next((i, j) for (i, j), outputs in sorted(columns.items())
                if i < j and any(degree[a] != degree[b] for a, b in outputs))
    eps = F(1, 10000)
    rho = np.eye(66)/66
    rho[i, j] = rho[j, i] = float(eps)
    heralded = np.zeros((66, 66))
    for (a, b), outputs in columns.items():
        for (x, y), value in outputs.items():
            heralded[x, y] += float(value)*rho[a, b]
    probability = np.trace(heralded)
    state = heralded/probability
    target = np.diag([float(x) for x in sigma])
    assert np.linalg.norm(state@target-target@state, 2) > 1e-10
    t = float(66*eps)
    assert np.linalg.eigvalsh(heralded/float(p.Q)-(1-t)*target).min() >= -1e-15
    assert np.linalg.eigvalsh((1+t)*target-heralded/float(p.Q)).min() >= -1e-15
    values, vectors = np.linalg.eigh(state)
    assert values.min() > 0
    log_state = (vectors*np.log(values))@vectors.T
    centered = log_state+np.log(probability/float(p.Q))*np.eye(66)
    log_target = np.diag(np.log(np.diag(target)))
    assert np.linalg.norm(centered-log_target, 2) <= t/(1-t)
    spectrum, vectors = np.linalg.eigh(centered)
    s = 16
    actual = (vectors*np.exp(-1j*s*spectrum))@vectors.T
    ideal = np.diag(np.exp(-1j*s*np.diag(log_target)))
    channel_bound = min(1, s*t/(1-t))
    assert np.linalg.norm(actual-ideal, 2) <= channel_bound
    # Keep an arbitrary probe reference, here an entangled two-dimensional one.
    probe = np.zeros((66, 2), complex)
    probe[i, 0] = probe[j, 1] = 1/np.sqrt(2)
    left, right = (actual@probe).ravel(), (ideal@probe).ravel()
    output_difference = np.outer(left, left.conj())-np.outer(right, right.conj())
    assert np.abs(np.linalg.eigvalsh(output_difference)).sum()/2 <= channel_bound


def test_fixed_seam_and_zero_warmup_do_not_supply_the_claimed_full_reference():
    states = p.relational.quantum.STATES
    support = {p.native.native(state, (0, 1), coin)[0]
               for state in states for coin in (0, 1)
               if p.native.native(state, (0, 1), coin)[1][1]}
    assert len(support) == 20 < 66
    initial = tuple(int(i in (0, 2)) for i in range(12))
    assert min(p.diagonal_trial(initial, 0)["conditional"]) == 0


@pytest.mark.parametrize("attempts,window,coverage", [
    (-1, F(16), F(1)), (True, F(16), F(1)), (1., F(16), F(1)),
    (3, -F(1), F(1)), (3, 16., F(1)), (3, F(16), 0),
    (3, F(16), 1.), (3, F(16), F(11, 10)),
])
def test_budget_rejects_undeclared_counts_windows_and_coverage(attempts, window, coverage):
    with pytest.raises(ValueError):
        p.budget(attempts, window, coverage)


def test_retained_preparation_receipt_and_outward_error_bounds():
    assert json.loads(p.OUTPUT.read_text(encoding="utf-8")) == p.run()
    exact = p.budget(369)
    compact = p._compact(exact)
    assert F(compact["minimum_work_eigenvalue_lower"]) <= exact["minimum_work_eigenvalue_lower"]
    assert F(compact["modular_channel_half_diamond_error_upper"]) >= exact["modular_channel_half_diamond_error_upper"]
