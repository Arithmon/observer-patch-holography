"""Independent local-action and original-vacuum checks for split evolution.

The oracle builds the edge form directly, propagates a dense phase-space
matrix, and carries the original action vacuum through that matrix. It does
not use the production recurrence or its tensor spectral calculation.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations, product
import json
from math import comb, sqrt
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import eigh

from .local_evolution import LocalEvolution, temporal_source, terminal_branch, transfer_bound


def dense_model(q):
    with localcontext() as context:
        context.prec = 70
        phi = (1 + Decimal(5).sqrt()) / 2
        nodes = sorted((b * phi) % 1 for b in range(q)) + [Decimal(1)]
        gaps = np.array([float(b - a) for a, b in zip(nodes, nodes[1:])])
        x = np.array([float(a) for a in nodes[1:-1]])
    mass = (gaps[:-1] + gaps[1:]) / 2
    n = q - 1
    incidence = np.zeros((q, n))
    for edge in range(q):
        if edge:
            incidence[edge, edge - 1] = -1
        if edge < n:
            incidence[edge, edge] = 1
    stiffness = incidence.T @ np.diag(1 / gaps) @ incidence
    whitened = stiffness / np.sqrt(np.outer(mass, mass))
    identity = np.eye(n)
    action = (np.eye(n**3) + np.kron(whitened, np.eye(n*n))
              + np.kron(identity, np.kron(whitened, identity))
              + np.kron(np.eye(n*n), whitened))

    # Direct polynomial profile and its beta-integral normalization, with
    # neither the producer's sine-coefficient recurrence nor its eigenvectors.
    normalization = sqrt(comb(34, 17) / (18 * 2**17))
    def profile(z):
        return np.sin(np.pi*z) * (1 + np.cos(np.pi*z))**8 / normalization

    transverse = np.sqrt(2 * mass) * np.sin(np.pi*x)
    preparation = np.kron(np.sqrt(mass)*profile(x)*(x < .45),
                          np.kron(transverse, transverse))
    detector = np.kron(np.sqrt(mass)*profile(1-x)*(x > .55),
                       np.kron(transverse, transverse))
    values, vectors = eigh(action)
    omega = np.sqrt(values)
    position_covariance = (vectors / (2 * omega)) @ vectors.T
    momentum_covariance = (vectors * (omega / 2)) @ vectors.T
    zero = np.zeros_like(action)
    covariance = np.block([[position_covariance, zero],
                           [zero, momentum_covariance]])
    return action, preparation, detector, covariance, omega, vectors


def dense_split(q, step, steps):
    action, preparation, detector, covariance, _, _ = dense_model(q)
    size = len(action)
    identity, zero = np.eye(size), np.zeros_like(action)
    tau = float(step)
    kick = np.block([[identity, zero], [-tau * action / 2, identity]])
    drift = np.block([[identity, tau * identity], [zero, identity]])
    evolution = np.linalg.matrix_power(kick @ drift @ kick, steps)
    initial_mean = np.concatenate((np.zeros(size), 4 * preparation))
    final_mean = evolution @ initial_mean
    final_covariance = evolution @ covariance @ evolution.T
    mean = detector @ final_mean[:size]
    variance = detector @ final_covariance[:size, :size] @ detector
    probability = (1 + np.exp(-variance / 2) * np.sin(mean)) / 2
    stationary_variance = detector @ covariance[:size, :size] @ detector
    return mean, variance, probability, stationary_variance


@pytest.mark.parametrize("q", (3, 5))
def test_stencil_and_profiles_match_independent_edge_energy(q):
    action, preparation, detector, *_ = dense_model(q)
    local = LocalEvolution(q)
    assert tuple(local.shape) == (q-1, q-1, q-1)
    np.testing.assert_allclose(local.preparation.reshape(-1), preparation,
                               rtol=0, atol=2e-14)
    np.testing.assert_allclose(local.detector.reshape(-1), detector,
                               rtol=0, atol=2e-14)
    values = np.random.default_rng(4821).normal(size=(3, *local.shape))
    expected = (values.reshape(3, -1) @ action.T).reshape(values.shape)
    np.testing.assert_allclose(local.apply_action(values), expected,
                               rtol=2e-14, atol=3e-12)


def test_one_action_has_only_onsite_and_adjacent_axis_support():
    local = LocalEvolution(5)
    source = (1, 1, 1)
    delta = np.zeros(local.shape)
    delta[source] = 1
    result = local.apply_action(delta)
    touched = {tuple(point) for point in np.argwhere(result != 0)}
    expected = {source}
    for axis, direction in product(range(3), (-1, 1)):
        point = list(source)
        point[axis] += direction
        expected.add(tuple(point))
    assert touched == expected
    assert result[source] > 0
    assert all(result[point] < 0 for point in touched - {source})


@pytest.mark.parametrize("q,step,steps", ((3, Fraction(1, 32), 0),
                                         (3, Fraction(1, 32), 9),
                                         (5, Fraction(1, 16), 12)))
def test_local_trajectory_matches_dense_canonical_matrix_and_vacuum(q, step, steps):
    expected_mean, expected_variance, expected_probability, _ = dense_split(q, step, steps)
    actual = LocalEvolution(q).evolve(step, steps)
    assert actual["mean"] == pytest.approx(expected_mean, rel=0, abs=2e-12)
    assert actual["variance"] == pytest.approx(expected_variance, rel=0, abs=2e-12)
    assert actual["plus_probability"] == pytest.approx(expected_probability, rel=0, abs=2e-12)
    assert Fraction(actual["step"]) == step
    assert actual["steps"] == steps
    assert Fraction(actual["model_time"]) == steps * step
    assert actual["oscillator_count"] == (q - 1)**3


def test_original_vacuum_is_evolved_not_replaced_by_a_stationary_covariance():
    q, step, steps = 5, Fraction(1, 16), 12
    mean, variance, probability, original_variance = dense_split(q, step, steps)
    actual = LocalEvolution(q).evolve(step, steps)
    assert variance - original_variance > .001
    unchanged_vacuum_probability = (1 + np.exp(-original_variance/2) * np.sin(mean)) / 2
    assert abs(probability - unchanged_vacuum_probability) > 1e-5
    assert actual["variance"] == pytest.approx(variance, rel=0, abs=2e-12)
    assert actual["plus_probability"] == pytest.approx(probability, rel=0, abs=2e-12)


def test_halving_step_converges_to_same_action_continuous_detector():
    q, time = 3, Fraction(3, 4)
    _, preparation, detector, _, omega, vectors = dense_model(q)
    f, g = vectors.T @ preparation, vectors.T @ detector
    exact_mean = 4 * np.dot(f*g, np.sin(float(time)*omega)/omega)
    exact_variance = np.dot(g*g, 1/(2*omega))
    exact_probability = (1 + np.exp(-exact_variance/2) * np.sin(exact_mean)) / 2
    errors = []
    for denominator in (16, 32, 64):
        step = Fraction(1, denominator)
        result = LocalEvolution(q).evolve(step, int(time / step))
        errors.append(abs(result["plus_probability"] - exact_probability))
        assert result["variance"] >= exact_variance - 2e-14
    assert all(3.5 < earlier/later < 4.6 for earlier, later in zip(errors, errors[1:]))


def history_law(response):
    """Independent native integer updates composed with any signed bias."""
    rows = {}
    for occupied in combinations(range(12), 2):
        initial = tuple(int(k in occupied) for k in range(12))
        for first_coin, bit, second_coin in product((0, 1), (-1, 1), (0, 1)):
            first = list(initial)
            total = first[0] + first[1]
            first[0] = (total + first_coin) // 2
            first[1] = total - first[0]
            current = initial[0] - first[0]
            other = 5 if bit == 1 else 1
            final = list(first)
            total = final[0] + final[other]
            final[0] = (total + second_coin) // 2
            final[other] = total - final[0]
            second_current = first[0] - final[0]
            key = (initial, first_coin, current, bit, second_coin, second_current)
            rows[key] = (Fraction(1, 2) + bit * current * response) / 264
    return rows


def test_split_bias_keeps_exact_native_history_error_and_control_contract():
    local = LocalEvolution(3)
    results = [local.evolve(Fraction(1, denominator), 3*denominator//4)
               for denominator in (16, 32)]
    biases = [Fraction.from_float(result["plus_probability"] - .5) for result in results]
    laws = [history_law(bias) for bias in biases]
    for bias, law in zip(biases, laws):
        assert len(law) == 528
        assert sum(law.values()) == 1
        assert sum(weight * key[-1] for key, weight in law.items()) == Fraction(5, 132)*bias
        assert sum(weight for key, weight in law.items() if key[3] == 1) == Fraction(1, 2)
        reversed_law = history_law(-bias)
        assert sum(weight * key[-1] for key, weight in reversed_law.items()) == -Fraction(5, 132)*bias
    variation = sum(abs(laws[0][key] - laws[1][key]) for key in laws[0]) / 2
    assert variation == Fraction(5, 33) * abs(biases[0] - biases[1])


def test_exact_step_rounding_and_error_transfer_at_window_edges_and_ties():
    parent = temporal_source()
    for level in parent["levels"]:
        step = Fraction(level["time_step"])
        tie_time = (level["first_window_step"] + Fraction(1, 2)) * step
        for reference in (Fraction(19, 20), Fraction(39, 40), Fraction(1), tie_time):
            result = transfer_bound(level, reference)
            quotient = reference / step
            rounded = (2*quotient.numerator + quotient.denominator) // (2*quotient.denominator)
            expected_steps = max(level["first_window_step"], min(level["last_window_step"], rounded))
            actual_time = expected_steps * step
            assert result["steps"] == expected_steps
            assert Fraction(result["model_time"]) == actual_time
            assert Fraction(19, 20) <= actual_time <= 1
            assert abs(actual_time - reference) <= step
            if reference == tie_time:
                assert expected_steps == level["first_window_step"] + 1
            same_time = Fraction(level["same_grid_total_error_upper"])
            nearby = min(Fraction(1), same_time + 2*abs(actual_time-reference))
            assert Fraction(result["reference_time_distance"]) == abs(actual_time-reference)
            assert Fraction(result["reference_probability_error_upper"]) == nearby
            assert Fraction(result["same_time_history_tv_upper"]) == 5*same_time/33
            assert Fraction(result["same_time_current_error_upper"]) == 5*same_time/132
            assert Fraction(result["reference_history_tv_upper"]) == 5*nearby/33
            assert Fraction(result["reference_current_error_upper"]) == 5*nearby/132
            assert Fraction(result["split_native_current_magnitude_lower"]) == (
                5*Fraction(level["split_response_lower"])/132)
    finest = transfer_bound(parent["levels"][-1])
    assert Fraction(finest["reference_current_error_upper"]) < Fraction(19, 10000)
    # The actually executed coarse grids must retain the parent's trivial
    # spatial bound; successful software execution cannot tighten it.
    for level in parent["levels"][:2]:
        assert Fraction(transfer_bound(level)["same_time_probability_error_upper"]) == 1


def test_terminal_squeezed_instrument_preserves_normalization_and_signed_preparation():
    value = LocalEvolution(5).evolve(Fraction(1, 16), 12)
    for current in (-1, 0, 1):
        branches = [terminal_branch(value, current, outcome) for outcome in (-1, 1)]
        assert sum(branch["probability"] for branch in branches) == pytest.approx(1, abs=2e-15)
        for branch in branches:
            raw = np.asarray(branch["gram"])
            gram = raw[:, :, 0] + 1j*raw[:, :, 1]
            raw = np.asarray(branch["coefficients"])
            coefficients = raw[:, 0] + 1j*raw[:, 1]
            np.testing.assert_allclose(gram, gram.conj().T, rtol=0, atol=1e-15)
            assert np.linalg.eigvalsh(gram).min() >= 0
            assert coefficients.conj() @ gram @ coefficients == pytest.approx(1, abs=2e-15)
            probability = .5 + branch["outcome"] * current * (value["plus_probability"]-.5)
            assert branch["probability"] == pytest.approx(probability, abs=2e-15)
            assert branch["q"] == value["q"] and branch["steps"] == value["steps"]


@pytest.mark.parametrize("q", (True, np.bool_(False), 1, 35, 3.0))
def test_invalid_local_grid_rejected(q):
    with pytest.raises(ValueError):
        LocalEvolution(q)


@pytest.mark.parametrize("step,steps", ((.01, 2), (Fraction(0), 2),
                                        (Fraction(-1, 32), 2), (Fraction(1), 2),
                                        (Fraction(1, 32), True), (Fraction(1, 32), 2.0),
                                        (Fraction(1, 32), -1), (Fraction(1, 32), 32769)))
def test_invalid_or_unstable_local_execution_rejected(step, steps):
    with pytest.raises(ValueError):
        LocalEvolution(3).evolve(step, steps)


def test_grid_step_resource_budget_rejected_before_advance():
    with pytest.raises(ValueError):
        LocalEvolution(34).evolve(Fraction(1, 32768), 32768)


def test_local_action_rejects_invalid_field_data():
    local = LocalEvolution(3)
    bad_fields = (np.zeros((2, 2)), np.zeros(local.shape, complex),
                  np.full(local.shape, np.nan), np.full(local.shape, np.inf),
                  np.ma.array(np.zeros(local.shape), mask=True))
    for values in bad_fields:
        with pytest.raises(ValueError):
            local.apply_action(values)


def test_retained_local_execution_and_analytic_refinements_have_distinct_scope():
    receipt = json.loads(Path(__file__).with_name("local_results.json").read_text())
    parent = temporal_source()
    assert receipt["all_pass"] and all(receipt["checks"].values())
    assert receipt["locally_executed_grids"] == [13, 34]
    assert receipt["refinement_bounds"] == [transfer_bound(level) for level in parent["levels"]]
    assert [row["q"] for row in receipt["refinement_bounds"]] == [13, 34, 89, 233]
    assert [row["local"]["q"] for row in receipt["executions"]] == [13, 34]
    for execution in receipt["executions"]:
        local = execution["local"]
        assert not local["full_per_site_event_log_retained"]
        assert local["executed_position_column_updates"] == 3*local["steps"]*local["oscillator_count"]
        assert local["checkpoints"][-1]["step"] == local["steps"]
        assert local["variance"] > execution["continuous_same_time"]["vacuum_variance"]
        for key in ("mean", "variance", "plus_probability"):
            assert local[key] == pytest.approx(execution["modal_oracle"][key], rel=0, abs=2e-10)
        sample = execution["sample"]
        assert sum(row["count"] for row in sample["record_counts"]) == sample["trials"] == 200000
        assert len(sample["record_counts"]) == 24
        assert sample["statistics"]["current"]["forecast"] == pytest.approx(
            5*(local["plus_probability"]-.5)/132, rel=0, abs=2e-16)
        assert len(execution["terminal_field_branches"]) == 6
    coarse = receipt["executions"][0]
    assert coarse["half_step_refinement"]["local"]["model_time"] == coarse["local"]["model_time"]
    assert coarse["half_step_refinement"]["absolute_probability_error"] < coarse["absolute_probability_error"]/3
