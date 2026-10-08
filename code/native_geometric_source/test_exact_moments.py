"""Original-input controls for exact finite stationary readout moments.

The linear oracle evolves a twelve-dimensional one-particle kernel.  The
nonlinear oracle enumerates the 66 original two-occupant configurations and
forms their full transition matrix.  Neither uses the producer's moment
assembly or the independent receipt verifier.
"""
from __future__ import annotations

import ast
from fractions import Fraction as F
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("exact_moment_producer", HERE / "build.py")
producer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)


@pytest.fixture(scope="module")
def carrier():
    configuration = json.loads((producer.EVIDENCE / "spec.json").read_text())
    native, _, _ = producer.native_modules(configuration)
    return native


@pytest.fixture(scope="module")
def edges():
    """Read original incidence, independently of the producer's carrier API."""
    source = producer.VENDOR / "oph_fpe/dynamics/self_readback_repair_closure.py"
    for node in ast.parse(source.read_text()).body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "ORIENTED_BASE_FACES"
            for target in node.targets
        ):
            faces = ast.literal_eval(node.value)
            return sorted({tuple(sorted((face[i], face[(i + 1) % 3])))
                           for face in faces for i in range(3)})
    raise AssertionError("original carrier incidence is missing")


def exact(value):
    return np.array([[F(entry) for entry in row] for row in value], dtype=object)


def linear_oracle(raised, edges, last_lag, name):
    """Uniform occupancy gives Cov(n)=kappa*(I-J/12).

    A fair endpoint coin replaces either endpoint by the seam mean in
    conditional expectation, so P=I-L/60.  Load is n-r/12, drive is L*n.
    At r=1 or 11 mismatch is H*n (or H*holes), H=5I+adjacency.
    """
    identity = np.eye(12, dtype=object)
    laplacian = np.zeros((12, 12), dtype=object)
    for a, b in edges:
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    transition = identity - laplacian / F(60)
    covariance = F(raised * (12 - raised), 132) * (
        identity - np.ones((12, 12), dtype=object) / F(12)
    )
    if name == "local_drive":
        covariance = laplacian @ covariance @ laplacian
    elif name == "local_mismatch":
        assert raised in (1, 11)
        readout = 10 * identity - laplacian
        covariance = readout @ covariance @ readout
    else:
        assert name == "load"
    result = []
    for _ in range(last_lag + 1):
        result.append(covariance)
        covariance = covariance @ transition
    return result


def averaged(moments, count, stride):
    total = count * moments[0]
    for lag in range(1, count):
        covariance = moments[stride * lag]
        total = total + (count - lag) * (covariance + covariance.T)
    return total / F(count * count)


def assert_record_moments(record, moments, windows):
    assert len(record["lag_covariance"]) == len(moments)
    for actual, expected in zip(record["lag_covariance"], moments, strict=True):
        np.testing.assert_array_equal(exact(actual), expected)
    np.testing.assert_array_equal(exact(record["instant_covariance"]), moments[0])
    for stride, key in ((1, "record_average_covariance"),
                        (2, "two_attempts_per_record_covariance")):
        assert set(record[key]) == {str(count) for count in windows}
        for count in windows:
            np.testing.assert_array_equal(
                exact(record[key][str(count)]), averaged(moments, count, stride)
            )


def test_original_half_filling_lag_eight_is_exact_and_positive(carrier, edges):
    row = producer.finite_control(6, carrier, [1, 2, 4, 5])
    actual = exact(row["readouts"]["local_drive"]["lag_covariance"][8])
    expected = linear_oracle(6, edges, 8, "local_drive")[8]
    assert expected[0, 0] == F(90750535351, 26730000000)
    # For every stationary reversible chain, C_8 is the Gram matrix of P^4 f.
    # The defective implementation gives a negative antipodal contrast.
    assert actual[0, 0] + actual[3, 3] - actual[0, 3] - actual[3, 0] >= 0
    np.testing.assert_array_equal(actual, expected)


def test_first_extended_window_cannot_have_negative_variance(carrier, edges):
    row = producer.finite_control(3, carrier, [1, 2, 4, 5])
    actual = exact(row["readouts"]["local_drive"]
                   ["two_attempts_per_record_covariance"]["5"])
    # The sum of all twelve drive records is identically zero in every state.
    # Before repair this returned variance is -35184372088832/10571923828125.
    assert sum(actual.flat) == 0
    expected = averaged(linear_oracle(3, edges, 8, "local_drive"), 5, 2)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("raised,windows", [(1, (1, 2, 4)), (6, (1, 2, 4)),
                                            (1, (1, 5, 10)), (6, (1, 5, 10))])
def test_all_linear_covariances_against_one_particle_kernel(carrier, edges, raised, windows):
    row = producer.finite_control(raised, carrier, windows)
    for name in ("load", "local_drive"):
        moments = linear_oracle(raised, edges, 2 * (max(windows) - 1), name)
        assert_record_moments(row["readouts"][name], moments, windows)
        np.testing.assert_array_equal(
            exact(row["readouts"][name]["residual_covariance"]),
            np.zeros((12, 12), dtype=object),
        )


@pytest.mark.parametrize("raised", [1, 11])
def test_mismatch_one_particle_and_one_hole_closed_form(carrier, edges, raised):
    windows = (1, 5, 10)
    row = producer.finite_control(raised, carrier, windows)
    moments = linear_oracle(raised, edges, 18, "local_mismatch")
    assert_record_moments(row["readouts"]["local_mismatch"], moments, windows)


def nonlinear_two_occupant_oracle(edges, last_lag):
    """Full original-state kernel: thirty stays plus thirty transpositions.

    This uses no polynomial closure.  Mismatch is evaluated directly on each
    binary state, and its stationary mean is independently 2*r*(12-r)/132
    per incident edge.  Matrix multiplication counts all length-t paths.
    """
    states = [frozenset(pair) for pair in itertools.combinations(range(12), 2)]
    index = {state: i for i, state in enumerate(states)}
    transition = np.zeros((66, 66), dtype=object)
    values = np.zeros((66, 12), dtype=object)
    for i, state in enumerate(states):
        transition[i, i] = 30
        for a, b in edges:
            target = state.symmetric_difference((a, b)) if (a in state) != (b in state) else state
            transition[i, index[target]] += 1
            if (a in state) != (b in state):
                values[i, a] += 1
                values[i, b] += 1
    assert all(sum(row) == 60 for row in transition)
    assert np.array_equal(transition, transition.T)
    mean = np.full(12, F(50, 33), dtype=object)
    assert all(F(sum(values[:, i]), 66) == mean[i] for i in range(12))
    propagated = values
    result = []
    for lag in range(last_lag + 1):
        result.append((values.T @ propagated) / F(66 * 60**lag) - np.outer(mean, mean))
        propagated = transition @ propagated
    return result


def test_nonlinear_mismatch_against_full_two_occupant_path_kernel(carrier, edges):
    windows = (1, 2, 5, 10)
    row = producer.finite_control(2, carrier, windows)
    moments = nonlinear_two_occupant_oracle(edges, 18)
    assert_record_moments(row["readouts"]["local_mismatch"], moments, windows)
    assert any(F(entry) != 0 for line in row["readouts"]["local_mismatch"]
               ["residual_covariance"] for entry in line)


def test_complement_preserves_covariances_and_flips_mismatch_projection(carrier):
    low = producer.finite_control(2, carrier, [1, 5, 10])
    high = producer.finite_control(10, carrier, [1, 5, 10])
    for name in ("load", "local_drive", "local_mismatch"):
        first, second = low["readouts"][name], high["readouts"][name]
        for key in ("mean", "lag_covariance", "instant_covariance", "residual_covariance",
                    "record_average_covariance", "two_attempts_per_record_covariance",
                    "cross_scale_covariance", "coarse_covariance"):
            assert first[key] == second[key]
        sign = -1 if name == "local_mismatch" else 1
        np.testing.assert_array_equal(exact(first["density_projection_B"]),
                                      sign * exact(second["density_projection_B"]))


def test_centered_covariance_does_not_overflow_before_fraction_conversion():
    values = np.array([[2**32, 0], [0, 2**32]], dtype=np.int64)
    expected = 2**62 * np.array([[1, -1], [-1, 1]], dtype=object)
    np.testing.assert_array_equal(producer.centered_cov(values, values), expected)


def test_centering_preserves_original_mixed_integer_values():
    values = np.array([[2**60 + 1], [np.int64(2**60)]], dtype=object)
    np.testing.assert_array_equal(producer.centered_cov(values, values),
                                  np.array([[F(1, 4)]], dtype=object))
