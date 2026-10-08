"""Original-input controls for exact finite stationary readout moments.

The linear oracle evolves a twelve-dimensional one-particle kernel.  The
nonlinear oracle enumerates the 66 original two-occupant configurations and
forms their full transition matrix.  Neither uses the producer's moment
assembly or the independent receipt verifier.
"""
from __future__ import annotations

import ast
import copy
from fractions import Fraction as F
import importlib.util
import itertools
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("exact_moment_producer", HERE / "build.py")
producer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)
verify_spec = importlib.util.spec_from_file_location("exact_moment_verifier", HERE / "verify.py")
verifier = importlib.util.module_from_spec(verify_spec)
verify_spec.loader.exec_module(verifier)


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


@pytest.mark.parametrize("stride,key", [(1, "record_average_covariance"),
                                       (2, "two_attempts_per_record_covariance")])
def test_record_covariance_from_native_path_sums(carrier, edges, stride, key):
    """Accumulate actual record sums, without a lag multiplicity formula."""
    states = [tuple(int(i in pair) for i in range(12))
              for pair in itertools.combinations(range(12), 2)]
    index = {state: i for i, state in enumerate(states)}
    transition = np.zeros((66, 66), dtype=object)
    value = np.zeros((66, 12), dtype=object)
    for i, state in enumerate(states):
        for a, b in edges:
            if state[a] != state[b]:
                value[i, a] += 1
                value[i, b] += 1
            for coin in (False, True):
                target = list(state)
                target[a], target[b] = carrier.integer_nearest_agreement(
                    state[a], state[b], ceiling_to_first=coin)
                transition[i, index[tuple(target)]] += 1
    assert all(sum(row) == 60 for row in transition)
    assert all(sum(column) == 60 for column in transition.T)
    windows = [1, 5, 11]
    record = producer.finite_control(2, carrier, windows)["readouts"]["local_mismatch"]
    mean = value.sum(axis=0) / F(66)
    gram = value.T @ value
    # For paths ending at each state, retain their summed record sums;
    # separately retain the global sum of squared record sums.
    first, second, paths_per_state = value.copy(), gram.copy(), 1
    for count in range(1, max(windows) + 1):
        if count in windows:
            covariance = second / F(66 * paths_per_state * count**2) - np.outer(mean, mean)
            np.testing.assert_array_equal(exact(record[key][str(count)]), covariance)
        if count == max(windows):
            break
        for _ in range(stride):
            first = transition.T @ first
        cross = first.T @ value
        paths_per_state *= 60**stride
        second = 60**stride * second + cross + cross.T + paths_per_state * gram
        first += paths_per_state * value


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


@pytest.mark.parametrize("raised", range(1, 12))
def test_mismatch_projection_has_independent_occupancy_factor(carrier, edges, raised):
    row = producer.finite_control(raised, carrier, [1])
    # Cov(n_a*n_b,z)/kappa = (r-1)/10*(e_a+e_b-2*1/12).
    # Thus projecting sum_j(n_i+n_j-2*n_i*n_j) multiplies H*Pi by (6-r)/5.
    h = 5 * np.eye(12, dtype=object)
    for a, b in edges:
        h[a, b] += 1
        h[b, a] += 1
    projection = F(6 - raised, 5) * h @ (
        np.eye(12, dtype=object) - np.ones((12, 12), dtype=object) / F(12)
    )
    readout = row["readouts"]["local_mismatch"]
    np.testing.assert_array_equal(exact(readout["density_projection_B"]), projection)
    assert [F(value) for value in readout["mean"]] == [F(5 * raised * (12 - raised), 66)] * 12
    if raised == 6:
        assert readout["residual_covariance"] == readout["instant_covariance"]


def test_cross_covariance_respects_large_origins_and_integer_coordinate_maps():
    x = np.array([[0, 3], [7, -2], [-1, 6]], dtype=object)
    y = np.array([[2], [-3], [5]], dtype=object)
    # Pairwise differences are an independent centering formula.
    expected = sum((np.outer(x[i] - x[j], y[i] - y[j])
                    for i in range(3) for j in range(3)),
                   np.zeros((2, 1), dtype=object)) / F(18)
    transform = np.array([[1, 2], [-3, 1]], dtype=object)
    x = (x @ transform) * 2**70 + np.array([2**120 + 1, -2**130], dtype=object)
    y = -2**65 * y + 2**150
    np.testing.assert_array_equal(producer.centered_cov(x, y),
                                  -(2**135) * transform.T @ expected)


@pytest.mark.parametrize("bad", [True, np.bool_(False), 1.0, 0.5, F(1, 2),
                                 float("nan"), float("inf"), None])
def test_original_noninteger_scalars_are_rejected(bad):
    values = [[2**60 + 1], [bad]]
    with pytest.raises(ValueError):
        producer.centered_cov(values, [[1], [2]])
    with pytest.raises(ValueError):
        producer.matrix(values)


def test_masked_and_malformed_populations_are_rejected():
    masked = np.ma.array([[1], [2]], mask=[[False], [True]])
    for values in (masked, np.array([1, 2]), np.zeros((0, 1), dtype=int)):
        with pytest.raises(ValueError):
            producer.centered_cov(values, values)
    with pytest.raises(ValueError):
        producer.centered_cov([[1], [2]], [[1]])
    with pytest.raises(ValueError):
        producer.matrix(masked)


@pytest.mark.parametrize("denominator", [True, 0, -1, 2.0, F(1, 2)])
def test_fraction_matrix_never_truncates_its_denominator(denominator):
    with pytest.raises(ValueError):
        producer.matrix([[3]], denominator)


@pytest.mark.parametrize("raised", [True, 6.0, np.int64(6), -1, 0, 12])
def test_finite_occupancy_requires_a_supported_exact_json_integer(carrier, edges, raised):
    with pytest.raises(ValueError):
        producer.finite_control(raised, carrier, [1])
    with pytest.raises(ValueError):
        verifier.check_finite({"raised": raised}, edges, carrier.antipode(), [1])


@pytest.mark.parametrize("windows", [[], [0], [-1], [1, 1], [True], [1.0],
                                    [np.int64(1)], range(1, 3), [[1]]])
def test_record_windows_reject_changed_numeric_meanings(carrier, edges, windows):
    with pytest.raises(ValueError):
        producer.finite_control(6, carrier, windows)
    with pytest.raises(ValueError):
        verifier.check_finite({"raised": 6}, edges, carrier.antipode(), windows)


def test_compressed_dynamics_reject_an_incompatible_native_primitive(carrier):
    original = dict(seams=carrier.seams, laplacian=carrier.laplacian,
                    antipode=carrier.antipode,
                    integer_nearest_agreement=carrier.integer_nearest_agreement)
    changed_repair = SimpleNamespace(**{**original, "integer_nearest_agreement":
                                       lambda a, b, **kwargs: (a, b)})
    with pytest.raises(ValueError, match="repair"):
        producer.finite_control(6, changed_repair, [1])
    changed_laplacian = SimpleNamespace(**{**original, "laplacian":
                                          lambda: np.zeros((12, 12), dtype=int)})
    with pytest.raises(ValueError, match="Laplacian"):
        producer.finite_control(6, changed_laplacian, [1])


@pytest.fixture(scope="module")
def extended_row(carrier):
    return producer.finite_control(2, carrier, [1, 5, 10])


@pytest.mark.parametrize("raised", range(1, 12))
def test_long_windows_pass_independent_replay_after_json_roundtrip(carrier, edges, extended_row, raised):
    row = extended_row if raised == 2 else producer.finite_control(raised, carrier, [1, 5, 10])
    restored = json.loads(json.dumps(row, allow_nan=False))
    checked = verifier.check_finite(restored, edges, carrier.antipode(), [1, 5, 10])
    assert checked["raised"] == row["raised"]
    assert checked["configurations"] == row["configurations"]


@pytest.mark.parametrize("kind", ["late_lag", "late_window", "projection", "missing_lag"])
def test_independent_long_window_replay_rejects_false_moments(carrier, edges, extended_row, kind):
    row = copy.deepcopy(extended_row)
    readout = row["readouts"]["local_mismatch"]
    if kind == "late_lag":
        readout["lag_covariance"][18][0][3] = "0"
    elif kind == "late_window":
        readout["two_attempts_per_record_covariance"]["10"][0][0] = "0"
    elif kind == "projection":
        readout["density_projection_B"][0][0] = "0"
    else:
        readout["lag_covariance"].pop()
    with pytest.raises(ValueError):
        verifier.check_finite(row, edges, carrier.antipode(), [1, 5, 10])


def test_detached_verifier_replays_and_rejects_without_the_producer(tmp_path, carrier, extended_row):
    detached = tmp_path / "code/native_geometric_source/verify.py"
    detached.parent.mkdir(parents=True)
    detached.write_bytes((HERE / "verify.py").read_bytes())
    incidence = producer.VENDOR / "oph_fpe/dynamics/self_readback_repair_closure.py"
    copied_incidence = tmp_path / incidence.relative_to(producer.RER)
    copied_incidence.parent.mkdir(parents=True)
    copied_incidence.write_bytes(incidence.read_bytes())
    (tmp_path / "input.json").write_text(json.dumps({
        "row": extended_row, "antipodes": list(carrier.antipode()), "windows": [1, 5, 10]
    }, allow_nan=False), encoding="utf-8")
    runner = tmp_path / "replay.py"
    runner.write_text('''from copy import deepcopy
import importlib.util
import json
from pathlib import Path
root = Path(__file__).parent
assert not (root / "code/native_geometric_source/build.py").exists()
spec = importlib.util.spec_from_file_location("detached", root / "code/native_geometric_source/verify.py")
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)
packet = json.loads((root / "input.json").read_text())
edges = v.source_edges()
def check(row):
    return v.check_finite(row, edges, packet["antipodes"], packet["windows"])
good = check(packet["row"])
rejected = []
for kind in ("lag18", "projection", "clock", "census"):
    row = deepcopy(packet["row"])
    drive = row["readouts"]["local_drive"]
    if kind == "lag18": drive["lag_covariance"][18][0][0] = "0"
    elif kind == "projection": row["readouts"]["local_mismatch"]["density_projection_B"][0][0] = "0"
    elif kind == "clock": drive["two_attempts_per_record_covariance"]["10"] = drive["record_average_covariance"]["10"]
    else: drive["lag_covariance"].pop()
    try: check(row)
    except ValueError: rejected.append(kind)
print(json.dumps({"accepted_occupancy": good["raised"], "rejected": rejected}))
''', encoding="utf-8")
    result = subprocess.run([sys.executable, "-E", "-P", "-W", "error", str(runner)],
                            cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == {"accepted_occupancy": 2,
                                        "rejected": ["lag18", "projection", "clock", "census"]}
