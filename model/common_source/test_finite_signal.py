"""Independent target construction and interval checks for the exact signal."""
from fractions import Fraction as Q
import json

import numpy as np
import pytest
from scipy.linalg import expm

from . import finite_signal as fs
from . import gauge_matter as gauge
from . import scalar_coupling as scalar


def numeric(vector):
    return np.array([complex(float(real), float(imag)) for real, imag in vector])


def reference_branches(enabled):
    # Independent field quadratures and projected oscillator squares, not
    # the certificate's rational stencil, construct this target.
    h = scalar.ideal_hamiltonian(2)
    initial = np.eye(16, dtype=complex)[:, 0]
    first, second = expm(-1j*h), expm(-.5j*h)
    x0 = np.kron(np.eye(8), np.array([[0, 1], [1, 0]]))
    x2 = np.kron(np.eye(2), np.kron(np.array([[0, 1], [1, 0]]), np.eye(4)))
    feedback = .6*np.eye(16)-.8j*x0
    result = {}
    for a in (0, 1):
        p = np.diag([int((i & 1) == a) for i in range(16)])
        v = second@(feedback if enabled and a else np.eye(16))@p@first@initial
        for b in (0, 1):
            result[a, b] = (np.eye(16)+(2*b-1)*x2)@v/2
    return result


def test_exact_stencil_equals_independent_projected_full_hamiltonian():
    h = fs.hamiltonian_matrix()
    assert gauge.dagger(h) == h
    numerical = gauge.numerical_matrix(h)
    assert np.linalg.norm(numerical-scalar.ideal_hamiltonian(2), 2) < 2e-15
    assert np.linalg.norm(numerical, 2) <= float(fs.HAMILTONIAN_NORM_UPPER)
    assert sum(map(Q, ("2", "1", "219/64", "9/64", "1/4"))) == Q(109, 16)
    # The corner retains the projected-square correction; using Q_d² and
    # P_d² without it gives a different scalar action and must not pass.
    assert h[12][12] == (Q(219, 64), Q(0))
    assert h[4][0] == (Q(1, 4), Q(0))
    assert h[5][1] == (-Q(1, 4), Q(0))


@pytest.mark.parametrize("degree", (12, 24, 40))
@pytest.mark.parametrize("time", (Q(1), Q(1, 2)))
def test_taylor_remainder_against_independent_matrix_exponential(degree, time):
    approximate, epsilon = fs.taylor_action(fs.initial_vector(), time, degree)
    exact_numerical = expm(-1j*float(time)*scalar.ideal_hamiltonian(2))[:, 0]
    # The explicit numerical floor belongs only to this corroborating test.
    assert np.linalg.norm(numeric(approximate)-exact_numerical) <= float(epsilon)+3e-14


@pytest.mark.parametrize("enabled", (False, True))
def test_every_retained_branch_matches_independent_quantum_instrument(enabled):
    actual, _, _, eta = fs.branch_vectors(enabled)
    expected = reference_branches(enabled)
    assert [record for record, _ in actual] == [(0, 0), (0, 1), (1, 0), (1, 1)]
    for record, vector in actual:
        assert np.linalg.norm(numeric(vector)-expected[record]) < float(eta)+3e-14
    assert abs(sum(float(fs.norm_squared(v)) for _, v in actual)-1) < 3e-14


def test_feedback_preserves_exact_record_probability_before_second_evolution():
    first, _ = fs.taylor_action(fs.initial_vector(), Q(1))
    for outcome in (0, 1):
        selected = fs.first_projection(first, outcome)
        controlled = fs.controlled_feedback(selected, outcome, True)
        assert fs.norm_squared(controlled) == fs.norm_squared(selected)
    # The second projection is quantum coordinate sign, not bare occupation.
    assert fs.norm_squared(fs.scalar_sign_projection(fs.initial_vector(), 1)) == Q(1, 2)


def test_rational_intervals_contain_independent_probabilities_and_exclude_zero():
    packet = fs.certificate()
    json.loads(json.dumps(packet))
    values = []
    for row in packet["configurations"]:
        expected = reference_branches(row["feedback_enabled"])
        p = sum(float(np.vdot(vector, vector).real)
                for (_, b), vector in expected.items() if b == 1)
        interval = list(map(Q, row["positive_scalar_continuous_probability_interval"]))
        assert float(interval[0])-3e-14 <= p <= float(interval[1])+3e-14
        values.append(p)
        for branch in row["retained_branches"]:
            vector = expected[tuple(branch["record"])]
            actual_probability = float(np.vdot(vector, vector).real)
            low, high = map(Q, branch["continuous_probability_interval"])
            assert float(low)-3e-14 <= actual_probability <= float(high)+3e-14
    low, high = map(Q, packet["continuous_contrast_interval"])
    assert low < high < 0
    assert float(low)-3e-14 <= values[0]-values[1] <= float(high)+3e-14
    source_low, source_high = map(Q, packet["ideal_source_word_contrast_interval"])
    assert source_low < low < high < source_high < 0
    assert packet["continuous_nonzero_response_certified"]
    assert packet["ideal_source_word_nonzero_response_certified"]
    assert not packet["scope"]["comparison_with_nature"]


def test_exact_refinement_comparison_requires_no_floating_oracle():
    # Compare exact polynomial centers at two degrees using the proved
    # remainder bounds. This check itself is entirely rational.
    for enabled in (False, True):
        lower, _, _, eta32 = fs.branch_vectors(enabled, 32)
        finer, _, _, eta64 = fs.branch_vectors(enabled, 64)
        centers = [sum((fs.norm_squared(vector) for (_, b), vector in collection if b == 1), Q(0))
                   for collection in (lower, finer)]
        radius = 4*(eta32+eta64)+2*(eta32*eta32+eta64*eta64)
        assert abs(centers[0]-centers[1]) <= radius


def test_zero_taylor_degree_and_coarse_source_bounds_do_not_claim_resolution():
    coarse = fs.certificate(degree=0, steps=1)
    assert coarse["continuous_contrast_interval"] == ["-1", "1"]
    assert not coarse["continuous_nonzero_response_certified"]
    assert not coarse["ideal_source_word_nonzero_response_certified"]
    tight_continuous_only = fs.certificate(degree=64, steps=1)
    assert tight_continuous_only["continuous_nonzero_response_certified"]
    assert not tight_continuous_only["ideal_source_word_nonzero_response_certified"]


@pytest.mark.parametrize("call", (lambda: fs.certificate(degree=True),
                                 lambda: fs.certificate(steps=True),
                                 lambda: fs.certificate(steps=0),
                                 lambda: fs.taylor_action(fs.initial_vector(), .5),
                                 lambda: fs.taylor_action(fs.initial_vector(), Q(-1)),
                                 lambda: fs.taylor_action(fs.initial_vector(), Q(1), -1),
                                 lambda: fs.branch_vectors(1),
                                 lambda: fs.branch_vectors(1, 64),
                                 lambda: fs.controlled_feedback(fs.initial_vector(), 2, True),
                                 lambda: fs.scalar_sign_projection(fs.initial_vector(), True)))
def test_invalid_certificate_inputs_rejected(call):
    with pytest.raises(ValueError):
        call()
