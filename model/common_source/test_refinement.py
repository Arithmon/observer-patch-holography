"""Independent algebra, operator and directed-arithmetic refinement checks."""
from fractions import Fraction as F
import json
import math

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from . import refinement as r
from . import scalar_coupling


def _mp_fraction(value):
    value = F(value)
    return mp.mpf(value.numerator) / value.denominator


def _independent_hamiltonian(d):
    """Physical infinite-oscillator matrix elements, then rectangular cutoff."""
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    q = (a + a.T) / np.sqrt(2)
    q2 = np.diag(np.arange(d) + .5)
    for n in range(d - 2):
        q2[n, n + 2] = q2[n + 2, n] = np.sqrt((n + 1) * (n + 2)) / 2
    n = np.diag(np.arange(d))
    identity = np.eye(d)
    total_n = np.kron(identity, n) + np.kron(n, identity)
    spring = 9 / 64 * (np.kron(identity, q2) + np.kron(q2, identity)
                       - 2 * np.kron(q, q))
    gauge = np.diag([0., 1., 1., 0.])
    for i in range(4):
        gauge[i ^ 1, i] = gauge[i ^ 2, i] = 1
    return (np.kron(total_n + np.eye(d * d) + spring, np.eye(4))
            + np.kron(np.eye(d * d), gauge)
            + np.sqrt(2) / 4 * np.kron(np.kron(identity, q), np.diag([1., -1., 1., -1.])))


def test_exact_vacuum_moment_by_independent_generating_function():
    z = sp.Symbol("z")
    generating = sp.sqrt(80 / (81 - z**2))
    value = generating
    for _ in range(4):
        value = z * sp.diff(value, z) + 2 * value
    exact = sp.simplify(value.subs(z, 1))
    assert r.initial_fourth_moment() == F(int(exact.p), int(exact.q))
    assert r.initial_fourth_moment() == F(143748821, 8192000)


def test_actual_normalized_vacuum_truncation_decreases_moment_and_tail_bound():
    with mp.workdps(80):
        previous = F(0)
        for d in (2, 4, 8, 16, 32):
            weights = [F(math.comb(2 * k, k), 4**k * 81**k) for k in range(d // 2)]
            truncated = sum(w * (2 * k + 2)**4 for k, w in enumerate(weights)) / sum(weights)
            assert previous < truncated < r.initial_fourth_moment()
            previous = truncated
            actual_tail = 1 - mp.sqrt(mp.mpf(80) / 81) * sum(_mp_fraction(w) for w in weights)
            assert 0 < actual_tail < _mp_fraction(r.preparation_tail_upper(d))


def test_moment_commutator_coefficient_inequalities_have_exact_positive_certificates():
    x, y = sp.symbols("x y")
    squeezing = 8 * x * (x + 2)**2 - ((x + 2)**4 - x**4)
    assert sp.expand(squeezing) == 8 * x**2 - 16
    linear = 4 * x**4 * (x + 1)**4 - (x - 1) * ((x + 1)**4 - x**4)**2
    coefficients = sp.Poly(sp.expand(linear.subs(x, y + 2)), y).all_coeffs()
    assert coefficients == [4, 64, 440, 1692, 3976, 5884, 5500, 3175, 959]
    assert all(c > 0 for c in coefficients)


@pytest.mark.parametrize("d", (2, 4, 6))
def test_interacting_number_moment_commutator_bound_on_every_finite_input(d):
    h = _independent_hamiltonian(d)
    if d in (2, 4):
        assert np.linalg.norm(h - scalar_coupling.ideal_hamiltonian(d), 2) < 3e-15
    weight = np.array([(i // 4 % d + i // (4 * d) + 2)**4 for i in range(4 * d * d)])
    normalized = 1j * h * (weight[None, :] - weight[:, None])
    normalized /= np.sqrt(weight[:, None] * weight[None, :])
    # A matrix inequality, checked on the whole compressed carrier rather
    # than only the initial vacuum. The general proof is in the companion doc.
    assert max(abs(np.linalg.eigvalsh(normalized))) < float(r.MOMENT_RATE)


def test_retained_gauge_branches_and_conditional_rotation_preserve_total_moment():
    d = 4
    rng = np.random.default_rng(604)
    state = rng.normal(size=(d * d, 4, 3)) + 1j * rng.normal(size=(d * d, 4, 3))
    state /= np.linalg.norm(state)
    weight = np.array([(i % d + i // d + 2)**4 for i in range(d * d)])
    before = np.sum(weight[:, None, None] * abs(state)**2)
    rotation = (3 * np.eye(4) - 4j * np.eye(4)[[1, 0, 3, 2]]) / 5
    after = 0
    for a in (0, 1):
        branch = state.copy()
        branch[:, [i for i in range(4) if i % 2 != a], :] = 0
        if a:
            branch = np.einsum("ij,sjr->sir", rotation, branch)
        after += np.sum(weight[:, None, None] * abs(branch)**2)
    assert after == pytest.approx(before, abs=2e-13)


@pytest.mark.parametrize("d", (2, 16, 64, 256))
def test_reported_directed_bound_dominates_independent_high_precision_expression(d):
    certificate = r.occupation_evolution_bound(d)
    with mp.workdps(90):
        moment = mp.mpf(143748821) / 8192000
        amplitude = mp.sqrt(moment)
        growth = mp.exp(mp.mpf(13) / 4 * mp.mpf(3) / 4)
        integral = amplitude * (growth - 1) / (mp.mpf(13) / 8)
        pair = d // 2
        omitted = mp.mpf(81) / 80 * math.comb(2 * pair, pair) / (mp.mpf(4)**pair * 81**pair)
        preparation = mp.sqrt(2 * omitted)
        leakage = mp.mpf(9) / (32 * (d + 2)) + 1 / (4 * mp.mpf(d + 2)**mp.mpf("1.5"))
        projected = preparation + integral * leakage
        outside = amplitude * growth / (d + 2)**2
        ideal_bound = mp.sqrt(projected**2 + outside**2)
        reported = _mp_fraction(certificate["pre_readout_stinespring_vector_error_upper"])
        assert ideal_bound <= reported
        assert reported - ideal_bound < mp.mpf("1e-28")


def test_unclipped_cutoff_bounds_decrease_without_a_rounding_floor():
    values = [F(r.occupation_evolution_bound(d)["pre_readout_stinespring_vector_error_upper"])
              for d in (4, 16, 64, 256, 1024)]
    assert all(a > b for a, b in zip(values, values[1:]))
    assert values[-1] < F(1, 100)
    zero_time = r.occupation_evolution_bound(4, F(0))
    assert F(zero_time["projected_duhamel_residual_integral_upper"]) == 0


def test_actual_sharp_finite_readout_is_not_compression_of_infinite_sign():
    from scipy.integrate import quad
    q = np.array([[0., 1.], [1., 0.]]) / np.sqrt(2)
    eigenvalues, eigenvectors = np.linalg.eigh(q)
    finite_sign = (eigenvectors * np.sign(eigenvalues)) @ eigenvectors.T
    # First two normalized Hermite functions: phi1(x)=sqrt(2)*x*phi0(x).
    compression_entry = quad(lambda x: 2 * np.sqrt(2 / np.pi) * x * np.exp(-x * x), 0, np.inf)[0]
    compressed_sign = compression_entry * np.array([[0., 1.], [1., 0.]])
    assert compression_entry == pytest.approx(np.sqrt(2 / np.pi))
    assert np.linalg.norm(finite_sign - compressed_sign, 2) > .2


def test_actual_sharp_readout_strong_convergence_diagnostic_on_vacuum_column():
    errors = []
    for d in (2, 4, 8, 16, 32, 64):
        a = np.diag(np.sqrt(np.arange(1, d)), 1)
        eigenvalues, eigenvectors = np.linalg.eigh((a + a.T) / np.sqrt(2))
        finite_column = ((eigenvectors * np.sign(eigenvalues)) @ eigenvectors.T)[:, 0]
        # Exact Hermite matrix elements <2k+1|sign(Q)|0>, from integrating
        # the Hermite recurrence at zero. The infinite column has norm one.
        infinite_entries = np.zeros(d)
        for k in range(d // 2):
            central = math.comb(2 * k, k) / 4**k
            infinite_entries[2 * k + 1] = (-1)**k * np.sqrt(2 * central / (np.pi * (2 * k + 1)))
        errors.append(np.sqrt(2 - 2 * np.dot(finite_column, infinite_entries)))
    assert all(a > b for a, b in zip(errors, errors[1:]))
    # This finite numerical diagnostic is not promoted to a uniform rate.
    assert not r.sharp_readout_status()["complete_history_finite_error_bound_supplied"]


def test_scope_cannot_transfer_pre_readout_bound_to_final_detector():
    report = r.report((4, 16))
    json.dumps(report)
    status = report["status"]
    assert status["complete_history_occupation_limit_proved"]
    assert status["sharp_readout_rate"] == "not supplied"
    assert not status["spatial_continuum_or_empirical_result"]
    assert not status["post_sharp_readout_fourth_moment_preservation_claimed"]
    for item in report["pre_readout_bounds"]:
        assert not item["includes_final_sharp_field_readout"]
        assert not item["includes_source_splitting_or_floating_execution"]


@pytest.mark.parametrize("d,time", ((True, F(1)), (3, F(1)), (0, F(1)),
                                    (2.0, F(1)), (4, True), (4, .5),
                                    (4, F(-1)), (4, F(2))))
def test_rejects_undeclared_dimension_and_model_time(d, time):
    with pytest.raises(ValueError):
        r.occupation_evolution_bound(d, time)
