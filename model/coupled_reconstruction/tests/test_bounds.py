"""Independent checks of the fixed real-product derivative certificate."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from model.coupled_reconstruction.bounds import (
    product_coefficient_radius, product_generator_cap_certificate,
)
from model.interacting_reconstruction.bounds import coefficient_radius


def _coherence(factors, time, j=0, k=1):
    state = np.eye(len(factors[0]), dtype=complex)[:, j]
    for generator in factors:
        state = expm(-1j*time*generator)@state
    return -float(np.imag(np.conj(state[j])*state[k]))


def _estimate(factors, u, v, j=0, k=1):
    return (v*v*_coherence(factors, u, j, k)/u
            - u*u*_coherence(factors, v, j, k)/v)/(v*v-u*u)


@pytest.mark.parametrize("dimension,initial,other", [(2, 0, 1), (4, 1, 3), (8, 0, 5)])
def test_noncommuting_real_products_obey_derivative_bound(dimension, initial, other):
    random = np.random.default_rng(710+dimension)
    factors = []
    cap = 0.
    for _ in range(4):
        raw = random.normal(size=(dimension, dimension))/4
        real = raw+raw.T
        spectrum = np.linalg.eigvalsh(real)
        centered = real-np.eye(dimension)*(spectrum[0]+spectrum[-1])/2
        cap += np.linalg.norm(centered, 2)+1e-12
        factors.append(real)
    u, v = .07, .19
    estimate = _estimate(factors, u, v, initial, other)
    target = sum(factors)[other, initial]
    certificate = product_coefficient_radius(u, v, cap)
    assert abs(estimate-target) <= certificate["radius"]
    assert _coherence(factors, -v, initial, other) == pytest.approx(
        -_coherence(factors, v, initial, other), abs=1e-14)
    assert certificate["continuous_semigroup_assumed"] is False


def test_global_generator_cancellation_does_not_bound_product_derivatives():
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    factors = [x, z, -(x+z)]
    assert np.array_equal(sum(factors), np.zeros((2, 2)))
    estimate = _estimate(factors, .1, .2)
    # The global Hamiltonian is zero, but this product is not its semigroup.
    # Using its norm would falsely give a zero derivative-estimation radius.
    assert abs(estimate) > 1e-4
    assert coefficient_radius(.2, 0) == 0
    sum_cap = sum(np.linalg.norm(h, 2) for h in factors)+1e-12
    assert abs(estimate) <= product_coefficient_radius(.1, .2, sum_cap)["radius"]


def test_centering_individual_generators_removes_only_global_phase():
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    factors = [x, z, x/2]
    shifted = [h+c*np.eye(2) for h, c in zip(factors, (100, -70, 21))]
    assert _estimate(factors, .05, .11) == pytest.approx(
        _estimate(shifted, .05, .11), abs=2e-13)


def test_changing_step_count_between_durations_invalidates_cancellation():
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    one_step = [x/2, z, x/2]
    two_steps = [x/4, z/2, x/4]*2
    u, v = .01, .02
    # Both schedules have the same first derivative and the same norm sum,
    # but their cubic coefficients differ.  The fixed-family bound cannot
    # be applied to this mixed pair of measurements.
    mixed = (v*v*_coherence(two_steps, u)/u
             - u*u*_coherence(one_step, v)/v)/(v*v-u*u)
    certificate = product_coefficient_radius(u, v, 2)
    assert abs(mixed-1) > 50*certificate["radius"]
    assert abs(_estimate(one_step, u, v)-1) <= certificate["radius"]
    assert abs(_estimate(two_steps, u, v)-1) <= certificate["radius"]


def test_product_radius_matches_exact_rational_dyadic_formula():
    t = F(1, 8192)
    result = product_coefficient_radius(t/2, t)
    assert F(result["finite_time_radius_exact"]) == 96**5*t**4/18
    assert result["radius"] == coefficient_radius(t, 96)
    assert result["fixed_scaled_factor_sequence_required"] is True
    assert result["independent_splitting_error_added"] is False


def test_probability_uncertainty_is_separate_from_product_remainder():
    t, error = F(1, 8192), F(1, 10**9)
    result = product_coefficient_radius(t/2, t, 96, error)
    assert F(result["probability_radius_exact"]) == 3*error/t
    assert result["clock_floating_and_physical_service_error_included"] is False


def test_default_backend_cap_is_an_independent_exact_audit():
    certificate = product_generator_cap_certificate()
    exact = F(certificate["exact_schedule_centered_norm_sum_upper"])
    assert 86 < exact < 87
    assert certificate["declared_cap_covers_source_schedule"] is True
    assert F(certificate["schedule_centered_norm_sum_upper"]) >= exact
    assert certificate["backend_audit_only"] is True
    assert product_generator_cap_certificate(declared_cap=80)[
        "declared_cap_covers_source_schedule"] is False


@pytest.mark.parametrize("first,second,cap", [(0, 1, 96), (1, 1, 96), (2, 1, 96),
                                              (1, 2, -1), (1, 2, np.nan)])
def test_invalid_product_bounds(first, second, cap):
    with pytest.raises(ValueError):
        product_coefficient_radius(first, second, cap)
