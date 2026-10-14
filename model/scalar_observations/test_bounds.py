"""Independent operator and Gaussian-form checks of the observation bounds."""
from fractions import Fraction as F
from math import pi

import numpy as np
import pytest
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import expm

from . import bounds
from .lattice import GalerkinParameters, SevenModeLattice, stiffness
from ..relativistic_scalar.lattice import local_operators


def _hermitian_norm(matrix):
    return float(np.max(np.abs(np.linalg.eigvalsh(matrix))))


def _split_matrices(model):
    a = sum(model._tensor_product({mode: model.onsite(mode)}) for mode in range(7)).toarray()
    return a, model.hamiltonian.toarray() - a


@pytest.mark.parametrize("dimension", [2, 4, 8])
def test_full_coordinate_square_norm_against_independent_fock_matrix(dimension):
    q2 = local_operators(dimension, 2)["q2"]
    assert _hermitian_norm(q2) <= float(bounds.coordinate_square_norm_upper(dimension)) + 2e-14
    # The full compression includes the boundary excursion absent from q_d^2.
    q = local_operators(dimension, 2)["q"]
    assert np.linalg.norm(q2 - q @ q) > .4


@pytest.mark.parametrize("dimension", [2, 4])
@pytest.mark.parametrize("grid", [5, 9, 33, None])
def test_onsite_bounds_cover_actual_family(dimension, grid):
    model = SevenModeLattice((dimension,) * 7, GalerkinParameters(grid=grid))
    for mode in (0, 1):
        onsite = model.onsite(mode)
        centered = onsite - onsite[0, 0] * np.eye(dimension)
        assert _hermitian_norm(centered) <= float(bounds.onsite_norm_upper(dimension, mode)) + 2e-14


@pytest.mark.parametrize("grid", [5, 17, None])
def test_global_nested_commutators_against_independent_matrix_products(grid):
    # One four-level mode retains noncommuting kinetic/potential terms;
    # the other six modes still contribute their interaction branches.
    dimensions = (2, 4, 2, 2, 2, 2, 2)
    model = SevenModeLattice(dimensions, GalerkinParameters(grid=grid))
    a, b = _split_matrices(model)
    ab = a @ b - b @ a
    aab = a @ ab - ab @ a
    bab = b @ ab - ab @ b
    certificate = bounds.commutator_bounds(dimensions)
    assert _hermitian_norm(aab) <= float(certificate["A_A_B_norm_upper"]) + 2e-12
    assert _hermitian_norm(bab) <= float(certificate["B_A_B_norm_upper"]) + 2e-12
    assert _hermitian_norm(aab) > .1


def test_product_bound_covers_actual_noncommuting_unitary():
    dimensions = (2, 4, 2, 2, 2, 2, 2)
    model = SevenModeLattice(dimensions, GalerkinParameters(grid=5))
    a, b = _split_matrices(model)
    time, steps = F(1, 8), 2
    half = expm(-.5j * float(time) / steps * a)
    step = half @ expm(-1j * float(time) / steps * b) @ half
    product = np.linalg.matrix_power(step, steps)
    actual = float(np.linalg.norm(product - expm(-1j * float(time) * (a + b)), 2))
    assert actual > 1e-9
    assert actual <= float(bounds.product_error(dimensions, time, steps))


def test_two_level_degeneracy_is_exact_splitting_not_an_occupation_claim():
    assert bounds.commutator_bounds(2) == {"A_A_B_norm_upper": F(0), "B_A_B_norm_upper": F(0)}
    assert bounds.splitting_history_bound(2, (F(1, 8),) * 2, 1) == 0
    assert bounds.spatial_history_bound(2, 5, F(1, 4)) > 0


@pytest.mark.parametrize("grid", [5, 9, 17, 33, 1000])
def test_rational_spatial_envelope_covers_actual_stiffness(grid):
    assert 0 < 1 - stiffness(grid) < float(bounds.kappa_error_upper(grid))


def test_spatial_bound_covers_independent_operator_difference():
    dimensions = (2, 4, 2, 2, 2, 2, 2)
    finite = SevenModeLattice(dimensions, GalerkinParameters(grid=5))
    target = SevenModeLattice(dimensions)
    delta = finite.hamiltonian.toarray() - target.hamiltonian.toarray()
    time = F(1, 4)
    assert float(time) * _hermitian_norm(delta) <= float(bounds.spatial_history_bound(dimensions, 5, time))
    assert bounds.spatial_history_bound(dimensions, None, time) == 0


def test_detector_energy_increment_by_independent_gaussian_quadrature():
    # Integrate actual derivatives of K_b(q) exp(-Omega*q^2/2+i*k*q).
    # Potential contributions cancel pointwise when both branches are kept.
    nodes, weights = hermgauss(64)
    omega, eta, theta, kick = 2., 1., pi / 2, .25
    q = nodes / np.sqrt(omega)
    weights = weights / np.sqrt(pi)
    angle = (eta * q + theta) / 2
    maps = (np.cos(angle), np.sin(angle))
    derivatives = (-eta * np.sin(angle) / 2, eta * np.cos(angle) / 2)
    before = omega / 4 + kick * kick / 2
    after = sum(.5 * np.dot(weights, np.abs(derivative + field * (-omega * q + 1j * kick))**2)
                for field, derivative in zip(maps, derivatives, strict=True))
    assert after - before == pytest.approx(eta * eta / 8, abs=2e-14)
    assert np.dot(weights, maps[0]**2 + maps[1]**2) == pytest.approx(1., abs=2e-14)


def test_blank_energy_and_conditional_feedback_envelope():
    # Full polynomial compression retains the unrestricted blank moments.
    model = SevenModeLattice((2,) * 7)
    assert model.hamiltonian[0, 0] == pytest.approx(float(F(705, 128)), abs=1e-14)
    assert bounds.blank_energy(initial_kick=1) == F(769, 128)
    principal = bounds.unrestricted_spatial_history(
        33, (F(1, 8),) * 2, eta=1, initial_kick=F(1, 4), feedback_kick=F(1, 4))
    assert principal["initial_energy_exact"] == F(709, 128)
    assert principal["after_first_read_energy_exact"] == F(725, 128)
    # theta=pi/2 gives each first bit probability 1/2 in the initial
    # Gaussian.  For an immediate detector, its real Kraus functions keep
    # conditional mean momentum k0; one branch receives the feedback kick.
    k0, k = F(1, 4), F(1, 4)
    immediate_feedback_energy = F(725, 128) + (k * k0 + k * k / 2) / 2
    assert immediate_feedback_energy <= principal["after_conditional_feedback_energy_upper"]


def test_unrestricted_fixed_band_rate_and_scope():
    settings = dict(times=(F(1, 8),) * 2, eta=1,
                    initial_kick=F(1, 4), feedback_kick=F(1, 4))
    coarse = bounds.unrestricted_spatial_history(9, **settings)
    fine = bounds.unrestricted_spatial_history(18, **settings)
    assert coarse["complete_history_tv_upper"] == 4 * fine["complete_history_tv_upper"]
    assert not fine["finite_occupation_error_included"]
    assert not fine["mode_band_removal_included"]
    assert fine["every_record_and_disturbed_state_retained"]
    assert bounds.unrestricted_spatial_history(None, **settings)["complete_history_tv_upper"] == 0
    with pytest.raises(ValueError, match="coupling"):
        bounds.unrestricted_spatial_history(9, coupling=0, **settings)
    # The finite-code comparison also covers the interaction-deletion control.
    assert bounds.generator_sum_bound(4, coupling=0) > 0
    assert bounds.product_error(4, F(1, 8), 2, coupling=0) == 0


@pytest.mark.parametrize("bad", [True, 4, 5.0, -1])
def test_invalid_spatial_grid_is_rejected(bad):
    with pytest.raises(ValueError):
        bounds.kappa_error_upper(bad)


def test_certificates_reject_float_times_and_verify_root_direction():
    with pytest.raises(ValueError):
        bounds.product_error(4, .125, 32)
    with pytest.raises(ValueError):
        bounds.unrestricted_spatial_history(9, (F(1, 8),))
    with pytest.raises(ValueError):
        bounds.product_error(4, F(1, 8), True)
    for value in (F(0), F(1, 2), F(2), F(81, 16), F(10**90, 7)):
        upper = bounds.sqrt_upper(value)
        assert upper * upper >= value
        if upper:
            assert (upper - F(1, 10**36))**2 < value or upper * upper == value
