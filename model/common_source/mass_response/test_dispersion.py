"""Independent spectral and refinement checks for the dimer continuation."""
from itertools import product
import math

import numpy as np
import pytest

from . import dispersion as d


def test_clifford_generators_retain_eight_components():
    gamma = d.gamma_matrices()
    assert len(gamma) == 6
    for i, a in enumerate(gamma):
        np.testing.assert_array_equal(a, a.conj().T)
        assert a.shape == (8, 8)
        for j, b in enumerate(gamma):
            np.testing.assert_array_equal(a @ b + b @ a,
                                          2 * np.eye(8) if i == j else np.zeros((8, 8)))


@pytest.mark.parametrize("k", [(0., 0., 0.), (.2, -.4, .8),
                              (math.pi, .3, -1.2), (2.1, -2.5, .7)])
def test_free_bloch_matrix_has_the_full_dispersion(k):
    h = d.bloch_hamiltonian(k)
    expected = d.dispersion_squared(k)
    np.testing.assert_allclose(h @ h, expected * np.eye(8), atol=2e-14, rtol=1e-14)
    values = np.linalg.eigvalsh(h)
    energy = math.sqrt(expected)
    np.testing.assert_allclose(values, [-energy] * 4 + [energy] * 4, atol=3e-14)


@pytest.mark.parametrize("side", [1, 2, 3, 4, 5])
def test_structural_gap_is_not_a_finite_volume_artifact(side):
    minimum = min(d.dispersion_squared(tuple(2 * math.pi * n / side for n in point))
                  for point in product(range(side), repeat=3))
    assert minimum == pytest.approx(d.finite_torus_gap(side) ** 2, abs=2e-15)


def test_factorized_reference_control_has_zero_structural_gap_at_fixed_clock():
    omega = d.SOURCE_OMEGA
    for level in [1, 2, 8]:
        a, b = d.refinement_coefficients(level, omega=omega, imbalance=0.)
        assert a == b
        assert d.finite_torus_gap(3, intra=a, inter=b) == 0.
        assert d.refined_dispersion_squared((0., 0., 0.), level,
                                            omega=omega, imbalance=0.) == 0.
    assert omega == pytest.approx(math.log(5), abs=3e-16)


@pytest.mark.parametrize("p", [(.1, .2, -.3), (1., -.5, .7), (2., 1.2, -.4)])
def test_refinement_error_has_the_claimed_sign_and_uniform_bound(p):
    target = d.continuum_dispersion_squared(p)
    errors = []
    for level in [1, 2, 4, 8, 16, 32]:
        a, b = d.refinement_coefficients(level)
        k = tuple(x / level for x in p)
        h = d.bloch_hamiltonian(k, intra=a, inter=b)
        actual = d.refined_dispersion_squared(p, level)
        np.testing.assert_allclose(h @ h, actual * np.eye(8), atol=5e-13, rtol=2e-13)
        error = target - actual
        assert error >= -2e-13
        assert error <= d.refinement_remainder_bound(p, level) + 2e-13
        errors.append(error)
    assert errors[-1] < errors[0] / 500


def test_interaction_dresses_finite_gap_and_neutral_gap_is_a_different_observable():
    free = d.ring_gap_predictions(coupling=0.)
    interacting = d.ring_gap_predictions()
    assert free["charged_addition_gap"] == pytest.approx(d.SOURCE_J)
    assert free["charge_transfer_pair_gap"] == pytest.approx(2 * d.SOURCE_J)
    assert interacting["charged_addition_gap"] == pytest.approx(1.028595953717031, abs=2e-14)
    assert interacting["single_ring_neutral_gap"] == pytest.approx(1.305496913774672, abs=2e-14)
    assert interacting["charge_transfer_pair_gap"] > free["charge_transfer_pair_gap"]
    assert interacting["charge_transfer_pair_gap"] != pytest.approx(interacting["single_ring_neutral_gap"])


@pytest.mark.parametrize("a,b,g", [(d.SOURCE_A, d.SOURCE_B, d.SOURCE_J),
                                  (1.8, .2, 0.), (1.8, .2, 5.)])
def test_analytic_ring_matches_independent_half_filled_occupation_matrix(a, b, g):
    spectra = d.ring_sector_energies(intra=a, inter=b, coupling=g)
    assert [len(spectra[n]) for n in range(5)] == [1, 4, 6, 4, 1]
    # Direct occupation basis 3,5,6,9,10,12. Boundary-hop CAR signs are
    # retained; this matrix is independent of the analytic root formulas.
    h = np.array(((0, -b, 0, 0, b, 0),
                  (-b, -g, -a, -a, 0, b),
                  (0, -a, 0, 0, -a, 0),
                  (0, -a, 0, 0, -a, 0),
                  (b, 0, -a, -a, -g, -b),
                  (0, b, 0, 0, -b, 0)))
    np.testing.assert_allclose(np.linalg.eigvalsh(h), spectra[2], atol=3e-14)
    assert spectra[2][0] < min(spectra[n][0] for n in [0, 1, 3, 4])


@pytest.mark.parametrize("mu", [-3., -.25, 0., .7, 5.])
def test_same_total_charge_pair_frequency_cancels_common_chemical_potential(mu):
    spectra = d.ring_sector_energies(chemical_potential=mu)
    transfer = spectra[3][0] + spectra[1][0] - 2 * spectra[2][0]
    assert transfer == pytest.approx(d.ring_gap_predictions()["charge_transfer_pair_gap"], abs=5e-15)


@pytest.mark.parametrize("bad", [0, -1, True, 1.2])
def test_refinement_and_torus_reject_nonpositive_noninteger_levels(bad):
    with pytest.raises(ValueError):
        d.refinement_coefficients(bad)
    with pytest.raises(ValueError):
        d.finite_torus_gap(bad)


@pytest.mark.parametrize("bad", [(1., 2.), (1., 2., float("nan")),
                                (1., 2., True), (1., 2., 1j)])
def test_momentum_requires_three_finite_real_coordinates(bad):
    with pytest.raises(ValueError):
        d.bloch_hamiltonian(bad)


def test_ring_branch_and_refinement_require_their_displayed_domains():
    with pytest.raises(ValueError):
        d.ring_gap_predictions(intra=.5, inter=1.)
    with pytest.raises(ValueError):
        d.ring_gap_predictions(coupling=-.1)
    with pytest.raises(ValueError):
        d.refinement_coefficients(1, omega=1., imbalance=1.)
