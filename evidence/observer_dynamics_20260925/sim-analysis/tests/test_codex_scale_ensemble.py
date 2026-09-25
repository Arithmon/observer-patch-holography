"""Independent small-data checks of the new scale/ensemble estimator.

These tests do not read or modify simulation data. In particular, the stationary
null is checked by exhaustive configurations and harmonics by direct quadrature.
"""

import importlib.util
import itertools
from pathlib import Path

import numpy as np
import pytest
from scipy.special import sph_harm_y


SPEC = importlib.util.spec_from_file_location(
    "codex_scale_ensemble", Path(__file__).resolve().parents[1] / "scripts" / "codex_scale_ensemble.py"
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_weighted_regression_has_intercept_and_exact_orthogonal_decomposition():
    x = np.array([-2., 1., 4.])
    weights = np.array([1., 2., 3.])
    # This residual has weighted mean zero and weighted covariance zero with x.
    orthogonal = np.array([3., -3., 1.])
    y = 2 * x + 7 + orthogonal
    result = AUDIT.decomposition(x, y, weights)
    assert result["response"] == pytest.approx(2.)
    assert result["initial_variance"] == pytest.approx(5.)
    assert result["terminal_variance"] == pytest.approx(25.)
    assert result["power_ratio"] == pytest.approx(5.)
    assert result["residual_power_ratio"] == pytest.approx(1.)
    assert result["coherence_squared"] == pytest.approx(.8)
    shifted = AUDIT.decomposition(x + 1000, y - 300, weights * 17)
    assert shifted == pytest.approx(result)
    # An unweighted calculation is different; this catches lost weights.
    assert AUDIT.decomposition(x, y, np.ones(3))["response"] != pytest.approx(2.)


@pytest.mark.parametrize("groups,weights", [
    ([2, 3, 3], [1., 1., 1.]),
    ([2, 3, 3], [2., 3., 3.]),
    ([1, 3], [4., 1.]),  # Disjoint groups need not cover the whole population.
])
def test_fixed_total_exchangeable_null_against_complete_enumeration(groups, weights):
    total, raised = 8, 3
    groups, weights = np.array(groups), np.array(weights)
    offsets = np.concatenate(([0], np.cumsum(groups)))
    variances = []
    for occupied in itertools.combinations(range(total), raised):
        field = np.zeros(total)
        field[list(occupied)] = 1.
        means = np.array([field[a:b].mean() for a, b in zip(offsets[:-1], offsets[1:])])
        mean = sum(weights * means) / sum(weights)
        variances.append(sum(weights * (means - mean) ** 2) / sum(weights))
    actual = AUDIT.exchangeable_null_variance(groups, weights, total, raised / total)
    assert actual == pytest.approx(np.mean(variances), rel=1e-13)


def test_iid_null_against_complete_independent_bit_enumeration():
    groups, weights = np.array([1, 2, 3]), np.array([1., 4., 2.])
    offsets = np.concatenate(([0], np.cumsum(groups)))
    variances = []
    for state in itertools.product((0., 1.), repeat=6):
        means = np.array([np.mean(state[a:b]) for a, b in zip(offsets[:-1], offsets[1:])])
        mean = np.average(means, weights=weights)
        variances.append(np.average((means - mean) ** 2, weights=weights))
    assert AUDIT.iid_null_variance(groups, weights, variance=.25) == pytest.approx(np.mean(variances))


def test_schedule_variance_and_cross_schedule_power_against_pairwise_definition():
    x, weights = np.array([-2., 1., 4.]), np.array([1., 2., 3.])
    residual = np.array([3., -3., 1.])
    fields = np.array([2 * x + 7 + residual, .5 * x - 3 - residual, -x + 9 + .3 * residual])
    scale = AUDIT.Scale("test", np.array([1, 2, 3]), x, weights, 72, .5)
    for seed, y in enumerate(fields):
        scale.add(y, seed, float(x.mean()), float(y.mean()))
    result = scale.result()
    means = np.average(fields, weights=weights, axis=1)
    centered = fields - means[:, None]
    a = weights / weights.sum()
    total = np.mean([np.sum(a * y * y) for y in centered])
    noise = np.sum(a * np.var(centered, axis=0, ddof=1))
    common = np.mean([np.sum(a * y * z) for i, y in enumerate(centered) for j, z in enumerate(centered) if i != j])
    assert result["total_terminal_variance"] == pytest.approx(total)
    assert result["conditional_schedule_variance_unbiased"] == pytest.approx(noise)
    assert result["cross_schedule_common_power_unbiased"] == pytest.approx(common)
    assert total == pytest.approx(noise + common)
    assert result["response_schedule_se"] == pytest.approx(result["response_schedule_sd"] / np.sqrt(3))


def test_negative_unbiased_common_power_is_preserved_not_clamped():
    x = np.array([-1., 0., 2.])
    scale = AUDIT.Scale("opposites", np.ones(3, dtype=int), x, np.ones(3), 36, .5)
    scale.add(x, 0, x.mean(), x.mean())
    scale.add(-x, 1, x.mean(), -x.mean())
    result = scale.result()
    assert result["cross_schedule_common_power_unbiased"] == pytest.approx(-x.var())
    assert result["conditional_schedule_variance_unbiased"] == pytest.approx(2 * x.var())
    assert result["schedule_noise_fraction"] == pytest.approx(2.)


def test_fft_harmonics_match_full_direct_midpoint_quadrature():
    nz, nphi, lmax = 13, 28, 6
    values = np.random.default_rng(419).normal(size=(nz, nphi)) + 12.
    theta = np.arccos(-1 + (np.arange(nz) + .5) * 2 / nz)
    phi = (np.arange(nphi) + .5) * 2 * np.pi / nphi
    field = values - values.mean()
    coefficients = AUDIT.angular_coefficients(values.ravel(), nz, nphi, lmax)
    for ell in range(1, lmax + 1):
        for m in range(ell + 1):
            harmonic = sph_harm_y(ell, m, theta[:, None], phi[None, :])
            direct = 4 * np.pi / (nz * nphi) * np.sum(field * harmonic.conj())
            assert coefficients[ell - 1][m] == pytest.approx(direct, abs=1e-13)
    # Monopole removal must make coefficients insensitive to an added constant.
    shifted = AUDIT.angular_coefficients((values + 100).ravel(), nz, nphi, lmax)
    for a, b in zip(coefficients, shifted):
        assert np.allclose(a, b, atol=1e-13)


def test_angular_summary_matches_full_negative_m_real_field_coefficients():
    # A real spherical field has a_l,-m = (-1)^m conj(a_l,m).
    x = np.array([1.2 + 0j, .4 + .3j, -.2 + .5j])
    ys = np.array([2 * x + [1., .5j, -.3], .7 * x + [-.2, .2, .4j], -x + [1., -.1j, .2]])

    def full(a):
        return np.concatenate([np.array([(-1) ** m * a[m].conjugate() for m in (2, 1)]), a])

    # Supply a harmless ell=1 row so the tested coefficients are assigned ell=2.
    dummy = np.array([1. + 0j, .5 + .1j])
    result = AUDIT.angular_summary([[dummy, x]] + [[dummy, y] for y in ys])[1]
    xx = np.mean(abs(full(x)) ** 2)
    power = np.array([np.mean(abs(full(y)) ** 2) for y in ys])
    cross = np.array([np.mean((full(y) * full(x).conjugate()).real) for y in ys])
    noise = np.mean(np.var(np.array([full(y) for y in ys]), axis=0, ddof=1))
    assert result["C_initial"] == pytest.approx(xx)
    assert result["C_terminal_mean"] == pytest.approx(power.mean())
    assert result["per_schedule_response"] == pytest.approx(cross / xx)
    assert result["conditional_schedule_noise_power"] == pytest.approx(noise)
    assert result["cross_schedule_common_power"] == pytest.approx(power.mean() - noise)
    assert result["coherence_squared_mean"] == pytest.approx(np.mean(cross ** 2 / (xx * power)))


def test_angular_grid_midpoints_roundtrip_and_poles_do_not_overflow():
    nz, nphi = 4, 8
    z = -1 + (np.arange(nz) + .5) * 2 / nz
    phi = (np.arange(nphi) + .5) * 2 * np.pi / nphi
    zz, pp = np.meshgrid(z, phi, indexing="ij")
    radius = np.sqrt(1 - zz ** 2)
    directions = np.stack([radius * np.cos(pp), radius * np.sin(pp), zz], axis=-1).reshape(-1, 3)
    directions = np.concatenate([directions, [[0., 0., -1.], [0., 0., 1.]]])
    labels, counts = AUDIT.angular_partition(directions, nz, nphi)
    assert labels[:nz * nphi].tolist() == list(range(nz * nphi))
    assert counts.sum() == nz * nphi + 2
    assert np.all(counts >= 1)
    with pytest.raises(ValueError, match="empty pixels"):
        AUDIT.angular_partition(np.array([[0., 0., 1.]]), nz, nphi)
