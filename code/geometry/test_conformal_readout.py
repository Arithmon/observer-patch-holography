"""Independent analytic, exact algebraic and adversarial geometry checks."""

import itertools
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import conformal_readout as cr
from quotient_cap_readout import kms_receipt


def circle(alpha, n=12):
    t = np.arange(n) * 2 * np.pi / n
    return np.column_stack((np.sin(alpha) * np.cos(t), np.sin(alpha) * np.sin(t),
                            np.full(n, np.cos(alpha))))


def inverse_stereo(z):
    if np.isinf(z):
        return np.array([0., 0., 1.])
    r = abs(z) ** 2
    return np.array([2 * z.real, 2 * z.imag, r - 1]) / (1 + r)


@pytest.mark.parametrize('alpha', [0.03, 0.4, np.pi / 2, 2.4, np.pi - 0.03])
@pytest.mark.parametrize('n', [3, 5, 30])
def test_small_great_and_major_caps_with_explicit_orientation(alpha, n):
    points = circle(alpha, n)
    expected = np.array([np.cos(alpha), 0, 0, 1]) / np.sin(alpha)
    fit = cr.fit_cap(points, interior_point=[0, 0, 1])
    np.testing.assert_allclose(fit.normal, expected, atol=2e-11)
    opposite = cr.produced_cap_normal(points, interior_point=[0, 0, -1])
    np.testing.assert_allclose(opposite, -expected, atol=2e-11)
    assert abs(-fit.normal[0] ** 2 + np.dot(fit.normal[1:], fit.normal[1:]) - 1) < 2e-11
    assert fit.max_incidence_residual <= 1e-10


def test_great_circle_rational_nullspace_is_spacelike_and_unoriented():
    rows = sp.Matrix([[-1, 1, 0, 0], [-1, 0, 1, 0], [-1, -1, 0, 0]])
    assert rows.rank() == 3
    (normal,) = rows.nullspace()
    eta = sp.diag(-1, 1, 1, 1)
    assert (normal.T * eta * normal)[0] == 1
    assert rows * normal == rows * (-normal) == sp.zeros(3, 1)
    points = np.array(rows[:, 1:]).astype(float)
    np.testing.assert_allclose(cr.produced_cap_normal(points, interior_point=[0, 0, 1]),
                               [0, 0, 0, 1], atol=1e-14)
    with pytest.raises(TypeError):
        cr.produced_cap_normal(points)  # no sign chosen by the implementation


@pytest.mark.parametrize('points', [[], [[1, 0, 0]], [[1, 0, 0]] * 3,
                                   [[1, 0, 0], [-1, 0, 0], [1, 0, 0]],
                                   [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                                   [[np.nan, 0, 1]] * 3,
                                   [[np.inf, 0, 1]] * 3,
                                   [[True, False, False]] * 3,
                                   [[1j, 0, 0]] * 3])
def test_cap_rejects_invalid_or_rank_deficient_data(points):
    with pytest.raises(ValueError):
        cr.fit_cap(points, interior_point=[0, 0, 1])


@pytest.mark.parametrize('budget', [True, np.nan, np.inf, -1, 1, 1j, '0.01'])
def test_residual_budget_is_validated(budget):
    with pytest.raises(ValueError):
        cr.fit_cap(circle(0.7), interior_point=[0, 0, 1], max_residual=budget)


def test_unsupported_noisy_circle_and_ambiguous_side_fail():
    points = circle(0.7)
    points[0] = [0, 0, -1]
    with pytest.raises(ValueError, match='residual'):
        cr.fit_cap(points, interior_point=[0, 0, 1])
    with pytest.raises(ValueError, match='side witness'):
        cr.fit_cap(circle(np.pi / 2), interior_point=[1, 0, 0])
    with pytest.raises(ValueError, match='unit sphere'):
        cr.fit_cap(circle(0.7), interior_point=[0, 0, 2])
    # A near-point circle has no resolved spacelike normalization.
    with pytest.raises(ValueError):
        cr.fit_cap(circle(1e-9), interior_point=[0, 0, 1])


def test_noise_budget_and_independent_direction_bound():
    points = circle(0.7)
    rng = np.random.default_rng(445)
    perturbed = points + 1e-4 * rng.normal(size=points.shape)
    perturbed /= np.linalg.norm(perturbed, axis=1)[:, None]
    fit = cr.fit_cap(perturbed, interior_point=[0, 0, 1], max_residual=0.001)
    with pytest.raises(ValueError, match='residual'):
        cr.fit_cap(perturbed, interior_point=[0, 0, 1])
    true_direction = np.array([np.cos(0.7), 0, 0, 1])
    true_direction /= np.linalg.norm(true_direction)
    computed = fit.normal / np.linalg.norm(fit.normal)
    sine = np.linalg.norm(computed - np.dot(computed, true_direction) * true_direction)
    error = np.linalg.norm(perturbed - points, ord=2)
    assert sine <= fit.direction_error_bound(error) < 0.01
    for invalid in (np.nan, np.inf, -1, True, fit.third_singular_value):
        with pytest.raises(ValueError):
            fit.direction_error_bound(invalid)


def test_infinity_safe_cross_ratios_all_four_positions():
    z = sp.symbols('z', real=True)
    for position in range(4):
        values = [sp.Rational(1, 3), sp.Rational(2, 5), sp.Rational(3, 2), sp.Rational(-4, 3)]
        values[position] = z
        a, b, c, d = values
        expected = sp.limit((a-c)*(b-d)/((a-d)*(b-c)), z, sp.oo)
        numeric = [float(v) if v != z else complex(np.inf, 0) for v in values]
        assert abs(cr.cross_ratio(*numeric) - float(expected)) < 1e-13


@pytest.mark.parametrize('bad', [True, '1', np.nan, complex(1, np.inf), -np.inf, [1]])
def test_cross_ratio_rejects_malformed_scalar(bad):
    with pytest.raises(ValueError):
        cr.cross_ratio(bad, 2, 3, 4)


def test_near_pole_and_near_coincidence_are_not_rounded_to_exact_data():
    with pytest.raises(ValueError):
        cr.cross_ratio(1e308, -1e308, 1, 2)
    with np.errstate(over='ignore', invalid='ignore'):
        with pytest.raises(ValueError):
            cr.stereographic([1e-310, 0, 1])
    assert cr.stereographic([0, 0, 1]) == complex(np.inf, 0)
    with pytest.raises(ValueError):
        cr.cross_ratio(1, 2, 1+1e-14, 3)


def test_exact_mobius_determinant_identity():
    a, b, c, d, x, y = sp.symbols('a b c d x y')
    assert sp.expand((a*x+b)*(c*y+d)-(a*y+b)*(c*x+d) - (a*d-b*c)*(x-y)) == 0


def test_poles_and_gauge_changes_reconstruct_from_scalar_data_only(monkeypatch):
    points = np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0], [0, 1, 0], [-1, 0, 0]])
    for gauge in itertools.permutations(range(3)):
        receipts = cr.cross_ratio_receipts(points, gauge)
        reconstructed = cr.reconstruct_from_cross_ratios(receipts, gauge)
        np.testing.assert_array_equal(reconstructed, receipts)
        assert not np.isnan(reconstructed).any()
    data = cr.cross_ratio_receipts(points, (0, 1, 2))
    def forbidden(*args, **kwargs):
        raise AssertionError('receipt consumer accessed source coordinates')
    monkeypatch.setattr(cr, 'stereographic', forbidden)
    monkeypatch.setattr(cr, 'sphere_points', forbidden)
    np.testing.assert_array_equal(cr.reconstruct_from_cross_ratios(data, (0, 1, 2)), data)


@pytest.mark.parametrize('alpha', [0.4, np.pi/2, 2.4])
def test_complete_scalar_receipt_to_oriented_cap_pipeline(alpha, monkeypatch):
    source = np.vstack((circle(alpha, 6), [0, 0, 1], [0, 0, -1]))
    gauge = (0, 1, 2)
    data = cr.cross_ratio_receipts(source, gauge)
    def forbidden(*args, **kwargs):
        raise AssertionError('consumer tried to generate its own source data')
    monkeypatch.setattr(cr, 'cross_ratio_receipts', forbidden)
    monkeypatch.setattr(cr, 'stereographic', forbidden)
    fit = cr.cap_from_cross_ratios(data, gauge, boundary_indices=list(range(6)), interior_index=6)
    complement = cr.cap_from_cross_ratios(data, gauge, boundary_indices=list(range(6)), interior_index=7)
    np.testing.assert_allclose(fit.normal, -complement.normal, atol=1e-13)
    # Three boundary gauge points are 0,1,infinity, so this representative is a great circle.
    assert abs(fit.normal[0]) < 1e-13
    actual_points = cr.sphere_from_cross_ratios(data, gauge)
    independent_points = np.array([inverse_stereo(z) for z in data])
    np.testing.assert_allclose(actual_points, independent_points, atol=1e-14)
    assert np.dot(fit.normal[1:], actual_points[6]) - fit.normal[0] > 0
    assert np.dot(fit.normal[1:], actual_points[7]) - fit.normal[0] < 0


@pytest.mark.parametrize('boundary,inside', [([0,1],3), ([0,0,2],3), ([0,1,8],3),
                                            ([0,1,2],True), ([0,1,2],1),
                                            ([0,1,2],9), ([False,1,2],3)])
def test_cap_receipt_pipeline_rejects_wrong_semantic_labels(boundary, inside):
    with pytest.raises(ValueError):
        cr.cap_from_cross_ratios([0,1,np.inf,1j], (0,1,2),
                                 boundary_indices=boundary, interior_index=inside)


@pytest.mark.parametrize('data,gauge', [
    ([0, 1, np.inf, np.nan], (0, 1, 2)),
    ([0, 1, np.inf, complex(0, np.inf)], (0, 1, 2)),
    ([0, 1, np.inf, -np.inf], (0, 1, 2)),
    ([0, 1, np.inf, 0], (0, 1, 2)),
    ([0, 1, np.inf, 1e-20], (0, 1, 2)),
    ([0, 1, np.inf], (0, 0, 2)),
    ([0, 1, np.inf], (0, 1, True)),
    ([0, 1, np.inf], (0, 1, 3)),
    ([0, 2, np.inf], (0, 1, 2)),
    (np.eye(3), (0, 1, 2)),
    (['0', '1', 'inf'], (0, 1, 2)),
])
def test_receipt_consumer_rejects_malformed_inputs(data, gauge):
    with pytest.raises(ValueError):
        cr.reconstruct_from_cross_ratios(data, gauge)


def test_lorentz_boost_covariance_and_conformal_metric_boundary():
    alpha, rapidity = 1.0, 0.7
    points = circle(alpha)
    c, s = np.cosh(rapidity), np.sinh(rapidity)
    boost = np.array([[c, 0, 0, s], [0, 1, 0, 0], [0, 0, 1, 0], [s, 0, 0, c]])
    eta = np.diag([-1., 1, 1, 1])
    np.testing.assert_allclose(boost.T @ eta @ boost, eta, atol=1e-14)
    def transform(points):
        rays = np.column_stack((np.ones(len(points)), points)) @ boost.T
        return rays[:, 1:] / rays[:, :1]
    transformed = transform(points)
    witness = transform(np.array([[0., 0, 1]]))[0]
    n = cr.produced_cap_normal(points, interior_point=[0, 0, 1])
    n2 = cr.produced_cap_normal(transformed, interior_point=witness)
    np.testing.assert_allclose(n2, boost @ n, atol=1e-13)
    before = cr.cross_ratio_receipts(points, (0, 1, 2))
    after = cr.cross_ratio_receipts(transformed, (0, 1, 2))
    finite = np.isfinite(before)
    np.testing.assert_allclose(after[finite], before[finite], atol=1e-12)
    # A boost changes angular size while preserving the complete supplied ratios.
    alpha2 = np.arccos(n2[0] / np.linalg.norm(n2[1:]))
    assert abs(alpha2 - alpha) > 0.1


def test_rotation_covariance_and_input_permutation():
    rng = np.random.default_rng(919)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    for alpha in (0.2, np.pi / 2, 2.5):
        points = circle(alpha) @ rotation.T
        witness = rotation[:, 2]
        expected = np.r_[np.cos(alpha), witness] / np.sin(alpha)
        for _ in range(4):
            points = points[rng.permutation(len(points))]
            np.testing.assert_allclose(cr.produced_cap_normal(points, interior_point=witness),
                                       expected, atol=1e-12)


def test_guards_survive_optimized_python():
    script = '''
from finite_incidence import incidence_complex
from quotient_cap_readout import spherical_incidence_receipt
from conformal_readout import reconstruct_from_cross_ratios
if spherical_incidence_receipt(incidence_complex([(0,1,2,3)])):
    raise RuntimeError('solid simplex accepted as sphere')
try:
    reconstruct_from_cross_ratios([0,1,float('inf'),float('nan')], (0,1,2))
except ValueError:
    pass
else:
    raise RuntimeError('NaN receipt accepted')
'''
    subprocess.run([sys.executable, '-O', '-c', script], cwd=Path(__file__).resolve().parent,
                   check=True, capture_output=True, text=True)


@pytest.mark.parametrize('kwargs', [dict(clock_scale=np.nan), dict(clock_scale=np.inf),
                                    dict(clock_scale=True), dict(clock_scale=-1),
                                    dict(clock_scale=1, tol=np.inf),
                                    dict(clock_scale=1, tol=100),
                                    dict(clock_scale=1, beta_target=1j)])
def test_scalar_clock_check_rejects_invalid_comparisons(kwargs):
    with pytest.raises(ValueError):
        kms_receipt(**kwargs)


def test_scalar_kms_agreement_does_not_supply_a_physical_clock():
    for supplied_beta in (1., 2., 2*np.pi, 10.):
        assert kms_receipt(supplied_beta, beta_target=supplied_beta)
    assert not kms_receipt(1.)
