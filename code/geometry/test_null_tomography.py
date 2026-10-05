"""Adversarial and independent controls for finite null tomography."""

import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from einstein_closure_receipts import (
    ETA, charges_of, generic_null_directions, null_vector,
    reconstruct_from_charges,
)
from geometry.null_tomography import (
    design_matrix, eta_project_out, fit_null_charges, tomography_directions,
    tracefree_basis,
)


@pytest.mark.parametrize("magnitude", [1e200, 1e-200])
def test_direction_units_cannot_turn_a_null_ray_timelike(magnitude):
    k = null_vector(np.array([magnitude, 0.0, 0.0]))
    np.testing.assert_array_equal(k, [1.0, 1.0, 0.0, 0.0])
    assert k @ ETA @ k == 0.0


@pytest.mark.parametrize("count", [1, 8, 12])
def test_zero_residual_does_not_certify_incomplete_tomography(count):
    rays = [np.array([1.0, 1.0, 0.0, 0.0])] * count
    hidden = np.diag([0.0, 0.0, 1.0, -1.0])
    np.testing.assert_array_equal(charges_of(hidden, rays), np.zeros(count))
    with pytest.raises(ValueError, match="rank|nine"):
        reconstruct_from_charges(np.zeros(count), rays)


def test_a_timelike_design_cannot_be_used_as_null_tomography():
    with pytest.raises(ValueError, match="null"):
        reconstruct_from_charges(np.zeros(12), [np.array([1., 0., 0., 0.])] * 12)


@pytest.mark.parametrize("bad", [np.zeros(3), [1., 0.], [np.nan, 1., 0.],
                                  [np.inf, 1., 0.], [1.+1j, 0., 0.]])
def test_malformed_directions_fail_closed(bad):
    with pytest.raises(ValueError):
        null_vector(bad)


def test_nonsymmetric_tensor_is_not_silently_projected_to_a_source():
    t = np.eye(4)
    t[0, 1] = 1.0
    with pytest.raises(ValueError, match="symmetric"):
        charges_of(t, generic_null_directions(12))


@pytest.fixture(scope="module")
def exact_frame():
    """Independent symbolic coordinate map, with no producer design calls."""
    x = sp.symbols("x:9")
    t = sp.Matrix([[x[0], x[1], x[2], x[3]],
                   [x[1], x[4], x[5], x[6]],
                   [x[2], x[5], x[7], x[8]],
                   [x[3], x[6], x[8], x[0] - x[4] - x[7]]])
    r = sp.sqrt(3) / 3
    directions = [sp.Matrix(k) for k in (
        (1, 1, 0, 0), (1, -1, 0, 0), (1, 0, 1, 0), (1, 0, -1, 0),
        (1, 0, 0, 1), (1, 0, 0, -1),
        (1, r, r, r), (1, r, r, -r), (1, r, -r, r))]
    design = sp.Matrix([(k.T * t * k)[0].expand() for k in directions]).jacobian(x)
    return x, t, directions, design


def test_exact_frame_reuses_existing_lean_theorem_and_recovers_every_coordinate(exact_frame):
    x, t, directions, a = exact_frame
    eta = sp.diag(-1, 1, 1, 1)
    assert all((k.T * eta * k)[0] == 0 for k in directions)
    assert sp.simplify(a.det()) == sp.Rational(8192, 27)
    # Exact symbolic inverse independently computed from tensor evaluation.
    inverse = a.inv()
    assert (inverse * a).applyfunc(sp.simplify) == sp.eye(9)
    rays = np.array([list(k) for k in directions], dtype=float)
    np.testing.assert_allclose(tomography_directions(), rays, rtol=0, atol=2e-16)
    for j in range(9):
        data = {c: int(i == j) for i, c in enumerate(x)}
        tensor = np.array(t.subs(data), dtype=float)
        exact_charges = np.array(a[:, j], dtype=float).ravel()
        fit = fit_null_charges(exact_charges, rays)
        np.testing.assert_allclose(fit.tensor, tensor, rtol=0, atol=8e-15)
        np.testing.assert_allclose(charges_of(tensor, rays), exact_charges,
                                   rtol=0, atol=1e-15)


def test_exact_cubic_countermodel_and_dependent_family_witness(exact_frame):
    x, t, directions, a = exact_frame
    # Degree-three spherical harmonic f(n)=n_x(n_y^2-n_z^2).
    # It vanishes on ALL nine frame directions, but is not a quadratic charge.
    extra = sp.Matrix([1, sp.Rational(2, 3), sp.Rational(1, 3), sp.Rational(2, 3)])
    directions = directions + [extra]
    cubic = sp.Matrix([sp.expand(k[1] * (k[2]**2 - k[3]**2)) for k in directions])
    assert cubic[:9, :] == sp.zeros(9, 1)
    assert cubic[9] == -sp.Rational(2, 9)
    row = sp.Matrix([(extra.T * t * extra)[0].expand()]).jacobian(x)
    extended = a.col_join(row)
    witness = extended.T.nullspace()[0]
    assert (witness.T * extended).applyfunc(sp.simplify) == sp.zeros(1, 9)
    violation = sp.simplify((witness.T * cubic)[0])
    assert violation != 0
    distance_squared = sp.simplify(violation**2 / (witness.T * witness)[0])
    rays = np.array([list(k) for k in directions], dtype=float)
    charges = np.array(cubic, dtype=float).ravel()
    short = fit_null_charges(charges[:9], rays[:9])
    np.testing.assert_array_equal(short.require_consistent(error_budget=0), np.zeros((4, 4)))
    fit = fit_null_charges(charges, rays)
    assert fit.residual_norm == pytest.approx(float(sp.sqrt(distance_squared)), rel=2e-14)
    expected_witness = np.array(witness, dtype=float).ravel()
    expected_witness /= np.linalg.norm(expected_witness)
    assert abs(fit.witness @ expected_witness) == pytest.approx(1., abs=2e-14)
    assert fit.witness_design_defect < 2e-14
    assert fit.witness_charge == pytest.approx(fit.residual_norm, rel=2e-14)
    with pytest.raises(ValueError, match="dependent-family"):
        fit.require_consistent(error_budget=fit.residual_norm / 2)
    fit.require_consistent(error_budget=fit.residual_norm * 1.01)


def test_frobenius_basis_and_sharp_noise_bound():
    basis = tracefree_basis()
    np.testing.assert_allclose(np.einsum("aij,bij->ab", basis, basis), np.eye(9), atol=4e-16)
    np.testing.assert_allclose(np.einsum("aij,ij->a", basis, ETA), 0., atol=3e-16)
    rays = np.array(generic_null_directions(17, seed=17))
    a = np.array([[k @ b @ k for b in basis] for k in rays])
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    t = np.array([[2., 1., 3., -2.], [1., 0., 2., 1.],
                  [3., 2., -1., 4.], [-2., 1., 4., 3.]])
    delta = 1e-3
    # Worst-case noise realizes equality, not just an arbitrary easy sample.
    q = charges_of(t, rays) + delta * u[:, -1]
    fit = fit_null_charges(q, rays)
    error = np.linalg.norm(fit.tensor - eta_project_out(t))
    assert fit.noise_amplification == pytest.approx(1 / s[-1], rel=2e-14)
    assert error == pytest.approx(fit.tensor_error_bound(delta), rel=1e-10)


@pytest.mark.parametrize("scale", [1e-200, 1., 1e200])
def test_tiny_and_large_inconsistent_charges_cannot_pass_by_norm_underflow(scale):
    rays = np.vstack([tomography_directions(), [1., 1., 0., 0.]])
    q = np.zeros(10)
    q[-1] = scale
    fit = fit_null_charges(q, rays)
    # Duplicated first ray: exact distance to image is |q10-q1|/sqrt(2).
    assert fit.residual_norm / scale == pytest.approx(1 / np.sqrt(2), rel=2e-14)
    assert fit.witness_charge / scale == pytest.approx(1 / np.sqrt(2), rel=2e-14)
    assert fit.witness_design_defect < 2e-14
    with pytest.raises(ValueError, match="dependent-family"):
        fit.require_consistent(error_budget=scale / 4)


def test_ray_representatives_and_permutation_do_not_change_fit_or_consistency():
    rays = np.array(generic_null_directions(15, seed=42))
    t = np.diag([2., 3., -1., 4.])
    q = charges_of(t, rays)
    q[0] += .25  # Also check invariance of an inconsistent family.
    reference = fit_null_charges(q, rays)
    weights = np.logspace(-140, 140, len(rays)) * (-1.)**np.arange(len(rays))
    permutation = np.random.default_rng(4).permutation(len(rays))
    fit = fit_null_charges((q * weights * weights)[permutation],
                           (rays * weights[:, None])[permutation])
    np.testing.assert_allclose(fit.tensor, reference.tensor, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(fit.singular_values, reference.singular_values, rtol=2e-14)
    assert fit.residual_norm == pytest.approx(reference.residual_norm, rel=1e-12)


@pytest.mark.parametrize("ray_scale,tensor_scale", [(1e200, 1e-200), (1e-200, 1e200)])
def test_compensating_tensor_and_ray_units_need_no_unrepresentable_square(ray_scale, tensor_scale):
    rays = tomography_directions()
    t = np.ones((4, 4)) + np.eye(4)
    raw = charges_of(t * tensor_scale, rays * ray_scale)
    expected_scale = ray_scale  # tensor_scale * ray_scale**2, without overflow
    np.testing.assert_allclose(raw / expected_scale, charges_of(t, rays), rtol=1e-14)
    fit = fit_null_charges(raw, rays * ray_scale)
    np.testing.assert_allclose(fit.tensor / tensor_scale, eta_project_out(t), atol=2e-14)
    with pytest.raises(ValueError, match="range"):
        design_matrix(rays * ray_scale)


def test_lorentz_covariance_of_consistent_reconstruction_modulo_metric():
    rapidity = .8
    c, s = np.cosh(rapidity), np.sinh(rapidity)
    boost = np.array([[c, s, 0., 0.], [s, c, 0., 0.],
                      [0., 0., 1., 0.], [0., 0., 0., 1.]])
    inverse = np.linalg.inv(boost)
    rays = np.array(generic_null_directions(16, seed=8))
    rng = np.random.default_rng(9)
    m = rng.normal(size=(4, 4))
    t = m + m.T
    moved_t = inverse.T @ t @ inverse
    moved_t = (moved_t + moved_t.T) / 2
    moved_rays = rays @ boost.T
    q = charges_of(t, rays)
    np.testing.assert_allclose(charges_of(moved_t, moved_rays), q, atol=1e-14)
    fit = fit_null_charges(q, moved_rays)
    np.testing.assert_allclose(fit.tensor, eta_project_out(moved_t), atol=2e-13)


def test_distinct_planar_and_numerically_collapsed_directions_do_not_supply_rank_nine():
    theta = np.linspace(0, 2*np.pi, 12, endpoint=False)
    planar = np.array([null_vector([np.cos(t), np.sin(t), 0.]) for t in theta])
    collapsed = np.array([null_vector([1., 1e-9*np.cos(t), 1e-9*np.sin(t)]) for t in theta])
    for rays in (planar, collapsed):
        with pytest.raises(ValueError, match="rank"):
            fit_null_charges(np.zeros(len(rays)), rays)


@pytest.mark.parametrize("data", [np.full(9, np.nan), np.full(9, np.inf), np.ones((9, 1)),
                                   np.ones(9, dtype=complex), np.ones(9, dtype=bool),
                                   ["1"]*9, [0.]*8, []])
def test_invalid_charges_are_not_coerced_into_valid_receipts(data):
    with pytest.raises(ValueError):
        fit_null_charges(data, tomography_directions())


@pytest.mark.parametrize("data", [[], np.zeros((9, 4)), np.ones((9, 3)), np.ones((9, 5)),
                                   np.full((9, 4), np.inf), np.eye(4),
                                   np.ones((9, 4), dtype=complex)])
def test_invalid_rays_fail_before_reconstruction(data):
    with pytest.raises(ValueError):
        fit_null_charges(np.zeros(len(data)), data)


@pytest.mark.parametrize("budget", [-1., np.nan, np.inf, True, 1j, [1.]])
def test_consistency_and_noise_budgets_cannot_be_nan_negative_or_coerced(budget):
    fit = fit_null_charges(np.zeros(9), tomography_directions())
    with pytest.raises(ValueError):
        fit.require_consistent(error_budget=budget)
    with pytest.raises(ValueError):
        fit.tensor_error_bound(budget)


def test_cancellation_does_not_erase_a_representable_quadratic_charge():
    t = np.diag([1e200, 0., -1e200, 0.])
    t[0, 2] = t[2, 0] = 1e-200
    # Independent analytic cancellation: T00 + 2*T02 + T22 = 2e-200.
    assert charges_of(t, [[1., 0., 1., 0.]])[0] == 2e-200


def test_unresolvable_charge_dynamic_range_and_output_units_fail_closed():
    rays = tomography_directions()
    q = np.ones(9)
    q[0], q[1] = 1e308, 1e-308
    with pytest.raises(ValueError, match="range"):
        fit_null_charges(q, rays)
    for scale in (1e200, 1e-200):
        with pytest.raises(ValueError, match="range"):
            charges_of(np.eye(4), rays * scale)
        with pytest.raises(ValueError, match="range"):
            fit_null_charges(np.ones(9), rays * scale)


def test_extended_exponent_inputs_cannot_be_silently_narrowed():
    if np.finfo(np.longdouble).maxexp == np.finfo(float).maxexp:
        pytest.skip("longdouble has no wider exponent range on this platform")
    for value in ("1e400", "1e-400"):
        with pytest.raises(ValueError, match="binary64"):
            fit_null_charges(np.full(9, np.longdouble(value)), tomography_directions())


def test_integer_conversion_cannot_erase_a_dependent_family_violation():
    rays = np.vstack([tomography_directions(), [1., 1., 0., 0.]])
    q = np.full(10, 2**60, dtype=np.int64)
    q[-1] += 1  # These different readings of the same ray round to one float.
    with pytest.raises(ValueError, match="representable"):
        fit_null_charges(q, rays)


def test_direction_conversion_cannot_erase_a_nonzero_component():
    with pytest.raises(ValueError, match="range"):
        null_vector([1e308, 1e-308, 0.])
    with pytest.raises(ValueError, match="range"):
        charges_of(np.eye(4), [[1e308, 1e308, 1e-308, 0.]])


def test_nonzero_noise_bound_cannot_underflow_to_zero():
    rays = generic_null_directions(200, seed=20)
    fit = fit_null_charges(np.zeros(200), rays)
    assert fit.noise_amplification < .5
    with pytest.raises(ValueError, match="range"):
        fit.tensor_error_bound(np.nextafter(0., 1.))
