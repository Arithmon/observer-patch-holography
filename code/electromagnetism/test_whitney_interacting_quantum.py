"""Independent finite checks for the interacting metric and Gauss reduction.

Exact simplex-moment formulae check the full kinetic metric at zero gauge
potential with nonuniform complex matter. Generic nonzero-potential tests
exercise covariance, the actual Schur complement, scalar-potential signs,
residual charge, and the zero-matter stabilizer. These numerical regressions
do not substitute for the paper's closed-form and operator-domain proofs.
"""
from functools import lru_cache
from dataclasses import replace
from math import factorial
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_interacting_quantum as quantum


@lru_cache(maxsize=1)
def fixture_data():
    mesh = quantum.geometry(4)
    rng = np.random.default_rng(48321)
    a = mesh.slice[:42, :30]@rng.normal(size=30)/5
    psi = rng.normal(size=13)+1j*rng.normal(size=13)
    return mesh, a, psi


def direct_zero_potential_metric(mesh, psi, charge, scalar_factor=2):
    """Integrate monomials explicitly; no quadrature or production fields."""
    expected = np.zeros((68, 68))
    for tet in mesh.tetrahedra:
        xyz = mesh.vertices[list(tet)]
        volume = abs(np.linalg.det(xyz[1:]-xyz[0]))/6
        gradients = np.linalg.inv(np.column_stack((np.ones(4), xyz)))[1:].T
        edge_columns, scalar_columns = [], []
        for edge, (u, v) in enumerate(mesh.edges):
            if u not in tet or v not in tet:
                continue
            i, j = tet.index(u), tet.index(v)
            edge_columns.append((edge, i, j))
            # d_a Psi = ie lambda_i lambda_j (psi_i - psi_j) at a=0.
            scalar_columns.append((edge, 1j*charge*(psi[u]-psi[v]), (i, j)))
        for i, vertex in enumerate(tet):
            scalar_columns.extend([(42+vertex, 1, (i,)), (55+vertex, 1j, (i,))])
        for e, i, j in edge_columns:
            for f, k, l in edge_columns:
                expected[e, f] += volume/20*(
                    (1+int(i == k))*(gradients[j]@gradients[l])
                    -(1+int(i == l))*(gradients[j]@gradients[k])
                    -(1+int(j == k))*(gradients[i]@gradients[l])
                    +(1+int(j == l))*(gradients[i]@gradients[k]))
        for col, value, powers in scalar_columns:
            for other, other_value, other_powers in scalar_columns:
                degree = powers+other_powers
                integral = 6*volume/factorial(len(degree)+3)
                for i in range(4):
                    integral *= factorial(degree.count(i))
                expected[col, other] += scalar_factor*np.real(np.conj(value)*other_value)*integral
    return expected


def test_actual_cone_full_metric_against_independent_simplex_moments():
    mesh, _, psi = fixture_data()
    assert (len(mesh.vertices), len(mesh.edges), len(mesh.tetrahedra)) == (13, 42, 20)
    charge = 0.7
    metric, _ = quantum.coefficients(np.zeros(42), psi, charge, mesh=mesh)
    independent = direct_zero_potential_metric(mesh, psi, charge)
    np.testing.assert_allclose(metric, independent, atol=2e-13, rtol=2e-13)
    # Dropping the complex scalar normalization is detected on all blocks.
    wrong = direct_zero_potential_metric(mesh, psi, charge, scalar_factor=1)
    assert np.linalg.norm(metric-wrong) > 1


def test_configuration_jacobian_and_vertical_identity_with_chain_rule_mutation():
    mesh, a, psi = fixture_data()
    charge = 0.7
    value, jacobian, _ = quantum.scalar_fields(a, psi, charge, mesh)
    rng = np.random.default_rng(5147)
    direction = rng.normal(size=68)
    eps = 1e-5
    def displaced(h):
        return quantum.scalar_fields(a+h*direction[:42],
            psi+h*(direction[42:55]+1j*direction[55:]), charge, mesh)[0]
    # Independent five-point derivative of the actual interpolated values.
    derivative = (displaced(-2*eps)-8*displaced(-eps)+8*displaced(eps)-displaced(2*eps))/(12*eps)
    np.testing.assert_allclose(jacobian@direction, derivative, atol=2e-10, rtol=2e-10)
    phi = rng.normal(size=13)
    expected = 1j*charge*(mesh.nodal@phi)*value
    actual = jacobian@quantum.vertical(psi, charge, mesh)@phi
    np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=2e-14)
    undressed = jacobian.copy()
    undressed[:, :42] = 0
    assert np.linalg.norm(undressed@quantum.vertical(psi, charge, mesh)@phi-expected) > 1
    wrong_sign = quantum.vertical(psi, charge, mesh).copy()
    wrong_sign[:42] *= -1
    assert np.linalg.norm(jacobian@wrong_sign@phi-expected) > 1


def test_generic_local_gauge_covariance_of_metric_and_potential():
    mesh, a, psi = fixture_data()
    charge = 0.7
    xi = mesh.mean_zero@np.linspace(-0.6, 0.9, 12)
    phase = np.exp(1j*charge*xi)
    metric, potential = quantum.coefficients(a, psi, charge, mesh=mesh)
    moved, moved_potential = quantum.coefficients(a+mesh.d@xi, phase*psi, charge, mesh=mesh)
    tangent = np.eye(68)
    tangent[42:55, 42:55] = np.diag(phase.real)
    tangent[42:55, 55:] = -np.diag(phase.imag)
    tangent[55:, 42:55] = np.diag(phase.imag)
    tangent[55:, 55:] = np.diag(phase.real)
    np.testing.assert_allclose(tangent.T@moved@tangent, metric, atol=2e-13, rtol=2e-13)
    np.testing.assert_allclose(moved_potential, potential, atol=3e-13, rtol=3e-13)


def test_coulomb_representative_is_unique_along_mean_zero_gauge_orbit():
    mesh, a, psi = fixture_data()
    xi = mesh.mean_zero@np.linspace(-1.0, 1.0, 12)
    recovered_a, recovered_psi, correction = quantum.coulomb_representative(
        a+mesh.d@xi, np.exp(0.7j*xi)*psi, 0.7, mesh)
    np.testing.assert_allclose(recovered_a, a, atol=2e-14, rtol=2e-14)
    np.testing.assert_allclose(recovered_psi, psi, atol=2e-14, rtol=2e-14)
    np.testing.assert_allclose(correction, -xi, atol=2e-14, rtol=2e-14)
    with pytest.raises(ValueError, match="Coulomb"):
        quantum.reduced_coefficients(a+mesh.d@xi, psi, mesh=mesh)


def test_schur_elimination_minimizes_actual_metric_and_rejects_principal_block():
    mesh, a, psi = fixture_data()
    metric, _ = quantum.coefficients(a, psi, 0.7, mesh=mesh)
    gamma, _, eta_map, inertia = quantum.reduced_coefficients(a, psi, 0.7, mesh=mesh)
    assert gamma.shape == (56, 56) and inertia.shape == (12, 12)
    assert np.linalg.eigvalsh(gamma).min() > 0.01
    r = quantum.vertical(psi, 0.7, mesh)@mesh.mean_zero
    rng = np.random.default_rng(9511)
    v, displacement = rng.normal(size=56), rng.normal(size=12)
    eta = eta_map@v
    horizontal = mesh.slice@v+r@eta
    np.testing.assert_allclose(r.T@metric@horizontal, 0, atol=3e-13)
    np.testing.assert_allclose(horizontal@metric@horizontal, v@gamma@v, atol=3e-12)
    trial = mesh.slice@v+r@(eta+displacement)
    np.testing.assert_allclose(trial@metric@trial-v@gamma@v,
                               displacement@inertia@displacement, atol=3e-12)
    # A Coulomb tangent is not horizontal for the scalar-coupled metric.
    assert np.linalg.norm(r.T@metric@mesh.slice@v) > 0.1
    assert np.linalg.norm(mesh.slice.T@metric@mesh.slice-gamma) > 0.1
    assert np.linalg.norm(r.T@metric@(mesh.slice@v-r@eta)) > 0.1


def test_residual_constant_gauss_equation_and_neutral_velocity():
    mesh, a, psi = fixture_data()
    charge, constant_phi = 0.7, 0.4
    metric, _ = quantum.coefficients(a, psi, charge, mesh=mesh)
    gamma, _, eta_map, _ = quantum.reduced_coefficients(a, psi, charge, mesh=mesh)
    j = np.concatenate((np.zeros(30), -psi.imag, psi.real))
    w = np.linspace(-0.4, 0.8, 56)+charge*constant_phi*j
    r = quantum.vertical(psi, charge, mesh)
    horizontal = mesh.slice@w+r@mesh.mean_zero@(eta_map@w)
    gauss = r.T@metric@horizontal
    np.testing.assert_allclose(mesh.mean_zero.T@gauss, 0, atol=1e-13)
    np.testing.assert_allclose(gauss.sum(), charge*(j@gamma@w), atol=2e-13)
    assert abs(gauss.sum()) > 0.1  # Twelve equations do not imply the thirteenth.
    neutral = w-j*(j@gamma@w)/(j@gamma@j)
    horizontal_neutral = mesh.slice@neutral+r@mesh.mean_zero@(eta_map@neutral)
    np.testing.assert_allclose(r.T@metric@horizontal_neutral, 0, atol=2e-13)


def test_residual_circle_preserves_reduced_metric_and_density():
    mesh, a, psi = fixture_data()
    alpha = 0.61
    gamma, potential, _, _ = quantum.reduced_coefficients(a, psi, 0.7, mesh=mesh)
    moved, moved_potential, _, _ = quantum.reduced_coefficients(a, np.exp(1j*alpha)*psi, 0.7, mesh=mesh)
    rotation = np.eye(56)
    rotation[30:43, 30:43] = np.cos(alpha)*np.eye(13)
    rotation[30:43, 43:] = -np.sin(alpha)*np.eye(13)
    rotation[43:, 30:43] = np.sin(alpha)*np.eye(13)
    rotation[43:, 43:] = np.cos(alpha)*np.eye(13)
    np.testing.assert_allclose(rotation.T@moved@rotation, gamma, atol=3e-13, rtol=3e-13)
    np.testing.assert_allclose(moved_potential, potential, atol=3e-13, rtol=3e-13)
    np.testing.assert_allclose(np.linalg.slogdet(moved)[1], np.linalg.slogdet(gamma)[1], atol=3e-13)


def test_zero_matter_metric_is_regular_while_full_gauge_orbit_has_stabilizer():
    mesh, a, _ = fixture_data()
    psi = np.zeros(13, dtype=complex)
    metric, _ = quantum.coefficients(a, psi, 0.7, mesh=mesh)
    gamma, _, eta_map, inertia = quantum.reduced_coefficients(a, psi, 0.7, mesh=mesh)
    assert np.linalg.eigvalsh(metric).min() > 0.01
    assert np.linalg.eigvalsh(gamma).min() > 0.01
    assert np.linalg.eigvalsh(inertia).min() > 0.01
    np.testing.assert_allclose(eta_map, 0, atol=1e-14)
    np.testing.assert_allclose(gamma, mesh.slice.T@metric@mesh.slice, atol=1e-13)
    r = quantum.vertical(psi, 0.7, mesh)
    assert np.linalg.matrix_rank(r.T@metric@r, tol=1e-10) == 12
    np.testing.assert_allclose(r@np.ones(13), 0, atol=1e-14)


def test_generic_coefficients_have_controlled_quadrature_order_sensitivity():
    mesh, a, psi = fixture_data()
    results = []
    for order in (4, 5, 6):
        # Hold the physical coordinate frame fixed; separately computed SVD
        # null-space frames need not use the same transverse coordinates.
        refined = replace(quantum.geometry(order), slice=mesh.slice, mean_zero=mesh.mean_zero)
        metric, potential = quantum.coefficients(a, psi, 0.7, mesh=refined)
        gamma, _, _, _ = quantum.reduced_coefficients(a, psi, 0.7, mesh=refined)
        results.append((metric, potential, gamma))
    coarse = np.array([np.max(abs(a-b)) for a, b in zip(results[0], results[1], strict=True)])
    fine = np.array([np.max(abs(a-b)) for a, b in zip(results[1], results[2], strict=True)])
    assert np.all(coarse > 1e-8)  # Gauge invariance alone hides these errors.
    assert np.all(coarse < 2e-6)
    assert np.all(fine < 1e-8) and np.all(fine < coarse/100)
    # These are finite order-comparison diagnostics, not interval error bounds.


@pytest.mark.parametrize("parameter,value", [("mass_squared", -1), ("quartic", -1),
                                            ("mass_squared", np.inf), ("quartic", np.nan)])
def test_reject_negative_or_nonfinite_potential_inputs(parameter, value):
    mesh, a, psi = fixture_data()
    with pytest.raises(ValueError, match="nonnegative finite"):
        quantum.coefficients(a, psi, mesh=mesh, **{parameter: value})


@pytest.mark.parametrize("charge", [np.nan, np.inf, -np.inf, 1j, [0.7]])
def test_reject_nonfinite_nonreal_or_nonscalar_charge(charge):
    mesh, a, psi = fixture_data()
    with pytest.raises(ValueError, match="finite real charge"):
        quantum.coefficients(a, psi, charge, mesh=mesh)


@pytest.mark.parametrize("part", ["edge_shape", "node_shape", "edge_nan", "node_inf", "complex_edge"])
def test_reject_malformed_configuration(part):
    mesh, a, psi = fixture_data()
    a, psi = a.copy(), psi.copy()
    if part == "edge_shape":
        a = a[:-1]
    elif part == "node_shape":
        psi = psi[:-1]
    elif part == "edge_nan":
        a[0] = np.nan
    elif part == "node_inf":
        psi[0] = np.inf
    else:
        a = a.astype(complex)+1j
    with pytest.raises(ValueError, match="coefficients required"):
        quantum.coefficients(a, psi, mesh=mesh)


def test_large_unwrapped_phases_obey_pointwise_dominant_node_bound():
    mesh, a, psi = fixture_data()
    for scale in (0, 1e3, 1e9):
        _, jacobian, _ = quantum.scalar_fields(scale*a, psi, 0.7, mesh)
        dressed = jacobian[:, 42:55]
        for i in range(13):
            v = 0.7*np.exp(1j*np.arange(13))
            v[i] = 1
            corner = mesh.nodal[:, i] >= 0.75
            assert corner.any()
            actual = abs(dressed[corner]@v)
            bound = 2*mesh.nodal[corner, i]-1
            assert np.all(actual >= bound-2e-14)
    # These are pointwise checks; fixed quadrature does not resolve the
    # oscillatory integrals at these large unwrapped configurations.
    lam = np.array([3/4, 1/12, 1/12, 1/12])
    phases = np.array([1, -1, -1, -1])
    assert lam@phases == pytest.approx(1/2)
    assert abs(lam[1:]@phases[1:]) < 1/2  # deleting the dominant term fails
    assert Fraction(1, 4)**3*Fraction(1, 4)*Fraction(1, 4) == Fraction(1, 1024)


def test_gaussian_normalization_by_independent_56_dimensional_radial_integral():
    from scipy.integrate import quad
    from scipy.special import gammaln
    sigma = 0.5
    log_area = np.log(2)+28*np.log(np.pi)-gammaln(28)
    def radial(r):
        if r == 0:
            return 0
        q = np.zeros(56)
        q[0] = r
        return np.exp(log_area+55*np.log(r)+2*quantum.gaussian_half_density_log(q, sigma))
    integral, error = quad(radial, 0, np.inf, epsabs=2e-12, epsrel=2e-12)
    assert integral == pytest.approx(1, abs=3e-12)
    assert error < 3e-12
    # An amplitude normalized as if dimension were 52 fails the radial check.
    assert abs(integral*(2*np.pi*sigma**2)**2-1) > 1


def test_gaussian_log_amplitude_density_factor_and_neutrality():
    mesh, a, psi = fixture_data()
    q = np.concatenate((mesh.slice[:42, :30].T@a, psi.real, psi.imag))
    logg = quantum.gaussian_half_density_log(q, Fraction(1, 2))
    logf = quantum.gaussian_state_log_amplitude(q, 0.5, 0.7, mesh)
    gamma, _, _, _ = quantum.reduced_coefficients(a, psi, 0.7, mesh=mesh)
    logdet = np.linalg.slogdet(gamma)[1]
    assert logf+logdet/4 == pytest.approx(logg, abs=3e-13)
    assert abs(logf+logdet/2-logg) > 1  # density versus metric determinant
    moved = np.exp(0.63j)*psi
    rotated = np.concatenate((q[:30], moved.real, moved.imag))
    assert quantum.gaussian_half_density_log(rotated, 0.5) == pytest.approx(logg, abs=3e-13)
    assert quantum.gaussian_state_log_amplitude(rotated, 0.5, 0.7, mesh) == pytest.approx(logf, abs=3e-13)


def test_gaussian_matter_moments_against_exact_wick_and_simplex_integration():
    sigma, volume, trace = Fraction(2, 3), Fraction(7, 5), Fraction(11, 7)
    def simplex_mean(powers):
        result = Fraction(factorial(3), factorial(3+sum(powers)))
        for power in powers:
            result *= factorial(power)
        return result
    quadratic = sum(simplex_mean(tuple(2*int(k == i) for k in range(4))) for i in range(4))
    quartic = sum(simplex_mean(tuple(2*int(k == i)+2*int(k == j) for k in range(4)))
                  for i in range(4) for j in range(4))
    moments = quantum.gaussian_initial_moments(sigma, volume, trace)
    # Independent real/imaginary variance and circular-Gaussian Wick pairing.
    assert moments['matter_l2'] == 2*sigma**2*quadratic*volume
    assert moments['matter_l4'] == 2*(2*sigma**2)**2*quartic*volume
    assert moments['magnetic_energy'] == sigma**2*trace/2
    assert all(isinstance(value, Fraction) for value in moments.values())
    assert moments['matter_l2'] != sigma**2*quadratic*volume
    assert moments['matter_l4'] != (2*sigma**2)**2*quartic*volume


def test_gaussian_state_rejects_wrong_coordinate_frame():
    mesh, _, _ = fixture_data()
    section = mesh.slice.copy()
    section[:, 0] *= 2
    with pytest.raises(ValueError, match='orthonormal Coulomb frame'):
        quantum.gaussian_state_log_amplitude(np.zeros(56), mesh=replace(mesh, slice=section))
    # An orthogonal exchange of field and scalar columns is also invalid.
    section = mesh.slice.copy()
    section[:, [0, 30]] = section[:, [30, 0]]
    with pytest.raises(ValueError, match='orthonormal Coulomb frame'):
        quantum.gaussian_state_log_amplitude(np.zeros(56), mesh=replace(mesh, slice=section))


@pytest.mark.parametrize('sigma', [0, -1, np.nan, np.inf, 1j, [1], True])
def test_gaussian_rejects_invalid_width(sigma):
    with pytest.raises(ValueError, match='sigma'):
        quantum.gaussian_half_density_log(np.zeros(56), sigma)


@pytest.mark.parametrize('q', [np.zeros(55), np.zeros((56, 1)), np.full(56, np.nan),
                             np.full(56, np.inf), np.ones(56, dtype=complex)*1j])
def test_gaussian_rejects_invalid_coordinates(q):
    with pytest.raises(ValueError, match='56-dimensional'):
        quantum.gaussian_half_density_log(q)


@pytest.mark.parametrize('arguments', [(0, 1, 1), (1, 0, 1), (1, 1, -1),
                                     (1, np.inf, 1), (1, 1, 1j), (True, 1, 1)])
def test_gaussian_moments_reject_invalid_parameters(arguments):
    with pytest.raises(ValueError, match='required sign'):
        quantum.gaussian_initial_moments(*arguments)
