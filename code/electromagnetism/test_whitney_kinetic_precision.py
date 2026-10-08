"""Original-input controls for positive Whitney kinetic reduction.

Reviewed baseline: 7521d5c3. At a=0 and uniform real matter T, the full
kinetic form is diag(M, 2N, 2N). Eliminating the mean-zero gauge velocity
must leave a positive form. These controls independently assemble the exact
degree-two Whitney mass, use simplex/Beta moments for each action component,
and diagonalize the scalar sector analytically. No producer quadrature,
thermal state, rounded Schur matrix or producer support decision is an oracle.

The API may refuse unresolved dense binary64 geometry with a 'kinetic
precision' ValueError. Ordinary resolved inputs must remain evaluable; the
factor-based Gaussian observable is required at T=1e6 even when materializing
a sufficiently accurate dense reduced metric is inappropriate. These finite
precision controls are not interval certificates or new physical claims.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_interacting_quantum as quantum
import whitney_quantum_packet as packet


@lru_cache(maxsize=1)
def original_geometry():
    """Exact-degree original cone forms, evaluated in a private 80-digit context.

    The twelve supplied golden vertices fix 30 edges of squared length four
    and 20 triangular faces. The barycentric gradient Gram matrix and simplex
    moment integral(lambda_i lambda_j)=vol*(1+delta_ij)/20 determine M.
    Edge orientation/order is immaterial to the determinant quotient below.
    """
    mp = mpmath.mp.clone()
    mp.dps = 80
    root = mp.sqrt(5)
    phi = (1 + root) / 2
    vertices = [(0, -1, -phi), (-1, -phi, 0), (-phi, 0, -1),
                (1, -phi, 0), (0, 1, -phi), (-phi, 0, 1),
                (phi, 0, -1), (0, -1, phi), (-1, phi, 0),
                (phi, 0, 1), (1, phi, 0), (0, 1, phi)]
    boundary = [(i, j) for i, j in combinations(range(12), 2)
                if abs(sum((vertices[i][k] - vertices[j][k]) ** 2 for k in range(3)) - 4)
                < mp.mpf("1e-70")]
    faces = [face for face in combinations(range(12), 3)
             if all(edge in boundary for edge in combinations(face, 2))]
    assert len(boundary) == 30 and len(faces) == 20
    edges = [(0, j) for j in range(1, 13)] + [(i + 1, j + 1) for i, j in boundary]
    volume, kappa = 10 + 10 * root / 3, 30 - 10 * root
    element_volume = (3 + root) / 6
    gradients = mp.matrix(4)
    for i in range(4):
        for j in range(4):
            gradients[i, j] = (kappa / volume if i == j == 0 else
                               (-7 + 3 * root) / 2 if i == 0 or j == 0 else
                               (3 - root) / 2 if i == j else 1 - root / 2)
    mass, incidence = mp.matrix(42), mp.matrix(42, 13)
    for e, (left, right) in enumerate(edges):
        incidence[e, left], incidence[e, right] = -1, 1
    for face in faces:
        tet = (0, *(j + 1 for j in face))
        local = [(e, tet.index(i), tet.index(j)) for e, (i, j) in enumerate(edges)
                 if i in tet and j in tet]
        for e, i, j in local:
            for f, k, l in local:
                mass[e, f] += element_volume / 20 * (
                    (1 + (i == k)) * gradients[j, l] - (1 + (i == l)) * gradients[j, k]
                    - (1 + (j == k)) * gradients[i, l] + (1 + (j == l)) * gradients[i, k])
    nodal = mp.matrix(13)
    nodal[0, 0] = volume / 10
    for j in range(1, 13):
        nodal[0, j] = nodal[j, 0] = volume / 80
        nodal[j, j] = volume / 40
    for i, j in boundary:
        nodal[i + 1, j + 1] = nodal[j + 1, i + 1] = volume / 200
    # If U spans ker(D.T M) orthonormally and C spans im(D), orthogonal
    # determinant complementation gives det(U.T M U)=det(M)det(C.T M C)
    # /det((M C).T(M C)). Pinning one vertex supplies a convenient basis C.
    gauge = incidence[:, 1:]
    mg = mass * gauge
    logdet_zero = (mp.log(mp.det(mass)) + mp.log(mp.det(gauge.T * mg))
                   - mp.log(mp.det(mg.T * mg)) + 2 * mp.log(mp.det(2 * nodal)))
    return mp, volume, kappa, nodal, logdet_zero


def uniform_oracle(amplitude, charge):
    """Absolute determinant, all 13 scalar eigenvalues and radial components.

    For x=(12,-1,...,-1), invariance and strict convexity reduce the full
    minimization to A*z^2+B*(1+e*T*z)^2, A=169*kappa, B=114*V/5.
    The other 11 scalar modes belong to the icosahedron adjacency eigenvalues
    sqrt(5), -sqrt(5), -1, with multiplicities 3,3,5. The remaining radial
    2x2 block uses its positive trace/determinant, rationalizing its small root.
    """
    mp, volume, kappa, nodal, logdet = original_geometry()
    amplitude, charge = mp.mpf(amplitude), mp.mpf(charge)
    coupling = amplitude * charge
    squared = coupling * coupling
    a, b = 169 * kappa, 114 * volume / 5
    remainder = a / (a + b * squared)
    eta = -b * coupling / (a + b * squared)
    momentum = [21 * volume * remainder / 10] + [volume * remainder / 5] * 12
    radial_b, radial_l = 19 * volume / 130, 13 * kappa / 12
    determinant = volume ** 2 / 80 * radial_l / (radial_l + squared * radial_b)
    trace = ((3 * volume / 10) * radial_l + squared * volume ** 2 / 80) / (
        radial_l + squared * radial_b)
    large = (trace + mp.sqrt(trace * trace - 4 * determinant)) / 2
    scalar_eigenvalues = [determinant / large, large]
    inertia_eigenvalues = [radial_l + squared * radial_b]
    logdet -= mp.log1p(squared * radial_b / radial_l)
    for adjacency, count in ((mp.sqrt(5), 3), (-mp.sqrt(5), 3), (-1, 5)):
        mass = volume * (5 + adjacency) / 100
        stiffness = mp.mpf(5) / 3 + (1 - mp.sqrt(5)) * adjacency / 6
        scalar_eigenvalues += [mass * stiffness / (stiffness + squared * mass)] * count
        inertia_eigenvalues += [stiffness + squared * mass] * count
        logdet -= count * mp.log1p(squared * mass / stiffness)
    return {
        "energy_twice": float(b * remainder), "eta_scale": float(eta),
        "momentum": np.array(list(map(float, momentum))),
        "moment_map": float(amplitude * sum(momentum)),
        "scalar_eigenvalues": np.array(sorted(map(float, scalar_eigenvalues))),
        "inertia_eigenvalues": np.array(sorted(map(float, inertia_eigenvalues))),
        "logdet": float(logdet), "scalar_mass": np.array((2 * nodal).tolist(), dtype=float),
        "potential": float(volume * (amplitude ** 2 / 2 + amplitude ** 4 / 8)),
    }


def state(amplitude, phase=1):
    psi = np.full(13, amplitude * phase, dtype=complex)
    x = np.r_[12., -np.ones(12)]
    velocity = np.r_[np.zeros(30), (1j * phase * x).real, (1j * phase * x).imag]
    return psi, x, velocity


def check_dense_components(result, mesh, amplitude, charge, phase=1, rtol=1e-8):
    gamma, potential, eta_map, inertia = result
    expected = uniform_oracle(amplitude, charge)
    _, x, velocity = state(amplitude, phase)
    # The unrotated, uncoupled scalar block must survive; checking a ratio or
    # determinant alone would miss correlated errors in the two sectors.
    uncoupled, reduced = ((slice(30, 43), slice(43, 56)) if phase == 1
                          else (slice(43, 56), slice(30, 43)))
    np.testing.assert_allclose(gamma[uncoupled, uncoupled], expected["scalar_mass"], rtol=2e-12, atol=2e-13)
    np.testing.assert_allclose(np.linalg.eigvalsh(gamma[reduced, reduced]),
                               expected["scalar_eigenvalues"], rtol=rtol, atol=0)
    np.testing.assert_allclose(np.linalg.eigvalsh(inertia), expected["inertia_eigenvalues"], rtol=2e-12, atol=0)
    np.testing.assert_allclose(mesh.mean_zero @ eta_map @ velocity, expected["eta_scale"] * x,
                               rtol=2e-10, atol=1e-13)
    observed = float(velocity @ gamma @ velocity)
    assert observed > 0
    assert observed == pytest.approx(expected["energy_twice"], rel=rtol, abs=0)
    expected_momentum = expected["momentum"] * (1 if phase == 1 else -1)
    np.testing.assert_allclose((gamma @ velocity)[reduced], expected_momentum, rtol=rtol, atol=0)
    assert potential == pytest.approx(expected["potential"], rel=2e-12, abs=0)
    sign, logdet = np.linalg.slogdet(gamma)
    assert sign == 1
    assert logdet == pytest.approx(expected["logdet"], abs=1e-7, rel=0)


@pytest.mark.parametrize("amplitude,charge,phase", [
    (0., .25, 1), (1., .25, 1), (1000., .25, 1), (1000., -.25, 1),
    (1000., .25, 1j), (2.**30, 0., 1), (2.**30, 2.**-32, 1),
])
def test_resolved_original_input_components_remain_evaluable(amplitude, charge, phase):
    mesh = quantum.geometry()
    psi, _, _ = state(amplitude, phase)
    result = quantum.reduced_coefficients(np.zeros(42), psi, charge, .5, .25, mesh)
    check_dense_components(result, mesh, amplitude, charge, phase)


@pytest.mark.parametrize("amplitude", [1e4, 1e5, 1e6])
def test_dense_precision_boundary_is_accurate_or_explicit(amplitude):
    mesh = quantum.geometry()
    psi, _, _ = state(amplitude)
    try:
        result = quantum.reduced_coefficients(np.zeros(42), psi, .25, .5, .25, mesh)
    except ValueError as error:
        assert "kinetic precision" in str(error)
    else:
        check_dense_components(result, mesh, amplitude, .25)


def test_dense_negative_kinetic_mode_is_refused_at_original_reproduction():
    # Baseline returns -1.6542323066914832e-13 for x.T gamma x instead of the
    # original-input positive minimum 2.0656721888405686e-14. A dense matrix
    # with a mode below its binary64 component precision must not be returned.
    mesh = quantum.geometry()
    psi, _, _ = state(1e9)
    assert uniform_oracle(1e9, .25)["energy_twice"] > 0
    with pytest.raises(ValueError, match="kinetic precision"):
        quantum.reduced_coefficients(np.zeros(42), psi, .25, .5, .25, mesh)


@pytest.mark.parametrize("amplitude,phase", [(1., 1), (1000., 1), (1e6, 1), (1e6, 1j)])
def test_gaussian_log_amplitude_uses_resolved_absolute_factor_density(amplitude, phase):
    psi, _, _ = state(amplitude, phase)
    coordinates = np.r_[np.zeros(30), psi.real, psi.imag]
    width = max(amplitude, 1.)
    mp, _, _, _, _ = original_geometry()
    # Select sigma=T at large T so the Gaussian term cannot drown out a wrong
    # density correction. Compute both terms independently from original inputs.
    logg = -14 * (mp.log(2 * mp.pi) + 2 * mp.log(width)) - mp.mpf(13) * amplitude ** 2 / (4 * mp.mpf(width) ** 2)
    expected = float(logg - uniform_oracle(amplitude, .25)["logdet"] / 4)
    actual = quantum.gaussian_state_log_amplitude(coordinates, width, .25)
    assert actual == pytest.approx(expected, abs=1e-7, rel=0)


def test_gaussian_large_density_is_accurate_or_precision_refusal():
    amplitude = 1e9
    coordinates = np.r_[np.zeros(30), np.full(13, amplitude), np.zeros(13)]
    mp, _, _, _, _ = original_geometry()
    logg = -14 * (mp.log(2 * mp.pi) + 2 * mp.log(amplitude)) - mp.mpf(13) / 4
    expected = float(logg - uniform_oracle(amplitude, .25)["logdet"] / 4)
    try:
        actual = quantum.gaussian_state_log_amplitude(coordinates, amplitude, .25)
    except ValueError as error:
        assert "kinetic precision" in str(error)
    else:
        assert actual == pytest.approx(expected, abs=1e-7, rel=0)


def test_phase_space_cotangent_components_are_independently_correct():
    amplitude = 1000.
    q68 = np.r_[np.zeros(42), np.full(13, amplitude), np.zeros(13)]
    velocity = np.r_[np.zeros(55), 12., -np.ones(12)]
    result = packet.phase_space(q68, velocity)
    expected = uniform_oracle(amplitude, .25)
    np.testing.assert_allclose(result["momentum"][43:], expected["momentum"], rtol=1e-8, atol=0)
    assert result["constant_moment_map"] == pytest.approx(expected["moment_map"], rel=1e-8, abs=0)
    assert result["momentum"] @ result["velocity"] == pytest.approx(expected["energy_twice"], rel=1e-8, abs=0)
    assert result["log_rho"] == pytest.approx(expected["logdet"] / 2, abs=1e-7, rel=0)


def test_phase_space_refuses_unresolved_original_cotangent():
    q68 = np.r_[np.zeros(42), np.full(13, 1e9), np.zeros(13)]
    velocity = np.r_[np.zeros(55), 12., -np.ones(12)]
    with pytest.raises(ValueError, match="kinetic precision"):
        packet.phase_space(q68, velocity)
