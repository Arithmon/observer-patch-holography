"""Finite coefficients for the declared interacting Whitney Hamiltonian.

The analytic construction uses exact element integrals. This small numerical
implementation evaluates its full kinetic metric, potential and mean-zero
gauge Schur complement by positive tensor Gauss quadrature. It does not
discretize the 56-dimensional quantum wavefunction or prove self-adjointness.
Geometry is reconstructed from the canonical Lean coordinates/incidences.
The coefficient evaluator also permits zero charge; the paper's neutral
interacting-quantization theorem assumes nonzero charge.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
from itertools import combinations, product
from math import isfinite

import numpy as np
from scipy.linalg import null_space

from verify_cone_whitney_bridge import source_mesh


@dataclass(frozen=True)
class Geometry:
    vertices: np.ndarray
    edges: tuple
    tetrahedra: tuple
    d: np.ndarray
    mass: np.ndarray
    stiffness: np.ndarray
    weights: np.ndarray
    nodal: np.ndarray
    nodal_grad: np.ndarray
    edge_forms: np.ndarray
    paths: np.ndarray
    path_grad: np.ndarray
    mean_zero: np.ndarray
    slice: np.ndarray


@lru_cache(maxsize=4)
def geometry(order=4):
    """The supplied cone with positive Duffy/Gauss element quadrature."""
    if type(order) is not int or not 3 <= order <= 8:
        raise ValueError("quadrature order must be an integer from 3 to 8")
    boundary_vertices, boundary_edges, boundary_faces = source_mesh()
    vertices = np.asarray([(0, 0, 0)] + boundary_vertices, dtype=float)
    edges = tuple([(0, j) for j in range(1, 13)] +
                  [(i+1, j+1) for i, j in boundary_edges])
    tets = tuple((0, *(j+1 for j in face)) for face in boundary_faces)
    d = np.zeros((42, 13))
    for row, (i, j) in enumerate(edges):
        d[row, i], d[row, j] = -1, 1
    points, gauss_weights = np.polynomial.legendre.leggauss(order)
    points, gauss_weights = (points+1)/2, gauss_weights/2
    barycentric, reference_weights = [], []
    for i, j, k in product(range(order), repeat=3):
        r, s, t = points[[i, j, k]]
        barycentric.append([(1-r)*(1-s)*(1-t), r, (1-r)*s, (1-r)*(1-s)*t])
        reference_weights.append(gauss_weights[i]*gauss_weights[j]*gauss_weights[k]
                                 *(1-r)**2*(1-s))
    barycentric, reference_weights = np.array(barycentric), np.array(reference_weights)
    weights, nodal, nodal_grad, edge_forms, paths, path_grad = [], [], [], [], [], []
    stiffness = np.zeros((42, 42))
    for tet in tets:
        xyz = vertices[list(tet)]
        gradients = np.linalg.inv(np.column_stack((np.ones(4), xyz)))[1:].T
        determinant = abs(np.linalg.det(xyz[1:]-xyz[0]))
        local_edges = [(n, tet.index(i), tet.index(j)) for n, (i, j) in enumerate(edges)
                       if i in tet and j in tet]
        curls = np.zeros((42, 3))
        for edge, i, j in local_edges:
            curls[edge] = 2*np.cross(gradients[i], gradients[j])
        stiffness += determinant/6*(curls@curls.T)
        for lam, ref_weight in zip(barycentric, reference_weights, strict=True):
            n, dn = np.zeros(13), np.zeros((13, 3))
            one, path, dp = np.zeros((42, 3)), np.zeros((13, 42)), np.zeros((13, 42, 3))
            n[list(tet)], dn[list(tet)] = lam, gradients
            for edge, i, j in local_edges:
                one[edge] = lam[i]*gradients[j]-lam[j]*gradients[i]
                path[tet[i], edge], path[tet[j], edge] = lam[j], -lam[i]
                dp[tet[i], edge], dp[tet[j], edge] = gradients[j], -gradients[i]
            weights.append(determinant*ref_weight)
            nodal.append(n)
            nodal_grad.append(dn)
            edge_forms.append(one)
            paths.append(path)
            path_grad.append(dp)
    weights, nodal, nodal_grad, edge_forms, paths, path_grad = map(
        np.asarray, (weights, nodal, nodal_grad, edge_forms, paths, path_grad))
    mass = np.einsum("q,qic,qjc->ij", weights, edge_forms, edge_forms)
    mean_zero = null_space(np.ones((1, 13)))
    transverse = null_space(d.T@mass)
    if transverse.shape != (42, 30):
        raise ValueError("unexpected transverse dimension")
    section = np.zeros((68, 56))
    section[:42, :30] = transverse
    section[42:, 30:] = np.eye(26)
    return Geometry(vertices, edges, tets, d, mass, stiffness, weights, nodal,
                    nodal_grad, edge_forms, paths, path_grad, mean_zero, section)


def _configuration(a, psi):
    if not np.isrealobj(a):
        raise ValueError("real edge coefficients required")
    a, psi = np.asarray(a, dtype=float), np.asarray(psi, dtype=complex)
    if a.shape != (42,) or psi.shape != (13,) or not np.isfinite(a).all() or not np.isfinite(psi).all():
        raise ValueError("finite 42-edge and 13-complex-node coefficients required")
    return a, psi


def _real_scalar(value, name):
    if np.ndim(value) != 0 or not np.isrealobj(value) or not np.isfinite(value):
        raise ValueError("finite real "+name+" required")
    return float(value)


def scalar_fields(a, psi, charge, mesh):
    """Value, real configuration Jacobian and covariant spatial gradient."""
    a, psi = _configuration(a, psi)
    charge = _real_scalar(charge, "charge")
    phase = np.exp(1j*charge*np.einsum("qie,e->qi", mesh.paths, a))
    w = mesh.nodal*phase
    value = w@psi
    dressing = 1j*charge*np.einsum("qi,i,qie->qe", w, psi, mesh.paths)
    jacobian = np.column_stack((dressing, w, 1j*w))
    grad_phase = np.einsum("qiec,e->qic", mesh.path_grad, a)
    grad = np.einsum("qi,i,qic->qc", phase, psi, mesh.nodal_grad)
    grad += 1j*charge*np.einsum("qi,i,qic->qc", w, psi, grad_phase)
    potential = np.einsum("qec,e->qc", mesh.edge_forms, a)
    covariant = grad-1j*charge*potential*value[:, None]
    return value, jacobian, covariant


def coefficients(a, psi, charge=1.0, mass_squared=1.0, quartic=1.0, mesh=None):
    """Full 68-real-coordinate kinetic G and the coupled potential V."""
    mesh = geometry() if mesh is None else mesh
    a, psi = _configuration(a, psi)
    if not all(np.ndim(x) == 0 and np.isrealobj(x) and np.isfinite(x) and x >= 0
               for x in (mass_squared, quartic)):
        raise ValueError("nonnegative finite mass-squared and quartic coupling required")
    value, jacobian, covariant = scalar_fields(a, psi, charge, mesh)
    weighted = np.sqrt(mesh.weights)[:, None]*jacobian
    metric = 2*np.real(weighted.conj().T@weighted)
    metric[:42, :42] += mesh.mass
    density = np.sum(abs(covariant)**2, axis=1)+mass_squared*abs(value)**2+quartic*abs(value)**4/2
    potential = a@mesh.stiffness@a/2+mesh.weights@density
    return metric, float(potential)


def vertical(psi, charge, mesh):
    """All thirteen infinitesimal gauge columns in [a, Re psi, Im psi]."""
    psi = np.asarray(psi, dtype=complex)
    charge = _real_scalar(charge, "charge")
    if psi.shape != (13,) or not np.isfinite(psi).all():
        raise ValueError("finite 13-complex-node coefficients required")
    return np.vstack((mesh.d, -charge*np.diag(psi.imag), charge*np.diag(psi.real)))


def reduced_coefficients(a, psi, charge=1.0, mass_squared=1.0, quartic=1.0, mesh=None):
    """Schur metric in the Euclidean-orthonormal slice coordinates.

    Return gamma, V, eta_map, inertia. For a slice velocity v (56 reals),
    the minimizing scalar potential is mean_zero @ eta_map @ v.
    The residual constant potential is retained separately; add e*c*J
    to the slice velocity before applying this map.
    """
    mesh = geometry() if mesh is None else mesh
    a, psi = _configuration(a, psi)
    if np.max(abs(mesh.d.T@mesh.mass@a)) > 1e-9*(1+np.linalg.norm(a)):
        raise ValueError("configuration must lie on the Coulomb slice")
    metric, potential = coefficients(a, psi, charge, mass_squared, quartic, mesh)
    r = vertical(psi, charge, mesh)@mesh.mean_zero
    inertia = r.T@metric@r
    coupling = r.T@metric@mesh.slice
    eta_map = -np.linalg.solve(inertia, coupling)
    gamma = mesh.slice.T@metric@mesh.slice+coupling.T@eta_map
    return gamma, potential, eta_map, inertia


def coulomb_representative(a, psi, charge=1.0, mesh=None):
    """The unique representative modulo mean-zero real gauge shifts."""
    mesh = geometry() if mesh is None else mesh
    a, psi = _configuration(a, psi)
    charge = _real_scalar(charge, "charge")
    b = mesh.mean_zero
    xi = -b@np.linalg.solve(b.T@mesh.d.T@mesh.mass@mesh.d@b,
                           b.T@mesh.d.T@mesh.mass@a)
    return a+mesh.d@xi, np.exp(1j*charge*xi)*psi, xi


def _state_arguments(q, sigma):
    if not np.isrealobj(q):
        raise ValueError("finite real 56-dimensional state coordinates required")
    q = np.asarray(q, dtype=float)
    if q.shape != (56,) or not np.isfinite(q).all():
        raise ValueError("finite real 56-dimensional state coordinates required")
    try:
        valid_sigma = (not isinstance(sigma, (bool, np.bool_)) and np.ndim(sigma) == 0
                       and np.isrealobj(sigma) and isfinite(sigma) and sigma > 0)
    except (TypeError, ValueError, OverflowError):
        valid_sigma = False
    if not valid_sigma:
        raise ValueError("positive finite real sigma required")
    return q, float(sigma)


def gaussian_half_density_log(q, sigma=1.0):
    """Log of the exactly normalized selected Gaussian in Lebesgue L2(R56).

    This is g=sqrt(rho)*f, not the curved-measure wavefunction f. Its squared
    modulus is the density of N(0, sigma**2 I56). No state evolution is run.
    """
    q, sigma = _state_arguments(q, sigma)
    return float(-14*(np.log(2*np.pi)+2*np.log(sigma))-(q/sigma)@(q/sigma)/4)


def gaussian_state_log_amplitude(q, sigma=1.0, charge=1.0, mesh=None):
    """Approximate log f=log g-log(det gamma)/4 at orthonormal slice q.

    Exact normalization belongs to the analytic state using the EXACT
    Riemannian density. Here gamma uses numerical element quadrature, so this
    pointwise amplitude is an approximation, not a normalization certificate.
    """
    q, sigma = _state_arguments(q, sigma)
    mesh = geometry() if mesh is None else mesh
    section = mesh.slice
    if (section.shape != (68, 56)
            or not np.allclose(section.T@section, np.eye(56), atol=1e-12, rtol=0)
            or not np.allclose(section[:42, 30:], 0, atol=1e-12, rtol=0)
            or not np.allclose(section[42:, :30], 0, atol=1e-12, rtol=0)
            or not np.allclose(section[42:, 30:], np.eye(26), atol=1e-12, rtol=0)
            or not np.allclose(mesh.d.T@mesh.mass@section[:42, :30], 0,
                               atol=1e-11, rtol=0)):
        raise ValueError("orthonormal Coulomb frame with standard scalar coordinates required")
    a = mesh.slice[:42, :30]@q[:30]
    psi = q[30:43]+1j*q[43:]
    gamma, _, _, _ = reduced_coefficients(a, psi, charge=charge, mesh=mesh)
    sign, logdet = np.linalg.slogdet(gamma)
    if sign <= 0 or not np.isfinite(logdet):
        raise ValueError("numerically positive reduced metric required")
    return float(gaussian_half_density_log(q, sigma)-logdet/4)


def gaussian_initial_moments(sigma, volume, transverse_stiffness_trace):
    """Exact algebraic selected-state moments, preserving rational inputs.

    Geometry/trace values supplied as floats retain their numerical precision;
    this function does not certify their geometric evaluation. Scalar real
    and imaginary coordinates each have variance sigma**2.
    """
    for name, value, positive in (("sigma", sigma, True), ("volume", volume, True),
                                 ("trace", transverse_stiffness_trace, False)):
        try:
            valid = (not isinstance(value, (bool, np.bool_)) and np.ndim(value) == 0
                     and np.isrealobj(value) and isfinite(value)
                     and (value > 0 if positive else value >= 0))
        except (TypeError, ValueError, OverflowError):
            valid = False
        if not valid:
            raise ValueError("finite real "+name+" with the required sign expected")
    return {
        "matter_l2": Fraction(4, 5)*sigma**2*volume,
        "matter_l4": Fraction(48, 35)*sigma**4*volume,
        "magnetic_energy": Fraction(1, 2)*sigma**2*transverse_stiffness_trace,
    }
