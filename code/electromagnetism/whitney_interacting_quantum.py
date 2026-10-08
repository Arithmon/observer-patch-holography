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
from decimal import Decimal
from functools import lru_cache
from fractions import Fraction
from itertools import combinations, product
from math import frexp, fsum, isfinite, ldexp

import numpy as np
from scipy.linalg import null_space, qr, solve_triangular

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
    if np.ma.isMaskedArray(a) or np.ma.isMaskedArray(psi):
        raise ValueError("masked configuration coefficients are not supplied values")
    a, psi = np.asarray(a, dtype=object), np.asarray(psi, dtype=object)
    if a.shape != (42,) or psi.shape != (13,):
        raise ValueError("finite 42-edge and 13-complex-node coefficients required")
    try:
        real_a = np.array([_real_scalar(x, "edge coefficient") for x in a])
        complex_psi = []
        for value in psi:
            if isinstance(value, (complex, np.complexfloating)):
                real, imag = value.real, value.imag
            else:
                real, imag = value, 0
            complex_psi.append(complex(_real_scalar(real, "scalar coefficient"),
                                       _real_scalar(imag, "scalar coefficient")))
    except ValueError as exc:
        raise ValueError("finite 42-edge and 13-complex-node coefficients required without input precision loss") from exc
    return real_a, np.array(complex_psi)


def _real_scalar(value, name):
    if (np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, float, np.integer, np.floating, Fraction, Decimal))):
        raise ValueError("finite real "+name+" required")
    try:
        exact = (Fraction(int(value)) if isinstance(value, (int, np.integer))
                 else Fraction(value) if isinstance(value, Fraction)
                 else Fraction(*value.as_integer_ratio()))
        result = float(exact)
        if not isfinite(result) or Fraction(result) != exact:
            raise ValueError("original scalar is not representable in binary64")
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError("finite real "+name+" required without input precision loss") from exc
    return result


def scalar_fields(a, psi, charge, mesh, *, include_covariant=True):
    """Value, Jacobian and optional covariant spatial gradient.

    Density-only callers do not evaluate unused spatial products, which may
    exceed the reporting range even when their kinetic factors are resolved.
    """
    a, psi = _configuration(a, psi)
    charge = _real_scalar(charge, "charge")
    phase = (np.exp(1j*charge*np.einsum("qie,e->qi", mesh.paths, a)) if charge
             else np.ones(mesh.nodal.shape, dtype=complex))
    w = mesh.nodal*phase
    value = w@psi
    # Pair the two endpoints before multiplying by barycentric weights.
    # Uniform matter at a=0 has exactly zero edge derivative; summing the
    # two large endpoint contributions separately loses that identity.
    left, right = np.asarray(mesh.edges).T
    transported = phase*psi
    dressing = ((1j*charge*mesh.nodal[:, left]*mesh.nodal[:, right]
                 * (transported[:, left]-transported[:, right])) if charge
                else np.zeros((len(mesh.weights), 42), dtype=complex))
    jacobian = np.column_stack((dressing, w, 1j*w))
    if not include_covariant:
        return value, jacobian, None
    # Partition of unity removes the constant mode before differentiation.
    # Node zero belongs to every tetrahedron of this cone.
    grad = np.einsum("qi,qic->qc", transported-transported[:, :1], mesh.nodal_grad)
    if charge:
        grad_phase = np.einsum("qiec,e->qic", mesh.path_grad, charge*a)
        grad += 1j*np.einsum("qi,i,qic->qc", w, psi, grad_phase)
        potential = np.einsum("qec,e->qc", mesh.edge_forms, charge*a)
        grad -= 1j*potential*value[:, None]
    covariant = grad
    return value, jacobian, covariant


def _potential_parameters(mass_squared, quartic):
    try:
        mass_squared = _real_scalar(mass_squared, "mass-squared")
        quartic = _real_scalar(quartic, "quartic coupling")
        if mass_squared < 0 or quartic < 0:
            raise ValueError("negative potential coupling")
    except ValueError as exc:
        raise ValueError("nonnegative finite mass-squared and quartic coupling required without input precision loss") from exc
    return mass_squared, quartic


def _positive_product_sum(groups):
    """Sum nonnegative products without first overflowing/underflowing powers.

    A group consists of (factor, integer power) pairs. Mantissas stay near
    unity until the final sum; the powers of two are tracked separately.
    This evaluates the supplied quadrature, not a bound on integration error.
    """
    mantissas, exponents = [], []
    for group in groups:
        arrays = np.broadcast_arrays(*(np.asarray(value, dtype=float) for value, _ in group))
        mantissa = np.ones(arrays[0].shape)
        exponent = np.zeros(arrays[0].shape, dtype=np.int64)
        for array, (_, power) in zip(arrays, group, strict=True):
            if not np.isfinite(array).all() or (power % 2 and np.any(array < 0)):
                raise ValueError("potential precision requires finite positive factors")
            fraction, binary_exponent = np.frexp(abs(array))
            mantissa *= fraction**power
            exponent += binary_exponent*power
        nonzero = mantissa != 0
        mantissas.extend(mantissa[nonzero].flat)
        exponents.extend(exponent[nonzero].flat)
    if not mantissas:
        return 0.
    maximum = int(max(exponents))
    # Terms too far below the largest term cannot change the reported sum.
    # Their underflow here is relative to that sum, not loss of its scale.
    with np.errstate(under="ignore"):
        scaled = np.ldexp(np.asarray(mantissas), np.asarray(exponents)-maximum)
    mantissa, adjustment = frexp(fsum(scaled))
    exponent = maximum+adjustment
    try:
        result = ldexp(mantissa, exponent)
        recovered = ldexp(result, -exponent)
    except OverflowError as exc:
        raise ValueError("potential precision outside finite reporting range") from exc
    if not isfinite(result) or result <= 0 or abs(recovered/mantissa-1) > 1e-12:
        raise ValueError("potential precision outside reliable reporting range")
    return result


def _magnetic_energy_terms(a, mesh):
    """Positive curl factors from original oriented face circulations.

    Three original binary64 edge values determine each circulation. Exact
    rational addition retains a small curl behind a large pure gradient;
    it is not a zero floor on the already damaged stiffness quadratic.
    """
    if not np.any(a):
        return []
    oriented = {}
    for index, (left, right) in enumerate(mesh.edges):
        oriented[left, right] = Fraction(float(a[index]))
        oriented[right, left] = -oriented[left, right]
    groups = []
    for tet in mesh.tetrahedra:
        pairs = list(combinations(range(1, 4), 2))
        try:
            circulation = np.array([float(oriented[tet[0], tet[i]]+oriented[tet[i], tet[j]]
                                          -oriented[tet[0], tet[j]]) for i, j in pairs])
        except OverflowError as exc:
            raise ValueError("magnetic potential precision outside finite reporting range") from exc
        scale = np.max(abs(circulation))
        if scale:
            points = mesh.vertices[list(tet)]
            gradients = np.linalg.inv(np.column_stack((np.ones(4), points)))[1:].T
            volume = abs(np.linalg.det(points[1:]-points[0]))/6
            curls = np.array([2*np.cross(gradients[i], gradients[j]) for i, j in pairs])
            normalized = (circulation/scale)@curls
            groups.append(((volume/2, 1), (scale, 2), (normalized, 2)))
    return groups


def _potential_energy(a, value, covariant, mass_squared, quartic, mesh):
    groups = _magnetic_energy_terms(a, mesh)
    for component in (covariant.real, covariant.imag):
        groups.append(((mesh.weights[:, None], 1), (component, 2)))
    if mass_squared:
        for component in (value.real, value.imag):
            groups.append(((mesh.weights, 1), (mass_squared, 1), (component, 2)))
    if quartic:
        for component in (value.real, value.imag):
            groups.append(((mesh.weights, 1), (quartic, 1), (.5, 1), (component, 4)))
        groups.append(((mesh.weights, 1), (quartic, 1), (value.real, 2), (value.imag, 2)))
    return _positive_product_sum(groups)


def coefficients(a, psi, charge=1.0, mass_squared=1.0, quartic=1.0, mesh=None):
    """Full 68-real-coordinate kinetic G and the coupled potential V."""
    mesh = geometry() if mesh is None else mesh
    a, psi = _configuration(a, psi)
    mass_squared, quartic = _potential_parameters(mass_squared, quartic)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            value, jacobian, covariant = scalar_fields(a, psi, charge, mesh)
            weighted = np.sqrt(mesh.weights)[:, None]*jacobian
            metric = 2*np.real(weighted.conj().T@weighted)
            metric[:42, :42] += mesh.mass
            potential = _potential_energy(a, value, covariant, mass_squared, quartic, mesh)
    except FloatingPointError as exc:
        raise ValueError("coefficient precision outside finite reporting range") from exc
    if not np.isfinite(metric).all():
        raise ValueError("coefficient precision outside finite reporting range")
    return metric, potential


def vertical(psi, charge, mesh):
    """All thirteen infinitesimal gauge columns in [a, Re psi, Im psi]."""
    _, psi = _configuration(np.zeros(42), psi)
    charge = _real_scalar(charge, "charge")
    if psi.shape != (13,) or not np.isfinite(psi).all():
        raise ValueError("finite 13-complex-node coefficients required")
    return np.vstack((mesh.d, -charge*np.diag(psi.imag), charge*np.diag(psi.real)))


KINETIC_RESOLUTION = 1e-7


@dataclass(frozen=True)
class KineticReduction:
    """Numerical square root of the reduced positive quadratic form.

    The resolution gates concern binary64 algebra on the supplied element
    factors. They are not interval bounds on quadrature or geometry error.
    A usable square root need not admit a resolved dense Gram matrix.
    """
    factor: np.ndarray
    potential: float | None
    eta_map: np.ndarray
    inertia: np.ndarray
    singular_values: np.ndarray
    resolution: float

    def logdet(self):
        return float(2*np.log(abs(np.diag(self.factor))).sum())

    def dense_metric(self):
        absolute_gram = abs(self.factor).T@abs(self.factor)
        smallest = self.singular_values[-1]
        gram_resolution = (56*np.finfo(float).eps
                           * (np.linalg.norm(absolute_gram, 2)/smallest)/smallest)
        if not np.isfinite(gram_resolution) or gram_resolution > KINETIC_RESOLUTION:
            raise ValueError("kinetic precision insufficient for a dense reduced metric")
        gamma = self.factor.T@self.factor
        try:
            np.linalg.cholesky(gamma)
        except np.linalg.LinAlgError as exc:
            raise ValueError("kinetic precision does not resolve a positive dense metric") from exc
        return gamma


def reduced_kinetic(a, psi, charge=1.0, mass_squared=1.0, quartic=1.0, mesh=None,
                    *, include_potential=True):
    """Minimize the original positive action using scaled Householder QR.

    If [V,H]=Q R are vertical and slice kinetic factors, the trailing block
    of R is the reduced factor and the first block determines the minimizer.
    No nearly equal Gram matrices are subtracted. The direct vertical scalar
    factor is the Ward identity i*e*Psi*lambda, before any rounded J@R product.
    Density-only callers may omit the potential, returned as None: an unused
    quartic energy can overflow even when the kinetic density is resolved.
    """
    mesh = geometry() if mesh is None else mesh
    a, psi = _configuration(a, psi)
    charge = _real_scalar(charge, "charge")
    mass_squared, quartic = _potential_parameters(mass_squared, quartic)
    if type(include_potential) is not bool:
        raise ValueError("include_potential must be a Boolean")
    # Compare in common units before any norm squares can overflow. This is
    # the same absolute-plus-relative slice tolerance at every finite scale.
    a_scale = max(1., float(np.max(abs(a))))
    normalized_a = a/a_scale
    if np.max(abs(mesh.d.T@mesh.mass@normalized_a)) > 1e-9*(1/a_scale+np.linalg.norm(normalized_a)):
        raise ValueError("configuration must lie on the Coulomb slice")
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            value, jacobian, covariant = scalar_fields(a, psi, charge, mesh, include_covariant=include_potential)
            electric = np.linalg.cholesky(mesh.mass).T
            weights = np.sqrt(2*mesh.weights)[:, None]
            scalar_h = weights*(jacobian@mesh.slice)
            scalar_v = weights*(1j*charge*value[:, None]*(mesh.nodal@mesh.mean_zero))
            h = np.vstack((electric@mesh.slice[:42], scalar_h.real, scalar_h.imag))
            v = np.vstack((electric@mesh.d@mesh.mean_zero, scalar_v.real, scalar_v.imag))
            joint = np.column_stack((v, h))
            scales = np.linalg.norm(joint, axis=0)
            if not np.isfinite(joint).all() or np.any(scales <= 0) or not np.isfinite(scales).all():
                raise ValueError("kinetic precision cannot resolve the source factors")
            normalized = joint/scales
            vertical_spectrum = np.linalg.svd(normalized[:, :12], compute_uv=False)
            triangular = qr(normalized, mode="r", check_finite=False)[0][:68]
            eta_map = (-solve_triangular(triangular[:12, :12], triangular[:12, 12:])
                       * scales[12:][None, :]/scales[:12, None])
            factor = triangular[12:, 12:]*scales[12:][None, :]
            singular = np.linalg.svd(factor, compute_uv=False)
            resolution = (68*np.finfo(float).eps*np.linalg.norm(h, 2)/singular[-1]
                          * vertical_spectrum[0]/vertical_spectrum[-1])
            if not np.isfinite(resolution) or resolution > KINETIC_RESOLUTION:
                raise ValueError("kinetic precision insufficient for the reduced factor")
            inertia = v.T@v
            potential = None
            if include_potential:
                potential = _potential_energy(a, value, covariant, mass_squared, quartic, mesh)
            if ((potential is not None and not np.isfinite(potential))
                    or not np.isfinite(eta_map).all() or not np.isfinite(inertia).all()):
                raise ValueError("kinetic precision outside finite reporting range")
    except (FloatingPointError, np.linalg.LinAlgError) as exc:
        raise ValueError("kinetic precision insufficient for the supplied configuration") from exc
    return KineticReduction(factor, potential, eta_map, inertia, singular, float(resolution))


def reduced_coefficients(a, psi, charge=1.0, mass_squared=1.0, quartic=1.0, mesh=None):
    """Return resolved dense gamma, V, eta_map and gauge inertia.

    The minimizing potential is mean_zero @ eta_map @ velocity. A finite,
    positive action whose weak modes cannot survive dense binary64 assembly
    is explicitly refused; density-only callers can use the square-root path.
    """
    reduction = reduced_kinetic(a, psi, charge, mass_squared, quartic, mesh)
    return reduction.dense_metric(), reduction.potential, reduction.eta_map, reduction.inertia


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
    if not np.isrealobj(q) or np.ma.isMaskedArray(q):
        raise ValueError("finite real 56-dimensional state coordinates required")
    q = np.asarray(q, dtype=object)
    if q.shape != (56,):
        raise ValueError("finite real 56-dimensional state coordinates required")
    try:
        q = np.array([_real_scalar(value, "state coordinate") for value in q])
    except ValueError as exc:
        raise ValueError("finite real 56-dimensional state coordinates required without input precision loss") from exc
    sigma = _real_scalar(sigma, "sigma")
    if sigma <= 0:
        raise ValueError("positive finite real sigma required")
    return q, float(sigma)


def gaussian_half_density_log(q, sigma=1.0):
    """Log of the exactly normalized selected Gaussian in Lebesgue L2(R56).

    This is g=sqrt(rho)*f, not the curved-measure wavefunction f. Its squared
    modulus is the density of N(0, sigma**2 I56). No state evolution is run.
    """
    q, sigma = _state_arguments(q, sigma)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            half_scaled = (q/sigma)/2
            result = float(-14*(np.log(2*np.pi)+2*np.log(sigma))-half_scaled@half_scaled)
    except FloatingPointError as exc:
        raise ValueError("Gaussian log amplitude outside finite reporting range") from exc
    if not isfinite(result):
        raise ValueError("Gaussian log amplitude outside finite reporting range")
    return result


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
    logdet = reduced_kinetic(a, psi, charge=charge, mesh=mesh, include_potential=False).logdet()
    result = float(gaussian_half_density_log(q, sigma)-logdet/4)
    if not isfinite(result):
        raise ValueError("Gaussian state log amplitude outside finite reporting range")
    return result


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
