"""Infer a finite free-field description from declared binary readouts alone.

The estimator imports no source executor, Hamiltonian or carrier coordinates.
Its inputs are pointer probabilities, intervention settings and an explicit
frequency alias bound.  Time is read in the dimensionless calibration clock's
units.  This is inverse identification within a supplied quadratic hypothesis;
it does not derive quantum structure, the chosen apparatus or a continuum.
"""
from __future__ import annotations

from math import pi
from numbers import Real

import numpy as np


def _real(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite admitted real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} exceeds the numerical representation") from exc
    if (not np.isfinite(result) or (positive and result <= 0)
            or (nonnegative and result < 0)):
        raise ValueError(f"{name} must be a finite admitted real number")
    return result


def _array(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} must not be masked")
    try:
        raw = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real array") from exc
    if raw.dtype.kind not in "fiu" or not np.all(np.isfinite(raw)):
        raise ValueError(f"{name} must be a finite real array")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        result = np.array(raw, dtype=float, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} exceeds the numerical representation")
    return result


def characteristic(probabilities):
    """Read chi from bit-zero probabilities at phases (0, pi/2).

    For K0=cos((eta Q+theta)/2), 2 p0(theta)-1 equals
    Re(exp(i theta) chi).  Both binary outcomes belong to each experiment;
    these probabilities are not conditional on a preferred field outcome.

    Floating source contractions may cross zero or one by at most 1e-12.
    Only that fixed endpoint dust is clamped in a copy; raw source receipts
    remain unchanged. This admission tolerance is not an error certificate.
    """
    p = _array(probabilities, "Ramsey probabilities")
    if p.ndim < 1 or p.shape[-1] != 2 or np.any((p < -1e-12) | (p > 1+1e-12)):
        raise ValueError("two phase probabilities in [0,1], up to endpoint roundoff, are required")
    p = np.clip(p, 0., 1.)
    chi = 2*p[..., 0]-1 + 1j*(1-2*p[..., 1])
    if np.any(np.abs(chi) > 1+1e-10):
        raise ValueError("phase probabilities are incompatible with a characteristic function")
    return chi


def clock_time(probabilities, *, visibility_tolerance=1e-8):
    """Read the clock's positive phase on the declared open window (0, pi).

    Clock calibrations are p0(0)=(1+cos(t))/2 and
    p0(pi/2)=(1-sin(t))/2.  The phase window is a protocol premise:
    wrapped observations alone cannot reveal an unbounded cycle count.
    """
    tolerance = _real(visibility_tolerance, "clock visibility tolerance", positive=True)
    if tolerance >= 1:
        raise ValueError("clock visibility tolerance must be less than one")
    p = _array(probabilities, "clock probabilities")
    if p.shape != (2,):
        raise ValueError("one pair of clock phase probabilities is required")
    value = complex(characteristic(p))
    if abs(abs(value)-1) > tolerance:
        raise ValueError("clock probabilities violate the unit-visibility calibration")
    phase = float(np.angle(value))
    if not 0 < phase < pi:
        raise ValueError("clock phase must lie in the declared positive window (0, pi)")
    return phase


def response_from_records(probabilities, *, gamma, eta, minimum_visibility=.05,
                          phase_margin=.05):
    """Return R[read,kick] from (+gamma,-gamma) complete Ramsey experiments.

    Input axes are (kick sign, kicked site, read site, Ramsey phase).
    Opposite kicks remove the zero-mean Gaussian visibility.  The declared
    phase window excludes observed branch-cut crossings; no finite sample
    can independently rule out an additional unseen phase winding.
    """
    gamma = _real(gamma, "kick strength", positive=True)
    eta = _real(eta, "detector strength", positive=True)
    visibility = _real(minimum_visibility, "minimum visibility", positive=True)
    margin = _real(phase_margin, "phase margin", positive=True)
    if visibility >= 1 or margin >= pi:
        raise ValueError("visibility and phase margin must be inside their admitted windows")
    p = _array(probabilities, "field probabilities")
    if p.ndim != 4 or p.shape[0] != 2 or p.shape[-1] != 2 or p.shape[1] != p.shape[2]:
        raise ValueError("field records need shape (2, sites, sites, 2)")
    if not 1 <= p.shape[1] <= 64:
        raise ValueError("one to sixty-four field channels are required")
    chi = characteristic(p)
    if np.any(np.abs(chi) < visibility):
        raise ValueError("field phase cannot be resolved at the observed visibility")
    phases = np.angle(chi)
    if np.any(np.abs(phases) >= pi-margin):
        raise ValueError("field phase is outside the declared unwrapped window")
    scale = 2*gamma*eta
    if not np.isfinite(scale) or scale == 0:
        raise ValueError("intervention scale cannot be represented")
    response = ((phases[0]-phases[1])/scale).T
    if not np.all(np.isfinite(response)):
        raise ValueError("response exceeds the numerical representation")
    return response


def _stiffness(value):
    matrix = _array(value, "stiffness")
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]
            or not 1 <= len(matrix) <= 64
            or not np.allclose(matrix, matrix.T, atol=1e-10, rtol=1e-10)):
        raise ValueError("a bounded real symmetric stiffness matrix is required")
    values, basis = np.linalg.eigh((matrix+matrix.T)/2)
    if values[0] <= 0:
        raise ValueError("a strictly positive stiffness matrix is required")
    return matrix, values, basis


def infer_stiffness(first, second, first_time, second_time, *,
                    frequency_upper_bound, relative_tolerance=2e-4,
                    maximum_condition=1e6):
    """Invert R(2t)=2 cos(t sqrt(K)) R(t) in an explicit nonaliasing band.

    For a free canonical field R(t)=K^(-1/2) sin(t sqrt(K)).  The supplied
    bound must obey t*frequency_upper_bound < pi.  It is necessary for a
    unique branch of matrix arccos and is not inferred from these two times.
    Residual checks test this quadratic model; they are floating diagnostics,
    not propagated confidence intervals or outward error certificates.
    """
    r1, r2 = _array(first, "first response"), _array(second, "second response")
    if (r1.ndim != 2 or r1.shape != r2.shape or r1.shape[0] != r1.shape[1]
            or not 1 <= len(r1) <= 64):
        raise ValueError("two square response matrices of the same size are required")
    time = _real(first_time, "first clock reading", positive=True)
    later = _real(second_time, "second clock reading", positive=True)
    upper = _real(frequency_upper_bound, "frequency upper bound", positive=True)
    tolerance = _real(relative_tolerance, "relative tolerance", positive=True)
    condition_limit = _real(maximum_condition, "condition limit", positive=True)
    if tolerance >= .1 or condition_limit < 1:
        raise ValueError("tolerance must be below .1 and condition limit at least one")
    if not np.isclose(later, 2*time, rtol=tolerance, atol=1e-12):
        raise ValueError("observed clock readings must have ratio two")
    if not np.isfinite(time*upper) or time*upper >= pi:
        raise ValueError("frequency upper bound violates the declared nonaliasing band")
    condition = float(np.linalg.cond(r1))
    if not np.isfinite(condition) or condition > condition_limit:
        raise ValueError("first response is singular or too ill-conditioned")
    response_scale = max(float(np.linalg.norm(r1, 2)), np.finfo(float).tiny)
    symmetry = max(float(np.linalg.norm(r-r.T, 2)) for r in (r1, r2))/response_scale
    if symmetry > tolerance:
        raise ValueError("responses violate reciprocity of the quadratic hypothesis")
    cosine = .5*np.linalg.solve(r1.T, r2.T).T
    cosine_asymmetry = float(np.linalg.norm(cosine-cosine.T, 2))
    if cosine_asymmetry > tolerance:
        raise ValueError("response ratio has no symmetric free-field generator")
    cosine = (cosine+cosine.T)/2
    values, basis = np.linalg.eigh(cosine)
    # Do not clip an inadmissible response onto an apparently physical model.
    if values[0] <= -1 or values[-1] >= 1:
        raise ValueError("response ratio has no strictly positive nonaliased spectrum")
    frequencies = np.arccos(values)/time
    if frequencies.max() > upper*(1+tolerance):
        raise ValueError("recovered spectrum exceeds the supplied frequency band")
    matrix = (basis*frequencies**2)@basis.T
    predicted_first = (basis*(np.sin(time*frequencies)/frequencies))@basis.T
    residual = float(np.linalg.norm(predicted_first-r1, 2))/response_scale
    if residual > tolerance:
        raise ValueError("response amplitude is inconsistent with canonical quadratic dynamics")
    return {"stiffness": matrix, "frequencies": np.sort(frequencies),
            "response_condition": condition, "response_reciprocity_residual": symmetry,
            "cosine_asymmetry": cosine_asymmetry,
            "canonical_response_relative_residual": residual,
            "frequency_upper_bound": upper,
            "frequency_alias_bound_is_supplied": True}


def cube_geometry(stiffness, *, relative_tolerance=2e-4):
    """Recover an unlabeled cube and its Gram matrix from eight field channels.

    Adjacency is inferred from the negative off-diagonal stiffness entries.
    Equal edge weights and the massive graph-Laplacian form are tested, not
    installed in the recovered matrix.  Edge length one is a declared unit;
    coordinates are determined only up to a Euclidean isometry.
    """
    matrix, _, _ = _stiffness(stiffness)
    tolerance = _real(relative_tolerance, "relative tolerance", positive=True)
    if matrix.shape != (8, 8) or tolerance >= .1:
        raise ValueError("cube recovery needs eight channels and tolerance below .1")
    weights = -matrix.copy()
    np.fill_diagonal(weights, -np.inf)
    neighbors = np.argsort(weights, axis=1)[:, -3:]
    adjacency = np.zeros((8, 8), dtype=int)
    adjacency[np.arange(8)[:, None], neighbors] = 1
    if not np.array_equal(adjacency, adjacency.T):
        raise ValueError("three strongest couplings do not define a reciprocal cube")
    coupling = float(np.mean(-matrix[adjacency.astype(bool)]))
    mass_squared = float(np.mean(matrix.sum(axis=1)))
    if coupling <= 0 or mass_squared <= 0:
        raise ValueError("positive mass and edge coupling are required")
    laplacian = 3*np.eye(8)-adjacency
    fitted = mass_squared*np.eye(8)+coupling*laplacian
    fit_residual = float(np.linalg.norm(matrix-fitted, 2))/max(float(np.linalg.norm(matrix, 2)), 1.)
    if fit_residual > tolerance:
        raise ValueError("stiffness does not have a common massive cube-Laplacian form")
    distances = np.where(adjacency, 1., np.inf)
    np.fill_diagonal(distances, 0.)
    for intermediate in range(8):
        distances = np.minimum(distances, distances[:, intermediate, None]+distances[None, intermediate, :])
    if (not np.all(np.isfinite(distances))
            or any(sorted(row.tolist()) != [0., 1., 1., 1., 2., 2., 2., 3.]
                   for row in distances)):
        raise ValueError("the inferred graph is not a connected cube")
    centering = np.eye(8)-np.ones((8, 8))/8
    # On a cube the graph distance equals the squared Euclidean distance
    # in the unit-edge embedding, not the Euclidean distance itself.
    gram = -.5*centering@distances@centering
    eigenvalues, basis = np.linalg.eigh(gram)
    if eigenvalues[0] < -1e-10 or np.count_nonzero(eigenvalues > 1e-10) != 3:
        raise ValueError("graph distances do not admit the declared rank-three cube embedding")
    coordinates = basis[:, -3:]*np.sqrt(eigenvalues[-3:])
    return {"adjacency": adjacency, "graph_distances": distances.astype(int),
            "gram": gram, "coordinates": coordinates,
            "spatial_dimension": 3, "mass": float(np.sqrt(mass_squared)),
            "coupling": coupling, "speed_per_unit_edge": float(np.sqrt(coupling)),
            "cube_fit_relative_residual": fit_residual,
            "edge_length_unit_is_supplied": True}


def recover(records_at_tau, records_at_2tau, clock_tau, clock_2tau, *,
            gamma, eta, frequency_upper_bound, relative_tolerance=2e-4):
    """Recover clock units, quadratic generator and cube from readout records."""
    first_time, second_time = clock_time(clock_tau), clock_time(clock_2tau)
    first = response_from_records(records_at_tau, gamma=gamma, eta=eta)
    second = response_from_records(records_at_2tau, gamma=gamma, eta=eta)
    inferred = infer_stiffness(first, second, first_time, second_time,
                               frequency_upper_bound=frequency_upper_bound,
                               relative_tolerance=relative_tolerance)
    geometry = cube_geometry(inferred["stiffness"], relative_tolerance=relative_tolerance)
    return {**inferred, **geometry, "clock_readings": (first_time, second_time),
            "response_at_tau": first, "response_at_2tau": second}


def predict_chi(stiffness, time, kick_site, read_site, *, gamma, eta, temperature=0):
    """Predict a withheld Ramsey characteristic within the fitted free model.

    The Gaussian vacuum (temperature zero) or thermal preparation is an
    additional state hypothesis.  Susceptibility alone neither determines
    this covariance nor establishes quantum mechanics.  Physical time and
    energy units are not supplied by this formula.
    """
    matrix, values, basis = _stiffness(stiffness)
    time = _real(time, "clock interval", nonnegative=True)
    gamma, eta = _real(gamma, "kick strength"), _real(eta, "detector strength")
    temperature = _real(temperature, "temperature", nonnegative=True)
    for index in (kick_site, read_site):
        if isinstance(index, (bool, np.bool_)) or not isinstance(index, (int, np.integer)) or not 0 <= index < len(matrix):
            raise ValueError("kick and read addresses must be admitted integer channel indices")
    frequencies = np.sqrt(values)
    thermal = np.ones_like(frequencies) if temperature == 0 else 1/np.tanh(frequencies/(2*temperature))
    variance = float(np.sum(basis[read_site]**2*thermal/(2*frequencies)))
    response = float(np.sum(basis[read_site]*basis[kick_site]*np.sin(time*frequencies)/frequencies))
    try:
        exponent = -.5*eta**2*variance+1j*eta*gamma*response
    except OverflowError as exc:
        raise ValueError("prediction exceeds the numerical representation") from exc
    if not np.isfinite(exponent):
        raise ValueError("prediction exceeds the numerical representation")
    return complex(np.exp(exponent))


def free_wightman(stiffness, time):
    """Vacuum two-point function under the independently declared CCR model.

    W_xy(t)=<q_x(t)q_y(0)>=[exp(-it sqrt(K))/(2 sqrt(K))]_xy.
    The imaginary commutator part is -R(t)/2.  This quantum prediction is
    not inferred from a retarded response alone, and is not a measured
    continuum Wightman distribution.
    """
    _, values, basis = _stiffness(stiffness)
    time = _real(time, "clock interval")
    frequencies = np.sqrt(values)
    return (basis*(np.exp(-1j*time*frequencies)/(2*frequencies)))@basis.T
