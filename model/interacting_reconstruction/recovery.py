"""Identify a finite local quartic action from observer-accessible readouts.

Only number-preparation labels, binary projector laws and calibrated clock
laws enter this inverse.  The oscillator basis, isolation protocols, real
generator hypothesis, norm caps and zero-or-separated-positive edge family
are supplied.  The bounds below control finite probe duration for ideal
probabilities; they do not certify source compilation, floating arithmetic,
clock calibration, occupation removal or an interacting continuum limit.
"""
from __future__ import annotations

from itertools import combinations
from math import sqrt
from numbers import Real

import numpy as np

from model.observer_reconstruction.recovery import clock_time, cube_geometry
from .bounds import richardson_radius


def _real(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite admitted real")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite admitted real") from exc
    if (not np.isfinite(result) or (positive and result <= 0)
            or (nonnegative and result < 0)):
        raise ValueError(f"{name} must be a finite admitted real")
    return result


def _probabilities(value, shape):
    if np.ma.isMaskedArray(value):
        raise ValueError("binary probabilities cannot be masked")
    try:
        raw = np.asarray(value)
    except (ValueError, TypeError) as exc:
        raise ValueError("binary probabilities must be a finite real array") from exc
    if (raw.shape != shape or raw.dtype.kind not in "fiu"
            or not np.isfinite(raw).all()):
        raise ValueError(f"binary probabilities must have shape {shape}")
    values = np.array(raw, dtype=float, copy=True)
    if (np.any(values < -1e-12) or np.any(values > 1+1e-12)
            or not np.allclose(values.sum(axis=-1), 1., atol=2e-10, rtol=0)):
        raise ValueError("both binary probabilities must be normalized and in [0,1]")
    return np.clip(values, 0., 1.)


def coefficient_from_records(probabilities, first_time, second_time, *,
                             norm_bound, probability_error=0.):
    """Recover one off-diagonal element of a real isolated generator.

    Array axes are (two increasing times, minus-i/plus-i projector, outcome).
    Outcome one projects onto (|a> +/- i|b>)/sqrt(2), after preparing |a>.
    Thus f(t)=p1(minus)-p1(plus) is odd and f'(0)=2 H[b,a].
    Richardson cancellation removes its cubic term without a gap or energy
    fit.  If ||H-cI|| <= B for some scalar c, |f^(5)| <= (2B)^5.

    For times u<v, the resulting coefficient error is at most
      (2B)^5 u^2 v^2 (u^2+v^2) / (240 (v^2-u^2)).
    This uses two separate Taylor remainders and is conservative.  Each
    supplied absolute probability error epsilon additionally contributes
      epsilon (v^2/u + u^2/v)/(v^2-u^2).
    Clock readings are treated as exact calibrated model-coordinate values.
    """
    p = _probabilities(probabilities, (2, 2, 2))
    u = _real(first_time, "first clock reading", positive=True)
    v = _real(second_time, "second clock reading", positive=True)
    bound = _real(norm_bound, "generator norm bound", positive=True)
    error = _real(probability_error, "probability error", nonnegative=True)
    if not u < v or error >= .5:
        raise ValueError("increasing clock readings and probability error below 1/2 are required")
    # Rational evaluation preserves the analytical bound for represented
    # inputs.  The inferred center still uses ordinary floating arithmetic.
    try:
        budget = richardson_radius(u, v, bound, probability_error=error)
        radius = budget["radius"]
        f = p[:, 0, 1]-p[:, 1, 1]
        element = (v*v*f[0]/(2*u)-u*u*f[1]/(2*v))/(v*v-u*u)
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("inverse coefficient exceeds the numerical representation") from exc
    if not np.isfinite((element, radius)).all():
        raise ValueError("inverse coefficient exceeds the numerical representation")
    return {**budget, "matrix_element": float(element),
            "lower": float(np.nextafter(element-radius, -np.inf)),
            "upper": float(np.nextafter(element+radius, np.inf)),
            "generator_norm_bound": bound}


def onsite_parameters(h02, h13, radius02, radius13, *, omega):
    """Resolve A and g in .5 p^2 + .5 A q^2 + g q^4.

    q^2 and q^4 are infinite oscillator polynomials compressed afterwards.
    With u=H20/sqrt(2), v=H31/sqrt(6), g=omega^2(v-u) and
    A=omega^2+omega(10u-6v).  A contains the incident bond contributions.
    Both transitions require at least four local occupation levels.
    """
    h02, h13 = _real(h02, "H20"), _real(h13, "H31")
    r02 = _real(radius02, "H20 radius", nonnegative=True)
    r13 = _real(radius13, "H31 radius", nonnegative=True)
    omega = _real(omega, "basis frequency", positive=True)
    u, v = h02/sqrt(2), h13/sqrt(6)
    # Strict rational upper approximations keep the propagation factors
    # conservative; the centers themselves have no floating certificate.
    ru, rv = r02*.707106781186548, r13*.408248290463864
    with np.errstate(over="ignore", invalid="ignore"):
        quartic = omega**2*(v-u)
        quadratic = omega**2+omega*(10*u-6*v)
        quartic_radius = np.nextafter(omega**2*(ru+rv), np.inf)
        quadratic_radius = np.nextafter(omega*(10*ru+6*rv), np.inf)
    if not np.isfinite((quartic, quadratic, quartic_radius, quadratic_radius)).all():
        raise ValueError("onsite action exceeds the numerical representation")
    return {"onsite_quadratic": float(quadratic),
            "onsite_quadratic_radius": float(quadratic_radius),
            "quartic": float(quartic), "quartic_radius": float(quartic_radius)}


def _common_interval(points, radii, name, *, positive=True):
    lower = float(np.max(points-radii))
    upper = float(np.min(points+radii))
    if lower > upper:
        raise ValueError(f"{name} has no common homogeneous interval")
    if positive and lower <= 0:
        raise ValueError(f"{name} is not resolved as strictly positive")
    center = float(np.clip(np.mean(points), lower, upper))
    return {"value": center, "lower": lower, "upper": upper,
            "radius": float(max(center-lower, upper-center))}


def _records(records):
    if not isinstance(records, dict):
        raise ValueError("calibration records must be a dictionary")
    expected = list(combinations(range(8), 2))
    try:
        pairs = np.asarray(records["pairs"])
        onsite = _probabilities(records["onsite_probabilities"], (8, 2, 2, 2))
        bonds = _probabilities(records["bond_probabilities"], (28, 2, 2))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("complete onsite and all-pair binary calibration records are required") from exc
    if pairs.dtype.kind not in "iu" or pairs.shape != (28, 2) or pairs.tolist() != [list(p) for p in expected]:
        raise ValueError("pairs must list every opaque channel pair once in lexicographic order")
    return onsite, bonds


def recover(records_first, records_second, clock_first, clock_second, *,
            omega, norm_bound_onsite, norm_bound_bond, probability_error=0.,
            edge_threshold=.1, require_positive_interaction=True):
    """Recover a positive homogeneous quartic cube action from probe records.

    A topology gap is an explicit action-family premise: each pair weight is
    zero or at least ``edge_threshold``.  An interval containing zero and
    lying below that threshold identifies absence only within this family.
    Intervals crossing the decision boundary are rejected.  All 28 pairs
    are assessed before any cube geometry or common coefficient is fitted.

    This returns the finite compressed action, its supplied-basis parameter
    intervals, a tested cube graph and a unit-edge rank-three Gram matrix.
    It does not identify the measured interacting mass gap with bare mass.
    """
    omega = _real(omega, "basis frequency", positive=True)
    if not isinstance(require_positive_interaction, (bool, np.bool_)):
        raise ValueError("require_positive_interaction must be boolean")
    threshold = _real(edge_threshold, "edge separation threshold", positive=True)
    first, bonds_first = _records(records_first)
    second, bonds_second = _records(records_second)
    times = (clock_time(clock_first), clock_time(clock_second))
    onsite_elements = np.empty((8, 2))
    onsite_radii = np.empty((8, 2))
    onsite_details = []
    for site in range(8):
        details = []
        for transition in range(2):
            result = coefficient_from_records(
                np.stack((first[site, transition], second[site, transition])),
                *times, norm_bound=norm_bound_onsite, probability_error=probability_error)
            onsite_elements[site, transition] = result["matrix_element"]
            onsite_radii[site, transition] = result["radius"]
            details.append(result)
        onsite_details.append(details)
    onsite = [onsite_parameters(*onsite_elements[i], *onsite_radii[i], omega=omega)
              for i in range(8)]
    a = np.array([item["onsite_quadratic"] for item in onsite])
    ar = np.array([item["onsite_quadratic_radius"] for item in onsite])
    g = np.array([item["quartic"] for item in onsite])
    gr = np.array([item["quartic_radius"] for item in onsite])
    adjacency = np.zeros((8, 8), dtype=int)
    weights, weight_radii = np.zeros((8, 8)), np.zeros((8, 8))
    raw_weights, raw_radii = np.zeros((8, 8)), np.zeros((8, 8))
    bond_details = []
    for index, (i, j) in enumerate(combinations(range(8), 2)):
        result = coefficient_from_records(
            np.stack((bonds_first[index], bonds_second[index])), *times,
            norm_bound=norm_bound_bond, probability_error=probability_error)
        weight = -2*omega*result["matrix_element"]
        radius = float(np.nextafter(2*omega*result["radius"], np.inf))
        lower, upper = weight-radius, weight+radius
        raw_weights[i, j] = raw_weights[j, i] = weight
        raw_radii[i, j] = raw_radii[j, i] = radius
        if lower >= threshold:
            adjacency[i, j] = adjacency[j, i] = 1
            weights[i, j] = weights[j, i] = weight
            weight_radii[i, j] = weight_radii[j, i] = radius
        elif not (lower <= 0 <= upper < threshold):
            raise ValueError("pair coupling is unresolved or incompatible with the supplied edge gap")
        bond_details.append(result)
    active = adjacency.astype(bool)
    if not active.any():
        raise ValueError("no spatial interaction edges were resolved")
    mass_squared = a-weights.sum(axis=1)
    mass_radius = ar+weight_radii.sum(axis=1)
    common_g = _common_interval(g, gr, "quartic interaction", positive=require_positive_interaction)
    if common_g["upper"] < 0:
        raise ValueError("quartic interaction is incompatible with the repulsive action family")
    interaction_resolved = common_g["lower"] > 0
    if not require_positive_interaction:
        # Intersect with the declared repulsive family for a matched free
        # control.  An interval touching zero is never a positive detection.
        common_g["lower"] = max(0., common_g["lower"])
        common_g["value"] = max(common_g["lower"], common_g["value"])
        common_g["radius"] = max(common_g["value"]-common_g["lower"],
                                 common_g["upper"]-common_g["value"])
    common_mass = _common_interval(mass_squared, mass_radius, "bare mass squared")
    common_kappa = _common_interval(weights[active], weight_radii[active], "spatial edge weight")
    laplacian = np.diag(adjacency.sum(axis=1))-adjacency
    stiffness = common_mass["value"]*np.eye(8)+common_kappa["value"]*laplacian
    geometry = cube_geometry(stiffness)
    if not np.array_equal(geometry["adjacency"], adjacency):
        raise ValueError("recovered geometry does not preserve the probed topology")
    # The geometry helper is shared with the free-field inverse.  Here its
    # mass and speed parameters describe the supplied bare action; neither
    # is an interacting pole mass, gap or measured signal velocity.
    geometry["bare_mass"] = geometry.pop("mass")
    geometry["bare_speed_per_unit_edge"] = geometry.pop("speed_per_unit_edge")
    geometry["bond_kappa"] = geometry.pop("coupling")
    return {**geometry, "bare_quadratic_stiffness": stiffness, "clock_readings": times,
            "basis_frequency": omega, "quartic": common_g["value"],
            "bare_mass_squared": common_mass["value"],
            "homogeneous_intervals": {"quartic": common_g,
                                      "bare_mass_squared": common_mass,
                                      "edge_weight": common_kappa},
            "onsite_quadratic": a, "onsite_quadratic_radius": ar,
            "onsite_quartic": g, "onsite_quartic_radius": gr,
            "site_bare_mass_squared": mass_squared, "site_bare_mass_squared_radius": mass_radius,
            "pair_weights": weights, "pair_weight_radii": weight_radii,
            "raw_pair_weights": raw_weights, "raw_pair_weight_radii": raw_radii,
            "onsite_coefficient_details": onsite_details, "bond_coefficient_details": bond_details,
            "edge_separation_threshold": threshold,
            "edge_gap_is_supplied": True, "isolated_generator_family_is_supplied": True,
            "positive_quartic_interaction_resolved": interaction_resolved,
            "finite_time_bounds_exclude_clock_and_floating_error": True,
            "bare_mass_is_not_measured_interacting_gap": True}
