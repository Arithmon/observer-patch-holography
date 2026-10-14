"""Exact oscillator moments and short-time comparison certificates.

This extends the factorial-weighted Bargmann polynomial representation and
spectral Taylor remainder used in ``model.common_source.infinite_history``
to the declared eight-site scalar cube.  Polynomials have finite support;
the Hamiltonian acting on them has no artificial occupation cutoff.
Here m=c=a=1, lambda is zero or one, and the oscillator basis has Omega=2.
Every duration is a supplied Hamiltonian model coordinate. Certificates
concern ideal arithmetic, not physical service errors.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt, lcm, prod


SITES = tuple((x, y, z) for x in range(2) for y in range(2) for z in range(2))
EDGES = tuple((i, j) for i, x in enumerate(SITES) for j, y in enumerate(SITES)
              if i < j and sum(abs(a-b) for a, b in zip(x, y)) == 1)
VACUUM = (0,) * 8
SparseState = dict[tuple[int, ...], F]


def _fraction(value, name, *, nonnegative=True):
    if type(value) not in (int, F):
        raise ValueError(f"{name} must be an exact rational")
    value = F(value)
    if nonnegative and value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _coupling(value):
    value = _fraction(value, "coupling")
    if value not in (F(0), F(1)):
        raise ValueError("the declared coupling is exactly zero or one")
    return int(value)


def cutoffs(local_cutoff=None):
    """Normalize a uniform cutoff or eight local dimensions; None is infinite."""
    if local_cutoff is None:
        return None
    if type(local_cutoff) is int:
        local_cutoff = (local_cutoff,) * 8
    if (not isinstance(local_cutoff, (tuple, list)) or len(local_cutoff) != 8
            or any(type(d) is not int or d < 1 for d in local_cutoff)):
        raise ValueError("cutoff must be a positive integer or eight positive integers")
    return tuple(local_cutoff)


def tensor_dimension(local_cutoff):
    """Dimension of the occupation tensor box, without constructing its basis."""
    dims = cutoffs(local_cutoff)
    if dims is None:
        raise ValueError("a finite occupation cutoff is required")
    return prod(dims)


def vacuum():
    return {VACUUM: F(1)}


def _validate_state(state):
    if not isinstance(state, dict):
        raise ValueError("state must map eight-component occupation tuples to rationals")
    for key, value in state.items():
        if (not isinstance(key, tuple) or len(key) != 8
                or any(type(n) is not int or n < 0 for n in key)):
            raise ValueError("invalid occupation tuple")
        _fraction(value, "coefficient", nonnegative=False)


def project(state, local_cutoff):
    _validate_state(state)
    dims = cutoffs(local_cutoff)
    return {key: F(value) for key, value in state.items() if value
            and (dims is None or all(n < d for n, d in zip(key, dims)))}


@lru_cache(maxsize=None)
def _weight(key):
    return prod(factorial(n) for n in key)


def norm_squared(state):
    """Bargmann norm: z^n denotes (a†)^n|0>, with squared norm prod(n!)."""
    _validate_state(state)
    return sum((F(c) ** 2 * _weight(key) for key, c in state.items()), F(0))


def sqrt_upper(value, digits=18):
    """Outward-rounded rational upper enclosure, verified by exact squaring."""
    value = _fraction(value, "squared norm")
    if type(digits) is not int or not 0 <= digits <= 120:
        raise ValueError("root precision must be an integer between zero and 120")
    scale = 10 ** digits
    numerator = isqrt(value.numerator * scale * scale // value.denominator)
    if F(numerator, scale) ** 2 < value:
        numerator += 1
    return F(numerator, scale)


def _add(out, key, value, dims):
    if not value or (dims is not None and any(n >= d for n, d in zip(key, dims))):
        return
    total = out.get(key, 0) + value
    if total:
        out[key] = total
    else:
        out.pop(key, None)


def _numerator(vector, coupling, dims=None, center=True):
    """Apply 384 times H, or H-<vac|H|vac>, to integer coefficients.

    H = 2N+8 - sum_edges Bx By/4 + coupling sum_x Bx^4/384,
    B=a+a†.  A finite dims tuple projects only AFTER each complete H;
    in particular it never replaces P B^4 P with (P B P)^4.
    """
    out = {}
    for key, coefficient in vector.items():
        diagonal = 768 * sum(key) + 6 * coupling * sum(n*(n+1) for n in key)
        if not center:
            diagonal += 3072 + 24 * coupling
        _add(out, key, coefficient * diagonal, dims)
        for x, y in EDGES:
            for dx, fx in ((1, 1), (-1, key[x])):
                for dy, fy in ((1, 1), (-1, key[y])):
                    if fx and fy:
                        dest = list(key)
                        dest[x] += dx
                        dest[y] += dy
                        _add(out, tuple(dest), -96 * coefficient * fx * fy, dims)
        if coupling:
            for site, n in enumerate(key):
                terms = ((4, 1), (2, 4*n+6),
                         (-2, n*(n-1)*(4*n-2)),
                         (-4, n*(n-1)*(n-2)*(n-3)))
                for delta, factor in terms:
                    if factor and n+delta >= 0:
                        dest = list(key)
                        dest[site] += delta
                        _add(out, tuple(dest), coefficient * factor, dims)
    return out


def apply_hamiltonian(state, coupling=F(1), local_cutoff=None, center=True):
    """Apply the full polynomial H or its exact tensor compression P H P.

    The common scalar subtracted when center=True is 8+coupling/16.
    Inputs to a compressed Hamiltonian must already lie in its tensor box.
    """
    _validate_state(state)
    lam, dims = _coupling(coupling), cutoffs(local_cutoff)
    if type(center) is not bool:
        raise ValueError("center must be Boolean")
    if dims is not None and any(value and any(n >= d for n, d in zip(key, dims))
                                for key, value in state.items()):
        raise ValueError("compressed Hamiltonian input lies outside its tensor box")
    denominator = lcm(*(F(c).denominator for c in state.values())) if state else 1
    integers = {key: int(F(value) * denominator) for key, value in state.items() if value}
    return {key: F(value, 384 * denominator)
            for key, value in _numerator(integers, lam, dims, center).items()}


def _order(value, minimum=0):
    if type(value) is not int or not minimum <= value <= 4:
        raise ValueError(f"moment order must be an integer between {minimum} and four")
    return value


def vacuum_moments(coupling=F(1), max_order=4, local_cutoff=None):
    """Exact squared norms of (H-<H>vac)^r|vac>, including r=0."""
    lam, order, dims = _coupling(coupling), _order(max_order), cutoffs(local_cutoff)
    vector, denominator = {VACUUM: 1}, 1
    rows = []
    for power in range(order + 1):
        squared = F(sum(value*value*_weight(key) for key, value in vector.items()),
                    denominator*denominator)
        rows.append({"order": power, "norm_squared": squared,
                     "support_size": len(vector),
                     "maximum_total_occupation": max((sum(key) for key in vector), default=0),
                     "maximum_local_occupation": max((max(key) for key in vector), default=0)})
        if power < order:
            vector = _numerator(vector, lam, dims)
            denominator *= 384
    return rows


def coefficient_comparison(local_cutoff, order=2, coupling=F(1)):
    """Compare full and compressed Hamiltonian powers, before Taylor remainders."""
    dims, order, lam = cutoffs(local_cutoff), _order(order), _coupling(coupling)
    if dims is None:
        raise ValueError("a finite occupation cutoff is required")
    full, finite, denominator = {VACUUM: 1}, {VACUUM: 1}, 1
    rows = []
    for power in range(order + 1):
        difference = {key: full.get(key, 0)-finite.get(key, 0)
                      for key in full.keys() | finite.keys()}
        squared = lambda state: F(sum(c*c*_weight(k) for k, c in state.items()),
                                   denominator*denominator)
        rows.append({"order": power, "full_norm_squared": squared(full),
                     "finite_norm_squared": squared(finite),
                     "difference_norm_squared": squared(difference),
                     "coefficients_match": not any(difference.values())})
        if power < order:
            full = _numerator(full, lam)
            finite = _numerator(finite, lam, dims)
            denominator *= 384
    return rows


def truncation_certificate(time, local_cutoff, order=2, coupling=F(1), digits=18):
    """State-specific unitary/Taylor bound from the supplied product vacuum.

    Retain differences through degree order-1, then bound both order-th
    spectral remainders.  The same scalar shift is used on both sides.
    Record total variation is at most this vector bound for identical
    bounded readouts; arbitrary interventions require their own bounds.
    """
    time, order = _fraction(time, "model time"), _order(order, 1)
    rows = coefficient_comparison(local_cutoff, order, coupling)
    differences = sum((time**r/factorial(r)
                       * sqrt_upper(rows[r]["difference_norm_squared"], digits)
                       for r in range(1, order)), F(0))
    remainder = time**order/factorial(order) * (
        sqrt_upper(rows[order]["full_norm_squared"], digits)
        + sqrt_upper(rows[order]["finite_norm_squared"], digits))
    return {"local_cutoffs": list(cutoffs(local_cutoff)), "order": order,
            "model_time": time, "coefficient_rows": rows,
            "coefficient_difference_upper": differences,
            "spectral_remainder_upper": remainder,
            "vector_error_upper": min(F(2), differences + remainder),
            "readout_tv_upper": min(F(1), differences + remainder),
            "floating_or_service_error_included": False}


def interaction_witness_certificate(time, local_cutoff=None, digits=18):
    """Positive interacting-minus-free probability for e=(1-i z0^4)/5.

    The input is the Omega=2 product vacuum.  For the local rank-one projector
    onto e, tensored with identity elsewhere, p(0)=1/25 and
    p'_1(0)-p'_0(0)=1/200.  Each centered Hamiltonian
    contributes probability remainder <=t^2 (||H psi||^2+||H^2 psi||).
    A finite execution must retain local levels zero through four at site0.
    """
    time, dims = _fraction(time, "model time"), cutoffs(local_cutoff)
    if dims is not None and dims[0] < 5:
        raise ValueError("the witness requires at least five levels at the origin")
    moments = [vacuum_moments(F(lam), 2, dims) for lam in (0, 1)]
    coefficient = sum((row[1]["norm_squared"]
                       + sqrt_upper(row[2]["norm_squared"], digits)
                       for row in moments), F(0))
    leading, remainder = time/F(200), time*time*coefficient
    return {"local_cutoffs": None if dims is None else list(dims),
            "initial_probability": F(1, 25), "contrast_derivative": F(1, 200),
            "model_time": time, "probability_remainder_coefficient_upper": coefficient,
            "first_order_contrast": leading, "remainder_upper": remainder,
            "contrast_lower": leading-remainder,
            "strictly_positive": leading > remainder,
            "moment_rows": moments, "floating_or_service_error_included": False}
