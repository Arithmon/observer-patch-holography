"""Finite-matrix error budgets for record-based interaction identification.

The bounds concern the supplied compressed canonical action and Hamiltonian
model coordinate.  They do not bound occupation truncation, clock calibration,
floating-point execution or physical source service errors.  Exact rational
arithmetic and outward square roots precede the exported upper floats.
"""
from __future__ import annotations

from collections.abc import Mapping
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
import math
from numbers import Integral, Real

from ..relativistic_scalar.exact import sqrt_upper


def _rational(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (Real, F)):
        raise ValueError(f"{name} must be a finite nonnegative real")
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"{name} must be a finite nonnegative real, positive where required")
    # Fractions retain exact supplied rationals; a float means its exact stored
    # binary value.  This is not an uncertainty certificate for that input.
    return F(value)


def _upper(value):
    """Smallest nearby float no smaller than a nonnegative exact rational."""
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("the upper bound is too large to represent")
    return math.nextafter(result, math.inf) if F(result) < value else result


def coefficient_radius(time, norm_cap, probability_error=0):
    """Bound a Richardson estimate of one real off-diagonal matrix entry.

    For real symmetric H, real initial |j>, and e_±=(|j>±i|k>)/sqrt(2),
    let s(t)=(p_-(t)-p_+(t))/2.  Then s is odd and s'(0)=H[k,j].
    For B >= ||H-c I||, the estimate [8s(t/2)-s(t)]/(3t) has error at
    most B^5 t^4/18.  A simultaneous absolute error delta in each of the
    four input probabilities adds 3 delta/t.  The latter must be supplied;
    exact ideal record laws have delta zero, not a finite-shot guarantee.
    """
    t = _rational(time, "model duration", positive=True)
    b = _rational(norm_cap, "centered generator norm cap")
    delta = _rational(probability_error, "per-probability error")
    return _upper(b**5*t**4/18 + 3*delta/t)


def richardson_radius(first_time, second_time, norm_cap, probability_error=0):
    """The same coefficient certificate for arbitrary increasing model times.

    For 0<u<v, the estimator is
    [v² s(u)/u-u² s(v)/v]/(v²-u²).  This cancels the cubic term without
    replacing independently decoded clock readings by an exact ratio two.
    Clock-reading uncertainty itself is not included in this certificate.
    """
    u = _rational(first_time, "first model duration", positive=True)
    v = _rational(second_time, "second model duration", positive=True)
    b = _rational(norm_cap, "centered generator norm cap")
    delta = _rational(probability_error, "per-probability error")
    if u >= v:
        raise ValueError("the second model duration must exceed the first")
    gap = v*v-u*u
    finite = (2*b)**5*u*u*v*v*(u*u+v*v)/(240*gap)
    probability = delta*(v*v/u+u*u/v)/gap
    return {"finite_time_radius": _upper(finite),
            "finite_time_radius_exact": str(finite),
            "probability_radius": _upper(probability),
            "probability_radius_exact": str(probability),
            "radius": _upper(finite+probability),
            "radius_exact": str(finite+probability)}


def simultaneous_shot_count(probability_error, laws, failure_probability=F(1, 20)):
    """Conservative fresh Bernoulli trials per law for a simultaneous margin.

    The union bound and Hoeffding inequality give failure probability at most
    2 M exp(-2 N delta^2).  Trials are independent fresh preparations within
    each law.  Correlated systematic error is not covered.  The returned
    resource count is a proposed sampling budget, not an executed experiment.
    """
    delta = _rational(probability_error, "probability margin", positive=True)
    alpha = _rational(failure_probability, "failure probability", positive=True)
    if delta >= 1 or alpha >= 1:
        raise ValueError("probability margin and failure probability must be below one")
    if isinstance(laws, bool) or not isinstance(laws, Integral) or laws < 1:
        raise ValueError("the number of probability laws must be a positive integer")
    with localcontext() as context:
        context.prec = 80
        context.rounding = ROUND_CEILING
        ratio = Decimal(2*int(laws)*alpha.denominator)/Decimal(alpha.numerator)
        # Decimal.ln is correctly rounded to nearest; one outward ulp makes
        # its upper bound explicit, including for extremely small margins.
        logarithm_upper = ratio.ln().next_plus()
        context.rounding = ROUND_FLOOR
        variance = Decimal(2*delta.numerator**2)/Decimal(delta.denominator**2)
        context.rounding = ROUND_CEILING
        trials = int((logarithm_upper/variance).to_integral_value(rounding=ROUND_CEILING))
    return {
        "probability_laws": int(laws),
        "trials_per_law": trials,
        "total_fresh_trials": int(laws)*trials,
        "simultaneous_probability_error": _upper(delta),
        "failure_probability_upper": _upper(alpha),
        "sampling_executed": False,
        "systematic_error_included": False,
    }


def oscillator_norm_bounds(dimension, omega=2):
    """Bounds for true compressed infinite-oscillator q, q² and q⁴.

    q² <= 2(N+1/2)/omega as quadratic forms.  Applying that estimate before
    and after one q gives ||P q^4 P|| <= (2d-1)(2d+1)/omega².  For PqP,
    the triangle bound on the truncated lowering and raising operators is
    ||PqP|| <= sqrt(2(d-1)/omega).  No power of truncated q is substituted
    for the compressed infinite polynomial.
    """
    if (isinstance(dimension, bool) or not isinstance(dimension, Integral)
            or dimension < 1):
        raise ValueError("occupation dimension must be a positive integer")
    d = int(dimension)
    omega = _rational(omega, "basis frequency", positive=True)
    return {"q": sqrt_upper(F(2*(d-1))/omega),
            "q2": F(2*d-1)/omega,
            "q4": F((2*d-1)*(2*d+1))/omega**2}


def parameter_uncertainty_bound(dimensions, onsite_quadratic_radii,
                                quartic_radii, bond_radii, duration_total,
                                omega=2):
    """Propagate coefficient intervals through a complete retained history.

    H=sum_i[p_i²/2+A_i q_i²/2+g_i q_i⁴]-sum_ij kappa_ij q_i q_j.
    Supply absolute radii for A_i, g_i and kappa_ij around the predicted
    generator.  ``bond_radii`` is a mapping (i,j)->radius, or a sequence of
    (i,j,radius) triples; all uncertain pairs, including inferred nonedges,
    must be included.  No cancellation or statistical independence is used.

    With identical initial state, instruments and feedback, Duhamel's
    inequality and retained-instrument contraction give complete-history
    TV <= sum(segment durations)*||delta H||, capped by one.  This also
    bounds the retained quantum-record trace distance.  It applies to the
    exact finite compressed evolutions; splitting error is a separate term.
    """
    dims = tuple(dimensions)
    if not dims:
        raise ValueError("at least one occupation dimension is required")
    norms = [oscillator_norm_bounds(d, omega) for d in dims]
    quadratic, quartic = tuple(onsite_quadratic_radii), tuple(quartic_radii)
    if len(quadratic) != len(dims) or len(quartic) != len(dims):
        raise ValueError("each site needs a quadratic and a quartic radius")
    quadratic = [_rational(x, "onsite quadratic radius") for x in quadratic]
    quartic = [_rational(x, "quartic radius") for x in quartic]
    total = sum((a*op["q2"]/2+g*op["q4"]
                 for a, g, op in zip(quadratic, quartic, norms)), F(0))
    pairs = ((i, j, radius) for (i, j), radius in bond_radii.items()) if isinstance(
        bond_radii, Mapping) else iter(bond_radii)
    visited = set()
    for i, j, radius in pairs:
        if (isinstance(i, bool) or isinstance(j, bool)
                or not isinstance(i, Integral) or not isinstance(j, Integral)
                or not 0 <= i < j < len(dims) or (i, j) in visited):
            raise ValueError("bond uncertainty requires unique pairs 0 <= i < j < sites")
        visited.add((i, j))
        total += _rational(radius, "bond radius")*norms[i]["q"]*norms[j]["q"]
    duration = _rational(duration_total, "total Hamiltonian model duration")
    tv = min(F(1), duration*total)
    return {
        "hamiltonian_norm_error_upper": _upper(total),
        "complete_history_tv_upper": _upper(tv),
        "complete_retained_trace_distance_upper": _upper(tv),
        "total_model_duration": _upper(duration),
        "exact_hamiltonian_norm_bound": str(total),
        "exact_complete_history_tv_bound": str(tv),
        "scope": {
            "finite_occupation_compression": True,
            "identical_initial_state_instruments_and_feedback": True,
            "all_record_branches_retained": True,
            "occupation_error_included": False,
            "splitting_error_included": False,
            "floating_or_physical_service_error_included": False,
        },
    }
