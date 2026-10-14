"""Exact polynomial certificates for a continued scalar-field record history.

The integer-pair/common-denominator representation and spectral Taylor
remainder generalize ``model.common_source.infinite_history``.  The full
eight-oscillator Hamiltonian oracle is in ``exact``.  Every branch remains
unnormalized.  No computed target probability or floating residual enters
the certificate.  Time is the supplied Hamiltonian model coordinate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import lcm

from . import exact


@dataclass
class Polynomial:
    coefficients: dict[tuple[int, ...], tuple[int, int]]
    denominator: int = 1


def _add(out, key, real, imaginary):
    old = out.get(key, (0, 0))
    value = old[0]+real, old[1]+imaginary
    if value != (0, 0):
        out[key] = value
    else:
        out.pop(key, None)


def vacuum():
    return Polynomial({exact.VACUUM: (1, 0)})


def norm_squared(vector):
    return F(sum((real*real+imaginary*imaginary)*exact._weight(key)
                 for key, (real, imaginary) in vector.coefficients.items()),
             vector.denominator**2)


def difference(left, right):
    denominator = lcm(left.denominator, right.denominator)
    out = {}
    for vector, sign in ((left, 1), (right, -1)):
        factor = sign*denominator//vector.denominator
        for key, (real, imaginary) in vector.coefficients.items():
            _add(out, key, factor*real, factor*imaginary)
    return Polynomial(out, denominator)


def _hamiltonian_numerator(coefficients, coupling, dims):
    real = {key: value[0] for key, value in coefficients.items() if value[0]}
    imaginary = {key: value[1] for key, value in coefficients.items() if value[1]}
    real = exact._numerator(real, coupling, dims)
    imaginary = exact._numerator(imaginary, coupling, dims)
    return {key: (real.get(key, 0), imaginary.get(key, 0))
            for key in real.keys() | imaginary.keys()}


def advance(vector, time, degree, coupling, dims=None):
    """Taylor polynomial and the norm of its next full Hamiltonian term.

    The spectral theorem gives a remainder at most
    |t|^(r+1)||H^(r+1)psi||/(r+1)! for any finite-support input.
    Full support is retained throughout every Hamiltonian power.
    """
    time = exact._fraction(time, "model time")
    if type(degree) is not int or not 1 <= degree <= 4:
        raise ValueError("Taylor degree must be an integer between one and four")
    coupling, dims = exact._coupling(coupling), exact.cutoffs(dims)
    term, total = dict(vector.coefficients), dict(vector.coefficients)
    denominator = vector.denominator
    maximum_support = len(term)
    for power in range(1, degree+2):
        term = {key: (time.numerator*imag, -time.numerator*real)
                for key, (real, imag) in _hamiltonian_numerator(term, coupling, dims).items()}
        factor = 384*time.denominator*power
        denominator *= factor
        maximum_support = max(maximum_support, len(term))
        if power == degree+1:
            remainder = Polynomial(term, denominator)
            return (Polynomial(total, denominator//factor),
                    exact.sqrt_upper(norm_squared(remainder)), maximum_support)
        total = {key: (real*factor, imag*factor) for key, (real, imag) in total.items()}
        for key, value in term.items():
            _add(total, key, *value)
    raise AssertionError("unreachable")


def project(vector, bit, instrument_sign=-1):
    """Local P=|e><e|, e=(1+sign*i z0^4)/5, and its complete complement."""
    if type(bit) is not int or bit not in (0, 1):
        raise ValueError("pointer bit must be zero or one")
    if type(instrument_sign) is not int or instrument_sign not in (-1, 1):
        raise ValueError("instrument sign must be minus or plus one")
    positive = {}
    for key, (real, imaginary) in vector.coefficients.items():
        if key[0] == 0:
            inner = real, imaginary
        elif key[0] == 4:
            inner = instrument_sign*24*imaginary, -instrument_sign*24*real
        else:
            continue
        zero, four = (0,)+key[1:], (4,)+key[1:]
        _add(positive, zero, *inner)
        _add(positive, four, -instrument_sign*inner[1], instrument_sign*inner[0])
    projected = Polynomial(positive, 25*vector.denominator)
    return projected if bit == 1 else difference(vector, projected)


def apply_feedback(vector, first_bit, enabled=True, phase=1):
    """The committed one-bit record phases the origin's fourth level by ±i."""
    if type(first_bit) is not int or first_bit not in (0, 1):
        raise ValueError("first pointer bit must be zero or one")
    if type(enabled) is not bool or type(phase) is not int or phase not in (-1, 1):
        raise ValueError("feedback requires a Boolean switch and phase sign ±1")
    if not enabled or first_bit == 0:
        return vector
    return Polynomial({key: (-phase*imag, phase*real) if key[0] == 4 else (real, imag)
                       for key, (real, imag) in vector.coefficients.items()}, vector.denominator)


def _run(time, degree, local_cutoff, coupling, feedback, feedback_phase, instrument_sign):
    time, dims = exact._fraction(time, "model time"), exact.cutoffs(local_cutoff)
    coupling = exact._coupling(coupling)
    if dims is not None and dims[0] < 5:
        raise ValueError("the origin must retain at least five occupation levels")
    if type(feedback) is not bool or type(feedback_phase) is not int or feedback_phase not in (-1, 1):
        raise ValueError("feedback requires a Boolean switch and phase sign ±1")
    first, first_error, maximum_support = advance(vacuum(), time, degree, coupling, dims)
    branches, vectors, second_errors = [], {}, []
    for bit in (0, 1):
        selected = apply_feedback(project(first, bit, instrument_sign), bit, feedback, feedback_phase)
        second, second_error, support = advance(selected, time, degree, coupling, dims)
        maximum_support = max(maximum_support, support)
        second_errors.append(second_error)
        branch_error = first_error+second_error
        for final_bit in (0, 1):
            vector = project(second, final_bit, instrument_sign)
            probability = norm_squared(vector)
            event_error = 2*exact.sqrt_upper(probability)*branch_error+branch_error**2
            branches.append({"records": [bit, final_bit],
                             "raw_polynomial_probability": probability,
                             "probability_interval": [max(F(0), probability-event_error),
                                                      min(F(1), probability+event_error)],
                             "event_probability_error_upper": event_error,
                             "polynomial_support_size": len(vector.coefficients)})
            vectors[bit, final_bit] = vector
    error = first_error+exact.sqrt_upper(sum((x*x for x in second_errors), F(0)))
    mass = sum((row["raw_polynomial_probability"] for row in branches), F(0))
    result = {"schema": "oph.relativistic_scalar.exact_history.v1",
              "model_time_per_segment": time, "taylor_degree": degree,
              "coupling": F(coupling), "feedback_enabled": feedback,
              "feedback_phase_sign": feedback_phase,
              "instrument_sign": instrument_sign,
              "local_cutoffs": None if dims is None else list(dims),
              "branches": branches, "raw_polynomial_mass": mass,
              "normalization_interval": [sum((row["probability_interval"][0] for row in branches), F(0)),
                                         sum((row["probability_interval"][1] for row in branches), F(0))],
              "first_segment_vector_error_upper": first_error,
              "second_segment_branch_errors_upper": second_errors,
              "complete_retained_vector_error_upper": error,
              "raw_weight_classical_half_l1_error_upper": error+error*error/2,
              "maximum_expanded_polynomial_support": maximum_support,
              "scope": {"initial_state": "Omega=2 product vacuum, not the interacting vacuum",
                        "time_parameter": "supplied Hamiltonian model coordinate",
                        "spatial_sites": 8,
                        "oscillator_space": "infinite" if dims is None else "exact P H P compression",
                        "all_four_record_branches_retained": True,
                        "branch_normalization_or_postselection": False,
                        "floating_or_service_error_included": False,
                        "infinite_spatial_volume_claimed": False,
                        "relativistic_continuum_limit_claimed": False}}
    return result, vectors


def certificate(time=F(1, 128), degree=2, local_cutoff=None, coupling=F(1),
                feedback=True, feedback_phase=1, instrument_sign=-1):
    """All four exact-law probability intervals, derived from rational polynomials.

    First-stage error propagates through the retained instrument as one
    vector. The second-stage remainder vectors have orthogonal first-record
    labels, hence combine by their Euclidean norm. This gives
    epsilon=epsilon_first+sqrt(sum_b epsilon_second[b]^2). The polynomial
    output need not have unit norm: its half-L1 probability error is at
    most epsilon+epsilon^2/2, rather than a normalization of its branches.
    """
    return _run(time, degree, local_cutoff, coupling, feedback, feedback_phase, instrument_sign)[0]


def comparison_certificate(local_cutoff, time=F(1, 128), degree=2, coupling=F(1),
                           feedback=True, feedback_phase=1, instrument_sign=-1):
    """Complete retained-output comparison: finite tensor box versus infinite oscillators."""
    dims = exact.cutoffs(local_cutoff)
    if dims is None:
        raise ValueError("a finite occupation cutoff is required")
    infinite, left = _run(time, degree, None, coupling, feedback, feedback_phase, instrument_sign)
    finite, right = _run(time, degree, dims, coupling, feedback, feedback_phase, instrument_sign)
    squared = sum((norm_squared(difference(left[key], right[key])) for key in left), F(0))
    polynomial_difference = exact.sqrt_upper(squared)
    error = (infinite["complete_retained_vector_error_upper"] + polynomial_difference
             + finite["complete_retained_vector_error_upper"])
    return {"schema": "oph.relativistic_scalar.exact_cutoff_comparison.v1",
            "infinite": infinite, "finite": finite,
            "polynomial_vector_difference_norm_squared": squared,
            "polynomial_vector_difference_upper": polynomial_difference,
            "complete_retained_vector_error_upper": min(F(2), error),
            "complete_history_tv_upper": min(F(1), error),
            "complete_retained_trace_distance_upper": min(F(1), error),
            "floating_or_service_error_included": False}
