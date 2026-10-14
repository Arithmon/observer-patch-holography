"""Exact finite-Hermite residual certificates for seven-mode observation laws.

Every polynomial curve is embedded in unrestricted oscillator space. Duhamel's
identity bounds its error by the integral of its full Schrodinger residual;
finite and exterior residuals are orthogonal. All arithmetic determining an
error budget is rational, including endpoint rounding. These are ideal-model
certificates, not floating-point or physical hardware error certificates.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import factorial, lcm, prod

from .bounds import sqrt_upper


@dataclass
class Vector:
    """Bargmann coefficients; trailing tuple entries are retained record bits."""
    coefficients: dict[tuple[int, ...], tuple[int, int]]
    denominator: int = 1


def _add(out, key, real, imaginary):
    old = out.get(key, (0, 0))
    value = old[0] + real, old[1] + imaginary
    if value != (0, 0):
        out[key] = value
    else:
        out.pop(key, None)


@lru_cache(maxsize=None)
def _weight(key):
    return prod(factorial(n) for n in key[:7])


def inner_real(left, right):
    if len(left.coefficients) > len(right.coefficients):
        left, right = right, left
    numerator = 0
    for key, (real, imaginary) in left.coefficients.items():
        other = right.coefficients.get(key)
        if other is not None:
            numerator += (real*other[0] + imaginary*other[1])*_weight(key)
    return F(numerator, left.denominator*right.denominator)


def norm_squared(vector):
    return inner_real(vector, vector)


def combine(vectors):
    denominator = lcm(*(v.denominator for v in vectors))
    out = {}
    for vector in vectors:
        scale = denominator // vector.denominator
        for key, (real, imaginary) in vector.coefficients.items():
            _add(out, key, scale*real, scale*imaginary)
    return Vector(out, denominator)


def _round_integer(numerator, denominator):
    # Nearest integer; ties away from zero. The certificate counts the exact
    # resulting displacement, so the tie convention has no analytic role.
    sign = -1 if numerator < 0 else 1
    return sign*((2*abs(numerator) + denominator)//(2*denominator))


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _signed_time(value):
    if type(value) not in (int, F):
        raise ValueError("model duration must be an exact rational")
    return F(value)


def quantize(vector, bits):
    bits = _positive_integer(bits, "rounding bits")
    denominator = 1 << bits
    out, displacement = {}, {}
    for key, (real, imaginary) in vector.coefficients.items():
        rounded = (_round_integer(real*denominator, vector.denominator),
                   _round_integer(imaginary*denominator, vector.denominator))
        if rounded != (0, 0):
            out[key] = rounded
        _add(displacement, key,
             rounded[0]*vector.denominator-real*denominator,
             rounded[1]*vector.denominator-imaginary*denominator)
    error = sqrt_upper(norm_squared(Vector(displacement, denominator*vector.denominator)))
    return Vector(out, denominator), error


@lru_cache(maxsize=None)
def _bpower(n, power):
    out = {n: 1}
    for _ in range(power):
        following = {}
        for k, coefficient in out.items():
            following[k+1] = following.get(k+1, 0)+coefficient
            if k:
                following[k-1] = following.get(k-1, 0)+k*coefficient
        out = following
    return tuple(out.items())


class Generator:
    """Exact rational H or Q, with P applied only after the full polynomial.

    The Hamiltonian uses kappa=mass=c=1, Omega=2, and g=lambda/V in {0,1}.
    Its common scalar <blank|H|blank> is subtracted; this changes no record
    probabilities. Q denotes the normalized cos(x) smeared coordinate.
    """

    def __init__(self, dimensions, kind="H", coupling=1):
        self.dimensions = _dimensions(dimensions)
        if kind not in ("H", "Q"):
            raise ValueError("generator kind must be H or Q")
        if type(coupling) is not int or coupling not in (0, 1):
            raise ValueError("coupling must be exactly zero or one")
        self.kind, self.coupling = kind, coupling
        self.denominator = 768 if kind == "H" else 2
        self._cache = {}

    def transitions(self, occupations):
        cached = self._cache.get(occupations)
        if cached is not None:
            return cached
        out = {}

        def add(changes, coefficient):
            dest = list(occupations)
            for mode, n in changes:
                dest[mode] = n
            dest = tuple(dest)
            out[dest] = out.get(dest, 0)+coefficient

        if self.kind == "Q":
            n = occupations[1]
            add(((1, n+1),), 1)
            if n:
                add(((1, n-1),), n)
        else:
            g = self.coupling
            # 768*<blank|H|blank> = 3936 + 294*g.
            add((), 768*sum(occupations)+384*7-3936-294*g)
            for mode, n in enumerate(occupations):
                for dest, coefficient in _bpower(n, 2):
                    add(((mode, dest),), 96*(1 if mode == 0 else 2)*coefficient)
                add(((mode, n+2),), -384)
                if n >= 2:
                    add(((mode, n-2),), -384*n*(n-1))
                if g:
                    for dest, coefficient in _bpower(n, 4):
                        add(((mode, dest),), g*(2 if mode == 0 else 3)*coefficient)
            if g:
                for left, right in combinations(range(7), 2):
                    coefficient = 6 if left > 0 and (left-1)//2 == (right-1)//2 else 12
                    for ni, ci in _bpower(occupations[left], 2):
                        for nj, cj in _bpower(occupations[right], 2):
                            add(((left, ni), (right, nj)), coefficient*ci*cj)
        result = tuple((dest, coefficient,
                        all(n < d for n, d in zip(dest, self.dimensions)))
                       for dest, coefficient in out.items() if coefficient)
        self._cache[occupations] = result
        return result

    def apply(self, vector, exterior=False):
        out = {}
        for key, (real, imaginary) in vector.coefficients.items():
            for dest, coefficient, inside in self.transitions(key[:7]):
                if inside == exterior:
                    continue
                _add(out, dest+key[7:], coefficient*real, coefficient*imaginary)
        return Vector(out, vector.denominator*self.denominator)


def _time_term(generator, vector, time, power):
    applied = generator.apply(vector)
    return Vector({key: (time.numerator*imaginary, -time.numerator*real)
                   for key, (real, imaginary) in applied.coefficients.items()},
                  applied.denominator*time.denominator*power)


def polynomial_step(generator, vector, time, degree=8, bits=100):
    """Advance one curve and certify both unrestricted and compressed endpoints.

    With p(s)=sum_j term_j*(s/t)^j, the integrated squared exterior residual
    is t*sum_jk Re<L term_j,L term_k>/(j+k+1). The interior residual has just
    its degree-r term. Cauchy--Schwarz bounds integral ||residual|| by the
    square root of |t| times that exact integral. The curve need not converge
    to the unrestricted exponential as degree tends to infinity.
    """
    time = _signed_time(time)
    degree = _positive_integer(degree, "Taylor degree")
    bits = _positive_integer(bits, "rounding bits")
    terms = [vector]
    for power in range(1, degree+2):
        terms.append(_time_term(generator, terms[-1], time, power))
    exterior = [generator.apply(term, exterior=True) for term in terms[:-1]]
    exterior_squared = F(0)
    for j, left in enumerate(exterior):
        exterior_squared += inner_real(left, left)/F(2*j+1)
        for k in range(j):
            exterior_squared += 2*inner_real(left, exterior[k])/F(j+k+1)
    exterior_squared *= time*time
    interior_squared = F((degree+1)**2, 2*degree+1)*norm_squared(terms[-1])
    if exterior_squared < 0:
        raise ArithmeticError("exact squared residual integral is negative")
    endpoint, rounding = quantize(combine(terms[:-1]), bits)
    return endpoint, {
        "full_residual_error_upper": sqrt_upper(exterior_squared+interior_squared),
        "interior_residual_error_upper": sqrt_upper(interior_squared),
        "exterior_residual_error_upper": sqrt_upper(exterior_squared),
        "rounding_error_upper": rounding,
        "full_integrated_squared_residual_times_duration": exterior_squared+interior_squared,
    }


def flow(generator, vector, time, degree=8, substeps=8, bits=100):
    time = _signed_time(time)
    degree = _positive_integer(degree, "Taylor degree")
    substeps = _positive_integer(substeps, "flow substeps")
    bits = _positive_integer(bits, "rounding bits")
    out = vector
    totals = {key: F(0) for key in ("unrestricted_vector_error_upper",
              "finite_vector_error_upper", "exterior_residual_error_upper",
              "rounding_error_upper")}
    for _ in range(substeps):
        out, step = polynomial_step(generator, out, F(time)/substeps, degree, bits)
        totals["unrestricted_vector_error_upper"] += step["full_residual_error_upper"]+step["rounding_error_upper"]
        totals["finite_vector_error_upper"] += step["interior_residual_error_upper"]+step["rounding_error_upper"]
        totals["exterior_residual_error_upper"] += step["exterior_residual_error_upper"]
        totals["rounding_error_upper"] += step["rounding_error_upper"]
    return out, totals


def _read_combine(positive, negative):
    # (cos(Q/2)-sin(Q/2), cos(Q/2)+sin(Q/2)); the common 1/sqrt(2)
    # is represented by halving the squared state weight outside this function.
    denominator = 2*lcm(positive.denominator, negative.denominator)
    out = {}
    for vector, phase in ((positive, 1), (negative, -1)):
        scale = denominator//(2*vector.denominator)
        for key, (real, imaginary) in vector.coefficients.items():
            _add(out, key+(0,), scale*(real-phase*imaginary), scale*(imaginary+phase*real))
            _add(out, key+(1,), scale*(real+phase*imaginary), scale*(imaginary-phase*real))
    return Vector(out, denominator)


def _dimensions(value):
    try:
        dims = (value,)*7 if type(value) is int else tuple(value)
    except TypeError as exc:
        raise ValueError("seven integer occupation dimensions of at least two required") from exc
    if len(dims) != 7 or any(type(d) is not int or d < 2 for d in dims):
        raise ValueError("seven integer occupation dimensions of at least two required")
    return dims


def certificate(dimensions=4, coupling=1, sign=1, interval=F(1, 8), *,
                degree=8, substeps=8, bits=100, prefinal_only=False, return_terminal_state=False):
    """Certify the principal complete four-record history at exact gradient.

    Preparation/feedback kicks are 1/4, read strength is one, read angle pi/2.
    No measured probability is used to construct the interval or error budget.
    With prefinal_only=True, return the exact polynomial state and error ledger
    immediately before the final read, for independently certified effects.
    """
    dims = _dimensions(dimensions)
    if type(prefinal_only) is not bool or type(return_terminal_state) is not bool:
        raise ValueError("state return flags must be Boolean")
    if type(coupling) is not int or coupling not in (0, 1):
        raise ValueError("coupling must be exactly zero or one")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("native preparation sign must be minus or plus one")
    if type(interval) not in (int, F) or interval < 0:
        raise ValueError("nonnegative exact Hamiltonian model interval required")
    if any(type(n) is not int or n < 1 for n in (degree, substeps, bits)):
        raise ValueError("degree, substeps and rounding bits must be positive integers")
    hamiltonian, coordinate = Generator(dims, "H", coupling), Generator(dims, "Q", coupling)
    vector, weight = Vector({(0,)*7: (1, 0)}), F(1)
    errors = {"unrestricted_vector_error_upper": F(0), "finite_vector_error_upper": F(0)}
    stages = []

    def advance(name, generator, time, selected=None):
        nonlocal vector
        current = vector if selected is None else Vector(
            {key: value for key, value in vector.coefficients.items() if selected(key)}, vector.denominator)
        updated, budget = flow(generator, current, time, degree, substeps, bits)
        if selected is not None:
            kept = Vector({key: value for key, value in vector.coefficients.items() if not selected(key)},
                          vector.denominator)
            updated = combine([updated, kept])
        vector = updated
        for key in errors:
            errors[key] += sqrt_upper(weight)*budget[key]
        stages.append({"operation": name, "raw_state_squared_weight": weight,
                       **{key: sqrt_upper(weight)*value for key, value in budget.items()}})

    def read(name):
        nonlocal vector, weight
        positive, ep = flow(coordinate, vector, -F(1, 2), degree, substeps, bits)
        negative, em = flow(coordinate, vector, F(1, 2), degree, substeps, bits)
        # flow uses exp(-i G t); positive denotes exp(+i Q/2).
        vector, rounding = quantize(_read_combine(positive, negative), bits)
        weight /= 2
        budget = {key: sqrt_upper(weight*(ep[key]**2+em[key]**2))
                  +sqrt_upper(weight)*rounding for key in errors}
        for key in errors:
            errors[key] += budget[key]
        stages.append({"operation": name, "raw_state_squared_weight": weight,
                       **budget, "coherent_weyl_branches_retained": True,
                       "recombination_rounding_error_upper": sqrt_upper(weight)*rounding})

    advance("initial_kick", coordinate, -F(sign, 4))
    advance("first_evolution", hamiltonian, interval)
    read("first_read")
    advance("conditional_feedback_kick", coordinate, -F(1, 4), lambda key: key[7] == 1)
    advance("second_evolution", hamiltonian, interval)
    if prefinal_only:
        return {"state": vector, "raw_state_squared_weight": weight,
                "stages": stages, **errors, "local_dimensions": list(dims),
                "coupling": coupling, "preparation_sign": sign}
    from .occupation_readout import terminal_probabilities
    terminal_effect = terminal_probabilities(vector, weight, errors["unrestricted_vector_error_upper"])
    read("second_read")
    rows = []
    for first in (0, 1):
        for last in (0, 1):
            branch = Vector({key: value for key, value in vector.coefficients.items()
                             if key[7:] == (first, last)}, vector.denominator)
            probability = weight*norm_squared(branch)
            bounds = {}
            for label, error in (("unrestricted", errors["unrestricted_vector_error_upper"]),
                                 ("finite", errors["finite_vector_error_upper"])):
                radius = 2*sqrt_upper(probability)*error+error*error
                bounds[label+"_probability_interval"] = [max(F(0), probability-radius), min(F(1), probability+radius)]
            rows.append({"bits": [first, last], "raw_polynomial_probability": probability, **bounds})
    eu, ed = errors["unrestricted_vector_error_upper"], errors["finite_vector_error_upper"]
    result = {"schema": "oph.scalar_observations.occupation_residual.v1",
            "local_dimensions": list(dims), "coupling": coupling, "preparation_sign": sign,
            "interval": F(interval), "degree": degree, "substeps_per_flow": substeps,
            "rounding_bits": bits, "raw_state_squared_weight": weight,
            "raw_polynomial_mass": sum((row["raw_polynomial_probability"] for row in rows), F(0)),
            "branches": rows, "stages": stages, **errors,
            "terminal_effect_certificate": terminal_effect,
            "finite_to_unrestricted_record_tv_upper": min(F(1), eu+ed),
            "finite_to_unrestricted_retained_trace_distance_upper": min(F(1), eu+ed),
            "raw_probability_half_l1_error_upper": eu+eu*eu/2,
            "finite_sector_dimension": dims[1]*prod((d+1)//2 for i, d in enumerate(dims) if i != 1),
            "occupied_polynomial_coefficients": len(vector.coefficients),
            "scope": {"gradient_stiffness": "1", "same_supplied_seven_mode_action": True,
                      "every_record_retained": True, "branch_normalization_or_postselection": False,
                      "rational_residual_and_rounding_arithmetic": True,
                      "source_splitting_error_included": False,
                      "floating_or_physical_service_error_included": False,
                      "momentum_band_removal_or_empirical_agreement_claimed": False}}
    if return_terminal_state:
        result["terminal_state"] = vector
    return result
