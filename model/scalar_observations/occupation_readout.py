"""Exact unrestricted terminal effects on finite Hermite polynomials.

The final recorded probabilities can be evaluated without approximating the
disturbed terminal quantum state.  For the principal eta=1, theta=pi/2 read,
E0=(I-sin(Q))/2 and E1=(I+sin(Q))/2.  Every polynomial probability is an
affine rational expression in exp(-1/8).  This is a classical-output bound;
it does not replace the full-output residual certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt, prod

from .bounds import sqrt_upper


def exp_minus_eighth(terms=32):
    """Directed alternating-series enclosure, using only exact rationals."""
    if type(terms) is not int or terms < 2 or terms % 2:
        raise ValueError('a positive even Taylor order of at least two required')
    term = F(1)
    upper = term
    for n in range(1, terms+1):
        term *= -F(1, 8*n)
        upper += term
    lower = upper+term*(-F(1, 8*(terms+1)))
    return lower, upper


@lru_cache(maxsize=None)
def weyl_element(m, n):
    """Return rational real/imaginary coefficients of <m|exp(iQ)|n>.

    |n> here means (a*)**n|0>, with norm sqrt(n!), and Q=(a+a*)/2.
    The common Gaussian factor exp(-1/8) is not included.
    """
    if any(type(k) is not int or k < 0 for k in (m, n)):
        raise ValueError('nonnegative integer occupations required')
    real = imaginary = F(0)
    for lowered in range(max(0, n-m), n+1):
        raised = m-n+lowered
        power = lowered+raised
        coefficient = F(factorial(m)*factorial(n),
                        factorial(lowered)*factorial(n-lowered)*factorial(raised)*2**power)
        if power % 4 == 0:
            real += coefficient
        elif power % 4 == 1:
            imaginary += coefficient
        elif power % 4 == 2:
            real -= coefficient
        else:
            imaginary -= coefficient
    return real, imaginary


def _root_lower(value, digits=36):
    if value < 0:
        raise ValueError('nonnegative squared norm required')
    scale = 10**digits
    return F(isqrt(value.numerator*scale*scale//value.denominator), scale)


def _affine_interval(constant, coefficient, exponential):
    lo, hi = exponential
    return ((constant+coefficient*lo, constant+coefficient*hi) if coefficient >= 0
            else (constant+coefficient*hi, constant+coefficient*lo))


def terminal_probabilities(vector, weight, vector_error):
    """Enclose all four unrestricted records from a pre-final approximate state.

    The state has seven occupation coordinates and one retained first bit.
    Its physical amplitudes are sqrt(weight) times the raw coefficients.
    ``vector_error`` bounds this entire flagged state against its unrestricted
    target; it is not divided among record branches. No normalization occurs.
    """
    if type(weight) not in (int, F) or weight <= 0:
        raise ValueError('positive rational squared state weight required')
    if type(vector_error) not in (int, F) or vector_error < 0:
        raise ValueError('nonnegative rational vector error required')
    if type(vector.denominator) is not int or vector.denominator < 1:
        raise ValueError('positive integer coefficient denominator required')
    groups = {}
    norms = [F(0), F(0)]
    for key, value in vector.coefficients.items():
        if (not isinstance(key, tuple) or len(key) != 8 or key[7] not in (0, 1)
                or any(type(n) is not int or n < 0 for n in key[:7])
                or not isinstance(value, tuple) or len(value) != 2
                or any(type(n) is not int for n in value)):
            raise ValueError('seven occupations, one bit and two integer coefficients required')
        norms[key[7]] += F(weight*(value[0]**2+value[1]**2)*prod(factorial(n) for n in key[:7]),
                           vector.denominator**2)
        spectator = key[:1]+key[2:7]
        groups.setdefault((key[7], spectator), {})[key[1]] = value
    imaginary = [F(0), F(0)]
    for (bit, spectator), amplitudes in groups.items():
        factor = F(weight*prod(factorial(n) for n in spectator), vector.denominator**2)
        for m, (mr, mi) in amplitudes.items():
            for n, (nr, ni) in amplitudes.items():
                wr, wi = weyl_element(m, n)
                # Im(conj(c_m) <m|exp(iQ)|n> c_n), before the Gaussian.
                imaginary[bit] += factor*((mr*nr+mi*ni)*wi+(mr*ni-mi*nr)*wr)
    exponential = exp_minus_eighth()
    rows = []
    for first in (0, 1):
        for last in (0, 1):
            constant = norms[first]/2
            coefficient = (1 if last else -1)*imaginary[first]/2
            raw_lower, raw_upper = _affine_interval(constant, coefficient, exponential)
            # Positivity follows from the full unrestricted effect; intersect
            # the analytic enclosure with that property before taking roots.
            lower_root = _root_lower(max(F(0), raw_lower))
            upper_root = sqrt_upper(max(F(0), raw_upper))
            probability = [max(F(0), lower_root-vector_error)**2,
                           min(F(1), (upper_root+vector_error)**2)]
            rows.append({'bits': [first, last], 'affine_constant': constant,
                'exp_minus_eighth_coefficient': coefficient,
                'polynomial_effect_probability_interval': [max(F(0), raw_lower), max(F(0), raw_upper)],
                'unrestricted_probability_interval': probability})
    return {'branches': rows, 'exponential_interval': list(exponential),
            'raw_polynomial_mass': sum(norms, F(0)),
            'prefinal_vector_error_upper': F(vector_error),
            'raw_probability_half_l1_error_upper': vector_error+vector_error**2/2,
            'terminal_quantum_output_certified': False,
            'final_effect_is_unrestricted': True,
            'branches_renormalized': False}
