"""Decimal decision arithmetic for the existing FZ-10 rule, not a new freeze."""
from decimal import Decimal, localcontext
import re


def number(value):
    if type(value) is not str or len(value) > 90 or re.fullmatch(r'-?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]{1,3})?', value) is None:
        raise ValueError('bounded decimal string required')
    result = Decimal(value)
    # At most 13 integer and 80 fractional places: the 110-digit context then
    # computes every subtraction and threshold product exactly. Division is
    # display-only and does not decide the verdict.
    if not result.is_finite() or abs(result) > Decimal('1e12') or result.as_tuple().exponent < -80:
        raise ValueError('finite bounded observable required')
    return result


def verdict(value, sigma, center=Decimal('1776.969027')):
    distance = abs(value-center)
    # The numbered FAIL bullet is unconditional; the precision gate belongs
    # to COMPATIBLE. "Every other outcome" is the third, residual category.
    if distance > 3*sigma:
        return 'FAIL'
    if distance <= 2*sigma and sigma <= Decimal('.045'):
        return 'COMPATIBLE'
    return 'INCONCLUSIVE'


def compare(payload):
    if type(payload) is not dict or payload.keys() != {'value', 'sigma', 'unit'} or payload['unit'] != 'MeV':
        raise ValueError('exact measured-value schema and MeV required')
    with localcontext() as c:
        c.prec = 110
        value, sigma = number(payload['value']), number(payload['sigma'])
        if value <= 0 or not Decimal('1e-30') <= sigma:
            raise ValueError('positive physical mass and bounded nonzero standard uncertainty required')
        center = Decimal('1776.969027')
        low, high = Decimal('1776.968991'), Decimal('1776.969063')
        d = abs(value-center)
        dmin = max(low-value, value-high, Decimal(0))
        dmax = max(abs(value-low), abs(value-high))
        interval_verdict = ('FAIL' if dmin > 3*sigma else
                            'COMPATIBLE' if dmax <= 2*sigma and sigma <= Decimal('.045') else
                            'INCONCLUSIVE')
        return dict(frozen_center_verdict=verdict(value, sigma),
                    distance_mev=str(d), distance_over_reported_sigma=str(d/sigma),
                    window_robust_verdict=interval_verdict,
                    interpretation='conditional balanced-mass relation; no OPH-versus-Koide discrimination')


def reference_log_likelihood_ratios(value, sigma, prediction):
    # Normal-error illustration only. No null tail area, fitted selector
    # penalty or Bayes factor against the free-mass model is asserted.
    with localcontext() as c:
        c.prec = 110
        y, s, p = map(number, (value, sigma, prediction))
        if y <= 0 or p <= 0 or not Decimal('1e-30') <= s:
            raise ValueError('positive masses and bounded nonzero standard uncertainty')
        return dict(oph_to_koide='0', twice_log_free_mass_to_balanced=str(((y-p)/s)**2))
