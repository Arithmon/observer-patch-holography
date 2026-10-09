"""Positive secular reduction of the declared one-plaquette S3 Hamiltonian.

K=1, four electric links, and character ordering (triv, sign, std) are fixed.
The rational root bracket is exact on the supplied coupling. Transcendental
readouts are high-precision numerical evaluations, not interval certificates.
See README.md for the reduction, reporting range and heat-kernel comparison.
"""

from decimal import Decimal
from fractions import Fraction
import math
from numbers import Integral, Real

import mpmath
import numpy as np


def _coupling(value):
    """Preserve original real scalar values before any binary64 conversion."""
    if isinstance(value, (bool, np.bool_)) or np.ma.isMaskedArray(value):
        raise ValueError("coupling must be a finite nonnegative real scalar")
    if isinstance(value, Integral):
        result = Fraction(int(value))
    elif isinstance(value, Fraction):
        result = Fraction(int(value.numerator), int(value.denominator))
    elif isinstance(value, (Real, Decimal)):
        try:
            numerator, denominator = value.as_integer_ratio()
            result = Fraction(int(numerator), int(denominator))
        except (ValueError, OverflowError, AttributeError) as exc:
            raise ValueError("coupling must be a finite nonnegative real scalar") from exc
    else:
        raise ValueError("coupling must be a finite nonnegative real scalar")
    if result < 0:
        raise ValueError("coupling must be nonnegative")
    return result


def ground_root_bracket(h, *, bits=256):
    """Enclose x=-E0>0 by rational bisection of the monotone secular equation.

    F(x)=x+12h-1-1/x-1/(24h+x) has positive derivative. The
    initial bracket has upper/lower ratio at most six for every h>=0;
    resolving a tiny root therefore costs no additional exponent steps.
    """
    h = _coupling(h)
    if type(bits) is not int or not 64 <= bits <= 4096:
        raise ValueError("root precision must be an integer from 64 to 4096 bits")
    if not h:
        return Fraction(2), Fraction(2)
    a, b = 24*h, 12*h-1
    lower = Fraction(1) if b <= 0 else 1/(b+2)
    upper = Fraction(2) if b <= 0 else min(Fraction(2), 2/b)
    for _ in range(bits):
        x = (lower+upper)/2
        # Multiplication by x*(a+x)>0 preserves the sign of F(x).
        value = x*(a+x)*(b+x)-(a+2*x)
        if value > 0:
            upper = x
        elif value < 0:
            lower = x
        else:
            return x, x
    return lower, upper


def _mp_fraction(ctx, value):
    return ctx.mpf(value.numerator)/value.denominator


def _readout(value, name, *, zero_allowed=False):
    result = float(value)
    if (not math.isfinite(result) or (result == 0 and not zero_allowed)
            or (value and abs((result-value)/value) > value.context.mpf("1e-12"))):
        raise ValueError(f"{name} exceeds the resolved binary64 reporting range")
    return result


def _state(h):
    h = _coupling(h)
    bounds = ground_root_bracket(h)
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    low, high = (_mp_fraction(ctx, bound) for bound in bounds)
    x = (low+high)/2
    a = 24*_mp_fraction(ctx, h)
    ratio = x/(a+x)
    norm = 1+ratio**2+x*x
    probabilities = {name: _readout(value/norm, f"{name} probability")
                     for name, value in (("triv", ctx.mpf(1)),
                                         ("sign", ratio**2), ("std", x*x))}
    return ctx, x, a, norm, probabilities, bounds


def s3_diagnostics(h):
    """Full positive sector law and cancellation-free held-out comparison.

    The fitted standard-sector time may be negative; `diffusion_fit` makes
    that distinction explicit. No finite coupling in this family satisfies
    the exact d_R heat-kernel law. Near t=0 the normalized readouts are only
    returned if the root bracket resolves them; the log discrepancy has no
    such division. An unresolved or unrepresentable readout raises.
    """
    ctx, x, a, norm, probabilities, bounds = _state(h)
    low, high = (_mp_fraction(ctx, bound) for bound in bounds)
    log_std_ratio = 2*ctx.log(x)-ctx.log(2)
    fit_time = -log_std_ratio/3
    # The bracket, not the rounded fit-time sign, decides if normalization
    # is resolved. Fitted time decreases with x.
    fit_low = (ctx.log(2)-2*ctx.log(high))/3
    fit_high = (ctx.log(2)-2*ctx.log(low))/3
    if fit_low*fit_high <= 0 or abs(fit_high-fit_low) > abs(fit_time)*ctx.mpf("1e-20"):
        raise ValueError("diffusion fit is unresolved at the root precision")
    z = x+a/2
    # y=x*(x+24h)=2z/(z-1). Avoid subtracting almost equal logarithms,
    # or rounding the sign/predicted ratio to one before taking its log.
    log_discrepancy = -2*ctx.log1p(1/(z-1))
    predicted = x**4/(4*norm)
    normalized = log_discrepancy/abs(2*log_std_ratio)
    ratio_excess = log_discrepancy/log_std_ratio
    return {
        "probabilities": probabilities,
        "fit_time": _readout(fit_time, "fitted time"),
        "diffusion_fit": fit_time > 0,
        "predicted_sign": _readout(predicted, "predicted sign probability"),
        "log_discrepancy": _readout(log_discrepancy, "log discrepancy"),
        "normalized_log_residual": _readout(normalized, "normalized log residual"),
        "log_ratio_excess": _readout(ratio_excess, "log-ratio excess"),
        "root_bracket": bounds,
    }


def s3_edge_distribution(h):
    """The three probabilities; every positive component must remain resolved."""
    # Keep the probability observable independent of ill-conditioned fit
    # normalization at the zero-diffusion boundary.
    return _state(h)[4]
