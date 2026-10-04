"""Independent full-source ray quadrature; a control of the separate exact bounds."""
from functools import lru_cache
from fractions import Fraction as F
import mpmath
from .format import need, rational


@lru_cache(maxsize=4, typed=True)
def bending_difference(epsilon, dps=40):
    """Integrate the actual common lapse and both radial metrics on the bump.

    Returns a decimal numerical control, not a certified interval. The exact
    analytic bounds are checked separately in check.verify_ray_enclosure.
    """
    need(type(epsilon) is F and 0 < epsilon <= F(1, 100), 'ray strength domain')
    need(type(dps) is int and 30 <= dps <= 60, 'bounded ray precision')
    ctx = mpmath.mp.clone()
    ctx.dps = dps
    e = ctx.mpf(epsilon.numerator)/epsilon.denominator
    r0 = ctx.mpf(1)/2
    def integrand(r):
        # nu(r)-nu(r0), reconstructed from the baseline Einstein-cluster mass.
        phase = ctx.quad(lambda x:e*x/(1+x**3-2*e*x*x), [r0, 1, r])
        b_root_A = r0*ctx.exp(phase)
        m0 = e*r**3/(1+r**3)
        d = e*(r-1)**4*(2-r)**4/10
        f_plus = 1-2*(m0+d)/r
        f_minus = 1-2*(m0-d)/r
        a, b = ctx.sqrt(f_plus), ctx.sqrt(f_minus)
        # Rationalize 1/a-1/b: avoids subtracting nearly equal weak-field roots.
        root_difference = (4*d/r)/(a*b*(a+b))
        return 2*b_root_A*root_difference/(r*r*ctx.sqrt(1-(b_root_A/r)**2))
    return ctx.nstr(ctx.quad(integrand, [1, ctx.mpf(3)/2, 2]), dps-8)


def verify_enclosures(rows):
    ctx = mpmath.mp.clone()
    ctx.dps = 50
    for row in rows:
        e = rational(row['epsilon'])
        lo, hi = rational(row['lower']), rational(row['upper'])
        actual = ctx.mpf(bending_difference(e))
        need(ctx.mpf(lo.numerator)/lo.denominator < actual < ctx.mpf(hi.numerator)/hi.denominator,
             'actual full-source null integral lies inside the exact enclosure')
