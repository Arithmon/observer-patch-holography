"""Outward rational intervals on a fixed decimal lattice; no floating arithmetic."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import factorial, isqrt

SCALE = 10**50
EXPORT_SCALE = 10**12


def ceil_div(a, b):
    return -((-a)//b)


@dataclass(frozen=True)
class I:
    lo: int
    hi: int

    def __post_init__(self):
        if type(self.lo) is not int or type(self.hi) is not int or self.lo > self.hi:
            raise ValueError("invalid integer interval")

    @staticmethod
    def of(x):
        if isinstance(x, I):
            return x
        if type(x) not in (int, Fraction):
            raise ValueError("only exact integers and Fractions are accepted")
        x = Fraction(x)*SCALE
        return I(x.numerator//x.denominator, ceil_div(x.numerator, x.denominator))

    def __add__(self, other):
        other = I.of(other)
        return I(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-I.of(other)

    def __rsub__(self, other):
        return I.of(other)+-self

    def __mul__(self, other):
        other = I.of(other)
        v = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return I(min(v)//SCALE, ceil_div(max(v), SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I.of(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError("division interval includes zero")
        v = [Fraction(a*SCALE, b) for a in (self.lo, self.hi)
             for b in (other.lo, other.hi)]
        a, b = min(v), max(v)
        return I(a.numerator//a.denominator, ceil_div(b.numerator, b.denominator))

    def __rtruediv__(self, other):
        return I.of(other)/self

    def __pow__(self, n):
        if type(n) is not int or n < 0:
            raise ValueError("nonnegative integer powers required")
        out = I.of(1)
        for _ in range(n):
            out = out*self
        return out

    def sqrt(self):
        if self.lo < 0:
            raise ValueError("negative square root")
        lo, hi = isqrt(self.lo*SCALE), isqrt(self.hi*SCALE)
        return I(lo, hi+(hi*hi != self.hi*SCALE))

    def abs(self):
        return I(0 if self.lo <= 0 <= self.hi else min(abs(self.lo), abs(self.hi)),
                 max(abs(self.lo), abs(self.hi)))

    def upper(self):
        return I(self.hi, self.hi)

    def encode(self):
        unit = SCALE//EXPORT_SCALE
        return [str(Fraction(self.lo//unit, EXPORT_SCALE)),
                str(Fraction(ceil_div(self.hi, unit), EXPORT_SCALE))]


def arctan_inverse(n, terms):
    a = sum((Fraction((-1)**j, (2*j+1)*n**(2*j+1)) for j in range(terms)), Fraction())
    b = a+Fraction((-1)**terms, (2*terms+1)*n**(2*terms+1))
    return I(I.of(min(a, b)).lo, I.of(max(a, b)).hi)


PI = 16*arctan_inverse(5, 90)-4*arctan_inverse(239, 30)


def sin(x):
    x = I.of(x)
    if x.lo != x.hi:
        mid = (x.lo+x.hi)//2
        rad = max(mid-x.lo, x.hi-mid)
        out = sin(I(mid, mid))
        return I(max(-SCALE, out.lo-rad), min(SCALE, out.hi+rad))
    # Subtract an integer number of EXACT 2*pi periods, enclosed by PI.
    periods = (x.lo+x.hi)//(4*PI.lo)
    x = x-2*periods*PI
    if x.abs().hi > 7*SCALE:
        raise ValueError("trigonometric argument budget exceeded")
    out, term = x, x
    xx = x*x
    for j in range(1, 80):
        term = -term*xx/((2*j)*(2*j+1))
        out = out+term
    # Taylor polynomial through degree159; derivative remainder <= |x|^160/160!.
    rem = (x.abs()**160/factorial(160)).hi
    return I(max(-SCALE, out.lo-rem), min(SCALE, out.hi+rem))


def exp(x):
    x = I.of(x)
    if x.abs().hi > 4*SCALE:
        raise ValueError("exponential argument budget exceeded")
    out, term = I.of(1), I.of(1)
    for j in range(1, 100):
        term = term*x/j
        out = out+term
    # exp(4)<81, using e<3, bounds the order100 Lagrange remainder.
    rem = (81*x.abs()**100/factorial(100)).hi
    return I(out.lo-rem, out.hi+rem)


def dot(x, y):
    if len(x) != len(y):
        raise ValueError("dot-product dimensions differ")
    return sum((a*b for a, b in zip(x, y)), I.of(0))


def nonnegative(x):
    """Intersect with a separately proved nonnegative quantity's range."""
    if x.hi < 0:
        raise ValueError("nonnegative quantity has negative upper bound")
    return I(max(0, x.lo), x.hi)
