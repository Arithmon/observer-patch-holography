"""Exact rational intervals; no binary floating arithmetic in the bounds."""
from fractions import Fraction as Q


class I:
    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = Q(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError('reversed interval')

    def __add__(self, other):
        y = other if isinstance(other, I) else I(other)
        return I(self.lo + y.lo, self.hi + y.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -(other if isinstance(other, I) else I(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        y = other if isinstance(other, I) else I(other)
        p = [a*b for a in (self.lo, self.hi) for b in (y.lo, y.hi)]
        return I(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(self, other):
        y = other if isinstance(other, I) else I(other)
        if y.lo <= 0 <= y.hi:
            raise ValueError('division by an interval containing zero')
        return self * I(1/y.hi, 1/y.lo)

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError('nonnegative integer power required')
        if n == 0:
            return I(1)
        if n % 2 == 0 and self.lo <= 0 <= self.hi:
            return I(0, max(abs(self.lo), abs(self.hi))**n)
        p = (self.lo**n, self.hi**n)
        return I(min(p), max(p))

    def pair(self):
        return [str(self.lo), str(self.hi)]


def phi():
    lower = Q('1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374')
    upper = lower + Q(1, 10**100)
    if not (1 < lower and lower*lower-lower-1 < 0 < upper*upper-upper-1):
        raise ValueError('golden root bracket')
    return I(lower, upper)


def magnitude(x):
    return max(abs(x.lo), abs(x.hi))
