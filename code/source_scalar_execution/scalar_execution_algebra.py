"""Exact arithmetic in Q(phi), with canonical bounded receipt parsing."""
from dataclasses import dataclass
from fractions import Fraction as F


@dataclass(frozen=True)
class Q:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        if type(self.a) not in (int, F) or type(self.b) not in (int, F):
            raise ValueError('exact rational algebraic coefficients required')
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def of(x):
        if isinstance(x, Q):
            return x
        if type(x) not in (int, F):
            raise ValueError('exact scalar required')
        return Q(x)

    def __add__(self, other):
        other = Q.of(other)
        return Q(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.a, -self.b)

    def __sub__(self, other):
        return self+-Q.of(other)

    def __rsub__(self, other):
        return Q.of(other)+-self

    def __mul__(self, other):
        other = Q.of(other)
        return Q(self.a*other.a+self.b*other.b,
                 self.a*other.b+self.b*other.a+self.b*other.b)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Q.of(other)
        norm = other.a*other.a+other.a*other.b-other.b*other.b
        if not norm:
            raise ValueError('division by zero')
        return self*Q((other.a+other.b)/norm, -other.b/norm)

    def __rtruediv__(self, other):
        return Q.of(other)/self

    def __pow__(self, n):
        if type(n) is not int or n < 0:
            raise ValueError('nonnegative integer power required')
        result = Q(1)
        for _ in range(n):
            result *= self
        return result

    def sign(self):
        c, b = 2*self.a+self.b, self.b
        if b == 0:
            return (c > 0)-(c < 0)
        if c >= 0 and b > 0:
            return 1
        if c <= 0 and b < 0:
            return -1
        d = c*c-5*b*b
        return ((d > 0)-(d < 0))*(1 if c > 0 else -1)

    def encode(self):
        return [str(self.a), str(self.b)]


def parse(value):
    if type(value) is not list or len(value) != 2:
        raise ValueError('algebraic value must have two rational strings')
    coeff = []
    for item in value:
        if type(item) is not str or len(item) > 3000:
            raise ValueError('invalid algebraic coefficient type or size')
        a = F(item)
        if str(a) != item or a.numerator.bit_length() > 10000 or a.denominator.bit_length() > 10000:
            raise ValueError('noncanonical or oversized algebraic coefficient')
        coeff.append(a)
    return Q(*coeff)


def rational_bounds(value, denominator=10**12):
    """Integer bisection using exact signs, never a floating conversion."""
    value = Q.of(value)
    if type(denominator) is not int or denominator <= 0:
        raise ValueError('positive exact denominator required')
    lo, hi = -denominator, denominator
    while (value-Q(F(lo, denominator))).sign() < 0:
        lo *= 2
    while (value-Q(F(hi, denominator))).sign() > 0:
        hi *= 2
    while hi-lo > 1:
        mid = (hi+lo)//2
        if (value-Q(F(mid, denominator))).sign() >= 0:
            lo = mid
        else:
            hi = mid
    return [str(F(lo, denominator)), str(F(hi, denominator))]
