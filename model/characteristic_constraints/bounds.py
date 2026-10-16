"""Exact rational witnesses for the standard two-unitary Gram constraint."""
from fractions import Fraction as F
from math import factorial

import numpy as np


# A complex rational is (real, imaginary). No rounded datum enters a certificate.
def add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def mul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def conj(a):
    return (a[0], -a[1])


def norm2(a):
    return a[0]**2+a[1]**2


def trig_interval(x, cosine=False):
    """Degree-120 Taylor polynomial with derivative-bounded remainder.

    This is an outward rational interval, not a floating accuracy estimate.
    The declared laboratory phases have |x| <= 25.
    """
    if type(x) is not F or abs(x) > 25 or type(cosine) is not bool:
        raise ValueError('rational phase in [-25,25] and Boolean parity required')
    terms = range(0 if cosine else 1, 121, 2)
    center = sum(((-1)**(k//2)*x**k/F(factorial(k)) for k in terms), F(0))
    error = abs(x)**121/F(factorial(121))
    # Compact exact outward endpoints; the rounding is also part of the bound.
    grid = 10**40
    lo, hi = center-error, center+error
    return F((lo*grid).__floor__(), grid), F((hi*grid).__ceil__(), grid)


def scale_interval(interval, factor):
    values = [factor*endpoint for endpoint in interval]
    return min(values), max(values)


def phase_real_interval(value, phase):
    c = scale_interval(trig_interval(phase, True), value[0])
    s = scale_interval(trig_interval(phase), -value[1])
    return c[0]+s[0], c[1]+s[1]


def gram(a, b, z, phase):
    a, b, z = (complex(float(x[0]), float(x[1])) for x in (a, b, z))
    w = np.exp(1j*float(phase))*z
    return np.array([[1, a.conjugate(), b], [a, 1, w],
                     [b.conjugate(), w.conjugate(), 1]], complex)


def witness_floor(a, b, z, phase):
    """Certify an entrywise complex-mean correction floor using a rational v.

    Float arithmetic chooses a witness only. Its quadratic form and bound
    are independently evaluated as rational intervals. The Schur witness
    avoids platform-dependent eigenvectors; it is not necessarily optimal.
    """
    for value in (a, b, z):
        if type(value) is not tuple or len(value) != 2 or any(type(x) is not F for x in value):
            raise ValueError('exact complex rational means required')
    if norm2(a) >= 1 or norm2(b) >= 1:
        raise ValueError('this witness requires strictly interior axis means')
    g = gram(a, b, z, phase)
    af, bf = g[1, 0], g[0, 2]
    t = -(g[1, 2]-af*bf).conjugate()/(1-abs(bf)**2)
    vector = [-af.conjugate()-bf*t, 1+0j, t]
    v = tuple((F(round(x.real*10**6), 10**6), F(round(x.imag*10**6), 10**6)) for x in vector)
    n = sum(map(norm2, v), F(0))
    constant = n+2*(mul(mul(conj(v[0]), conj(a)), v[1])[0]
                    +mul(mul(conj(v[0]), b), v[2])[0])
    rotating = mul(mul(conj(v[1]), z), v[2])
    lo, hi = phase_real_interval(rotating, phase)
    qlo, qhi = constant+2*lo, constant+2*hi
    # |v* DeltaG v| <= 2 epsilon sum_i<j |vi||vj| <= 2 epsilon ||v||².
    floor = max(F(0), -qhi/(2*n))
    return {'vector': [[str(x), str(y)] for x, y in v],
            'norm_squared': str(n), 'quadratic_form_interval': [str(qlo), str(qhi)],
            'complex_entry_correction_lower_bound': str(floor),
            'complex_entry_correction_lower_bound_float': float(floor)}
