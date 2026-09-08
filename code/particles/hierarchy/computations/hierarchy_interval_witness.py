#!/usr/bin/env python3
"""Independent interval witness for the R_U hierarchy root (issue #518).

Recomputes the hierarchy zero-crossing certificate with genuine interval
arithmetic (mpmath.iv, outward rounding) and a certified derivative enclosure
obtained by forward-mode dual numbers over intervals. Independent of the
existing hierarchy_recompute.py float evaluation and of the shipped
certificate JSON: same declared formula stack, different arithmetic and
different code path.

Witness content:
  W1 (bracket):    Phi(lo) and Phi(hi) are enclosed by intervals of opposite
                   definite sign, evaluated with outward rounding.
  W2 (derivative): Phi'(I) over each subdivision piece of I_U is enclosed by
                   an interval with definite negative sign, via interval dual
                   numbers, so Phi is strictly decreasing on I_U and the root
                   is unique.
  W3 (perturbation controls): perturbing a defining antecedent (beta
   coefficient, weak multiplicity, pixel value) must break the certificate.
"""
import json
import sys
from fractions import Fraction
import hashlib
from pathlib import Path
import mpmath as mp
from mpmath import iv
from mpmath.libmp import to_rational

CODE_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(CODE_ROOT))
from interval_decimal import interval_json  # noqa: E402

iv.dps = 60

# Fixed antecedents from the frozen packet; no values fitted in this computation.
P_STR = '1.630968209403959324879279847782648941'
I_LO = '0.041123336195630494'
I_HI = '0.041125336195630496'
N2, N3 = 128, 64


class Dual:
    """Forward-mode dual number over mpmath intervals: value + derivative."""

    __slots__ = ('v', 'd')

    def __init__(self, v, d=None):
        self.v = iv.mpf(v) if not isinstance(v, type(iv.mpf(0))) else v
        zero = iv.mpf(0)
        self.d = zero if d is None else (iv.mpf(d) if not isinstance(d, type(iv.mpf(0))) else d)

    @staticmethod
    def const(v):
        return Dual(v, 0)

    @staticmethod
    def var(v):
        return Dual(v, 1)

    def __add__(self, o):
        o = o if isinstance(o, Dual) else Dual.const(o)
        return Dual(self.v + o.v, self.d + o.d)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, Dual) else Dual.const(o)
        return Dual(self.v - o.v, self.d - o.d)

    def __rsub__(self, o):
        return Dual.const(o).__sub__(self)

    def __mul__(self, o):
        o = o if isinstance(o, Dual) else Dual.const(o)
        return Dual(self.v * o.v, self.d * o.v + self.v * o.d)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, Dual) else Dual.const(o)
        return Dual(self.v / o.v, (self.d * o.v - self.v * o.d) / (o.v * o.v))

    def __rtruediv__(self, o):
        return Dual.const(o).__truediv__(self)

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def exp(self):
        e = iv.exp(self.v)
        return Dual(e, e * self.d)

    def log(self):
        return Dual(iv.log(self.v), self.d / self.v)

    def sqrt(self):
        s = iv.sqrt(self.v)
        return Dual(s, self.d / (2 * s))


def dexp(x):
    return x.exp() if isinstance(x, Dual) else iv.exp(x)


def dlog(x):
    return x.log() if isinstance(x, Dual) else iv.log(x)


def dsqrt(x):
    return x.sqrt() if isinstance(x, Dual) else iv.sqrt(x)


def phi(a, P, b1, b2, b3, weak_mult=4, mu_iters=64):
    """The declared hierarchy readback Phi(a): heat-kernel edge entropies minus P/4.

    `a` may be a Dual (interval dual number) or an interval. All constants are
    promoted through C() so no mpmath-side reflected operator ever fires on a
    Dual operand.
    """
    dual_mode = isinstance(a, Dual)

    def C(x):
        x = iv.mpf(x) if not hasattr(x, 'a') else x
        return Dual.const(x) if dual_mode else x

    one = iv.mpf(1)
    pi = iv.pi
    MU = iv.exp(-2 * pi) * iv.mpf(P) ** (one / 6)
    v = C(iv.mpf(P) ** (-one / 2)) * dexp(C(-2 * pi) / (C(weak_mult) * a))

    def alpha(mu, b):
        return C(1) / (C(1) / a + C(b / (2 * pi)) * dlog(C(MU) / mu))

    def fmu(mu):
        return v / C(2) * dsqrt(
            C(4 * pi) * alpha(mu, b2) + C(4 * pi) * C(iv.mpf(3) / 5) * alpha(mu, b1))

    mu = v * C(iv.mpf('0.3714'))
    for _ in range(mu_iters):
        mu = fmu(mu)

    t2 = C(4 * pi ** 2) * alpha(mu, b2)
    t3 = C(4 * pi ** 2) * alpha(mu, b3)

    def ell_su2(t):
        Z = C(0)
        S = C(0)
        for n in range(N2 + 1):
            j = iv.mpf(n) / 2
            d = 2 * j + 1
            Cas = j * (j + 1)
            w = C(d) * dexp(-(t * C(Cas)))
            Z = Z + w
            S = S + w * C(iv.log(d))
        return S / Z

    def ell_su3(t):
        Z = C(0)
        S = C(0)
        for p in range(N3 + 1):
            for q in range(N3 + 1):
                d = iv.mpf((p + 1) * (q + 1) * (p + q + 2)) / 2
                Cas = iv.mpf(p * p + q * q + p * q + 3 * p + 3 * q) / 3
                w = C(d) * dexp(-(t * C(Cas)))
                Z = Z + w
                S = S + w * C(iv.log(d))
        return S / Z

    return ell_su2(t2) + ell_su3(t3) - C(iv.mpf(P) / 4)


def interval_str(x):
    """Export genuine outer endpoints, independently of scalar mp precision."""
    bounds = interval_json(x, 40)
    return f"[{bounds['lo']}, {bounds['hi']}]"


def subdivision(pieces=8):
    """Exact rational knots and outward cells covering the declared interval."""
    if type(pieces) is not int or pieces < 1:
        raise ValueError("pieces must be a positive integer")
    lo, hi = Fraction(I_LO), Fraction(I_HI)
    knots = [lo + (hi - lo) * k / pieces for k in range(pieces + 1)]
    enclosed = [iv.mpf(q.numerator) / iv.mpf(q.denominator) for q in knots]
    cells = [iv.mpf([a.a, b.b]) for a, b in zip(enclosed, enclosed[1:])]
    for k, cell in enumerate(cells):
        left, right = (Fraction(*to_rational(x)) for x in cell._mpi_)
        if not left <= knots[k] < knots[k + 1] <= right:
            raise ValueError("outward subdivision lost exact rational coverage")
    return knots, cells


def sign_of(x):
    if x.b < 0:
        return -1
    if x.a > 0:
        return 1
    return 0


def _produce_at_working_precision():
    b1, b2, b3 = iv.mpf(33) / 5, iv.mpf(1), iv.mpf(-3)
    lo, hi = iv.mpf(I_LO), iv.mpf(I_HI)

    # W1: outward-rounded endpoint enclosures with opposite definite signs.
    f_lo = phi(lo, P_STR, b1, b2, b3)
    f_hi = phi(hi, P_STR, b1, b2, b3)
    w1_ok = sign_of(f_lo) * sign_of(f_hi) == -1

    # W2: certified derivative enclosure over a subdivision of I_U.
    pieces = 8
    knots, cells = subdivision(pieces)
    derivs = []
    w2_ok = True
    for a_piece in cells:
        res = phi(Dual.var(a_piece), P_STR, b1, b2, b3)
        derivs.append(res.d)
        if not (res.d.b < 0):
            w2_ok = False

    # W3: perturbation controls — each must break W1 or W2.
    controls = {}
    f_lo_p = phi(lo, P_STR, b1, b2 + iv.mpf('0.05'), b3)
    f_hi_p = phi(hi, P_STR, b1, b2 + iv.mpf('0.05'), b3)
    controls['beta2_perturbed_bracket_broken'] = not (sign_of(f_lo_p) * sign_of(f_hi_p) == -1)
    f_lo_w = phi(lo, P_STR, b1, b2, b3, weak_mult=5)
    f_hi_w = phi(hi, P_STR, b1, b2, b3, weak_mult=5)
    controls['weak_mult_perturbed_bracket_broken'] = not (sign_of(f_lo_w) * sign_of(f_hi_w) == -1)
    P_pert = '1.7309682094039593'
    f_lo_P = phi(lo, P_pert, b1, b2, b3)
    f_hi_P = phi(hi, P_pert, b1, b2, b3)
    controls['pixel_perturbed_bracket_broken'] = not (sign_of(f_lo_P) * sign_of(f_hi_P) == -1)

    derivative_union = iv.mpf([min(d.a for d in derivs), max(d.b for d in derivs)])

    out = {
        'witness': 'independent_hierarchy_interval_witness',
        'issue': 518,
        'arithmetic': f'mpmath.iv outward-rounded intervals, dps={iv.dps}, forward-mode dual derivative',
        'interval_I_U': [I_LO, I_HI],
        'phi_at_lower': interval_str(f_lo),
        'phi_at_upper': interval_str(f_hi),
        'bracket_sign_change': bool(w1_ok),
        'derivative_enclosure_pieces': pieces,
        'derivative_subdivision': {
            'exact_rational_knots': [str(k) for k in knots],
            'outward_cells': [interval_str(cell) for cell in cells],
            'derivative_enclosures': [interval_str(d) for d in derivs],
            'entire_declared_interval_covered': True,
        },
        'derivative_enclosure_union': interval_str(derivative_union),
        'derivative_strictly_negative': bool(w2_ok),
        'unique_root_certified': bool(w1_ok and w2_ok),
        'perturbation_controls': controls,
        'controls_all_fail_closed': all(controls.values()),
        'formula_scope': {
            'pixel_decimal': P_STR,
            'su2_max_n': N2,
            'su3_max_p_and_q': N3,
            'mu_iterations': 64,
            'infinite_sum_or_fixed_point_limit_certified_here': False,
            'physical_hierarchy_attachment': False,
        },
        'endpoint_serialization': 'exact binary endpoints rounded down/up to 40 decimal significant digits',
        'source_pins': {
            path.relative_to(CODE_ROOT.parent).as_posix(): {
                'bytes': path.stat().st_size,
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in (Path(__file__).resolve(), CODE_ROOT / 'interval_decimal.py')
        },
        'historical_correction': {
            'prior_receipt_git_object': '5e1816846d14dbd73cdeeb3e8ab596099cc06076:code/particles/hierarchy/certificates/independent_interval_witness_receipt.json',
            'prior_receipt_sha256': 'ae6262c86b592796b1574f3c58ae4ca44df1516fa169c2b1c6ef6f735a5b2291',
            'prior_producer_sha256': '8823c747ace9a45737604df5830d5d31b682c1b86fec88250aeae94e6f36929b',
            'reason': 'Scalar rounding did not preserve outer endpoints, and scalar-rounded subdivision omitted part of the declared interval; this is a fresh finite-formula witness.',
        },
    }
    return out


def produce():
    """Use the declared interval precision without changing caller settings."""
    previous = iv.dps
    try:
        iv.dps = 60
        return _produce_at_working_precision()
    finally:
        iv.dps = previous


def main():
    out = produce()
    print(json.dumps(out, indent=2))
    return 0 if (out['unique_root_certified'] and out['controls_all_fail_closed']) else 1


if __name__ == '__main__':
    sys.exit(main())
