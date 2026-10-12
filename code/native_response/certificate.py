"""Directed-interval certificate of N4 mean nonclosure at two fixed angles.

This independently reconstructs the source over intervals containing
Q(sqrt(5),i); it imports neither source.py nor the numerical exponential.
Protected 3+3 endpoint blocks reduce each calculation to dimension 18.

For the source chart there are two absolute off-diagonal Gram entries in
each row, with four occurrences each. Proper carrier covariance and the
pentagonal link give these classes; representatives are the common-port
pairs ((0,2),(0,4)) and ((0,2),(0,8)). Thus a certified positive value of
b0-4|b1|-4|b2| makes the entire thirty-seam Gram positive definite.
See README.md for the locality reduction and symmetry argument.
"""
from __future__ import annotations

from math import factorial
from numbers import Integral

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext


DEGREE = 24
TAIL_DENOMINATOR = 2**50
REPRESENTATIVES = ((0, 2), (0, 4), (0, 8))


class IntervalCalculation:
    """Own the precision contexts so callers' mpmath settings cannot alter it."""

    def __init__(self):
        self.iv = MPIntervalContext()
        self.iv.dps = 35
        self.mp = MPContext()
        self.mp.dps = 80
        self.generators = self._source()
        self.center = [sum((d[s] for d in self.generators), self.iv.zeros(3))
                       for s in (0, 1)]

    def _source(self):
        iv = self.iv
        root = iv.sqrt(5)
        frames = []
        for phi in ((1 + root) / 2, (1 - root) / 2):
            frame = []
            for offset in range(3):
                for s in (1, -1):
                    for t in (1, -1):
                        v = [iv.mpf(0), iv.mpf(s), t * phi]
                        frame.append(v[offset:] + v[:offset])
            frames.append(frame)

        def hat(v):
            x, y, z = v
            return iv.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])

        result = []
        for v, vc in zip(*frames):
            norm = sum(x * x for x in v)
            quadrupole = iv.matrix([
                [v[i] * v[j] / 2 +
                 ((iv.mpf(1) / 12 - norm / 6) if i == j else 0)
                 for j in range(3)] for i in range(3)])
            result.append((hat([x / 2 for x in v]) + iv.j * quadrupole,
                           hat([x / 2 for x in vc])))
        return result

    def star(self, matrix):
        # Explicit components avoid reliance on the interval-complex
        # conjugate implementation and retain directed interval endpoints.
        iv = self.iv
        return iv.matrix([
            [iv.mpc(matrix[j, i].real, -matrix[j, i].imag)
             for j in range(matrix.rows)] for i in range(matrix.cols)])

    @staticmethod
    def trace(matrix):
        return sum(matrix[i, i] for i in range(matrix.rows))

    def kron(self, left, right):
        return self.iv.matrix([
            [left[i // right.rows, j // right.cols]
             * right[i % right.rows, j % right.cols]
             for j in range(left.cols * right.cols)]
            for i in range(left.rows * right.rows)])

    def generator(self, edge, first_block, second_block):
        iv = self.iv
        a, b = edge
        selected = (first_block, second_block)
        g = [self.generators[b][r] - self.generators[a][r] for r in selected]
        h = [self.generators[a][r] + self.generators[b][r] - self.center[r] / 6
             for r in selected]
        G = self.kron(g[0], iv.eye(3)) + self.kron(iv.eye(3), g[1])
        H = self.kron(h[0], iv.eye(3)) + self.kron(iv.eye(3), h[1])
        amplitude = H - iv.j * G
        opposite = -self.star(amplitude)
        result = iv.zeros(18)
        for i in range(9):
            for j in range(9):
                result[i, j + 9] = amplitude[i, j]
                result[i + 9, j] = opposite[i, j]
        return result

    def exponential(self, generator, theta):
        """Enclose exp(theta T), including source, rounding and Taylor error."""
        iv = self.iv
        argument = theta * generator
        # Both induced norms enclose the operator norm through
        # ||A||operator <= sqrt(||A||1 ||A||infinity). Checking both also
        # keeps this bound valid independently of the skew-Hermitian proof.
        for i in range(18):
            row = sum(abs(argument[i, j]) for j in range(18))
            column = sum(abs(argument[j, i]) for j in range(18))
            if not (row.b < 2 and column.b < 2):
                raise ArithmeticError("Taylor norm bound was not proved")
        # Sum_{k>24}2^k/k! <= exp(2)2^25/25! <8*2^25/25! <2^-50.
        if not 8 * 2**(DEGREE + 1) * TAIL_DENOMINATOR < factorial(DEGREE + 1):
            raise ArithmeticError("Taylor tail bound was not proved")
        result, term = iv.eye(18), iv.eye(18)
        for k in range(1, DEGREE + 1):
            term = (term * argument) / k
            result += term
        # An operator-norm error bounds the absolute error of every complex
        # entry, and therefore each real and imaginary component separately.
        error = iv.mpf([-1, 1]) / TAIL_DENOMINATOR
        for i in range(18):
            for j in range(18):
                x = result[i, j]
                result[i, j] = iv.mpc(x.real + error, x.imag + error)
        return result

    def local(self, edge, theta):
        iv = self.iv
        contrast = iv.diag([1] * 9 + [-1] * 9)
        outputs = {}
        for s in (0, 1):
            for t in (0, 1):
                unitary = self.exponential(self.generator(edge, s, t), theta)
                pullback = self.star(unitary) * contrast * unitary
                private = iv.matrix([
                    [(pullback[i, j] + pullback[i + 9, j + 9]) / 2
                     for j in range(9)] for i in range(9)])
                outputs[s, t] = (unitary * self.kron(iv.eye(2), private)
                                 * self.star(unitary))
        beta = sum(self.trace(contrast * q) for q in outputs.values()).real / 72
        reduced = [[iv.zeros(6), iv.zeros(6)] for _ in (0, 1)]
        for (s, t), output in outputs.items():
            for sign in (0, 1):
                residual = iv.matrix([
                    [output[i + sign * 9, j + sign * 9]
                     - (beta * (1 - 2 * sign) if i == j else 0)
                     for j in range(9)] for i in range(9)])
                for i in range(3):
                    for j in range(3):
                        reduced[0][sign][3 * s + i, 3 * s + j] += sum(
                            residual[3 * i + k, 3 * j + k] for k in range(3)) / 6
                        reduced[1][sign][3 * t + i, 3 * t + j] += sum(
                            residual[3 * k + i, 3 * k + j] for k in range(3)) / 6
        return beta, reduced

    def rational_enclosure(self, interval):
        """Serialize broad rational endpoints, checking outward containment."""
        denominator = 10**12
        lower = int(self.mp.floor(self.mp.mpf(interval.a) * denominator))
        upper = int(self.mp.ceil(self.mp.mpf(interval.b) * denominator))
        iv = self.iv
        if not (iv.mpf(lower) / denominator <= interval.a and
                iv.mpf(upper) / denominator >= interval.b):
            raise ArithmeticError("serialized endpoints do not enclose result")
        return {"lower_numerator": lower, "upper_numerator": upper,
                "denominator": denominator}


def certify(theta_denominator: int) -> dict:
    """Certify the fixed source at theta=1/16 or theta=1/8 in about 15 seconds."""
    if (isinstance(theta_denominator, bool) or
            not isinstance(theta_denominator, Integral) or
            theta_denominator not in (8, 16)):
        raise ValueError("certificate covers only theta=1/8 and theta=1/16")
    calculation = IntervalCalculation()
    iv = calculation.iv
    theta = iv.mpf(1) / int(theta_denominator)
    values = [calculation.local(edge, theta) for edge in REPRESENTATIVES]
    beta, reduced = values[0]
    diagonal = (iv.mpf(10) / 33) * beta * (1 - beta)
    cross = []
    for _, other in values[1:]:
        value = (9 * calculation.trace(reduced[0][0] * other[0][0])
                 + calculation.trace(reduced[0][1] * other[0][1])).real / (66 * 6)
        cross.append(value)
    margin = diagonal - 4 * abs(cross[0]) - 4 * abs(cross[1])
    if not margin.a > 0:
        raise ArithmeticError("positive full-network Gram margin not certified")
    names = ("beta", "diagonal", "cross_first", "cross_second", "gershgorin_margin")
    enclosures = {name: calculation.rational_enclosure(value)
                  for name, value in zip(names, (beta, diagonal, *cross, margin))}
    return {"theta_numerator": 1, "theta_denominator": int(theta_denominator),
            "taylor_degree": DEGREE, "entry_tail_denominator": TAIL_DENOMINATOR,
            "interval_decimal_precision": 35,
            "representative_edges": [list(edge) for edge in REPRESENTATIVES],
            "off_diagonal_multiplicities": [4, 4], "enclosures": enclosures}


def certify_fixed_angles() -> list[dict]:
    return [certify(denominator) for denominator in (16, 8)]


if __name__ == "__main__":
    import json
    print(json.dumps(certify_fixed_angles(), indent=2, sort_keys=True))
