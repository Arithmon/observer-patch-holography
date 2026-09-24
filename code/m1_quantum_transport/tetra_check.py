"""Independent exact Q(sqrt(3),i) and Grassmann reconstruction of the optimizer."""

from dataclasses import dataclass
from fractions import Fraction as F
import itertools

import mpmath as mp

from .check import need, digest, decimal


@dataclass(frozen=True)
class K:
    a: F = F(0)
    b: F = F(0)
    c: F = F(0)
    d: F = F(0)

    @staticmethod
    def lift(x):
        return x if isinstance(x, K) else K(F(x))

    def __add__(self, other):
        other = K.lift(other)
        return K(*(x+y for x, y in zip(self.parts(), other.parts())))

    __radd__ = __add__

    def __neg__(self):
        return K(*(-x for x in self.parts()))

    def __sub__(self, other):
        return self+-K.lift(other)

    def __rsub__(self, other):
        return K.lift(other)+-self

    def __mul__(self, other):
        e, f, g, h = K.lift(other).parts()
        a, b, c, d = self.parts()
        return K(a*e+3*b*f-c*g-3*d*h, a*f+b*e-c*h-d*g,
                 a*g+3*b*h+c*e+3*d*f, a*h+b*g+c*f+d*e)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = F(scalar)
        return K(*(x/scalar for x in self.parts()))

    def __bool__(self):
        return any(self.parts())

    def parts(self):
        return self.a, self.b, self.c, self.d

    def conjugate(self):
        return K(self.a, self.b, -self.c, -self.d)


J, R = K(c=F(1)), K(b=F(1))


def eye(n):
    return [[K(int(i == j)) for j in range(n)] for i in range(n)]


def add(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, z):
    return [[z*x for x in row] for row in a]


def mul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), K())
             for j in range(len(b[0]))] for i in range(len(a))]


def dagger(a):
    return [[a[j][i].conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def diagonal(values):
    return [[K.lift(x) if i == j else K() for j in range(len(values))] for i, x in enumerate(values)]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), K())


def encoded(a):
    return [[[str(x) for x in z.parts()] for z in row] for row in a]


def frames():
    # Complex conference matrix; no candidate matrices or symbolic solver are imported.
    s = [[K.lift(x) for x in row] for row in
         [[0, 1, 1, 1], [1, 0, J, -J], [1, -J, 0, J], [1, J, -J, 0]]]
    need(mul(s, s) == scale(eye(4), 3) and dagger(s) == s, "conference matrix identity")
    return scale(s, R/3), [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]


def bands(coin, points):
    answer = {}
    for sign in (1, -1):
        p = scale(add(eye(4), scale(coin, sign)), F(1, 2))
        need(mul(p, p) == p and trace(p) == K(2), "two-dimensional native band")
        x = [scale(mul(mul(p, diagonal([point[i] for point in points])), p), R) for i in range(3)]
        for i in range(3):
            need(dagger(x[i]) == x[i] and mul(x[i], x[i]) == p, "compressed Pauli square")
            need(mul(x[i], x[(i+1) % 3]) == scale(x[(i+2) % 3], sign*J), "compressed chirality")
        effects = [mul(mul(p, diagonal([int(i == j) for i in range(4)])), p) for j in range(4)]
        total = [[K()]*4 for _ in range(4)]
        for point, effect in zip(points, effects):
            expected = p
            for i in range(3):
                expected = add(expected, scale(x[i], point[i]*R/3))
            need(effect == scale(expected, F(1, 4)), "minimal tetrahedral positive effect")
            total = add(total, effect)
        need(total == p, "all tetrahedral outcomes retained")
        for i, a in enumerate(effects):
            for j, b in enumerate(effects):
                need(trace(mul(a, b)) == K(F(1, 4) if i == j else F(1, 12)), "tetrahedral overlap")
        answer[str(sign)] = dict(projector=encoded(p), rank=2, chirality=sign,
                                 pauli=[encoded(y) for y in x], effects=[encoded(e) for e in effects])
    return answer


def spectrum(coin, points):
    rows, zeros, pi_points, middle = [], [], [], []
    phases = [K(1), -J, K(-1), J]
    for k in itertools.product(range(4), repeat=3):
        u = mul(diagonal([phases[sum(x*y for x, y in zip(point, k)) % 4] for point in points]), coin)
        need(mul(dagger(u), u) == eye(4), "complete tetrahedral momentum unitarity")
        u2 = mul(u, u)
        middle_coefficient = F(-2, 3)*sum((-1)**x for x in k)
        # Verify the full quartic identity, trace and eigenvalue degeneracies.
        need(add(add(mul(u2, u2), scale(u2, middle_coefficient)), eye(4)) == [[K()]*4 for _ in range(4)],
             "tetrahedral spectral polynomial")
        need(trace(u) == K() and trace(u2) == K(-2*middle_coefficient), "all four spectral roots")
        coefficients = [K(1), K(), K(middle_coefficient), K(), K(1)]
        rows.append([list(k), encoded(u), [[str(x) for x in z.parts()] for z in coefficients]])
        if all(x % 2 == 0 for x in k):
            zeros.append(list(k))
            pi_points.append(list(k))
        if all(x % 2 == 1 for x in k):
            middle.append(list(k))
    need(len(zeros) == len(pi_points) == len(middle) == 8, "complete tetrahedral node catalog")
    return dict(q=4, momenta=64, modes=256, rows_sha256=digest(rows), zero_points=zeros,
                pi_points=pi_points, middle_points=middle)


def fock(coin):
    # Expand transformed creation operators, with anticommutation at each insertion.
    full = [[K()]*16 for _ in range(16)]
    for source in range(16):
        polynomial = {0: K(1)}
        for column in range(4):
            if not source & (1 << column):
                continue
            changed = {}
            for state, coefficient in polynomial.items():
                for row in range(4):
                    if state & (1 << row):
                        continue
                    sign = (-1)**sum(bool(state & (1 << i)) for i in range(row+1, 4))
                    target = state | (1 << row)
                    changed[target] = changed.get(target, K())+sign*coefficient*coin[row][column]
            polynomial = changed
        for target, coefficient in polynomial.items():
            full[target][source] = coefficient
    need(mul(dagger(full), full) == eye(16), "all-sector minimal coin unitarity")
    return dict(dimension=16, sector_dimensions=[1, 4, 6, 4, 1], matrix_sha256=digest(encoded(full)))


def thermal():
    result = {}
    with mp.workdps(75):
        for q in (4, 8, 12, 16):
            cosine = [mp.mpf((-1)**(4*k//q)) if 4*k % q == 0 else mp.cospi(mp.mpf(4)*k/q)
                      for k in range(q)]
            energies, zeros = [], 0
            for k in itertools.product(range(q), repeat=3):
                arg = sum(cosine[i] for i in k)/3
                need(-1 <= arg <= 1, "tetrahedral spectral angle domain")
                epsilon = mp.acos(arg)/2
                zeros += epsilon == 0
                energies.extend([epsilon*q/mp.sqrt(3), (mp.pi-epsilon)*q/mp.sqrt(3)])
            cases = {}
            for beta in (1, 2):
                partition = mp.fsum([2*mp.log1p(mp.exp(-beta*e)) for e in energies])
                energy = mp.fsum([2*e*mp.exp(-beta*e)/(1+mp.exp(-beta*e)) for e in energies])
                cases[str(beta)] = [decimal(x) for x in (partition, energy, partition+beta*energy)]
            result[str(q)] = dict(momenta=q**3, modes=4*q**3, zero_momenta=zeros, cases=cases)
    return result


def reconstruct():
    coin, points = frames()
    trajectories = []
    for point in points:
        for f in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
            time = "0" if f == 0 else "sqrt(3)" if f == 1 else f"{f.numerator}*sqrt(3)/{f.denominator}"
            if f.numerator == 1 and f not in (0, 1):
                time = f"sqrt(3)/{f.denominator}"
            trajectories.append([list(point), str(f), [str(f*x) for x in point], time])
    signals = []
    for channel, point in enumerate(points):
        prepared = [[row[channel]] for row in coin]
        output = mul(coin, prepared)
        need(output == [[K(int(i == channel))] for i in range(4)], "lossless native occupation record")
        signals.append(dict(channel=channel, prepared=encoded(prepared), output=encoded(output),
                            probabilities=[str(int(i == channel)) for i in range(4)],
                            empty_probabilities=["0"]*4, destination=list(point), time="sqrt(3)"))
    return dict(points=[list(x) for x in points], coin=encoded(coin), bands=bands(coin, points),
                spectrum=spectrum(coin, points), flight_duration="sqrt(3)", native_speed="1", field_speed="1/3",
                trajectories_sha256=digest(trajectories), signals=signals, fock=fock(coin), thermal=thermal())
