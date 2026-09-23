"""Independent rational-field, path, exterior-power and spectral reconstruction."""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math

import mpmath as mp


def need(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def same(a, b, message):
    need(canonical(a) == canonical(b), message)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


@dataclass(frozen=True)
class Q:
    """Exact a+b*sqrt(2)."""
    a: F = F(0)
    b: F = F(0)

    @staticmethod
    def lift(x):
        return x if isinstance(x, Q) else Q(F(x))

    def __add__(self, other):
        other = self.lift(other)
        return Q(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.a, -self.b)

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        return Q(self.a*other.a+2*self.b*other.b, self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.lift(other)
        denominator = other.a**2-2*other.b**2
        need(denominator != 0, "division by zero in exact field")
        return self*Q(other.a/denominator, -other.b/denominator)

    def __bool__(self):
        return bool(self.a or self.b)


@dataclass(frozen=True)
class Z:
    """Exact Gaussian extension of the real quadratic field."""
    re: Q = Q()
    im: Q = Q()

    @staticmethod
    def lift(x):
        return x if isinstance(x, Z) else Z(Q.lift(x))

    def __add__(self, other):
        other = self.lift(other)
        return Z(self.re+other.re, self.im+other.im)

    __radd__ = __add__

    def __neg__(self):
        return Z(-self.re, -self.im)

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        return Z(self.re*other.re-self.im*other.im, self.re*other.im+self.im*other.re)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.lift(other)
        denominator = other.re*other.re+other.im*other.im
        result = self*other.conjugate()
        return Z(result.re/denominator, result.im/denominator)

    def conjugate(self):
        return Z(self.re, -self.im)

    def __bool__(self):
        return bool(self.re or self.im)


ZERO, ONE, J = Z(), Z.lift(1), Z(Q(), Q(1))


def eye(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def mat(rows):
    return [[Z.lift(x) for x in row] for row in rows]


SIGMA = (mat([[0, 1], [1, 0]]), mat([[0, -J], [J, 0]]), mat([[1, 0], [0, -1]]))


def add(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, k):
    return [[k*x for x in row] for row in a]


def mul(a, b):
    out = [[ZERO]*len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(b[k]):
                    if y:
                        out[i][j] += x*y
    return out


def dagger(a):
    return [[a[j][i].conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), ZERO)


def encoded(a):
    need(all(x.re.b == 0 and x.im.b == 0 for row in a for x in row), "expected rational matrix")
    return [[[str(x.re.a), str(x.im.a)] for x in row] for row in a]


@lru_cache(None)
def determinant(a):
    if not a:
        return ONE
    if len(a) == 1:
        return a[0][0]
    # Choose the sparsest row, unlike the producer's gate-state updates.
    r = min(range(len(a)), key=lambda i: sum(bool(x) for x in a[i]))
    total = ZERO
    for c, x in enumerate(a[r]):
        if x:
            minor = tuple(tuple(v for j, v in enumerate(row) if j != c)
                          for i, row in enumerate(a) if i != r)
            total += (-1)**(r+c)*x*determinant(minor)
    return total


def native_word(order):
    projectors = {(axis, sign): scale(add(eye(2), scale(SIGMA[axis], sign)), F(1, 2))
                  for axis in range(3) for sign in (-1, 1)}
    history = []
    final = {}
    for length in range(1, 4):
        prefixes = []
        for signs in itertools.product((-1, 1), repeat=length-1):
            source = [0, 0, 0]
            amplitude = eye(2)
            for axis, sign in zip(order, signs):
                source[axis] += sign
                amplitude = mul(projectors[axis, sign], amplitude)
            prefixes.append((source, amplitude))
        for source, amplitude in sorted(prefixes):
            axis = order[length-1]
            for sign in (-1, 1):
                target = list(source)
                target[axis] += sign
                out = mul(projectors[axis, sign], amplitude)
                positions = []
                for numerator in range(5):
                    fraction = F(numerator, 4)
                    pos = [F(x) for x in source]
                    pos[axis] += sign*fraction
                    positions.append([str(length-1+fraction), [str(x) for x in pos]])
                history.append(dict(stage=length, source=source, target=target, sign=sign,
                                    amplitude=encoded(out), positions=positions))
                if length == 3:
                    final[tuple(target)] = out
    # Exact convolution identity proves the whole translation-invariant unitary.
    convolution = {}
    for x, a in final.items():
        for y, b in final.items():
            displacement = tuple(y[i]-x[i] for i in range(3))
            convolution[displacement] = add(convolution.get(displacement, mat([[0, 0], [0, 0]])),
                                           mul(dagger(a), b))
    need(convolution.pop((0, 0, 0)) == eye(2), "native word lost norm")
    need(all(not z for matrix in convolution.values() for row in matrix for z in row),
         "native word is not a unitary convolution")
    return dict(order=list(order), duration="3", speed_cap="1", effective_speed="1/3",
                intermediate_events=len(history), history_sha256=digest(history),
                kernel=[[list(k), encoded(v)] for k, v in sorted(final.items())])


def spectral_catalog(q):
    cosine = ([Q(1), Q(), Q(-1), Q()] if q == 4 else
              [Q(1), Q(0, F(1, 2)), Q(), Q(0, F(-1, 2)), Q(-1),
               Q(0, F(-1, 2)), Q(), Q(0, F(1, 2))])
    sine = [cosine[(i-q//4) % q] for i in range(q)]
    rows, zeros, pi_modes = [], [], []
    for k in itertools.product(range(q), repeat=3):
        cx, cy, cz = (cosine[i] for i in k)
        sx, sy, sz = (sine[i] for i in k)
        coefficients = [cx*cy*cz+sx*sy*sz, sx*cy*cz-cx*sy*sz,
                        cx*sy*cz+sx*cy*sz, cx*cy*sz-sx*sy*cz]
        need(sum((v*v for v in coefficients), Q()) == Q(1), "complete Bloch spectrum not unitary")
        if coefficients[0] == Q(1):
            zeros.append(list(k))
        if coefficients[0] == Q(-1):
            pi_modes.append(list(k))
        rows.append([list(k), [[str(x.a), str(x.b)] for x in coefficients]])
    need(len(zeros) == len(pi_modes) == 8, "all eight cones at each eigenphase are required")
    return dict(q=q, momenta=q**3, modes=2*q**3, rows=rows,
                zero_points=zeros, pi_points=pi_modes)


def nodes():
    out = []
    for index in itertools.product(range(4), repeat=3):
        if all(i % 2 == 0 for i in index) and sum(i//2 for i in index) % 2 == 0:
            sign = 1
        elif all(i % 2 == 1 for i in index) and sum(i//2 for i in index) % 2 == 0:
            sign = -1
        else:
            continue
        out.append(dict(index=list(index), chirality=sign,
                        derivative=[["1", "0", "0"], ["0", str(sign), "0"], ["0", "0", "1"]]))
    return out


def decimal(x):
    return mp.nstr(x, 21, strip_zeros=False)


def thermal():
    answer = {}
    with mp.workdps(75):
        for q in (4, 8, 12, 16):
            cosine, sine = [], []
            for k in range(q):
                if 4*k % q == 0:
                    quarter = (4*k//q) % 4
                    cosine.append(mp.mpf([1, 0, -1, 0][quarter]))
                    sine.append(mp.mpf([0, 1, 0, -1][quarter]))
                else:
                    cosine.append(mp.cospi(mp.mpf(2)*k/q))
                    sine.append(mp.sinpi(mp.mpf(2)*k/q))
            energies = []
            for x, y, z in itertools.product(range(q), repeat=3):
                scalar = cosine[x]*cosine[y]*cosine[z]+sine[x]*sine[y]*sine[z]
                need(-1 <= scalar <= 1, "spectral angle domain")
                energies.append(mp.acos(scalar)*q/3)
            cases = {}
            for beta in (1, 2):
                partition = mp.fsum([2*mp.log1p(mp.exp(-beta*e)) for e in energies])
                energy = mp.fsum([2*e*mp.exp(-beta*e)/(1+mp.exp(-beta*e)) for e in energies])
                cases[str(beta)] = [decimal(partition), decimal(energy), decimal(partition+beta*energy)]
            answer[str(q)] = dict(momenta=q**3, zero_momenta=sum(e == 0 for e in energies), cases=cases)
    return answer


def exterior(one):
    modes = len(one)
    occupied = [tuple(i for i in range(modes) if state & (1 << i)) for state in range(1 << modes)]
    return [[determinant(tuple(tuple(one[i][j] for j in source) for i in target))
             if len(target) == len(source) else ZERO for source in occupied] for target in occupied]


def fock(sites):
    modes, dimension = 2*sites, 1 << (2*sites)
    rotation = mat([[F(3, 5), F(4, 5)], [F(-4, 5), F(3, 5)]])
    matrices = []
    for inverse in (False, True):
        for site in range(sites):
            step = eye(modes)
            for i in range(2):
                for j in range(2):
                    step[2*site+i][2*site+j] = rotation[j][i] if inverse else rotation[i][j]
            matrices.append(step)
        flight = [[ZERO]*modes for _ in range(modes)]
        for target in range(modes):
            site, spin = divmod(target, 2)
            source = 2*((site-1+2*spin) % sites)+spin
            flight[target][source] = ONE
        matrices.append(flight)
    phase = eye(modes)
    phase[1][1] = J
    matrices.append(phase)
    one = eye(modes)
    snapshots = []
    for gate in matrices:
        one = mul(gate, one)
        snapshots.append(exterior(one))
    full = snapshots[-1]
    need(mul(dagger(one), one) == eye(modes), "one-particle inverse failed")
    need(mul(dagger(full), full) == eye(dimension), "all-sector Fock inverse failed")
    history = []
    for initial in range(dimension):
        for stage, snapshot in enumerate(snapshots):
            column = [[i, encoded([[row[initial]]])[0][0]]
                      for i, row in enumerate(snapshot) if row[initial]]
            history.append([initial, stage, column])
    return dict(sites=sites, modes=modes, fock_dimension=dimension, local_gates=2*sites+1,
                flight_stages=2, duration="2", basis_stage_evaluations=len(matrices)*dimension,
                one_particle=encoded(one), full_fock_sha256=digest(encoded(full)),
                full_history_sha256=digest(history),
                sector_dimensions=[math.comb(modes, n) for n in range(modes+1)])


def validate_resolution(stages):
    """Check actual paired positive effects and native velocities, then derive M."""
    need(type(stages) is list, "resolution must be a list")
    total = F(0)
    velocity = [[F(0)]*3 for _ in range(3)]
    for row in stages:
        need(type(row) is list and len(row) == 3, "malformed flight resolution")
        w, v, r = row
        need(type(w) is str and type(v) is list and type(r) is list
             and len(v) == len(r) == 3 and all(type(x) is str for x in v+r),
             "resolution requires exact rational strings")
        w, v, r = F(w), [F(x) for x in v], [F(x) for x in r]
        need(w >= 0, "negative elapsed time")
        need(sum(x*x for x in v) <= 1, "native velocity exceeds its cap")
        need(sum(x*x for x in r) <= 1, "velocity effect is not positive")
        total += w
        for i in range(3):
            for j in range(3):
                velocity[i][j] += w*v[i]*r[j]
    need(total <= 1, "uncharged flight time exceeds the whole period")
    return velocity, total


def budget():
    cases = {
        "isotropic": [["1/3", ["1", "0", "0"], ["1", "0", "0"]],
                      ["1/3", ["0", "1", "0"], ["0", "1", "0"]],
                      ["1/3", ["0", "0", "1"], ["0", "0", "1"]]],
        "anisotropic": [["1/2", ["1", "0", "0"], ["1", "0", "0"]],
                        ["1/3", ["0", "1", "0"], ["0", "1", "0"]],
                        ["1/6", ["0", "0", "1"], ["0", "0", "1"]]],
        "waiting": [["1/4", ["1", "0", "0"], ["1", "0", "0"]],
                    ["1/4", ["0", "1", "0"], ["0", "1", "0"]],
                    ["1/4", ["0", "0", "1"], ["0", "0", "1"]]],
        "rotated": [["1/2", ["3/5", "4/5", "0"], ["1", "0", "0"]],
                    ["1/3", ["-4/5", "3/5", "0"], ["0", "1", "0"]],
                    ["1/6", ["0", "0", "1"], ["0", "0", "-1"]]],
    }
    out = {}
    for name, stages in cases.items():
        matrix, charged = validate_resolution(stages)
        out[name] = dict(stages=stages, velocity=[[str(x) for x in row] for row in matrix],
                         charged_flight_time=str(charged), wait_time=str(1-charged))
    return out


def finite_delay():
    rows = []
    for a, b in itertools.product(range(4), repeat=2):
        defects = ["64" if a and b and a != x and b != y else "0"
                   for x in range(1, 4) for y in range(1, 4)]
        rows.append([[a, b], defects])
    total = 256*sum(F(2*a+b, 17)**2 for a in range(1, 4) for b in range(1, 4))
    probabilities = []
    with mp.workdps(75):
        for denominator in (64, 16, 4):
            response = mp.sin(mp.mpf(2)/denominator)
            probabilities.append([f"1/{denominator}", decimal((1+response)/2),
                                  decimal((1-response)/2), decimal(response)])
    return dict(pauli_basis=rows, mixed_interaction_defect=str(total), signals=probabilities)


def current():
    rows = []
    for a, b, d, e in itertools.product(range(-1, 2), repeat=4):
        rho = a*a+b*b+d*d+e*e
        if rho == 0:
            continue
        j = [2*(a*d+b*e), 2*(a*e-b*d), a*a+b*b-d*d-e*e]
        need(sum(x*x for x in j) == rho*rho, "Pauli current flux identity")
        rows.append([[a, b, d, e], str(rho), [str(x) for x in j]])
    return dict(spinors=len(rows), rows=rows)


@lru_cache(None)
def reconstruct():
    from . import tetra_check, stencil_check
    return dict(words={"".join(map(str, p)): native_word(p) for p in itertools.permutations(range(3))},
                spectra={str(q): spectral_catalog(q) for q in (4, 8)}, nodes=nodes(), thermal=thermal(),
                fock={str(q): fock(q) for q in (2, 3)}, budget=budget(),
                finite_delay=finite_delay(), current=current(), tetra=tetra_check.reconstruct(),
                stencil=stencil_check.reconstruct())


def verify_evidence(candidate):
    same(candidate, reconstruct(), "evidence differs from complete independent reconstruction")
