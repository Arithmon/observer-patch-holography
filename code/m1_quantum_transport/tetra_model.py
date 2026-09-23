"""Executed minimal tetrahedral optimizer, with its complete four-channel spectrum."""

import itertools
import json
import hashlib

import mpmath as mp
import sympy as s


I, R = s.I, s.sqrt(3)
S = s.Matrix([[0, 1, 1, 1], [1, 0, I, -I], [1, -I, 0, I], [1, I, -I, 0]])
COIN = S/R
POINTS = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def entry(x):
    x = s.expand(x)
    re, im = s.expand(s.re(x)), s.expand(s.im(x))
    result = [s.expand(re-re.coeff(R)*R), re.coeff(R),
              s.expand(im-im.coeff(R)*R), im.coeff(R)]
    if not all(v.is_Rational for v in result):
        raise ValueError("tetrahedral coefficient escaped its exact field")
    return [str(v) for v in result]


def encode(a):
    return [[entry(x) for x in row] for row in a.tolist()]


def bands():
    out = {}
    for sign in (1, -1):
        p = (s.eye(4)+sign*COIN)/2
        velocity = [R*p*s.diag(*(point[i] for point in POINTS))*p for i in range(3)]
        effects = [p*s.diag(*(int(i == j) for i in range(4)))*p for j in range(4)]
        out[str(sign)] = dict(projector=encode(p), rank=int(p.rank()), chirality=sign,
                              pauli=[encode(x) for x in velocity], effects=[encode(x) for x in effects])
    return out


def spectrum():
    rows, zeros, pi_points, middle = [], [], [], []
    for k in itertools.product(range(4), repeat=3):
        phases = [(-I)**sum(x*y for x, y in zip(point, k)) for point in POINTS]
        u = s.diag(*phases)*COIN
        coefficients = [entry(x) for x in u.charpoly().all_coeffs()]
        rows.append([list(k), encode(u), coefficients])
        if u.charpoly().eval(1) == 0:
            zeros.append(list(k))
        if u.charpoly().eval(-1) == 0:
            pi_points.append(list(k))
        if s.expand(u.charpoly().eval(I)) == 0:
            middle.append(list(k))
    return dict(q=4, momenta=64, modes=256, rows_sha256=digest(rows),
                zero_points=zeros, pi_points=pi_points, middle_points=middle)


def fock():
    occupied = [tuple(i for i in range(4) if state & (1 << i)) for state in range(16)]
    rows = [[COIN.extract(target, source).det() if len(target) == len(source) else s.Integer(0)
             for source in occupied] for target in occupied]
    full = s.Matrix(rows)
    return dict(dimension=16, sector_dimensions=[1, 4, 6, 4, 1], matrix_sha256=digest(encode(full)))


def thermal():
    out = {}
    with mp.workdps(65):
        c = mp.matrix([[0, 1, 1, 1], [1, 0, mp.j, -mp.j],
                       [1, -mp.j, 0, mp.j], [1, mp.j, -mp.j, 0]])/mp.sqrt(3)
        for q in (4, 8, 12, 16):
            phases = [(-mp.j)**(4*k//q) if 4*k % q == 0 else mp.exp(-2*mp.pi*mp.j*k/q)
                      for k in range(q)]
            energies, zeros = [], 0
            for k in itertools.product(range(q), repeat=3):
                u = mp.diag([phases[sum(x*y for x, y in zip(point, k)) % q] for point in POINTS])*c
                # Exact endpoint classifications avoid rounding an acos argument to its boundary.
                if all(2*x % q == 0 for x in k):
                    epsilon = mp.mpf(0)
                    zeros += 1
                elif all(4*x % q == 0 and 2*x % q != 0 for x in k):
                    epsilon = mp.pi/2
                else:
                    squared = u*u
                    arg = mp.re(sum(squared[i, i] for i in range(4)))/4
                    if not -1 < arg < 1:
                        raise ValueError("tetrahedral eigenphase escaped its interior")
                    epsilon = mp.acos(arg)/2
                energies.extend([epsilon*q/mp.sqrt(3), (mp.pi-epsilon)*q/mp.sqrt(3)])
            cases = {}
            for beta in (1, 2):
                logz = mp.fsum(2*mp.log1p(mp.exp(-beta*e)) for e in energies)
                energy = mp.fsum(2*e/(1+mp.exp(beta*e)) for e in energies)
                cases[str(beta)] = [mp.nstr(x, 21, strip_zeros=False) for x in (logz, energy, logz+beta*energy)]
            out[str(q)] = dict(momenta=q**3, modes=4*q**3, zero_momenta=zeros, cases=cases)
    return out


def candidate():
    trajectories = [[list(point), str(f), [str(f*x) for x in point], str(R*f)]
                    for point in POINTS for f in (s.Rational(0), s.Rational(1, 4), s.Rational(1, 2),
                                                 s.Rational(3, 4), s.Rational(1))]
    signals = []
    for channel, point in enumerate(POINTS):
        prepared = COIN[:, channel]
        output = (COIN*prepared).applyfunc(s.simplify)
        signals.append(dict(channel=channel, prepared=encode(prepared), output=encode(output),
                            probabilities=[str(s.expand(s.conjugate(x)*x)) for x in output],
                            empty_probabilities=["0"]*4, destination=list(point), time="sqrt(3)"))
    return dict(points=[list(x) for x in POINTS], coin=encode(COIN), bands=bands(), spectrum=spectrum(),
                flight_duration="sqrt(3)", native_speed="1", field_speed="1/3",
                trajectories_sha256=digest(trajectories), signals=signals, fock=fock(), thermal=thermal())
