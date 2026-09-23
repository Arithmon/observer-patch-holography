"""Candidate evidence from exact matrix products and executed fermionic gates."""

from functools import lru_cache
import hashlib
import itertools
import json

import mpmath as mp
import sympy as s


I = s.I
ROOT = s.sqrt(2)
PAULI = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -I], [I, 0]]),
         s.Matrix([[1, 0], [0, -1]]))
IDENTITY = s.eye(2)
PERMUTATIONS = tuple(itertools.permutations(range(3)))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def pair(value):
    value = s.expand(value)
    return [str(s.re(value)), str(s.im(value))]


def matrix(value):
    return [[pair(x) for x in row] for row in value.tolist()]


def quadratic(value):
    value = s.expand(value)
    b = value.coeff(ROOT)
    a = s.expand(value-b*ROOT)
    if not (a.is_Rational and b.is_Rational):
        raise ValueError("coefficient outside the declared exact quadratic field")
    return [str(a), str(b)]


def native_word(order):
    kernel = {(0, 0, 0): IDENTITY}
    history = []
    for stage, axis in enumerate(order):
        nxt = {}
        for source, amplitude in sorted(kernel.items()):
            for sign in (-1, 1):
                target = list(source)
                target[axis] += sign
                target = tuple(target)
                projector = (IDENTITY+sign*PAULI[axis])/2
                outgoing = projector*amplitude
                nxt[target] = nxt.get(target, s.zeros(2))+outgoing
                positions = []
                for fraction in (s.Rational(0), s.Rational(1, 4), s.Rational(1, 2),
                                 s.Rational(3, 4), s.Rational(1)):
                    pos = [s.Rational(x) for x in source]
                    pos[axis] += sign*fraction
                    positions.append([str(stage+fraction), [str(x) for x in pos]])
                history.append(dict(stage=stage+1, source=list(source), target=list(target),
                                    sign=sign, amplitude=matrix(outgoing), positions=positions))
        kernel = nxt
    rows = [[list(k), matrix(v)] for k, v in sorted(kernel.items())]
    return dict(order=list(order), duration="3", speed_cap="1", effective_speed="1/3",
                intermediate_events=len(history), history_sha256=digest(history),
                kernel=rows)


@lru_cache(None)
def spectral_catalog(q):
    if q not in (4, 8):
        raise ValueError("complete exact catalog is q=4 or q=8")
    cosines = ([1, 0, -1, 0] if q == 4 else
               [1, ROOT/2, 0, -ROOT/2, -1, -ROOT/2, 0, ROOT/2])
    sinuses = [cosines[(i-q//4) % q] for i in range(q)]
    steps = [[cosines[k]*IDENTITY-I*sinuses[k]*PAULI[axis] for k in range(q)]
             for axis in range(3)]
    rows = []
    zeros, pi_modes = [], []
    for k in itertools.product(range(q), repeat=3):
        u = steps[2][k[2]]*steps[1][k[1]]*steps[0][k[0]]
        coefficients = [s.expand(s.trace(u)/2)]
        coefficients += [s.expand(I*s.trace(p*u)/2) for p in PAULI]
        encoded = [quadratic(x) for x in coefficients]
        if coefficients[0] == 1:
            zeros.append(list(k))
        if coefficients[0] == -1:
            pi_modes.append(list(k))
        rows.append([list(k), encoded])
    return dict(q=q, momenta=q**3, modes=2*q**3, rows=rows,
                zero_points=zeros, pi_points=pi_modes)


def nodal_derivatives():
    rows = []
    angles = (s.Integer(0), s.pi/2, s.pi, 3*s.pi/2)
    for index in itertools.product(range(4), repeat=3):
        steps = [s.cos(angles[index[a]])*IDENTITY-I*s.sin(angles[index[a]])*PAULI[a]
                 for a in range(3)]
        u = steps[2]*steps[1]*steps[0]
        if u != IDENTITY:
            continue
        velocity = []
        for a in range(3):
            derivative = -s.sin(angles[index[a]])*IDENTITY-I*s.cos(angles[index[a]])*PAULI[a]
            word = list(steps)
            word[a] = derivative
            du = word[2]*word[1]*word[0]
            velocity.append([str(s.expand(s.trace(p*I*du)/2)) for p in PAULI])
        rows.append(dict(index=list(index), derivative=velocity,
                         chirality=int(s.Matrix(velocity).det())))
    return rows


def number(value):
    return mp.nstr(value, 21, strip_zeros=False)


def trig(q, k):
    if (4*k) % q == 0:
        j = (4*k//q) % 4
        return mp.mpf((1, 0, -1, 0)[j]), mp.mpf((0, 1, 0, -1)[j])
    angle = 2*mp.pi*k/q
    return mp.cos(angle), mp.sin(angle)


def thermal_catalog():
    out = {}
    with mp.workdps(65):
        pauli = (mp.matrix([[0, 1], [1, 0]]), mp.matrix([[0, -1j], [1j, 0]]),
                 mp.matrix([[1, 0], [0, -1]]))
        for q in (4, 8, 12, 16):
            steps = [[trig(q, k)[0]*mp.eye(2)-mp.j*trig(q, k)[1]*pauli[a]
                      for k in range(q)] for a in range(3)]
            energies = []
            zeros = 0
            for k in itertools.product(range(q), repeat=3):
                u = steps[2][k[2]]*steps[1][k[1]]*steps[0][k[0]]
                half_trace = mp.re(u[0, 0]+u[1, 1])/2
                if not -1 <= half_trace <= 1:
                    raise ValueError("unitary eigenphase outside its domain")
                energy = mp.acos(half_trace)*q/3
                energies.append(energy)
                zeros += energy == 0
            cases = {}
            for beta in (1, 2):
                logz = 2*mp.fsum(mp.log1p(mp.exp(-beta*e)) for e in energies)
                energy = 2*mp.fsum(e/(mp.exp(beta*e)+1) for e in energies)
                cases[str(beta)] = [number(logz), number(energy), number(logz+beta*energy)]
            out[str(q)] = dict(momenta=q**3, zero_momenta=zeros, cases=cases)
    return out


def basis_gate(state, site, coin):
    a, b = 2*site, 2*site+1
    occupancy = ((state >> a) & 1)+2*((state >> b) & 1)
    if occupancy == 0:
        return [(state, s.Integer(1))]
    if occupancy == 3:
        return [(state, s.expand(coin.det()))]
    source = 0 if occupancy == 1 else 1
    cleared = state & ~(1 << a) & ~(1 << b)
    return [(cleared | (1 << (a+target)), coin[target, source]) for target in range(2)
            if coin[target, source] != 0]


def flight_basis(state, sites, retain_sign=True):
    occupied = [i for i in range(2*sites) if state & (1 << i)]
    moved = [2*((i//2+(1 if i % 2 == 0 else -1)) % sites)+i % 2 for i in occupied]
    inversions = sum(moved[i] > moved[j] for i in range(len(moved)) for j in range(i+1, len(moved)))
    return sum(1 << i for i in moved), (-1)**inversions if retain_sign else 1


def fock_program(sites, retain_sign=True):
    r = s.Matrix([[s.Rational(3, 5), s.Rational(4, 5)],
                  [-s.Rational(4, 5), s.Rational(3, 5)]])
    phase = s.diag(1, I)
    instructions = [("coin", site, r) for site in range(sites)]
    instructions += [("flight", None, None)]
    instructions += [("coin", site, r.T) for site in range(sites)]
    instructions += [("flight", None, None), ("coin", 0, phase)]
    dimension = 1 << (2*sites)
    columns = []
    traces = []
    for initial in range(dimension):
        state = {initial: s.Integer(1)}
        for stage, (kind, site, coin) in enumerate(instructions):
            out = {}
            for basis, amplitude in state.items():
                terms = ([flight_basis(basis, sites, retain_sign)] if kind == "flight"
                         else basis_gate(basis, site, coin))
                for target, coefficient in terms:
                    out[target] = s.expand(out.get(target, 0)+coefficient*amplitude)
            state = {i: x for i, x in out.items() if x != 0}
            traces.append([initial, stage, [[i, pair(x)] for i, x in sorted(state.items())]])
        columns.append([state.get(i, s.Integer(0)) for i in range(dimension)])
    unitary = s.Matrix.hstack(*(s.Matrix(x) for x in columns))
    single = [1 << i for i in range(2*sites)]
    one = unitary.extract(single, single)
    return dict(sites=sites, modes=2*sites, fock_dimension=dimension,
                local_gates=2*sites+1, flight_stages=2, duration="2",
                basis_stage_evaluations=len(instructions)*dimension,
                one_particle=matrix(one), full_fock_sha256=digest(matrix(unitary)),
                full_history_sha256=digest(traces),
                sector_dimensions=[s.binomial(2*sites, n).__int__() for n in range(2*sites+1)])


def budget_examples():
    # Each pair has effects (I +/- r.sigma)/2 and is charged its time fraction.
    cases = {
        "isotropic": [(s.Rational(1, 3), (1, 0, 0), (1, 0, 0)),
                      (s.Rational(1, 3), (0, 1, 0), (0, 1, 0)),
                      (s.Rational(1, 3), (0, 0, 1), (0, 0, 1))],
        "anisotropic": [(s.Rational(1, 2), (1, 0, 0), (1, 0, 0)),
                        (s.Rational(1, 3), (0, 1, 0), (0, 1, 0)),
                        (s.Rational(1, 6), (0, 0, 1), (0, 0, 1))],
        "waiting": [(s.Rational(1, 4), (1, 0, 0), (1, 0, 0)),
                    (s.Rational(1, 4), (0, 1, 0), (0, 1, 0)),
                    (s.Rational(1, 4), (0, 0, 1), (0, 0, 1))],
        "rotated": [(s.Rational(1, 2), (s.Rational(3, 5), s.Rational(4, 5), 0), (1, 0, 0)),
                    (s.Rational(1, 3), (-s.Rational(4, 5), s.Rational(3, 5), 0), (0, 1, 0)),
                    (s.Rational(1, 6), (0, 0, 1), (0, 0, -1))],
    }
    out = {}
    for name, stages in cases.items():
        velocity = s.zeros(3)
        for weight, spatial, spin in stages:
            velocity += weight*s.Matrix(spatial)*s.Matrix(spin).T
        out[name] = dict(stages=[[str(w), [str(x) for x in u], [str(x) for x in r]] for w, u, r in stages],
                         velocity=[[str(x) for x in row] for row in velocity.tolist()],
                         charged_flight_time=str(sum(w for w, _, _ in stages)),
                         wait_time=str(1-sum(w for w, _, _ in stages)))
    return out


def finite_delay():
    local = [IDENTITY, *PAULI]
    probes = [(s.kronecker_product(a, IDENTITY), s.kronecker_product(IDENTITY, b))
              for a in PAULI for b in PAULI]
    rows = []
    interaction = s.zeros(4)
    for a, b in itertools.product(range(4), repeat=2):
        h = s.kronecker_product(local[a], local[b])
        defects = []
        for x, y in probes:
            first = h*x-x*h
            defect = first*y-y*first
            defects.append(str(s.trace(defect.conjugate().T*defect)))
        rows.append([[a, b], defects])
        if a and b:
            interaction += s.Rational(2*a+b, 17)*h
    total = 0
    for x, y in probes:
        first = interaction*x-x*interaction
        defect = first*y-y*first
        total += s.trace(defect.conjugate().T*defect)
    probabilities = []
    with mp.workdps(65):
        h = mp.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
        initial = [mp.matrix([1, 1, mp.j, mp.j])/2,
                   mp.matrix([1, -1, mp.j, -mp.j])/2]
        for denominator in (64, 16, 4):
            t = mp.mpf(1)/denominator
            u = mp.cos(t)*mp.eye(4)-mp.j*mp.sin(t)*h
            states = [u*x for x in initial]
            p = [abs(x[0])**2+abs(x[1])**2 for x in states]
            probabilities.append([f"1/{denominator}", *[number(x) for x in p], number(p[0]-p[1])])
    return dict(pauli_basis=rows, mixed_interaction_defect=str(s.expand(total)), signals=probabilities)


def coherent_current():
    rows = []
    for raw in itertools.product((-1, 0, 1), repeat=4):
        if raw == (0, 0, 0, 0):
            continue
        psi = s.Matrix([raw[0]+I*raw[1], raw[2]+I*raw[3]])
        rho = s.expand((psi.conjugate().T*psi)[0])
        current = [s.expand((psi.conjugate().T*p*psi)[0]) for p in PAULI]
        rows.append([list(raw), str(rho), [str(x) for x in current]])
    return dict(spinors=len(rows), rows=rows)


def candidate():
    from . import tetra_model, stencil_model
    return dict(words={"".join(map(str, p)): native_word(p) for p in PERMUTATIONS},
                spectra={str(q): spectral_catalog(q) for q in (4, 8)},
                nodes=nodal_derivatives(), thermal=thermal_catalog(),
                fock={str(q): fock_program(q) for q in (2, 3)},
                budget=budget_examples(), finite_delay=finite_delay(), current=coherent_current(),
                tetra=tetra_model.candidate(), stencil=stencil_model.candidate())
