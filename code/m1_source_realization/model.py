"""Producer for source controls, code channels and explicit gate executions."""

import itertools
import json
import math
from numbers import Real
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as sp


ROOT = Path(__file__).resolve().parents[2]


def source_generators():
    packet = json.loads((ROOT/'code/source_selection_model/response.json').read_text(encoding='utf-8'))
    out = []
    for blocks in packet['generators']:
        matrix = sp.zeros(6)
        for block, rows in enumerate(blocks):
            for i, row in enumerate(rows):
                for j, value in enumerate(row):
                    a, b, c, d = map(sp.Rational, value)
                    matrix[3*block+i, 3*block+j] = a+b*sp.sqrt(5)+sp.I*(c+d*sp.sqrt(5))
        out.append(matrix)
    return out


def field(x):
    x = sp.expand(sp.simplify(x))
    return [str(x.coeff(sp.sqrt(5), 0)), str(x.coeff(sp.sqrt(5), 1))]


def controls():
    generators = source_generators()
    gram = sp.Matrix(12, 12, lambda i, j: sp.simplify(-sp.trace(generators[i]*generators[j])))
    inverse = gram.inv().applyfunc(sp.simplify)
    targets = []
    for i in range(3):
        h = sp.zeros(6)
        h[i, i] = sp.I
        targets.append((f'phase_{i}', h))
    for i, j in itertools.combinations(range(3), 2):
        for kind in ('real', 'imaginary'):
            h = sp.zeros(6)
            h[i, j], h[j, i] = (1, -1) if kind == 'real' else (sp.I, sp.I)
            targets.append((f'{kind}_{i}{j}', h))
    centre = sp.diag(sp.I, sp.I, sp.I, 0, 0, 0)
    targets.append(('relative_phase', centre))
    rows = []
    for name, h in targets:
        rhs = sp.Matrix([-sp.trace(d*h) for d in generators])
        coefficients = (inverse*rhs).applyfunc(sp.simplify)
        if sum((x*d for x, d in zip(coefficients, generators)), sp.zeros(6)).applyfunc(sp.simplify) != h:
            raise ValueError('source control reconstruction')
        rows.append(dict(name=name, coefficients=[field(x) for x in coefficients]))
    return dict(gram_determinant=field(gram.det()), rows=rows)


def choi(kraus):
    d = len(kraus[0])
    # Input first, output second: vectorization by columns.
    return sum(np.outer(k.T.reshape(-1), k.T.reshape(-1).conj()) for k in kraus)/d


def fourier(d):
    return np.exp(2j*np.pi*np.outer(np.arange(d), np.arange(d))/d)/np.sqrt(d)


def weyl(d, shift, phase):
    out = np.zeros((d, d), complex)
    for j in range(d):
        out[(j+shift) % d, j] = np.exp(2j*np.pi*phase*j/d)
    return out


def classical_fidelity(kraus, basis):
    return sum(sum(abs(np.vdot(v, k@v))**2 for k in kraus) for v in basis.T)/len(basis)


def selected_channel(d, ez, ex):
    if type(d) is not int or d < 2:
        raise ValueError('channel dimension must be an integer at least two')
    if any(isinstance(e, (bool, np.bool_)) or not isinstance(e, Real) or not math.isfinite(e)
           or not 0 <= e <= 1 for e in (ez, ex)):
        raise ValueError('agreement errors must be finite probabilities')
    # Once a fidelity bound is weaker than uniform guessing, entropy is
    # maximized at uniformity rather than by saturating the loose bound.
    ez, ex = min(ez, (d-1)/d), min(ex, (d-1)/d)
    p = np.array([1-ez]+[ez/(d-1)]*(d-1))
    q = np.array([1-ex]+[ex/(d-1)]*(d-1))
    return [np.sqrt(p[a]*q[b])*weyl(d, a, b) for a in range(d) for b in range(d)]


def channels():
    rows = []
    for d in (2, 3, 4):
        for ez, ex in ((0., 0.), (0., (d-1)/d), (.02, .03), (0., .05)):
            operators = selected_channel(d, ez, ex)
            j = choi(operators)
            eigenvalues = np.maximum(np.linalg.eigvalsh(j), 0)
            omega = np.eye(d).reshape(-1)/np.sqrt(d)
            delta = j-np.outer(omega, omega)
            rows.append(dict(d=d, ez=ez, ex=ex,
                             choi_eigenvalues=sorted(eigenvalues.tolist()),
                             z_fidelity=float(classical_fidelity(operators, np.eye(d))),
                             x_fidelity=float(classical_fidelity(operators, fourier(d))),
                             entanglement_fidelity=float(np.vdot(omega, j@omega).real),
                             entropy=float(-sum(x*np.log(x) for x in eigenvalues if x > 1e-14)),
                             choi_half_trace_distance=float(sum(abs(np.linalg.eigvalsh(delta)))/2)))
    return rows


def code_transfer_kraus(d, keep_phase=True):
    """Move a proper code of a six-level source into a fresh output.

    The complement is an explicit failure flag. It is neither postselected
    away nor represented as a reversible map of the whole six-level source.
    """
    if type(d) is not int or d not in (2, 4) or type(keep_phase) is not bool:
        raise ValueError('the source interface admits only qubit and pair proper codes')
    out = []
    if keep_phase:
        k = np.zeros((d+1, 6), complex)
        k[:d, :d] = np.eye(d)
        out.append(k)
    else:
        for i in range(d):
            k = np.zeros((d+1, 6), complex)
            k[i, i] = 1
            out.append(k)
    for i in range(d, 6):
        k = np.zeros((d+1, 6), complex)
        k[d, i] = 1
        out.append(k)
    return out


def act(kraus, rho):
    return sum(k@rho@k.conj().T for k in kraus)


def instruments():
    rows = []
    for d in (2, 4):
        for coherent in (False, True):
            ks = code_transfer_kraus(d, coherent)
            images = []
            for i, j in itertools.product(range(6), repeat=2):
                e = np.zeros((6, 6), complex)
                e[i, j] = 1
                z = act(ks, e)
                images.append([[float(x.real), float(x.imag)] for x in z.reshape(-1)])
            rows.append(dict(d=d, coherent=coherent,
                             completeness_error=float(np.linalg.norm(sum(k.conj().T@k for k in ks)-np.eye(6))),
                             matrix_unit_images=images))
    return rows


def entangler():
    h = np.diag([1., 1., 1., 0., 0., 0.])
    rows = []
    for theta in (0., math.pi/7, math.pi/2, math.pi):
        pulse = expm(-1j*theta*h)
        code = np.exp(1j*theta)*pulse[:4, :4]
        plus = np.ones(4)/2
        psi = (code@plus).reshape(2, 2)
        reduced = psi@psi.conj().T
        rows.append(dict(theta=theta,
                         code_diagonal=[[float(z.real), float(z.imag)] for z in np.diag(code)],
                         concurrence=float(2*abs(np.linalg.det(psi))),
                         marginal_purity=float(np.trace(reduced@reduced).real)))
    return rows


def decompose(unitary):
    """Complete adjacent-row QR; no entry or final diagonal is dropped."""
    a = unitary.copy()
    moves = []
    size = len(a)
    for col in range(size-1):
        for lower in range(size-1, col, -1):
            upper = lower-1
            x, y = a[upper, col], a[lower, col]
            radius = np.hypot(abs(x), abs(y))
            g = np.eye(2, dtype=complex) if radius == 0 else np.array([[x.conjugate(), y.conjugate()], [-y, x]])/radius
            a[[upper, lower], :] = g@a[[upper, lower], :]
            moves.append((upper, lower, g.conj().T))
    diagonal = np.diag(a).copy()
    return diagonal, list(reversed(moves))


def gate_catalog():
    c = np.array([[0, 1, 1, 1], [1, 0, 1j, -1j], [1, -1j, 0, 1j], [1, 1j, -1j, 0]])/np.sqrt(3)
    diagonal, moves = decompose(c)
    actual = np.diag(diagonal)
    stages = []
    for upper, lower, g in moves:
        # An SU(2) rotation of the one-occupation code, realized inside U(3).
        native = np.eye(6, dtype=complex)
        native[1:3, 1:3] = g
        stages.append(dict(pair=[upper, lower], matrix=[[[float(z.real), float(z.imag)] for z in row] for row in g],
                           native_unitarity_error=float(np.linalg.norm(native.conj().T@native-np.eye(6)))))
        actual[[upper, lower], :] = g@actual[[upper, lower], :]
    return dict(diagonal=[[float(z.real), float(z.imag)] for z in diagonal], stages=stages,
                reconstruction_error=float(np.linalg.norm(actual-c)))


def candidate():
    return dict(controls=controls(), channels=channels(), instruments=instruments(),
                entangler=entangler(), coin=gate_catalog())
