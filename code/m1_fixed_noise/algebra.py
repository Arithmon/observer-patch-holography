"""Complete Steane decoder, native interfaces and exact memory controls."""

import math

import numpy as np

from m1_fermionic_source.circuits import native_one, native_cz
from .recovery import CHECKS, correction, parity, syndrome


def decoder():
    words = [0]
    for h in CHECKS:
        words += [w ^ h for w in words]
    rows = []
    for s in range(64):
        x, z = correction(s)
        branches = []
        for logical in range(2):
            branches.append(sorted([[w ^ (127*logical) ^ x,
                                     (-1)**parity(z & (w ^ (127*logical)))] for w in words]))
        rows.append(dict(syndrome=s, correction=[x, z], kraus_rows=branches, denominator_squared=8))
    return rows


def interfaces():
    h, t = native_one('h')[:2, :2], np.exp(1j*math.pi/8)*native_one('rz', -math.pi/8)[:2, :2]
    target_h = np.kron(h, np.eye(2))
    cx = target_h@native_cz()[:4, :4]@target_h
    # Native CZ differs by an immaterial scalar; keep the actual matrix.
    return dict(h=encode(h), t=encode(t), cx=encode(cx))


def encode(a):
    a = np.asarray(a, dtype=complex)
    return np.stack((a.real, a.imag), axis=-1).tolist()


def memory():
    bad, unchanged = [0]*8, [0]*8
    rows = decoder()
    code = np.zeros((128, 2))
    for j, entries in enumerate(rows[0]['kraus_rows']):
        for i, sign in entries:
            code[i, j] = sign/math.sqrt(8)
    recovery = []
    for row in rows:
        k = np.zeros((2, 128))
        for j, entries in enumerate(row['kraus_rows']):
            for i, sign in entries:
                k[j, i] = sign/math.sqrt(8)
        recovery.append(k)
    branch_checks = []
    for z in range(128):
        phase = np.array([(-1)**parity(i & z) for i in range(128)])
        out = phase[:, None]*code
        s = syndrome(0, z)
        logical = recovery[s]@out
        is_bad = abs(np.trace(logical))/2 < .5
        bad[z.bit_count()] += int(is_bad)
        unchanged[z.bit_count()] += int(abs(np.trace(code.T@out))/2 > .5)
        if z == 0 or z & (z-1) == 0:
            branch_checks.append(dict(z=z, syndrome=s, logical=encode(logical)))
    return dict(logical_error_by_weight=bad, raw_entangled_fidelity_by_weight=unchanged,
                single_phase_branches=branch_checks)


def accounting():
    operators = []
    for row in decoder():
        k = np.zeros((2, 128), dtype=complex)
        for j, entries in enumerate(row['kraus_rows']):
            for i, sign in entries:
                k[j, i] = sign/math.sqrt(8)
        operators.append(k)
    code = operators[0].conj().T
    effect = np.array([[.35, .12+.05j], [.12-.05j, .75]])
    physical = sum(k.conj().T@effect@k for k in operators)
    psi = np.array([math.sqrt(.6), 1j*math.sqrt(.4)])
    cases = []
    for phase in (0, 1, 2, 4):
        state = np.array([(-1)**parity(i & phase) for i in range(128)])*(code@psi)
        cases.append(dict(phase=phase, decoded_expectation=float(np.vdot(state, physical@state).real),
                          raw_code_penalty=float(1-np.linalg.norm(code.conj().T@state)**2)))
    return dict(logical_effect=encode(effect), logical_state=encode(psi), cases=cases,
                abort_value=1., before_encoding_phase_error=.125,
                raw_input_entangled_fidelity=.875, after_encoding_corrected_fidelity=1.)
