"""Continuous native M6 noise, not an endpoint-only Pauli surrogate."""

import numpy as np
from scipy.linalg import expm

from .algebra import encode


def choi(channel, d):
    answer = np.zeros((d*d, d*d), dtype=complex)
    for i in range(d):
        for j in range(d):
            basis = np.zeros((d, d))
            basis[i, j] = 1
            out = (channel@basis.reshape(-1, order='F')).reshape((d, d), order='F')
            answer[i*d:(i+1)*d, j*d:(j+1)*d] = out/d
    return answer


def examples():
    result = []
    for kind, rate, duration in [('mix', .31, .7), ('pair_phase', .23, .41)]:
        d = 6
        h = np.zeros((d, d), dtype=complex)
        if kind == 'mix':
            h[0, 1], h[1, 0] = -.7j, .7j
        else:
            h[:3, :3] = np.pi*np.eye(3)
        dephase = np.diag([float(i == j) for j in range(d) for i in range(d)])
        generator = -1j*(np.kron(np.eye(d), h)-np.kron(h.T, np.eye(d)))
        channel = expm(duration*(generator+rate*(dephase-np.eye(d*d))))
        unitary = expm(-1j*duration*h)
        ideal = np.kron(unitary.conj(), unitary)
        zero = np.exp(-rate*duration)
        remainder = (channel-zero*ideal)/(1-zero)
        result.append(dict(kind=kind, rate=rate, duration=duration, zero_jump_probability=float(zero),
                           channel_choi=encode(choi(channel, d)), fault_choi=encode(choi(remainder, d))))
    return result
