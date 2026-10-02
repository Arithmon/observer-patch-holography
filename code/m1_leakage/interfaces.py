"""Producer: continue actual source transfers into a fresh code blank."""

import numpy as np
from m1_source_realization.model import code_transfer_kraus
from m1_fermionic_source.circuits import native_cz, native_one
from .format import pack, need


def reduction(d):
    need(type(d) is int and d in (2, 4), 'native code dimension')
    transfer = code_transfer_kraus(d)
    # A transfer has a code plus an orthogonal failure symbol. These two
    # continuation branches do not coherently add code and failure amplitudes.
    ok = np.eye(d, d+1, dtype=complex)
    failure = np.zeros((d, d+1), complex)
    failure[0, d] = 1
    groups = [[ok@transfer[0]], [failure@k for k in transfer[1:]]]
    return groups


def gate(name):
    if name == 'cz':
        return native_cz()[:4, :4]
    if name == 'h_t':
        h = native_one('h')[:2, :2]
        t = np.exp(1j*np.pi/8)*native_one('rz', -np.pi/8)[:2, :2]
        return np.kron(h, t)
    raise ValueError('frozen native gate')


def candidate():
    reductions = [dict(d=d, groups=[[pack(k) for k in group] for group in reduction(d)])
                  for d in (2, 4)]
    pair = []
    groups = reduction(2)
    for name in ('cz', 'h_t'):
        u = gate(name)
        branches = []
        for first in range(2):
            for second in range(2):
                ks = [u@np.kron(a, b) for a in groups[first] for b in groups[second]]
                branches.append(dict(flags=[first, second], maps=[pack(k) for k in ks]))
        pair.append(dict(name=name, branches=branches))
    return dict(reductions=reductions, pair_gates=pair)
