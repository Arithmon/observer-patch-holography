"""Actual two-bank fermionic clock instruction catalogue and executions."""

import itertools
import math

import numpy as np

from m1_source_realization.model import decompose

from .model import encode
from .spatial import graph


COIN = np.array([[0, 1, 1, 1], [1, 0, 1j, -1j], [1, -1j, 0, 1j], [1, 1j, -1j, 0]])/np.sqrt(3)


def layers(spacing, masses=(.7, 1.1), speed=3.):
    diagonal, moves = decompose(COIN)
    def coin(offset):
        return [('phase', [offset+j], np.angle(z)) for j, z in enumerate(diagonal)]+[
            ('mix', [offset+i, offset+j], g) for i, j, g in moves]
    tape = coin(0)+coin(8)
    tape += [('ferry', [j], None) for j in range(16)]
    tape += [('restore', [j], None) for j in range(16)]
    tape += coin(4)+coin(12)
    for mass_index, mass in enumerate(masses):
        mu = mass*math.sqrt(3)*spacing/speed
        g = np.array([[math.cos(mu), -1j*math.sin(mu)], [-1j*math.sin(mu), math.cos(mu)]])
        tape += [('mix', [8*mass_index+j, 8*mass_index+j+4], g) for j in range(4)]
    return tape


def execute(side, spacing, columns=None):
    vertices, _, flights = graph(side)
    active = [i for i, (_, bank, _) in enumerate(vertices) if bank == 0]
    size = len(active)
    seed = np.eye(size, dtype=complex) if columns is None else np.asarray(columns, complex)
    if seed.ndim == 1:
        seed = seed[:, None]
    state = np.zeros((len(vertices), seed.shape[1]), complex)
    state[active] = seed
    ids = {v: i for i, v in enumerate(vertices)}
    for kind, modes, value in layers(spacing):
        if kind == 'ferry':
            for i, j in flights:
                if vertices[i][2] == modes[0]:
                    state[[i, j]] = state[[j, i]]
        else:
            for cell in itertools.product(range(side), repeat=3):
                if kind == 'restore':
                    pair = [ids[(cell, bank, modes[0])] for bank in (0, 1)]
                    state[pair] = state[pair[::-1]]
                elif kind == 'phase':
                    state[ids[(cell, 0, modes[0])]] *= np.exp(1j*value)
                else:
                    pair = [ids[(cell, 0, j)] for j in modes]
                    state[pair] = value@state[pair]
    return state[active], state[[i for i in range(len(vertices)) if i not in set(active)]]


def evidence():
    rows = []
    for side, spacing in ((1, .07), (2, .13), (3, .03)):
        n = 16*side**3
        probe = np.exp(1j*np.arange(n)**2/17)/np.sqrt(n)
        output, banks = execute(side, spacing, probe)
        tape = layers(spacing)
        counts = {kind: sum(k == kind for k, _, _ in tape) for kind in ('phase', 'mix', 'ferry', 'restore')}
        rows.append(dict(side=side, spacing=spacing, layers=counts,
                         pauli_rotations=counts['phase']+6*counts['mix']+8*(counts['ferry']+counts['restore']),
                         program=[dict(kind=k, modes=m, parameter=(None if v is None else
                                  float(v) if k == 'phase' else encode(v))) for k, m, v in tape],
                         output=encode(output[:, 0]), banks=encode(banks[:, 0])))
    return rows
