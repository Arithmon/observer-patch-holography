"""Predict disturbed field histories using inferred coefficients only.

The oscillator basis, tensor-factor interpretation, field family and detector
are supplied. No backend graph, parameter, source matrix or coordinate enters
this module. Sparse matrix evolution is observer-side numerical inference,
not a source-compiled implementation of that inference algorithm.
"""
from __future__ import annotations

from math import prod
from numbers import Integral, Real

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..relativistic_scalar.lattice import local_operators


def _array(value, shape, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot be masked")
    result = np.asarray(value)
    if (result.shape != shape or result.dtype.kind not in 'fiu'
            or not np.isfinite(result).all()):
        raise ValueError(f"{name} must have shape {shape} and finite real entries")
    return result.astype(float)


def _real(value, name):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not np.isfinite(value)):
        raise ValueError(f"{name} must be a finite real")
    return float(value)


class RecoveredField:
    """Finite local scalar action on opaque labels; topology comes from records."""

    def __init__(self, dimensions, quadratic, quartic, bonds, *, omega=2.):
        if (not isinstance(dimensions, (tuple, list)) or not dimensions
                or any(isinstance(d, (bool, np.bool_)) or not isinstance(d, Integral)
                       or d < 2 for d in dimensions)):
            raise ValueError('occupation dimensions must be integers at least two')
        self.dimensions = tuple(int(d) for d in dimensions)
        n = len(dimensions)
        self.dimension = prod(dimensions)
        self.quadratic = _array(quadratic, (n,), 'onsite quadratic coefficient')
        self.quartic = _array(quartic, (n,), 'onsite quartic coefficient')
        self.bonds = _array(bonds, (n, n), 'pair coefficients')
        if (not np.allclose(self.bonds, self.bonds.T, rtol=0, atol=1e-13)
                or np.any(np.diag(self.bonds) != 0)):
            raise ValueError('pair coefficients must be symmetric with zero diagonal')
        self.operators = tuple(local_operators(d, omega) for d in dimensions)
        h = sparse.csr_matrix((self.dimension, self.dimension), dtype=float)
        for site, op in enumerate(self.operators):
            local = op['p2']/2+self.quadratic[site]*op['q2']/2+self.quartic[site]*op['q4']
            h += self.embed({site: local})
        for i in range(n):
            for j in range(i+1, n):
                if self.bonds[i, j] != 0:
                    h -= self.bonds[i, j]*self.embed(
                        {i: self.operators[i]['q'], j: self.operators[j]['q']})
        h.eliminate_zeros()
        self.hamiltonian = h.tocsr()

    def embed(self, operators):
        out = sparse.csr_matrix([[1.]])
        for site in reversed(range(len(self.dimensions))):
            factor = operators.get(site, sparse.eye(self.dimensions[site], format='csr'))
            out = sparse.kron(out, sparse.csr_matrix(factor), format='csr')
        return out

    def evolve(self, state, duration):
        duration = _real(duration, 'calibrated clock interval')
        scale = -1j*duration
        return expm_multiply(scale*self.hamiltonian, state,
                             traceA=scale*self.hamiltonian.diagonal().sum())

    def history(self, duration, *, read_site=0, sign=-1, feedback_phase=np.pi/2):
        """Two complete binary reads on one disturbed state, with first-bit feedback."""
        duration = _real(duration, 'calibrated clock interval')
        feedback_phase = _real(feedback_phase, 'feedback phase')
        if duration <= 0:
            raise ValueError('positive calibrated clock interval required')
        if (isinstance(read_site, (bool, np.bool_)) or not isinstance(read_site, Integral)
                or not 0 <= read_site < len(self.dimensions)):
            raise ValueError('read site must name an accessible tensor factor')
        if isinstance(sign, (bool, np.bool_)) or sign not in (-1, 1):
            raise ValueError('pointer sign must be +/-1')
        d = self.dimensions[read_site]
        if d < 5:
            raise ValueError('the supplied readout needs at least five local levels')
        vector = np.zeros(d, complex)
        vector[0], vector[4] = 1/5, sign*1j*np.sqrt(24)/5
        yes = np.outer(vector, vector.conj())
        kraus = tuple(self.embed({read_site: k}) for k in (np.eye(d)-yes, yes))
        feedback = np.eye(d, dtype=complex)
        feedback[4, 4] = np.exp(1j*feedback_phase)
        feedback = self.embed({read_site: feedback})
        blank = np.zeros(self.dimension, complex)
        blank[0] = 1
        state = self.evolve(blank, duration)
        branches = np.column_stack([k@state for k in kraus])
        branches[:, 1] = feedback@branches[:, 1]
        branches = self.evolve(branches, duration)
        rows = []
        for first in (0, 1):
            for last, k in enumerate(kraus):
                final = k@branches[:, first]
                rows.append({'bits': [first, last], 'probability': float(np.vdot(final, final).real)})
        return rows
