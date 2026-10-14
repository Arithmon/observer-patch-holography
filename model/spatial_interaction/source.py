"""Local operations emitted and executed by the retained twelve-port compiler.

Each site uses two occupation bits; the two-site code is nx+4*ny.
Only number-conserving operations restricted to total occupation <=2 are
exported. Padding is retained in full local source checks, not postselected.
"""
from functools import lru_cache
import math

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c

SECTORS = ((0,), (1, 4), (2, 5, 8))
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'identity_phase_pulses', 'proper_code_load_events',
              'proper_code_export_events', 'proper_code_service_events')


def bond_generator():
    """Unit-hopping graph Laplacian, zero on unused local sectors."""
    out = np.zeros((16, 16))
    out[np.ix_(SECTORS[1], SECTORS[1])] = [[1, -1], [-1, 1]]
    out[np.ix_(SECTORS[2], SECTORS[2])] = [
        [2, -math.sqrt(2), 0], [-math.sqrt(2), 2, -math.sqrt(2)],
        [0, -math.sqrt(2), 2]]
    return out


def rotation_basis():
    out = np.eye(16)
    out[np.ix_((2, 5), (2, 5))] = np.array([[1, -1], [1, 1]])/math.sqrt(2)
    return out


def _executed(word, qubits, target):
    raw = c.execute(word, qubits)
    matrix, stability = c.stabilize_executed_unitary(raw)
    residual = float(np.linalg.norm(raw-target, 2))
    if residual > 2e-10:
        raise ArithmeticError('full-phase local source word disagrees with independent target')
    resources = c.resources(word)
    return matrix, {
        'resources': {key: resources[key] for key in COUNT_KEYS},
        'full_phase_target_residual': residual,
        'numerical_stabilization': stability,
    }


@lru_cache(maxsize=64)
def bond(angle):
    """Executed exp(-i angle L_bond), independently checked on all 16 inputs."""
    generator = bond_generator()
    # Fix the analytic eigenbasis, including padding, rather than depending
    # on an eigensolver's arbitrary choice in degenerate eigenspaces.
    basis = np.eye(16)
    values = np.zeros(16)
    basis[np.ix_(SECTORS[1], SECTORS[1])] = np.array([[1, 1], [1, -1]])/math.sqrt(2)
    values[list(SECTORS[1])] = [0, 2]
    basis[np.ix_(SECTORS[2], SECTORS[2])] = np.array([
        [0.5, 1/math.sqrt(2), 0.5], [1/math.sqrt(2), 0, -1/math.sqrt(2)],
        [0.5, -1/math.sqrt(2), 0.5]])
    values[list(SECTORS[2])] = [0, 2, 4]
    change = c.compile_real_orthogonal(basis)
    word = c.inverse(change)+c.compile_diagonal(-float(angle)*values)+change
    return _executed(word, 4, expm(-1j*float(angle)*generator))


@lru_cache(maxsize=32)
def onsite(angle):
    """Executed contact phase for n=0,1,2; unused n=3 padding is fixed."""
    phases = -float(angle)*np.array([0., 0., 1., 0.])
    matrix, info = _executed(c.compile_diagonal(phases), 2, np.diag(np.exp(1j*phases)))
    info['numerical_offdiagonal_compression_residual'] = float(np.linalg.norm(
        matrix-np.diag(np.diag(matrix)), 2))
    return matrix, info


@lru_cache(maxsize=4)
def preparation():
    basis = rotation_basis()
    return _executed(c.compile_real_orthogonal(basis), 4, basis)


@lru_cache(maxsize=8)
def feedback(angle):
    # Phase the split configuration relative to the two local doublons.
    phases = np.zeros(16)
    phases[5] = float(angle)
    return _executed(c.compile_diagonal(phases), 4, np.diag(np.exp(1j*phases)))


def _pointer_word(mask):
    phases = np.zeros(32)
    phases[16:] = np.pi*np.asarray(mask, dtype=float)
    hadamard = [c.pulse('H', math.pi/2, [4])]
    return hadamard+c.compile_diagonal(phases)+c.inverse(hadamard)


@lru_cache(maxsize=8)
def instrument(kind):
    """Retained binary pointer, optional Ramsey basis, and disturbed state.

    The initial blank pointer and its computational readout are apparatus.
    Both Kraus outcomes are returned; no branch is renormalized or erased.
    """
    if kind not in ('occupation', 'ramsey'):
        raise ValueError('unknown instrument')
    mask = np.array([i % 4 == 2 for i in range(16)]) if kind == 'occupation' else np.arange(16) == 2
    change = []
    basis = np.eye(16, dtype=complex)
    if kind == 'ramsey':
        phases = np.zeros(16)
        phases[5] = math.pi/2
        change = c.compile_real_orthogonal(rotation_basis())+c.compile_diagonal(phases)
        basis = np.diag(np.exp(1j*phases))@rotation_basis()
    word = c.inverse(change)+_pointer_word(mask)+change
    yes = basis@np.diag(mask.astype(float))@basis.conj().T
    no = np.eye(16)-yes
    target = np.block([[no, yes], [yes, no]])
    matrix, info = _executed(word, 5, target)
    kraus = (matrix[:16, :16], matrix[16:, :16])
    info['completeness_residual'] = float(np.linalg.norm(
        sum(k.conj().T@k for k in kraus)-np.eye(16), 2))
    return kraus, info
