"""Small native Ramsey witness for the four-level-invisible B(N) action.

All eight computational outputs are retained.  The preparation and read use
the high bit of the existing three-qubit proper code; evolution compiles its
complete eight-level diagonal, including unused occupation levels.  No dense
eigenbasis compiler, postselection or target matrix is inserted in evolution.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c
from .response import finite


DIMENSION = 8
TIME = math.pi/2
BASE_COEFFICIENTS = (Fraction(0), Fraction(1), Fraction(1, 4), Fraction(1, 8))


def energies(extra=0.):
    extra = finite(extra, 'B(N) coefficient')
    return np.array([float(sum(coefficient*math.comb(n, k)
                              for k, coefficient in enumerate(BASE_COEFFICIENTS)))
                     +extra*math.comb(n, 4) for n in range(DIMENSION)])


def words(extra, time=TIME):
    time = finite(time, 'model-coordinate duration')
    if time < 0:
        raise ValueError('nonnegative model-coordinate duration required')
    # Fixed read phase uses the locked pair c4=0 versus c4=1, never actual
    # outputs or the tested candidate's value.  It is the midpoint phase.
    phase = time*(energies(0)[4]+energies(1)[4])/2
    return {
        'prepare': [c.pulse('Y', -math.pi/4, [2])],
        'evolve': c.compile_diagonal((-time*energies(extra)).tolist()),
        'read': [c.pulse('Z', -phase/2, [2]), c.pulse('X', math.pi/4, [2])],
    }


def execute(extra, time=TIME, *, commit_before_read=False):
    if type(commit_before_read) is not bool:
        raise ValueError('explicit Boolean commitment switch required')
    stages = words(extra, time)
    blank = np.zeros(DIMENSION, complex)
    blank[0] = 1
    prepared = c.execute(stages['prepare'], 3, blank)
    evolved = c.execute(stages['evolve'], 3, prepared)
    if commit_before_read:
        # Retain every pre-read occupation and every later read code.  These
        # are complete sequential record probabilities, not postselection.
        reader = c.execute(stages['read'], 3)
        joint = np.abs(reader)**2*np.abs(evolved)[None, :]**2
        probabilities = joint.sum(axis=1)
    else:
        probabilities = np.abs(c.execute(stages['read'], 3, evolved))**2
        joint = None
    tape = sum(stages.values(), [])
    matrix = c.execute(stages['evolve'], 3)
    target = np.diag(np.exp(-1j*float(time)*energies(extra)))
    return {
        'B_coefficient': float(extra), 'model_time': float(time),
        'committed_before_read': commit_before_read,
        'complete_read_probabilities': probabilities.tolist(),
        'joint_old_occupation_then_read_code': None if joint is None else joint.T.tolist(),
        'outcome_zero_probability': float(probabilities[0]),
        'complete_probability_mass': float(probabilities.sum()),
        'word_sha256': hashlib.sha256(json.dumps(tape, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        'resources': c.resources(tape),
        'stage_pulses': {name: len(word) for name, word in stages.items()},
        'full_phase_evolution_target_residual': float(np.linalg.norm(matrix-target, 2)),
        'numerical_residual_is_certificate': False,
    }


def independent_reference(extra, time=TIME):
    """Direct dense Hamiltonian exponential and analytic measurement vector.

    Neither source words nor the source execution matrices enter this route.
    The supplied preparation, Hamiltonian and detector convention are shared.
    """
    time = finite(time, 'model-coordinate duration')
    if time < 0:
        raise ValueError('nonnegative model-coordinate duration required')
    prepared = np.zeros(DIMENSION, complex)
    prepared[[0, 4]] = 1/math.sqrt(2)
    state = expm(-1j*time*np.diag(energies(extra)))@prepared
    midpoint = time*(energies(0)[4]+energies(1)[4])/2
    # R = exp(i*pi*X/4) exp(-i*midpoint*Z/2) on the high bit.
    x = np.array([[0, 1], [1, 0]], complex)
    z = np.diag([1., -1.])
    reader = np.kron(expm(1j*math.pi*x/4)@expm(-1j*midpoint*z/2), np.eye(4))
    return np.abs(reader@state)**2


def background_slope_check():
    """The B(N) derivative response persists with noncommuting quartic H0.

    This semantic check does not execute the complete coupled field.  It
    verifies the general commutator identity on a noncommuting local example.
    """
    a = np.diag(np.sqrt(np.arange(1, DIMENSION)), 1)
    q = (a+a.T)/math.sqrt(2)
    base = np.diag(np.arange(DIMENSION))+q@q@q@q/24
    extra = np.diag([math.comb(n, 4) for n in range(DIMENSION)])
    rows = []
    for n in range(1, 6):
        psi = np.zeros(DIMENSION, complex)
        detector = np.zeros(DIMENSION, complex)
        psi[[0, n]] = 1/math.sqrt(2)
        detector[0], detector[n] = 1/math.sqrt(2), -1j/math.sqrt(2)
        density = np.outer(psi, psi.conj())
        effect = np.outer(detector, detector.conj())
        slopes = [float(np.trace(effect@(-1j*(h@density-density@h))).real)
                  for h in (base, base+extra)]
        rows.append({'occupation': n, 'base_derivative': slopes[0],
                     'changed_derivative': slopes[1],
                     'twice_derivative_difference': 2*(slopes[1]-slopes[0]),
                     'exact_binomial_response': math.comb(n, 4)})
    return {'background': 'N+Q^4/24 with Q=(a+a_dagger)/sqrt(2), dimension 8',
            'background_commutes_with_added_term': False,
            'native_source_execution': False, 'records': rows}
