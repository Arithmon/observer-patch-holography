"""Local canonical scalar operations executed by the twelve-port compiler.

All exported matrices are contractions of emitted source pulses, checked on
every input column against independent finite targets.  The targets never
replace an executed matrix.  The supplied proper-code load/export services,
blank pointers and computational pointer readout have the inherited scope;
floating diagnostics are not physical gate-error certificates.

Occupation codes are little endian.  A bond on (left, right) uses the index
n_left + d_left*n_right, including when the two dimensions differ.  There
is no number-sector projection or branch normalization in this module.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import math
from numbers import Integral, Real

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c
from .lattice import ScalarParameters, local_operators


COUNT_KEYS = (
    "native_pulses", "single_code_pulses", "pair_code_pulses",
    "identity_phase_pulses", "proper_code_load_events",
    "proper_code_export_events", "proper_code_service_events",
)


def _dimension(d):
    if (isinstance(d, (bool, np.bool_)) or not isinstance(d, Integral)
            or d not in (2, 4, 8)):
        raise ValueError("source occupation dimension must be 2, 4 or 8")
    return int(d)


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return value


def _parameters(parameters):
    if parameters is None:
        return ScalarParameters()
    if not isinstance(parameters, ScalarParameters):
        raise ValueError("parameters must be a ScalarParameters instance")
    return parameters


def _basis(matrix):
    values, vectors = np.linalg.eigh(matrix)
    # Fix column signs by the first supported component, avoiding argmax
    # ties that can change sign across LAPACK implementations.
    vectors[abs(vectors) < 1e-14] = 0
    for column in vectors.T:
        pivot = int(np.flatnonzero(abs(column) > 1e-12)[0])
        if column[pivot] < 0:
            column *= -1
    return values, vectors


def _counts(word):
    resources = c.resources(word)
    # The emitted schedule retains every zero-angle request.  Its exact
    # floating zero diagnostic is not pinned, since eigenvector dust can
    # differ across backends without changing the pulse/service count.
    return {key: resources[key] for key in COUNT_KEYS}


def _execute(stages, qubits, target):
    word = [event for _, stage in stages for event in stage]
    raw = c.execute(word, qubits)
    residual = float(np.linalg.norm(raw-target, 2))
    if not math.isfinite(residual) or residual > 2e-10:
        raise ArithmeticError("full-phase local source word disagrees with independent target")
    matrix, guard = c.stabilize_executed_unitary(raw)
    matrix.flags.writeable = False
    stage_counts = [{"stage": name, "resources": _counts(stage)} for name, stage in stages]
    resources = {key: sum(stage["resources"][key] for stage in stage_counts) for key in COUNT_KEYS}
    return matrix, {
        "resources": resources,
        "source_stages": stage_counts,
        "local_register_qubits": qubits,
        "full_phase_target_residual": residual,
        "numerical_stabilization": guard,
        "service_failure_outcomes_retained": True,
        "valid_code_failure_probability_in_exact_arithmetic": 0,
        "physical_service_error_certified": False,
    }


@lru_cache(maxsize=64)
def _onsite(d, degree, dt, parameters):
    operators = local_operators(d, parameters.omega)
    hamiltonian = (operators["p2"]/2
                   + (parameters.mass**2+parameters.kappa*degree)*operators["q2"]/2
                   + parameters.g*operators["q4"])
    values, basis = _basis(hamiltonian)
    change = c.compile_real_orthogonal(basis)
    stages = [
        ("onsite_inverse_basis", c.inverse(change)),
        ("onsite_eigenvalue_phase", c.compile_diagonal(-dt*values)),
        ("onsite_basis", change),
    ]
    return _execute(stages, d.bit_length()-1, expm(-1j*dt*hamiltonian))


def onsite(d, degree, dt, parameters=None):
    """Execute exp(-i dt h_site), preserving its scalar phase and top levels.

    The degree diagonal belongs to the open-graph gradient term.  The
    matrices of q^2, p^2 and q^4 are compressed infinite-oscillator powers.
    ``dt`` is a supplied Hamiltonian model-coordinate duration.
    """
    d = _dimension(d)
    if (isinstance(degree, (bool, np.bool_)) or not isinstance(degree, Integral)
            or degree < 0):
        raise ValueError("a nonnegative integer site degree required")
    matrix, metadata = _onsite(d, int(degree), _real(dt, "model duration"), _parameters(parameters))
    return matrix.copy(), deepcopy(metadata)


@lru_cache(maxsize=64)
def _bond(dleft, dright, dt, parameters):
    left = local_operators(dleft, parameters.omega)["q"]
    right = local_operators(dright, parameters.omega)["q"]
    left_values, left_basis = _basis(left)
    right_values, right_basis = _basis(right)
    left_bits, right_bits = dleft.bit_length()-1, dright.bit_length()-1
    left_change = c.compile_real_orthogonal(left_basis)
    right_change = c.remap(c.compile_real_orthogonal(right_basis),
                           list(range(left_bits, left_bits+right_bits)))
    # Separate local basis words avoid compiling a generic orthogonal matrix
    # of size dleft*dright, which can exceed the compiler's four-qubit bound.
    phases = dt*parameters.kappa*np.kron(right_values, left_values)
    stages = [
        ("bond_inverse_left_basis", c.inverse(left_change)),
        ("bond_inverse_right_basis", c.inverse(right_change)),
        ("bond_product_phase", c.compile_diagonal(phases)),
        ("bond_right_basis", right_change),
        ("bond_left_basis", left_change),
    ]
    target = expm(1j*dt*parameters.kappa*np.kron(right, left))
    return _execute(stages, left_bits+right_bits, target)


def bond(dleft, dright, dt, parameters=None):
    """Execute the bond unitary exp(+i dt kappa q_left q_right).

    The local basis index is n_left + dleft*n_right.  Both unequal code
    dimensions and the complete local occupation boundary are retained.
    """
    matrix, metadata = _bond(_dimension(dleft), _dimension(dright),
                             _real(dt, "model duration"), _parameters(parameters))
    return matrix.copy(), deepcopy(metadata)


def _origin_dimension(d):
    if _dimension(d) != 8:
        raise ValueError("the declared origin instrument uses the eight-level code")
    return 8


def _pointer_word(d):
    pointer = d.bit_length()-1
    phases = np.zeros(2*d)
    phases[d] = math.pi
    hadamard = [c.pulse("H", math.pi/2, [pointer])]
    return hadamard+c.compile_diagonal(phases)+c.inverse(hadamard)


@lru_cache(maxsize=2)
def _instrument(d, sign):
    real_basis = np.eye(d)
    cosine, sine = 1/5, math.sqrt(24)/5
    real_basis[np.ix_((0, 4), (0, 4))] = [[cosine, -sine], [sine, cosine]]
    phases = np.zeros(d)
    phases[4] = sign*math.pi/2
    change = c.compile_real_orthogonal(real_basis)+c.compile_diagonal(phases)
    stages = [
        ("pointer_inverse_field_basis", c.inverse(change)),
        ("pointer_controlled_flip", _pointer_word(d)),
        ("pointer_restore_field_basis", change),
    ]
    vector = np.zeros(d, complex)
    vector[0], vector[4] = cosine, sign*1j*sine
    yes = np.outer(vector, vector.conj())
    no = np.eye(d)-yes
    target = np.block([[no, yes], [yes, no]])
    matrix, metadata = _execute(stages, d.bit_length(), target)
    kraus = (matrix[:d, :d], matrix[d:, :d])
    metadata["pointer_bit_one"] = "projector onto (|0> + sign*i*sqrt(24)|4>)/5"
    metadata["projector_phase_sign"] = sign
    metadata["completeness_residual"] = float(np.linalg.norm(
        sum(k.conj().T@k for k in kraus)-np.eye(d), 2))
    return kraus, metadata


def instrument(d=8, sign=-1):
    """Return both disturbed field branches of a blank binary pointer read.

    Bit one reads P_e, e=(|0>+sign*i*sqrt(24)|4>)/5; bit zero reads I-P_e.
    The circuit changes basis, flips the pointer only for field code zero,
    and restores the field basis.  No outcome branch is normalized.
    """
    d = _origin_dimension(d)
    if isinstance(sign, (bool, np.bool_)) or not isinstance(sign, Integral) or sign not in (-1, 1):
        raise ValueError("projector phase sign must be -1 or +1")
    kraus, metadata = _instrument(d, int(sign))
    return tuple(k.copy() for k in kraus), deepcopy(metadata)


@lru_cache(maxsize=16)
def _feedback(d, phase):
    phases = np.zeros(d)
    phases[4] = phase
    return _execute([("feedback_occupation_four_phase", c.compile_diagonal(phases))],
                    d.bit_length()-1, np.diag(np.exp(1j*phases)))


def feedback(d=8, phase=math.pi/2):
    """Execute exp(i phase |4><4|) for record-controlled continuation."""
    matrix, metadata = _feedback(_origin_dimension(d), _real(phase, "feedback phase"))
    return matrix.copy(), deepcopy(metadata)
