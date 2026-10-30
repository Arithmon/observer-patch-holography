"""Complete source-operated two-point spectroscopy on one/two field sites.

The ground or Gibbs preparation, both Ramsey outputs, the purifier and all
first export failures are retained. Coherent probe controls and logical
encoding are declared apparatus; CAR parity is not claimed to prohibit their
encoded intermediate operations or to identify a laboratory fermion source.
"""
from __future__ import annotations

from fractions import Fraction as F
from dataclasses import dataclass
import math
from types import MappingProxyType
import numpy as np
from scipy.linalg import expm

from .spectrum import FieldSpectrum
from ..full_fock_field import field
from ..full_fock_field.local_simulation import split_bound
from ..field_thermal.thermal import execute_with_exports, failure_probability
from ..observer.source_scalar import compiler
from ..reference_observation.observation import frozen
from ..interacting_field.lattice import real


HADAMARD = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
PAULI = {'I': np.eye(2), 'X': np.array([[0,1],[1,0]]),
         'Y': np.array([[0,-1j],[1j,0]]), 'Z': np.diag([1,-1])}


@dataclass(frozen=True)
class _Preparation:
    state: np.ndarray
    target: np.ndarray
    failures: tuple
    error: float
    resources: object


def pauli_matrix(label):
    if type(label) is not str or not label or any(x not in PAULI for x in label):
        raise ValueError('nonempty tensor Pauli string required')
    out = np.ones((1, 1), complex)
    for letter in reversed(label): out = np.kron(out, PAULI[letter])
    return out


def _admit(spectrum):
    if not isinstance(spectrum, FieldSpectrum) or spectrum.geometry.volume not in (1, 2):
        raise ValueError('source spectroscopy supports one/two complete field sites')


def prepare(spectrum):
    _admit(spectrum)
    state = spectrum.reference_root.T.reshape(-1)
    word = compiler.compile_real_state(np.abs(state))+compiler.compile_diagonal(np.angle(state))
    initial = np.zeros_like(state); initial[0] = 1
    actual, failures = execute_with_exports(word, 2*spectrum.geometry.modes, initial)
    error = math.sqrt(float(np.linalg.norm(actual-state))**2+failure_probability(failures))
    if not math.isfinite(error) or error > 2e-9:
        raise ArithmeticError('source reference preparation failed numerical admission')
    return _Preparation(frozen(actual), frozen(state), tuple(failures), error,
                        MappingProxyType(compiler.resources(word)))


def _pointer_phase(pointer, phase):
    return [compiler.pulse('I', phase/2, [pointer]), compiler.pulse('Z', -phase/2, [pointer])]


def controlled_pauli(label, pointer):
    """C(P) from commuting native Pauli rotations, including its phase."""
    if type(pointer) is not int or pointer < len(label) or pointer > 8:
        raise ValueError('distinct bounded pointer bit required')
    pauli_matrix(label)
    base = list(label)+['I']*(pointer+1-len(label))
    coupled = base.copy(); coupled[pointer] = 'Z'
    return (compiler.pauli_rotation(base, math.pi/4)
            +compiler.pauli_rotation(coupled, -math.pi/4)
            +_pointer_phase(pointer, -math.pi/2))


def forward_word(spectrum, duration, steps):
    _admit(spectrum)
    if type(duration) not in (int, F) or type(steps) is not int or steps < 1:
        raise ValueError('exact rational model duration and positive split count required')
    word = list(field.source_word(spectrum.geometry, float(duration), steps))
    # N commutes with H and with every ideal split factor. This exact number
    # phase implements K=H-mu*N without fitting or changing the interaction.
    for mode in range(spectrum.geometry.modes):
        word += _pointer_phase(mode, spectrum.chemical_potential*float(duration))
    return word


def correlator(spectrum, first, second, duration=F(1, 4), steps=4, phase=0.):
    """Measure Re(exp(i*phase) Tr(rho P(t)Q)) with all pointer/abort outputs.

    P and Q are full tensor-Pauli observables. Chronological controlled Q,
    forward K, controlled P and reverse K give the branch operator U* P U Q.
    The joint reference is prepared by the same source interface.
    """
    return _correlator(spectrum, first, second, duration, steps, phase, prepare(spectrum))


def _correlator(spectrum, first, second, duration, steps, phase, prepared):
    """Internal reuse of a preparation executed in this call's source path.

    Public callers cannot supply a state, target or claimed resource counter.
    Each independently counted setting starts with the same literal preparation.
    """
    _admit(spectrum)
    modes = spectrum.geometry.modes
    if not isinstance(first, str) or not isinstance(second, str) or len(first) != modes or len(second) != modes:
        raise ValueError('Pauli strings must cover the complete field')
    p, q = pauli_matrix(first), pauli_matrix(second)
    phase = real(phase, 'pointer phase')
    if not math.isfinite(phase):
        raise ValueError('finite pointer phase required')
    target = spectrum.reference_root.T.reshape(-1)
    if not isinstance(prepared, _Preparation) or not np.array_equal(prepared.target, target):
        raise ValueError('prepared reference does not match this spectral model')
    pointer = 2*modes
    hadamard = compiler.remap(compiler.compile_real_orthogonal(HADAMARD), [pointer])
    forward = forward_word(spectrum, duration, steps)
    tape = (hadamard+controlled_pauli(second, pointer)+forward
            +controlled_pauli(first, pointer)+compiler.inverse(forward)
            +_pointer_phase(pointer, phase)+hadamard)
    initial = np.zeros(2*spectrum.dimension**2, complex)
    initial[:spectrum.dimension**2] = prepared.state
    actual, failures = execute_with_exports(tape, pointer+1, initial)
    outputs = tuple(frozen(x) for x in actual.reshape(2, -1))
    probabilities = [float(np.vdot(x, x).real) for x in outputs]
    prefix_abort = failure_probability(prepared.failures)
    current_abort = failure_probability(failures)
    full_mass = sum(probabilities)+prefix_abort+current_abort
    unitary = expm(-1j*float(duration)*spectrum.generator)
    operator = unitary.conj().T@p@unitary@q
    columns = target.reshape(spectrum.dimension, spectrum.dimension).T
    changed = (operator@columns).T.reshape(-1)
    exact = ((target+np.exp(1j*phase)*changed)/2, (target-np.exp(1j*phase)*changed)/2)
    error = math.sqrt(sum(float(np.linalg.norm(a-b))**2 for a,b in zip(outputs, exact, strict=True))
                      +prefix_abort+current_abort)
    eta = split_bound(spectrum.geometry.sides, duration, steps,
                      F(spectrum.geometry.speed), F(spectrum.geometry.interaction), refined=True)
    expected = float((np.exp(1j*phase)*np.trace(spectrum.reference@operator)).real)
    mean = probabilities[0]-probabilities[1]
    if (not all(math.isfinite(x) for x in [full_mass, error, expected, mean, float(eta)])
            or abs(full_mass-1) > 2e-9 or error > 2*float(eta)+4e-9):
        raise ArithmeticError('complete source spectroscopy failed its numerical/error admission')
    return dict(pointer_outputs=outputs, first_export_failures=failures,
                preparation_first_export_failures=prepared.failures,
                pointer_probabilities=probabilities, preparation_abort_probability=prefix_abort,
                continuation_abort_probability=current_abort, complete_probability=full_mass,
                measured_mean=mean, ideal_mean=expected,
                complete_vector_error=error, ideal_operator_error_bound=float(2*eta),
                preparation_error=prepared.error, resources=compiler.resources(tape),
                normalized_success_or_pointer_outcome=False)


def components(kind, index, modes):
    """A=(X+iY)/2 for a lowering probe, expressed in unitary Pauli words."""
    if kind == 'single' and type(index) is int and 0 <= index < modes:
        x = 'Z'*index+'X'+'I'*(modes-index-1)
        y = 'Z'*index+'Y'+'I'*(modes-index-1)
        return [(1., x)], [(1., y)]
    if kind == 'pair' and type(index) is int and 0 <= 2*index < modes:
        left, right = 'I'*(2*index), 'I'*(modes-2*index-2)
        return [(.5,left+'XX'+right),(-.5,left+'YY'+right)], [(.5,left+'XY'+right),(.5,left+'YX'+right)]
    raise ValueError('a valid single-mode or onsite-pair probe is required')


def probe_readout(spectrum, kind='single', index=0, duration=F(1,4), steps=4):
    """Reconstruct complete addition/removal correlations from counted records.

    U(1) invariance gives C_XX=C_YY and C_YX=-C_XY. Thus Gplus=(C_XX-iC_XY)/2
    and Gminus=conj((C_XX+iC_XY)/2). Every constituent Ramsey run is complete;
    this linear reconstruction is not conditioning on successful particles.
    """
    operator, charge = spectrum.probe(kind, index)
    measure = spectrum.measure(operator, charge)
    x, y = components(kind, index, spectrum.geometry.modes)
    prepared = prepare(spectrum)
    records, correlations = [], []
    total_pulses, worst_error = 0, 0.
    for second_terms in (x, y):
        value = 0j
        for a, first in x:
            for b, second in second_terms:
                quadratures = []
                for phase in (0., -math.pi/2):
                    result = _correlator(spectrum, first, second, duration, steps, phase, prepared)
                    quadratures.append(result['measured_mean'])
                    total_pulses += result['resources']['native_pulses']
                    worst_error = max(worst_error, result['complete_vector_error'])
                    records.append(dict(first=first,second=second,phase=phase,
                                        pointer_probabilities=result['pointer_probabilities'],
                                        preparation_abort_probability=result['preparation_abort_probability'],
                                        continuation_abort_probability=result['continuation_abort_probability'],
                                        complete_probability=result['complete_probability'],
                                        source_pulses=result['resources']['native_pulses']))
                value += a*b*complex(*quadratures)
        correlations.append(value)
    forward = (correlations[0]-1j*correlations[1])/2
    reverse = ((correlations[0]+1j*correlations[1])/2).conjugate()
    exact_forward, exact_reverse = (measure.correlation(float(duration), name) for name in ('forward','reverse'))
    residual = max(abs(forward-exact_forward), abs(reverse-exact_reverse))
    if not math.isfinite(residual) or residual > 4*worst_error+2e-9:
        raise ArithmeticError('source reconstructed Green function lost complete-measure agreement')
    encode = lambda z: [float(z.real), float(z.imag)]
    return dict(probe=kind,index=index,CAR_charge=charge,duration=str(duration),split_steps=steps,
                reference=spectrum.reference_kind,chemical_potential=spectrum.chemical_potential,
                observed_forward=encode(forward),observed_reverse=encode(reverse),
                complete_Lehmann_forward=encode(exact_forward),complete_Lehmann_reverse=encode(exact_reverse),
                maximum_correlation_residual=float(residual),maximum_complete_vector_error=worst_error,
                preparation_source_pulses_per_run=prepared.resources['native_pulses'],
                spectroscopy_source_pulses=total_pulses,
                independent_preparations_required=len(records),
                total_source_pulses_with_preparation=total_pulses+len(records)*prepared.resources['native_pulses'],
                initialized_bits=2*spectrum.geometry.modes+1,
                retained_purifier_bits=spectrum.geometry.modes,records=records,
                physical_apparatus_error_or_clock_calibrated=False)
