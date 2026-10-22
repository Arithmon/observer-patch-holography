"""Source-word realization of a declared charged relational readout.

The protocol uses the existing six-level processor controls and C/E repair
words.  It is an ideal encoded execution: proper-code load/export, reset and
pointer acquisition retain their separately declared apparatus contracts.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction as F
import math

import numpy as np

from ..federated_quantum_sector import construction as quantum
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler as source


TOLERANCE = 3e-11
CURRENT_CODES = {0: 0, 1: 1, -1: 2}
LOCAL_CODES = (0, 2, 1, 4)  # (0,0), (1,0), (0,1), (1,1) in the existing seam code


def preparation_word():
    """Blank -> i (|one_a,0_ref> + |vac,1_ref>)/sqrt(2)."""
    return ([source.pulse('X', math.pi/2, [1]),
             source.pulse('Y', -math.pi/4, [3])]
            + source.cnot(3, 1))


def analyzer_word():
    """Compile the effect (I+X_a,ref)/2, fair on the orthogonal space.

    Rotate the selected pair to addresses zero/one, then perform a pointer
    Ramsey test.  The controlled phase is 0 on plus, pi on minus and pi/2
    elsewhere, so the pointer's real interference is exactly X_a,ref.
    The emitted instrument is specified by this word, not by only its POVM.
    """
    basis = circuit.basis_swap_word(8, 0, 4) + [source.pulse('Y', math.pi/4, [1])]
    phases = [0. if address < 16 or address == 16 else
              math.pi if address == 18 else math.pi/2
              for address in range(32)]
    return (basis + [source.pulse('Y', -math.pi/4, [4])]
            + source.compile_diagonal(phases)
            + [source.pulse('Y', math.pi/4, [4])]
            + source.inverse(basis))


@lru_cache(maxsize=1)
def _preparation():
    blank = np.zeros(16, complex)
    blank[0] = 1
    result = source.execute(preparation_word(), 4, blank)
    result.flags.writeable = False
    return result


@lru_cache(maxsize=1)
def _coin():
    result = source.execute([source.pulse('Y', -math.pi/4, [0])], 1, [1., 0.])
    result.flags.writeable = False
    return result


def _local_image(mask, coin):
    target, current, imbalance = quantum.native_image(mask, (0, 1), coin)
    return target, current, imbalance


def _record_code(law, current, imbalance):
    return CURRENT_CODES[current] + (4*(abs(imbalance) % 2) if law == 'E' else 0)


def _permutation_image(law, coin, code):
    """Exact basis-address action of the existing compiled repair swaps."""
    if law not in ('C', 'E') or coin not in (0, 1) or code not in LOCAL_CODES:
        raise ValueError('an admitted local repair column required')
    address = code
    for left, right in circuit.REPAIR_SWAPS['C', coin]:
        if address == left:
            address = right
        elif address == right:
            address = left
    if law == 'E' and ((address >> coin) & 1):
        address ^= 32
    return address


@lru_cache(maxsize=4)
def _repair_columns(law, coin):
    """Execute the existing C/E source word on every local binary column."""
    if law not in ('C', 'E') or type(coin) is not int or coin not in (0, 1):
        raise ValueError('a C/E law and exact binary placement coin required')
    bits = 3+circuit.RECORD_BITS[law]
    blank = np.zeros((1 << bits, 4), complex)
    for mask, code in enumerate(LOCAL_CODES):
        blank[code, mask] = 1
    result = source.execute(circuit.repair_word(law, coin), bits, blank)
    for mask in range(4):
        target, current, imbalance = _local_image(mask, coin)
        expected = LOCAL_CODES[target]+(_record_code(law, current, imbalance) << 3)
        if _permutation_image(law, coin, LOCAL_CODES[mask]) != expected:
            raise ArithmeticError('exact repair swaps disagree with native vacuum/one columns')
        column = np.zeros(1 << bits, complex)
        column[expected] = 1
        if np.max(abs(result[:, mask]-column)) > TOLERANCE:
            raise ArithmeticError('source repair word disagrees with native vacuum/one/binary columns')
    result.flags.writeable = False
    return result


@lru_cache(maxsize=1)
def _analyzer_columns():
    blank = np.zeros((32, 16), complex)
    blank[:16] = np.eye(16)
    result = source.execute(analyzer_word(), 5, blank)
    plus, minus = result[:16], result[16:]
    effect = np.eye(16, dtype=complex)/2
    effect[2, 8] = effect[8, 2] = 1/2
    if (np.linalg.norm(plus.conj().T@plus-effect, 2) > TOLERANCE
            or np.linalg.norm(minus.conj().T@minus-(np.eye(16)-effect), 2) > TOLERANCE
            or np.linalg.norm(result.conj().T@result-np.eye(16), 2) > TOLERANCE):
        raise ArithmeticError('compiled analyzer has the wrong complete two-outcome effect')
    result.flags.writeable = False
    return result


def _expected_density(law, coin, code):
    """Independent exact matrix-unit branch, with source/reference coherence."""
    result = np.zeros((16, 16), complex)
    source_ref = ((1, 0), (0, 1))
    for left_mask, left_ref in source_ref:
        left, left_current, left_imbalance = _local_image(left_mask, coin)
        if _record_code(law, left_current, left_imbalance) != code:
            continue
        for right_mask, right_ref in source_ref:
            right, right_current, right_imbalance = _local_image(right_mask, coin)
            if _record_code(law, right_current, right_imbalance) == code:
                result[LOCAL_CODES[left]+8*left_ref, LOCAL_CODES[right]+8*right_ref] += .25
    return result


def _branches(law):
    """All coin, raw-record and pointer branches; no successful subset."""
    bits = circuit.RECORD_BITS[law]
    output = []
    prepared = _preparation()
    analyzer = _analyzer_columns()
    for coin in (0, 1):
        local = _repair_columns(law, coin)
        repaired = np.zeros((1 << (3+bits+1)), complex)
        for ref in (0, 1):
            for source_mask in range(4):
                repaired[(1 << (3+bits))*ref:(1 << (3+bits))*(ref+1)] += (
                    _coin()[coin]*prepared[LOCAL_CODES[source_mask]+8*ref]*local[:, source_mask])
        for code in range(1 << bits):
            vector = np.zeros(16, complex)
            for ref in (0, 1):
                for source_mask in range(8):
                    vector[source_mask+8*ref] = repaired[source_mask+(code << 3)+(ref << (3+bits))]
            if max(abs(vector[index]) for index in range(16) if index not in (1, 2, 8)) > TOLERANCE:
                raise ArithmeticError('source/reference total-one sector leaked under native repair')
            expected = _expected_density(law, coin, code)
            error = np.linalg.norm(np.outer(vector, vector.conj())-expected, 'fro')
            if error > TOLERANCE:
                raise ArithmeticError('an executed raw record branch differs from the exact native unit')
            final = analyzer@vector
            inverse = circuit.record_values(law, coin)
            for pointer, amplitude in enumerate((final[:16], final[16:])):
                output.append({'coin': coin, 'raw_record': code,
                               'current_or_contrast': inverse.get(code),
                               'abort': code not in inverse,
                               'pointer': pointer,
                               'probability': float(np.vdot(amplitude, amplitude).real),
                               'density_error': float(error)})
    return output


def source_word_readout(seams):
    """Source-executed one-seam gate and exact uniform-seam aggregate."""
    seams = tuple(seams)
    if (len(seams) != 30 or len(set(seams)) != 30
            or any(type(a) is not int or type(b) is not int or not 0 <= a < b < 12
                   for a, b in seams)):
        raise ValueError('thirty distinct oriented twelve-port seams required')
    incident = sum(a == 0 for a, _ in seams)
    if incident != 5 or {port for seam in seams for port in seam} != set(range(12)):
        raise ValueError('the native degree-five port-zero carrier required')
    rates = {}
    for law in ('C', 'E'):
        exact = quantum.native_unit(law, 1, 0, (0, 1))
        if any(key != (1, 0) for key in exact):
            raise ArithmeticError('incident-seam charge coherence moved or split')
        rates[law] = exact.get((1, 0), F(0))
    remote_source_word_columns = 0
    remote_exact_units = 0
    for law in ('C', 'E'):
        for coin in (0, 1):
            vacuum = np.zeros(1 << (3+circuit.RECORD_BITS[law]), complex)
            vacuum[0] = 1
            if np.linalg.norm(_repair_columns(law, coin)[:, 0]-vacuum) > TOLERANCE:
                raise ArithmeticError('a remote-seam source word changed its vacuum factor')
            remote_source_word_columns += 1
        for seam in seams:
            if 0 not in seam:
                if quantum.native_unit(law, 1, 0, seam) != {(1, 0): F(1)}:
                    raise ArithmeticError('a remote native seam changed source/reference coherence')
                remote_exact_units += 1
    rows = {}
    maximum_error = 0.
    for law in ('C', 'E'):
        branches = _branches(law)
        total = sum(row['probability'] for row in branches)
        plus = sum(row['probability'] for row in branches if row['pointer'] == 0)
        expected = float((1+rates[law])/2)
        if abs(total-1) > TOLERANCE or abs(plus-expected) > TOLERANCE:
            raise ArithmeticError('source-executed relational readout lost its native rate')
        maximum_error = max(maximum_error, max(row['density_error'] for row in branches),
                            abs(total-1), abs(plus-expected))
        rows[law] = {
            'all_coin_raw_record_pointer_branches': len(branches),
            'zero_probability_abort_branches': sum(row['abort'] and abs(row['probability']) < TOLERANCE
                                                   for row in branches),
            'one_incident_seam_coherence_rate_exact': str(rates[law]),
            'one_incident_seam_plus_probability_exact': str((1+rates[law])/2),
            'source_executed_one_incident_seam_plus_probability': plus,
            'uniform_thirty_seam_plus_probability_exact': str((incident*(1+rates[law])/2+len(seams)-incident)/len(seams)),
        }
    return {
        'source_word_preparation': 'declared total-one source/reference coherent state',
        'source_word_analyzer_effect': '(I+X_port,ref)/2, fair on orthogonal subspace',
        'native_law': 'existing compiled C/E local repair with every raw record and coin retained',
        'source_code_qubits': 4,
        'pointer_qubits': 1,
        'source_preparation_pulses': len(preparation_word()),
        'source_analyzer_pulses': len(analyzer_word()),
        'source_repair_pulses': {law: len(circuit.repair_word(law, 0)) for law in ('C', 'E')},
        'exact_repair_column_maps_checked': 16,
        'incident_seam_count_per_port': incident,
        'carrier_seam_count': len(seams),
        'remote_seam_source_word_vacuum_columns_checked': remote_source_word_columns,
        'remote_seam_exact_coherence_units_checked': remote_exact_units,
        'rows': rows,
        'maximum_floating_branch_or_probability_error': maximum_error,
        'numerical_acceptance_tolerance': TOLERANCE,
        'proper_code_service_failures_physically_bounded': False,
        'preparation_analyzer_or_C_E_law_selected_by_native_records': False,
        'uniform_seam_scheduler_or_physical_clock_selected': False,
        'particle_pole_or_rest_mass_identified': False,
    }
