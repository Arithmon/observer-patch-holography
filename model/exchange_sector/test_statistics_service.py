"""Independent endpoint/phase oracles and hostile source-service admissions."""
from copy import deepcopy
import math
import warnings

import numpy as np
import pytest

from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from . import statistics as exact, statistics_service as service


def address(pair):
    return circuit.encode(tuple(int(p in pair) for p in range(12)))


def endpoint_hop(state, seam, character, controlled):
    """Independent direct generator exponential, used only as a target."""
    a, b = seam
    target = state.copy()
    offset = 512 if controlled else 0
    for spectator in range(12):
        if spectator in seam:
            continue
        left, right = address((a, spectator)), address((b, spectator))
        phase = -1j*(character if a < spectator < b else 1)
        target[offset+left] = phase*state[offset+right]
        target[offset+right] = phase*state[offset+left]
    return target


@pytest.fixture(scope='module')
def certificate():
    return service.certificate()


def test_fresh_certificate_complete_source_scope(certificate):
    assert certificate['layout']['record_selected_diagonals_checked'] == 180
    assert certificate['layout']['full_register_hopping_layouts_checked'] == 120
    assert certificate['layout']['full_address_integer_Walsh_inversion']
    for law in ('C', 'E'):
        for character in (-1, 1):
            witness = certificate['one_attempt_witnesses'][law][str(character)]
            assert witness['combined_raw_record_Stinespring_error_checked']
            assert witness['untouched_reference_dimension'] == 'arbitrary'
            for row in certificate['triangular_exchange']['observations'][law][str(character)].values():
                assert row['combined_complex_record_and_archive_amplitudes_checked']
    assert certificate['source_admission']['interval_or_hardware_certificate'] is False


@pytest.mark.parametrize('seam', [(0, 5), (0, 11), (5, 11)])
@pytest.mark.parametrize('character', [-1, 1])
@pytest.mark.parametrize('controlled', [False, True])
def test_complete_hop_on_unknown_complex_inputs_padding_and_reference(seam, character, controlled):
    rng = np.random.default_rng(5821)
    size = 1024 if controlled else 512
    state = rng.normal(size=(size, 3))+1j*rng.normal(size=(size, 3))
    state /= np.linalg.norm(state)
    actual = service.apply_hop(seam, state, character, controlled=controlled)
    target = endpoint_hop(state, seam, character, controlled)
    error = np.linalg.norm(actual-target)
    assert error < 2e-8
    assert error <= service.hop_error_bound(seam, character, controlled=controlled)+1e-12
    # All references participate in the same linear action.
    reference = np.array([[1, 1j], [.3, -.4j], [.2j, -.6]])
    replay = service.apply_hop(seam, state@reference, character, controlled=controlled)
    assert np.linalg.norm(replay-actual@reference) < 2e-11


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_output_phase_source_execution_has_exact_identity_complement(current):
    rng = np.random.default_rng(408)
    state = rng.normal(size=(512, 2))+1j*rng.normal(size=(512, 2))
    state /= np.linalg.norm(state)
    actual = service.apply_phase((0, 5), current, state, -1)
    target = state.copy()
    destination = 5 if current == 1 else 0
    if current:
        for spectator in (1, 2, 3, 4):
            target[address((destination, spectator))] *= -1
    assert np.linalg.norm(actual-target) < 2e-10
    assert np.linalg.norm(service.apply_phase((0, 5), current, state, 1)-state) == 0


@pytest.mark.parametrize('law', ['C', 'E'])
def test_complete_record_reference_instrument_with_independent_endpoint_oracle(law):
    rng = np.random.default_rng(9218)
    inputs = rng.normal(size=(66, 3))+1j*rng.normal(size=(66, 3))
    inputs /= np.linalg.norm(inputs)
    source = np.zeros((512, 3), complex)
    for pair, row in zip(exact.PAIRS, inputs, strict=True):
        source[address(pair)] = row
    actual = service.repair_branches(law, (0, 5), source, -1)
    target = {(coin, raw): np.zeros_like(source)
              for coin in (0, 1) for raw in range(4 if law == 'C' else 8)}
    for coin in (0, 1):
        for pair, row in zip(exact.PAIRS, inputs, strict=True):
            occupied = set(pair)
            total = int(0 in occupied)+int(5 in occupied)
            left = total//2+coin*(total % 2)
            current = int(0 in occupied)-left
            raw = {0: 0, 1: 1, -1: 2}[current]+(4*(total % 2) if law == 'E' else 0)
            final = occupied-{0, 5}
            if left:
                final.add(0)
            if total-left:
                final.add(5)
            sign = -1 if current and any(0 < s < 5 for s in pair if s not in (0, 5)) else 1
            target[coin, raw][address(final)] += sign*row/math.sqrt(2)
    assert len(actual) == (8 if law == 'C' else 16)
    error = np.linalg.norm(np.concatenate([row['amplitude']-target[row['coin'], row['record_code']]
                                         for row in actual]))
    assert error < 2e-8
    assert sum(row['abort'] for row in actual) == (2 if law == 'C' else 6)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_retained_loop_archive_destroys_visibility_without_changing_native_history(law):
    distributions = {}
    for character in (-1, 1):
        for dephase in (False, True):
            rows = service.loop_events(law, character, dephase=dephase)
            service._validate_grid(law, rows, loop=True, dephase=dephase)
            path = np.zeros(2)
            records = np.zeros((2, 4 if law == 'C' else 8))
            for row in rows:
                block = abs(row['amplitude'].reshape(2, 512))**2
                path += block.sum(axis=1)
                records[row['coin'], row['record_code']] += block.sum()
            wanted = [.5, .5] if dephase else [1, 0] if character == 1 else [0, 1]
            np.testing.assert_allclose(path, wanted, atol=2e-8)
            distributions[character, dephase] = records
    for records in distributions.values():
        np.testing.assert_allclose(records, distributions[1, False], atol=2e-8)


def test_returned_array_metadata_and_words_do_not_poison_execution():
    source = np.zeros(512, complex)
    source[address((0, 11))] = 1
    first = service.apply_hop((0, 5), source)
    assert not first.flags.writeable
    # NumPy >=2.4 deprecates this caller attack; exercise it deliberately
    # without suppressing warnings from any production source execution.
    with warnings.catch_warnings(record=True) as notices:
        warnings.simplefilter('always', DeprecationWarning)
        first.strides = (0,)
    assert all(issubclass(w.category, DeprecationWarning) and 'strides' in str(w.message)
               for w in notices)
    fresh = service.apply_hop((0, 5), source)
    assert abs(fresh[address((5, 11))]+1j) < 2e-8
    word = service.phase_word((0, 5))
    word[0]['angle'] = 991.
    assert service.phase_word((0, 5))[0]['angle'] != 991.


@pytest.mark.parametrize('bad', [True, False, None, -2, 0, 2, 1., '-1', np.int64(-1)])
def test_character_requires_explicit_exact_integer(bad):
    for call in (lambda: service.phase_word((0, 5), character=bad),
                 lambda: service.witness_columns('C', bad),
                 lambda: service.loop_resources('C', bad),
                 lambda: service.hop_layout((0, 5), bad)):
        with pytest.raises(ValueError):
            call()


@pytest.mark.parametrize('bad', [True, None, 1., '1', 2, np.int64(1)])
def test_current_requires_exact_signed_integer(bad):
    with pytest.raises(ValueError):
        service.phase_word((0, 5), bad)


@pytest.mark.parametrize('bad', [(5, 0), [0, 5], (False, 5), (0, 2), (0, 0), None])
def test_seam_requires_native_ordered_integer_pair(bad):
    with pytest.raises(ValueError):
        service.hop_layout(bad)


@pytest.mark.parametrize('bad', [0, 1, None, 'false', np.bool_(False)])
def test_control_switch_requires_boolean(bad):
    with pytest.raises(ValueError):
        service.hop_layout((0, 5), controlled=bad)
    with pytest.raises(ValueError):
        service.loop_resources('C', dephase=bad)


@pytest.mark.parametrize('bad', [np.zeros(511), np.zeros((512, 0)), np.ones(512, bool),
                               np.full(512, np.nan), np.full(512, np.inf),
                               np.ma.array(np.zeros(512)), [False]*512, [None]*512])
def test_state_admission_rejects_missing_coerced_or_nonfinite_input(bad):
    with pytest.raises(ValueError):
        service.apply_phase((0, 5), 1, bad)


@pytest.mark.parametrize('fault', ['drop', 'duplicate', 'bool_coin', 'bool_abort', 'null_current', 'nan', 'shape'])
def test_public_repair_rejects_corrupt_raw_grid_before_feedback(monkeypatch, fault):
    source = np.zeros(512, complex)
    source[address((0, 11))] = 1
    producer = service.recorded.repair_branches

    def broken(law, state):
        rows = deepcopy(producer(law, state))
        if fault == 'drop':
            rows.pop()
        elif fault == 'duplicate':
            rows[-1] = deepcopy(rows[0])
        elif fault == 'bool_coin':
            rows[0]['coin'] = False
        elif fault == 'bool_abort':
            rows[0]['abort'] = 0
        elif fault == 'null_current':
            rows[0]['current'] = None
        elif fault == 'nan':
            rows[0]['amplitude'][0] = np.nan
        else:
            rows[0]['amplitude'] = np.zeros(511)
        return rows

    monkeypatch.setattr(service.recorded, 'repair_branches', broken)
    with pytest.raises(ArithmeticError):
        service.repair_branches('C', (0, 5), source)


@pytest.mark.parametrize('fault', ['scalar', 'angle', 'wire', 'cnot'])
def test_emitted_walsh_word_audit_detects_phase_and_parity_tampering(fault):
    quarters = service.phase_quarters((0, 5))
    word = service.diagonal_word(quarters)
    if fault == 'scalar':
        assert word[0]['kind'] == 'I'
        word.pop(0)
    elif fault == 'angle':
        word[0]['angle'] *= -1
    elif fault == 'wire':
        phase = next(x for x in word if x['kind'] == 'Z')
        phase['targets'] = [(phase['targets'][0]+1) % 9]
    else:
        word.pop(next(i for i, x in enumerate(word) if x['kind'] == 'J'))
    with pytest.raises(ArithmeticError):
        service._audit_diagonal_word(quarters, word)


@pytest.mark.parametrize('fault', ['phase', 'nan', 'missing'])
def test_source_gray_guard_rejects_finite_wrong_global_phase_and_nonfinite(monkeypatch, fault):
    def broken(word, bits):
        if fault == 'nan':
            return np.full((512, 512), np.nan)
        if fault == 'missing':
            return np.eye(511)
        target = np.eye(512, dtype=complex)
        target[[0, 256]] = target[[256, 0]]
        return 1j*target
    monkeypatch.setattr(compiler, 'execute', broken)
    # Unique word key avoids touching a prior valid cached execution.
    key = (('I', 913.751, (0,)),)
    with pytest.raises(ArithmeticError):
        service._executed_gray(key)


def test_resource_accounts_include_all_paths_and_interventions(certificate):
    for law in ('C', 'E'):
        for character in (-1, 1):
            resource = service.witness_resources(law, character)
            assert resource['blank_source_qubits'] == 9
            assert len(resource['branches']) == (8 if law == 'C' else 16)
            for row in resource['branches']:
                assert row['native_pulses'] == row['proper_code_load_events']
                assert row['native_pulses'] == row['proper_code_export_events']
                assert row['proper_code_service_events'] == 2*row['native_pulses']
                assert row['proper_code_complement_outcome_slots'] == (
                    4*row['single_code_pulses']+2*row['pair_code_pulses'])
            ordinary = service.loop_resources(law, character)
            dephased = service.loop_resources(law, character, dephase=True)
            assert ordinary['native_pulses'] == dephased['native_pulses']
            assert ordinary['blank_coherent_path_qubits'] == 1
            assert dephased['retained_intermediate_path_measurements'] == 1
            assert dephased['dephasing_archive_slots'] == 2


def test_literal_carrier_address_binding_rejects_shared_chart_permutation(monkeypatch):
    encode = circuit.encode
    left, right = address((4, 8)), address((4, 9))
    def poisoned(state):
        value = encode(state)
        return right if value == left else left if value == right else value
    monkeypatch.setattr(circuit, 'encode', poisoned)
    with pytest.raises(ArithmeticError, match='literal carrier binding'):
        service._source_address_certificate()


def test_complete_export_admission_rejects_missing_zero_complement(monkeypatch):
    transfer = compiler.transfer_kraus
    monkeypatch.setattr(compiler, 'transfer_kraus', lambda dimension: transfer(dimension)[:-1])
    with pytest.raises(ArithmeticError, match='complete proper-code export'):
        service._proper_code_certificate()
