"""Exact action replay and resealed-history falsifiers for the q5 execution."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import scalar_execution_algebra as arithmetic
import source_scalar_execution as producer
import verify_source_scalar_execution as verifier


@pytest.fixture(scope='module')
def packet():
    return verifier.load()


@pytest.fixture(scope='module')
def model():
    return verifier.model()


def reseal(trace):
    """Repair every hash after a forged semantic edit, including the tail."""
    chain = '0'*64
    for event in trace['events']:
        event['ledger_parent'] = chain
        material = {k: event[k] for k in ('id', 'reads', 'write', 'ledger_parent')}
        chain = sha256(verifier.canonical(material)).hexdigest()
        event['hash'] = chain
    trace['final_hash'] = chain


def test_full_independent_replay(packet):
    result = verifier.verify(packet)
    assert result['events_replayed'] == 5888
    assert result['dynamic_field_reads_replayed'] == 34944
    lo, hi = map(F, result['full_mass_norm_error_final'])
    assert F(95, 10000) < lo < hi < F(96, 10000)
    assert not result['old_count_clock_applied']


def test_default_loader_receipt(packet):
    assert verifier.load(verifier.OUTPUT) == packet


def test_independent_exact_field_operations():
    for a in range(-4, 5):
        for b in range(-4, 5):
            x = arithmetic.Q(F(a, 3), F(b, 7))
            y = verifier.parse(x.encode())
            assert x.sign() == y.sign()
            assert (x*x).encode() == (y*y).encode()
            assert arithmetic.parse(y.encode()) == x
            if a or b:
                assert (1/x).encode() == (1/y).encode()
                assert x*(1/x) == arithmetic.Q(1)


@pytest.mark.parametrize('value', [True, 0.1, 1j, '1', None],
                         ids=['bool', 'float', 'complex', 'string', 'none'])
def test_exact_input_types(value):
    with pytest.raises(ValueError):
        arithmetic.Q(value)
    with pytest.raises(ValueError):
        verifier.R(value)


@pytest.mark.parametrize('value', [[True, '0'], ['1.0', '0'], ['2/2', '0'],
    ['1/0', '0'], ['0', '1', '2'], ['0'*3001, '0'], ['nan', '0']],
    ids=['bool', 'decimal', 'unreduced', 'zero-denominator', 'arity', 'oversize', 'nan'])
def test_strict_algebraic_loader(value):
    with pytest.raises((ValueError, ZeroDivisionError)):
        arithmetic.parse(value)
    with pytest.raises(ValueError):
        verifier.parse(value)


@pytest.mark.parametrize('raw', [b'{"a":1,"a":2}', b'{"x":1.0}', b'{"x":NaN}',
    b'{"x":Infinity}', b'\xff', b' '*12_000_001],
    ids=['duplicate', 'float', 'nan', 'infinity', 'encoding', 'size-budget'])
def test_strict_json_loader(tmp_path, raw):
    path = tmp_path/'bad.json'
    path.write_bytes(raw)
    with pytest.raises((ValueError, UnicodeDecodeError)):
        verifier.load(path)


@pytest.mark.parametrize('kind', ['missing-read', 'wrong-writer', 'wrong-version', 'wrong-resource',
    'wrong-read-value', 'wrong-write-value', 'late-write-value', 'extra-read', 'swapped-read',
    'old-memory-read', 'bool-version', 'false-count', 'extra-event-field', 'coherent-zero-intervention'])
def test_resealed_false_histories_rejected(packet, model, kind):
    trace = deepcopy(packet['traces']['ascending_intervention'])
    event = trace['events'][200]
    if kind == 'missing-read':
        event['reads'].pop()
    elif kind == 'wrong-writer':
        event['reads'][-1][2][1] = (event['reads'][-1][2][1]+1) % 64
    elif kind == 'wrong-version':
        event['reads'][-1][1] += 1
    elif kind == 'wrong-resource':
        event['reads'][-1][0] = (event['reads'][-1][0]+1) % 128
    elif kind == 'wrong-read-value':
        event['reads'][-1][3] = ['123', '0']
    elif kind == 'wrong-write-value':
        event['write'][3] = ['123', '0']
    elif kind == 'late-write-value':
        trace['events'][-1]['write'][3] = ['123', '0']
    elif kind == 'extra-read':
        event['reads'].append(deepcopy(event['reads'][-1]))
    elif kind == 'swapped-read':
        event['reads'][-1], event['reads'][-2] = event['reads'][-2], event['reads'][-1]
    elif kind == 'old-memory-read':
        event['reads'][0][1] = 1
        event['reads'][0][2][0] = 0
    elif kind == 'bool-version':
        event['reads'][-1][1] = True
    elif kind == 'false-count':
        trace['dynamic_field_reads'] -= 1
    elif kind == 'extra-event-field':
        event['hidden_read'] = event['reads'][0]
    else:
        trace = deepcopy(packet['traces']['ascending_baseline'])
    reseal(trace)
    with pytest.raises(ValueError):
        verifier.verify_trace(trace, model[1], model[2], list(range(64)))


@pytest.mark.parametrize('kind', ['clock', 'vacuum', 'probability', 'substeps', 'mass', 'source', 'pin-bool'])
def test_false_scope_or_source_rejected(packet, kind):
    fake = deepcopy(packet)
    if kind == 'clock':
        fake['scope']['old_count_clock_theorem_applied'] = True
    elif kind == 'vacuum':
        fake['scope']['modified_hamiltonian_vacuum_substituted'] = True
    elif kind == 'probability':
        fake['scope']['quantum_covariance_or_probability_queried'] = True
    elif kind == 'substeps':
        fake['model']['steps'] = 3
    elif kind == 'mass':
        fake['mass_Qphi'][0] = ['1', '0']
    elif kind == 'source':
        fake['source_addresses'][0]['D6_control'][0] += 1
    else:
        fake['source_pins'][verifier.PINS[0]]['bytes'] = True
    with pytest.raises(ValueError):
        verifier.verify(fake)


@pytest.mark.parametrize('kind', ['zero-tail', 'last-error', 'commutator', 'ancestry'])
def test_false_full_comparison_rejected(packet, kind):
    fake = deepcopy(packet)
    if kind == 'zero-tail':
        fake['diagnostics']['taylor_remainder_upper'] = '0'
    elif kind == 'last-error':
        fake['diagnostics']['comparisons'][-1]['full_mass_norm_discretization_error'] = ['0', '0']
    elif kind == 'commutator':
        fake['diagnostics']['comparisons'][-1]['canonical_commutator_over_i_hbar_Qphi'] = ['0', '0']
    else:
        fake['diagnostics']['actual_field_interval_counts'][-1]['inclusive_count'] -= 1
    with pytest.raises(ValueError):
        verifier.verify(fake)


def test_reordered_execution_is_value_equal_not_hash_equal(packet):
    for preparation in ('baseline', 'intervention'):
        a = packet['traces']['ascending_'+preparation]
        b = packet['traces']['descending_'+preparation]
        assert a['layers'] == b['layers']
        assert a['final_hash'] != b['final_hash']


def test_cfl_energy_and_nontrivial_intervention(packet, model):
    mass, rows, velocity, _, _ = model
    assert F(packet['model']['cfl_upper']) < 3 < 4
    # The primitive action itself is symmetric after mass weighting.
    for i, row in enumerate(rows):
        for j, coefficient in row:
            assert mass[i]*coefficient == mass[j]*dict(rows[j])[i]
    expected = verifier.inner(mass, velocity, velocity)/2
    assert verifier.parse(packet['diagnostics']['modified_energy_Qphi']) == expected
    witness = packet['diagnostics']['first_coarse_layer_outside_radius']
    assert verifier.parse(witness['field_response_Qphi']).sign() > 0
    assert (verifier.parse(witness['x_displacement_squared_Qphi'])-F(1, 5)).sign() > 0
    assert [r['inclusive_count'] for r in packet['diagnostics']['actual_field_interval_counts']] == [144, 576, 1024]
    assert len(packet['diagnostics']['comparisons']) == 21


def test_invalid_producer_schedule_fails_before_execution():
    with pytest.raises(ValueError):
        producer.execute([], [], [True]+list(range(1, 64)))
    with pytest.raises(ValueError):
        producer.execute([], [], list(range(63)))
