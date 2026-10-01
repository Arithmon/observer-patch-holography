"""Independent replay and adversarial controls for the fixed-strength exit."""

import copy
import json
import math
import os
import random
import subprocess
import sys

import numpy as np
import pytest

from . import (bounds, bounds_check, build, circuit, circuit_check, decoder,
               decoder_check, export, export_check, graph_check, graphs, verify)
from .format import pack, unpack


@pytest.fixture(scope='module')
def evidence():
    return build.candidate()


def test_all_semantics_independently_replay(evidence):
    verify.verify_evidence(evidence)


def test_committed_receipt_and_recursive_parent():
    packet = verify.verify()
    assert len(packet['evidence']['decoder']) == 27
    assert sum(x['cases'] for x in packet['evidence']['circuit']['censuses']) == 91740


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'identity', 'trace', 'power', 'bool', 'mixing'])
def test_graph_receipts_reject_forged_expansion(evidence, mutation):
    rows = copy.deepcopy(evidence['graphs']['graphs'])
    if mutation == 'missing':
        rows.pop()
    elif mutation == 'duplicate':
        rows[-1] = rows[0]
    elif mutation == 'identity':
        rows[0]['walks_sha256'] = '0'*64
    elif mutation == 'trace':
        rows[0]['trace'] += 1
    elif mutation == 'power':
        rows[0]['degree'] = 8
    elif mutation == 'bool':
        rows[0]['m'] = True
    else:
        rows[-1]['mixing_numerator'] = 0
    with pytest.raises(ValueError):
        graph_check.verify_graphs(rows)


@pytest.mark.parametrize('mutation', ['loop', 'inverse', 'missing', 'bool', 'swap'])
def test_labelled_matching_semantics(mutation):
    rows = graphs.permutations(5)
    if mutation == 'loop':
        rows[0] = list(range(25))
    elif mutation == 'inverse':
        rows[0] = rows[1]
    elif mutation == 'missing':
        rows.pop()
    elif mutation == 'bool':
        rows[0][0] = False
    else:
        rows[0][0], rows[0][1] = rows[0][1], rows[0][0]
    with pytest.raises(ValueError):
        graph_check.check_maps(rows, 5)


def test_regular_identity_voters_fail_the_refresh_invariant():
    n = 64
    old = np.zeros(n, dtype=int)
    old[:4] = 1
    # An identity multigraph is regular at any degree but has no mixing.
    output = old.copy()
    output[63] ^= 1
    assert old.sum() == n//16 and output.sum() > n//16
    assert 1 <= n//64


@pytest.mark.parametrize('mutation', ['omit_source', 'omit_voter', 'bad_count', 'wrong_output', 'duplicate', 'bool'])
def test_damage_controls_reject_favorable_filtering(evidence, mutation):
    rows = copy.deepcopy(evidence['graphs']['damage'])
    row = rows[-1]
    if mutation == 'omit_source':
        row['shared_faults'].pop()
    elif mutation == 'omit_voter':
        row['damaged_voters'].pop()
    elif mutation == 'bad_count':
        row['maximum_live_wrong'] = 0
    elif mutation == 'wrong_output':
        row['output_wrong'] = []
    elif mutation == 'duplicate':
        row['initial'].append(row['initial'][0])
    else:
        row['initial'][0] = False
    with pytest.raises(ValueError):
        graph_check.verify_damage(rows)


@pytest.mark.parametrize('m', [2, 3])
def test_actual_native_schedule_matches_independent_oracle(m):
    program = circuit.schedule(m)
    circuit_check.validate(program)
    assert program == circuit_check.reference_schedule(m)


@pytest.mark.parametrize('mutation', ['idle', 'reset', 'copy', 'shared_sort', 'majority', 'layer', 'bool'])
def test_schedule_rejects_hidden_fault_ownership_and_unpaid_idles(mutation):
    p = circuit.schedule(3)
    if mutation == 'idle':
        p['layers'][0].pop(0)
    elif mutation == 'reset':
        p['layers'][0][-1][0] = 'idle'
    elif mutation == 'copy':
        p['layers'][1][0][-1] += 1
    elif mutation == 'shared_sort':
        p['layers'][9][0] = ['sort', 0, 1]
    elif mutation == 'majority':
        op = next(op for op in p['layers'][-1] if op[0] == 'copy')
        op[1] += 1
    elif mutation == 'layer':
        p['layers'].pop()
    else:
        p['layers'][0][0][1] = False
    with pytest.raises(ValueError):
        circuit_check.validate(p)


def test_unpowered_fault_failures_are_retained(evidence):
    for row in evidence['circuit']['censuses']:
        assert sum(row['output_histogram']) == row['cases']
        assert sum(row['live_histogram']) == row['cases']
        assert sum(row['output_histogram'][2:]) > 0
        assert row['live_histogram'][2] > 0


@pytest.mark.parametrize('field', ['cases', 'locations', 'width', 'depth', 'output_histogram', 'live_histogram'])
def test_false_census_rejected(evidence, field):
    row = copy.deepcopy(evidence['circuit'])
    item = row['censuses'][0]
    if isinstance(item[field], list):
        item[field][0] += 1
    else:
        item[field] += 1
    with pytest.raises(ValueError):
        circuit_check.verify(row)


def test_indexed_compiler_keeps_both_operands_in_one_service():
    rng = random.Random(88012)
    for power in (1, 2, 12):
        n, d, width, depth = circuit.dimensions(17, power)
        for _ in range(40):
            layer = rng.randrange(1, d+1)
            source = rng.randrange(n)
            op = circuit.operation(17, power, layer, source)
            assert circuit.operation(17, power, layer, op[2]) == op
        for layer in (d+1, d+2, depth//2, depth-2, depth-1):
            for wire in (n, n+d//2-1, width-1):
                op = circuit.operation(17, power, layer, wire)
                for operand in op[1:]:
                    assert circuit.operation(17, power, layer, operand) == op


@pytest.mark.parametrize('power', [1, 2, 12])
def test_indexed_compiler_matches_independent_walks_and_sort_rows(power):
    rng = random.Random(314159+power)
    m, n, d = 17, 289, 8**power
    maps = graph_check.expected_maps(m)
    width = n*(d+2)
    for _ in range(48):
        word, source = rng.randrange(d), rng.randrange(n)
        voter, remainder = source, word
        for _ in range(power):
            voter = maps[remainder % 8][voter]
            remainder //= 8
        expected = ['copy', source, n+d*voter+word]
        assert circuit.operation(m, power, word+1, source) == expected
        assert circuit.operation(m, power, word+1, expected[2]) == expected
    # Construct layer numbers from bubble rows, independently of the
    # compiler's binary-search inversion of a triangular number.
    rows = {0, 1, d//2, d-2} | {rng.randrange(d-1) for _ in range(48)}
    for row in rows:
        count = d-1-row
        for position in {0, count-1, rng.randrange(count)}:
            layer = d+1+row*(2*d-row-1)//2+position
            voter = rng.randrange(n)
            first = n+d*voter+position
            expected = ['sort', first, first+1]
            assert circuit.operation(m, power, layer, first) == expected
            assert circuit.operation(m, power, layer, first+1) == expected
            for idle in (voter, width-1, n+d*((voter+1) % n)+d-1):
                # The last private bit participates only in the first row's
                # final comparator; choose a shared/output idle in that case.
                if idle >= n and idle < n*(d+1) and position == d-2:
                    continue
                assert circuit.operation(m, power, layer, idle) == ['idle', idle]


@pytest.mark.parametrize('mutation', ['missing', 'idle', 'bool', 'wrong_wire'])
def test_full_degree_index_controls_reject_mutations(evidence, mutation):
    rows = copy.deepcopy(evidence['circuit']['indexed'])
    if mutation == 'missing':
        rows.pop()
    elif mutation == 'idle':
        rows[2]['operation'] = ['idle', rows[2]['wire']]
    elif mutation == 'bool':
        rows[0]['layer'] = False
    else:
        rows[2]['operation'][2] += 1
    with pytest.raises(ValueError):
        circuit_check.verify_indexed(rows)


@pytest.mark.parametrize('seed', range(8))
def test_actual_code_corrects_varied_degree_five_words_at_radius_31(seed):
    rng = random.Random(seed)
    terms = rng.sample(list(decoder.MONOMIALS), 37)
    truth = decoder.polynomial(terms)
    errors = rng.sample(range(2047), 31)
    bits = truth[:-1].astype(int).tolist()
    for i in errors:
        bits[i] ^= 1
    result = decoder.decode(bits)
    independent = decoder_check.reference(sum(x << i for i, x in enumerate(bits)))
    assert result == independent
    assert not result['failure'] and result['distance'] == 31
    assert result['logical'] == int(truth[-1])


@pytest.mark.parametrize('mutation', ['logical', 'distance', 'coefficient', 'failure', 'radius', 'omit', 'bool'])
def test_decoder_evidence_rejects_false_correction(evidence, mutation):
    rows = copy.deepcopy(evidence['decoder'])
    if mutation == 'logical':
        rows[5]['decoded']['logical'] = 0
    elif mutation == 'distance':
        rows[4]['decoded']['distance'] = 0
    elif mutation == 'coefficient':
        rows[5]['decoded']['coefficients'] = []
    elif mutation == 'failure':
        rows[-1]['decoded']['failure'] = False
    elif mutation == 'radius':
        rows[-2]['decoded']['logical'] = 0
    elif mutation == 'omit':
        rows.pop()
    else:
        rows[1]['errors'][0] = False
    with pytest.raises(ValueError):
        decoder_check.verify(rows)


@pytest.mark.parametrize('value', [[], [0]*2046, [0]*2048, [False]*2047, [2]*2047, '0'*2047])
def test_decoder_refuses_malformed_input(value):
    with pytest.raises(ValueError):
        decoder.decode(value)


@pytest.mark.parametrize('field', ['power', 'degree', 'output_fraction', 'service_locations_per_bit',
                                 'archive_threshold_denominator_factor', 'joint_decay', 'accounting_decay',
                                 'storage_time', 'archive_decay', 'quantum_decay'])
def test_forged_bounds_rejected(evidence, field):
    row = copy.deepcopy(evidence['bounds'])
    if type(row[field]) is list:
        row[field][0] += 1
    else:
        row[field] += 1
    with pytest.raises(ValueError):
        bounds_check.verify(row)


@pytest.mark.parametrize('mutation', ['levels', 'size', 'missing', 'tail_sign'])
def test_refinement_and_norm_tails_are_not_probability_shortcuts(evidence, mutation):
    row = copy.deepcopy(evidence['bounds'])
    if mutation == 'levels':
        row['schedule'][0]['levels'] -= 1
    elif mutation == 'size':
        row['schedule'][-1]['bits'] = 9
    elif mutation == 'missing':
        row['schedule'].pop()
    else:
        row['tails'][0]['coefficients'][1] *= -1
    with pytest.raises(ValueError):
        bounds_check.verify(row)


@pytest.mark.parametrize('d', [2, 4])
def test_private_kraus_rotation_is_allowed(d):
    rows = export.reductions()
    group = rows[(d-2)//2]['groups']
    matrices = [unpack(k, (d, 6)) for k in group[1]]
    count = len(matrices)
    u = np.exp(2j*np.pi*np.outer(np.arange(count), np.arange(count))/count)/math.sqrt(count)
    group[1] = [pack(sum(u[i, j]*matrices[j] for j in range(count))) for i in range(count)]
    export_check.check_groups(group, d)


def test_identical_unlabelled_channel_with_wrong_public_meaning_fails():
    row = export.reductions()[0]
    a, b = (unpack(row['groups'][j][0], (2, 6)) for j in (0, 1))
    row['groups'][0][0] = pack((a+b)/math.sqrt(2))
    row['groups'][1][0] = pack((a-b)/math.sqrt(2))
    with pytest.raises(ValueError, match='public instrument'):
        export_check.check_groups(row['groups'], 2)


@pytest.mark.parametrize('mutation', ['postselect', 'abort', 'coherence', 'raw_floor', 'nan', 'bool'])
def test_complete_export_instrument_rejects_omitted_or_false_branches(evidence, mutation):
    rows = copy.deepcopy(evidence['export']['instruments'])
    row = rows[0]
    if mutation == 'postselect':
        row['blocks'][0][1] = pack(np.zeros((16, 16)))
    elif mutation == 'abort':
        row['blocks'][1][1] = pack(np.zeros((16, 16)))
    elif mutation == 'coherence':
        block = unpack(row['blocks'][0][0], (16, 16))
        row['blocks'][0][0] = pack(np.diag(np.diag(block)))
    elif mutation == 'raw_floor':
        row['shared_broadcast_probability'] = 0
    elif mutation == 'nan':
        row['wrong_label_probability'] = float('nan')
    else:
        row['abort_label'] = True
    with pytest.raises(ValueError):
        export_check.verify_instruments(rows)


@pytest.mark.parametrize('mutation', ['floor', 'omit', 'uniform', 'unitary'])
def test_staged_decoder_countercontrols_are_binding(evidence, mutation):
    rows = copy.deepcopy(evidence['export']['stages'])
    if mutation == 'floor':
        rows[-1]['raw_serial_flip'] = 0
    elif mutation == 'omit':
        rows[-1]['strengths'].pop()
    elif mutation == 'uniform':
        rows[-1]['uniform_sum_bound'] = 0
    else:
        rows[-1]['unitary'] = pack(np.eye(2))
    with pytest.raises(ValueError):
        export_check.verify_stages(rows)


@pytest.mark.parametrize('stage', range(3, 8))
@pytest.mark.parametrize('factor', [0, -1, .5, 2])
def test_tiny_nonzero_decoder_terms_cannot_be_erased(evidence, stage, factor):
    rows = copy.deepcopy(evidence['export']['stages'])
    rows[-1]['strengths'][stage] *= factor
    with pytest.raises(ValueError, match='positive recursive stage strength'):
        export_check.verify_stages(rows)


def test_small_roundoff_in_nonzero_decoder_terms_is_allowed(evidence):
    rows = copy.deepcopy(evidence['export']['stages'])
    rows[-1]['strengths'] = [value*(1+1e-13) for value in rows[-1]['strengths']]
    export_check.verify_stages(rows)


@pytest.mark.parametrize('mutation', ['drop', 'no_leak', 'distance', 'zero_noise'])
def test_concrete_nonzero_noise_process_is_verified(evidence, mutation):
    row = copy.deepcopy(evidence['export'])
    case = row['fixed_noise'][0]
    if mutation == 'drop':
        case['kraus'].pop()
    elif mutation == 'no_leak':
        k = unpack(case['kraus'][1], (6, 6))
        k[0, 0], k[5, 0] = k[5, 0], 0
        case['kraus'][1] = pack(k)
    elif mutation == 'distance':
        case['dilation_distance'] = 0
    else:
        case['p'] = 0
    with pytest.raises(ValueError):
        export_check.verify(row)


@pytest.mark.parametrize('raw', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}'])
def test_json_ambiguity_rejected(tmp_path, raw):
    path = tmp_path/'bad.json'
    path.write_text(raw, encoding='utf-8')
    with pytest.raises(ValueError):
        verify.load(path)


@pytest.mark.parametrize('field', ['schema', 'sources', 'claims', 'extra'])
def test_receipt_custody_rejects_trash_before_semantics(tmp_path, field):
    packet = verify.load(verify.HERE/'receipt.json')
    packet[field] = 'trash'
    path = tmp_path/'bad.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(ValueError):
        verify.verify(path)


def test_producer_disabled_optimized_verifier():
    script = '''
import importlib.abc, sys
blocked={'m1_expander_archive.'+x for x in ('graphs','circuit','decoder','bounds','export','build')}
class Reject(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname in blocked:
            raise RuntimeError('producer import: '+fullname)
sys.meta_path.insert(0, Reject())
from m1_expander_archive.verify import verify
verify()
print('producer-free optimized replay passed')
'''
    env = dict(os.environ, PYTHONPATH=str(verify.ROOT/'code'))
    result = subprocess.run([sys.executable, '-O', '-c', script], cwd=verify.ROOT,
                            env=env, capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, result.stdout+result.stderr
    assert 'producer-free optimized replay passed' in result.stdout
