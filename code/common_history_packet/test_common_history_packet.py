from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

import pytest

spec = importlib.util.spec_from_file_location('common_history_packet_verifier', Path(__file__).with_name('verify.py'))
verifier = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = verifier
spec.loader.exec_module(verifier)


@pytest.fixture(scope='module')
def context():
    # Full parent replay is exercised by the standalone scientific verifier.
    # Mutation tests use the independently assembled finite model directly.
    return verifier.context(full_parent=False)


@pytest.fixture(scope='module')
def packet():
    return verifier.load()


def reseal(trace):
    chain = verifier.digest({'model': trace['model_sha256'], 'case': trace['case']})
    for event in trace['events']:
        event['parent'] = chain
        event['hash'] = verifier.digest({k: v for k, v in event.items() if k != 'hash'})
        chain = event['hash']
    trace['final_hash'] = chain


@pytest.mark.parametrize('mutation', ['writer', 'value', 'missing_route', 'disabled_feedback', 'foreign_record'])
def test_resealed_event_forgeries_fail(packet, context, mutation):
    trace = deepcopy(packet['cases'][0])
    if mutation == 'missing_route':
        trace['events'].pop(next(i for i, e in enumerate(trace['events']) if e['op'] == 'route'))
    elif mutation == 'disabled_feedback':
        event = next(e for e in trace['events'] if e['op'] == 'feedback')
        event['op'] = 'feedback_disabled'
        del event['writes'][next(k for k in event['writes'] if k.startswith('u/'))]
    elif mutation == 'foreign_record':
        event = next(e for e in trace['events'] if e['op'] == 'route')
        entry = next(iter(event['reads'].values()))
        entry['value'] = ['1/1000000', '0']
    else:
        event = next(e for e in trace['events'] if e['op'] == 'evolve')
        if mutation == 'writer':
            # Another real writer with the same zero value is still wrong.
            key = next(k for k in event['reads'] if k.startswith('u/'))
            event['reads'][key]['writer'] = 0
        else:
            event['writes'][next(iter(event['writes']))] = ['1', '0']
    reseal(trace)
    with pytest.raises(ValueError):
        verifier.replay(trace, context)


@pytest.mark.parametrize('mutation', ['clock', 'action', 'readout', 'population', 'error_budget', 'quantum_promotion'])
def test_resealed_scientific_substitutions_fail(packet, context, mutation):
    changed = dict(packet)
    if mutation in ('error_budget', 'quantum_promotion'):
        changed['contract'] = dict(packet['contract'])
        changed['contract']['record_error' if mutation == 'error_budget' else 'quantum_claim'] = '0' if mutation == 'error_budget' else True
    else:
        changed['model'] = deepcopy(packet['model'])
        if mutation == 'clock':
            changed['model']['tau'] = ['-1/70', '1/35']
        elif mutation == 'action':
            changed['model']['action'][0][0][1][0] = str(F(changed['model']['action'][0][0][1][0])+1)
        elif mutation == 'readout':
            changed['model']['detector'][32] = ['0', '0']
        else:
            changed['model']['addresses'][0]['coordinate_Qphi'][0] = ['0', '0']
    with pytest.raises(ValueError):
        verifier.verify(changed, context)


def test_positive_receipt_has_useful_full_window_and_real_feedback_consumption():
    root = Path(__file__).resolve().parents[2]
    verifier.check_manifest()
    receipt = json.loads((root/'evidence/common_history_packet_20260925/verification.json').read_text())
    assert receipt['packet_sha256'] == sha256((root/'evidence/common_history_packet_20260925/packet.json').read_bytes()).hexdigest()
    assert receipt['verdict'] == 'PASS' and receipt['parent_full_independent_replay']
    assert 0 < F(receipt['comparison_total_error_upper']) < F(receipt['paired_response_lower'])
    assert F(receipt['strict_margin']) > F(1, 1000)
    for name, resources in receipt['resources'].items():
        assert resources['operation_counts']['route'] == 17*32
        if name != 'disabled_feedback':
            assert resources['reconstruction_mismatches'] == 0
            assert resources['restored_writer_reads_in_future_action'] > 0
    assert receipt['resources']['disabled_feedback']['reconstruction_mismatches'] > 0


def test_error_bound_cannot_be_made_small_by_discarding_control_histories(packet, context):
    changed = dict(packet)
    changed['cases'] = packet['cases'][:3]
    with pytest.raises(ValueError, match='complete control cases'):
        verifier.verify(changed, context)


def test_closed_custody_rejects_unlisted_nested_payload(tmp_path):
    folder = tmp_path/'evidence/common_history_packet_20260925'
    folder.mkdir(parents=True)
    for name in ('packet.json', 'verification.json', 'manifest.json'):
        (folder/name).write_text('{}')
    (folder/'extra').mkdir()
    (folder/'extra/unlisted.json').write_text('{}')
    with pytest.raises(ValueError, match='closed evidence file inventory'):
        verifier.check_manifest(tmp_path)
