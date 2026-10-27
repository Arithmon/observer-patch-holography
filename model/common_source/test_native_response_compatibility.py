"""Receipt custody and strict-claim controls for the composed current result."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest

from . import native_response_compatibility as packet


@pytest.fixture(scope='module')
def execution():
    return packet.execute()


def test_real_execution_roundtrips_original_json_types(execution):
    report, arrays = execution
    restored = json.loads(json.dumps(report, allow_nan=False))
    packet.compare_execution(restored, arrays, report, arrays)


@pytest.mark.parametrize('mutation', ['empty', 'missing', 'promote', 'dimension'])
def test_exact_claim_mutations_fail(execution, mutation):
    report, arrays = execution
    bad = deepcopy(report)
    if mutation == 'empty':
        bad = {}
    elif mutation == 'missing':
        del bad['controller']
    elif mutation == 'promote':
        bad['comparison']['physical_standard_model_derived'] = True
    else:
        bad['reduction'] = {'dimension': 12}
    with pytest.raises(ValueError):
        packet.compare_execution(bad, arrays, report, arrays)


@pytest.mark.parametrize('bad', [True, 1.0, '1', None])
def test_exact_scalar_aliases_fail(bad):
    with pytest.raises(ValueError):
        packet._exact(1, bad)


def test_checker_owns_fresh_replay_and_preserves_existing_outputs(tmp_path, monkeypatch, execution):
    calls = []
    def replay():
        calls.append(True)
        return execution
    monkeypatch.setattr(packet, 'execute', replay)
    path = tmp_path / 'response.json'
    packet.write(path)
    packet.check(path)
    assert len(calls) == 2
    original = path.read_bytes(), path.with_suffix('.npz').read_bytes()
    with pytest.raises(ValueError, match='fresh'):
        packet.write(path)
    assert original == (path.read_bytes(), path.with_suffix('.npz').read_bytes())
    assert len(calls) == 2


def test_forged_binding_and_duplicate_key_rejected_before_replay(tmp_path, monkeypatch):
    def forbidden():
        raise AssertionError('invalid custody reached source execution')
    monkeypatch.setattr(packet, 'execute', forbidden)
    path = tmp_path / 'bad.json'
    path.write_text('{"source_bindings":{}}', encoding='utf-8')
    with pytest.raises(ValueError, match='bindings'):
        packet.check(path)
    path.write_text('{"schema":1,"schema":1}', encoding='utf-8')
    with pytest.raises(ValueError):
        packet.check(path)


@pytest.mark.parametrize('mutation', ['reference_unitary', 'controller_unitary', 'controller_coefficient', 'cooling_unitary'])
def test_rehashed_wrong_source_artifact_still_fails(tmp_path, monkeypatch, execution, mutation):
    monkeypatch.setattr(packet, 'execute', lambda: execution)
    path = tmp_path / 'response.json'
    packet.write(path)
    artifact = path.with_suffix('.npz')
    with np.load(artifact, allow_pickle=False) as stored:
        arrays = {key: stored[key] for key in stored.files}
    if mutation == 'reference_unitary':
        arrays['reference__controls'][15] = np.eye(8)
    elif mutation == 'controller_unitary':
        # A global phase preserves every Gram matrix but violates the actual
        # complete emitted-column action required by the source receipt.
        arrays['controller__controls'][0, 0] *= 1j
    elif mutation == 'controller_coefficient':
        coefficients = arrays['controller__coefficients']
        index = tuple(np.argwhere(coefficients != 0)[0])
        coefficients[index] *= -1
    else:
        arrays['cooling__replaced'][:, 7] *= 1j
    with artifact.open('wb') as output:
        np.savez_compressed(output, **arrays)
    report = json.loads(path.read_text(encoding='utf-8'))
    report['artifact']['sha256'] = hashlib.sha256(artifact.read_bytes()).hexdigest()
    path.write_text(json.dumps(report, allow_nan=False), encoding='utf-8')
    with pytest.raises(ValueError):
        packet.check(path)


def test_all_components_and_array_owners_share_the_retained_sources(execution):
    report, arrays = execution
    assert report['schema'] == 'oph.native_response_compatibility.v2'
    assert set(report) == {'schema', 'comparison', 'reduction', 'biased', 'dynamics',
                           'reference', 'cooling', 'controller', 'source_composition'}
    assert report['dynamics']['bias'] == '16'
    split = packet.split_arrays(arrays)
    assert set(split) == {'reference', 'cooling', 'controller'}
    assert split['reference']['controls'] is arrays['reference__controls']
    assert split['cooling']['counter'] is arrays['cooling__counter']
    assert split['controller']['controls'].shape == (17, 18, 8, 8)
    assert split['controller']['coefficients'].shape == (17, 18, 18)
    assert not report['source_composition']['independent_surrogate_clock_samples_used']


@pytest.mark.parametrize('mutation', ['unprefixed', 'unknown_owner', 'missing_owner', 'ambiguous_name'])
def test_array_namespaces_fail_closed_before_owner_validation(execution, mutation):
    _, arrays = execution
    bad = dict(arrays)
    if mutation == 'unprefixed':
        bad['controls'] = bad.pop('reference__controls')
    elif mutation == 'unknown_owner':
        bad['foreign__controls'] = bad.pop('reference__controls')
    elif mutation == 'missing_owner':
        bad = {key: value for key, value in bad.items() if not key.startswith('cooling__')}
    else:
        bad['reference__nested__controls'] = bad.pop('reference__controls')
    with pytest.raises(ValueError, match='owner'):
        packet.split_arrays(bad)


@pytest.mark.parametrize('value', [0., 1e-20, 1e-5, True])
def test_small_fixed_admissions_cannot_use_absolute_replay_slack(value):
    with pytest.raises(ValueError):
        packet._fixed_numbers({'nested': [2e-11]}, {'nested': [value]}, 'budget')


def test_controller_compare_receives_the_same_validated_reference_and_clock(execution, monkeypatch):
    report, arrays = execution
    called = []
    controller = packet._controller()
    original = controller.compare_execution

    def checked(*args, **kwargs):
        assert kwargs['reference_arrays']['native_C'] is arrays['reference__native_C']
        assert kwargs['cooling_arrays']['counter'] is arrays['cooling__counter']
        assert kwargs['fresh_reference_arrays']['native_C'] is arrays['reference__native_C']
        assert kwargs['fresh_cooling_arrays']['counter'] is arrays['cooling__counter']
        called.append(True)
        return original(*args, **kwargs)

    monkeypatch.setattr(controller, 'compare_execution', checked)
    packet.compare_execution(report, arrays, report, arrays)
    assert called == [True]


def test_exact_biased_weight_diagnostic_cannot_be_forged_to_zero(execution):
    report, arrays = execution
    bad = deepcopy(report)
    bad['biased']['examples'][1]['first_triplet_weight_diagnostic'] = 0.
    with pytest.raises(ValueError, match='biased'):
        packet.compare_execution(bad, arrays, report, arrays)
