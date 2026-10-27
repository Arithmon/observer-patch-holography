"""Bounded receipt replay must reject valid replacement experiments too."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest

from . import native_family_protocol as protocol


@pytest.fixture(scope='module')
def execution():
    return protocol.build()


def test_complete_array_schema_and_process_boundaries(execution):
    report, arrays = execution
    protocol.validate_arrays(report, arrays)
    assert arrays['memory_probabilities'].shape == (2, 9, 1536, 3)
    assert arrays['field_spectrum'].shape == (180,)
    assert len(report['interaction_records']['all_catalog_entries']) == 3
    assert len(report['interaction_records']['duty_catalog']) == 4
    np.testing.assert_array_equal(arrays['duty_ranks'][2:], np.full((2, 3), 2))
    assert not any(row['pulse_duty_renormalized'] for row in report['interaction_records']['duty_catalog'])
    assert not report['record_processes']['static_mean_service_identified_with_finite_CQ_memory']
    assert not report['scope']['dimensionless_coefficients_are_physical_pole_masses']
    assert report['filtered_current']['exact']['logical_current_dimension'] == 12
    assert report['filtered_current']['theorem_interfaces']['Higgs_tensor_substitution_checks'] == 72
    assert len(report['filtered_current']['exact']['declared_from_filtered']) == 12


@pytest.mark.parametrize('change', ['missing_zero_labels', 'shape', 'nonfinite', 'extra'])
def test_complete_artifact_schema_fails_closed(execution, change):
    report, arrays = execution
    altered = dict(arrays)
    if change == 'missing_zero_labels':
        altered['memory_labels'] = arrays['memory_labels'][1:]
    elif change == 'shape':
        altered['field_spectrum'] = arrays['field_spectrum'][:-1]
    elif change == 'nonfinite':
        altered['catalog_coefficients'] = arrays['catalog_coefficients'].copy()
        altered['catalog_coefficients'].flat[0] = np.nan
    else:
        altered['untracked'] = np.zeros(1)
    with pytest.raises(ValueError):
        protocol.validate_arrays(report, altered)


def test_static_records_from_another_valid_coefficient_fail_source_binding(execution):
    report, arrays = execution
    replacement = protocol.field.forward_records(np.diag([.1, .2, .7]))
    protocol.field.decode_records(replacement['probabilities'], replacement['survival'])
    altered = dict(arrays)
    altered['static_probabilities'] = arrays['static_probabilities'].copy()
    altered['static_survival'] = arrays['static_survival'].copy()
    altered['static_probabilities'][0] = replacement['probabilities']
    altered['static_survival'][0] = replacement['survival']
    with pytest.raises(ValueError, match='differ from'):
        protocol.validate_arrays(report, altered)


def _retained_pair(tmp_path, report, arrays):
    path = tmp_path/'receipt.json'
    artifact = path.with_suffix('.npz')
    np.savez_compressed(artifact, **arrays)
    retained = deepcopy(report)
    retained.update(source_bindings=protocol.source_bindings(),
                    artifact={'path': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()})
    path.write_text(json.dumps(retained, allow_nan=False), encoding='utf-8')
    return path


def test_retained_receipt_replays_every_complete_array(execution, tmp_path):
    report, arrays = execution
    assert protocol.check(_retained_pair(tmp_path, report, arrays))


def test_rehashing_and_decoding_a_valid_replacement_does_not_admit_it(execution, tmp_path):
    report, arrays = execution
    altered = {key: value.copy() for key, value in arrays.items()}
    # Another complete, valid current instrument is still the wrong protocol.
    altered['current_kraus'] = protocol.gates.instrument(0)['kraus']
    assert protocol.gates.validate_kraus(altered['current_kraus']) == 0
    path = _retained_pair(tmp_path, report, altered)
    with pytest.raises(ValueError, match='fresh source'):
        protocol.check(path)


def test_changed_exact_claim_and_physical_promotion_are_rejected(execution, tmp_path):
    report, arrays = execution
    for changed in ('rank', 'mass'):
        altered = deepcopy(report)
        if changed == 'rank':
            altered['filtered_current']['exact']['logical_current_dimension'] = 13
        else:
            altered['scope']['dimensionless_coefficients_are_physical_pole_masses'] = True
        with pytest.raises(ValueError, match='exact claim'):
            protocol.check(_retained_pair(tmp_path, altered, arrays))


def test_fresh_sources_bind_all_four_producers_and_receipt():
    sources = protocol.source_bindings()
    for name in ('native_filtered_current_control', 'native_family_record_memory',
                 'native_gauge_family_records', 'native_record_yukawa_field', 'native_family_protocol'):
        assert 'model/common_source/'+name+'.py' in sources


def test_diagnostic_ratios_allow_conditioned_roundoff_without_relaxing_exact_claims():
    protocol._compare(1179523.1449, 1179523.14487, 'receipt.family_eigenvalue_ratios[2]')
    with pytest.raises(ValueError, match='numerical value'):
        protocol._compare(1179524.1449, 1179523.14487, 'receipt.family_eigenvalue_ratios[2]')
    with pytest.raises(ValueError, match='exact claim'):
        protocol._compare(True, False, 'receipt.physical_mass_claim')
    with pytest.raises(ValueError, match='original type'):
        protocol._compare(12., 12, 'receipt.current_dimension')


def test_canonical_verifier_discovers_the_retained_protocol():
    from model import verify
    assert ('model.common_source.native_family_protocol', None,
            'model/common_source/native_family_protocol_results.json') in verify.REPLAYS


def test_receipt_duplicate_json_keys_fail_closed(execution, tmp_path):
    path = _retained_pair(tmp_path, *execution)
    path.write_text('{"schema": "first", '+path.read_text()[1:], encoding='utf-8')
    with pytest.raises(ValueError, match='duplicate'):
        protocol.check(path)
