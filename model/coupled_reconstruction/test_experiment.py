"""Fresh whole-field calibration, restricted inverse access, and transfer tests."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest

from ..interacting_reconstruction import source as isolated_source
from ..interacting_reconstruction.prediction import RecoveredField
from . import experiment, recovery, source


@pytest.fixture(scope='module')
def executed():
    events, records, inputs = [], [], []
    original_calibration, original_recover = source.calibration_records, recovery.recover
    original_history, original_audit = RecoveredField.history, experiment.target_audit
    original_cap = experiment.bounds.product_generator_cap_certificate

    def forbidden(*args, **kwargs):
        raise AssertionError('isolated probe or hidden backend queried during identification')

    def calibration(*args, **kwargs):
        events.append('whole_field_calibration')
        value = original_calibration(*args, **kwargs)
        records.append(deepcopy(value))
        return value

    def infer(*args, **kwargs):
        events.append('record_only_inverse')
        inputs.append(deepcopy(args))
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(experiment, 'target_audit', forbidden)
            patch.setattr(experiment.bounds, 'product_generator_cap_certificate', forbidden)
            patch.setattr(isolated_source, 'supplied_target', forbidden)
            patch.setattr(isolated_source.c, 'execute', forbidden)
            return original_recover(*args, **kwargs)

    def predict(field, *args, **kwargs):
        events.append('prediction')
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(experiment, 'target_audit', forbidden)
            patch.setattr(source, 'calibration_records', forbidden)
            patch.setattr(isolated_source.c, 'execute', forbidden)
            return original_history(field, *args, **kwargs)

    def audit():
        events.append('target_audit')
        return original_audit()

    def cap(*args, **kwargs):
        events.append('cap_audit')
        return original_cap(*args, **kwargs)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(source, 'calibration_records', calibration)
        patch.setattr(recovery, 'recover', infer)
        patch.setattr(RecoveredField, 'history', predict)
        patch.setattr(experiment, 'target_audit', audit)
        patch.setattr(experiment.bounds, 'product_generator_cap_certificate', cap)
        patch.setattr(isolated_source, 'onsite_probe', forbidden)
        patch.setattr(isolated_source, 'bond_probe', forbidden)
        result = experiment.run()
    return result, events, records, inputs


def canonical(value):
    return json.dumps(experiment.serializable(value), sort_keys=True, allow_nan=False)


def test_only_whole_field_records_and_clock_laws_reach_inverse(executed):
    result, events, records, inputs = executed
    assert events[:3] == ['whole_field_calibration']*2+['record_only_inverse']
    assert events[-2:] == ['target_audit', 'cap_audit']
    assert events.count('prediction') == 8
    assert len(inputs) == 1 and len(inputs[0]) == 4
    assert canonical(inputs[0][:2]) == canonical(records)
    assert canonical(records) == canonical(result['calibration_records'])
    for actual, retained in zip(inputs[0][2:], result['calibration_clock_probabilities'], strict=True):
        np.testing.assert_array_equal(actual, np.asarray(retained)[:, 0])


def test_recovered_interaction_graph_and_all_coupled_histories_pass(executed):
    result = executed[0]
    inferred = result['recovered']
    assert inferred['bare_mass_squared'] == pytest.approx(1., abs=2e-7)
    assert inferred['bond_kappa'] == pytest.approx(1., abs=2e-7)
    assert inferred['quartic'] == pytest.approx(1/24, abs=2e-7)
    assert inferred['homogeneous_intervals']['quartic']['lower'] > 0
    assert sum(map(sum, inferred['adjacency'])) == 24
    assert inferred['spatial_dimension'] == 3
    for row in result['levels']:
        assert row['all_pass']
        assert max(row['prediction_history_tv_by_phase']) < 2e-8
        assert min(row['quadratic_deletion_history_tv_by_phase']) > row['ideal_inferred_to_source_history_tv_upper']
        assert row['native_history_count'] == 1056
        assert abs(sum(x['probability'] for x in row['native_executed_history'])-1) < 2e-9


def test_schedule_certificate_and_source_bindings_remain_explicit(executed):
    result = executed[0]
    cap = result['post_prediction_backend_audit']['product_generator_cap_certificate']
    assert cap['declared_cap_covers_source_schedule']
    assert cap['schedule_centered_norm_sum_upper'] < experiment.PRODUCT_NORM_CAP
    protocol = result['protocol']
    assert protocol['every_onsite_and_bond_term_participates']
    assert not protocol['generator_isolation_required']
    assert not protocol['product_family_is_exact_continuous_semigroup']
    assert not protocol['separate_calibration_splitting_error_needed']
    assert not protocol['inference_arithmetic_compiled_into_source']
    assert not result['floating_clock_and_physical_service_error_certified']
    assert all(result['checks'].values()) and result['all_pass']
    for path, digest in result['source_binding_sha256'].items():
        assert digest == hashlib.sha256((experiment.ROOT/path).read_bytes()).hexdigest()
