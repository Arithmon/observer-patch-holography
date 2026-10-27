"""Fixed native feedback records, blind inversion and complement accounting."""
from copy import deepcopy
from unittest.mock import patch

import numpy as np
import pytest

from . import native_repair_mass_readout as h


@pytest.fixture(scope='module')
def execution():
    return h.execute()


def test_fixed_detectors_are_complete_and_keep_four_dimensional_complement():
    for f in range(3):
        for g in range(3):
            effects = h.fixed_effects(f, g)
            np.testing.assert_allclose(effects.sum(axis=0), np.eye(6), atol=3e-16, rtol=0)
            assert [round(np.trace(e).real) for e in effects] == [1, 1, 4]
            assert min(np.linalg.eigvalsh(e).min() for e in effects) >= -1e-16


def test_source_forward_path_never_reads_the_averaged_kernel_or_its_spectrum(execution):
    architecture = {k: v for k, v in h.source.architecture().items()
                    if k not in ('Y', 'e2', 'determinant')}
    with patch.object(h.source, 'architecture', return_value=architecture):
        _, arrays = h.execute()
    np.testing.assert_array_equal(arrays['readout_probabilities'], execution[1]['readout_probabilities'])


def test_blind_inverse_contains_source_kernel_with_explicit_finite_channel_error(execution):
    report, arrays = execution
    bounds = report['error_budget']
    with patch.object(h.source, 'architecture', side_effect=AssertionError('inverse accessed source')):
        decoded = h.decode_records(arrays['readout_probabilities'], time=h.TIME,
                                   entry_error_upper=bounds['entry_error_upper'])
    actual = np.asarray(h.source.architecture()['Y'], float)
    assert np.max(abs(decoded['estimate']-actual)) <= bounds['entry_error_upper']
    for value, (lo, hi) in zip(np.linalg.eigvalsh(actual), decoded['eigenvalue_intervals'], strict=True):
        assert lo <= value <= hi
    assert decoded['strictly_positive_kernel_resolved']
    assert decoded['particle_names_or_calibrated_mass_inputs'] is None


def test_final_memory_and_all_pointer_outcomes_are_unconditionally_retained(execution):
    report, arrays = execution
    assert arrays['readout_probabilities'].shape == (9, 5, 3)
    np.testing.assert_allclose(arrays['readout_probabilities'].sum(axis=(1, 2)), np.ones(9), atol=1e-13)
    assert np.any(arrays['readout_probabilities'][:, :, 2] > 0)
    assert report['error_budget']['stationary_product_initialization_required']
    assert not report['scope']['forward_evolution_uses_averaged_kernel_exponential']


@pytest.mark.parametrize('change', ['missing', 'nonfinite', 'extra', 'probability', 'estimate'])
def test_retained_readout_arrays_fail_closed(execution, change):
    arrays = deepcopy(execution[1])
    if change == 'missing': arrays.pop('readout_terminal_blocks')
    elif change == 'nonfinite': arrays['readout_terminal_blocks'][0, 0, 0, 0] = np.nan
    elif change == 'extra': arrays['untracked'] = np.zeros(1)
    elif change == 'probability': arrays['readout_probabilities'][0, 0, 0] += 1e-5
    else: arrays['readout_estimate'][0, 0] += 1e-5
    with pytest.raises(ValueError): h.validate_arrays(arrays)


def test_inverse_never_renormalizes_a_lost_terminal_record(execution):
    probabilities = execution[1]['readout_probabilities'].copy()
    probabilities[0, 0] = 0
    with pytest.raises(ValueError, match='complete'):
        h.decode_records(probabilities, time=h.TIME, entry_error_upper=1e-5)
