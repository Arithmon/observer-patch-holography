"""The large-source producer must survive its serialized public replay."""
import json
from fractions import Fraction as F

import pytest

from ..observer.scalar_response.experiment import canonical, compare_receipt
from .receipt import run


def test_source_screen_receipt_roundtrip_and_complete_budgets():
    actual = run()
    serialized = json.loads(canonical(actual))
    # Catch non-JSON tuple fields, which otherwise fail only in replay CI.
    compare_receipt(serialized, actual)
    assert actual['geometry']['mask_counts'] == [23, 5, 4, 1]
    assert actual['contract']['population_dual_dimension'] == 24754
    assert len(actual['shared_population_blocks']) == 21
    assert F(actual['budgets']['two_completed_screen_allowance']) < F(2731, 10**13)
    for channel in actual['channels']:
        assert len(channel['coherences']) == 21
        assert F(channel['full_screen_to_export_midpoint_allowance']) < F(7779, 10**13)
    witness = actual['quantum_witnesses']
    assert witness['completed_bell_C_E_distance'] == '101/968'
    assert F(witness['C_word_1_0_witness_violation_lower']) > F(26, 1000)
    assert witness['word_1_0_divided_by_its_probability'] is False
    for control in actual['runtime_controls']['executions']:
        assert control['checkpoint_replayed']
        assert control['full_output_distance_numeric'] < 2e-12
    # Rational budgets are exact receipt fields, not tolerance-based floats.
    serialized['budgets']['two_completed_screen_allowance'] = '0'
    with pytest.raises(ValueError):
        compare_receipt(serialized, actual)
