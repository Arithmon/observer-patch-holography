from fractions import Fraction as F

import numpy as np
import pytest

from . import bounds, protocol
from .lattice import GalerkinParameters, SevenModeLattice


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('feedback', [False, True])
def test_complete_source_history_matches_independent_continuous_reference(sign, feedback):
    # Dimension two has commuting evolution, so this also tests the exact-zero
    # temporal bound instead of allowing a generous error budget to hide bugs.
    lattice = SevenModeLattice((2,)*7, GalerkinParameters(grid=5))
    actual = protocol.history(lattice, steps=2, sign=sign, feedback=feedback)
    reference = protocol.history(lattice, source_execution=False, sign=sign, feedback=feedback)
    assert [row['bits'] for row in actual] == [[0, 0], [0, 1], [1, 0], [1, 1]]
    assert abs(sum(row['probability'] for row in actual)-1) < 2e-12
    assert protocol.tv(actual, reference) < 2e-12
    assert bounds.splitting_history_bound(2, (F(1, 8),)*2, 2) == 0


def test_interacting_noncommuting_history_and_spatial_map():
    lattice = SevenModeLattice((4,)*7, GalerkinParameters(grid=5))
    reference = protocol.history(lattice, source_execution=False)
    actual = protocol.history(lattice, steps=8)
    assert protocol.tv(actual, reference) < float(bounds.splitting_history_bound(4, (F(1, 8),)*2, 8))
    exact_gradient = protocol.history(SevenModeLattice((4,)*7), source_execution=False)
    assert protocol.tv(reference, exact_gradient) < float(bounds.spatial_history_bound(4, 5, F(1, 4)))


def test_feedback_changes_future_but_not_the_already_committed_record():
    lattice = SevenModeLattice((4,)*7)
    connected = protocol.history(lattice, source_execution=False)
    no_feedback = protocol.history(lattice, source_execution=False, feedback=False)
    reset = protocol.history(lattice, source_execution=False, reset=True)
    assert protocol.tv(connected, no_feedback) > .001
    assert protocol.tv(connected, reset) > .01
    for bit in (0, 1):
        first = lambda law: sum(row['probability'] for row in law if row['bits'][0] == bit)
        assert abs(first(connected)-first(no_feedback)) < 2e-13
        assert abs(first(connected)-first(reset)) < 2e-13


def test_native_records_have_the_complete_joint_alphabet_and_feedback():
    lattice = SevenModeLattice((2,)*7)
    laws = [protocol.history(lattice, steps=1, sign=sign) for sign in (-1, 1)]
    actual = protocol.native_records(laws)
    detached = protocol.native_records(laws, feedback=False)
    assert len(actual) == 1056
    assert len({(row['initial'], *row['records']) for row in actual}) == 1056
    assert abs(sum(row['probability'] for row in actual)-1) < 2e-12
    assert any(a['final_port_zero'] != b['final_port_zero'] for a, b in zip(actual, detached))
    assert [row['probability'] for row in actual] == [row['probability'] for row in detached]


def test_resource_ledger_counts_emitted_factors_without_runtime_shortcuts():
    lattice = SevenModeLattice((2,)*7)
    one = protocol.source_cost(lattice, 1)
    two = protocol.source_cost(lattice, 2)
    assert 0 < one['native_pulses'] < two['native_pulses'] < 2*one['native_pulses']
    assert two['proper_code_service_events'] > one['proper_code_service_events']


def test_blank_calibration_does_not_refit_or_certify_later_records():
    report = protocol.detector_calibration(SevenModeLattice((4,)*7))
    assert report['fitted_parameters'] == []
    assert report['experimental_data_used'] is False
    assert all(0 < row['absolute_occupation_difference'] < .001 for row in report['settings'])


def test_strong_pointer_exposes_a_history_cutoff_problem():
    settings = protocol.Settings(eta=2., initial_kick=1., feedback_kick=1.)
    four = protocol.history(SevenModeLattice((4,)*7), settings=settings, source_execution=False)
    wider = protocol.history(SevenModeLattice((4, 8, 4, 4, 4, 4, 4)), settings=settings, source_execution=False)
    # A tiny blank-state calibration error cannot stand in for the complete
    # history's error after preparation, disturbance, feedback and evolution.
    assert protocol.tv(four, wider) > .01


@pytest.mark.parametrize('arguments', [dict(interval=-1), dict(interval=.1),
    dict(interval=True), dict(mode=True), dict(mode=7), dict(eta=float('nan')),
    dict(theta=float('inf')), dict(initial_kick=True)])
def test_invalid_observation_contract_is_rejected(arguments):
    with pytest.raises(ValueError):
        protocol.Settings(**arguments)


def test_record_alphabet_mismatch_is_rejected():
    with pytest.raises(ValueError):
        protocol.tv([{'bits': [0, 1], 'probability': 1}], [{'bits': [1, 0], 'probability': 1}])
