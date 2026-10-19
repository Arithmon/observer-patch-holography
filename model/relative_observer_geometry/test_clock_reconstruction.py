from copy import deepcopy
from fractions import Fraction as F

import pytest

from . import clock_reconstruction as clock


@pytest.fixture(scope='module')
def schedules():
    seams = tuple(sorted(tuple(sorted(e)) for e in clock.event._source().seams))
    return ({e: F(1, 30) for e in seams},
            {e: F(i+1, 465) for i, e in enumerate(seams)})


def test_nonuniform_scheduler_is_identified_from_clock_records_only(schedules):
    uniform, biased = schedules
    for weights in schedules:
        recovered, exposure = clock.recover(clock.clock_records(weights))
        assert recovered == weights
        assert sum(exposure) == 2
    assert clock.structural_gram(clock.clock_records(uniform)) == clock.structural_gram(clock.clock_records(biased))
    assert len(set(clock.recover(clock.clock_records(biased))[1])) > 1


@pytest.mark.parametrize('z,ticks', [(F(1, 2), 1), (F(1), 2)])
def test_inferred_weights_predict_full_native_stopped_joint_law(schedules, z, ticks):
    for weights in schedules:
        records = clock.clock_records(weights)
        predicted = clock.predict_field(records, 0, 1, z, ticks)
        assert predicted == clock.native_field(weights, 0, 1, z, ticks)
        assert sum(predicted.values()) == clock.renewal_pgf(records, 0, 1, z, ticks)


def test_same_structural_geometry_does_not_license_uniform_kinetics(schedules):
    uniform, biased = (clock.clock_records(w) for w in schedules)
    assert clock.predict_field(uniform, 0, 1) != clock.predict_field(biased, 0, 1)
    assert clock.predict_field(uniform, 0, 1) == clock.event.predicted_field(0, 1)


@pytest.mark.parametrize('mutation', ['duplicate', 'missing', 'hidden', 'mean', 'coincidence', 'noncanonical', 'bool'])
def test_inconsistent_or_undeclared_pair_evidence_is_rejected(schedules, mutation):
    records = deepcopy(clock.clock_records(schedules[1]))
    if mutation == 'duplicate':
        records[-1] = deepcopy(records[0])
    elif mutation == 'missing':
        records.pop()
    elif mutation == 'hidden':
        records[0]['source_seam'] = [0, 1]
    elif mutation == 'mean':
        records[0]['mean_B_per_A'] = '1/999'
    elif mutation == 'coincidence':
        records[0]['terminal_coincidence'] = '1/999'
    elif mutation == 'noncanonical':
        records[0]['terminal_coincidence'] = '0/2'
    else:
        records[0]['ports'][0] = False
    with pytest.raises(ValueError):
        clock.recover(records)


def test_disconnected_positive_clock_support_does_not_certify_native_geometry():
    weights = {(i, i+1): F(1, 6) for i in range(0, 12, 2)}
    with pytest.raises(ValueError, match='disconnected'):
        clock.recover(clock.clock_records(weights))


def test_other_connected_graph_can_be_identified_without_becoming_native_rank_three():
    weights = {(i, i+1): F(1, 11) for i in range(11)}
    records = clock.clock_records(weights)
    assert clock.recover(records)[0] == weights
    with pytest.raises(ValueError, match='rank-three'):
        clock.structural_gram(records)


def test_port_relabelling_transports_the_inferred_scheduler(schedules):
    permutation = tuple(reversed(range(12)))
    records = clock.clock_records(schedules[1])
    for row in records:
        row['ports'] = [permutation[x] for x in row['ports']]
    recovered, _ = clock.recover(records)
    assert recovered == {tuple(sorted(permutation[x] for x in e)): p
                         for e, p in schedules[1].items()}


def test_zero_exposure_and_nonprobability_schedules_are_rejected():
    with pytest.raises(ValueError, match='positive probability'):
        clock.clock_records({(0, 1): F(1)})
    with pytest.raises(ValueError, match='sum to one'):
        clock.clock_records({(0, 1): F(1, 2)})


def test_temporal_dependence_is_not_identified_by_pair_means_and_coincidences(schedules):
    uniform = schedules[0]
    records = clock.clock_records(uniform)
    for row in records:
        control = clock.periodic_control(tuple(uniform), *row['ports'])
        assert str(control['mean_B_per_A']) == row['mean_B_per_A']
        assert str(control['terminal_coincidence']) == row['terminal_coincidence']
    control = clock.periodic_control(tuple(uniform))
    assert control['zero_B_count_probability'] == F(4, 5)
    assert clock.renewal_pgf(records, 0, 1, F(0)) == F(4, 9)


def test_oracle_cache_does_not_expose_mutable_fraction_objects(schedules):
    result = clock.native_field(schedules[0], 0, 1)
    result[1, 1]._numerator = 999
    assert clock.native_field(schedules[0], 0, 1)[1, 1] == F(14, 407)
