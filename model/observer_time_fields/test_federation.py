"""Independent source histories, full-sector error budgets and hostile inputs."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from ..federated_observer import source as native
from . import federation as f


def _source_step(population, seams, denominator):
    result = {}
    for occupied, weight in population.items():
        for seam in seams:
            for coin in (0, 1):
                target, _ = native.binary_step(occupied, seam, coin)
                result[target] = result.get(target, F(0))+weight/denominator
    return result


def _direct_one_read(initial, cap, feedback):
    data = f.source()
    label = tuple(int(port in initial) for port in data.read_ports)
    waiting, stopped = {initial: F(1)}, {}
    for _ in range(cap):
        for occupied, value in _source_step(waiting, data.incident, data.denominator).items():
            stopped[occupied] = stopped.get(occupied, F(0))+value
        waiting = _source_step(waiting, data.remote, data.denominator)
    result = {record: {} for record in f.records(1)}
    result[label, (), 'timeout'] = waiting
    for bit in (0, 1):
        selected = {occupied: value for occupied, value in stopped.items()
                    if int(data.read_ports[0] in occupied) == bit}
        if feedback and bit:
            selected = _source_step(selected, (data.read_ports,), 2)
        result[label, (bit,), 'complete'] = selected
    return result


@pytest.mark.parametrize('feedback', (False, True))
def test_every_output_density_is_below_independent_complete_native_history(feedback):
    data = f.source()
    initial = tuple(sorted((data.read_ports[0], data.remote[0][0])))
    expected = _direct_one_read(initial, 2, feedback)
    actual = f.execute(cap=2, ticks=1, grid_bits=40, feedback=feedback, initial=initial)
    scale = 1 << actual.grid_bits
    exact_mass = F(0)
    index = {occupied: i for i, occupied in enumerate(data.states)}
    for column, record in enumerate(actual.records):
        targets = expected[record]
        exact_mass += sum(targets.values(), F(0))
        for row in np.flatnonzero(actual.counts[:, column]):
            assert F(int(actual.counts[row, column]), scale) <= targets.get(data.states[row], F(0))
        for occupied, value in targets.items():
            assert F(int(actual.counts[index[occupied], column]), scale) <= value
    assert exact_mass == 1
    proof = f.budget(actual)
    assert proof['complete_l1_error'] == exact_mass-F(sum(map(int, actual.counts.ravel())), scale)
    assert proof['finite_half_nuclear_error'] == proof['complete_l1_error']/2
    assert proof['exact_timeout_probability'] == F(565, 630)**2
    assert proof['complete_l1_error'] < F(1, 10**7)
    if feedback:
        # The control genuinely places a particle across the inter-carrier seam.
        assert any(data.read_ports[1] in occupied for occupied in expected[
            tuple(int(p in initial) for p in data.read_ports), (1,), 'complete'])


def test_remote_waits_change_actual_states_and_exclusion_is_not_independent_walks():
    data = f.source()
    remote, incident, repair = f.operators()
    initial = tuple(sorted((data.read_ports[0], data.remote[0][0])))
    index = data.states.index(initial)
    assert remote[:, index].nnz > 1
    # Both occupied endpoints cannot hop onto each other. The complete fair
    # native repair fixes this state, retaining both particles and total mass.
    blocked = data.states.index(data.read_ports)
    column = repair[:, blocked].toarray().ravel()
    assert np.count_nonzero(column) == 1 and column[blocked] == 2
    assert np.all(np.asarray(remote.sum(axis=0)).ravel() == 1130)
    assert np.all(np.asarray(incident.sum(axis=0)).ravel() == 130)
    assert remote.shape == (28680, 28680)


def test_finer_integer_grid_improves_each_complete_original_entry():
    coarse = f.execute(cap=4, ticks=2, grid_bits=30, feedback=True)
    fine = f.execute(cap=4, ticks=2, grid_bits=40, feedback=True)
    assert coarse.records == fine.records
    assert np.all(fine.counts//(1 << 10) >= coarse.counts)
    assert f.budget(fine)['complete_l1_error'] < f.budget(coarse)['complete_l1_error']


def test_complete_menus_include_impossible_initial_labels_and_every_timeout():
    data = f.source()
    initial = data.read_ports
    result = f.execute(cap=1, ticks=4, initial=initial)
    assert len(result.records) == 124
    assert len(set(result.records)) == 124
    assert sum(status == 'complete' for _, _, status in result.records) == 64
    assert sum(status == 'timeout' for _, _, status in result.records) == 60
    for column, (label, _, _) in enumerate(result.records):
        if label != (1, 1):
            assert not np.any(result.counts[:, column])
    proof = f.budget(result)
    assert sum(proof['per_record_mass_lower']) <= 1
    assert proof['exact_timeout_probability'] == 1-F(65, 630)**4


def test_current_only_tracker_matches_full_source_paths_with_remote_and_feedback_moves():
    data = f.source()
    internal = next(edge for edge in data.incident if data.read_ports[0] in edge
                    and edge != data.read_ports)
    tape = (data.remote[0], data.read_ports, internal, data.remote[1], data.read_ports)
    starts = (data.read_ports, tuple(sorted((data.read_ports[0], data.remote[0][0]))),
              tuple(sorted(data.remote[0])))
    for initial in starts:
        label = tuple(int(port in initial) for port in data.read_ports)
        for coins in product((0, 1), repeat=len(tape)):
            state, records = initial, ()
            for event, (seam, coin) in enumerate(zip(tape, coins, strict=True), 1):
                state, current = native.binary_step(state, seam, coin)
                if seam in data.incident:
                    records += ((event, *seam, coin, current),)
                assert f.readback(label, records) == tuple(int(port in state) for port in data.read_ports)
    audit = f.readback_audit()
    assert audit['all_global_state_seam_coin_addresses'] == 36136800
    assert audit['feedback_addresses_already_included'] == 57360


def test_resources_include_feedback_in_native_cost_and_incident_clock():
    data = f.source()
    initial = tuple(sorted((data.read_ports[0], data.remote[0][0])))
    actual = f.execute(cap=2, ticks=1, grid_bits=40, feedback=True, initial=initial)
    resource = f.resources(actual)
    expected = _direct_one_read(initial, 2, True)
    expected_feedback = sum((sum(population.values(), F(0)) for (_, word, status), population
                             in expected.items() if status == 'complete' and word == (1,)), F(0))
    assert F(resource['expected_feedback_native_repairs_exact']) == expected_feedback
    assert resource['maximum_total_native_attempts'] == 3
    assert resource['maximum_total_incident_opportunities'] == 2
    assert f.read_clock_values((1, 0, 1, 0), feedback=True) == (1, 3, 4, 5)
    assert f.read_clock_values((0, 1, 0, 1), feedback=True) == (1, 2, 3, 4)
    assert resource['retained_terminal_integer_bytes'] == actual.counts.nbytes


def test_actual_carrier_faces_supply_full_sector_transfer_floor():
    proof = f.transfer_certificate()
    assert proof['B_lower'] == proof['symmetric_T_lower'] == '1/21'
    assert len(proof['actual_carrier_triangle_covers']) == 2
    for row in proof['actual_carrier_triangle_covers']:
        assert row['internal_seams'] == 30 and row['triangles'] == 20
        assert row['edge_triangle_multiplicity'] == 2


def test_production_resolves_connected_cross_carrier_fields_and_causal_feedback():
    result = f.summary()
    assert result['all_pass'] is True
    for name in ('baseline', 'record_controlled_repair'):
        row = result[name]
        assert row['complete_native_state_dimension'] == math.comb(240, 2)
        assert row['complete_output_count'] == len(row['complete_outcome_mass_lower']) == 124
        assert F(row['complete_l1_rounding_error']) < F(1, 10**6)
        assert F(row['exact_timeout_probability']) < F(1, 10**11)
    low, high = map(F, result['baseline']['fourth_cumulant_interval'])
    assert low <= high < -F(1, 100000)
    low, high = map(F, result['second_field_feedback_minus_baseline_interval'])
    assert F(1, 1000) < low <= high
    for lower, upper in result['baseline']['field_means']:
        assert F(lower) <= F(2, 240) <= F(upper)


def test_fourth_cumulant_partition_bound_has_correct_independent_lipschitz_constant():
    partitions = f._partitions((0, 1, 2, 3))
    assert len(partitions) == len(set(partitions)) == 15
    assert sum(math.factorial(len(p)-1)*len(p) for p in partitions) == 75


@pytest.mark.parametrize('name,bad', (
    ('cap', True), ('cap', 0), ('cap', 513), ('cap', F(2)),
    ('ticks', True), ('ticks', 0), ('ticks', 5), ('ticks', 2.),
    ('grid_bits', True), ('grid_bits', 15), ('grid_bits', 52), ('grid_bits', 100),
    ('feedback', 1), ('feedback', None),
    ('initial', (1, True)), ('initial', (2, 1)), ('initial', (1, 1)),
    ('initial', (1, 240)), ('initial', [1, 2])))
def test_original_parameters_rejected_before_cached_execution(name, bad):
    f.execute(cap=1, ticks=1)
    options = {'cap': 1, 'ticks': 1, name: bad}
    with pytest.raises(ValueError):
        f.execute(**options)


def test_source_operators_and_cached_result_do_not_expose_cached_mutable_envelopes():
    data = f.source()
    seam = data.read_ports
    data.__dict__.update(read_ports=(0, 1), ports=1)
    assert f.source().read_ports == seam and f.source().ports == 240
    remote, _, _ = f.operators()
    remote.data[:] = 0
    assert f.operators()[0].nnz > 0 and np.any(f.operators()[0].data)
    result = f.execute(cap=1, ticks=1)
    seal = result.seal
    with pytest.raises(ValueError):
        result.counts.flags.writeable = True
    result.__dict__.update(cap=0, seal='wrong', counts=np.zeros((1, 1), dtype=np.int64))
    fresh = f.execute(cap=1, ticks=1)
    assert fresh.seal == seal
    f.validate_execution(fresh)


def test_resealed_positive_wrong_complete_law_is_rejected_by_cold_source_replay():
    original = f.execute(cap=1, ticks=1, initial=(1, 12))
    wrong = np.zeros_like(original.counts)
    wrong[0] = original.counts.sum(axis=0)
    assert wrong.sum() == original.counts.sum()
    changed = replace(original, counts=wrong, seal='')
    changed = replace(changed, seal=f._seal(changed))
    f._compute.cache_clear()
    for operation in (f.validate_execution, f.budget, f.field_statistics, f.resources, f.summary):
        with pytest.raises(ValueError):
            operation(changed)


def test_counter_guards_and_original_integer_array_dispatch():
    data = f.source()
    a, b = data.read_ports
    for label, records in (((True, 0), ()), ((0, 0), ((1, a, b, 0, 1),)),
                           ((1, 0), ((True, a, b, 0, 1),)),
                           ((1, 0), ((1, *data.remote[0], 0, 0),))):
        with pytest.raises(ValueError):
            f.readback(label, records)
    class DispatchTrap(np.ndarray):
        def __getattribute__(self, name):
            if name in ('shape', 'dtype', 'tobytes'):
                raise AssertionError('subclass inspected before rejection')
            return super().__getattribute__(name)
    actual = f.execute(cap=1, ticks=1)
    with pytest.raises(ValueError):
        f.validate_execution(replace(actual, counts=actual.counts.view(DispatchTrap)))


def test_status_subclass_cannot_relabel_the_original_history_during_validation():
    produced = f.execute(cap=1, ticks=1, grid_bits=40, initial=(1, 12))
    genuine_interval = f.field_statistics(produced)['moments'][(0,)]
    assert genuine_interval[0] <= F(62, 65) <= genuine_interval[1]

    class RelabelAtComparison(str):
        calls = 0
        victim = None
        replacement = None

        def __eq__(self, other):
            self.calls += 1
            self.victim.__dict__['records'] = self.replacement
            return str.__eq__(self, other)

        __hash__ = str.__hash__

    label, word, original_status = produced.records[0]
    status = RelabelAtComparison(original_status)
    declared = ((label, word, status), *produced.records[1:])
    claimed = replace(produced, records=declared)
    status.victim = claimed
    status.replacement = tuple((label, tuple(1-bit for bit in word), state)
                              if state == 'complete' else (label, word, state)
                              for label, word, state in produced.records)
    # JSON encodes the same string payload, so source hashes alone cannot
    # reject this object. On the old reader the callback relabelled the
    # caller after budget() validated its copy, producing a wrong field law.
    assert f._seal(claimed) == produced.seal
    with pytest.raises(ValueError, match='record'):
        f.field_statistics(claimed)
    assert status.calls == 0
    assert claimed.records is declared


@pytest.mark.parametrize('location', ('record', 'label', 'word', 'label_bit', 'word_bit',
                                     'status', 'seal', 'cap', 'initial'))
def test_original_record_grammar_and_header_types_reject_callbacks_before_dispatch(location):
    class TupleTrap(tuple):
        def __iter__(self):
            raise AssertionError('tuple subclass iterated before type rejection')

        def __eq__(self, other):
            raise AssertionError('tuple subclass compared before type rejection')

    class IntTrap(int):
        def __eq__(self, other):
            raise AssertionError('integer subclass compared before type rejection')

        def __bool__(self):
            raise AssertionError('integer subclass tested before type rejection')

    class StringTrap(str):
        def __eq__(self, other):
            raise AssertionError('string subclass compared before type rejection')

    original = f.execute(cap=1, ticks=1, grid_bits=40, initial=(1, 12))
    if location in ('seal', 'cap', 'initial'):
        changed = {'seal': StringTrap(original.seal), 'cap': IntTrap(original.cap),
                   'initial': TupleTrap(original.initial)}[location]
        supplied = replace(original, **{location: changed})
    else:
        label, word, status = original.records[-1]
        if location == 'label':
            label = TupleTrap(label)
        elif location == 'word':
            word = TupleTrap(word)
        elif location == 'label_bit':
            label = (IntTrap(label[0]), label[1])
        elif location == 'word_bit':
            word = (IntTrap(word[0]),)
        elif location == 'status':
            status = StringTrap(status)
        record = (label, word, status)
        if location == 'record':
            record = TupleTrap(record)
        supplied = replace(original, records=(*original.records[:-1], record))
    with pytest.raises(ValueError):
        f.validate_execution(supplied)
