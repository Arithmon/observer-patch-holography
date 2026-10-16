"""Independent full-configuration source oracle for density-compatible fields.

The oracle propagates all 2**N binary configurations using the original
endpoint rule and exact rational masses. It does not use exclusion duality,
independence after evolution, a stopped-kernel inverse, or production arrays.
"""
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

from ..observer import native


PATH = (4, ((0, 1), (1, 2), (2, 3)), ((0, 1), (1, 2)))
TRIANGLE = (3, ((0, 1), (0, 2), (1, 2)), ((0, 1), (0, 2), (1, 2)))


def _coin_step(population, seams, denominator):
    result = defaultdict(F)
    for state, mass in population.items():
        for seam in seams:
            for coin in (0, 1):
                final, (retained_coin, current) = native.native(state, seam, coin)
                assert retained_coin == coin and current == state[seam[0]]-final[seam[0]]
                result[final] += mass/denominator
    return dict(result)


def _interval(population, seams, incident, cap, *, remote_identity=False, reverse=False):
    """Enumerate first-hit alternatives and preserve the actual timeout state."""
    remote = tuple(edge for edge in seams if edge not in incident)
    denominator = 2*len(seams)

    def wait(state):
        if remote_identity:
            return {configuration: mass*F(len(remote), len(seams))
                    for configuration, mass in state.items()}
        return _coin_step(state, remote, denominator)

    waiting = population
    stopped = defaultdict(F)
    for _ in range(cap):
        for state, mass in _coin_step(waiting, incident, denominator).items():
            stopped[state] += mass
        waiting = wait(waiting)
    if reverse:
        # Deliberately wrong chronology, used only as a negative control.
        stopped = defaultdict(F)
        advanced = _coin_step(population, incident, denominator)
        for _ in range(cap):
            for state, mass in advanced.items():
                stopped[state] += mass
            advanced = wait(advanced)
    return dict(stopped), waiting


def _records(ticks):
    labels = tuple(product((0, 1), repeat=2))
    complete = tuple((label, first, last, ticks, 'complete')
                     for label in labels for first, last in product((0, 1), repeat=2))
    timeouts = tuple((label, first, None, completed, 'timeout')
                     for completed in range(ticks) for label in labels
                     for first in ((None,) if completed == 0 else (0, 1)))
    return (*complete, *timeouts)


def _forward_oracle(graph, *, ticks, cap, feedback, rho=F(1, 3),
                    remote_identity=False, reverse=False):
    vertices, seams, incident = graph
    initial = {state: rho**sum(state)*(1-rho)**(vertices-sum(state))
               for state in product((0, 1), repeat=vertices)}
    live = {(label, None): {state: mass for state, mass in initial.items()
                           if state[:2] == label}
            for label in product((0, 1), repeat=2)}
    terminal = {record: {} for record in _records(ticks)}
    for completed in range(ticks):
        following = {}
        for (label, first), population in live.items():
            arrived, waiting = _interval(population, seams, incident, cap,
                                         remote_identity=remote_identity, reverse=reverse)
            terminal[label, first, None, completed, 'timeout'] = waiting
            selected = {first: arrived} if completed else {
                bit: {state: mass for state, mass in arrived.items() if state[0] == bit}
                for bit in (0, 1)}
            for bit, state in selected.items():
                if completed == 0 and feedback and bit:
                    state = _coin_step(state, ((0, 1),), 2)
                if completed+1 == ticks:
                    for last in (0, 1):
                        terminal[label, bit, last, ticks, 'complete'] = {
                            configuration: mass for configuration, mass in state.items()
                            if configuration[1] == last}
                else:
                    following[label, bit] = state
        live = following
    return terminal


def _probabilities(terminal):
    return {record: sum(population.values(), F()) for record, population in terminal.items()}


def _production(graph, **kwargs):
    from . import refinement
    vertices, edges, incident = graph
    return refinement.law(tuple(range(vertices)), edges, incident, **kwargs)


@pytest.mark.parametrize('graph', (PATH, TRIANGLE))
@pytest.mark.parametrize('ticks', (2, 3))
@pytest.mark.parametrize('cap', (1, 2, 3, 4))
@pytest.mark.parametrize('feedback', (False, True))
def test_native_oracle_retains_initial_disclosure_and_every_timeout_prefix(graph, ticks, cap, feedback):
    law = _probabilities(_forward_oracle(graph, ticks=ticks, cap=cap, feedback=feedback))
    assert len(law) == 16+4+8*(ticks-1)
    assert sum(law.values(), F()) == 1
    q = F(len(graph[1])-len(graph[2]), len(graph[1]))**cap
    for completed in range(ticks):
        assert sum((mass for record, mass in law.items()
                    if record[4] == 'timeout' and record[3] == completed), F()) == (1-q)**completed*q
    for label in product((0, 1), repeat=2):
        assert sum((mass for record, mass in law.items() if record[0] == label), F()) == (
            F(1, 3)**sum(label)*F(2, 3)**(2-sum(label)))
    produced = _production(graph, ticks=ticks, cap=cap, feedback=feedback, bits=40)
    intervals = produced.intervals()
    assert len(produced.records) == len(set(produced.records)) == len(law)
    assert set(intervals) == set(law)
    for record, exact in law.items():
        low, high = intervals[record]
        assert type(low) is type(high) is F
        assert 0 <= low <= exact <= high <= 1
    assert produced.timeout == 1-(1-q)**ticks
    assert sum((law[record]-intervals[record][0] for record in law), F()) == produced.mass_deficit
    assert produced.mass_deficit < F(1, 10**6)
    intervals.clear()
    assert len(produced.intervals()) == len(law)


@pytest.mark.parametrize('bits', (8, 16, 44))
@pytest.mark.parametrize('rho', (F(1, 5), F(1, 3), F(4, 5)))
def test_entire_law_is_enclosed_even_at_coarse_grid_and_other_product_densities(bits, rho):
    oracle = _probabilities(_forward_oracle(PATH, ticks=3, cap=4, feedback=True, rho=rho))
    result = _production(PATH, ticks=3, cap=4, feedback=True, bits=bits, rho=rho)
    for record, (low, high) in result.intervals().items():
        assert low <= oracle[record] <= high
    assert sum(result.lower_mass) <= 1 <= sum(result.upper_mass)


def test_complete_law_resolves_a_small_exact_native_feedback_witness():
    record = ((1, 0), 1, 1, 2, 'complete')
    passive = _production(PATH, ticks=2, cap=1, feedback=False, bits=44)
    active = _production(PATH, ticks=2, cap=1, feedback=True, bits=44)
    assert passive.intervals()[record][0] <= F(13, 486) <= passive.intervals()[record][1]
    assert active.intervals()[record][0] <= F(17, 486) <= active.intervals()[record][1]
    assert passive.intervals()[record][1] < active.intervals()[record][0]
    assert passive.timeout == active.timeout == F(5, 9)


@pytest.mark.parametrize('changed', ({'ticks': True}, {'cap': True}, {'bits': True},
                                     {'rho': 1/3}, {'feedback': 1}))
def test_public_law_rejects_numeric_aliases_before_cached_production(changed):
    kwargs = dict(ticks=2, cap=1, bits=16, rho=F(1, 3), feedback=False)
    _production(PATH, **kwargs)
    with pytest.raises(ValueError):
        _production(PATH, **dict(kwargs, **changed))


def test_undisclosed_remote_changes_to_a_read_port_cannot_claim_the_local_dictionary():
    # Removing (1,2) from the clock would hide genuine changes of read port1.
    # A generic hidden-state measurement law is not this current-record observer.
    wrong_clock = (PATH[0], PATH[1], ((0, 1),))
    with pytest.raises(ValueError):
        _production(wrong_clock, ticks=2, cap=1, bits=16)


@pytest.mark.parametrize('reader', ('intervals', 'mass_deficit', 'statistics', 'validate', 'export'))
@pytest.mark.parametrize('change', ('lower', 'upper', 'timeout', 'policy', 'missing',
                                   'nested_boolean', 'float_grid', 'boolean_port'))
def test_all_field_readers_reject_fabricated_source_laws(reader, change):
    from . import refinement
    law = _production(PATH, ticks=2, cap=2, bits=20, feedback=True)
    if change == 'lower':
        forged = law._replace(lower_mass=(law.lower_mass[0]+F(1, 2**40), *law.lower_mass[1:]))
    elif change == 'upper':
        forged = law._replace(upper_mass=(law.upper_mass[0]+F(1, 2**40), *law.upper_mass[1:]))
    elif change == 'timeout':
        forged = law._replace(timeout=law.timeout+F(1, 2**40))
    elif change == 'policy':
        forged = law._replace(feedback=False)
    elif change == 'missing':
        forged = law._replace(records=law.records[:-1])
    elif change == 'nested_boolean':
        first = ((False, law.records[0][0][1]), *law.records[0][1:])
        assert first == law.records[0]  # Python equality would conceal the alias.
        forged = law._replace(records=(first, *law.records[1:]))
    elif change == 'float_grid':
        forged = law._replace(bits=float(law.bits))
    else:
        forged = law._replace(ports=(False, *law.ports[1:]))
    with pytest.raises(ValueError):
        if reader == 'intervals':
            forged.intervals()
        elif reader == 'mass_deficit':
            _ = forged.mass_deficit
        elif reader == 'export':
            refinement.export(forged, F(0))
        else:
            getattr(refinement, reader)(forged)


@pytest.mark.parametrize('feedback', (False, True))
def test_reported_four_point_cumulant_encloses_direct_full_configuration_moments(feedback):
    from . import refinement
    oracle = _probabilities(_forward_oracle(TRIANGLE, ticks=3, cap=1, feedback=feedback))
    result = _production(TRIANGLE, ticks=3, cap=1, feedback=feedback, bits=44)
    moments = {}
    for degree in range(1, 5):
        for subset in combinations(range(4), degree):
            moments[subset] = sum((mass for (initial, first, last, _, status), mass in oracle.items()
                if status == 'complete' and all((*initial, first, last)[i] for i in subset)), F())
    # The explicit four-variable cumulant formula is independent of the
    # production partition generator and its interval arithmetic.
    fourth = moments[0, 1, 2, 3]
    cumulant = fourth-sum((moments[i,]*moments[tuple(j for j in range(4) if j != i)]
                          for i in range(4)), F())
    cumulant -= sum((moments[0, i]*moments[tuple(j for j in range(1, 4) if j != i)]
                    for i in range(1, 4)), F())
    cumulant += 2*sum((moments[pair]*math.prod(moments[j,] for j in range(4) if j not in pair)
                      for pair in combinations(range(4), 2)), F())
    cumulant -= 6*math.prod(moments[i,] for i in range(4))
    reported = refinement.statistics(result)
    for i, (low, high) in enumerate(reported['means']):
        assert low <= moments[i,] <= high
    assert reported['joint_four'][0] <= fourth <= reported['joint_four'][1]
    assert reported['connected_four'][0] <= cumulant <= reported['connected_four'][1]


@pytest.mark.parametrize('value', (F(0), F(1), F(1, 3), F(-1, 3), F(1, 2**100), F(-1, 2**100)))
def test_compact_rational_export_rounds_outward_even_for_negative_cumulants(value):
    from . import refinement
    low, high = refinement.dyadic(value), refinement.dyadic(value, upper=True)
    assert low <= value <= high
    assert 0 <= value-low < F(1, 2**80)
    assert 0 <= high-value < F(1, 2**80)


def test_export_preserves_all_records_and_charges_its_additional_rounding():
    from . import refinement
    law = _production(PATH, ticks=3, cap=4, bits=40, feedback=True)
    exported = refinement.export(law, F(0))
    lower = {}
    for row in exported['records']:
        record = (tuple(row['initial_bits']), row['first_field'], row['final_field'],
                  row['scheduled_ticks'], row['status'])
        assert record not in lower
        low, high = map(F, row['probability_interval'])
        exact_low, exact_high = law.intervals()[record]
        assert 0 <= low <= exact_low <= exact_high <= high <= 1
        lower[record] = low
    assert set(lower) == set(law.records)
    assert F(exported['capped_complete_readout_l1_allowance']) >= 1-sum(lower.values())
    assert F(exported['interval_cap_timeout_allowance']) >= law.timeout
    assert F(exported['eventual_readout_missing_mass_allowance']) >= 1-sum(
        mass for record, mass in lower.items() if record[-1] == 'complete')


@pytest.mark.parametrize('seed,feedback,cap', ((5, True, 500000), (1, True, 500000), (1, False, 1)))
def test_full_source_projection_uses_both_actual_initial_bits_and_continuing_currents(seed, feedback, cap):
    from ..native_action import clock_histories as source
    from . import refinement
    result = refinement.history(seed, feedback=feedback, cap=cap)
    initial = source.execute(seed, ticks=4, cap=cap, feedback=feedback, stop_after_slots=0)
    final = source.execute(seed, ticks=4, cap=cap, feedback=feedback)
    assert result['external_union_slots_used'] == final.slots <= cap
    assert result['serialized_continuation_preserves_full_states_records_and_rng']
    assert result['serialized_checkpoint_bytes'] > sum(source.source().ports)
    for side, row in enumerate(result['screens']):
        addresses = dict(source.source().matched_addresses[side])
        expected_initial = [initial.states[side][addresses[port]] for port in (0, 1)]
        assert row['initial_bits'] == expected_initial
        scheduled = [record for record in final.records[side] if record[0] == 'scheduled']
        assert row['first_field'] == (scheduled[0][-1] if scheduled else None)
        assert row['scheduled_ticks'] == len(scheduled) == final.scheduled_hits
        assert row['native_feedback_attempts'] == int(feedback and bool(scheduled) and scheduled[0][-1] == 1)
        assert row['feedback_signed_currents'] == [record[7] for record in final.records[side]
                                                   if record[0] == 'feedback']
        assert row['total_incident_ticks'] == len(final.records[side]) == (
            row['scheduled_ticks']+row['native_feedback_attempts'])
        assert row['native_attempts'] == final.uniform_attempts[side]+final.feedback_attempts[side]
        complete = final.scheduled_hits == 4
        assert row['status'] == ('complete' if complete else 'timeout')
        assert row['final_field'] == (final.states[side][addresses[1]] if complete else None)
        assert row['terminal_state_sha256'] == hashlib.sha256(final.states[side]).hexdigest()
        reconstructed = [expected_initial[port]+sum(record[7]*(int(record[5] == port)-int(record[4] == port))
                             for record in final.records[side]) for port in (0, 1)]
        assert reconstructed == [final.states[side][addresses[port]] for port in (0, 1)]
        if seed == 5:
            assert row['feedback_signed_currents'] == [1]
    if cap == 1:
        assert all(row['status'] == 'timeout' and row['final_field'] is None for row in result['screens'])


@pytest.fixture(scope='module')
def refinement_receipt():
    from . import refinement
    return refinement.run()


def test_retained_refinement_receipt_replays_actual_source_and_rejects_forgery(refinement_receipt):
    from tools import strict_json
    from ..observer.scalar_response.experiment import compare_receipt
    retained = strict_json.load(Path(__file__).with_name('refinement_results.json'))
    compare_receipt(retained, refinement_receipt)

    def exact_types(value):
        assert type(value) in (dict, list, str, int, bool, type(None))
        if type(value) is dict:
            assert all(type(key) is str for key in value)
            for item in value.values():
                exact_types(item)
        elif type(value) is list:
            for item in value:
                exact_types(item)

    exact_types(refinement_receipt)
    assert [len(row['records']) for row in retained['field_protocols']] == [44, 44]
    for mutate in (
        lambda r: r['source'].__setitem__('matched_ports', 221),
        lambda r: r['field_protocols'][0]['records'].pop(),
        lambda r: r['field_protocols'][0]['records'][0].__setitem__('first_field', False),
        lambda r: r['field_protocols'][0]['records'][0].__setitem__('probability_interval', ['0', '0']),
        lambda r: r['field_protocols'][0].__setitem__('capped_complete_readout_l1_allowance', '0'),
        lambda r: r['field_protocols'][1].__setitem__('either_full_screen_eventual_connected_four', ['0', '0']),
        lambda r: r.__setitem__('source_bindings', {}),
    ):
        forged = deepcopy(retained)
        mutate(forged)
        with pytest.raises(ValueError):
            compare_receipt(forged, refinement_receipt)


def test_real_refinement_cli_rejects_a_changed_scientific_record(tmp_path):
    from tools import strict_json
    from ..observer.scalar_response.experiment import canonical
    forged = strict_json.load(Path(__file__).with_name('refinement_results.json'))
    forged['field_protocols'][0]['records'][0]['probability_interval'] = ['0', '0']
    path = tmp_path/'wrong-refinement.json'
    path.write_text(canonical(forged), encoding='utf-8')
    process = subprocess.run([sys.executable, '-m', 'model.observer_time_fields.refinement',
                              '--check', str(path)], text=True, capture_output=True)
    assert process.returncode != 0
    assert 'Receipt mismatch' in process.stderr and 'probability_interval' in process.stderr


@pytest.fixture(scope='module')
def weighted_effective_source():
    import sympy as sp
    states = tuple(product((0, 1), repeat=4))
    index = {state: i for i, state in enumerate(states)}
    operators = {}
    for edge in ((0, 1), (0, 2), (1, 2), (2, 3)):
        matrix = sp.zeros(len(states))
        for column, state in enumerate(states):
            for coin in (0, 1):
                final, _ = native.native(state, edge, coin)
                matrix[index[final], column] += sp.Rational(1, 2)
        assert matrix.T == matrix and matrix*matrix == matrix
        operators[edge] = matrix
    identity = sp.eye(len(states))
    remote = operators[2, 3]/4
    incident = sum((operators[edge] for edge in ((0, 1), (0, 2), (1, 2))), sp.zeros(len(states)))/4
    inverse_root = identity+(2/sp.sqrt(3)-1)*operators[2, 3]
    completed = incident*(identity+operators[2, 3]/3)
    effective = (inverse_root*incident*inverse_root).applyfunc(sp.simplify)
    assert completed*(identity-remote) == incident
    assert (inverse_root*(identity-remote)*inverse_root-identity).applyfunc(sp.simplify) == sp.zeros(len(states))
    # rho=1/5 makes the square-root Bernoulli boundary exactly rational,
    # while retaining every particle-number sector with nonuniform weights.
    omega = sp.Matrix([sp.Rational(2**(4-sum(state)), 25) for state in states])
    probability = omega.applyfunc(lambda value: value**2)
    assert (omega.T*omega)[0] == sum(probability) == 1
    projectors = {(port, bit): sp.diag(*(int(state[port] == bit) for state in states))
                  for port in range(3) for bit in (0, 1)}
    return states, operators[0, 1], completed, effective, omega, probability, projectors


def test_weighted_effective_transfer_matches_every_adaptive_full_density_field_word(weighted_effective_source):
    import sympy as sp
    states, repair, completed, effective, omega, probability, projectors = weighted_effective_source
    ones, identity = sp.ones(1, len(states)), sp.eye(len(states))
    for feedback in (False, True):
        total = sp.Rational(0)
        for a, b, first, last in product((0, 1), repeat=4):
            initial = projectors[0, a]*projectors[1, b]
            intervention = repair if feedback and first else identity
            direct = (ones*projectors[1, last]*completed**2*intervention
                      *projectors[0, first]*completed*initial*probability)[0]
            transformed = (omega.T*projectors[1, last]*effective**2*intervention
                           *projectors[0, first]*effective*initial*omega)[0]
            assert sp.simplify(direct-transformed) == 0
            total += direct
        assert total == 1


def test_effective_dictionary_does_not_accept_nonstationary_boundaries_or_remote_fields(weighted_effective_source):
    import sympy as sp
    states, _, completed, effective, omega, probability, projectors = weighted_effective_source
    ones = sp.ones(1, len(states))
    changed = sp.Matrix([sp.Rational(2**sum(state[i] == i % 2 for i in range(4)), 25)
                         for state in states])
    assert (changed.T*changed)[0] == 1
    direct = (ones*projectors[0, 1]*completed*projectors[1, 0]
              *changed.applyfunc(lambda value: value**2))[0]
    transformed = (changed.T*projectors[0, 1]*effective*projectors[1, 0]*changed)[0]
    assert sp.simplify(direct-transformed) == -sp.Rational(3, 10000)
    direct = (ones*projectors[2, 1]*completed*projectors[2, 1]*probability)[0]
    transformed = (omega.T*projectors[2, 1]*effective*projectors[2, 1]*omega)[0]
    assert sp.simplify(direct-transformed) == -sp.Rational(7, 300)+sp.sqrt(3)/75
    assert sp.simplify(direct-transformed) != 0


@pytest.fixture(scope='module')
def observer_prefixes():
    from ..native_action import clock_histories as source
    from . import refinement
    full = source.execute(5, ticks=4, feedback=True)
    first_slot = next(row[1] for row in full.records[0] if row[0] == 'scheduled')
    prefixes = tuple(source.execute(5, ticks=4, feedback=True, stop_after_slots=slot)
                     for slot in (first_slot-1, first_slot))
    encoded = tuple(refinement.checkpoint_to_json(prefix) for prefix in prefixes)
    return full, prefixes, encoded


@pytest.mark.parametrize('after_first', (False, True))
def test_observer_checkpoint_continues_saved_memory_using_only_the_new_suffix(observer_prefixes, after_first):
    from ..native_action import clock_histories as source
    from . import refinement
    full, prefixes, encoded = observer_prefixes
    restored = refinement.checkpoint_from_json(encoded[int(after_first)])
    assert restored.source == prefixes[int(after_first)]
    assert restored.first_fields == ((1, 1) if after_first else (None, None))
    continued = source.execute(5, ticks=4, feedback=True, checkpoint=restored.source)
    assert continued == full
    for side in (0, 1):
        suffix = continued.records[side][len(restored.source.records[side]):]
        actual = refinement._counter_view(restored.counters[side], suffix, restored.first_fields[side])
        uninterrupted = refinement._counter_view(restored.initial_bits[side], full.records[side])
        assert actual == uninterrupted
        if after_first:
            assert refinement._counter_view(restored.counters[side], suffix)[1] is None
            with pytest.raises(ArithmeticError):
                refinement._counter_view(restored.counters[side], full.records[side], restored.first_fields[side])


def test_observer_checkpoint_rejects_forged_memory_aliases_and_source_ancestry(observer_prefixes):
    from ..native_action import clock_histories as source
    from . import refinement
    _, _, encoded = observer_prefixes
    good = json.loads(encoded[1])

    def bad_ancestry(data):
        inner = json.loads(data['source_checkpoint'])
        inner['source_digest'] = '0'*64
        data['source_checkpoint'] = json.dumps(inner)

    for mutate in (
        lambda d: d.__setitem__('extra', 0),
        lambda d: d.pop('initial_bits'),
        lambda d: d.__setitem__('schema', 'other-schema'),
        lambda d: d['initial_bits'][0].__setitem__(0, 1-d['initial_bits'][0][0]),
        lambda d: d['initial_bits'][0].__setitem__(0, bool(d['initial_bits'][0][0])),
        lambda d: d['counters'][0].__setitem__(1, 1-d['counters'][0][1]),
        lambda d: d['counters'][0].__setitem__(1, float(d['counters'][0][1])),
        lambda d: d['first_fields'].__setitem__(0, 0),
        lambda d: d['first_fields'].__setitem__(0, True),
        bad_ancestry,
    ):
        forged = deepcopy(good)
        mutate(forged)
        with pytest.raises(ValueError):
            refinement.checkpoint_from_json(json.dumps(forged))
    # A warm replay cache must not be the reason the numeric alias is rejected.
    source._compute.cache_clear()
    forged = deepcopy(good)
    forged['first_fields'][0] = True
    with pytest.raises(ValueError):
        refinement.checkpoint_from_json(json.dumps(forged))
    with pytest.raises(ValueError):
        refinement.checkpoint_from_json('{"schema":"duplicate",'+encoded[1][1:])


def test_oracle_detects_removed_remote_motion_and_record_controlled_repair():
    kwargs = {'ticks': 3, 'cap': 4, 'feedback': True}
    actual = _probabilities(_forward_oracle(PATH, **kwargs))
    held = _probabilities(_forward_oracle(PATH, remote_identity=True, **kwargs))
    passive = _probabilities(_forward_oracle(PATH, **dict(kwargs, feedback=False)))
    assert actual != held
    assert actual != passive


def test_full_terminal_oracle_detects_reversed_stopping_even_when_local_law_cancels_it():
    kwargs = {'ticks': 3, 'cap': 4, 'feedback': True}
    actual = _forward_oracle(PATH, **kwargs)
    reversed_order = _forward_oracle(PATH, reverse=True, **kwargs)
    # Stationary local cylinder probabilities alone cannot certify chronology.
    # The similarity commutes with our reads and internal intervention, so a
    # genuine source-order test needs terminal states or a nonstationary input.
    assert _probabilities(actual) == _probabilities(reversed_order)
    assert actual != reversed_order
