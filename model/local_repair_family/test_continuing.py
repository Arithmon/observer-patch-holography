"""Independent endpoint/FIFO oracles and hostile stationary-archive controls."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from . import continuing as c


def params(mix=F(1), repair=F(1)):
    return c.local.Parameters(port_weights=(25, 25)+(8,)*10,
                              c_weight=mix, repair_probability=repair)


def index(pair):
    return c.STATES.index(tuple(int(i in pair) for i in range(12)))


def independent_map(record):
    if record.kind == 'hold':
        return dict(enumerate(range(66)))
    a, b = record.seam
    result = {}
    for i, state in enumerate(c.STATES):
        total = state[a]+state[b]
        left = total//2+record.coin*(total % 2)
        current = state[a]-left
        raw = {0: 0, 1: 1, -1: 2}[current]
        if record.law == 'E':
            raw += 4*(total % 2)
        if raw == record.raw:
            final = list(state); final[a], final[b] = left, total-left
            result[i] = c.STATES.index(tuple(final))
    return result


def trace(density):
    return sum((v for (i, a, j, b), v in density.items() if i == j and a == b), F())


def reference(density):
    result = {}
    for (i, a, j, b), value in density.items():
        if i == j:
            result[a, b] = result.get((a, b), 0)+value
    return {key: value for key, value in result.items() if value}


@pytest.mark.parametrize('mix,repair', [(F(1), F(1)), (F(0), F(1)), (F(1, 3), F(3, 4))])
def test_every_raw_slot_and_reference_column_has_independent_native_action(mix, repair):
    p = params(mix, repair)
    vector = {(index((0, 2)), 2): .5, (index((2, 3)), 5): .5j,
              (index((1, 7)), 2): -.5, (index((4, 8)), 5): .5}
    density = {(i, a, j, b): x*complex(y).conjugate()
               for (i, a), x in vector.items() for (j, b), y in vector.items()}
    outputs = list(c.step_blocks(density, (), p))
    assert len(outputs) == 721 and len({fifo[0] for fifo, _, _ in outputs}) == 721
    total = {}
    for fifo, exported, actual in outputs:
        record = fifo[0]; mapping = independent_map(record)
        assert exported is None
        if record.kind == 'hold':
            w = 1-repair
        else:
            a, b = record.seam
            w = repair*F(1, 30)*(mix if record.law == 'C' else 1-mix)
            w *= F(p.port_weights[a if record.coin else b], p.port_weights[a]+p.port_weights[b])
        expected = {(mapping[i], a, mapping[j], b): w*v for (i, a, j, b), v in density.items()
                    if w and i in mapping and j in mapping and v}
        assert actual == expected
        for key, value in actual.items():
            total[key] = total.get(key, 0)+value
    assert trace(total) == pytest.approx(1, abs=1e-15)
    for key, value in reference(density).items():
        assert reference(total)[key] == pytest.approx(value, abs=1e-15)


@pytest.mark.parametrize('mix,repair', [(F(1), F(1)), (F(0), F(1)), (F(1, 3), F(3, 4))])
def test_stationary_source_fifo_is_actual_supported_joint_state(mix, repair):
    p = params(mix, repair)
    sigma = c.local.reference(p)
    actual = c.stationary(p)
    independent = {}
    for branch in c.branches(p):
        row = {(j, 0, j, 0): branch.weight*sigma[i]
               for i, j in independent_map(branch.record).items() if branch.weight}
        if row:
            independent[(branch.record,)] = row
    assert actual == independent
    assert c.advance_marginal(actual, p) == actual
    assert trace(c.source_sum(actual)) == 1
    assert c.certificate(p)['exact_joint_source_one_slot_stationarity']


def test_complete_history_never_renormalizes_and_reconstructs_all_exports():
    p = params(F(1, 3), F(3, 4))
    tape = (c.Record('repair', (0, 1), 'C', 0, 1, 1, 1), c.Record('hold'),
            c.Record('repair', (2, 3), 'E', 0, 5, 1, 1), c.Record('hold'))
    initial = {(index((0, 2)), 4, index((0, 2)), 4): F(7, 9)}
    expected = initial
    catalogue = {b.record: b for b in c.branches(p)}
    for record in tape:
        mapping = independent_map(record); weight = catalogue[record].weight
        expected = {(mapping[i], a, mapping[j], b): weight*v for (i, a, j, b), v in expected.items()
                    if i in mapping and j in mapping}
    assert expected and trace(expected) < F(1, 1000)
    for capacity in (1, 2, 4):
        replay = c.history(initial, tape, p, capacity=capacity)
        assert replay['density'] == expected
        assert replay['exported']+replay['fifo'] == tape
        assert len(replay['exported']) == max(0, 4-capacity)


def test_two_slot_arbitrary_reference_identity_uses_only_unread_old_prefix():
    p = params()
    x, y = index((0, 2)), index((2, 3))
    # An entangled input with reference labels beyond a two-row implementation.
    initial = {(x, 2, x, 2): F(1, 2), (y, 5, y, 5): F(1, 2),
               (x, 2, y, 5): -.5j, (y, 5, x, 2): .5j}
    active = [b.record for b in c.branches(p) if b.record.law == 'C' and independent_map(b.record)]
    def manual(density, record):
        a, b = record.seam
        weight = F(p.port_weights[a if record.coin else b], 30*(p.port_weights[a]+p.port_weights[b]))
        mapping = independent_map(record)
        return {(mapping[i], r, mapping[j], s): weight*v for (i, r, j, s), v in density.items()
                if i in mapping and j in mapping and v}
    unread = {}
    for record in active:
        for key, value in manual(initial, record).items():
            unread[key] = unread.get(key, 0)+value
    expected = {}
    for first in active:
        middle = manual(unread, first)
        if not middle: continue
        for last in active:
            output = manual(middle, last)
            if output: expected[first, last] = output
    actual = {(): initial}
    for _ in range(3): actual = c.advance_marginal(actual, p, capacity=2)
    assert actual.keys() == expected.keys()
    for fifo in actual:
        assert actual[fifo].keys() == expected[fifo].keys()
        for key in actual[fifo]:
            assert actual[fifo][key] == pytest.approx(expected[fifo][key], abs=1e-18)
    for key, value in reference(initial).items():
        assert reference(c.source_sum(actual))[key] == pytest.approx(value, abs=2e-14)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_word_reference_is_exact_corner_and_modular_ratio_cancels_scalar(law):
    p = params(F(1, 2))
    label = 1
    first = c.Record('repair', (0, 5), law, 0, 1 if law == 'C' else 5, label, 1)
    second = c.Record('repair', (6, 7), law, 0, 0, 0, 0)
    tape = first, second
    row = c.word_reference(p, tape)
    direct = c.stationary_word(p, tape)
    assert direct == {(j, 0, j, 0): v for j, v in zip(row['support'], row['diagonal'], strict=True)}
    assert len(row['support']) > 1 and trace(direct) == row['mass']
    sigma = c.local.reference(p)
    for (i, j), (k, ell) in product(row['pairs'], repeat=2):
        assert sigma[j]/sigma[ell] == sigma[i]/sigma[k]
    normalized = tuple(v/row['mass'] for v in row['diagonal'])
    assert all(a/b == sigma[i]/sigma[j]
               for i, a in zip(row['support'], normalized) for j, b in zip(row['support'], normalized))


def test_zero_slots_stay_zero_but_cannot_supply_a_modular_corner():
    p = params()
    initial = {(0, 0, 0, 0): F(1)}
    zero = c.Record('repair', (0, 1), 'C', 0, 3)
    result = c.history(initial, (zero, c.Record('hold')), p)
    assert result['density'] == {} and result['exported'] == (zero,)
    with pytest.raises(ValueError, match='zero branches'):
        c.word_reference(p, (zero,))
    with pytest.raises(ValueError, match='zero branches'):
        c.word_reference(p, (c.Record('hold'),))  # Hold probability zero.


def test_stationary_record_cannot_be_replaced_by_an_independent_source_copy():
    p = params(); sigma = c.local.reference(p); omega = c.stationary(p)
    record = c.Record('repair', (0, 1), 'C', 0, 1, 1, 1)
    actual = omega[(record,)]; event = trace(actual)
    complement = [i for i, x in enumerate(c.STATES) if not (x[0] == 0 and x[1] == 1)]
    assert sum(actual.get((i, 0, i, 0), F()) for i in complement) == 0
    assert event*sum(sigma[i] for i in complement) == F(22020, 6759003)
    assert c.archive_counterexample()['source_times_FIFO_product_claim_is_false']


def test_first_archive_bit_keeps_input_reference_information_forever():
    p = params(); catalogue = c.branches(p)
    wanted = c.Record('repair', (0, 1), 'C', 0, 1, 1, 1)
    # Store the first-event indicator beside the untouched input tag. This is a
    # classical marginal of the full archive, not a modification of the source.
    initial = {(index((0, 2)), 0, index((0, 2)), 0): F(1, 2),
               (index((2, 3)), 1, index((2, 3)), 1): F(1, 2)}
    after = {}
    for branch in catalogue:
        flag = int(branch.record == wanted)
        for (i, a, j, b), value in c._apply(initial, branch).items():
            key = i, 2*flag+a, j, 2*flag+b
            after[key] = after.get(key, F())+value
    joint = reference(after)
    expected = {(0, 0): F(59, 120), (1, 1): F(1, 2), (2, 2): F(1, 120)}
    assert joint == expected
    for _ in range(3):
        after = c.source_sum(c.advance_marginal({(): after}, p))
        assert reference(after) == expected
    record_marginal = (F(119, 120), F(1, 120))
    distance = sum(abs(expected.get((2*a+r, 2*a+r), F())-record_marginal[a]/2)
                   for a, r in product(range(2), repeat=2))/2
    assert distance == F(1, 120)


def test_mixing_bound_applies_after_fifo_filling_and_not_to_full_archive(monkeypatch):
    from . import mixing
    observed = []
    monkeypatch.setattr(mixing, 'mixing_bound', lambda p, n: observed.append(n) or F(1, 9))
    assert c.preparation_bound(params(), 64, capacity=3) == F(1, 9)
    assert observed == [61]
    with pytest.raises(ValueError, match='filled'):
        c.preparation_bound(params(), 2, capacity=3)
    assert c.resource_counts(64, capacity=3)['exported_records'] == 61
    assert not c.certificate(params())['full_archive_stationarity_or_reference_factorization_claimed']


def test_scalar_generator_compression_does_not_imply_scalar_finite_tick():
    theta = (2, 1, 2)+(1,)*9
    V = np.zeros((66, 2))
    for column, pairs in enumerate((((0, 2), (1, 3)), ((0, 3), (1, 2)))):
        for pair in pairs: V[index(pair), column] = 1/math.sqrt(2)
    # Independently apply the floor/ceiling endpoint oracle to both codewords.
    # Flagged correction requires each complete branch Gram to be scalar;
    # this checks the premise rather than only its later phase consequence.
    for branch in c.branches(params(F(1, 2))):
        mapping = independent_map(branch.record)
        image = np.zeros_like(V)
        for i, j in mapping.items(): image[j] += V[i]
        gram = image.T@image
        np.testing.assert_allclose(gram, gram[0, 0]*np.eye(2), atol=0, rtol=0)
    powers = np.array([x[0]+x[2] for x in c.STATES])
    np.testing.assert_allclose(V.T@np.diag(powers)@V, np.eye(2), atol=3e-16)
    np.testing.assert_allclose(np.diag((-1.)**powers)@V, V@np.diag([1, -1]), atol=0)
    halfway = np.diag((1j)**powers)@V
    assert np.linalg.norm(halfway-V@(V.T@halfway)) == pytest.approx(1)
    row = c.correctable_clock_certificate()
    assert row['finite_tick_countercontrol']['theta'] == list(theta)
    assert not row['finite_compressed_phase_scalar_claimed_without_invariance']


@pytest.mark.parametrize('damage', ['current', 'coin_bool', 'raw_bool', 'law', 'hold_fields'])
def test_malformed_records_cannot_change_the_declared_controller(damage):
    r = c.Record('repair', (0, 1), 'C', 0, 1, 1, 1)
    edits = {'current': {'current': -1}, 'coin_bool': {'coin': False}, 'raw_bool': {'raw': True},
             'law': {'law': 'coherent'}, 'hold_fields': {'kind': 'hold'}}
    with pytest.raises(ValueError):
        c.history({(0, 0, 0, 0): F(1)}, (replace(r, **edits[damage]),), params())


@pytest.mark.parametrize('value', [float('nan'), complex(0, float('inf')), True, '1/2'])
def test_nonfinite_or_disguised_density_coefficients_are_refused(value):
    with pytest.raises(ValueError):
        c.history({(0, 0, 0, 0): value}, (), params())


def test_fraction_and_public_result_mutations_do_not_change_later_instruments():
    p = params(F(1, 2))
    first = c.branches(p)
    selected = next(b for b in first if b.weight)
    selected.weight._numerator = -1
    assert all(b.weight >= 0 for b in c.branches(p))
    original = c.stationary(p)
    next(iter(original.values())).clear()
    assert trace(c.source_sum(c.stationary(p))) == 1
    p.c_weight._numerator = 2; p.c_weight._denominator = 4
    with pytest.raises(ValueError, match='canonical Fraction'):
        c.branches(p)


def test_long_positive_history_refuses_float_underflow_and_keeps_exact_mass():
    p = params(F(1), F(1)); x = index((2, 3))
    record = c.Record('repair', (0, 1), 'C', 0, 0, 0, 0)
    tape = (record,)*256
    expected = F(1, 60)**256
    for scalar in (1, F(1)):
        actual = c.history({(x, 0, x, 0): scalar}, tape, p)
        assert actual['density'] == {(x, 0, x, 0): expected}
        assert len(actual['exported']) == 255
    with pytest.raises(ArithmeticError, match='underflow'):
        c.history({(x, 0, x, 0): 1.0}, tape, p)
    entangled_phase = {(x, 0, x, 0): .5, (x, 1, x, 1): .5,
                       (x, 0, x, 1): -.5j, (x, 1, x, 0): .5j}
    with pytest.raises(ArithmeticError, match='underflow'):
        c.history(entangled_phase, tape, p)


def test_numeric_component_precision_boundary_keeps_valid_neighbor_and_exact_tiny_mass():
    p = params(F(1), F(1)); x = index((2, 3)); key = x, 0, x, 0
    record = c.Record('repair', (0, 1), 'C', 0, 0, 0, 0)
    tiny = F(1, 10**400)
    assert c.apply_branch({key: tiny}, p, record) == {key: tiny/60}
    adjacent = c.apply_branch({key: 1e-321}, p, record)[key]
    assert adjacent > 0 and adjacent == F(1, 60)*1e-321
    for bad in (1e-323, complex(1e-323, .5), complex(.5, 1e-323)):
        with pytest.raises(ArithmeticError, match='underflow'):
            c.apply_branch({key: bad}, p, record)
    # A structurally absent outcome remains an exact zero, not an underflow.
    zero = c.Record('repair', (0, 1), 'C', 0, 3)
    assert c.apply_branch({key: 1e-323}, p, zero) == {}


def test_committed_history_snapshots_caller_records_and_repeated_occurrences():
    p = params(F(1), F(1)); x = index((2, 3))
    record = c.Record('repair', (0, 1), 'C', 0, 0, 0, 0)
    result = c.history({(x, 0, x, 0): F(1)}, (record, record), p)
    for name in ('raw', 'label', 'current'):
        object.__setattr__(record, name, 1)
    assert result['fifo'][0].raw == result['exported'][0].raw == 0
    assert result['fifo'][0] is not result['exported'][0]
    assert c.history({(x, 0, x, 0): F(1)},
                     result['exported']+result['fifo'], p)['density'] == result['density']
    object.__setattr__(result['exported'][0], 'raw', 1)
    assert result['fifo'][0].raw == 0


def test_step_generator_snapshots_inputs_and_detaches_each_yield():
    p = params(F(1), F(1)); x = index((2, 3)); key = x, 0, x, 0
    old = tuple(c.Record('repair', (0, 1), 'C', 0, 0, 0, 0) for _ in range(2))
    density = {key: F(1)}
    stream = c.step_blocks(density, old, p, capacity=2)
    p.c_weight._numerator = 0
    density[key]._numerator = 2
    for record in old:
        for name in ('raw', 'label', 'current'):
            object.__setattr__(record, name, 1)
    first, exported, _ = next(stream)
    assert first[0].raw == exported.raw == 0
    for record in (first[0], exported):
        for name in ('raw', 'label', 'current'):
            object.__setattr__(record, name, 1)
    second, exported_again, output = next(stream)
    assert second[0].raw == exported_again.raw == 0
    assert output == {key: F(1, 60)}
