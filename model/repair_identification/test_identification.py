"""Independent native endpoint and hostile-input controls for identification.

The endpoint oracle below does not use a coherence table, the probe prediction
formula, or the production channel application.  It groups the actual primitive
endpoint images by their committed native records before discarding labels.
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import hashlib
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from model.observer import native
from model.record_laws import family as f
from model.repair_identification import core


STATES = native.states()
INDEX = {state: index for index, state in enumerate(STATES)}


def _table(law):
    values = (1, F(1, 2), F(1, 2), 0, F(1, 2)) if law == 'C' else (1, 0, 0, 0, 0)
    return dict(zip(('k', 'g', 'h', 'm', 'n'), map(f.S, values), strict=True))


def _native_images(law):
    images = []
    for state in STATES:
        row = []
        for coin in (0, 1):
            output, (_, current) = native.native(state, (0, 1), coin)
            private = (coin, current if law == 'C' else state[0] - state[1])
            row.append((private, current, INDEX[output]))
        images.append(tuple(row))
    return tuple(images)


def _native_operator(law, density):
    """Exact common-current CP blocks of the original C/E endpoint rule."""
    images = _native_images(law)
    blocks = {current: {} for current in (-1, 0, 1)}
    for (left, right), value in density.items():
        for coin in (0, 1):
            label, current, a = images[left][coin]
            other, _, b = images[right][coin]
            if label == other:
                block = blocks[current]
                block[a, b] = block.get((a, b), f.S()) + value * F(1, 2)
    return {current: {unit: value for unit, value in block.items() if value}
            for current, block in blocks.items()}


def _control_operator(law, density):
    if law in ('C', 'E'):
        return _native_operator(law, density)
    if law == 'E_collision_half':
        first, second = (_native_operator(candidate, density) for candidate in ('C', 'E'))
        result = {0: second[0]}
        for current in (-1, 1):
            result[current] = {unit: (first[current].get(unit, f.S())
                                      + second[current].get(unit, f.S())) * F(1, 2)
                               for unit in first[current].keys() | second[current].keys()}
        return result
    assert law == 'C_phase_quarter'
    images = _native_images('C')
    result = {current: {} for current in (-1, 0, 1)}
    for (left, right), value in density.items():
        for coin in (0, 1):
            label, current, a = images[left][coin]
            other, _, b = images[right][coin]
            if label != other:
                continue
            phase_a = f.S(0, 0, 1) if current == 0 and sum(STATES[left][:2]) == 1 else f.S(1)
            phase_b = f.S(0, 0, 1) if current == 0 and sum(STATES[right][:2]) == 1 else f.S(1)
            result[current][a, b] = (result[current].get((a, b), f.S())
                                     + phase_a * phase_b.conj() * value * F(1, 2))
    return {current: {unit: value for unit, value in block.items() if value}
            for current, block in result.items()}


def _parts(operator):
    return {unit: value.parts() for unit, value in operator.items() if value}


def _probe_input(probe):
    x, y = probe.input_pair
    phase = f.S(1) if probe.phase == 0 else f.S(0, 0, 1)
    return {(x, x): f.S(F(1, 2)), (y, y): f.S(F(1, 2)),
            (x, y): phase.conj() * F(1, 2), (y, x): phase * F(1, 2)}


def _project(density, pair, plus):
    """Full rank-one Lueders branch, with the complement left coherent."""
    def column(source):
        result = {target: F(1, 2) for target in pair} if source in pair else {}
        if plus:
            return result
        result = {target: -value for target, value in result.items()}
        result[source] = result.get(source, F()) + 1
        return {target: value for target, value in result.items() if value}

    out = {}
    for (left, right), value in density.items():
        for a, c in column(left).items():
            for b, d in column(right).items():
                out[a, b] = out.get((a, b), f.S()) + value * c * d
    return {unit: value for unit, value in out.items() if value}


def _trace(density):
    return sum((value for (a, b), value in density.items() if a == b), f.S())


def _native_probe(law, probe):
    repaired = _native_operator(law, _probe_input(probe))
    return {(current, plus): _project(block, probe.output_pair, plus)
            for current, block in repaired.items() for plus in (False, True)}


@pytest.mark.parametrize('law', ('C', 'E'))
def test_current_maps_match_actual_native_endpoints_on_every_operator_entry(law):
    # At fixed current, x -> x+jD is injective. Distinct complex coefficients
    # therefore compare all 6,084 matrix units without cancellation between
    # source units or thousands of repeated admissibility validations.
    incoming = {(a, b): f.S(78*a+b+1, 0, 79*b-a+1)
                for a, b in product(range(78), repeat=2)}
    expected = _native_operator(law, incoming)
    for current in (-1, 0, 1):
        assert _parts(core.apply_current(_table(law), current, incoming)) == _parts(expected[current])


@pytest.mark.parametrize('law', ('C', 'E'))
def test_every_probe_probability_comes_from_a_complete_native_instrument(law):
    predicted = core.predictions(_table(law))
    assert len(predicted) == len(core.PROBES) == 9
    for probe, value in zip(core.PROBES, predicted, strict=True):
        outputs = _native_probe(law, probe)
        assert len(outputs) == 6
        assert sum((_trace(rho) for rho in outputs.values()), f.S()) == 1
        assert _trace(outputs[probe.current, True]) == value
        # A target current is an outcome, not a condition for accepting a shot.
        if probe.parameter != 'k':
            assert sum((_trace(rho) for (j, _), rho in outputs.items()
                        if j != probe.current), f.S()).sign() > 0


@pytest.mark.parametrize('law,control', (('C_phase_quarter', 'C_phase_quarter'),
                                       ('E_collision_half', 'E_n_quarter')))
def test_phase_and_collision_controls_follow_actual_endpoint_record_operations(law, control):
    predicted = core.predictions(core.controls()[control])
    for probe, value in zip(core.PROBES, predicted, strict=True):
        blocks = _control_operator(law, _probe_input(probe))
        outputs = {(current, plus): _project(block, probe.output_pair, plus)
                   for current, block in blocks.items() for plus in (False, True)}
        assert sum((_trace(rho) for rho in outputs.values()), f.S()) == 1
        assert _trace(outputs[probe.current, True]) == value


def _interior_table():
    # Strict diagonal dominance of all current blocks, with every one of the
    # nine real coordinates nonzero, prevents an accidentally rank-deficient
    # named-candidate test from concealing a phase-sign or inverse error.
    return f.table(f.S(F(1, 8), 0, F(1, 16)),
                   f.S(F(1, 16), 0, F(-1, 32)),
                   f.S(F(-1, 16), 0, F(1, 32)), F(1, 16),
                   f.S(F(1, 8), 0, F(-1, 16)))


def test_nine_real_coordinates_are_reconstructed_with_their_phase_signs():
    table = _interior_table()
    assert f.positive_semidefinite(table)
    reconstructed = core.reconstruct(core.predictions(table))
    assert {name: value.parts() for name, value in reconstructed.items()} == {
        name: value.parts() for name, value in table.items()}
    # The imaginary preparation has rho_xy=-i/2, so the positive imaginary
    # multiplier raises the joint plus probability above its diagonal baseline.
    imaginary = {probe.parameter: value for probe, value in
                 zip(core.PROBES, core.predictions(table), strict=True) if probe.phase == 1}
    assert (imaginary['k'] - f.S(F(1, 2))).sign() > 0
    assert (imaginary['g'] - f.S(F(3, 8))).sign() < 0
    assert (imaginary['h'] - f.S(F(3, 8))).sign() > 0
    assert (imaginary['n'] - f.S(F(3, 8))).sign() < 0


def test_all_nine_coordinates_vary_independently_inside_the_cp_domain():
    base = f.table(0, 0, 0, 0, 0)
    baseline = core.predictions(base)
    columns = []
    for parameter, imaginary in (('k', False), ('k', True), ('g', False), ('g', True),
                                  ('h', False), ('h', True), ('m', False),
                                  ('n', False), ('n', True)):
        changed = deepcopy(base)
        changed[parameter] = f.S(0, 0, F(1, 16)) if imaginary else f.S(F(1, 16))
        assert f.positive_semidefinite(changed)
        column = []
        for initial, value in zip(baseline, core.predictions(changed), strict=True):
            difference = (value-initial)*16
            assert difference.b == difference.c == difference.d == 0
            column.append(difference.a)
        columns.append(tuple(column))
    assert tuple(columns) == tuple(tuple(F(1, 2) if i == j else F(0)
                                         for i in range(9)) for j in range(9))


def test_collision_coherence_is_invisible_on_the_entire_binary_shell():
    e = _table('E')
    collision = _table('E')
    collision['n'] = f.S(F(1, 4))
    assert f.positive_semidefinite(collision)
    binary = tuple(INDEX[state] for state in native.states(binary=True))
    incoming = {(a, b): f.S(78*a+b+1, 0, b-a)
                for a, b in product(binary, repeat=2)}
    for current in (-1, 0, 1):
        assert _parts(core.apply_current(e, current, incoming)) == _parts(
            core.apply_current(collision, current, incoming))
    original, changed = core.predictions(e), core.predictions(collision)
    differences = {probe.name: b-a for probe, a, b in
                   zip(core.PROBES, original, changed, strict=True) if a != b}
    assert differences == {'n_re': f.S(F(1, 8))}
    assert core.reconstruct(changed)['n'] == F(1, 4)


@pytest.mark.parametrize('parameter', ('h', 'm', 'n'))
def test_existing_empty_seam_probes_do_not_identify_the_remaining_parameters(parameter):
    first = f.table(0, 0, 0, 0, 0)
    second = deepcopy(first)
    second[parameter] = f.S(F(1, 8))
    assert f.positive_semidefinite(second)
    values_a, values_b = core.predictions(first), core.predictions(second)
    for probe, a, b in zip(core.PROBES, values_a, values_b, strict=True):
        assert (a != b) == (probe.parameter == parameter and probe.phase == 0)


def test_complete_current_bound_detects_the_collision_only_channel_difference():
    first, second = _table('E'), _table('E')
    second['n'] = f.S(F(1, 4))
    # On the n probe the changed current block is [[0,1/8],[1/8,0]];
    # its half trace norm is 1/8. Thus this is also an attained lower bound,
    # despite equality on every binary input and equality of all diagonals.
    assert core.instrument_bound(first, second) == F(1, 8)
    assert core.instrument_bound(first, first) == 0
    assert core.history_bound(3, F(1, 8), F(1, 32), F(1, 64)) == F(27, 64)
    assert core.history_bound(10, F(1, 8)) == 1
    assert core.coordinate_error_bound(F(1, 100)) == F(7, 100)
    assert core.coordinate_error_bound(F(1, 7)) == 1


def test_common_current_interface_cannot_be_replaced_by_an_unread_channel():
    # Let L_0,L_1 be C's coin-zero wait/move branches. Replacing them by
    # (L_0 +/- L_1)/sqrt(2) preserves their unread CP sum exactly but changes
    # the record instrument. A basis input on (1,0) originally moves with
    # probability 1/2; labelling the new plus branch "move" gives 1/4.
    x = INDEX[tuple(1 if p in (0, 2) else 0 for p in range(12))]
    images = _native_images('C')
    old = {}
    for source, row in enumerate(images):
        for label, _, target in row:
            old.setdefault(label, {})[source] = (target, f.ROOT_HALF)
    wait, move = old[(0, 0)], old[(0, 1)]
    assert not wait.keys() & move.keys()
    mixed = {}
    for sign in (-1, 1):
        branch = {a: (b, amplitude * f.ROOT_HALF) for a, (b, amplitude) in wait.items()}
        branch.update({a: (b, sign * amplitude * f.ROOT_HALF)
                       for a, (b, amplitude) in move.items()})
        mixed[sign] = branch
    # All-input equality follows coefficient by coefficient, including the
    # cross-current terms whose signs cancel only after forgetting the label.
    for a, b in product(range(78), repeat=2):
        def image(branches):
            out = {}
            for branch in branches:
                if a in branch and b in branch:
                    u, left = branch[a]
                    v, right = branch[b]
                    out[u, v] = out.get((u, v), f.S()) + left * right.conj()
            return _parts(out)
        assert image((wait, move)) == image(tuple(mixed.values()))
    assert move[x][1] * move[x][1].conj() == F(1, 2)
    assert mixed[1][x][1] * mixed[1][x][1].conj() == F(1, 4)
    assert core.instrument_bound(_table('C'), _table('C')) == 0


class _IntAlias(int):
    pass


class _FractionAlias(F):
    pass


class _SurdAlias(f.S):
    pass


class _DictAlias(dict):
    pass


class _TupleAlias(tuple):
    pass


@pytest.mark.parametrize('bad', (True, 0.5, float('nan'), float('inf'),
                                complex(1, 0), _IntAlias(1), _FractionAlias(1),
                                _SurdAlias(1)))
def test_exact_inverse_rejects_original_numeric_aliases_and_inexact_values(bad):
    values = list(core.predictions(_table('E')))
    values[0] = bad
    with pytest.raises(ValueError):
        core.reconstruct(tuple(values))


@pytest.mark.parametrize('numerator,denominator',
                         ((-1, -1), (1, 0), (1.0, 1), (1, 1.0),
                          (True, 1), (1, True), (_IntAlias(1), 1),
                          (1, _IntAlias(1)), (2, 2), (0, 2)))
def test_exact_inverse_rejects_corrupted_fraction_internals(numerator, denominator):
    corrupted = F(1)
    corrupted._numerator = numerator
    corrupted._denominator = denominator
    values = list(core.predictions(_table('E')))
    values[0] = corrupted
    with pytest.raises(ValueError):
        core.reconstruct(tuple(values))


@pytest.mark.parametrize('component', ('a', 'b', 'c', 'd'))
def test_exact_inverse_validates_original_surd_components(component):
    value = f.S(F(1, 2))
    setattr(value, component, True)
    values = list(core.predictions(_table('E')))
    values[0] = value
    with pytest.raises(ValueError):
        core.reconstruct(tuple(values))


def test_inverse_rejects_an_impossible_table_without_projecting_to_cp():
    values = list(core.predictions(_table('E')))
    values[next(i for i, p in enumerate(core.PROBES) if p.name == 'g_re')] += F(1, 1000)
    # k=1 forces h=g. Every probability remains in [0,1], but this changed
    # g with unchanged h cannot define a completely positive repair law.
    with pytest.raises(ValueError):
        core.reconstruct(tuple(values))


@pytest.mark.parametrize('change', ('missing', 'extra', 'list', 'alias'))
def test_inverse_requires_the_complete_fixed_probe_tuple(change):
    values = core.predictions(_table('E'))
    malformed = {'missing': values[:-1], 'extra': (*values, f.S()),
                 'list': list(values), 'alias': _TupleAlias(values)}[change]
    with pytest.raises(ValueError):
        core.reconstruct(malformed)


def test_validated_tables_and_predictions_do_not_retain_mutable_surd_aliases():
    def snapshot(value):
        return tuple((part.numerator, part.denominator) for part in value.parts())

    original = _interior_table()
    expected = {name: snapshot(value) for name, value in original.items()}
    detached = core.validate_table(original)
    detached['g'].a._numerator = 100
    assert {name: snapshot(value) for name, value in original.items()} == expected
    prediction = core.predictions(original)
    honest = tuple(snapshot(value) for value in prediction)
    prediction[0].a._numerator = 100
    assert tuple(snapshot(value) for value in core.predictions(original)) == honest


@pytest.mark.parametrize('count', (True, -1, 1.0, _IntAlias(1)))
def test_history_bound_requires_a_genuine_nonnegative_call_count(count):
    with pytest.raises(ValueError):
        core.history_bound(count, F(1, 10))


@pytest.mark.parametrize('bad', (True, F(-1, 10), 0.1, float('nan'),
                                _IntAlias(0), _FractionAlias(0)))
def test_error_budgets_are_exact_nonnegative_original_values(bad):
    with pytest.raises(ValueError):
        core.coordinate_error_bound(bad)
    with pytest.raises(ValueError):
        core.history_bound(2, bad)


def test_sparse_operator_boundary_rejects_numeric_key_aliases_and_wrong_menus():
    for current in (True, 0.0, _IntAlias(0), 2):
        with pytest.raises(ValueError):
            core.apply_current(_table('C'), current, {(0, 0): f.S(1)})
    for density in ({(True, 0): f.S(1)}, {(0, 78): f.S(1)},
                    {_TupleAlias((0, 0)): f.S(1)}, _DictAlias({(0, 0): f.S(1)})):
        with pytest.raises(ValueError):
            core.apply_current(_table('C'), 0, density)


def _numeric(value):
    return complex(float(value.a) + np.sqrt(2)*float(value.b),
                   float(value.c) + np.sqrt(2)*float(value.d))


def _lifted_control(law, density, reference):
    blocks = {current: {} for current in (-1, 0, 1)}
    for ra, rb in product(range(reference), repeat=2):
        component = {(a//reference, b//reference): value
                     for (a, b), value in density.items()
                     if a % reference == ra and b % reference == rb}
        for current, output in _control_operator(law, component).items():
            blocks[current].update({(a*reference+ra, b*reference+rb): value
                                    for (a, b), value in output.items()})
    return blocks


def _lifted_project(density, pair, plus, reference):
    result = {}
    for ra, rb in product(range(reference), repeat=2):
        component = {(a//reference, b//reference): value
                     for (a, b), value in density.items()
                     if a % reference == ra and b % reference == rb}
        projected = _project(component, pair, plus)
        result.update({(a*reference+ra, b*reference+rb): value
                       for (a, b), value in projected.items()})
    return result


def _joint_probe_input(probe, reference):
    if reference == 1:
        return _probe_input(probe)
    x, y = probe.input_pair
    phase = f.S(1) if probe.phase == 0 else f.S(0, 0, 1)
    a, b = x*2, y*2+1
    return {(a, a): f.S(F(1, 2)), (b, b): f.S(F(1, 2)),
            (a, b): phase.conj()*F(1, 2), (b, a): phase*F(1, 2)}


def _native_read(law, probe, incoming, reference):
    return {(current, pointer): _lifted_project(block, probe.output_pair, pointer == 1, reference)
            for current, block in _lifted_control(law, incoming, reference).items()
            for pointer in (-1, 1)}


def _encoded_density(density, reference):
    from model.observer.source_discrimination import circuit
    out = np.zeros((8*reference, 8*reference), complex)
    for (a, b), value in density.items():
        ca, cb = circuit.SOURCE_CODES[a//reference], circuit.SOURCE_CODES[b//reference]
        assert ca < 8 and cb < 8
        out[ca*reference+a % reference, cb*reference+b % reference] = _numeric(value)
    return out


def _source_coarse(branches, reference, history=False):
    """Return active blocks AND a bound on all omitted full-code entries.

    This is a test comparison, not a normalized/projection-based state. For
    v=(a,b), ||vv*-aa*||_1 <= 2||a||||b||+||b||^2 includes both coherences
    with all 504 spectator rows and the complete outside population.
    """
    result, leakage = {}, 0.
    for branch in branches:
        full = branch.amplitude.reshape(512, reference)
        active, outside = full[:8].reshape(-1), full[8:].reshape(-1)
        leakage += 2*np.linalg.norm(active)*np.linalg.norm(outside) + np.linalg.norm(outside)**2
        if history:
            key = tuple(None if raw is None or raw[-1] is None else (current, raw[-1])
                        for current, raw in zip(branch.currents, branch.record, strict=True))
        else:
            key = None if branch.abort else (branch.current, branch.record[-1])
        result[key] = result.get(key, np.zeros((8*reference,)*2, complex)) + np.outer(active, active.conj())
    return result, leakage


@pytest.mark.parametrize('law', ('C', 'E', 'C_phase_quarter', 'E_collision_half'))
def test_all_emitted_probes_match_independent_complete_native_output_states(law):
    from model.repair_identification import source
    for probe in core.PROBES:
        executed = source.execute_probe(law, probe.name)
        assert len({branch.record for branch in executed.branches}) == len(executed.branches)
        actual, leakage = _source_coarse(executed.branches, 1)
        expected = _native_read(law, probe, _probe_input(probe), 1)
        assert set(actual) == {*expected, None}
        error = leakage
        for key, block in actual.items():
            target = np.zeros_like(block) if key is None else _encoded_density(expected[key], 1)
            error += np.sum(np.abs(block-target))
        assert error < 2e-9
        assert sum(float(np.vdot(branch.amplitude, branch.amplitude).real)
                   for branch in executed.branches) == pytest.approx(1, abs=2e-10)
        assert any(branch.abort for branch in executed.branches)
        for branch in executed.branches:
            assert dict(branch.resources)['declared_blank_data_preparations'] == 1
            assert dict(branch.resources)['native_repairs'] == 1
            assert branch.amplitude.flags.writeable is False


@pytest.mark.parametrize('law', ('C', 'E', 'C_phase_quarter', 'E_collision_half'))
def test_emitted_two_event_history_preserves_reference_and_disturbed_state(law):
    from model.repair_identification import source
    from model.observer.source_discrimination import circuit
    probe = next(p for p in core.PROBES if p.name == 'g_re')
    reference = 2
    prefix = _native_read(law, probe, _joint_probe_input(probe, reference), reference)
    expected = {}
    for first, density in prefix.items():
        controlled = {}
        for (a, b), value in density.items():
            phase = 1
            if first == (0, 1):
                phase = (-1)**((circuit.SOURCE_CODES[a//reference] & 1)
                              + (circuit.SOURCE_CODES[b//reference] & 1))
            controlled[a, b] = phase*value
        expected.update({(first, second): output for second, output in
                         _native_read(law, probe, controlled, reference).items()})
    executed = source.continue_history(law, probe.name, feedback=True, reference=True)
    actual, leakage = _source_coarse(executed.branches, reference, history=True)
    assert len({branch.record for branch in executed.branches}) == len(executed.branches)
    error = leakage
    for key, block in actual.items():
        target = _encoded_density(expected.get(key, {}), reference)
        error += np.sum(np.abs(block-target))
    assert set(expected) <= set(actual)
    assert error < 8e-9
    for first in prefix:
        probability = sum(np.trace(block).real for (a, _), block in actual.items() if a == first)
        assert probability == pytest.approx(_numeric(_trace(prefix[first])).real, abs=3e-10)
    for branch in executed.branches:
        costs = dict(branch.resources)
        assert costs['declared_blank_data_preparations'] == 1
        assert costs['joint_reference_preparations'] == 1
        assert costs['native_repairs'] == (1 if branch.record[1] is None else 2)
        first = branch.record[0]
        assert costs.get('feedback_phase_controls', 0) == int(
            first[-1] == 1 and branch.currents[0] == 0)


class _StringAlias(str):
    pass


class _HiddenArray(np.ndarray):
    def nonzero(self):
        return (np.array([], dtype=int),)


@pytest.mark.parametrize('bad', (np.zeros(511), np.zeros((512, 0)),
                                np.zeros(512, bool), np.zeros(512, object),
                                np.ma.array(np.zeros(512)),
                                np.zeros(512).view(_HiddenArray)))
def test_source_boundaries_validate_original_complete_arrays(bad):
    from model.repair_identification import source
    with pytest.raises(ValueError):
        source.repair_branches('C', bad)
    with pytest.raises(ValueError):
        source.readout_branches('g_re', bad)


@pytest.mark.parametrize('bad', (float('nan'), float('inf'), complex(0, float('inf'))))
def test_source_boundaries_do_not_hide_nonfinite_outside_entries(bad):
    from model.repair_identification import source
    value = np.zeros((512, 2), complex)
    value[511, 1] = bad
    with pytest.raises(ValueError):
        source.repair_branches('C', value)
    with pytest.raises(ValueError):
        source.readout_branches('g_re', value)


def test_source_protocol_names_and_choices_reject_original_type_aliases():
    from model.repair_identification import source
    with pytest.raises(ValueError):
        source.execute_probe(_StringAlias('C'), 'g_re')
    with pytest.raises(ValueError):
        source.execute_probe('C', _StringAlias('g_re'))
    with pytest.raises(ValueError):
        source.preparation('g_re', reference=1)
    with pytest.raises(ValueError):
        source.continue_history('C', 'g_re', feedback=1)


def test_source_pointer_preserves_the_unobserved_native_complement():
    from model.repair_identification import source
    from model.observer.source_discrimination import circuit
    # Two actual native states outside the probe's spectator-zero chart.
    # Both have the same deterministic minus read, but their relative phase
    # survives only a Lueders complement, not an occupation measurement.
    a, b = tuple(code for code in circuit.SOURCE_CODES if code >= 8)[:2]
    incoming = np.zeros(512, complex)
    incoming[a], incoming[b] = 1/np.sqrt(2), 1j/np.sqrt(2)
    outputs = source.readout_branches('g_re', incoming)
    assert set(outputs) == {-1, 1}
    assert np.linalg.norm(outputs[1]) < 2e-10
    assert np.linalg.norm(outputs[-1]-incoming) < 2e-10
    assert outputs[-1][a]*outputs[-1][b].conjugate() == pytest.approx(-0.5j, abs=2e-10)


def test_full_pointer_cache_agrees_with_literal_source_execution_on_reference_columns():
    from model.repair_identification import source
    from model.observer.source_discrimination import experiment as recorded
    incoming = np.zeros((512, 2), complex)
    incoming[0, 0], incoming[200, 1] = 0.5, 0.5j
    incoming[511, 0], incoming[8, 1] = -0.5, 0.5
    cached = source.readout_branches('g_re', incoming)
    literal = recorded.readout_branches(incoming)
    assert set(cached) == set(literal) == {-1, 1}
    for pointer in (-1, 1):
        assert np.linalg.norm(cached[pointer]-literal[pointer]) < 2e-13


def _complete_counts(law='C', attempted=65536):
    from model.repair_identification import calibration
    rows = []
    for probe in core.PROBES:
        blocks = _native_read(law, probe, _probe_input(probe), 1)
        counts = {name: 0 for name in calibration.OUTCOMES}
        for (current, pointer), block in blocks.items():
            mass = _trace(block)
            assert mass.b == mass.c == mass.d == 0
            count = mass.a*attempted
            assert count.denominator == 1
            counts[f'current_{current}_{"plus" if pointer == 1 else "minus"}'] = count.numerator
        rows.append({'probe': probe.name, 'attempted': attempted, 'counts': counts})
    return {'schema': calibration.SCHEMA, 'settings': rows}


def _analyze_counts(evidence):
    from model.repair_identification import calibration
    return calibration.analyze_counts(evidence, (F(1, 64),)*9,
                                      {'C': _table('C'), 'E': _table('E')}, 2, (F(0),)*9)


def test_count_interface_uses_all_attempts_and_charges_both_error_sources():
    from model.repair_identification import calibration
    evidence = _complete_counts()
    result = _analyze_counts(evidence)
    assert result['compatible_locked_candidates'] == ['C']
    assert result['physical_law_selected'] is False
    assert result['count_evidence']['simultaneous_sampling_failure_bound'] == '9/64'
    assert result['complete_history_diameter_bound'] == '7/16'
    # Make 1/8 of every setting explicit aborts by moving counts from an
    # occupied bin. The denominator remains attempted, including all aborts.
    changed = deepcopy(evidence)
    for row in changed['settings']:
        target = max((key for key in row['counts'] if key != 'abort'), key=row['counts'].get)
        row['counts'][target] -= row['attempted']//8
        row['counts']['abort'] += row['attempted']//8
    analyzed = calibration.analyze_counts(changed, (F(1, 64),)*9,
        {'C': _table('C')}, 2, (F(1, 8),)*9)
    assert analyzed['compatible_locked_candidates'] == ['C']
    assert analyzed['count_evidence']['abort_counts'] == [8192]*9
    assert analyzed['count_evidence']['premises_established_by_counts'] is False
    for row, center, probe in zip(changed['settings'], analyzed['centers'], core.PROBES, strict=True):
        assert F(center) == F(row['counts'][f'current_{probe.current}_plus'], row['attempted'])


@pytest.mark.parametrize('mutation', ('drop_abort', 'postselect', 'missing_probe', 'duplicate_probe',
                                     'reorder', 'extra_bin', 'bool_count', 'float_count',
                                     'alias_count', 'bool_attempted', 'alias_probe', 'alias_dict'))
def test_count_interface_rejects_incomplete_or_aliased_evidence(mutation):
    evidence = _complete_counts()
    row = evidence['settings'][0]
    if mutation == 'drop_abort':
        del row['counts']['abort']
    elif mutation == 'postselect':
        other = evidence['settings'][2]
        other['attempted'] = sum(value for key, value in other['counts'].items() if key.startswith('current_0'))
    elif mutation == 'missing_probe':
        evidence['settings'].pop()
    elif mutation == 'duplicate_probe':
        evidence['settings'][-1] = deepcopy(row)
    elif mutation == 'reorder':
        evidence['settings'][0], evidence['settings'][1] = evidence['settings'][1], row
    elif mutation == 'extra_bin':
        row['counts']['unread'] = 0
    elif mutation == 'bool_count':
        row['counts']['abort'] = False
    elif mutation == 'float_count':
        row['counts']['abort'] = 0.0
    elif mutation == 'alias_count':
        row['counts']['abort'] = _IntAlias(0)
    elif mutation == 'bool_attempted':
        row['attempted'] = True
    elif mutation == 'alias_probe':
        row['probe'] = _StringAlias(row['probe'])
    else:
        evidence = _DictAlias(evidence)
    with pytest.raises(ValueError):
        _analyze_counts(evidence)


def test_noisy_centers_are_not_cp_projected_or_mistaken_for_family_exclusion():
    from model.repair_identification import calibration
    centers = [value.a for value in core.predictions(_table('E'))]
    centers[2] += F(1, 1000)  # g != h at k=1 is not an admissible point.
    with pytest.raises(ValueError):
        core.reconstruct(tuple(centers))
    result = calibration.analyze_tubes(tuple(centers), (F(1, 500),)*9, {'E': _table('E')}, 3)
    assert result['compatible_locked_candidates'] == ['E']
    assert result['cp_projection_performed'] is False
    zero = tuple(value.a for value in core.predictions(f.table(0, 0, 0, 0, 0)))
    result = calibration.analyze_tubes(zero, (F(0),)*9, {'C': _table('C'), 'E': _table('E')}, 1)
    assert result['compatible_locked_candidates'] == []
    assert result['full_family_feasibility'] == 'not_decided'
    assert result['candidate_set_exhaustive'] is False


def test_complete_matrix_comparison_does_not_erase_antihermitian_errors():
    from model.repair_identification import bridge
    first = {record: np.zeros((8, 8), complex) for record in product(range(7), repeat=2)}
    second = {record: value.copy() for record, value in first.items()}
    second[0, 0][0, 1], second[0, 0][1, 0] = 1, -1
    # Arithmetic comparison must charge the entire returned error, even when
    # it is not Hermitian; symmetrizing first would incorrectly return zero.
    assert bridge.distance(first, second) == pytest.approx(1, abs=1e-14)


class _FalseDifference(np.ndarray):
    def __sub__(self, other):
        return np.zeros(self.shape, complex)


def test_complete_matrix_comparison_rejects_an_array_that_hides_its_difference():
    from model.repair_identification import bridge
    changed = np.eye(8, dtype=complex).view(_FalseDifference)
    with pytest.raises(ValueError):
        bridge.distance({(0,): changed}, {(0,): np.zeros((8, 8), complex)})
    with pytest.raises(ValueError):
        bridge.distance({(False,): np.eye(8)}, {(0,): np.eye(8)})


def test_complete_matrix_comparison_cannot_wrap_integer_subtraction():
    from model.repair_identification import bridge
    first = {(0,): np.array([[np.iinfo(np.int64).max]], dtype=np.int64)}
    second = {(0,): np.array([[np.iinfo(np.int64).min]], dtype=np.int64)}
    try:
        distance = bridge.distance(first, second)
    except ValueError:
        return  # Rejecting non-floating comparison matrices is also valid.
    assert distance > 2**62


def test_full_source_reduction_guard_charges_spectator_leakage():
    from model.repair_identification import bridge
    matrix = np.zeros((512, 8), complex)
    matrix[:8] = np.eye(8)
    matrix[511, 0] = 1e-10
    diagnostic = bridge.reduction_diagnostics([matrix])
    assert diagnostic['stacked_leakage_norm'] == pytest.approx(1e-10, rel=1e-12)
    assert diagnostic['projection_change_bound_numeric'] >= 0.999e-10
    matrix[511, 0] = 1e-4
    with pytest.raises(ArithmeticError):
        bridge.reduction_diagnostics([matrix])
    matrix[511, 0] = float('nan')
    with pytest.raises(ValueError):
        bridge.reduction_diagnostics([matrix])
    with pytest.raises(ValueError):
        bridge.reduction_diagnostics([np.zeros((512, 8)).view(_HiddenArray)])


@pytest.fixture(scope='module')
def raw_reference_history():
    from model.repair_identification import source
    return source.continue_history('C', 'g_re', reference=True)


@pytest.mark.parametrize('mutation', ('drop_zero', 'duplicate_zero', 'bool_current', 'bool_pointer',
                                     'float_reference', 'flatten_reference', 'array_subclass',
                                     'bad_cost', 'wrong_word', 'nan_outside'))
def test_common_adapter_rejects_incomplete_or_aliased_raw_histories(raw_reference_history, mutation):
    from model.repair_identification import bridge
    history = raw_reference_history
    branches = list(history.branches)
    zero = min(range(len(branches)), key=lambda i: np.linalg.norm(branches[i].amplitude))
    index = next(i for i, branch in enumerate(branches) if branch.currents[0] == 0
                 and branch.record[0][-1] == 1 and branch.record[1] is not None)
    selected = branches[index]
    if mutation == 'drop_zero':
        branches.pop(zero)
    elif mutation == 'duplicate_zero':
        branches.append(branches[zero])
    elif mutation == 'bool_current':
        branches[index] = selected._replace(currents=(False, selected.currents[1]))
    elif mutation == 'bool_pointer':
        first, second = selected.record
        branches[index] = selected._replace(record=((*first[:-1], True), second))
    elif mutation == 'float_reference':
        history = history._replace(reference_dimension=2.0)
    elif mutation == 'flatten_reference':
        branches[index] = selected._replace(amplitude=selected.amplitude.ravel())
    elif mutation == 'array_subclass':
        branches[index] = selected._replace(amplitude=selected.amplitude.view(_HiddenArray))
    elif mutation == 'bad_cost':
        costs = dict(selected.resources)
        costs['native_repairs'] -= 1
        branches[index] = selected._replace(resources=tuple(sorted(costs.items())))
    elif mutation == 'wrong_word':
        name, digest = history.word_hashes[0]
        history = history._replace(word_hashes=((name, '0'*len(digest)), *history.word_hashes[1:]))
    else:
        changed = selected.amplitude.copy()
        changed[511, 1] = float('nan')
        branches[index] = selected._replace(amplitude=changed)
    history = history._replace(branches=tuple(branches))
    with pytest.raises(ValueError):
        bridge.coarse_raw(history)


def test_producer_rejects_a_probability_invisible_reference_state_corruption(raw_reference_history, monkeypatch):
    from model.repair_identification import experiment, source
    branches = list(raw_reference_history.branches)
    index = max((i for i, branch in enumerate(branches)
                 if np.count_nonzero(np.abs(branch.amplitude) > 1e-4) >= 2),
                key=lambda i: np.linalg.norm(branches[i].amplitude))
    branch = branches[index]
    changed = branch.amplitude.copy()
    changed.flat[int(np.argmax(np.abs(changed)))] *= -1
    assert np.vdot(changed, changed) == pytest.approx(np.vdot(branch.amplitude, branch.amplitude), abs=1e-15)
    branches[index] = branch._replace(amplitude=changed)
    corrupted = raw_reference_history._replace(branches=tuple(branches))
    monkeypatch.setattr(source, 'continue_history', lambda *args, **kwargs: corrupted)
    with pytest.raises(ArithmeticError, match='continuing source history'):
        experiment._history('C', 'g_re', True, _table('C'))


@pytest.mark.parametrize('mutation', ('drop_zero', 'swap_candidate', 'nan_outside'))
def test_probe_statistic_consumer_checks_complete_source_provenance(monkeypatch, mutation):
    from model.repair_identification import experiment, source
    execution = source.execute_probe('C', 'g_re')
    if mutation == 'drop_zero':
        branches = list(execution.branches)
        index = min(range(len(branches)), key=lambda i: np.linalg.norm(branches[i].amplitude))
        branches.pop(index)
        execution = execution._replace(branches=tuple(branches))
    elif mutation == 'swap_candidate':
        execution = execution._replace(candidate='E')
    else:
        branches = list(execution.branches)
        changed = branches[0].amplitude.copy()
        changed[511] = float('nan')
        branches[0] = branches[0]._replace(amplitude=changed)
        execution = execution._replace(branches=tuple(branches))
    monkeypatch.setattr(source, 'execute_probe', lambda *args, **kwargs: execution)
    probe = next(p for p in core.PROBES if p.name == 'g_re')
    with pytest.raises(ValueError):
        experiment._probe('C', probe, _table('C'))


@pytest.fixture(scope='module')
def retained_receipt():
    from tools import strict_json
    return strict_json.load(Path(__file__).with_name('results.json'))


def test_retained_receipt_covers_the_complete_executed_scope_and_source_bindings(retained_receipt):
    receipt = retained_receipt
    assert set(receipt['candidates']) == {'C', 'E', 'C_phase_quarter', 'E_collision_half'}
    expected_raw = {'C': 14, 'E': 26, 'C_phase_quarter': 14, 'E_collision_half': 34}
    for law, candidate in receipt['candidates'].items():
        assert [row['probe'] for row in candidate['source_probes']] == [p.name for p in core.PROBES]
        for row in candidate['source_probes']:
            assert len(row['complete_outcome_probabilities']) == 7
            assert row['raw_outcomes'] == expected_raw[law]
            assert row['maximum_branch_resources']['declared_blank_data_preparations'] == 1
    histories = receipt['continuing_histories']
    assert len(histories) == 10
    assert sum(row['untouched_reference_dimension'] == 2 for row in histories) == 5
    for row in histories:
        assert row['complete_visible_histories'] == len(row['complete_visible_probabilities']) == 49
        assert row['maximum_raw_branch_resources']['declared_blank_data_preparations'] == 1
        assert row['canonical_runtime']['checkpoint_replay_distance'] == 0
        supplied = row['canonical_runtime']['supplied_preparation_cost']
        if row['untouched_reference_dimension'] == 2:
            assert supplied['declared_blank_data_preparations'] == supplied['joint_reference_preparations'] == 1
        else:
            assert supplied == {}
    assert receipt['error_budget']['machine_roundoff_certified'] is False
    assert receipt['scope']['physical_repair_law_selected'] is False
    assert 'model/common_source/protocol.py' in receipt['sources']
    root = Path(__file__).resolve().parents[2]
    for path, digest in receipt['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest() == digest


@pytest.mark.parametrize('mutation', ('drop_abort', 'drop_history', 'change_phase_probability',
                                     'bool_count', 'float_reference', 'nan', 'source_hash',
                                     'physical_claim', 'omit_preparation', 'exact_table'))
def test_receipt_comparator_rejects_scientific_and_structural_forgeries(retained_receipt, mutation):
    from model.repair_identification import experiment
    changed = deepcopy(retained_receipt)
    probe = changed['candidates']['C']['source_probes'][0]
    history = changed['continuing_histories'][0]
    if mutation == 'drop_abort':
        probe['complete_outcome_probabilities'].pop()
    elif mutation == 'drop_history':
        history['complete_visible_probabilities'].pop()
    elif mutation == 'change_phase_probability':
        values = history['complete_visible_probabilities']
        index = max(range(len(values)), key=values.__getitem__)
        values[index] -= 1/64
        values[(index+1) % len(values)] += 1/64
    elif mutation == 'bool_count':
        probe['raw_outcomes'] = True
    elif mutation == 'float_reference':
        history['untouched_reference_dimension'] = float(history['untouched_reference_dimension'])
    elif mutation == 'nan':
        history['full_probability_mass'] = float('nan')
    elif mutation == 'source_hash':
        path = next(iter(changed['sources']))
        changed['sources'][path] = '0'*64
    elif mutation == 'physical_claim':
        changed['scope']['physical_repair_law_selected'] = True
    elif mutation == 'omit_preparation':
        row = next(row for row in changed['continuing_histories'] if row['untouched_reference_dimension'] == 2)
        row['canonical_runtime']['supplied_preparation_cost'] = {}
    else:
        changed['candidates']['E_collision_half']['table']['n'] = '0'
    with pytest.raises(ValueError):
        experiment.compare_receipt(changed, retained_receipt)


def test_real_cli_rejects_duplicate_receipt_keys_before_source_execution(tmp_path):
    forged = tmp_path/'duplicate.json'
    forged.write_text('{"schema":"oph.repair_identification.results.v1","all_pass":true,"all_pass":false}',
                      encoding='utf-8')
    root = Path(__file__).resolve().parents[2]
    # The pinned Windows interpreter has an isolated ._pth; bootstrap the
    # actual checked-out tree explicitly rather than inheriting another tree.
    script = ('import runpy,sys; from pathlib import Path; '
              'p=Path(sys.argv.pop(1)); sys.path.insert(0,str(p/"code")); '
              'sys.path.insert(0,str(p)); '
              'runpy.run_module("model.repair_identification",run_name="__main__")')
    result = subprocess.run([sys.executable, '-c', script, str(root), '--check', str(forged)],
                            cwd=root, capture_output=True, text=True, timeout=30)
    assert result.returncode != 0
    assert 'duplicate JSON key' in result.stderr
    assert '"all_pass": true' not in result.stdout
