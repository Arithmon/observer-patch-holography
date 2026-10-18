"""Independent native-branch controls for the quantum incident-clock bridge.

The oracle starts with the original integer endpoint rule and its committed
record equality.  It never uses the producer's channel matrices, square root,
coherence formula, or matrix-composition helper to predict an output.
"""
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
from math import comb
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from ..observer import native, quantum_mixing
from . import quantum_clock as q


ROOT = (0, 1, 5)
REGIONS = ((0, 1), (0, 5))
MARKS = ((0, 1, 0), (1, 0, 0), (1, 1, -1), (1, 1, 0), (1, 1, 1))
TIMEOUT = 'timeout'


@lru_cache(maxsize=1)
def _native_images():
    """Only classical endpoint execution is cached, with both original labels."""
    result = []
    for seam in quantum_mixing.seams():
        seam = tuple(seam)
        remote = not set(seam).intersection(ROOT)
        ticks = tuple(int(bool(set(seam).intersection(region))) for region in REGIONS)
        columns = []
        for port in range(12):
            state = tuple(int(i == port) for i in range(12))
            images = []
            for coin in (0, 1):
                final, (_, current) = native.native(state, seam, coin)
                assert sum(final) == 1 and set(final) <= {0, 1}
                images.append((final.index(1), current,
                               state[seam[0]]-state[seam[1]],
                               final[0]-state[0]))
            columns.append(tuple(images))
        result.append((remote, ticks, tuple(columns)))
    assert len(result) == 30 and sum(row[0] for row in result) == 18
    return tuple(result)


def _add(target, source, weight=F(1)):
    for key, value in source.items():
        target[key] += weight*value


def _step(candidate, density, *, remote):
    """One actual uniform seam/coin call, with its original hidden record."""
    outputs = {None: defaultdict(F)} if remote else {
        mark: defaultdict(F) for mark in MARKS}
    label = 1 if candidate == 'C' else 2
    for is_remote, ticks, columns in _native_images():
        if remote != is_remote:
            continue
        for (left, right), value in density.items():
            for coin in (0, 1):
                a, b = columns[left][coin], columns[right][coin]
                if a[label] != b[label]:
                    continue
                # A definite C/E hidden record cannot mix visible currents.
                assert a[1] == b[1] and a[3] == b[3]
                mark = None if remote else (*ticks, a[3])
                outputs[mark][a[0], b[0]] += value/F(60)
    return {mark: {key: value for key, value in state.items() if value}
            for mark, state in outputs.items()}


def _remote(candidate, state, calls):
    for _ in range(calls):
        state = _step(candidate, state, remote=True)[None]
    return state


def _direct_interval(candidate, initial, cap, *, effective):
    """Execute all possible remote-before/incident/remote-after lengths."""
    outputs = {mark: defaultdict(F) for mark in MARKS}
    prefixes = [initial]
    for _ in range(cap):
        prefixes.append(_step(candidate, prefixes[-1], remote=True)[None])
    if not effective:
        for prefix in prefixes[:-1]:
            for mark, state in _step(candidate, prefix, remote=False).items():
                _add(outputs[mark], state)
    else:
        for total_remote in range(cap):
            for before in range(total_remote+1):
                after = total_remote-before
                weight = F(comb(2*before, before)*comb(2*after, after),
                           4**total_remote)
                for mark, state in _step(candidate, prefixes[before], remote=False).items():
                    _add(outputs[mark], _remote(candidate, state, after), weight)
    outputs[TIMEOUT] = prefixes[-1]
    return {mark: {key: value for key, value in state.items() if value}
            for mark, state in outputs.items()}


def _column(channel, left, right):
    if left == right:
        return {(row, row): channel.population[row][left]
                for row in range(12) if channel.population[row][left]}
    value = channel.coherence[left][right]
    return {(left, right): value} if value else {}


def _apply(channel, state):
    result = defaultdict(F)
    for (left, right), value in state.items():
        _add(result, _column(channel, left, right), value)
    return {key: value for key, value in result.items() if value}


@pytest.mark.parametrize('candidate', ('C', 'E'))
def test_actual_source_instrument_on_every_operator_input(candidate):
    source = q.source(candidate)
    assert source.root == ROOT and source.probability == F(2, 5)
    assert set(dict(source.incident)) == set(MARKS)
    for left, right in product(range(12), repeat=2):
        state = {(left, right): F(1)}
        expected = _step(candidate, state, remote=False)
        assert _column(source.remote, left, right) == _step(
            candidate, state, remote=True)[None]
        for mark, channel in source.incident:
            assert _column(channel, left, right) == expected[mark]


@pytest.mark.parametrize('candidate', ('C', 'E'))
@pytest.mark.parametrize('cap', (1, 2, 3))
@pytest.mark.parametrize('effective', (False, True))
def test_capped_maps_match_raw_native_histories_on_every_input(candidate, cap, effective):
    maps = q.capped_maps(candidate, cap, effective=effective)
    assert set(maps) == {*MARKS, TIMEOUT}
    for left, right in product(range(12), repeat=2):
        expected = _direct_interval(candidate, {(left, right): F(1)}, cap,
                                    effective=effective)
        for mark, channel in maps.items():
            assert _column(channel, left, right) == expected[mark]
        if left == right:
            assert sum(value for state in expected.values()
                       for (i, j), value in state.items() if i == j) == 1
            assert sum(expected[TIMEOUT].values()) == F(3, 5)**cap


@pytest.mark.parametrize('candidate,retention', (('C', F(7, 12)), ('E', F(1, 4))))
def test_native_clock_coherence_witness_is_not_a_classical_population_test(candidate, retention):
    unit = {(0, 1): F(1)}
    remote = _step(candidate, unit, remote=True)[None]
    incident = defaultdict(F)
    for state in _step(candidate, unit, remote=False).values():
        _add(incident, state)
    assert remote == {(0, 1): F(3, 5)}
    assert incident[0, 1]/(1-remote[0, 1]) == retention
    source = q.source(candidate)
    assert source.remote.coherence[0][1] == remote[0, 1]
    assert sum(channel.coherence[0][1] for _, channel in source.incident) == incident[0, 1]


def test_classical_triangle_floor_cannot_be_reused_for_the_quantum_E_channel():
    # The three triangle seams all erase |100><010| for E, while their
    # classical population sum has a strictly positive lower bound.
    left = (1, 0, 0)
    right = (0, 1, 0)
    for seam in ((0, 1), (0, 2), (1, 2)):
        for coin in (0, 1):
            native.native(left, seam, coin)
            native.native(right, seam, coin)
            assert left[seam[0]]-left[seam[1]] != right[seam[0]]-right[seam[1]]


def test_full_binary_triangle_floor_matches_independent_native_eigenspaces():
    import sympy as sp

    certificate = q.triangle_certificate()
    assert certificate['triangle_ports'] == ROOT
    assert set(certificate['triangle_seams']) <= {tuple(edge) for edge in quantum_mixing.seams()}
    assert certificate['full_binary_C_incident_hs_lower'] == F(1, 30)
    states = tuple(product((0, 1), repeat=3))
    index = {state: i for i, state in enumerate(states)}
    for candidate in ('C', 'E'):
        matrix = sp.zeros(64)
        for left, right in product(range(8), repeat=2):
            for seam, coin in product(((0, 1), (0, 2), (1, 2)), (0, 1)):
                a, b = states[left], states[right]
                x, (_, current_a) = native.native(a, seam, coin)
                y, (_, current_b) = native.native(b, seam, coin)
                label_a = current_a if candidate == 'C' else a[seam[0]]-a[seam[1]]
                label_b = current_b if candidate == 'C' else b[seam[0]]-b[seam[1]]
                if label_a == label_b:
                    matrix[8*index[x]+index[y], 8*left+right] += 1
        # Direct exact eigenspaces are independent of the producer's
        # annihilator and spectral-projector multiplicity calculation.
        spectrum = tuple(sorted((int(value), int(size)) for value, size in matrix.eigenvals().items()))
        assert spectrum == certificate['count_spectra_with_multiplicity'][candidate]
        assert sum(size for _, size in spectrum) == 64


class _IntAlias(int):
    pass


class _StringAlias(str):
    pass


@pytest.mark.parametrize('candidate', (True, 0, None, 'c', _StringAlias('C')))
def test_source_rejects_candidate_aliases_before_cache_lookup(candidate):
    q.source('C')
    with pytest.raises(ValueError):
        q.source(candidate)


@pytest.mark.parametrize('cap', (True, 0, 33, 2.0, F(2), _IntAlias(2)))
def test_interval_rejects_noncanonical_caps_before_cache_lookup(cap):
    q.capped_maps('C', 2)
    with pytest.raises(ValueError):
        q.capped_maps('C', cap)


@pytest.mark.parametrize('effective', (0, 1, 'False', None))
def test_interval_requires_an_actual_protocol_boolean(effective):
    with pytest.raises(ValueError):
        q.capped_maps('C', 2, effective=effective)


@pytest.mark.parametrize('root', ((False, 1, 5), (0, _IntAlias(1), 5),
                                (0, 1, 1), (0, 1, 12), [0, 1, 5]))
def test_root_support_rejects_nested_aliases_and_malformed_geometry(root):
    q.source('C')
    with pytest.raises(ValueError):
        q.source('C', root=root)


def _fraction_snapshot(channel):
    return tuple(tuple((x.numerator, x.denominator) for x in row)
                 for matrix in (channel.population, channel.coherence) for row in matrix)


def test_public_source_fractions_cannot_poison_cached_source():
    first = q.source('C')
    expected = _fraction_snapshot(first.remote)
    fraction = first.remote.population[0][0]
    assert type(fraction) is F
    fraction._numerator += fraction.denominator
    assert _fraction_snapshot(q.source('C').remote) == expected


def test_public_interval_maps_and_fractions_are_detached_from_cached_proof():
    first = q.capped_maps('C', 2, effective=True)
    expected = {mark: _fraction_snapshot(channel) for mark, channel in first.items()}
    fraction = first[TIMEOUT].population[0][0]
    fraction._numerator += fraction.denominator
    first.clear()
    assert {mark: _fraction_snapshot(channel)
            for mark, channel in q.capped_maps('C', 2, effective=True).items()} == expected


@pytest.mark.parametrize('candidate', ('C', 'E'))
@pytest.mark.parametrize('effective', (False, True))
def test_largest_cap_preserves_complete_mass_and_reference_positivity(candidate, effective):
    # Exact scalar mass plus a separately assembled Choi Gram diagnostic.
    # The all-matrix-unit raw branch oracle above supplies the exact CP proof.
    maps = q.capped_maps(candidate, 32, effective=effective)
    for column in range(12):
        assert sum(channel.population[row][column] for channel in maps.values()
                   for row in range(12)) == 1
        assert sum(maps[TIMEOUT].population[row][column] for row in range(12)) == F(3, 5)**32
    for channel in maps.values():
        assert all(value >= 0 for row in channel.population for value in row)
        gram = np.array([[float(channel.population[i][i] if i == j
                                else channel.coherence[i][j])
                          for j in range(12)] for i in range(12)])
        assert np.array_equal(gram, gram.T)
        assert np.linalg.eigvalsh(gram).min() >= -1e-13


@pytest.mark.parametrize('candidate', ('C', 'E'))
@pytest.mark.parametrize('effective', (False, True))
def test_same_marginal_clocks_do_not_replace_their_joint_record(candidate, effective):
    cap = 3
    maps = q.capped_maps(candidate, cap, effective=effective)
    success = 1-F(3, 5)**cap
    joint = {(0, 0): F(0), (0, 1): F(1, 4),
             (1, 0): F(1, 4), (1, 1): F(1, 2)}
    for column in range(12):
        for pair, conditional in joint.items():
            assert sum(channel.population[row][column]
                       for mark, channel in maps.items() if mark != TIMEOUT and mark[:2] == pair
                       for row in range(12)) == success*conditional
    independent = {pair: F(3, 4)**sum(pair)*F(1, 4)**(2-sum(pair)) for pair in joint}
    assert sum(joint[pair] for pair in joint if pair[0]) == F(3, 4)
    assert sum(independent[pair] for pair in joint if pair[0]) == F(3, 4)
    assert sum(abs(joint[pair]-independent[pair]) for pair in joint)/2 == F(1, 8)


@pytest.mark.parametrize('candidate,passive,active', (
    ('C', F(49, 144), F(7, 144)), ('E', F(1, 16), -F(1, 16))))
@pytest.mark.parametrize('cap', (1, 3, 32))
@pytest.mark.parametrize('effective', (False, True))
def test_record_selected_phase_changes_the_unconditional_coherent_field(candidate, passive, active,
                                                                         cap, effective):
    maps = q.capped_maps(candidate, cap, effective=effective)
    ordinary = sum(maps[mark].coherence[0][1] for mark in MARKS)
    selected = sum((-1 if mark[:2] == (0, 1) else 1)*maps[mark].coherence[0][1]
                   for mark in MARKS)
    retained = (1-F(3, 5)**cap)**2
    assert ordinary**2 == passive*retained
    assert selected*ordinary == active*retained
    assert selected*ordinary != ordinary**2


def _initial(reference):
    positions = (0, reference+(reference-1))
    return {(left, right): F(1, 2) for left, right in product(positions, repeat=2)}


def _reference_interval(candidate, state, cap, effective, reference):
    outputs = {mark: {} for mark in (*MARKS, TIMEOUT)}
    for r, s in product(range(reference), repeat=2):
        block = {(i//reference, j//reference): value for (i, j), value in state.items()
                 if i % reference == r and j % reference == s}
        if not block:
            continue
        for mark, final in _direct_interval(candidate, block, cap, effective=effective).items():
            outputs[mark].update({(i*reference+r, j*reference+s): value
                                  for (i, j), value in final.items()})
    return outputs


def _phase_measurement(state, outcome, reference):
    def column(port):
        if outcome == 2:
            return ((port, F(1)),) if port not in (0, 1) else ()
        if port not in (0, 1):
            return ()
        sign = -1 if outcome == 1 else 1
        return ((0, F(sign**port, 2)), (1, F(sign**(port+1), 2)))

    output = defaultdict(F)
    for (i, j), value in state.items():
        for left, a in column(i//reference):
            for right, b in column(j//reference):
                output[left*reference+i % reference, right*reference+j % reference] += a*b*value
    return {key: value for key, value in output.items() if value}


def _raw_history(candidate, cap, *, effective, feedback, reference):
    """All native field/clock/ref outputs, with physical feedback and reads."""
    outputs = {}
    for first, state in _reference_interval(candidate, _initial(reference), cap,
                                             effective, reference).items():
        if first == TIMEOUT:
            outputs[(TIMEOUT,)] = state
            continue
        if feedback and first[:2] == (0, 1):
            state = {key: (-1)**((key[0]//reference == 0)+(key[1]//reference == 0))*value
                     for key, value in state.items()}
        for second, final in _reference_interval(candidate, state, cap, effective, reference).items():
            if second == TIMEOUT:
                outputs[first, TIMEOUT] = final
            else:
                for outcome in range(3):
                    outputs[first, second, ('phase', outcome)] = _phase_measurement(final, outcome, reference)
    return outputs


def _sparse(array):
    return {(i, j): value for (i, j), value in np.ndenumerate(array) if value}


@pytest.mark.parametrize('candidate', ('C', 'E'))
@pytest.mark.parametrize('effective', (False, True))
@pytest.mark.parametrize('reference', (False, True))
@pytest.mark.parametrize('cap,feedback', ((1, False), (1, True), (2, True)))
def test_complete_quantum_history_matches_independent_endpoint_recursion(candidate, effective,
                                                                         reference, cap, feedback):
    from . import quantum_experiment as experiment

    dimension = 2 if reference else 1
    expected = _raw_history(candidate, cap, effective=effective, feedback=feedback,
                            reference=dimension)
    actual = experiment.exact_history(candidate, cap, effective=effective,
                                      feedback=feedback, reference=reference)
    assert set(actual) == set(expected) and len(actual) == 81
    assert sum(len(word) == 1 for word in actual) == 1
    assert sum(len(word) == 2 for word in actual) == 5
    assert sum(len(word) == 3 for word in actual) == 75
    for word, array in actual.items():
        assert type(array) is np.ndarray and array.shape == (12*dimension, 12*dimension)
        assert _sparse(array) == expected[word]
    assert sum(np.trace(array) for array in actual.values()) == 1
    assert sum(np.trace(array) for word, array in actual.items() if len(word) < 3) == (
        1-(1-F(3, 5)**cap)**2)


def _root_trace(state, reference):
    """Partial trace in the full twelve-qubit embedding, restricted to N=1."""
    local = {port: index+1 for index, port in enumerate(ROOT)}
    output = defaultdict(F)
    for (i, j), value in state.items():
        a, b = i//reference, j//reference
        remote_a = () if a in local else (a,)
        remote_b = () if b in local else (b,)
        if remote_a == remote_b:
            output[local.get(a, 0)*reference+i % reference,
                   local.get(b, 0)*reference+j % reference] += value
    return {key: value for key, value in output.items() if value}


@pytest.mark.parametrize('reference', (1, 2))
def test_root_restriction_traces_remote_labels_instead_of_collecting_amplitudes(reference):
    from . import quantum_experiment as experiment

    # Every source/reference matrix unit, including invisible exterior
    # coherences and root/reference coherences, tests the actual partial trace.
    for i, j in product(range(12*reference), repeat=2):
        matrix = np.empty((12*reference, 12*reference), dtype=object)
        matrix.fill(F(0))
        matrix[i, j] = F(1)
        actual = experiment.root_output(matrix, reference=reference)
        assert actual.shape == (4*reference, 4*reference)
        assert _sparse(actual) == _root_trace({(i, j): F(1)}, reference)


def test_untouched_reference_separates_quantum_output_from_equal_phase_probabilities():
    from . import quantum_experiment as experiment

    outputs = {law: experiment.exact_history(law, 2, feedback=True, reference=True)
               for law in ('C', 'E')}
    assert set(outputs['C']) == set(outputs['E'])
    assert all(np.trace(outputs['C'][word]) == np.trace(outputs['E'][word])
               for word in outputs['C'])
    assert any(not np.array_equal(outputs['C'][word], outputs['E'][word])
               for word in outputs['C'])
    for law in ('C', 'E'):
        # A local analyzer on one half of the Bell input is not a joint Bell
        # analyzer. Its equal outcome probabilities do not erase the retained
        # conditional correlations with the untouched reference.
        contrast = sum((1 if word[-1] == ('phase', 0) else -1)*np.trace(state)
                       for word, state in outputs[law].items()
                       if len(word) == 3 and word[-1] != ('phase', 2))
        assert contrast == 0


@pytest.mark.parametrize('candidate', ('C', 'E'))
@pytest.mark.parametrize('effective', (False, True))
def test_channel_application_preserves_reference_indices_and_nonhermitian_parts(candidate, effective):
    maps = q.capped_maps(candidate, 2, effective=effective)
    for (left, right), r, s in product(((0, 1), (1, 0), (2, 2), (0, 5), (2, 11)),
                                      range(2), range(2)):
        unit = {(2*left+r, 2*right+s): F(1)}
        expected = _reference_interval(candidate, unit, 2, effective, 2)
        density = np.empty((24, 24), dtype=object)
        density.fill(F(0))
        density[2*left+r, 2*right+s] = F(1)
        for mark, channel in maps.items():
            assert _sparse(q.apply(channel, density, reference_dimension=2)) == expected[mark]


class _ArrayAlias(np.ndarray):
    def tolist(self):
        raise AssertionError('an untrusted subclass method must not run')


class _FractionAlias(F):
    pass


@pytest.mark.parametrize('bad_entry', (True, 0.0, float('nan'), _IntAlias(0), _FractionAlias(0)))
def test_exact_apply_rejects_numeric_aliases_in_original_entries(bad_entry):
    density = np.full((12, 12), F(0), dtype=object)
    density[0, 0] = bad_entry
    with pytest.raises(ValueError):
        q.apply(q.source('C').remote, density)


def test_exact_apply_rejects_array_subclass_before_dispatching_user_methods():
    density = np.full((12, 12), F(0), dtype=object).view(_ArrayAlias)
    with pytest.raises(ValueError):
        q.apply(q.source('C').remote, density)


@pytest.mark.parametrize('reference', (True, 0, 5, 2.0, _IntAlias(2)))
def test_exact_apply_requires_canonical_reference_dimensions(reference):
    density = np.full((12, 12), F(0), dtype=object)
    with pytest.raises(ValueError):
        q.apply(q.source('C').remote, density, reference_dimension=reference)


class _HiddenNonzero(np.ndarray):
    @property
    def flat(self):
        return np.zeros(self.size, dtype=self.dtype).flat


@pytest.mark.parametrize('forged_side', ('actual', 'target'))
def test_binary_error_bound_rejects_subclasses_that_hide_original_entries(forged_side):
    from . import quantum_experiment as experiment

    actual = np.zeros((2, 2), dtype=complex)
    target = np.full((2, 2), F(0), dtype=object)
    if forged_side == 'actual':
        actual[0, 0] = 1
        actual = actual.view(_HiddenNonzero)
    else:
        target[0, 0] = F(1)
        target = target.view(_HiddenNonzero)
    with pytest.raises(ValueError):
        experiment.binary_output_bound({(TIMEOUT,): actual}, {(TIMEOUT,): target})


@pytest.mark.parametrize('bad', (True, 0.0, _IntAlias(0), _FractionAlias(0)))
def test_binary_error_bound_rejects_original_exact_target_aliases(bad):
    from . import quantum_experiment as experiment

    actual = np.zeros((2, 2), dtype=float)
    target = np.full((2, 2), F(0), dtype=object)
    target[0, 1] = bad
    with pytest.raises(ValueError):
        experiment.binary_output_bound({(TIMEOUT,): actual}, {(TIMEOUT,): target})


@pytest.mark.parametrize('actual', (np.zeros((2, 2), dtype=np.float32),
                                   np.zeros((2, 2), dtype=object),
                                   np.full((2, 2), float('nan')),
                                   np.full((2, 2), float('inf'))))
def test_binary_error_bound_requires_finite_original_binary64_arrays(actual):
    from . import quantum_experiment as experiment

    target = np.full((2, 2), F(0), dtype=object)
    with pytest.raises(ValueError):
        experiment.binary_output_bound({(TIMEOUT,): actual}, {(TIMEOUT,): target})


def test_binary_error_bound_preserves_imaginary_nonhermitian_error_and_missing_records():
    from . import quantum_experiment as experiment

    actual = np.zeros((2, 2), dtype=complex)
    actual[0, 1] = 1+2j
    target = np.full((2, 2), F(0), dtype=object)
    assert experiment.binary_output_bound({(TIMEOUT,): actual}, {(TIMEOUT,): target}) == F(3, 2)
    with pytest.raises(ValueError):
        experiment.binary_output_bound({}, {(TIMEOUT,): target})


@pytest.mark.parametrize('candidate,effective', (('C', False), ('E', True)))
def test_common_runtime_and_serialized_continuation_keep_every_reference_record(candidate, effective):
    from . import quantum_experiment as experiment

    result = experiment.execute_history(candidate, 2, effective=effective,
                                        feedback=True, reference=True)
    expected = _raw_history(candidate, 2, effective=effective, feedback=True, reference=2)
    assert set(result['outputs']) == set(expected)
    assert result['checkpoint_bytes'] > 0
    assert result['history'].reference_dimension == 2
    assert len(result['history'].branches) == 108
    discrepancy = 0.
    for word, actual in result['outputs'].items():
        target = np.zeros((24, 24), dtype=float)
        for (i, j), value in expected[word].items():
            target[i, j] = float(value)
        discrepancy += np.linalg.svd(actual-target, compute_uv=False).sum()/2
    assert discrepancy < 1e-12
    assert 0 <= result['error_bound'] <= experiment.NUMERIC_ALLOWANCE
    for branch in result['history'].branches:
        first, second, phase = branch.record
        if first == 5 and (second != 5 or phase != 2):
            assert not np.any(branch.density)
            continue
        if second == 5 and phase != 2:
            assert not np.any(branch.density)
            continue
        resources = dict(branch.resources)
        assert resources['native_attempt_budget'] == (2 if first == 5 else 4)
        assert resources['interval_requests'] == (1 if first == 5 else 2)
        assert resources.get('local_phase_controls', 0) == int(first == 0)
        assert resources.get('phase_analyzers', 0) == int(first != 5 and second != 5)


@pytest.mark.parametrize('candidate', ('C', 'E'))
@pytest.mark.parametrize('effective', (False, True))
def test_actual_preanalyzer_branch_retains_C_reference_entanglement(candidate, effective):
    from . import quantum_experiment as experiment
    from ..primitive_observer import runtime

    mark = (1, 1, 0)
    first = _reference_interval(candidate, _initial(2), 2, effective, 2)[mark]
    terminal = _reference_interval(candidate, first, 2, effective, 2)[mark]
    local = _root_trace(terminal, 2)
    coherence = F(25, 1152)*(1-F(3, 5)**2)**2 if candidate == 'C' else F(0)
    assert local.get((2, 5), F(0)) == local.get((5, 2), F(0)) == coherence
    assert local.get((3, 3), F(0)) == local.get((4, 4), F(0)) == 0
    # After reference transposition this exact principal minor is negative
    # for C. The final local phase measurement is not used for this witness.
    assert -coherence**2 < 0 if candidate == 'C' else coherence == 0

    program = experiment.runtime_program(candidate, 2, effective=effective, feedback=True)
    prepared = np.zeros((24, 24))
    for (i, j), value in _initial(2).items():
        prepared[i, j] = float(value)
    prefix = runtime.execute(program, initial=prepared, stop_after=3)
    selected = next(branch for branch in prefix.branches if branch.record == (3, 3))
    target = np.zeros((24, 24))
    for (i, j), value in terminal.items():
        target[i, j] = float(value)
    assert np.linalg.norm(selected.density-target, ord='fro') < 1e-13


def _conditional_record(candidate, state, seam, coin, label, current, reference):
    """Replay one supplied native audit record without the sampler's groups."""
    images = []
    for port in range(12):
        loads = tuple(int(i == port) for i in range(12))
        target, (_, moved) = native.native(loads, seam, coin)
        hidden = moved if candidate == 'C' else loads[seam[0]]-loads[seam[1]]
        images.append((target.index(1), hidden, moved))
    result = defaultdict(F)
    for (left, right), value in state.items():
        a, b = images[left//reference], images[right//reference]
        if a[1:] == b[1:] == (label, current):
            result[a[0]*reference+left % reference, b[0]*reference+right % reference] += value
    probability = sum(value for (i, j), value in result.items() if i == j)
    assert 0 < probability <= 1
    return {key: value/probability for key, value in result.items() if value}


def _audit_trajectory(row):
    from ..observer.scalar_response.experiment import canonical

    reference = row['reference_dimension']
    state = _initial(reference)
    public = []
    controls = 0
    seams = {tuple(edge) for edge in quantum_mixing.seams()}
    for edge, coin, label, current in row['native_repairs']:
        seam = tuple(edge)
        assert seam in seams and type(coin) is int and coin in (0, 1)
        state = _conditional_record(row['law'], state, seam, coin, label, current, reference)
        if set(seam).intersection(ROOT):
            mark = (*(int(bool(set(seam).intersection(region))) for region in REGIONS),
                    current*(int(seam[1] == 0)-int(seam[0] == 0)))
            public.append(mark)
            if len(public) == 1 and row['feedback'] and mark[:2] == (0, 1):
                # Every following remote selective branch commutes with this
                # phase, so it can be replayed at the actual incident boundary.
                state = {key: (-1)**((key[0]//reference == 0)+(key[1]//reference == 0))*value
                         for key, value in state.items()}
                controls += 1
    if row['public_word'][-1] == TIMEOUT:
        assert len(public) < 2
        public.append(TIMEOUT)
    else:
        assert len(public) == 2
        outcome = row['public_word'][-1]
        assert outcome[0] == 'phase' and outcome[1] in (0, 1, 2)
        state = _phase_measurement(state, outcome[1], reference)
        probability = sum(value for (i, j), value in state.items() if i == j)
        assert probability > 0
        state = {key: value/probability for key, value in state.items()}
        public.append(outcome)
    assert public == row['public_word']
    assert controls == row['local_phase_controls']
    assert row['native_calls'] == len(row['native_repairs']) <= 2*row['cap']
    reduced = _root_trace(state, reference)
    rows = [[str(reduced.get((i, j), F(0))) for j in range(4*reference)]
            for i in range(4*reference)]
    digest = hashlib.sha256(canonical([[('terminal',), rows]]).encode()).hexdigest()
    assert digest == row['terminal_root_reference_sha256']


@pytest.mark.parametrize('candidate', ('C', 'E'))
@pytest.mark.parametrize('effective', (False, True))
@pytest.mark.parametrize('reference', (False, True))
def test_seeded_source_calls_replay_private_records_without_replacing_quantum_state(candidate,
                                                                                   effective, reference):
    from . import quantum_experiment as experiment

    for cap in (1, 3):
        exact = experiment.exact_history(candidate, cap, effective=effective,
                                         feedback=True, reference=reference)
        for seed in range(8):
            row = experiment.source_trajectory(candidate, cap, seed, effective=effective,
                                                feedback=True, reference=reference)
            _audit_trajectory(row)
            assert np.trace(exact[tuple(row['public_word'])]) > 0
    first = experiment.source_trajectory(candidate, 3, 0, effective=effective, reference=reference)
    assert first == experiment.source_trajectory(candidate, 3, 0, effective=effective, reference=reference)


@pytest.mark.parametrize('effective', (False, True))
def test_unfiltered_short_cap_examples_include_both_timeout_prefixes_and_completion(effective):
    from . import quantum_experiment as experiment

    rows = [experiment.source_trajectory('C', 1, seed, effective=effective) for seed in range(8)]
    assert any(row['public_word'] == [TIMEOUT] and row['native_calls'] == 1 for row in rows)
    assert any(len(row['public_word']) == 2 and row['public_word'][-1] == TIMEOUT
               and row['native_calls'] == 2 for row in rows)
    assert any(len(row['public_word']) == 3 and row['native_calls'] == 2 for row in rows)


@pytest.mark.parametrize('seed', (True, -1, 1.0, _IntAlias(1)))
def test_native_trajectory_requires_an_original_nonnegative_seed(seed):
    from . import quantum_experiment as experiment

    with pytest.raises(ValueError):
        experiment.source_trajectory('C', 2, seed)


@pytest.mark.parametrize('arguments', ({'cap': True}, {'effective': 1},
                                      {'feedback': 0}, {'reference': 1}))
def test_native_trajectory_rejects_aliased_protocol_parameters(arguments):
    from . import quantum_experiment as experiment

    parameters = {'cap': 2, 'seed': 0, **arguments}
    with pytest.raises(ValueError):
        experiment.source_trajectory('C', **parameters)


def test_rational_sampler_uses_exact_mass_boundaries_including_zero_outcomes():
    from . import quantum_experiment as experiment

    class Draw:
        def __init__(self, selected):
            self.selected = selected

        def randrange(self, denominator):
            assert denominator == 6
            return self.selected

    weights = (F(0), F(1, 6), F(1, 3), F(1, 2), F(0))
    assert [experiment._draw(Draw(i), weights) for i in range(6)] == [1, 2, 2, 3, 3, 3]


def test_cap_sixteen_C_passive_runtime_preserves_exact_zero_outputs_and_full_mass():
    from . import quantum_experiment as experiment

    result = experiment.execute_history('C', 16, feedback=False, reference=False)
    assert result['error_bound'] <= experiment.NUMERIC_ALLOWANCE
    assert len(result['outputs']) == 81
    assert sum(np.trace(state) for state in result['exact'].values()) == 1
    for word, target in result['exact'].items():
        if not np.any(target):
            assert np.linalg.norm(result['outputs'][word], ord='fro') < 1e-13


@pytest.fixture(scope='module')
def replayed_quantum_receipt():
    from . import quantum_experiment as experiment
    from ..observer.scalar_response.experiment import canonical, compare_receipt
    from tools import strict_json

    path = Path(__file__).with_name('quantum_results.json')
    retained = strict_json.loads(path.read_text(encoding='utf-8'))
    recomputed = strict_json.loads(canonical(experiment.run()))
    compare_receipt(retained, recomputed)
    return retained


def test_retained_complete_receipt_replays_actual_source_and_boundaries(replayed_quantum_receipt):
    receipt = replayed_quantum_receipt
    assert len(receipt['experiments']) == 8
    assert all(len(row['experiments']) == 2 for row in receipt['experiments'])
    for row in receipt['experiments']:
        for result in row['experiments']:
            assert result['complete_output_labels'] == 81
            assert F(result['timeout_probability']) > 0
            assert sum(map(F, result['phase_outcomes_unconditioned']))+F(result['timeout_probability']) == 1
    assert any(row['public_word'][-1] == TIMEOUT for row in receipt['native_source_trajectories'])
    assert any(row['negative_minor_certifies_entanglement']
               for row in receipt['pre_analyzer_reference_witnesses'])


@pytest.mark.parametrize('mutation', ('drop_timeout', 'change_phase_with_same_mass',
                                     'change_reference_state_with_same_traces',
                                     'numeric_bool_alias', 'numeric_float_alias'))
def test_receipt_rejects_missing_outputs_and_quantum_changes_with_unchanged_source_bindings(
        replayed_quantum_receipt, mutation):
    from ..observer.scalar_response.experiment import compare_receipt

    original = replayed_quantum_receipt
    altered = deepcopy(original)
    result = altered['experiments'][0]['experiments'][0]
    if mutation == 'drop_timeout':
        del result['timeout_probability']
    elif mutation == 'change_phase_with_same_mass':
        field = result['phase_outcomes_unconditioned']
        field[0], field[1] = str(F(field[0])+F(1, 1000)), str(F(field[1])-F(1, 1000))
        assert sum(map(F, field)) == sum(map(F, original['experiments'][0]['experiments'][0][
            'phase_outcomes_unconditioned']))
    elif mutation == 'change_reference_state_with_same_traces':
        result = next(row for row in altered['experiments']
                      if row['untouched_reference_dimension'] == 2)['experiments'][0]
        result['root_reference_output_sha256'] = '0'*64
    elif mutation == 'numeric_bool_alias':
        altered['contract']['conserved_total'] = True
    else:
        altered['bounds']['maximum_native_calls'] = float(altered['bounds']['maximum_native_calls'])
    assert altered['source_sha256'] == original['source_sha256']
    with pytest.raises(ValueError, match='Receipt'):
        compare_receipt(altered, original)


def test_retained_receipt_pipeline_rejects_duplicate_keys_and_nonfinite_json(replayed_quantum_receipt):
    from ..observer.scalar_response.experiment import canonical, compare_receipt
    from tools import strict_json

    payload = canonical(replayed_quantum_receipt)
    duplicate = payload.replace('"schema":', '"schema": "forged", "schema":', 1)
    with pytest.raises(ValueError):
        strict_json.loads(duplicate)
    count = replayed_quantum_receipt['bounds']['maximum_native_calls']
    nonfinite = payload.replace(f'"maximum_native_calls": {count}',
                                '"maximum_native_calls": NaN', 1)
    assert nonfinite != payload
    # The shared parser owns duplicate-key rejection; the exact receipt
    # comparison owns scalar-type and nonfinite-value rejection.
    with pytest.raises(ValueError):
        compare_receipt(strict_json.loads(nonfinite), replayed_quantum_receipt)


def test_real_quantum_cli_rejects_changed_reference_output_with_original_mass_and_hashes(
        replayed_quantum_receipt, tmp_path):
    from ..observer.scalar_response.experiment import canonical

    altered = deepcopy(replayed_quantum_receipt)
    row = next(row for row in altered['experiments'] if row['untouched_reference_dimension'] == 2)
    row['experiments'][0]['root_reference_output_sha256'] = '0'*64
    path = tmp_path/'forged-quantum-receipt.json'
    path.write_text(canonical(altered), encoding='utf-8')
    root = Path(__file__).resolve().parents[2]
    bootstrap = (
        "import sys; from pathlib import Path; root=Path(sys.argv.pop(1)); "
        "sys.path.insert(0,str(root/'code')); sys.path.insert(0,str(root)); "
        "from runpy import run_module; "
        "run_module('model.observer_correspondence.quantum_experiment',run_name='__main__')")
    result = subprocess.run([sys.executable, '-c', bootstrap, str(root), '--check', str(path)],
                            cwd=root, capture_output=True, text=True, timeout=180,
                            env={**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
    assert result.returncode != 0
    assert 'Receipt mismatch' in result.stderr
    assert 'root_reference_output_sha256' in result.stderr
