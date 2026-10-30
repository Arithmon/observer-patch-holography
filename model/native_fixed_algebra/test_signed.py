"""Actual signed channels, exact algebra and independent reference controls."""
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest

from . import fixed, signed as s
from ..exchange_sector import statistics
from ..local_repair_family import local_bias as local, mixing
from ..observer import record_recovery as records


def choices(shift=0):
    return tuple(F((i+shift) % 3, 2) for i in range(30))


@pytest.mark.parametrize('character', [1, -1])
@pytest.mark.parametrize('weight', [F(0), F(1, 2), F(1)])
def test_both_complete_lifts_preserve_the_same_faithful_reference(character, weight):
    parameters = local.Parameters(port_weights=(25, 25)+(8,)*10, c_weight=weight,
                                  repair_probability=F(1, 16))
    result = s.faithful_channel_certificate(parameters, character)
    assert result['raw_slots_including_zero_and_hold'] == 721
    assert result['same_faithful_diagonal_reference']
    assert result['complete']
    assert F(result['smallest_reference_weight']) > 0


@pytest.mark.parametrize('character', [1, -1])
def test_effects_and_actual_signed_transition_tree_generate_full_matrix_algebra(character):
    rows = s.branches(mixing.canonical_parameters(), character, choices())
    result = s.algebra_certificate(rows)
    assert len(set(result['effect_signatures'])) == 66
    assert result['commutant_dimension'] == 1
    assert result['generated_star_algebra_dimension'] == 66**2
    assert len(result['transition_forest']) == 65
    lookup = {row.label: row for row in rows}
    for edge in result['transition_forest']:
        row = lookup[edge['branch']]
        assert row.probability > 0
        assert (edge['source'], edge['target'], edge['nonzero_phase']) in row.columns
        assert {edge['parent'], edge['child']} == {edge['source'], edge['target']}
    assert result['scalar_record_gauge_equivalence_claimed'] is False


@pytest.mark.parametrize('character', [1, -1])
@pytest.mark.parametrize('ports', [(1,)*12, (25, 25)+(8,)*10, (1, 4)*6])
def test_signed_seam_blocks_are_exact_positive_GNS_contractions(character, ports):
    result = s.seam_GNS_certificate(local.Parameters(port_weights=ports), character)
    assert result['matrix_units_checked'] == 60*66**2
    assert result['weighted_self_adjoint_positive_contractions']
    assert (result['negative_matrix_entries_retained'] > 0) == (character == -1)


@pytest.mark.parametrize('character', [1, -1])
def test_signed_dephasing_comparison_and_exact_uniform_gap(character):
    comparison = s.extra_dephasing_certificate(character)
    assert comparison['matrix_unit_coin_seam_checks'] == 60*66**2
    assert comparison['C_minus_E_diagonal_nonnegative_quadratic_form']
    assert comparison['C_minus_E_positive_map_claimed'] is False
    result = s.uniform_gap_certificate(character)
    assert result['absolute_coherence_row_sum_upper'] == '5/6'
    assert result['entrywise_absolute_domination_verified']
    assert result['signed_channel_replaced_by_absolute_values'] is False
    assert F(result['uniform_gap_lower']) == F(fixed.uniform_gap_certificate()['uniform_gap_lower'])


def _bitstring_phase(source, target):
    before = sum(1 << i for i, value in enumerate(source) if value)
    after = sum(1 << i for i, value in enumerate(target) if value)
    if before == after:
        return 1
    origin = (before & ~after).bit_length()-1
    destination = (after & ~before).bit_length()-1
    phase = (-1)**((before & ((1 << origin)-1)).bit_count())
    before ^= 1 << origin
    phase *= (-1)**((before & ((1 << destination)-1)).bit_count())
    return phase


def _oracle(parameters, matrix, character, law_choices):
    states = local.basis()
    index = {state: i for i, state in enumerate(states)}
    output = (1-float(parameters.repair_probability))*matrix.copy()
    for seam, sw, cw in zip(local.SEAMS, parameters.seam_weights, law_choices):
        for law, choice in (('C', cw), ('E', 1-cw)):
            for coin in (0, 1):
                probability = (parameters.repair_probability*F(sw, sum(parameters.seam_weights))
                               *choice*local.coin_probability(parameters, seam, coin))
                if not probability:
                    continue
                images = [records.recorded_step(law, state, seam, coin) for state in states]
                for i, (left, left_record) in enumerate(images):
                    si = _bitstring_phase(states[i], left) if character == -1 else 1
                    for j, (right, right_record) in enumerate(images):
                        if left_record == right_record:
                            sj = _bitstring_phase(states[j], right) if character == -1 else 1
                            output[i, j] += float(probability)*si*sj*matrix[index[left], index[right]]
    return output


@pytest.mark.parametrize('character', [1, -1])
def test_actual_signed_quantum_map_matches_independent_endpoint_and_Fock_oracle(character):
    parameters = local.Parameters(port_weights=(1, 4)*6, repair_probability=F(3, 4))
    rng = np.random.default_rng(672)
    matrix = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    actual = s.apply_heisenberg(parameters, matrix, character, choices())
    np.testing.assert_allclose(actual, _oracle(parameters, matrix, character, choices()), rtol=2e-14, atol=2e-14)


@pytest.mark.parametrize('character', [1, -1])
def test_signed_dirichlet_identity_and_near_conservation_bound(character):
    parameters = mixing.canonical_parameters()
    rng = np.random.default_rng(891)
    matrix = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    result = s.dirichlet_control(parameters, matrix, character, choices())
    assert result['dirichlet_trace_pairing'] == pytest.approx(result['half_weighted_commutator_sum'], rel=2e-13)
    gamma = float(s.gap_bound(parameters, character, choices()))
    assert gamma*result['weighted_variance'] <= result['dirichlet_trace_pairing']+1e-12
    assert gamma**2*result['weighted_variance'] <= result['weighted_fixed_point_residual_squared']+1e-12
    fixed_identity = s.dirichlet_control(parameters, np.eye(66), character, choices())
    assert fixed_identity['half_weighted_commutator_sum'] == 0


def _norm(matrix, pi):
    return np.sqrt(np.sum(pi*np.sum(abs(matrix)**2, axis=0)))


def test_exogenous_character_switching_with_untouched_reference_contracts_to_reference_algebra():
    ports = (25, 25)+(8,)*10
    parameters = [local.Parameters(port_weights=ports, seam_weights=tuple(1+(i+k) % 16 for i in range(30)),
                                   repair_probability=F(1, 2)) for k in range(3)]
    pi = np.array(list(map(float, local.reference(parameters[0]))))
    rng = np.random.default_rng(108)
    blocks = rng.normal(size=(2, 2, 66, 66))+1j*rng.normal(size=(2, 2, 66, 66))
    for a in range(2):
        for b in range(2):
            blocks[a, b] -= np.dot(pi, np.diag(blocks[a, b]))*np.eye(66)
    bound = sum(_norm(block, pi)**2 for row in blocks for block in row)
    for k, (parameter, character) in enumerate(zip(parameters, (-1, 1, -1))):
        mixture = choices(k)
        for a in range(2):
            for b in range(2):
                blocks[a, b] = s.apply_heisenberg(parameter, blocks[a, b], character, mixture)
                assert abs(np.dot(pi, np.diag(blocks[a, b]))) < 1e-12
        bound *= (1-float(s.gap_bound(parameter, character, mixture)))**2
        assert sum(_norm(block, pi)**2 for row in blocks for block in row) <= bound+1e-12
    np.testing.assert_allclose(s.apply_heisenberg(parameters[0], 2j*np.eye(66), -1), 2j*np.eye(66), atol=1e-14)


@pytest.mark.parametrize('character', [1, -1])
def test_removed_seams_restore_held_port_charge_outside_exploration_hypothesis(character):
    rows = tuple(replace(row, probability=row.probability*F(30, 25))
                 for row in s.branches(local.Parameters(repair_probability=F(1)), character)
                 if row.label == ('hold',) or 0 not in row.label[1])
    algebra = s.algebra_certificate(rows)
    assert sorted(map(len, algebra['components'])) == [11, 55]
    assert algebra['commutant_dimension'] == 2
    charge = np.array([state[0] for state in local.basis()])
    actual = np.zeros(66)
    for row in rows:
        for i, j, sign in row.columns:
            actual[i] += float(row.probability)*sign**2*charge[j]
    np.testing.assert_allclose(actual, charge, atol=1e-14)


def test_equal_gap_does_not_mean_the_two_quantum_channels_are_equal():
    matrix = np.zeros((66, 66), complex)
    matrix[statistics.INDEX[(0, 1)], statistics.INDEX[(0, 6)]] = 1
    parameters = local.Parameters(c_weight=F(1))
    assert np.linalg.norm(s.apply_heisenberg(parameters, matrix, 1)
                          -s.apply_heisenberg(parameters, matrix, -1)) > .01
    assert s.gap_bound(parameters, 1) == s.gap_bound(parameters, -1)


def test_exact_branchwise_similarity_is_not_a_positive_channel():
    result = s.similarity_certificate()
    assert result['complete_raw_slots'] == 720
    assert result['branch_matrix_units_checked'] > 100000
    assert result['same_spectrum_for_every_fixed_selector']
    assert result['S_is_positive_or_a_physical_channel'] is False
    assert result['negative_normalized_output_expectation'] == '-1/3'
    rng = np.random.default_rng(309)
    matrix = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    parameters = mixing.canonical_parameters()
    transformed = s.matrix_unit_transform(s.apply_heisenberg(
        parameters, s.matrix_unit_transform(matrix), 1, choices()))
    np.testing.assert_allclose(transformed, s.apply_heisenberg(parameters, matrix, -1, choices()), atol=2e-14)
    np.testing.assert_array_equal(s.matrix_unit_transform(s.matrix_unit_transform(matrix)), matrix)
    vector = np.zeros(66)
    addresses = [statistics.INDEX[pair] for pair in ((0, 1), (0, 2), (1, 2))]
    vector[addresses] = 1/np.sqrt(3)
    output = s.matrix_unit_transform(np.outer(vector, vector))
    probe = np.zeros(66)
    probe[addresses] = np.array([1, -1, 1])/np.sqrt(3)
    assert (probe@output@probe).real == pytest.approx(-1/3, abs=1e-15)


def test_changing_endpoint_weights_does_not_preserve_the_common_reference():
    initial = local.reference(local.Parameters())
    output = [F()]*66
    for row in s.branches(mixing.canonical_parameters(), -1):
        for i, j, sign in row.columns:
            output[j] += row.probability*sign**2*initial[i]
    assert tuple(output) != initial


def test_character_mutation_fails_even_when_absolute_effects_stay_the_same(monkeypatch):
    original = statistics.branch
    def changed(*args):
        row = original(*args)
        if row.inputs and any(i != j for i, j in zip(row.inputs, row.outputs)):
            row = replace(row, signs=(-row.signs[0], *row.signs[1:]))
        return row
    monkeypatch.setattr(statistics, 'branch', changed)
    with pytest.raises(ArithmeticError, match='selected exchange character'):
        s.faithful_channel_certificate(local.Parameters(), -1)


@pytest.mark.parametrize('damage', ['zero_phase', 'duplicate_output', 'negative_weight'])
def test_algebra_certificate_rejects_invalid_signed_partial_injections(damage):
    rows = list(s.branches(local.Parameters(), -1))
    index = next(i for i, row in enumerate(rows) if len(row.columns) >= 2 and row.probability)
    row = rows[index]
    columns = list(row.columns)
    if damage == 'zero_phase':
        i, j, _ = columns[0]
        columns[0] = (i, j, 0)
    elif damage == 'duplicate_output':
        i, _, sign = columns[1]
        columns[1] = (i, columns[0][1], sign)
    rows[index] = replace(row, columns=tuple(columns),
                          probability=F(-1) if damage == 'negative_weight' else row.probability)
    with pytest.raises(ValueError):
        s.algebra_certificate(rows)


@pytest.mark.parametrize('character', [True, False, 0, 1.0, '-1', None])
def test_invalid_character_is_not_silently_accepted(character):
    with pytest.raises(ValueError):
        s.branches(local.Parameters(), character)
    with pytest.raises(ValueError):
        s.gap_bound(local.Parameters(), character)
