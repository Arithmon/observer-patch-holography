from fractions import Fraction as F

import numpy as np
import pytest

from . import fixed
from ..local_repair_family import local_bias as local, mixing
from ..observer import record_recovery as records


def choices():
    return tuple(F(index % 3, 2) for index in range(30))


@pytest.mark.parametrize('weight', [F(0), F(1, 2), F(1)])
def test_actual_complete_channel_has_same_faithful_reference(weight):
    parameters = local.Parameters(port_weights=(25, 25)+(8,)*10,
                                  seam_weights=tuple(1+i % 16 for i in range(30)),
                                  c_weight=weight, repair_probability=F(1, 16))
    result = fixed.faithful_channel_certificate(parameters)
    assert result['dimension'] == 66
    assert result['complete'] and result['faithful_invariant_reference']
    assert F(result['smallest_reference_weight']) > 0


def test_independent_seam_law_choices_preserve_reference_and_positive_diagonal_span():
    parameters = mixing.canonical_parameters()
    result = fixed.faithful_channel_certificate(parameters, choices())
    assert result['faithful_invariant_reference']
    positive = fixed.positive_diagonal_Kraus_span(parameters, choices())
    assert positive['minimum_positive_diagonal'] > 0
    assert positive['coherent_sum_implemented'] is False


def test_extra_E_measurement_is_a_quadratic_form_comparison_on_all_matrix_units():
    result = fixed.extra_dephasing_certificate()
    assert result['matrix_unit_coin_seam_checks'] == 30*2*66**2
    assert result['additional_C_matrix_unit_terms'] > 0
    assert result['C_minus_E_positive_map_claimed'] is False


@pytest.mark.parametrize('ports', [(1,)*12, (25, 25)+(8,)*10, (1, 4)*6])
def test_all_actual_seam_GNS_blocks_are_exact_positive_contractions(ports):
    result = fixed.seam_GNS_certificate(local.Parameters(port_weights=ports))
    assert result['matrix_units_checked'] == 60*66**2
    assert result['two_dimensional_blocks'] > 0
    assert result['weighted_self_adjoint_positive_contractions']


def _native_heisenberg_oracle(parameters, observable, law_choices):
    states = local.basis()
    index = {state: i for i, state in enumerate(states)}
    output = (1-float(parameters.repair_probability))*observable.copy()
    for edge, seam_weight, c_weight in zip(local.SEAMS, parameters.seam_weights, law_choices, strict=True):
        for law, law_weight in (('C', c_weight), ('E', 1-c_weight)):
            for coin in (0, 1):
                weight = (parameters.repair_probability*F(seam_weight, sum(parameters.seam_weights))
                          * law_weight*local.coin_probability(parameters, edge, coin))
                if not weight:
                    continue
                # Independent endpoint calls; no compiled partial maps.
                images = [records.recorded_step(law, state, edge, coin) for state in states]
                for i, (left, left_record) in enumerate(images):
                    for j, (right, right_record) in enumerate(images):
                        if left_record == right_record:
                            output[i, j] += float(weight)*observable[index[left], index[right]]
    return output


def test_complex_observable_pullback_matches_actual_endpoint_oracle():
    parameters = local.Parameters(port_weights=(1, 4)*6, repair_probability=F(3, 4))
    rng = np.random.default_rng(197)
    observable = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    actual = fixed.apply_heisenberg(parameters, observable, choices())
    expected = _native_heisenberg_oracle(parameters, observable, choices())
    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-14)


def test_weighted_commutator_identity_and_near_fixed_residual_bound():
    parameters = mixing.canonical_parameters()
    rng = np.random.default_rng(31)
    observable = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    result = fixed.dirichlet_control(parameters, observable, choices())
    assert result['dirichlet_trace_pairing'] == pytest.approx(
        result['half_weighted_commutator_sum'], rel=2e-13)
    gamma = float(fixed.gap_bound(parameters, choices()))
    assert gamma*result['weighted_variance'] <= result['dirichlet_trace_pairing']+1e-12
    assert gamma**2*result['weighted_variance'] <= result['weighted_fixed_point_residual_squared']+1e-12
    identity = fixed.dirichlet_control(parameters, np.eye(66), choices())
    assert abs(identity['dirichlet_trace_pairing']) < 1e-14
    assert identity['half_weighted_commutator_sum'] == 0


def test_constant_expectation_in_one_stationary_state_is_not_operator_conservation():
    parameters = mixing.canonical_parameters()
    pi = np.array(list(map(float, local.reference(parameters))))
    observable = np.diag([state[0] for state in local.basis()])
    evolved = fixed.apply_heisenberg(parameters, observable, choices())
    assert np.dot(pi, np.diag(evolved)) == pytest.approx(np.dot(pi, np.diag(observable)))
    assert np.linalg.norm(evolved-observable) > .1


def test_removing_incident_seams_restores_a_nonconstant_conserved_charge():
    # This controlled exception is outside positive all-seam exploration.
    observable = np.diag([state[0] for state in local.basis()]).astype(complex)
    held = tuple(edge for edge in local.SEAMS if 0 not in edge)
    output = np.zeros_like(observable)
    for edge in held:
        for coin, _, pairs in local.raw_maps('C', edge):
            if pairs:
                inputs, outputs = zip(*pairs, strict=True)
                output[np.ix_(inputs, inputs)] += observable[np.ix_(outputs, outputs)]/(2*len(held))
    np.testing.assert_allclose(output, observable, atol=1e-14)
    assert len(set(np.diag(observable))) == 2


def _weighted_norm(matrix, pi):
    return np.sqrt(np.sum(pi*np.sum(abs(matrix)**2, axis=0)))


def test_exogenous_different_selectors_contract_one_common_reference():
    ports = (25, 25)+(8,)*10
    parameters = [local.Parameters(port_weights=ports, seam_weights=tuple(1+(i+j) % 16 for i in range(30)),
                                   repair_probability=F(1, 2)) for j in range(3)]
    rng = np.random.default_rng(61)
    matrix = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    pi = np.array(list(map(float, local.reference(parameters[0]))))
    matrix -= np.dot(pi, np.diag(matrix))*np.eye(66)
    bound = _weighted_norm(matrix, pi)
    for step, parameter in enumerate(parameters):
        mixture = tuple(F((i+step) % 3, 2) for i in range(30))
        matrix = fixed.apply_heisenberg(parameter, matrix, mixture)
        bound *= 1-float(fixed.gap_bound(parameter, mixture))
        assert _weighted_norm(matrix, pi) <= bound+1e-12
        assert abs(np.dot(pi, np.diag(matrix))) < 1e-12


def test_reference_operator_extension_centers_on_reference_algebra_not_scalars():
    parameters = mixing.canonical_parameters()
    pi = np.array(list(map(float, local.reference(parameters))))
    rng = np.random.default_rng(74)
    blocks = rng.normal(size=(66, 2, 66, 2))+1j*rng.normal(size=(66, 2, 66, 2))
    for a in range(2):
        for b in range(2):
            matrix = blocks[:, a, :, b]
            matrix -= np.dot(pi, np.diag(matrix))*np.eye(66)
    original = sum(_weighted_norm(blocks[:, a, :, b], pi)**2 for a in range(2) for b in range(2))
    final = sum(_weighted_norm(fixed.apply_heisenberg(parameters, blocks[:, a, :, b]), pi)**2
                for a in range(2) for b in range(2))
    gamma = float(fixed.gap_bound(parameters))
    assert final <= (1-gamma)**2*original+1e-12
    untouched = fixed.apply_heisenberg(parameters, 2j*np.eye(66))
    np.testing.assert_allclose(untouched, 2j*np.eye(66), atol=1e-14)


def test_rational_gap_is_positive_and_family_floor_is_conservative():
    base = fixed.uniform_gap_certificate()
    assert base['offdiagonal_matrix_units'] == 66*65
    assert F(base['uniform_gap_lower']) == min(F(base['exact_population_minorization'])/20, F(1, 6))
    extremal = local.Parameters(port_weights=(1, 4)*6, seam_weights=(1,)+(16,)*29,
                                repair_probability=F(1, 16))
    assert fixed.gap_bound(extremal) >= F(base['uniform_gap_lower'])/10240 > 0


@pytest.mark.parametrize('weights', [(F(0),)*29, (0,)*30, (F(-1),)*30, (F(2),)*30])
def test_invalid_per_seam_selection_fails_admission(weights):
    with pytest.raises(ValueError):
        fixed.branches(local.Parameters(), weights)


def test_mutated_original_fraction_does_not_enter_source_certificate():
    parameters = local.Parameters(repair_probability=F(3, 4))
    object.__setattr__(parameters.repair_probability, '_denominator', 0)
    with pytest.raises(ValueError):
        fixed.faithful_channel_certificate(parameters)
