"""Independent exchange algebra, log bounds and complete readout controls."""
from fractions import Fraction as F
from itertools import combinations
import math

import numpy as np
import pytest
import sympy as sp

from . import construction as c


def test_local_hard_core_operator_is_identity_minus_actual_endpoint_swap():
    kernel = c.source_kernel()
    index = {state: i for i, state in enumerate(kernel.states)}
    for seam in kernel.seams:
        swap = np.zeros((66, 66), dtype=np.int64)
        for i, state in enumerate(kernel.states):
            changed = list(state)
            a, b = seam
            changed[a], changed[b] = changed[b], changed[a]
            swap[index[tuple(changed)], i] = 1
        term = np.asarray(c.seam_term(seam), dtype=np.int64)
        np.testing.assert_array_equal(term, np.eye(66, dtype=np.int64)-swap)
        np.testing.assert_array_equal(term@term, 2*term)
        np.testing.assert_array_equal(term, term.T)
    g = np.asarray(c.generator(), dtype=np.int64)
    np.testing.assert_array_equal(g.sum(axis=1), np.zeros(66))
    np.testing.assert_array_equal(60*np.eye(66, dtype=np.int64)-g, kernel.counts)


def test_source_density_interaction_cannot_be_deleted_at_fixed_hopping():
    result = c.local_field_audit()
    assert result['deleted_density_one_attempt_total_mass_from_uniform'] == '65/66'
    assert result['overlapping_seam_commutator_nonzero']
    assert result['disjoint_seam_commutator'] == '0'
    assert not result['deleted_density_interaction_constant_vector_is_ground']


def test_native_field_terms_are_local_on_record_recovered_finite_metric():
    result = c.spatial_locality_audit()
    assert result['rank_three_unit_gram_certified']
    assert result['native_exchange_terms_checked'] == 30
    assert result['pair_shell_counts'] == {'edge': 30, 'other': 36}
    assert result['every_native_exchange_term_is_geometrically_edge_local']
    assert result['native_generator_is_exact_sum_of_local_edge_terms']
    assert result['logarithm_can_have_longer_range_terms']
    assert not result['physical_length_scale_or_continuum_claimed']


def test_exact_log_polynomial_against_independent_symbolic_matrix_calculus():
    d = sp.Matrix(c.generator())/60
    direct = d+d*d/2+d*d*d/3
    assert sp.Matrix(c.log_polynomial(3)) == direct
    assert c.log_remainder(0) == F(1, 2)
    assert c.log_remainder(16) == F(1, 1463588514)
    zero = c.log_polynomial(0)
    assert all(x == 0 for row in zero for x in row)
    assert all(type(x) is F for row in c.log_polynomial(16) for x in row)


def test_log_tail_bounds_scalar_spectrum_independently_at_high_precision():
    from decimal import Decimal, localcontext
    with localcontext() as context:
        context.prec = 100
        for denominator in (3, 4, 7):
            x = Decimal(1)/denominator
            for order in (0, 1, 2, 4, 8, 16, 24):
                polynomial = sum((x**k/k for k in range(1, order+1)), Decimal(0))
                remainder = c.log_remainder(order)
                difference = -(1-x).ln()-polynomial
                assert 0 <= difference <= Decimal(remainder.numerator)/remainder.denominator


def test_logarithm_certifies_longer_range_entry_without_floating_sign_decision():
    result = c.logarithm_audit()
    witness = result['nonlocal_logarithm_witness']
    i, j = witness['state_indices']
    assert i != j and c.generator()[i][j] == 0
    assert F(witness['H_entry_strictly_positive_lower']) > 0
    assert c.log_polynomial(16)[i][j]-c.log_remainder(16) == F(witness['H_entry_strictly_positive_lower'])
    assert result['faithful_transfer_reference']['trace_K'] == '56'


@pytest.mark.parametrize('ports,witness', [
    ((0, 0, 0, 0), '-200297/101616768000'),
    ((0, 1, 0, 1), '-3701/4064670720'),
    ((0, 3, 0, 3), '-4537/20323353600'),
])
def test_four_record_moments_match_independent_weighted_transfer_products(ports, witness):
    kernel = c.source_kernel()
    counts = np.asarray(kernel.counts, dtype=object)
    score = []
    for port in ports:
        marks = [np.asarray(m, dtype=object) for m in c.currents.marked_counts(kernel, port)]
        score.append(marks[2]-marks[0])
    def moment(indices):
        row = np.ones(66, dtype=object)
        for i in range(4):
            row = row@(score[i] if i in indices else counts)
        return F(int(sum(row)), 66*60**4)
    result = c.current_moments(ports)
    assert result['complete_word_count'] == 81
    assert result['connected_four_point'] == witness
    assert F(result['four_point']) == moment(range(4))
    for i, j in combinations(range(4), 2):
        assert F(result['pair_covariances'][f'{i},{j}']) == moment((i, j))
    assert result['four_point'] != result['Gaussian_matching_two_point_prediction']


def test_complete_native_CP_readout_keeps_all_states_records_costs_and_resume():
    result = c.execute_currents((0, 1, 0, 1), checkpoint=True)
    assert len(result['complete_record_law']) == 81
    assert sum(F(row['exact_probability']) for row in result['complete_record_law']) == 1
    assert result['full_terminal_dimension'] == 66
    assert result['full_terminal_density_residual_diagnostic'] < 1e-12
    assert result['serialized_checkpoint_complete_distance_diagnostic'] < 1e-12
    zero = [row for row in result['complete_record_law'] if F(row['exact_probability']) == 0]
    assert len(zero) == 10
    assert all(row['executed_probability_diagnostic'] == 0 for row in zero)


def test_port_integration_does_not_grant_initial_hidden_loads():
    law, denominator = c.exact_current_law((0, 0, 0, 0))
    assert len(law) == 81 and len(law[1, 1, 1, 1]) == 66
    assert sum(law[1, 1, 1, 1]) == 0
    assert denominator == 855360000
    # Current history (0,0,0,0) leaves both possible initial occupations.
    kernel = c.source_kernel()
    assert {state[0] for state, n in zip(kernel.states, law[0, 0, 0, 0], strict=True) if n} == {0, 1}


def test_one_read_checkpoint_preserves_supported_short_program():
    result = c.execute_currents((0,), checkpoint=True)
    assert len(result['complete_record_law']) == 3
    assert result['serialized_checkpoint_complete_distance_diagnostic'] == 0


def test_three_class_reconstruction_is_minimal_for_complete_current_bit_words():
    result = c.record_spectral_reduction()
    assert result['minimal_linear_record_dimension'] == 3
    assert result['all_one_word_Hankel_determinant'] == '1/9702990000'
    reduction = c.clock.lumping(c.source_kernel())
    p, d = sp.Matrix(reduction['matrix']), sp.diag(*reduction['stationary'])
    root = sp.diag(*(sp.sqrt(value) for value in reduction['stationary']))
    symmetric = root*p*root.inv()
    assert symmetric == symmetric.T
    assert symmetric.trace() == sp.Rational(8, 3)
    assert (d*p-p.T*d).is_zero_matrix
    # A marked transition insertion is not a positive measurement effect.
    marked = root*sp.Matrix(reduction['marked'][1])*root.inv()
    assert marked.det() < 0


@pytest.mark.parametrize('bad', [True, -1, 25, 2., '2', None])
def test_log_orders_reject_undeclared_values_and_alias_types(bad):
    for operation in (c.log_polynomial, c.log_remainder):
        with pytest.raises(ValueError):
            operation(bad)


@pytest.mark.parametrize('bad', [(), (True,), (-1,), (12,), (1.,), (0,)*5, [0, 1]])
def test_current_ports_fail_closed(bad):
    for operation in (c.exact_current_law, c.current_program, c.current_moments):
        with pytest.raises(ValueError):
            operation(bad)


def test_seams_and_four_point_domain_fail_closed():
    for bad in ((0, 0), (0, True), (0, 1.), [0, 1]):
        with pytest.raises(ValueError):
            c.seam_term(bad)
    with pytest.raises(ValueError):
        c.current_moments((0, 1))
    with pytest.raises(ValueError):
        c.execute_currents((0,), checkpoint=1)
