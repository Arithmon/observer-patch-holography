from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations

import pytest
import sympy as sp

from . import event_geometry as e


@pytest.fixture(scope='module')
def receipt():
    return e.run()


@pytest.mark.parametrize('delta,pair', [(1, (0, 1)), (0, (0, 2))])
def test_clock_renewal_factorization_and_exact_moments(delta, pair):
    parameters = e.renewal_parameters(delta)
    z = sp.Symbol('z')
    endpoint = parameters['terminal_B_probability']
    success = parameters['geometric_success_probability']
    pgf = (1-endpoint+endpoint*z)*success/(1-(1-success)*z)
    assert sp.simplify(pgf-(5-delta+delta*z)/(10-delta-(5-delta)*z)) == 0
    mean = sp.diff(pgf, z).subs(z, 1)
    variance = sp.diff(pgf, z, 2).subs(z, 1)+mean-mean**2
    assert mean == parameters['mean'] == 1
    assert variance == parameters['variance']
    for argument in (F(0), F(1, 3), F(1)):
        for n in (1, 2, 5):
            assert e.clock_pgf(*pair, argument, n) == F(pgf.subs(z, argument)**n)


def test_finite_bins_keep_the_infinite_count_tail():
    adjacent, overflow = e.count_law(0, 1, maximum=2)
    assert adjacent == (F(4, 9), F(25, 81), F(100, 729))
    assert overflow == F(80, 729)
    assert sum(adjacent)+overflow == 1
    law, tail = e.count_law(0, 2, count=2, maximum=4)
    assert law == tuple(F(k+1, 2**(k+2)) for k in range(5))
    assert tail > 0 and sum(law)+tail == 1


def test_clock_only_reconstruction_and_covariant_port_relabelling():
    records = e.clock_records()
    laplacian, gram = e.recover(records)
    assert gram.rank() == 3
    assert (gram*gram-4*gram).applyfunc(sp.simplify) == sp.zeros(12)
    permutation = (7, 2, 10, 0, 5, 11, 8, 4, 1, 9, 3, 6)
    renamed = [{'ports': [permutation[a], permutation[b]],
                'zero_B_before_next_A': row['zero_B_before_next_A']}
               for row in records for a, b in [row['ports']]]
    new_laplacian, new_gram = e.recover(renamed)
    for a, b in combinations(range(12), 2):
        assert new_laplacian[permutation[a], permutation[b]] == laplacian[a, b]
        assert sp.simplify(new_gram[permutation[a], permutation[b]]-gram[a, b]) == 0


def test_finite_clock_sample_classifier_has_a_joint_unconditional_budget():
    intervals = 427680
    rows = [{'ports': row['ports'], 'intervals': intervals,
             'zero_counts': int(intervals*F(row['zero_B_before_next_A']))}
            for row in e.clock_records()]
    inferred = e.recover_from_counts(rows)
    assert inferred['unit_gram'] == e.recover(e.clock_records())[1]
    assert inferred['unconditional_rejection_or_wrong_graph_probability_bound'] == F(1, 20)
    assert e.geometry_error_bound(intervals) == F(1, 20)
    rows[0]['zero_counts'] = intervals
    with pytest.raises(ValueError):
        e.recover_from_counts(rows)


def test_all_native_field_outcomes_and_waiting_tails_match(receipt):
    assert receipt['independent_full_state_comparisons'] == 9
    assert receipt['native_endpoint_branches_per_pair'] == 3960
    for row in receipt['joint_native_field_predictions']:
        mass = sum(F(outcome['weight']) for outcome in row['joint_bit_generating_weights'])
        assert mass == e.clock_pgf(*row['ports'], F(row['B_count_PGF_argument']), row['A_ticks'])


def test_two_local_bits_and_two_A_ticks_resolve_all_three_geometric_shells(receipt):
    control = receipt['same_clock_identity_repair_control']
    assert control['native_joint_11_by_shell'] == {'1': '14/407', '2': '7/407', '3': '58/3663'}
    assert control['two_tick_joint_11_by_shell'] == {'1': '1/66', '2': '1/66', '3': '1/66'}
    # Equal single-bit marginals and equal clock laws do not fix the joint field.
    for pair in ((0, 1), (0, 2), (0, 3)):
        native = e.predicted_field(*pair)
        identity = e.identity_field_control(*pair)
        assert sum(native.values()) == sum(identity.values()) == 1
        assert sum(native[1, bit] for bit in (0, 1)) == F(1, 6)
        assert sum(identity[1, bit] for bit in (0, 1)) == F(1, 6)
        assert native != identity


def test_clock_geometry_is_not_invariant_under_a_changed_scheduler():
    # Doubling seam(0,1) keeps independent scheduling but changes incidence law.
    # The resulting P(B=0 before A) is4/10, outside the uniform-source class.
    rows = deepcopy(e.clock_records())
    rows[0]['zero_B_before_next_A'] = '2/5'
    with pytest.raises(ValueError, match='outside'):
        e.recover(rows)


@pytest.mark.parametrize('change', [
    lambda rows: rows.pop(),
    lambda rows: rows.__setitem__(1, rows[0]),
    lambda rows: rows[0].__setitem__('seam', [0, 1]),
    lambda rows: rows[0].__setitem__('ports', [False, 1]),
    lambda rows: rows[0].__setitem__('zero_B_before_next_A', '8/18'),
    lambda rows: rows[0].__setitem__('zero_B_before_next_A', F(4, 9)),
])
def test_malformed_or_extra_access_clock_records_are_rejected(change):
    records = deepcopy(e.clock_records()); change(records)
    with pytest.raises(ValueError):
        e.recover(records)


def test_scalar_inputs_reject_boolean_float_and_corrupted_rationals():
    for value in (True, 1.0, -1, 12):
        with pytest.raises(ValueError):
            e.clock_pgf(value, 2, F(1))
    for z in (True, 0.5, F(-1), F(2)):
        with pytest.raises(ValueError):
            e.predicted_field(0, 1, z)
    corrupted = F(1, 2); corrupted._denominator = 4
    corrupted._numerator = 2
    with pytest.raises(ValueError):
        e.clock_pgf(0, 1, corrupted)
    for delta in (True, -1, 2, F(1)):
        with pytest.raises(ValueError):
            e.renewal_parameters(delta)
