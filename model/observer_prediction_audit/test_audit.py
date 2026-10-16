"""Independent controls for null identities, data custody and trial design."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from tools import strict_json
from . import readout, prospective
from .construction import FOLDER, compare, run


def test_original_binary_effect_enclosure_and_hostile_array_boundary():
    identity = np.eye(2, dtype=complex)
    zero = np.zeros((2, 2), complex)
    half = identity/math.sqrt(2)
    menu = {'real_zero': (identity, zero), 'imag_zero': (half, half),
            'imag_plus': (identity, zero), 'imag_minus': (zero, identity)}
    result = readout.operator_audit(menu, F(1, 10**10))
    assert result['original_binary_Kraus_entries_certified']
    corrupt = dict(menu); corrupt['imag_minus'] = (identity, zero)
    with pytest.raises(ArithmeticError):
        readout.operator_audit(corrupt, F(1, 10**10))
    class Hostile(np.ndarray):
        @property
        def ndim(self):
            raise AssertionError('subclass attribute dispatched')
    with pytest.raises(ValueError):
        readout._matrix(identity.view(Hostile))
    with pytest.raises(ValueError):
        readout._matrix(np.zeros((0, 0), complex))


def test_fresh_source_controls_include_nontrivial_signed_response():
    result = readout.execute(strict_json.load(FOLDER/'protocol.json'))
    assert result['dimension'] == 8
    row = result['executed_complete_laws']['real_neighbour_superposition']
    assert abs(row['imag_plus']['mean']) > .1
    assert abs(row['imag_plus']['mean']+row['imag_minus']['mean']) < 1e-10
    for control in result['executed_complete_laws'].values():
        for name, record in control.items():
            assert len(record['probabilities']) == 2
            assert abs(sum(record['probabilities'])-1) < 1e-10
            assert min(record['probabilities']) >= 0
        assert abs(control['real_zero']['mean']-1) < 1e-10
        assert abs(control['imag_zero']['mean']) < 1e-10


def test_all_empirical_rows_and_projection_floor_are_retained():
    result = readout.compare_data()
    assert len(result['all_origin_measurements']) == 46
    assert len(result['nominal_opposite_coordinate_pairs']) == 240
    assert result['origin_statistics']['real_origin_observed_opposite_frequencies'] == ['1/500', '0', '3/500']
    for row in result['nominal_opposite_coordinate_pairs']:
        a, b = F(row['negative_mean']), F(row['positive_mean'])
        projected = (a-b)/2
        squared = (a-projected)**2+(b+projected)**2
        assert squared == F(row['least_squared_spin_residual_if_exact_opposites'])
        # Any alternative odd prediction has exactly this larger squared error.
        alternate = F(3, 7)
        assert (a-alternate)**2+(b+alternate)**2 == squared+2*(alternate-projected)**2
    floors = result['all_707_previously_scored_positive_duration_means']
    assert .00629 < floors['minimum_binary_probability_rmse_from_origins_alone'] < .00630
    assert .03419 < floors['minimum_binary_probability_rmse_with_nominal_antisymmetry'] < .03420
    assert result['statistical_or_hardware_error_inferred'] is False


@pytest.mark.parametrize('bad', ('nan', 'inf', '1.01', '-1.001', '1/0'))
def test_bad_means_rejected(bad):
    with pytest.raises(ValueError):
        readout._mean(bad)


def test_exact_complete_source_laws_and_affinity_design():
    c, e = prospective.exact_law('C'), prospective.exact_law('E')
    assert sum(c) == sum(e) == 1
    assert c[0] == e[0] == F(28, 33)
    assert c[1] == e[1] == 0
    assert c[2]/F(5, 66) == F(1921, 2592)
    assert e[2]/F(5, 66) == F(12401, 20000)
    assert sum(abs(a-b) for a, b in zip(c, e))/2 == F(12259, 668250)
    design = prospective.trial_design(c, e)
    assert design['fixed_total_trials'] == 3579
    assert design['fixed_total_trials']*F(design['log_inverse_affinity_lower']) >= F(design['log_inverse_error_upper'])
    exact_affinity = sum(math.sqrt(float(a*b)) for a, b in zip(c, e))
    assert abs(exact_affinity-float(F(design['affinity_upper']))) < 1e-14
    assert exact_affinity**design['fixed_total_trials'] < .01
    assert float(F(design['affinity_upper']))**(design['fixed_total_trials']-1) > .01
    with pytest.raises(ValueError):
        prospective.trial_design(c, c)


def test_exact_future_decision_uses_failures_and_rejects_missing_or_impossible_records():
    c, e = prospective.exact_law('C'), prospective.exact_law('E')
    # Complete rounded count fixtures for software only, not sampled evidence.
    n = prospective.trial_design(c, e)['fixed_total_trials']
    for law, expected in ((c, 'C'), (e, 'E')):
        counts = [int(n*p) for p in law]
        counts[0] += n-sum(counts)
        assert prospective.classify(tuple(counts), c, e, n) == expected
    assert prospective.classify((n-1,1,0,0,0,0), c, e, n) == 'outside_both_declared_models'
    assert prospective.classify((n,0,0,0,0,0), c, e, n) == 'C'
    with pytest.raises(ValueError):
        prospective.classify((True,0,0,0,0,0), c, e, 1)
    with pytest.raises(ValueError):
        prospective.classify((0,0,1,0,0,0), c, e, n)


def test_two_trial_likelihood_errors_are_bounded_by_product_affinity():
    from itertools import product
    c, e = prospective.exact_law('C'), prospective.exact_law('E')
    errors = {'C': F(0), 'E': F(0)}
    for word in product(range(6), repeat=2):
        counts = tuple(word.count(index) for index in range(6))
        selected = prospective.classify(counts, c, e, 2)
        if selected == 'E':
            errors['C'] += c[word[0]]*c[word[1]]
        if selected == 'C':
            errors['E'] += e[word[0]]*e[word[1]]
    bound = F(prospective.trial_design(c, e)['affinity_upper'])**2
    assert all(error <= bound for error in errors.values())


def test_complete_run_is_json_safe_and_replay_keeps_exact_bound_and_data_strings():
    result = run()
    assert json.loads(json.dumps(result)) == result
    altered = json.loads(json.dumps(result))
    altered['prospective_native_discrimination']['trial_design']['individual_error_upper'] = '1/99'
    with pytest.raises((ValueError, AssertionError)):
        compare(altered, result)
    altered = json.loads(json.dumps(result))
    altered['empirical_null_comparison']['all_origin_measurements'][0]['observed_mean'] = '0'
    with pytest.raises((ValueError, AssertionError)):
        compare(altered, result)
