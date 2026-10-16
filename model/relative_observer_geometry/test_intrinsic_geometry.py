from copy import deepcopy
from fractions import Fraction as F

import pytest

from ..observer_reconstruction import native_geometry as ng
from . import intrinsic_geometry as g


@pytest.fixture(scope='module')
def calibrations():
    records = ng.execute_calibration()
    seams = tuple(sorted({tuple(row['seam']) for row in records}))
    return records, seams, g.execute_biased_calibration(seams)


def test_structural_geometry_survives_preparation_dependent_schedule_and_deterministic_coins(calibrations):
    uniform, _, biased = calibrations
    before = ng.recover(uniform)
    after = g.recover_structure(biased)
    assert len(biased) == 540
    assert {row['coin'] for row in biased if row['preparation_port'] == 0} == {0}
    assert after['seam_exposures_exact'][0] != after['seam_exposures_exact'][1]
    assert after['laplacian_exact'] == before['laplacian_exact']
    assert after['gram_qsqrt5'] == before['gram_qsqrt5']
    assert after['rank_three_structural_gram_certified']


def test_missing_exposure_does_not_inherit_the_icosahedral_certificate(calibrations):
    _, seams, _ = calibrations
    result = g.recover_structure(ng.execute_calibration({seams[0]: F(1)}))
    assert not result['rank_three_structural_gram_certified']


@pytest.mark.parametrize('coin_value', [None, 0, 1])
def test_one_native_shot_per_setting_recovers_gram_without_probabilities(calibrations, coin_value):
    records, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams, coin_value)
    assert len(shots) == 360
    assert all('probability' not in row for row in shots)
    result = g.recover_from_shots(shots)
    assert result['gram_qsqrt5'] == ng.recover(records)['gram_qsqrt5']
    assert not result['probability_law_oracle_used']
    assert result['conditional_current_variance_zero_under_native_model']
    assert not result['laboratory_noise_or_preparation_error_bounded']


def test_unequal_finite_repetition_counts_leave_structural_recovery_unchanged(calibrations):
    _, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams)
    repeated = [row for index, row in enumerate(shots) for _ in range(1+index % 5)]
    assert g.recover_from_shots(shots)['gram_qsqrt5'] == g.recover_from_shots(repeated)['gram_qsqrt5']


def test_finite_shots_reject_a_wrong_native_current_or_probability_oracle(calibrations):
    _, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams)
    changed = deepcopy(shots)
    changed[0]['current'] = 0
    with pytest.raises(ValueError, match='even-load'):
        g.recover_from_shots(changed)
    with pytest.raises(ValueError, match='not probability laws'):
        g.recover_from_shots(calibrations[0])
    with pytest.raises(ValueError, match='inconsistent currents'):
        g.recover_from_shots(shots+changed[:1])


def test_finite_shot_missing_coverage_cannot_recover_the_full_metric(calibrations):
    _, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams)
    with pytest.raises(ValueError, match='same positive seam support'):
        g.recover_from_shots(shots[1:])
    incomplete = [row for row in shots if tuple(row['seam']) != seams[0]]
    assert not g.recover_from_shots(incomplete)['rank_three_structural_gram_certified']


def test_relabeling_covariance_without_source_coordinates(calibrations):
    _, _, records = calibrations
    permutation = (7, 2, 10, 0, 5, 11, 8, 4, 1, 9, 3, 6)
    relabelled = deepcopy(records)
    for row in relabelled:
        row['preparation_port'] = permutation[row['preparation_port']]
        a, b = [permutation[x] for x in row['seam']]
        row['seam'] = sorted((a, b))
        if a > b:
            row['current'] *= -1
            row['coin'] = 1-row['coin']
    first, second = g.recover_structure(records), g.recover_structure(relabelled)
    assert second['rank_three_structural_gram_certified']
    for i in range(12):
        for j in range(12):
            assert first['gram_qsqrt5'][i][j] == second['gram_qsqrt5'][permutation[i]][permutation[j]]


@pytest.mark.parametrize('mutate', [
    lambda rows: rows.append(deepcopy(rows[0])),
    lambda rows: rows[0].update(probability='0'),
    lambda rows: rows[0].update(current=0),
    lambda rows: rows[0].update(coin=True),
    lambda rows: rows[0].update(preparation_port=False),
    lambda rows: rows[0].update(seam=[1, 0]),
    lambda rows: rows[0].update(hidden_loads=[]),
    lambda rows: rows.pop(),
])
def test_corrupted_or_incomplete_calibration_rejected(calibrations, mutate):
    rows = deepcopy(calibrations[0])
    mutate(rows)
    with pytest.raises(ValueError):
        g.recover_structure(rows)


def test_preparation_specific_missing_support_rejected(calibrations):
    rows = deepcopy(calibrations[0])
    seam = rows[0]['seam']
    rows = [row for row in rows if row['preparation_port'] != 0 or row['seam'] != seam]
    for row in rows:
        if row['preparation_port'] == 0:
            row['probability'] = '1/58'
    with pytest.raises(ValueError, match='same positive seam support'):
        g.recover_structure(rows)


def test_doubling_any_one_seam_splits_slowest_native_response_from_three_to_two(calibrations):
    records, seams, _ = calibrations
    laplacian = [[F(x) for x in row] for row in g.recover_structure(records)['laplacian_exact']]
    # Check all edges, without relying on an assumed symmetry action.
    for seam in seams:
        result = g.scheduling_control(laplacian, seam)
        assert result['smallest_positive_weighted_laplacian_multiplicity_exact'] == 2
        assert result['normalized_long_attempt_response_gram_rank'] == 2
        assert result['structural_gram_rank_after_exposure_removal'] == 3


def test_full_response_inverse_metric_has_exact_rank_eleven_and_matches_poisson_response(calibrations):
    records, _, _ = calibrations
    laplacian = [[F(x) for x in row] for row in g.recover_structure(records)['laplacian_exact']]
    result = g.inverse_response_metric(laplacian)
    assert result['centered_gram_rank_exact'] == 11
    assert result['squared_distance_shell_counts'] == {'11/30': 30, '7/15': 30, '1/2': 6}
    gram = [[F(x) for x in row] for row in result['gram_exact']]
    resistance = [[F(x) for x in row] for row in result['squared_euclidean_distance_exact']]
    for a in range(12):
        for b in range(a+1, 12):
            potential = [row[a]-row[b] for row in gram]
            forcing = [sum(laplacian[i][j]*potential[j] for j in range(12)) for i in range(12)]
            assert forcing == [F(i == a)-F(i == b) for i in range(12)]
            assert potential[a]-potential[b] == resistance[a][b]


def test_exact_quadratic_field_rank_handles_a_nonrational_pivot():
    z, r = (F(0), F(0)), (F(0), F(1))
    assert g.qrank([[r, z], [z, r]]) == 2
    assert g.qrank([[r, r], [r, r]]) == 1


def test_complete_single_carrier_receipt(calibrations):
    result = g.run()
    assert result['all_pass']
    assert result['same_structural_gram_as_uniform_protocol_exact']
    assert result['finite_shot_recovery']['native_branch_calls'] == 360
    assert not result['scope']['complete_ensemble_probability_law_needed_for_even_preparation_structure']
    assert not result['scope']['federation_needed_for_finite_carrier_result']
    assert not result['scope']['physical_length_time_or_lab_calibration_established']
