"""Independent checks of publication ingestion and the observation split."""
from __future__ import annotations

import hashlib
import json
import shutil

import numpy as np
import pytest

from .data import FOLDER, load_grid, load_observations, orient, split_quadrant
from .stats import calibrate_axes, predict, score


@pytest.fixture
def grids():
    return (np.loadtxt(FOLDER / 'data/cf_0ms_real.csv', delimiter=','),
            np.loadtxt(FOLDER / 'data/cf_0ms_imag.csv', delimiter=','))


def test_deposit_matches_independent_published_checksums():
    published = {
        'cf_0ms_real.csv': 'f130544185e52fd4763eaf7e6a8c379e',
        'cf_0ms_imag.csv': 'fa478867a6d579284fea4836c0b9f187',
        'SOURCE_README.md': '04e5a8d33f1e12c73ae6e9c012d6929e',
    }
    observations = load_observations()
    for name, checksum in published.items():
        assert hashlib.md5((FOLDER / 'data' / name).read_bytes()).hexdigest() == checksum
    assert observations['provenance']['accepted_counts_included'] is False
    assert observations['provenance']['raw_shot_records_included'] is False


def test_publication_redundancy_and_noisy_central_line(grids):
    real, imag = grids
    assert np.array_equal(real, real[::-1, :])
    assert np.array_equal(real, real[:, ::-1])
    assert np.array_equal(imag, imag[:, ::-1])
    assert np.array_equal(imag[:10], -imag[20:10:-1])
    assert np.max(np.abs(imag[10])) > .1
    metadata = orient(real, imag)
    assert metadata['imaginary_coordinate_array_axis'] == 0
    assert metadata['coordinate_direction_known'] is False
    assert metadata['coordinate_extent_resolved'] is False
    assert metadata['mirror_copies_are_additional_observations'] is False


def test_split_has_disjoint_unique_symmetry_representatives(grids):
    real, imag = grids
    result = split_quadrant(real, imag)
    calibration, validation = result['calibration'], result['validation']
    assert len(calibration) == 21
    assert len(validation) == 100
    assert set(calibration).isdisjoint(validation)
    assert set(calibration) | set(validation) == {
        (i, j) for i in range(11) for j in range(11)}
    assert all(i == 0 or j == 0 for i, j in calibration)
    assert all(i > 0 and j > 0 for i, j in validation)
    for (i, j), value in (calibration | validation).items():
        assert value == complex(real[10+j, 10+i], imag[10+j, 10+i])
    # A named off-axis value catches swapped coordinates and sign changes.
    assert validation[1, 1] == pytest.approx(.432 + .656j, abs=1e-15)
    assert result['reported_origin'] == pytest.approx(.984 + .112j, abs=1e-15)
    assert calibration[0, 0] == result['reported_origin']


def test_transposing_deposit_changes_metadata_but_not_observations(grids):
    real, imag = grids
    ordinary = split_quadrant(real, imag)
    transposed = split_quadrant(real.T, imag.T)
    assert transposed['orientation']['imaginary_coordinate_array_axis'] == 1
    assert transposed['calibration'] == ordinary['calibration']
    assert transposed['validation'] == ordinary['validation']


def test_central_imaginary_noise_is_never_replaced_with_exact_zero(grids):
    real, imag = grids
    changed = imag.copy()
    changed[10, 8] = changed[10, 12] = .35
    observations = split_quadrant(real, changed)
    assert observations['calibration'][2, 0].imag == .35
    assert observations['reported_origin'].imag == imag[10, 10]


def test_negative_row_convention_conjugates_off_axis_observations(grids):
    real, imag = grids
    ordinary = split_quadrant(real, imag)
    reversed_rows = split_quadrant(real[::-1], imag[::-1])
    for key, value in ordinary['validation'].items():
        assert reversed_rows['validation'][key] == value.conjugate()
    assert reversed_rows['orientation']['coordinate_direction_known'] is False


def test_ambiguous_real_only_symmetry_does_not_choose_an_orientation(grids):
    real, _ = grids
    with pytest.raises(ValueError, match='unique'):
        orient(real, np.zeros_like(real))


def test_a_broken_mirror_is_not_treated_as_an_independent_extra_trial(grids):
    real, imag = grids
    changed = real.copy()
    changed[0, 1] += .01
    with pytest.raises(ValueError, match='unique'):
        orient(changed, imag)


@pytest.mark.parametrize('invalid', [np.nan, np.inf, 1.001, -1.001])
def test_non_spin_entries_rejected_on_load_and_orientation(tmp_path, grids, invalid):
    real, imag = grids
    changed = real.copy()
    changed[0, 0] = invalid
    path = tmp_path / 'bad.csv'
    np.savetxt(path, changed, delimiter=',')
    with pytest.raises(ValueError):
        load_grid(path)
    with pytest.raises(ValueError):
        orient(changed, imag)


def test_wrong_shape_rejected(tmp_path, grids):
    real, imag = grids
    path = tmp_path / 'wrong_shape.csv'
    np.savetxt(path, real[:20], delimiter=',')
    with pytest.raises(ValueError):
        load_grid(path)
    with pytest.raises(ValueError):
        orient(real[:20], imag)


def test_mutated_bytes_fail_even_if_csv_values_are_unchanged(tmp_path):
    copied = tmp_path / 'published'
    shutil.copytree(FOLDER / 'data', copied / 'data')
    shutil.copyfile(FOLDER / 'provenance.json', copied / 'provenance.json')
    path = copied / 'data/cf_0ms_real.csv'
    path.write_bytes(path.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='binding changed'):
        load_observations(copied)


def test_changing_sha_alone_does_not_bypass_published_md5(tmp_path):
    copied = tmp_path / 'published'
    shutil.copytree(FOLDER / 'data', copied / 'data')
    provenance = json.loads((FOLDER / 'provenance.json').read_text())
    path = copied / 'data/cf_0ms_imag.csv'
    path.write_bytes(path.read_bytes() + b'\n')
    provenance['files'][path.name]['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    provenance['files'][path.name]['bytes'] = path.stat().st_size
    (copied / 'provenance.json').write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match='binding changed'):
        load_observations(copied)


def test_unexpected_file_inventory_rejected_before_ingestion(tmp_path):
    provenance = json.loads((FOLDER / 'provenance.json').read_text())
    provenance['files']['unpublished.csv'] = {}
    (tmp_path / 'provenance.json').write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match='inventory'):
        load_observations(tmp_path)


def test_off_axis_changes_cannot_enter_axial_calibration(grids):
    real, imag = grids
    original = split_quadrant(real, imag)
    changed = real.copy()
    for row in (8, 12):
        for column in (7, 13):
            changed[row, column] = .9
    altered = split_quadrant(changed, imag)
    assert altered['calibration'] == original['calibration']
    assert altered['validation'] != original['validation']
    assert predict(calibrate_axes(altered['calibration']), altered['validation']) == (
        predict(calibrate_axes(original['calibration']), original['validation']))


def test_product_predictions_preserve_origin_and_match_independent_array_calculation(grids):
    real, imag = grids
    z = real + 1j*imag
    expected = z[11:, 10, None] * z[10, None, 11:]
    residual = z[11:, 11:] - expected
    expected_rmse = float(np.sqrt(np.mean(np.abs(residual)**2)))
    observations = load_observations()
    calibration = calibrate_axes(observations['calibration'])
    predictions = predict(calibration, observations['validation'])
    scored = score(predictions, observations['validation'])
    assert predictions[1, 1]['mean'] == pytest.approx(expected[0, 0], abs=1e-15)
    assert abs(predictions[1, 1]['mean'] - expected[0, 0]/z[10, 10]) > .05
    assert calibration['axis_means'][0, 0] == observations['reported_origin']
    assert calibration['real_statistics'] == 42
    assert scored['withheld_points'] == 100
    assert scored['real_statistics'] == 200
    assert scored['complex_rmse'] == pytest.approx(expected_rmse, abs=1e-15)
    assert expected_rmse == pytest.approx(.10803244966731063, abs=1e-14)
    assert scored['processed_means_establish_sampling_premise'] is False
    assert scored['physical_service_error_budget_included'] is False
    assert scored['empirical_confirmation'] is False
