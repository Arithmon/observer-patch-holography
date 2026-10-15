"""Publication binding, coordinate orientation, and calibration isolation.

Positive-time parser tests use synthetic means. They do not inspect the retained
positive-time experimental outcomes to choose a model or protocol.
"""
import csv
import hashlib
import json
import shutil

import numpy as np
import pytest

from .data import (CHARACTERISTIC_FILE, FOLDER, POSITION_FILE, PUBLISHED_MD5,
                   load_calibration, load_metadata, parse_characteristic,
                   parse_position, verify_provenance)


def _write(path, rows):
    with path.open('w', newline='') as handle:
        csv.writer(handle).writerows(rows)
    return path


def _position_rows():
    rows = [['time (ms)', *list(np.linspace(-.4, .4, 13) / np.sqrt(2))]]
    for index in range(40):
        rows.append([index * .4, *[.01 * index + .001 * probe for probe in range(13)]])
    return rows


def _characteristic_rows():
    rows = [['', '0 ms', '0 ms', '2 ms', '2 ms', '4 ms', '4 ms'],
            ['', 'Real', 'Imag.', 'Real', 'Imag.', 'Real', 'Imag.']]
    for index, coordinate in enumerate(np.linspace(0, 5, 50)):
        rows.append([coordinate, *[.1 * column + .001 * index for column in range(6)]])
    return rows


def _copy_published(tmp_path):
    shutil.copytree(FOLDER / 'data', tmp_path / 'data')
    shutil.copyfile(FOLDER / 'provenance.json', tmp_path / 'provenance.json')
    return tmp_path


def test_all_retained_files_match_independently_published_bytes():
    provenance = verify_provenance()
    assert set(provenance['files']) == set(PUBLISHED_MD5)
    assert provenance['raw_shot_records_included'] is False
    assert provenance['accepted_counts_included'] is False
    assert provenance['license'] == 'CC-BY-4.0'


def test_actual_headers_determine_time_rows_and_explicit_probe_columns():
    metadata = load_metadata()
    position = metadata['position']
    assert len(position['times_ms']) == 40
    assert position['times_ms'][0] == 0
    assert position['times_ms'][-1] == pytest.approx(15.6)
    assert len(position['probe_coordinates']) == 13
    assert position['probe_coordinates'][0] == pytest.approx(-.4 / np.sqrt(2))
    characteristic = metadata['characteristic']
    assert characteristic['times_ms'].tolist() == [0., 2., 4.]
    assert len(characteristic['probe_coordinates']) == 50
    assert characteristic['probe_coordinates'][33] == 3.367346939
    assert metadata['two_dimensional_scans_used'] is False


def test_actual_calibration_contains_only_zero_time_and_no_reconstructed_counts():
    data = load_calibration()
    assert data['position']['times_ms'].tolist() == [0.]
    assert data['characteristic']['times_ms'].tolist() == [0.]
    assert data['position']['spin_imag'].shape == (1, 13)
    assert data['characteristic']['spin_real'].shape == (1, 50)
    assert data['characteristic']['spin_imag'].shape == (1, 50)
    assert set(data['position']) == {'times_ms', 'probe_coordinates', 'spin_imag'}


def test_position_parser_preserves_sign_time_axis_and_every_probe(tmp_path):
    rows = _position_rows()
    rows[1][1] = -.5
    path = _write(tmp_path / 'position.csv', rows)
    calibration = parse_position(path, calibration=True)
    validation = parse_position(path, calibration=False)
    assert calibration['spin_imag'].shape == (1, 13)
    assert validation['spin_imag'].shape == (39, 13)
    assert calibration['spin_imag'][0, 0] == -.5
    assert validation['spin_imag'][0, 12] == pytest.approx(.022)
    assert validation['spin_imag'][-1, 0] == pytest.approx(.39)
    assert validation['times_ms'][0] == .4
    assert calibration['probe_coordinates'][0] < 0


def test_characteristic_parser_does_not_swap_time_real_imag_or_probe_axes(tmp_path):
    path = _write(tmp_path / 'characteristic.csv', _characteristic_rows())
    calibration = parse_characteristic(path, calibration=True)
    validation = parse_characteristic(path, calibration=False)
    assert calibration['spin_real'].shape == calibration['spin_imag'].shape == (1, 50)
    assert validation['spin_real'].shape == validation['spin_imag'].shape == (2, 50)
    assert calibration['spin_real'][0, 7] == pytest.approx(.007)
    assert calibration['spin_imag'][0, 7] == pytest.approx(.107)
    assert validation['spin_real'][0, 7] == pytest.approx(.207)
    assert validation['spin_imag'][1, 49] == pytest.approx(.549)


@pytest.mark.parametrize('kind', ['position', 'characteristic'])
def test_calibration_never_parses_positive_time_values(tmp_path, kind):
    if kind == 'position':
        rows = _position_rows()
        for row in rows[2:]:
            row[1:] = ['unread held-out value'] * 13
        parser = parse_position
    else:
        rows = _characteristic_rows()
        for row in rows[2:]:
            row[3:] = ['unread held-out value'] * 4
        parser = parse_characteristic
    path = _write(tmp_path / 'separated.csv', rows)
    parser(path, calibration=True)
    with pytest.raises(ValueError, match='numeric'):
        parser(path, calibration=False)


@pytest.mark.parametrize('kind', ['position', 'characteristic'])
def test_validation_never_parses_calibration_values(tmp_path, kind):
    if kind == 'position':
        rows = _position_rows()
        rows[1][1:] = ['unread calibration value'] * 13
        parser = parse_position
    else:
        rows = _characteristic_rows()
        for row in rows[2:]:
            row[1:3] = ['unread calibration value'] * 2
        parser = parse_characteristic
    path = _write(tmp_path / 'separated.csv', rows)
    parser(path, calibration=False)
    with pytest.raises(ValueError, match='numeric'):
        parser(path, calibration=True)


@pytest.mark.parametrize('mutation', ['time', 'probe', 'transpose', 'header'])
def test_position_orientation_errors_are_rejected(tmp_path, mutation):
    rows = _position_rows()
    if mutation == 'time':
        rows[1:] = rows[:0:-1]
    elif mutation == 'probe':
        rows[0][1:] = rows[0][:0:-1]
    elif mutation == 'transpose':
        rows = list(zip(*rows))
    else:
        rows[0][0] = 'time (s)'
    with pytest.raises(ValueError):
        parse_position(_write(tmp_path / 'bad.csv', rows), calibration=True)


@pytest.mark.parametrize('mutation', ['time', 'component', 'probe', 'transpose'])
def test_characteristic_orientation_errors_are_rejected(tmp_path, mutation):
    rows = _characteristic_rows()
    if mutation == 'time':
        rows[0][3] = '4 ms'
    elif mutation == 'component':
        rows[1][1], rows[1][2] = rows[1][2], rows[1][1]
    elif mutation == 'probe':
        rows[2:] = rows[:1:-1]
    else:
        rows = list(zip(*rows))
    with pytest.raises(ValueError):
        parse_characteristic(_write(tmp_path / 'bad.csv', rows), calibration=True)


@pytest.mark.parametrize('value', [np.nan, np.inf, -np.inf, 1.0001, -1.0001, 'text'])
@pytest.mark.parametrize('kind', ['position', 'characteristic'])
def test_invalid_selected_spin_means_are_rejected(tmp_path, value, kind):
    rows = _position_rows() if kind == 'position' else _characteristic_rows()
    parser = parse_position if kind == 'position' else parse_characteristic
    rows[1 if kind == 'position' else 2][1] = value
    with pytest.raises(ValueError):
        parser(_write(tmp_path / 'bad.csv', rows), calibration=True)


def test_raw_byte_tampering_fails_even_for_whitespace_only_changes(tmp_path):
    folder = _copy_published(tmp_path)
    path = folder / 'data' / POSITION_FILE
    path.write_bytes(path.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='binding changed'):
        load_calibration(folder)


def test_rehashing_manifest_cannot_override_published_checksum(tmp_path):
    folder = _copy_published(tmp_path)
    path = folder / 'data' / CHARACTERISTIC_FILE
    content = path.read_bytes() + b'\n'
    path.write_bytes(content)
    provenance = json.loads((folder / 'provenance.json').read_text())
    entry = provenance['files'][path.name]
    entry.update(bytes=len(content), sha256=hashlib.sha256(content).hexdigest(),
                 published_md5=hashlib.md5(content).hexdigest())
    (folder / 'provenance.json').write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match='binding changed'):
        verify_provenance(folder)


def test_altering_unscored_two_dimensional_scan_is_detected(tmp_path):
    folder = _copy_published(tmp_path)
    path = folder / 'data/cf_4ms_imag.csv'
    path.write_bytes(path.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='binding changed'):
        verify_provenance(folder)


def test_inventory_is_validated_before_files_are_loaded(tmp_path):
    provenance = json.loads((FOLDER / 'provenance.json').read_text())
    del provenance['files'][POSITION_FILE]
    (tmp_path / 'provenance.json').write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match='inventory'):
        verify_provenance(tmp_path)


@pytest.mark.parametrize('key,value', [('accepted_counts_included', True),
                                      ('raw_shot_records_included', True),
                                      ('doi', '10.0000/unrelated')])
def test_provenance_cannot_claim_unavailable_trial_records(tmp_path, key, value):
    provenance = json.loads((FOLDER / 'provenance.json').read_text())
    provenance[key] = value
    (tmp_path / 'provenance.json').write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match='scope changed'):
        verify_provenance(tmp_path)
