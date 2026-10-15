"""Synthetic CSV controls; no later empirical values are inspected here."""
import csv
import hashlib
import json
import shutil

import pytest

from . import data


def synthetic_csv(path, *, zero='0.5', later='0.25'):
    rows = [['time (ms)', *'abcde']]
    rows.extend([f'{index*.4:.1f}', *([zero if index == 0 else later]*5)]
                 for index in range(26))
    with path.open('w', newline='') as handle:
        csv.writer(handle).writerows(rows)
    return path


def test_metadata_and_calibration_do_not_convert_later_response_strings(tmp_path):
    path = synthetic_csv(tmp_path/'poison.csv', later='DO NOT OPEN LATER RESPONSE')
    assert data.parse_metadata(path)['labels'] == list('abcde')
    assert data.parse_metadata(path)['times_ms'] == [index*2/5 for index in range(26)]
    assert data.parse_zero(path) == dict.fromkeys('abcde', .5)
    with pytest.raises(ValueError):
        data.parse_validation(path)


def test_validation_does_not_read_zero_time_responses(tmp_path):
    path = synthetic_csv(tmp_path/'zero_poison.csv', zero='DO NOT FIT')
    observations = data.parse_validation(path)
    assert len(observations) == 125
    assert len({(row['case'], row['time_ms']) for row in observations}) == 125
    assert {row['case'] for row in observations} == set('abcde')
    assert all(row['time_ms'] > 0 and row['observed_position_estimator'] == .25
               for row in observations)


@pytest.mark.parametrize('target', ['zero', 'later'])
@pytest.mark.parametrize('value', ['nan', 'inf', '-inf'])
def test_only_selected_nonfinite_response_subset_is_rejected(tmp_path, target, value):
    path = synthetic_csv(tmp_path/'nonfinite.csv', **{target: value})
    data.parse_metadata(path)
    with pytest.raises(ValueError):
        (data.parse_zero if target == 'zero' else data.parse_validation)(path)


@pytest.mark.parametrize('mutation', ['header', 'missing_row', 'extra_column', 'off_grid', 'nan_time'])
def test_coordinate_or_shape_changes_are_rejected(tmp_path, mutation):
    path = synthetic_csv(tmp_path/'bad.csv')
    with path.open() as handle:
        rows = list(csv.reader(handle))
    if mutation == 'header':
        rows[0][1] = 'wrong'
    elif mutation == 'missing_row':
        rows.pop()
    elif mutation == 'extra_column':
        rows[2].append('0')
    elif mutation == 'off_grid':
        rows[3][0] = '.81'
    else:
        rows[2][0] = 'nan'
    with path.open('w', newline='') as handle:
        csv.writer(handle).writerows(rows)
    for parser in (data.parse_metadata, data.parse_zero, data.parse_validation):
        with pytest.raises(ValueError):
            parser(path)


@pytest.fixture
def bound_folder(tmp_path):
    # Copy uninterpreted bytes for binding tests, never parse empirical outcomes.
    shutil.copytree(data.FOLDER/'data', tmp_path/'data')
    shutil.copyfile(data.FOLDER/'provenance.json', tmp_path/'provenance.json')
    return tmp_path


def test_frozen_published_bytes_match_independent_identifiers(bound_folder):
    assert set(data.verify_provenance(bound_folder)['files']) == set(data.PUBLISHED_MD5)


@pytest.mark.parametrize('mutation', ['empty', 'omit', 'foreign', 'doi', 'license', 'raw', 'counts', 'rewritten'])
def test_changed_provenance_cannot_reauthorize_different_data(bound_folder, mutation):
    path = bound_folder/'provenance.json'
    provenance = json.loads(path.read_text())
    if mutation == 'empty':
        provenance['files'] = {}
    elif mutation == 'omit':
        del provenance['files'][data.NAME]
    elif mutation == 'foreign':
        provenance['files']['foreign.csv'] = provenance['files'][data.NAME]
    elif mutation in ('doi', 'license'):
        provenance[mutation] = 'changed'
    elif mutation in ('raw', 'counts'):
        provenance['raw_shot_records_included' if mutation == 'raw' else 'accepted_counts_included'] = True
    else:
        raw = b'fabricated replacement'
        (bound_folder/'data'/data.NAME).write_bytes(raw)
        provenance['files'][data.NAME].update(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(),
                                             published_md5=hashlib.md5(raw).hexdigest())
    path.write_text(json.dumps(provenance))
    with pytest.raises(ValueError):
        data.verify_provenance(bound_folder)
