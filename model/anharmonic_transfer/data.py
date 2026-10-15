"""Byte-bound data with separate zero-time and later-outcome access paths."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path

FOLDER = Path(__file__).resolve().parent
NAME = 'asymmetric_wells_expectation_x.csv'
LABELS = tuple('abcde')
PUBLISHED_MD5 = {'SOURCE_README.md': '04e5a8d33f1e12c73ae6e9c012d6929e',
                 NAME: 'b3313dde8b217f5dfd397e8ea9608d78'}


def verify_provenance(folder=FOLDER):
    folder = Path(folder)
    provenance = json.loads((folder/'provenance.json').read_text())
    if (set(provenance['files']) != set(PUBLISHED_MD5)
            or provenance['doi'] != '10.5281/zenodo.18752287'
            or provenance['license'] != 'CC-BY-4.0'
            or provenance['raw_shot_records_included'] is not False
            or provenance['accepted_counts_included'] is not False):
        raise ValueError('Published-data inventory or provenance scope changed')
    for name, item in provenance['files'].items():
        raw = (folder/'data'/name).read_bytes()
        if (len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']
                or item['published_md5'] != PUBLISHED_MD5[name]
                or hashlib.md5(raw).hexdigest() != PUBLISHED_MD5[name]):
            raise ValueError(f'Published data bytes changed: {name}')
    return provenance


def _rows(path):
    rows = list(csv.reader(io.StringIO(Path(path).read_text())))
    if len(rows) != 27 or rows[0] != ['time (ms)', *LABELS] or any(len(r) != 6 for r in rows):
        raise ValueError('Expected time and five cases at 26 durations')
    times = [float(r[0]) for r in rows[1:]]
    if any(not math.isfinite(t) or abs(t-.4*i) > 1e-10 for i, t in enumerate(times)):
        raise ValueError('Expected published times 0,.4,...,10 ms')
    return rows, times


def parse_metadata(path):
    """Convert times only; positive-time response strings remain uninterpreted."""
    _, times = _rows(path)
    return {'labels': list(LABELS), 'times_ms': times}


def parse_zero(path):
    rows, _ = _rows(path)
    values = [float(v) for v in rows[1][1:]]
    if not all(math.isfinite(v) for v in values):
        raise ValueError('Finite zero-time observations required')
    return dict(zip(LABELS, values, strict=True))


def parse_validation(path):
    rows, times = _rows(path)
    result = []
    for time, row in zip(times[1:], rows[2:], strict=True):
        for case, value in zip(LABELS, row[1:], strict=True):
            value = float(value)
            if not math.isfinite(value):
                raise ValueError('Finite derived observations required')
            result.append({'case': case, 'time_ms': time, 'observed_position_estimator': value})
    return result


def load_metadata():
    provenance = verify_provenance()
    return {**parse_metadata(FOLDER/'data'/NAME), 'provenance': provenance}


def load_zero():
    verify_provenance()
    return parse_zero(FOLDER/'data'/NAME)


def load_validation():
    verify_provenance()
    return parse_validation(FOLDER/'data'/NAME)
