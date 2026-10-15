"""Read immutable published CF means, keeping calibration separate from validation.

CSV coordinates and imaginary signs are preserved verbatim. Their identification
with a physical displacement is a protocol premise, not an ingestion fit. The
downloaded two-dimensional scans are provenance-bound but are not scored here.
No accepted shot counts or independent shot records are inferred from means.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


FOLDER = Path(__file__).resolve().parent
POSITION_FILE = 'symmetric_well_characteristic_functions_for_x.csv'
CHARACTERISTIC_FILE = 'symmetric_well_characteristic_functions_for_P.csv'
PUBLISHED_MD5 = {
    'SOURCE_README.md': '04e5a8d33f1e12c73ae6e9c012d6929e',
    POSITION_FILE: '9bb5fcec7441c62517cbb264d9ff89d0',
    CHARACTERISTIC_FILE: 'f4f3da41d68cf631d1c323bfeec7c466',
    'cf_2ms_real.csv': 'a199954171691f29ad0162334d249c93',
    'cf_2ms_imag.csv': 'bf381902bcdbc59489cee63b79e40f4a',
    'cf_4ms_real.csv': '20228fafdee578898a27635ee6a0bf98',
    'cf_4ms_imag.csv': '532151858f13958fb119a96748104594',
}


def verify_provenance(folder=FOLDER):
    """Check byte counts and both local SHA-256 and independently pinned MD5."""
    folder = Path(folder)
    provenance = json.loads((folder / 'provenance.json').read_text())
    if set(provenance['files']) != set(PUBLISHED_MD5):
        raise ValueError('Unexpected published-data file inventory')
    if (provenance['doi'] != '10.5281/zenodo.18752287'
            or provenance['license'] != 'CC-BY-4.0'
            or provenance['raw_shot_records_included'] is not False
            or provenance['accepted_counts_included'] is not False):
        raise ValueError('Published-data provenance scope changed')
    for name, expected_md5 in PUBLISHED_MD5.items():
        entry = provenance['files'][name]
        content = (folder / 'data' / name).read_bytes()
        if (entry['published_md5'] != expected_md5
                or hashlib.md5(content).hexdigest() != expected_md5
                or hashlib.sha256(content).hexdigest() != entry['sha256']
                or len(content) != entry['bytes']):
            raise ValueError(f'Published-data binding changed: {name}')
    return provenance


def _rows(path, shape):
    with Path(path).open(newline='') as handle:
        rows = list(csv.reader(handle))
    if len(rows) != shape[0] or any(len(row) != shape[1] for row in rows):
        raise ValueError(f'Expected CSV shape {shape}, including headers')
    return rows


def _numbers(values, label):
    try:
        result = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f'Expected numeric {label}') from error
    if not np.isfinite(result).all():
        raise ValueError(f'Expected finite {label}')
    return result


def _means(values):
    result = _numbers(values, 'spin means')
    if np.any(np.abs(result) > 1):
        raise ValueError('Expected spin means in [-1,1]')
    return result


def _position_metadata(rows):
    if rows[0][0] != 'time (ms)':
        raise ValueError('Expected time (ms) position header')
    coordinates = _numbers(rows[0][1:], 'position probe coordinates')
    times = _numbers([row[0] for row in rows[1:]], 'position times')
    if not np.allclose(coordinates, np.linspace(-.4, .4, 13) / np.sqrt(2),
                       rtol=0, atol=2e-15):
        raise ValueError('Unexpected position probe coordinates or orientation')
    if not np.allclose(times, np.arange(40) * .4, rtol=0, atol=5e-15):
        raise ValueError('Unexpected position time grid or orientation')
    return times, coordinates


def _characteristic_metadata(rows):
    if rows[0] != ['', '0 ms', '0 ms', '2 ms', '2 ms', '4 ms', '4 ms']:
        raise ValueError('Unexpected characteristic time headers or orientation')
    if rows[1] != ['', 'Real', 'Imag.', 'Real', 'Imag.', 'Real', 'Imag.']:
        raise ValueError('Unexpected characteristic component headers or orientation')
    coordinates = _numbers([row[0] for row in rows[2:]], 'characteristic probes')
    # One deposited coordinate is printed to nine decimal places. Preserve it;
    # do not silently replace the column with an ideal linspace.
    if not np.allclose(coordinates, np.linspace(0, 5, 50), rtol=0, atol=3e-10):
        raise ValueError('Unexpected characteristic probe coordinates or orientation')
    return np.array([0., 2., 4.]), coordinates


def load_metadata(folder=FOLDER):
    """Expose only coordinate/time headers, without converting any spin means."""
    folder = Path(folder)
    provenance = verify_provenance(folder)
    position = _rows(folder / 'data' / POSITION_FILE, (41, 14))
    characteristic = _rows(folder / 'data' / CHARACTERISTIC_FILE, (52, 7))
    position_times, position_coordinates = _position_metadata(position)
    characteristic_times, characteristic_coordinates = _characteristic_metadata(characteristic)
    return {
        'position': {'times_ms': position_times, 'probe_coordinates': position_coordinates},
        'characteristic': {'times_ms': characteristic_times,
                           'probe_coordinates': characteristic_coordinates},
        'coordinates_and_imaginary_signs_preserved': True,
        'two_dimensional_scans_used': False,
        'provenance': provenance,
    }


def parse_position(path, *, calibration):
    """Parse only t=0 or only t>0 means; keep each time's probes in row order."""
    rows = _rows(path, (41, 14))
    times, coordinates = _position_metadata(rows)
    indices = [0] if calibration else list(range(1, len(times)))
    return {'times_ms': times[indices], 'probe_coordinates': coordinates,
            'spin_imag': _means([rows[index + 1][1:] for index in indices])}


def parse_characteristic(path, *, calibration):
    """Parse only requested time columns, explicitly alternating real and imag."""
    rows = _rows(path, (52, 7))
    times, coordinates = _characteristic_metadata(rows)
    indices = [0] if calibration else [1, 2]
    return {
        'times_ms': times[indices], 'probe_coordinates': coordinates,
        'spin_real': _means([[row[1 + 2 * index] for row in rows[2:]]
                             for index in indices]),
        'spin_imag': _means([[row[2 + 2 * index] for row in rows[2:]]
                             for index in indices]),
    }


def _load(folder, *, calibration):
    folder = Path(folder)
    provenance = verify_provenance(folder)
    return {
        'position': parse_position(folder / 'data' / POSITION_FILE,
                                   calibration=calibration),
        'characteristic': parse_characteristic(folder / 'data' / CHARACTERISTIC_FILE,
                                               calibration=calibration),
        'provenance': provenance,
    }


def load_calibration(folder=FOLDER):
    """Load the t=0 subset without parsing positive-time observation values."""
    return _load(folder, calibration=True)


def load_validation(folder=FOLDER):
    """Load strictly positive times, for scoring after the protocol is fixed."""
    return _load(folder, calibration=False)
