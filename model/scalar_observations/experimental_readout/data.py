"""Read frozen published spin means without reconstructing trials or successes."""
from pathlib import Path
import hashlib
import json

import numpy as np


FOLDER = Path(__file__).resolve().parent


def load_grid(path):
    grid = np.loadtxt(path, delimiter=',')
    if grid.shape != (21, 21) or not np.isfinite(grid).all() or np.max(np.abs(grid)) > 1:
        raise ValueError('Expected a finite 21 by 21 spin-mean grid in [-1,1]')
    return grid


def orient(real, imag):
    """Identify the conjugated coordinate from exact deposited mirror copies.

    The central line is an actual noisy measurement, not an odd-symmetry copy.
    Its imaginary values must not be forced to zero. This establishes which
    array coordinate is the imaginary displacement, but not its sign or scale.
    """
    for grid in (real, imag):
        if grid.shape != (21, 21) or not np.isfinite(grid).all() or np.max(np.abs(grid)) > 1:
            raise ValueError('Expected a finite 21 by 21 spin-mean grid in [-1,1]')
    z = real+1j*imag
    candidates = []
    for axis in (0, 1):
        negative = np.take(z, range(10), axis=axis)
        positive = np.take(z, range(20, 10, -1), axis=axis)
        if (np.array_equal(negative, positive.conj())
                and np.array_equal(z, np.flip(z, axis=1-axis))):
            candidates.append(axis)
    if len(candidates) != 1:
        raise ValueError('Deposited mirrors do not determine a unique imaginary-coordinate axis')
    return {'imaginary_coordinate_array_axis': candidates[0],
            'real_coordinate_array_axis': 1-candidates[0],
            'coordinate_direction_known': False,
            'coordinate_extent_resolved': False,
            'mirror_copies_are_additional_observations': False,
            'central_line_imaginary_values_preserved': True}


def split_quadrant(real, imag):
    orientation = orient(real, imag)
    z = real+1j*imag
    # Canonical array order is (imaginary-coordinate index, real index).
    if orientation['imaginary_coordinate_array_axis'] == 1:
        z = z.T
    calibration, validation = {}, {}
    for j in range(11):
        for i in range(11):
            destination = calibration if i == 0 or j == 0 else validation
            destination[i, j] = complex(z[10+j, 10+i])
    return {'calibration': calibration, 'validation': validation,
            'orientation': orientation, 'reported_origin': complex(z[10, 10])}


def load_observations(folder=FOLDER):
    folder = Path(folder)
    provenance = json.loads((folder/'provenance.json').read_text())
    expected = {'cf_0ms_real.csv', 'cf_0ms_imag.csv', 'SOURCE_README.md'}
    if set(provenance['files']) != expected:
        raise ValueError('Unexpected published-data file inventory')
    for name, entry in provenance['files'].items():
        content = (folder/'data'/name).read_bytes()
        if (hashlib.sha256(content).hexdigest() != entry['sha256']
                or hashlib.md5(content).hexdigest() != entry['published_md5']
                or len(content) != entry['bytes']):
            raise ValueError(f'Published-data binding changed: {name}')
    result = split_quadrant(load_grid(folder/'data/cf_0ms_real.csv'),
                            load_grid(folder/'data/cf_0ms_imag.csv'))
    result['provenance'] = provenance
    return result
