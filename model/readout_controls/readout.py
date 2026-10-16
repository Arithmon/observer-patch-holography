"""Fresh source null instruments and exact arithmetic on deposited means."""
from __future__ import annotations

import csv
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np

from ..anharmonic_observations import data, source
from ..calibrated_measurement.bounds import root_interval
from ..primitive_observer import runtime


def _matrix(value):
    if type(value) is not np.ndarray:
        raise ValueError('an original base ndarray is required before dispatch')
    if (value.ndim != 2 or value.shape[0] < 1 or value.shape[0] != value.shape[1] or value.dtype.kind != 'c'
            or value.dtype.itemsize > 16 or not np.isfinite(value).all()):
        raise ValueError('a finite square binary complex matrix is required')
    return [[(F(z.real), F(z.imag)) for z in row] for row in value]


def _effect(matrix):
    matrix = _matrix(matrix)
    n = len(matrix)
    return [[(sum((matrix[k][i][0]*matrix[k][j][0]+matrix[k][i][1]*matrix[k][j][1]
                   for k in range(n)), F(0)),
              sum((matrix[k][i][0]*matrix[k][j][1]-matrix[k][i][1]*matrix[k][j][0]
                   for k in range(n)), F(0)))
             for j in range(n)] for i in range(n)]


def _combine(a, b, sign=1):
    return [[(x[0]+sign*y[0], x[1]+sign*y[1]) for x, y in zip(ar, br, strict=True)]
            for ar, br in zip(a, b, strict=True)]


def _norm(matrix):
    return root_interval(sum((a*a+b*b for row in matrix for a, b in row), F(0)))[1]


def operator_audit(operations, guarantee):
    if type(operations) is not dict or set(operations) != {'real_zero', 'imag_zero', 'imag_plus', 'imag_minus'}:
        raise ValueError('the four fixed complete detector settings are required')
    if type(guarantee) is not F or not 0 < guarantee < 1:
        raise ValueError('a positive exact guarantee smaller than one is required')
    moments, rows, size = {}, {}, None
    for name, matrices in operations.items():
        if type(matrices) is not tuple or len(matrices) != 2:
            raise ValueError('both pointer outcomes are required')
        effects = [_effect(matrix) for matrix in matrices]
        n = len(effects[0])
        if len(effects[1]) != n or size is not None and size != n:
            raise ValueError('all detector outcomes must share one carrier')
        size = n
        identity = [[(F(i == j), F(0)) for j in range(n)] for i in range(n)]
        moments[name] = _combine(*effects, -1)
        defect = _norm(_combine(_combine(*effects), identity, -1))
        if defect > guarantee:
            raise ArithmeticError('source effect sum exceeds its fixed guarantee')
        rows[name] = {'effect_sum_defect_upper': str(guarantee)}
    defects = {
        'real_origin_minus_identity': _norm(_combine(moments['real_zero'], identity, -1)),
        'imag_origin': _norm(moments['imag_zero']),
        'opposite_imaginary_sum': _norm(_combine(moments['imag_plus'], moments['imag_minus'])),
    }
    if any(value > guarantee for value in defects.values()):
        raise ArithmeticError('source state-independent identity exceeds its guarantee')
    # Complete comparison K_b T^-1/2 is a mathematical instrument. The raw
    # source matrices and runtime states are never normalized or repaired here.
    completion = guarantee*(1+guarantee)/(1-guarantee)
    return {'all_input_effect_bounds': rows,
            'raw_identity_guarantees': {key: str(guarantee) for key in defects},
            'complete_comparison_null_mean_error_upper': str(guarantee+completion),
            'complete_comparison_opposite_sum_error_upper': str(guarantee+2*completion),
            'completion_is_mathematical_not_source_executed': True,
            'original_binary_Kraus_entries_certified': True,
            'inherited_polar_word_representation_error_certified': False,
            'floating_state_contraction_or_hardware_error_certified': False}


def execute(protocol):
    dimension = protocol['source_detector_dimension']
    h = protocol['source_probe_h']
    settings = {'real_zero': (0., 0.), 'imag_zero': (0., -math.pi/2),
                'imag_plus': (h, -math.pi/2), 'imag_minus': (-h, -math.pi/2)}
    operations, metadata = {}, {}
    for name, (v, phase) in settings.items():
        operations[name], metadata[name] = source.detector(dimension, 0., v, phase)
    certificate = operator_audit(operations, F(protocol['source_effect_guarantee']))
    vacuum = np.zeros((dimension, dimension), complex); vacuum[0, 0] = 1
    vector = np.zeros(dimension, complex); vector[0] = 1/math.sqrt(2); vector[-1] = 1j/math.sqrt(2)
    real_vector = np.zeros(dimension, complex); real_vector[:2] = 1/math.sqrt(2)
    states = {'vacuum': vacuum, 'maximally_mixed': np.eye(dimension)/dimension,
              'real_neighbour_superposition': np.outer(real_vector, real_vector.conj()),
              'complex_edge_superposition': np.outer(vector, vector.conj())}
    outputs = {}
    for state_name, density in states.items():
        rows = {}
        for name, matrices in operations.items():
            instrument = runtime.Instrument(kraus=tuple((matrix,) for matrix in matrices), visible=True,
                binding={'source_setting': metadata[name], 'postselection': False})
            program = runtime.Program('source_readout_null', dimension, density,
                (runtime.Event('readout', instruments=(instrument,)),),
                binding={'input_control': state_name, 'evolution_action_is_irrelevant': True})
            history = runtime.execute(program)
            probabilities = [branch.probability for branch in history.branches]
            if len(probabilities) != 2:
                raise ArithmeticError('the complete pointer law lost an outcome')
            rows[name] = {'probabilities': probabilities, 'mean': probabilities[0]-probabilities[1]}
        if (abs(rows['real_zero']['mean']-1) > 1e-10 or abs(rows['imag_zero']['mean']) > 1e-10
                or abs(rows['imag_plus']['mean']+rows['imag_minus']['mean']) > 1e-10):
            raise ArithmeticError('fresh source execution contradicts the null identity')
        outputs[state_name] = rows
    return {'dimension': dimension, 'certificate': certificate,
            'executed_complete_laws': outputs, 'source_metadata': metadata,
            'ideal_origin_law_independent_of_occupation_cutoff': True,
            'physical_preparation_or_action_fit': False}


def _read(path):
    with path.open(newline='') as stream:
        return list(csv.reader(stream))


def _mean(value):
    try:
        number = F(value)
    except (TypeError, ValueError, ZeroDivisionError) as error:
        raise ValueError('a finite printed rational scalar mean is required') from error
    if not -1 <= number <= 1:
        raise ValueError('a printed mean must lie in [-1,1]')
    return number


def compare_data(folder=data.FOLDER):
    """All rows retained, no count inference, correction, normalization or fit."""
    provenance = data.verify_provenance(folder)
    folder = Path(folder)
    position = _read(folder/'data'/data.POSITION_FILE)
    characteristic = _read(folder/'data'/data.CHARACTERISTIC_FILE)
    data._position_metadata(position)
    data._characteristic_metadata(characteristic)
    coordinates = list(map(F, position[0][1:]))
    if len(position) != 41 or len(characteristic) != 52 or coordinates[6] != 0:
        raise ValueError('the complete declared empirical tables are required')
    if F(characteristic[2][0]) != 0:
        raise ValueError('the characteristic origin must be present')
    origins, pairs = [], []
    for row in position[1:]:
        time, means = F(row[0]), list(map(_mean, row[1:]))
        origins.append({'dataset': 'position_scan', 'duration_ms': str(time),
                        'component': 'imag', 'observed_mean': str(means[6]),
                        'predicted_mean': '0', 'residual': str(means[6])})
        for j in range(6):
            k = 12-j
            pairs.append({'duration_ms': str(time), 'negative_coordinate': str(coordinates[j]),
                          'positive_coordinate': str(coordinates[k]),
                          'printed_coordinate_sum': str(coordinates[j]+coordinates[k]),
                          'negative_mean': str(means[j]), 'positive_mean': str(means[k]),
                          'predicted_sum': '0', 'observed_sum': str(means[j]+means[k]),
                          'least_squared_spin_residual_if_exact_opposites': str((means[j]+means[k])**2/2)})
    for index, time in enumerate((0, 2, 4)):
        for component, offset, wanted in (('real', 1, F(1)), ('imag', 2, F(0))):
            mean = _mean(characteristic[2][2*index+offset])
            origins.append({'dataset': 'characteristic_scan', 'duration_ms': str(time),
                            'component': component, 'observed_mean': str(mean),
                            'predicted_mean': str(wanted), 'residual': str(mean-wanted)})
    zero_sse = sum((F(row['residual'])**2 for row in origins), F(0))
    positive_null_sse = sum((F(row['residual'])**2 for row in origins if F(row['duration_ms']) > 0), F(0))
    odd_sse = sum((F(row['least_squared_spin_residual_if_exact_opposites'])
                   for row in pairs if F(row['duration_ms']) > 0), F(0))
    real_deficits = [(1-F(row['observed_mean']))/2 for row in origins if row['component'] == 'real']
    return {'provenance': provenance, 'all_origin_measurements': origins,
            'nominal_opposite_coordinate_pairs': pairs,
            'origin_statistics': {'scalar_means': len(origins), 'squared_residual_sum': str(zero_sse),
                'rmse_diagnostic': math.sqrt(float(zero_sse/len(origins))),
                'real_origin_observed_opposite_frequencies': list(map(str, real_deficits))},
            'all_707_previously_scored_positive_duration_means': {
                'action_independent_origin_squared_residual_floor': str(positive_null_sse),
                'minimum_binary_probability_rmse_from_origins_alone': math.sqrt(float(positive_null_sse/707))/2,
                'conditional_opposite_pair_squared_residual_floor': str(odd_sse),
                'minimum_binary_probability_rmse_with_nominal_antisymmetry': math.sqrt(float((positive_null_sse+odd_sse)/707))/2,
                'floor_is_descriptive_not_a_statistical_rejection': True},
            'maximum_printed_opposite_coordinate_sum': str(max(abs(F(row['printed_coordinate_sum'])) for row in pairs)),
            'statistical_or_hardware_error_inferred': False,
            'missing_metadata': ['accepted trial counts', 'raw outcomes', 'same-state and opposite-coordinate apparatus alignment',
                                 'independently calibrated readout error budget'],
            'action_change_can_remove_ideal_origin_residuals': False}
