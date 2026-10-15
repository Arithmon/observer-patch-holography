"""Native load covariance, recorded increments and a controlled Gaussian limit.

Replication is an explicit product-history premise. Neither a physical clock
nor a spatial continuum is inferred from a growing number of copies.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
import math

import numpy as np

from ..common_source.native_fields import mean_operator
from ..observer import native
from ..primitive_observer import runtime
from ..primitive_observer.transfer import native_counts, exact_audit


VALUES = (-1, 0, 1)
WEIGHTS = (1, -2, 1)
TOLERANCE = 2e-10


def matmul(left, right):
    return tuple(tuple(sum((a*b for a, b in zip(row, col, strict=True)), F(0))
                       for col in zip(*right, strict=True)) for row in left)


def transpose(matrix):
    return tuple(zip(*matrix, strict=True))


def bounded(value, name):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= TOLERANCE:
        raise ArithmeticError(f'{name} exceeds the numerical comparison tolerance')
    return value


def source_moments(kernel):
    """Check the closure and covariance on the actual finite transition law."""
    exact_audit(kernel)
    centered = tuple(tuple(F(x)-F(1, 6) for x in state) for state in kernel.states)
    a = mean_operator()
    c = tuple(tuple(sum((x[i]*x[j] for x in centered), F(0))/66
                    for j in range(12)) for i in range(12))
    expected = tuple(tuple(F(5, 33)*(F(i == j)-F(1, 12)) for j in range(12))
                     for i in range(12))
    if c != expected:
        raise ArithmeticError('the stationary native load covariance changed')
    conditional = matmul(tuple(tuple(F(v, 60) for v in row) for row in kernel.counts), centered)
    if conditional != matmul(centered, transpose(a)):
        raise ArithmeticError('the native conditional mean is not the claimed closed response')
    aca = matmul(matmul(a, c), transpose(a))
    d = tuple(tuple(c[i][j]-aca[i][j] for j in range(12)) for i in range(12))
    # Compare several lag covariances with direct finite-state propagation.
    propagated = centered
    covariance = c
    for _ in range(5):
        direct = tuple(tuple(sum((centered[k][i]*propagated[k][j] for k in range(66)), F(0))/66
                             for j in range(12)) for i in range(12))
        if direct != covariance:
            raise ArithmeticError('finite-state and closed lag covariance disagree')
        propagated = matmul(tuple(tuple(F(v, 60) for v in row) for row in kernel.counts), propagated)
        covariance = matmul(covariance, transpose(a))
    return a, c, d


def marked_counts(kernel, port=0):
    if type(port) is not int or not 0 <= port < 12:
        raise ValueError('port must be an integer from zero through eleven')
    return tuple(tuple(tuple(value if kernel.states[j][port]-kernel.states[i][port] == delta else 0
                             for j, value in enumerate(row)) for i, row in enumerate(kernel.counts))
                 for delta in VALUES)


def word_law(kernel, steps=3, port=0):
    """Full signed-current record law, with all terminal native populations."""
    if type(steps) is not int or not 1 <= steps <= 4:
        raise ValueError('one through four retained native attempts are supported')
    marked = marked_counts(kernel, port)
    branches = {(): (F(1, 66),)*66}
    for _ in range(steps):
        following = {}
        for word, distribution in branches.items():
            for delta, matrix in zip(VALUES, marked, strict=True):
                row = [F(0)]*66
                for i, mass in enumerate(distribution):
                    for j, count in enumerate(matrix[i]):
                        if mass and count:
                            row[j] += mass*F(count, 60)
                following[word+(delta,)] = tuple(row)
        branches = following
    if sum(sum(row) for row in branches.values()) != 1:
        raise ArithmeticError('complete signed-current law lost probability mass')
    return branches


def current_instrument(kernel, port=0):
    """Coarsen actual native C records to the signed load change at one port."""
    marked_counts(kernel, port)  # Validate the public port domain first.
    index = {state: i for i, state in enumerate(kernel.states)}
    maps = {}
    deltas = {}
    for e, seam in enumerate(kernel.seams):
        for coin in (0, 1):
            for i, state in enumerate(kernel.states):
                output, record = native.native(state, seam, coin)
                key = (e, *record)
                delta = -(int(port == seam[0])-int(port == seam[1]))*record[1]
                if output[port]-state[port] != delta:
                    raise ArithmeticError('native record fails the pathwise continuity identity')
                if key in deltas and deltas[key] != delta:
                    raise ArithmeticError('signed change is not determined by the admitted record')
                deltas[key] = delta
                maps.setdefault(key, np.zeros((66, 66), complex))[index[output], i] = 1/math.sqrt(60)
    groups = tuple(tuple(maps[key] for key in sorted(maps) if deltas[key] == delta) for delta in VALUES)
    return runtime.Instrument(kraus=groups, visible=True,
        resources=(('native_attempts', 1), ('retained_signed_current_records', 1)),
        binding={'native_law': 'C', 'record_values_in_outcome_order': list(VALUES),
                 'port': port, 'record': 'delta x_port = -(1_port=a - 1_port=b)*j',
                 'scheduler': 'uniform thirty seams and fair coin',
                 'all_other_native_labels': 'probability summed',
                 'native_channel_semantics': True, 'native_pulse_words_emitted': False})


def current_program(kernel, steps=3, port=0):
    if type(steps) is not int or not 1 <= steps <= 4:
        raise ValueError('one through four retained native attempts are supported')
    instrument = current_instrument(kernel, port)
    return runtime.Program('native_signed_current_history', 66, np.eye(66)/66,
        tuple(runtime.Event(f'current_{i}', instruments=(instrument,)) for i in range(steps)),
        binding={'preparation': 'stationary mixture I66/66',
                 'initial_loads_are_not_observed': True, 'port': port,
                 'parameter': 'integer native attempts'})


def execute_current_check(kernel):
    program = current_program(kernel)
    history = runtime.execute(program)
    law = word_law(kernel)
    error = 0.
    rows = []
    for branch in history.branches:
        word = tuple(VALUES[i] for i in branch.record)
        terminal = law[word]
        error += float(np.linalg.norm(branch.density-np.diag(list(map(float, terminal))), 'fro'))
        if dict(branch.resources) != {'native_attempts': 3, 'retained_signed_current_records': 3}:
            raise ArithmeticError('native current resource accounting changed')
        rows.append({'signed_current_word': list(word), 'probability': str(sum(terminal)),
                     'executed_probability': branch.probability})
    if len(rows) != 27:
        raise ArithmeticError('the complete three-record dictionary changed')
    prefix = runtime.execute(program, stop_after=1)
    restored = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(program, checkpoint=restored)
    return {'records': rows, 'terminal_dimension_per_record': 66,
            'complete_terminal_density_residual': bounded(error, 'signed current terminal output'),
            'checkpoint_distance': bounded(runtime.complete_distance(history, continued), 'checkpoint')}


def score_distribution(law, weights=WEIGHTS):
    weights = tuple(weights)
    if not weights or any(type(w) is not int for w in weights) or any(len(word) != len(weights) for word in law):
        raise ValueError('one integer score weight per retained current is required')
    result = defaultdict(F)
    for word, terminal in law.items():
        result[sum(w*x for w, x in zip(weights, word, strict=True))] += sum(terminal)
    return {score: probability for score, probability in sorted(result.items()) if probability}


def score_moments(distribution):
    if (not distribution or any(type(k) is not int or type(p) not in (int, F) or p < 0
                                for k, p in distribution.items()) or sum(distribution.values()) != 1):
        raise ValueError('a complete exact normalized integer-score law is required')
    mean = sum((F(k)*p for k, p in distribution.items()), F(0))
    if mean != 0:
        raise ValueError('the Gaussian fluctuation comparison requires a centered score')
    variance = sum((k*k*p for k, p in distribution.items()), F(0))
    third_absolute = sum((abs(k)**3*p for k, p in distribution.items()), F(0))
    fourth = sum((k**4*p for k, p in distribution.items()), F(0))
    return variance, third_absolute, fourth-3*variance**2


def characteristic_comparison(distribution, replicas, theta=1.):
    if type(replicas) is not int or not 1 <= replicas <= 1000000:
        raise ValueError('one through one million independent replicas are supported')
    if isinstance(theta, bool) or not math.isfinite(float(theta)) or abs(float(theta)) > 10:
        raise ValueError('a finite Fourier coordinate with magnitude at most ten is required')
    theta = float(theta)
    variance, third, fourth_cumulant = score_moments(distribution)
    single = sum(float(p)*np.exp(1j*theta*k/math.sqrt(replicas)) for k, p in distribution.items())
    exact_product = single**replicas
    gaussian = math.exp(-float(variance)*theta*theta/2)
    # Taylor remainder + telescoping of M factors, each of modulus <= 1.
    analytic_bound = min(2., float(third)*abs(theta)**3/(6*math.sqrt(replicas))
                         +float(variance)**2*theta**4/(8*replicas))
    error = abs(exact_product-gaussian)
    if error > analytic_bound+TOLERANCE:
        raise ArithmeticError('the finite product violates its characteristic-function bound')
    return {'independent_replicas': replicas, 'Fourier_coordinate': theta,
            'native_product_characteristic_real': float(exact_product.real),
            'native_product_characteristic_imag': float(exact_product.imag),
            'Gaussian_limit_characteristic': gaussian, 'difference_diagnostic': float(error),
            'analytic_characteristic_error_bound': analytic_bound,
            'variance_exact': str(variance),
            'fourth_cumulant_exact': str(fourth_cumulant/replicas),
            'perfectly_correlated_replica_variance_control': str(replicas*variance)}


def gaussian_action_check(a, c, d, steps=3):
    if type(steps) is not int or not 1 <= steps <= 8:
        raise ValueError('one through eight history intervals are supported')
    # Orthonormal coordinates only remove the conserved total-load direction.
    basis = np.linalg.qr(np.eye(12)[:, :11]-np.eye(12)[:, [11]])[0]
    a_r = basis.T@np.asarray(a, float)@basis
    c_r = basis.T@np.asarray(c, float)@basis
    d_r = basis.T@np.asarray(d, float)@basis
    if np.linalg.eigvalsh(d_r)[0] <= 0 or not 0 < np.linalg.eigvalsh(a_r)[0] < 1:
        raise ArithmeticError('the native mean-free Gaussian action is not positive')
    size = (steps+1)*11
    transform = np.eye(size)
    precision_of_innovations = np.zeros((size, size))
    precision_of_innovations[:11, :11] = np.linalg.inv(c_r)
    for i in range(1, steps+1):
        transform[11*i:11*(i+1), 11*(i-1):11*i] = -a_r
        precision_of_innovations[11*i:11*(i+1), 11*i:11*(i+1)] = np.linalg.inv(d_r)
    precision = transform.T@precision_of_innovations@transform
    covariance = np.block([[c_r@np.linalg.matrix_power(a_r, abs(i-j))
                            for j in range(steps+1)] for i in range(steps+1)])
    residual = np.linalg.norm(precision@covariance-np.eye(size), 'fro')
    # Exact lag law predicts recorded increments without granting load readout.
    return {'mean_free_dimensions': 11, 'history_intervals': steps,
            'history_scalar_coordinates': size,
            'innovation_covariance_min_eigenvalue': float(np.linalg.eigvalsh(d_r)[0]),
            'mean_multiplier_spectrum': list(map(float, np.linalg.eigvalsh(a_r))),
            'precision_times_native_covariance_residual': bounded(residual, 'Gaussian action inverse'),
            'stationary_initial_term_included': True,
            'parameter': 'native attempt spacing one',
            'Gaussian_noise_supplied_at_single_carrier': False}
