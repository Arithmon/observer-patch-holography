"""Recover a finite detector Hilbert space from incident-clock correlations.

The recovery function consumes rational correlations only. Native source
enumeration is a separate certificate of those inputs and their continuation.
Exact rank is a law-level result; sampled records have a separate error budget.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path
import random
from typing import NamedTuple

import sympy as sp

from ..observer import native
from ..observer_time_fields import core
from ..native_action import incident_refinement as clocks
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS, *clocks.screen.SOURCE_PATHS,
    'model/native_action/incident_refinement.py',
    'model/relative_observer_geometry/temporal_records.py'))))


def _rational(value):
    return F(int(value.p), int(value.q))


def _matrix_rows(matrix):
    return tuple(tuple(_rational(value) for value in row) for row in matrix.tolist())


def _positive(matrix):
    # Sylvester's criterion is exact here; no numerical rank tolerance occurs.
    return all(matrix[:i, :i].det(method='domain-ge') > 0
               for i in range(1, matrix.rows+1))


class Representation(NamedTuple):
    moments: tuple
    dimension: int
    gram: tuple
    shifted_gram: tuple
    transfer: tuple
    polynomial: tuple  # ascending powers, leading coefficient one
    gram_lower: F


def recover(moments):
    """Recover the cyclic positive transfer from an exact flat moment sequence.

    The declared model class is c_k=<f,T**k f>, with 0<T<I on this cyclic
    space. A zero Schur remainder at lag 2r proves the degree-r recurrence
    under that class. Finite noisy moments do not establish this equality.
    """
    if (type(moments) is not tuple or not 3 <= len(moments) <= 31
            or len(moments) % 2 != 1
            or any(type(value) is not F for value in moments)):
        raise ValueError('an odd bounded tuple of exact rational correlations is required')
    values = [sp.Rational(value.numerator, value.denominator) for value in moments]
    if values[0] <= 0:
        raise ValueError('a nonzero positive field variance is required')
    for rank in range(1, (len(values)+1)//2):
        gram = sp.Matrix(rank, rank, lambda i, j: values[i+j])
        if not _positive(gram):
            raise ValueError('the leading temporal Gram is not positive definite')
        next_column = sp.Matrix(values[rank:2*rank])
        coefficients = gram.inv()*next_column
        remainder = values[2*rank]-(next_column.T*coefficients)[0]
        if remainder < 0:
            raise ValueError('the temporal moments have a negative Schur remainder')
        if remainder:
            continue
        shifted = sp.Matrix(rank, rank, lambda i, j: values[i+j+1])
        if not _positive(shifted) or not _positive(gram-shifted):
            raise ValueError('the recovered transfer is not a strict positive contraction')
        transfer = sp.zeros(rank)
        for i in range(rank-1):
            transfer[i+1, i] = 1
        transfer[:, rank-1] = coefficients
        if gram*transfer != shifted:
            raise ArithmeticError('the record metric does not make transfer self-adjoint')
        for index in range(rank, len(values)):
            if values[index] != sum(coefficients[j]*values[index-rank+j] for j in range(rank)):
                raise ValueError('later supplied correlations violate the flat recurrence')
        return Representation(tuple(F(v.numerator, v.denominator) for v in moments), rank,
            _matrix_rows(gram), _matrix_rows(shifted), _matrix_rows(transfer),
            tuple(_rational(-v) for v in coefficients)+(F(1),),
            _rational(1/sp.trace(gram.inv())))
    raise ValueError('the supplied correlations do not certify a finite flat extension')


def predict(representation, lag):
    """Predict a correlation from a previously checked record representation."""
    if type(representation) is not Representation or representation != recover(representation.moments):
        raise ValueError('an authentic recovered temporal representation is required')
    if type(lag) is not int or not 0 <= lag <= 128:
        raise ValueError('an exact correlation lag from zero through 128 is required')
    values = list(representation.moments)
    r = representation.dimension
    while len(values) <= lag:
        values.append(-sum(representation.polynomial[j]*values[-r+j] for j in range(r)))
    return values[lag]


@lru_cache(maxsize=1)
def _mean_kernel():
    source = core.source()
    identity, remote, incident = sp.eye(12), sp.zeros(12), sp.zeros(12)
    for a, b in source.seams:
        swap = sp.eye(12)
        swap[a, a] = swap[b, b] = 0
        swap[a, b] = swap[b, a] = 1
        if {a, b}.intersection(source.region):
            incident += (identity+swap)/source.denominator
        else:
            remote += (identity+swap)/source.denominator
    return incident*(identity-remote).inv()


def source_correlations(max_lag=20):
    """All four two-field joint laws, via the independent full 66-state kernel.

    Only initial local bits and incident signed currents are observer inputs.
    Full-state propagation certifies their exact laws; it is not observer data.
    """
    if type(max_lag) is not int or not 14 <= max_lag <= 32:
        raise ValueError('an exact source lag from fourteen through 32 is required')
    source, kernel = core.source(), core.completed()
    indicators = tuple(tuple(state[port] for state in source.states) for port in (0, 1))
    columns = indicators
    denominator = len(source.states)
    rows = []
    mean_kernel, mean_power = _mean_kernel(), sp.eye(12)
    for lag in range(max_lag+1):
        joint = tuple(tuple(F(sum(a*b for a, b in zip(indicators[i], columns[j], strict=True)),
                                  denominator) for j in (0, 1)) for i in (0, 1))
        for i in (0, 1):
            for j in (0, 1):
                independent = sp.Rational(1, 36)+sp.Rational(5, 33)*(mean_power[i, j]-sp.Rational(1, 12))
                if joint[i][j] != _rational(independent):
                    raise ArithmeticError('full native records and the occupation reduction disagree')
        rows.append(joint)
        if lag != max_lag:
            columns = tuple(tuple(sum(a*b for a, b in zip(row, column, strict=True))
                                  for row in kernel.counts) for column in columns)
            denominator *= kernel.denominator
            mean_power = mean_kernel*mean_power
    return tuple(rows)


def _outward(value, upper=False):
    scale = 10**12
    units = -((-value.numerator*scale)//value.denominator) if upper else value.numerator*scale//value.denominator
    return F(units, scale)


def _negative_log_bounds(value, terms=24):
    # -log(x)=2 sum z^(2j+1)/(2j+1), z=(1-x)/(1+x), 0<x<1.
    z = (1-value)/(1+value)
    lower = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    return lower, lower+2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))


def spectrum(representation):
    if type(representation) is not Representation or representation != recover(representation.moments):
        raise ValueError('an authentic recovered temporal representation is required')
    x = sp.Symbol('x')
    polynomial = sp.Poly(sum(sp.Rational(c.numerator, c.denominator)*x**i
                             for i, c in enumerate(representation.polynomial)), x)
    intervals = sp.polys.polytools.intervals(polynomial, eps=sp.Rational(1, 10**12))
    result = []
    for (left, right), multiplicity in intervals:
        left, right = _rational(left), _rational(right)
        if not 0 < left <= right < 1 or multiplicity != 1:
            raise ArithmeticError('the detector polynomial lost its distinct positive roots')
        energy = _negative_log_bounds(right)[0], _negative_log_bounds(left)[1]
        result.append({'transfer_interval': [str(_outward(left)), str(_outward(right, True))],
            'energy_interval_per_incident_tick': [str(_outward(energy[0])), str(_outward(energy[1], True))]})
    return result


def sampling_bound(samples, tolerance, settings=4):
    """Two-sided Hoeffding union bound, with rational exponential majorant."""
    if (type(samples) is not int or not 1 <= samples <= 10**9
            or type(settings) is not int or not 1 <= settings <= 64
            or type(tolerance) is not F or not 0 < tolerance < 1):
        raise ValueError('bounded sample count, settings and rational tolerance required')
    exponent = 2*samples*tolerance*tolerance
    denominator = sum((exponent**j/math.factorial(j) for j in range(65)), F(0))
    return min(F(1), 2*settings/denominator)


def perturbation_budget(representation, error):
    """Uniform correlation error gives a certified generalized-eigenvalue bound.

    If every used correlation changes by <=error, both Hankel errors have
    operator norm <=r*error. With g<=lambda_min(H0), the ordered generalized
    eigenvalues change by <=2*r*error/(g-r*error), when r*error<g.
    This does not certify exact rank or the infinite continuation from samples.
    """
    if type(representation) is not Representation or representation != recover(representation.moments):
        raise ValueError('an authentic recovered temporal representation is required')
    if type(error) is not F or error < 0:
        raise ValueError('a nonnegative exact correlation error is required')
    eta = representation.dimension*error
    if eta >= representation.gram_lower:
        return None
    return 2*eta/(representation.gram_lower-eta)


def sampled_records(seed=24017, samples=2048):
    """Independent reset histories; later fields use incident currents only.

    The seed reproduces executions, not an IID-entropy certificate. The bound
    is conditional on the stated independent uniform preparation/scheduler and
    fair native coin services. There is no postselection or capped timeout.
    """
    if type(seed) is not int or not 0 <= seed < 2**32 or type(samples) is not int or not 1 <= samples <= 8192:
        raise ValueError('a bounded exact seed and sample count required')
    source = core.source()
    rng = random.Random(seed)
    rows = []
    for lag in (1, 2, 3, 4):
        successes, attempts, digest = 0, 0, hashlib.sha256()
        for _ in range(samples):
            state = source.states[rng.randrange(66)]
            initial, counter, ticks = state[0], state[0], 0
            digest.update(bytes((initial, state[1])))
            while ticks < lag:
                seam = source.seams[rng.randrange(30)]
                state, (_, current) = native.native(state, seam, rng.randrange(2))
                attempts += 1
                if set(seam).intersection(source.region):
                    ticks += 1
                    counter += current*(int(seam[1] == 0)-int(seam[0] == 0))
                    if counter not in (0, 1):
                        raise ArithmeticError('authentic incident currents lost the local field')
                    digest.update(bytes((*seam, current+1)))
            # This hidden-state equality is an audit, never an observer input.
            if counter != state[0]:
                raise ArithmeticError('record counter disagrees with the native occupation')
            successes += initial*counter
        rows.append({'lag_in_incident_ticks': lag, 'samples': samples,
            'joint_ones_count': successes, 'native_attempts': attempts,
            'accessible_record_sha256': digest.hexdigest()})
    return rows


def screen_compression():
    """Certified finite temporal subspace of both actual screen processes.

    This uses their Bernoulli density and whole-root clock, not the isolated
    total-two source's spectrum. It is a compression, not a closed model.
    """
    ports, edges, incident, observed, exits, _ = clocks.geometry()
    counts = clocks.exit_counts(ports, edges, observed, exits, 48)
    intervals = [(F(2, 9), F(2, 9))]
    for row in clocks.density_response(ports, edges, incident, 3):
        low, high = map(F, row['return_probability_interval'])
        spatial, _ = clocks.exit_bound(counts, row['scheduled_ticks'], len(incident))
        intervals.append((F(2, 9)*low-spatial, F(2, 9)*high+spatial))
    # The certificate uses the exported intervals themselves, including their
    # additional outward rounding, so consumers can reconstruct its budget.
    intervals = [(_outward(lo), _outward(hi, True)) for lo, hi in intervals]
    centers = [(lo+hi)/2 for lo, hi in intervals]
    error = max((hi-lo)/2 for lo, hi in intervals)
    gram = sp.Matrix(2, 2, lambda i, j: sp.Rational(centers[i+j]))
    shifted = sp.Matrix(2, 2, lambda i, j: sp.Rational(centers[i+j+1]))
    lowers = []
    for matrix, allowance in ((gram, 2*error), (shifted, 2*error), (gram-shifted, 4*error)):
        if not _positive(matrix):
            raise ArithmeticError('the source temporal pencil has no positive midpoint certificate')
        lower = _rational(1/sp.trace(matrix.inv()))-allowance
        if lower <= 0:
            raise ArithmeticError('the full-screen temporal positivity is unresolved by its error budget')
        lowers.append(lower)
    # Min-max perturbation of the generalized Rayleigh quotient. Since both
    # pencils are positive contractions, numerator and denominator errors add.
    spectral_error = 4*error/lowers[0]
    x = sp.Symbol('x')
    polynomial = sp.Poly((shifted-x*gram).det(), x).monic()
    roots = sp.polys.polytools.intervals(polynomial, eps=sp.Rational(1, 10**12))
    spectral = []
    for (left, right), multiplicity in roots:
        if multiplicity != 1:
            raise ArithmeticError('the two temporal Ritz values merged')
        spectral.append([str(_outward(_rational(left)-spectral_error)),
                         str(_outward(_rational(right)+spectral_error, True))])
    return {'source_levels': [3, 5], 'preparation': 'Bernoulli product density 1/3',
        'clock_region_ports': 12, 'incident_seams': len(incident),
        'centered_root_field_covariance_intervals': [[str(_outward(lo)), str(_outward(hi, True))]
                                                    for lo, hi in intervals],
        'temporal_vectors': ['f', 'T f'], 'temporal_gram_dimension_lower': 2,
        'gram_positive_lower': str(lowers[0]), 'shifted_gram_positive_lower': str(lowers[1]),
        'gram_minus_shifted_positive_lower': str(lowers[2]),
        'uniform_covariance_midpoint_error': str(error),
        'ordered_compressed_transfer_eigenvalue_intervals': spectral,
        'generalized_eigenvalue_error_upper': str(spectral_error),
        'intervals_hold_separately_for_either_actual_full_screen': True,
        'compressed_transfer_is_closed_or_equals_isolated_carrier_spectrum': False,
        'all_level_continuum_or_physical_time_claimed': False}


def run():
    joint = source_correlations(20)
    centered = tuple(row[0][0]-F(1, 36) for row in joint)
    scalar = recover(centered[:15])
    # Both fields are already readable from the same incident record service.
    plus = tuple(sum(row[i][j]-F(1, 36) for i in (0, 1) for j in (0, 1)) for row in joint)
    minus = tuple(sum((1 if i == j else -1)*(row[i][j]-F(1, 36))
                       for i in (0, 1) for j in (0, 1)) for row in joint)
    symmetric, antisymmetric = recover(plus[:9]), recover(minus[:7])
    if (scalar.dimension, symmetric.dimension, antisymmetric.dimension) != (7, 4, 3):
        raise ArithmeticError('the exact detector-visible mode census changed')
    if any(predict(scalar, k) != centered[k] for k in range(15, 21)):
        raise ArithmeticError('the detector recurrence failed held-out native lags')
    if any(row[0][0] != row[1][1] or row[0][1] != row[1][0] for row in joint):
        raise ArithmeticError('the root exchange symmetry changed')
    # Exact spatial observability rank, independent of the correlation pencil.
    kernel = _mean_kernel()
    q = sp.eye(12)-sp.ones(12)/12
    rows = []
    for port in (0, 1):
        vector = q[:, port].T
        for _ in range(12):
            rows.append(vector)
            vector = vector*kernel
    observability = sp.Matrix.vstack(*rows)
    visible_rank = observability.rank()
    if visible_rank != 7:
        raise ArithmeticError('the calibrated occupation-space visibility changed')
    dark = sp.Matrix.hstack(*sp.Matrix.vstack(observability, sp.ones(1, 12)).nullspace())
    dark_transfer, parameters = dark.gauss_jordan_solve(kernel*dark)
    if dark.cols != 4 or parameters.rows or observability*dark != sp.zeros(24, 4):
        raise ArithmeticError('the centered dark occupation space is not invariant')
    sampled = sampled_records()
    tolerance = F(1, 16)
    for row in sampled:
        predicted = joint[row['lag_in_incident_ticks']][0][0]
        discrepancy = abs(F(row['joint_ones_count'], row['samples'])-predicted)
        row.update(exact_joint_probability=str(predicted), absolute_error=str(discrepancy),
                   within_declared_error=discrepancy <= tolerance)

    def representation(value):
        return {'dimension': value.dimension, 'recurrence_polynomial_ascending': list(map(str, value.polynomial)),
            'gram_minimum_eigenvalue_lower': str(value.gram_lower),
            'transfer_spectrum': spectrum(value)}

    return {'schema': 'oph.detector_temporal_records.v1',
        'source': {'binary_states': 66, 'ports': 12, 'clock_region': [0, 1], 'incident_seams': 9,
            'preparation': 'uniform binary total-two source; both local initial bits disclosed',
            'tick': 'each scheduled seam meeting port zero or one, including zero current',
            'physical_seconds_or_modular_time_identified': False},
        'exact_two_field_joint_ones_by_lag': [[[str(p) for p in row] for row in block] for block in joint],
        'centered_single_field': representation(scalar),
        'symmetric_two_field': representation(symmetric),
        'antisymmetric_two_field': representation(antisymmetric),
        'held_out_lags': list(range(15, 21)), 'held_out_recurrence_matches_full_native_source': True,
        'spatial_visibility': {'centered_occupation_dimension': 11, 'visible_dimension': int(visible_rank),
            'dark_dimension': 11-int(visible_rank), 'rank_three_carrier_band_identified_with_temporal_modes': False,
            'dark_transfer_polynomial_ascending': [str(_rational(c)) for c in
                                                  reversed(dark_transfer.charpoly().all_coeffs())],
            'dark_to_linear_two_field_correlations_not_all_record_protocols': True,
            'physical_spatial_coordinates_recovered': False},
        'finite_record_budget': {'histories': sampled, 'simultaneous_joint_probability_tolerance': str(tolerance),
            'conditional_failure_probability_upper': str(_outward(sampling_bound(2048, tolerance), True)),
            'independent_fair_source_services_assumed_not_certified_by_seed': True,
            'raw_seven_mode_spectral_error_certified_at_this_tolerance': perturbation_budget(scalar, tolerance),
            'certified_spectral_error_at_correlation_error_1e_minus_26': str(
                perturbation_budget(scalar, F(1, 10**26))),
            'finite_samples_prove_exact_rank_or_continuation': False},
        'actual_screen_temporal_compression': screen_compression(),
        'scope': {'exact_laws_reconstruct_finite_cyclic_positive_transfer': True,
            'positive_transfer_class_is_existing_native_theorem': True,
            'coherent_real_time_implementation_or_lab_duration_derived': False,
            'same_spectrum_on_refined_screens_asserted': False},
        'source_bindings': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Detector temporal Hilbert space and native correlation predictions replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
