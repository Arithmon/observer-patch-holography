"""Certified finite coordinate-sign projectors and complete-history errors.

The eigensolver supplies only a candidate basis. Exact integer residuals in
the factorial Bargmann metric certify its proximity to the source projector.
No eigensolver accuracy assumption, fitted convergence rate or output
probability enters the certificate. See READOUT.md for the polar-basis proof.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import math

from . import halfline


def _dimension(value):
    if type(value) is not int or value < 2 or value > 128 or value % 2:
        raise ValueError('an even integer pointer dimension from 2 through 128 required')
    return value


def _upper_root(value):
    return halfline._sqrt_upper(value)


@lru_cache(maxsize=8)
def _seed(d):
    """Reproducible approximate eigenvectors, rounded before certification."""
    import mpmath as mp
    d = _dimension(d)
    basis_scale, value_scale = 10**(d+30), 10**30
    with mp.workdps(d+60):
        matrix = mp.matrix(d)
        for n in range(1, d):
            matrix[n-1, n] = matrix[n, n-1] = mp.sqrt(n)
        values, vectors = mp.eigsy(matrix)
        basis = tuple(tuple(int(mp.nint(vectors[n, j]*basis_scale/mp.sqrt(math.factorial(n))))
                            for j in range(d)) for n in range(d))
        diagonal = tuple(int(mp.nint(values[j]*value_scale)) for j in range(d))
    return basis, basis_scale, diagonal, value_scale


def certify_eigenbasis(basis, basis_scale, diagonal, value_scale):
    """Independently certify a supplied rational candidate against a+a*.

    Columns use the unnormalized Bargmann basis. This deliberately accepts
    arbitrary candidate entries: a corrupted eigensolver is either enclosed
    with its real error or rejected, never trusted because it ran normally.
    """
    if type(basis) is not tuple:
        raise ValueError('immutable integer eigenbasis rows required')
    d = _dimension(len(basis))
    if (any(type(row) is not tuple or len(row) != d
            or any(type(x) is not int for x in row) for row in basis)
            or type(diagonal) is not tuple or len(diagonal) != d
            or any(type(x) is not int or x == 0 for x in diagonal)
            or type(basis_scale) is not int or basis_scale <= 0
            or type(value_scale) is not int or value_scale <= 0):
        raise ValueError('square integer basis, nonzero integer spectrum and positive integer scales required')
    weights = tuple(math.factorial(n) for n in range(d))
    scale2 = basis_scale**2
    gram_squared = 0
    for i in range(d):
        for j in range(i, d):
            entry = sum(weights[n]*basis[n][i]*basis[n][j] for n in range(d))
            entry -= scale2 if i == j else 0
            gram_squared += (1 if i == j else 2)*entry**2
    epsilon = _upper_root(F(gram_squared, scale2**2))
    if epsilon >= 1:
        raise ArithmeticError('candidate columns do not certify an invertible near-orthonormal basis')
    residual_squared = 0
    for n in range(d):
        for j in range(d):
            applied = ((basis[n-1][j] if n else 0)
                       + ((n+1)*basis[n+1][j] if n+1 < d else 0))
            error = applied*value_scale-basis[n][j]*diagonal[j]
            residual_squared += weights[n]*error**2
    residual = _upper_root(F(residual_squared, (basis_scale*value_scale)**2))
    maximum = F(max(map(abs, diagonal)), value_scale)
    gap = F(min(map(abs, diagonal)), value_scale)
    # ||a+a*|| <= 2 sqrt(d-1) <= 2d. Polar correction V of U satisfies
    # ||U-V||_F <= ||U*U-I||_F = epsilon < 1 in the physical metric.
    operator_residual = residual+(2*d+maximum)*epsilon
    # In the eigenbases of A and C=V Lambda V*, differing-sign entries
    # of P_A-P_C are bounded by entries of A-C divided by gap.
    projection_error = operator_residual/gap+epsilon*(2+epsilon)
    if projection_error >= F(1, 4):
        raise ArithmeticError('candidate does not certify a useful finite sign projector')
    return {
        'dimension': d,
        'basis_orthogonality_frobenius_error_upper': str(epsilon),
        'source_eigenvector_residual_frobenius_upper': str(residual),
        'candidate_spectral_gap_lower': str(gap),
        'polar_source_residual_frobenius_upper': str(operator_residual),
        'positive_projector_operator_error_upper': str(projection_error),
        'source_operator': 'a+a_dagger in the factorial Bargmann metric; same sign as Q_d',
        'eigensolver_accuracy_assumed': False,
    }


@lru_cache(maxsize=8)
def _certified_pointer(d):
    basis, scale, diagonal, value_scale = _seed(_dimension(d))
    certificate = certify_eigenbasis(basis, scale, diagonal, value_scale)
    positive = tuple(j for j, value in enumerate(diagonal) if value > 0)
    weights = tuple(math.factorial(n) for n in range(d))
    matrix = tuple(tuple(F(weights[column]*sum(basis[row][j]*basis[column][j]
                                              for j in positive), scale**2)
                         for column in range(d)) for row in range(d))
    return matrix, certificate


def pointer_certificate(d):
    """A fresh compact certificate; cached internal evidence cannot be edited."""
    return dict(_certified_pointer(_dimension(d))[1])


def apply_pointer(vector, d, outcome):
    """Apply the certified rational approximation without renormalizing output."""
    d, vector = _dimension(d), halfline._vector(vector)
    if type(outcome) is not int or outcome not in (0, 1):
        raise ValueError('binary integer pointer outcome required')
    if any(n0 >= d or n1 >= d for n0, n1, _ in vector):
        raise ValueError('input must belong to the declared rectangular occupation code')
    matrix = _certified_pointer(d)[0]
    groups = {}
    for (n0, n1, gauge), value in vector.items():
        groups.setdefault((n1, gauge), []).append((n0, value))
    output = {}
    for (n1, gauge), entries in groups.items():
        for row in range(d):
            real = imaginary = F(0)
            for column, (a, b) in entries:
                entry = matrix[row][column]
                if outcome == 0:
                    entry = (1 if row == column else 0)-entry
                real += entry*a
                imaginary += entry*b
            if real or imaginary:
                output[row, n1, gauge] = real, imaginary
    return output


def complete_cutoff_comparison(history, d):
    """Enclose ideal finite/infinite complete outputs for both feedback modes.

    ``history`` is an internally recomputed infinite_history.certificate with
    retained polynomial vectors, not a verifier for untrusted JSON evidence.
    The public receipt checker independently regenerates it and every field.
    The native preparation current is zero, matching those actual vectors.
    """
    from .occupation import occupation_cutoff_bound
    from .sharp_output import pure_output_distance_interval
    from .infinite_history import _interval
    d = _dimension(d)
    normalization = F(history['initial_polynomial_norm_squared'])
    pointer = pointer_certificate(d)
    eta = F(pointer['positive_projector_operator_error_upper'])
    evolution = occupation_cutoff_bound(d)
    evolution_error = F(evolution['pre_final_readout_trace_distance_upper'])
    results = []
    for configuration in history['configurations']:
        lower = upper = mass = F(0)
        records = []
        for branch in configuration['retained_first_records']:
            vector = {(n0, n1, gauge): (F(real), F(imag))
                      for n0, n1, gauge, real, imag in branch['pre_readout_polynomial']}
            mass += halfline.norm_squared(vector)/normalization
            for outcome in (0, 1):
                projected = apply_pointer(vector, d, outcome)
                lo, hi = pure_output_distance_interval(vector, projected, outcome, normalization, normalization)
                lower += lo
                upper += hi
                records.append({'first_record': branch['first_record'], 'sign_record': outcome,
                                'approximate_output_half_trace_norm_interval': _interval(lo, hi)})
        # Both record branches have vector error <= sqrt(2)*eta*sqrt(mass).
        # Half trace norm <= mass*(sqrt(2)*eta+eta**2), enlarged below.
        pointer_error = mass*(2*eta+eta*eta)
        polynomial_error = F(configuration['complete_output_half_trace_norm_error_upper'])
        comparison_error = evolution_error+2*polynomial_error+pointer_error
        results.append({
            'feedback_enabled': configuration['feedback_enabled'],
            'retained_record_comparisons': records,
            'polynomial_instrument_distance_interval': _interval(max(F(0), lower-pointer_error),
                                                                upper+pointer_error),
            'pointer_seed_density_error_upper': str(pointer_error),
            'infinite_history_polynomial_error_upper': str(polynomial_error),
            'complete_finite_to_infinite_trace_distance_interval': _interval(
                max(F(0), lower-comparison_error), min(F(1), upper+comparison_error)),
        })
    return {'dimension': d, 'pointer_certificate': pointer,
            'pre_readout_trace_distance_upper': str(evolution_error),
            'configurations': results,
            'finite_continuous_to_infinite_instrument': True,
            'source_execution_at_this_dimension_claimed': False,
            'includes_source_splitting_or_physical_calibration': False}
