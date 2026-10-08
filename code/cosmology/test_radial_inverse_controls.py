"""Original-input exact controls for finite minimum-prior radial inversion.

The oracle solves the small rational KKT equations by Gaussian elimination.
It never forms a normal matrix, whitens a metric, or uses a numerical rank or
null-space decomposition.  All supplied examples are exactly dyadic, so the
oracle and the public routine receive the same input values.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np
import pytest

from oph_radial_lift_330 import minimum_prior_continuation


def _solve_exact(matrix, rhs):
    """Small nonsingular rational solve, independent of LAPACK/SVD."""
    n = len(rhs)
    rows = [[F(x) for x in row] + [F(b)] for row, b in zip(matrix, rhs)]
    assert all(len(row) == n + 1 for row in rows)
    for column in range(n):
        pivot = next(i for i in range(column, n) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        divisor = rows[column][column]
        rows[column] = [x / divisor for x in rows[column]]
        for i in range(n):
            if i != column:
                multiple = rows[i][column]
                rows[i] = [x - multiple * y for x, y in zip(rows[i], rows[column])]
    return [row[-1] for row in rows]


def _dot(a, b):
    return sum((F(x) * F(y) for x, y in zip(a, b)), F(0))


def _exact_control(a, c, p0, q):
    """Use explicitly independent constraint rows to get p, R, N, objective."""
    n, m = len(p0), len(c)
    kkt = [list(row) + [a[j][i] for j in range(m)] for i, row in enumerate(q)]
    kkt += [list(row) + [0] * m for row in a]
    p = _solve_exact(kkt, [_dot(row, p0) for row in q] + list(c))[:n]
    columns = [
        _solve_exact(kkt, [0] * n + [a[i][j] for i in range(m)])[:n]
        for j in range(n)
    ]
    resolution = [[columns[j][i] for j in range(n)] for i in range(n)]
    null = [[F(i == j) - resolution[i][j] for j in range(n)] for i in range(n)]
    delta = [x - F(y) for x, y in zip(p, p0)]
    objective = _dot(delta, [_dot(row, delta) for row in q]) / 2
    assert [_dot(row, p) for row in a] == list(c)
    return p, resolution, null, objective


def _assert_components(actual, exact, zero_scales):
    """Nonzero entries get relative checks with no global absolute floor."""
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(exact, dtype=object)
    scales = np.broadcast_to(np.asarray(zero_scales), expected.shape)
    assert actual.shape == expected.shape
    assert np.all(np.isfinite(actual))
    for index in np.ndindex(expected.shape):
        value = float(expected[index])
        if value:
            assert actual[index] == pytest.approx(value, rel=5e-10, abs=0.0), index
        else:
            assert abs(actual[index]) <= 5e-11 * scales[index], index


def _check(a, c, p0, q, *, supplied_a=None, supplied_c=None):
    expected_p, expected_r, expected_n, objective = _exact_control(a, c, p0, q)
    original_a = a if supplied_a is None else supplied_a
    original_c = c if supplied_c is None else supplied_c
    # Also establish compatibility of the redundant original equations exactly.
    assert [_dot(row, expected_p) for row in original_a] == list(original_c)
    actual = minimum_prior_continuation(
        np.asarray(original_a, dtype=float),
        np.asarray(original_c, dtype=float),
        prior_center=np.asarray(p0, dtype=float),
        prior_precision=np.asarray(q, dtype=float),
    )
    qdiag = np.asarray([float(q[i][i]) for i in range(len(q))])
    p_scale = np.sqrt(float(2 * objective) / qdiag)
    projector_scale = np.sqrt(qdiag[None, :] / qdiag[:, None])
    _assert_components(actual.p, expected_p, p_scale)
    _assert_components(actual.resolution, expected_r, projector_scale)
    _assert_components(actual.null_projector, expected_n, projector_scale)
    assert actual.objective == pytest.approx(float(objective), rel=5e-10, abs=0.0)
    assert actual.effective_rank == len(a)

    # Componentwise feasibility is measured in each original equation's units.
    # Computing exact products of the returned floats avoids a second rounded
    # dot product masking a bad producer residual.
    residual = [F(rhs) - _dot(row, actual.p) for row, rhs in zip(original_a, original_c)]
    residual_scale = [
        abs(float(rhs)) + sum(abs(float(x) * y) for x, y in zip(row, actual.p))
        for row, rhs in zip(original_a, original_c)
    ]
    for error, scale in zip(residual, residual_scale):
        assert abs(float(error)) <= 5e-10 * scale
    assert math.isfinite(actual.residual_norm)

    # Exact null columns give independent feasible variations.  Their metric
    # inner product with the displacement must vanish, including for dense Q.
    delta = [F(float(x)) - F(y) for x, y in zip(actual.p, p0)]
    gradient = [_dot(row, delta) for row in q]
    for column in zip(*expected_n):
        terms = [F(x) * y for x, y in zip(column, gradient)]
        assert abs(float(sum(terms))) <= 5e-10 * sum(abs(float(x)) for x in terms)
    return actual


BASE_A = [[1, 2, -1, 0], [0, 1, 1, 2]]
BASE_C = [F(3, 8), F(-5, 4)]
BASE_P0 = [F(1, 4), F(-1, 2), F(1, 8), F(3, 4)]
BASE_Q = [[4, 1, 1, 0], [1, 3, 0, 1], [1, 0, 5, 1], [0, 1, 1, 4]]


@pytest.mark.parametrize("screen_exponent", [-40, 0, 40])
@pytest.mark.parametrize("objective_exponent", [-40, 0, 40])
def test_dense_prior_exact_kkt_in_independent_units(screen_exponent, objective_exponent):
    screen_unit, objective_unit = F(2) ** screen_exponent, F(2) ** objective_exponent
    a = [[screen_unit * x for x in row] for row in BASE_A]
    c = [screen_unit * x for x in BASE_C]
    q = [[objective_unit * x for x in row] for row in BASE_Q]
    _check(a, c, BASE_P0, q)


def test_redundant_constraints_preserve_exact_minimum_and_projectors():
    redundant_a = BASE_A + [[x + 2 * y for x, y in zip(*BASE_A)]]
    redundant_c = BASE_C + [BASE_C[0] + 2 * BASE_C[1]]
    _check(BASE_A, BASE_C, BASE_P0, BASE_Q, supplied_a=redundant_a, supplied_c=redundant_c)


def test_resolved_small_direction_with_dense_prior_matches_exact_kkt():
    small = F(1, 2**24)
    a = [[1, 1, 0], [0, small, small]]
    c = [F(1, 2**12), F(3, 2**40)]
    q = [[4, 1, 0], [1, 3, 1], [0, 1, 2]]
    _check(a, c, [0, 0, 0], q)


def test_anisotropic_prior_keeps_small_nonzero_projector_components():
    a = [[1, 1, 0], [0, 1, 1]]
    q = [[1, 0, 0], [0, 2**16, 0], [0, 0, 2**32]]
    _check(a, [1, -2], [F(1, 4), F(1, 8), F(1, 2)], q)


def test_source_coordinate_units_transform_the_minimum_and_both_projectors():
    # p=S*x changes A to A*S and Q to S^T*Q*S; this is the same
    # constrained problem, with different units for each source coefficient.
    units = [F(1, 2**10), F(1), F(2**10), F(4)]
    a = [[x * units[j] for j, x in enumerate(row)] for row in BASE_A]
    q = [[x * units[i] * units[j] for j, x in enumerate(row)] for i, row in enumerate(BASE_Q)]
    p0 = [x / unit for x, unit in zip(BASE_P0, units)]
    _check(a, BASE_C, p0, q)


@pytest.mark.parametrize("row_exponents", [(-12, 12), (12, -12)])
def test_each_constraint_keeps_its_own_units(row_exponents):
    units = [F(2) ** exponent for exponent in row_exponents]
    a = [[unit * x for x in row] for unit, row in zip(units, BASE_A)]
    c = [unit * x for unit, x in zip(units, BASE_C)]
    _check(a, c, BASE_P0, BASE_Q)
