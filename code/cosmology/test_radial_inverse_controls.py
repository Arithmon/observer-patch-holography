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

from oph_radial_lift_330 import RadialLiftInputError, minimum_prior_continuation


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


def _assert_components(actual, exact, zero_scales, *, rtol=5e-10):
    """Nonzero entries get relative checks with no global absolute floor."""
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(exact, dtype=object)
    scales = np.broadcast_to(np.asarray(zero_scales), expected.shape)
    assert actual.shape == expected.shape
    assert np.all(np.isfinite(actual))
    for index in np.ndindex(expected.shape):
        value = float(expected[index])
        if value:
            assert actual[index] == pytest.approx(value, rel=rtol, abs=0.0), index
        else:
            assert abs(actual[index]) <= 5e-11 * scales[index], index


def _check(a, c, p0, q, *, supplied_a=None, supplied_c=None, component_rtol=5e-10):
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
    _assert_components(actual.p, expected_p, p_scale, rtol=component_rtol)
    _assert_components(actual.resolution, expected_r, projector_scale, rtol=component_rtol)
    _assert_components(actual.null_projector, expected_n, projector_scale, rtol=component_rtol)
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
    assert len(actual.residual) == len(residual)
    for reported, exact, scale in zip(actual.residual, residual, residual_scale):
        assert math.isfinite(reported)
        assert abs(reported - float(exact)) <= 5e-15 * scale
    exact_norm = math.hypot(*(float(x) for x in residual))
    assert abs(actual.residual_norm - exact_norm) <= 5e-15 * math.hypot(*residual_scale)

    # Exact null columns give independent feasible variations.  Their metric
    # inner product with the displacement must vanish, including for dense Q.
    delta = [F(float(x)) - F(y) for x, y in zip(actual.p, p0)]
    gradient = [_dot(row, delta) for row in q]
    returned_objective = _dot(delta, gradient) / 2
    assert actual.objective == pytest.approx(float(returned_objective), rel=5e-13, abs=0.0)
    for column in zip(*expected_n):
        terms = [F(x) * y for x, y in zip(column, gradient)]
        assert abs(float(sum(terms))) <= 5e-10 * sum(abs(float(x)) for x in terms)
    for i, r_column in enumerate(zip(*actual.resolution)):
        for j, n_column in enumerate(zip(*actual.null_projector)):
            terms = [F(x) * _dot(row, n_column) for x, row in zip(r_column, q)]
            assert abs(float(sum(terms))) <= 5e-10 * math.sqrt(qdiag[i] * qdiag[j])
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


def test_full_column_rank_with_redundant_data_has_zero_null_projector():
    a = [[1, 2], [2, -1]]
    c = [F(5, 8), F(-5, 4)]
    _check(
        a, c, [F(1, 2), F(-1, 4)], [[3, 1], [1, 2]],
        supplied_a=a + [[4, 3]], supplied_c=c + [2 * c[0] + c[1]],
    )


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


def _accurate_or_resolution_refusal(a, c, p0, q):
    # These exact compatible problems may exceed the solver's reportable
    # precision.  A numerical-resolution refusal must not become a false
    # minimum or an original-constraint null-space claim.
    try:
        _check(a, c, p0, q)
    except RadialLiftInputError as error:
        assert any(word in str(error).lower() for word in ("resol", "round", "report", "represent", "numerical null"))


def test_large_cancelling_prior_cannot_hide_complete_loss_of_correction():
    _accurate_or_resolution_refusal(
        [[1, 1]], [F(1, 2**20)], [2**40, -(2**40)], [[2, 1], [1, 3]],
    )


def test_rounding_the_continuation_cannot_halve_the_minimum_objective():
    _accurate_or_resolution_refusal(
        [[1, 1]], [2**53 + 2], [2**53, 0], [[1, 0], [0, 1]],
    )


def test_resolved_correction_survives_a_large_cancelling_prior():
    # The optimum's two increments are exact multiples of the local ulp.
    # This nearby case must remain supported, including its small objective.
    _check([[1, 1]], [F(3, 2**20)], [2**20, -(2**20)], [[2, 1], [1, 3]])


def test_prior_weighting_cannot_create_an_original_operator_null_direction():
    _accurate_or_resolution_refusal(
        [[1, 0], [0, 1]], [1, 0], [0, 0], [[1, 0], [0, 2**96]],
    )


def test_equal_raw_and_weighted_rank_counts_do_not_validate_the_null_subspace():
    _accurate_or_resolution_refusal(
        [[1, 0], [0, F(1, 2**48)]], [0, F(1, 2**48)], [0, 0],
        [[2**96, 0], [0, F(1, 2**96)]],
    )


def test_prior_weighting_can_resolve_a_small_original_operator_direction():
    _check(
        [[1, 0], [0, F(1, 2**48)]], [1, F(1, 2**48)], [0, 0],
        [[1, 0], [0, F(1, 2**96)]],
    )


def test_parallel_whitened_rows_do_not_make_an_original_constraint_null():
    _accurate_or_resolution_refusal(
        [[1, 0], [1, 1]], [1, 1], [0, 0], [[1, 0], [0, 2**96]],
    )


@pytest.mark.parametrize("exponent", [40, 48])
def test_correlated_prior_weak_mode_requires_resolved_minimum(exponent):
    # Q is exactly (3,5)^T(3,5) plus a positive second-coordinate term.
    # Choosing the row perturbation as its square root gives objective 1/4.
    # A feasible p alone cannot certify that the dense prior was respected.
    weak = F(1, 2**exponent)
    row_step = F(1, 2 ** (exponent // 2))
    _accurate_or_resolution_refusal(
        [[3, 5 + row_step]], [1], [0, 0], [[9, 15], [15, 25 + weak]],
    )


def test_ordinary_correlated_prior_matches_exact_minimum():
    _check([[3, 5 + F(1, 2**5)]], [1], [0, 0], [[9, 15], [15, 25 + F(1, 2**10)]])


def test_resolved_correlated_prior_keeps_its_weak_mode():
    # This intermediate case remains inside the declared 1e-7 resolution
    # budget, while the exact KKT values permit a tighter component check.
    _check(
        [[3, 5 + F(1, 2**10)]], [1], [0, 0], [[9, 15], [15, 25 + F(1, 2**20)]],
        component_rtol=1e-8,
    )


def test_near_parallel_constraints_need_a_resolved_forward_solution():
    _accurate_or_resolution_refusal(
        [[1, 1], [1, 1 + F(1, 2**36)]], [1, 1 + F(1, 2**38)],
        [0, 0], [[1, 0], [0, 1]],
    )


def test_resolved_near_parallel_constraints_keep_both_directions():
    _check(
        [[1, 1], [1, 1 + F(1, 2**18)]], [1, 1 + F(1, 2**20)],
        [0, 0], [[1, 0], [0, 1]],
    )


def test_anisotropic_null_projector_preserves_a_tiny_positive_component():
    # The exact N[0,0] is 1/(1+2**96), despite R[0,0] rounding to 1.
    _check([[1, 1]], [1], [0, 0], [[1, 0], [0, 2**96]])
