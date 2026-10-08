"""Retained failures of the finite radial inverse at main 1ba8a011.

Controls use original supplied matrices, exact small solutions and scalar
norm identities. They do not infer correctness from the producer's fit flag.
"""

from __future__ import annotations

import math
from dataclasses import asdict

import numpy as np
import pytest

from oph_radial_lift_330 import (
    RadialLiftInputError,
    forward_residual,
    minimum_prior_continuation,
    radial_null_space_report,
)


def test_resolved_direction_is_not_lost_by_squaring_the_condition_number():
    matrix = np.diag([1.0, 1e-8])
    result = minimum_prior_continuation(
        matrix, [1e-3, 1e-11], prior_center=[0.0, 0.0], prior_precision=np.eye(2)
    )
    np.testing.assert_allclose(result.p, [1e-3, 1e-3], rtol=1e-12, atol=0)
    assert result.effective_rank == 2
    np.testing.assert_allclose(result.resolution, np.eye(2), rtol=0, atol=1e-13)
    np.testing.assert_allclose(result.null_projector, np.zeros((2, 2)), rtol=0, atol=1e-13)


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_constraint_units_preserve_resolved_solution_and_projectors(scale):
    result = minimum_prior_continuation(
        scale * np.eye(2), [scale, 2 * scale],
        prior_center=[0.0, 0.0], prior_precision=np.eye(2),
    )
    np.testing.assert_allclose(result.p, [1.0, 2.0], rtol=1e-12, atol=0)
    assert result.objective == pytest.approx(2.5, rel=1e-12, abs=0)
    assert result.effective_rank == 2
    np.testing.assert_allclose(result.null_projector, 0, rtol=0, atol=1e-13)


@pytest.mark.parametrize("units", [(1e-200, 1e200), (1e200, 1e-200), (0.125, 1024.0)])
def test_independent_equation_units_do_not_destroy_a_resolved_identity(units):
    result = minimum_prior_continuation(
        np.diag(units), [units[0], 2 * units[1]],
        prior_center=[0.0, 0.0], prior_precision=np.eye(2),
    )
    np.testing.assert_allclose(result.p, [1, 2], rtol=1e-12, atol=0)
    np.testing.assert_array_equal(result.resolution, np.eye(2))
    np.testing.assert_array_equal(result.null_projector, np.zeros((2, 2)))
    assert result.objective == pytest.approx(2.5, rel=1e-12, abs=0)


@pytest.mark.parametrize("scale", [1e-200, 1e-12, 1.0, 1e200])
def test_inconsistent_redundant_rows_are_not_accepted_at_small_units(scale):
    with pytest.raises(RadialLiftInputError):
        minimum_prior_continuation(
            [[1.0], [1.0]], [scale, 2 * scale],
            prior_center=[0.0], prior_precision=[[1.0]],
        )


def test_a_dominant_equation_cannot_hide_an_impossible_small_zero_row():
    # A global range residual misses the second equation by 100%, while its
    # Euclidean residual is tiny relative to the unrelated first equation.
    with pytest.raises(RadialLiftInputError):
        minimum_prior_continuation(
            [[1, 0], [0, 0]], [1, 1e-100],
            prior_center=[0, 0], prior_precision=np.eye(2),
        )


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_forward_norms_preserve_small_and_large_nonzero_residuals(scale):
    result = forward_residual(np.eye(2), [0.0, 0.0], [scale, 2 * scale])
    assert result["absolute_l2_residual"] == pytest.approx(math.sqrt(5) * scale, rel=1e-14, abs=0)
    assert result["relative_l2_residual"] == pytest.approx(1.0, rel=1e-14, abs=0)


def test_forward_prediction_preserves_a_small_term_between_cancelling_terms():
    result = forward_residual([[2**53, 1, -(2**53)]], [1, 1, 1], [1])
    assert result["predicted"] == [1.0]
    assert result["residual"] == [0.0]
    assert result["absolute_l2_residual"] == 0.0
    assert result["relative_l2_residual"] == 0.0


def test_cancelling_products_can_have_a_representable_forward_result():
    result = forward_residual([[1e200, 1e200]], [1e200, -1e200], [0.0])
    assert result["predicted"] == [0.0]
    assert result["absolute_l2_residual"] == result["relative_l2_residual"] == 0.0


def test_nonzero_residual_has_no_finite_relative_error_against_zero_target():
    with pytest.raises(RadialLiftInputError):
        forward_residual([[1.0]], [1.0], [0.0])


@pytest.mark.parametrize("consumer", ["null", "continuation", "forward"])
@pytest.mark.parametrize("kind", ["mask", "bool", "nan", "infinity", "complex", "rounded_integer"])
def test_original_invalid_operator_entries_are_not_silently_coerced(consumer, kind):
    matrix = [[1.0, 0.0], [0.0, 1.0]]
    if kind == "mask":
        matrix = np.ma.array(matrix, mask=[[False, True], [False, False]])
    elif kind == "bool":
        matrix[0][0] = True
    elif kind == "nan":
        matrix[0][0] = math.nan
    elif kind == "infinity":
        matrix[0][0] = math.inf
    elif kind == "complex":
        matrix[0][0] = 1 + 1j
    else:
        matrix[0][0] = 2**53 + 1
    with pytest.raises(RadialLiftInputError):
        if consumer == "null":
            radial_null_space_report(matrix)
        elif consumer == "continuation":
            minimum_prior_continuation(matrix, [1, 2], prior_center=[0, 0], prior_precision=np.eye(2))
        else:
            forward_residual(matrix, [1, 2], [1, 2])


@pytest.mark.parametrize("rtol", [0.0, -1e-12, 1.0, math.inf, math.nan, True])
@pytest.mark.parametrize("consumer", ["null", "continuation"])
def test_rank_cutoff_must_be_a_finite_fraction(rtol, consumer):
    with pytest.raises(RadialLiftInputError):
        if consumer == "null":
            radial_null_space_report(np.eye(2), rtol=rtol)
        else:
            minimum_prior_continuation(np.eye(2), [1, 2], prior_center=[0, 0],
                                       prior_precision=np.eye(2), rtol=rtol)


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_prior_symmetry_does_not_depend_on_objective_units(scale):
    with pytest.raises(RadialLiftInputError):
        minimum_prior_continuation(
            [[1, 0]], [1], prior_center=[0, 0],
            prior_precision=scale * np.array([[1.0, 0.1], [0.0, 1.0]]),
        )


@pytest.mark.parametrize("field", ["screen_cl", "prior_center", "prior_precision"])
@pytest.mark.parametrize("value", [math.nan, math.inf, True])
def test_all_original_inverse_data_are_validated(field, value):
    arguments = {"screen_cl": [1.0, 2.0], "prior_center": [0.0, 0.0],
                 "prior_precision": [[1.0, 0.0], [0.0, 1.0]]}
    if field == "prior_precision":
        arguments[field][0][0] = value
    else:
        arguments[field][0] = value
    with pytest.raises(RadialLiftInputError):
        minimum_prior_continuation(np.eye(2), **arguments)


def test_zero_map_preserves_the_prior_and_has_no_condition_number():
    import json

    matrix = np.zeros((2, 3))
    report = radial_null_space_report(matrix)
    assert report.rank == 0 and report.nullity == 3
    assert report.condition_number_nonzero is None
    json.dumps(asdict(report), allow_nan=False)
    result = minimum_prior_continuation(matrix, [0, 0], prior_center=[1, 2, 3],
                                        prior_precision=np.eye(3))
    np.testing.assert_array_equal(result.p, [1, 2, 3])
    np.testing.assert_array_equal(result.resolution, np.zeros((3, 3)))
    np.testing.assert_array_equal(result.null_projector, np.eye(3))
    assert result.objective == result.residual_norm == 0
    assert result.effective_rank == 0
    json.dumps(asdict(result), allow_nan=False)


def test_roundoff_is_not_a_second_independent_constraint_at_tiny_rtol():
    matrix = [[1.0, 1.0], [2.0, 2.0]]
    assert radial_null_space_report(matrix, rtol=1e-300).rank == 1
    result = minimum_prior_continuation(matrix, [2, 4], prior_center=[0, 0],
                                        prior_precision=np.eye(2), rtol=1e-300)
    np.testing.assert_allclose(result.p, [1, 1], rtol=1e-12, atol=0)
    assert result.effective_rank == 1
