"""Exact residual, canonical ancestry and native clock covariance controls."""
from fractions import Fraction as F
import json
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

from ..observer.scalar_response.experiment import compare_receipt
from . import ancestry_green as g


def test_level_zero_certificate_against_independent_exact_schur_inverse():
    graph, _, r, f_integer, denominator = g.elimination.matrices(0)
    r = sp.Matrix(r.tolist())
    f = sp.Matrix(f_integer.tolist())/denominator
    schur = r.T*f.inv(method='DM')*r/12
    exact = (schur+sp.ones(20)/20).inv(method='DM')-sp.ones(20)/20
    result = g.covariance(0)
    error_squared = sum((F(exact[i, j])-result['green'][i][j])**2
                         for i in range(20) for j in range(20))
    assert error_squared <= result['green_error']**2
    assert result['green_error'] < F(1, 1000000)
    # Independently recover the root-clock term from exact edge Schur data,
    # without the full-port numerical solution or its encoded differences.
    root_rows = sp.diag(*(int(a//12 == 0 or b//12 == 0)
                          for a, b in graph.inter_seams))
    pairing = r.T*root_rows*f.inv(method='DM')*r*exact/12
    projector = sp.eye(20)-sp.ones(20)/20
    clock = exact+(projector-pairing-pairing.T)/132
    clock_error_squared = sum((F(clock[i, j])-result['incident'][i][j])**2
                               for i in range(20) for j in range(20))
    assert clock_error_squared <= result['incident_error']**2


@pytest.mark.parametrize('level', range(4))
def test_fixed_base_face_dictionary_and_constant_free_covariances(level):
    result = g.covariance(level)
    assert result['ports'] == 240*4**level
    assert result['ancestry_scale'] == 4**level
    for name in ('green', 'attempt', 'incident'):
        matrix = result[name]
        assert all(sum(row) == 0 for row in matrix)
        assert all(matrix[i][j] == matrix[j][i] for i in range(20) for j in range(20))
    assert result['green_error'] < F(1, 10000)
    correction = np.asarray(result['incident'], dtype=float)-np.asarray(result['green'], dtype=float)
    assert np.linalg.norm(correction, 2) <= float(result['incident_clock_correction_bound']
        +result['incident_error']+result['green_error'])
    attempt_delta = np.asarray(result['attempt'], dtype=float)-np.asarray(result['green'], dtype=float)
    assert abs(np.linalg.norm(attempt_delta, 2)-float(result['attempt_clock_correction_bound'])) < 1e-12


def test_exact_residual_detects_a_changed_numerical_solution():
    _, laplacian, score, _, _, inverse = g.problem(0)
    candidate = np.zeros(score.shape)
    original, allowance, residual = g.certify_solution(laplacian, score, candidate, inverse, 1)
    candidate[1, 0] = 1000
    changed, bad_allowance, bad_residual = g.certify_solution(laplacian, score, candidate, inverse, 1)
    assert np.any(changed != original)
    assert bad_residual > residual and bad_allowance > allowance
    assert not np.any(changed.sum(axis=1))


def test_changed_ancestry_is_rejected_instead_of_relabelling_the_readout(monkeypatch):
    builder = g.elimination.source._source()
    raw = builder.build_federation(0, 'port_pair')
    fake = SimpleNamespace(build_federation=lambda level, mode: raw,
        build_geodesic_icosahedral_tower=lambda level: SimpleNamespace(
            cell_ancestor_ids=lambda **kwargs: np.arange(20)[::-1]))
    monkeypatch.setattr(g.elimination.source, '_source', lambda: fake)
    with pytest.raises(ArithmeticError, match='dictionary'):
        g.problem(0)


def test_exact_covariance_comparison_intervals_and_json_scalar_contract():
    result = g.run()
    compare_receipt(json.loads(json.dumps(result)), result)
    for row in result['consecutive_same_readout_comparisons']:
        assert float(F(row['operator_difference_lower'])) <= row['operator_difference_diagnostic']
        assert row['operator_difference_diagnostic'] <= float(F(row['operator_difference_upper']))
        assert row['finite_comparison_does_not_bound_the_unexecuted_refinement_tail']
    assert not result['acceptance']['all_level_Green_Cauchy_tail_or_limit_proved']


@pytest.mark.parametrize('level', [True, -1, 4, 1.0])
def test_source_level_is_a_declared_original_integer(level):
    with pytest.raises(ValueError):
        g.problem(level)


def test_unsafe_dyadic_integer_arithmetic_is_rejected():
    _, laplacian, score, _, _, inverse = g.problem(0)
    with pytest.raises(ValueError, match='allowance'):
        g.certify_solution(laplacian, score, np.full(score.shape, 1e20), inverse, 1)
    with pytest.raises(ValueError):
        g.certify_solution(laplacian, score, np.zeros(score.shape), inverse, 1, bits=True)
    corrupted = F(1); corrupted._numerator = 1.0
    with pytest.raises(ValueError):
        g.certify_solution(laplacian, score, np.zeros(score.shape), corrupted, 1)


@pytest.mark.parametrize('changed', ['operator', 'score', 'score_normalization'])
def test_exact_residual_rejects_noninteger_or_unnormalized_source_data(changed):
    _, laplacian, score, _, _, inverse = g.problem(0)
    if changed == 'operator':
        laplacian = laplacian.astype(float)
    elif changed == 'score':
        score = score.astype(float)
    else:
        score = score.copy(); score[0, 0] = 20
    with pytest.raises(ValueError, match='integer native'):
        g.certify_solution(laplacian, score, np.zeros(score.shape), inverse, 1)
