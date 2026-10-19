"""Independent triangular compactness and native action-range controls."""
from collections import Counter
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import local_continuum as c
from ..observer.scalar_response.experiment import compare_receipt


@pytest.mark.parametrize('depth', range(6))
def test_exact_reflection_multiplicities_and_dual_patch_census(depth):
    result = c.reflected_square_control(depth)
    m, vertices, edges = c.triangle_patch(depth)
    assert len(vertices) == m*m
    assert len(edges) == 3*m*(m-1)//2
    assert result['square_edges'] == 2*m*(m-1)
    degrees = Counter(i for edge in edges for i in edge)
    for i, (kind, _, _) in enumerate(vertices):
        if kind == 'd':
            assert degrees[i] == 3


@pytest.mark.parametrize('depth', range(4))
def test_actual_native_ancestry_has_the_patch_stencil(depth):
    assert c.source_patch_control(depth)['base_patches_checked'] == 20


@pytest.mark.parametrize('depth', range(3))
def test_patch_bound_has_independent_exact_positive_matrix_certificate(depth):
    row = c.exact_patch_energy_control(depth)
    assert row['rational_LDL_positive_pivots'] == 4**depth


def test_patch_variance_is_neumann_variance_without_pinned_boundary_values():
    # Nonzero, nonconstant boundary values exercise the actual free patch.
    m, vertices, edges = c.triangle_patch(3)
    values = [F(7*i*i-3*j+int(kind == 'd')*11) for kind, i, j in vertices]
    mean = sum(values)/len(values)
    variance = sum((value-mean)**2 for value in values)
    energy = sum((values[a]-values[b])**2 for a, b in edges)
    assert 0 < variance < 26*m*m*energy
    translated = [value+17 for value in values]
    shifted_mean = sum(translated)/len(translated)
    assert sum((value-shifted_mean)**2 for value in translated) == variance


def test_missing_internal_edge_cannot_pass_the_native_patch_binding(monkeypatch):
    original = c.uniform.dual_edges
    def omit(faces):
        edges = original(faces)
        # At level one this edge is internal to the first coarse triangle.
        victim = next(edge for edge in edges if edge[0] < 4 and edge[1] < 4)
        return tuple(edge for edge in edges if edge != victim)
    monkeypatch.setattr(c.uniform, 'dual_edges', omit)
    with pytest.raises(ArithmeticError, match='patch stencil'):
        c.source_patch_control(1)


@pytest.mark.parametrize('level', range(3))
def test_native_source_polynomial_range_and_compactness(level):
    row = c.native_operator_control(level)
    compare_receipt(json.loads(json.dumps(row)), row)
    assert row['absolute_rescaled_coupling_diagnostic'] <= float(F(row['absolute_rescaled_coupling_upper']))
    assert all(value['integer_polynomial_support_checked'] for value in row['exact_polynomial_support'])
    assert row['compactness_controls'][-1]['tail_norm_squared_diagnostic'] == 0


def test_nonlocal_internal_kernel_is_rejected_by_polynomial_support(monkeypatch):
    original = c.source.matrices(0)
    graph, laplacian, r, integer_f, denominator = original
    changed = integer_f.copy()
    # Find a pair of cross seams with disjoint endpoint cells.
    for i in range(len(r)):
        for j in range(i+1, len(r)):
            if not np.any((r[i] != 0) & (r[j] != 0)):
                changed[i, j] += 1; changed[j, i] += 1
                monkeypatch.setattr(c.source, 'matrices', lambda level:
                    (graph, laplacian, r, changed, denominator))
                with pytest.raises(ArithmeticError, match='long-range'):
                    c.native_operator_control(0)
                return
    raise AssertionError('the native seam graph has no disjoint seam pair')


def test_schur_mean_elimination_does_not_assume_a_Markov_matrix():
    row = c.non_markov_control()
    assert F(row['exact_positive_Schur_off_diagonal']) > 0
    assert not row['finite_cell_mean_generator_is_a_Markov_generator']
    assert not row['continuum_locality_or_nonlocality_decided_by_this_control']


def test_exact_integer_polynomial_rejects_powers_that_would_overflow(monkeypatch):
    graph, laplacian, r, integer_f, denominator = c.source.matrices(0)
    # The rational matrix is unchanged, but this unreduced representation
    # cannot safely execute the exact support powers in fixed-width integers.
    monkeypatch.setattr(c.source, 'matrices', lambda level:
        (graph, laplacian, r, integer_f*10**6, denominator*10**6))
    with pytest.raises(ArithmeticError, match='overflow int64'):
        c.native_operator_control(0)


@pytest.mark.parametrize('coarse', (0, 1, 10, 100))
def test_all_level_exact_tail_bounds_need_no_graph_allocation(coarse):
    assert c.compactness_tail_squared(coarse) == F(851968, 4**coarse)
    assert c.covariance_eigenvalue_tail(coarse) == F(624, 4**coarse)
    assert c.compactness_tail_squared(coarse+1) == c.compactness_tail_squared(coarse)/4


def test_rescaled_separated_action_and_macroscopic_decay():
    assert c.separated_action_bound(3, 5) == 2
    exponents = [c.macroscopic_action_exponent(k, F(1, 4)) for k in (5, 8, 12)]
    assert exponents == [6, -44, -996]
    # Exact ceiling for a noninteger separation, including the cell-size margin.
    assert c.macroscopic_action_exponent(5, F(1, 3)) == 3
    assert c.macroscopic_action_exponent(100, F(1, 4)) < -10**20


@pytest.mark.parametrize('value', (True, -1, 1., '2'))
def test_level_parameters_require_original_nonnegative_integers(value):
    for function in (c.compactness_tail_squared, c.covariance_eigenvalue_tail, c.triangle_patch):
        with pytest.raises(ValueError):
            function(value)


@pytest.mark.parametrize('distance', (True, 0, -1, 1., '2'))
def test_support_distance_requires_a_positive_integer(distance):
    with pytest.raises(ValueError):
        c.separated_action_bound(0, distance)


@pytest.mark.parametrize('separation', (0, .25, F(0), F(-1, 4), True))
def test_comparison_separation_requires_an_exact_positive_fraction(separation):
    with pytest.raises(ValueError):
        c.macroscopic_action_exponent(8, separation)


def test_coarse_cells_do_not_get_an_unjustified_separation_bound():
    with pytest.raises(ValueError, match='not yet certified disjoint'):
        c.macroscopic_action_exponent(0, F(1, 4))


def test_finite_control_allocation_caps_are_explicit():
    for function, value in ((c.triangle_patch, 6), (c.exact_patch_energy_control, 4),
                            (c.native_operator_control, 3), (c.source_patch_control, 4)):
        with pytest.raises(ValueError):
            function(value)
