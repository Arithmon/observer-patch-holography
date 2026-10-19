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
                            (c.native_operator_control, 3), (c.source_patch_control, 4),
                            (c.tent_refinement_control, 5)):
        with pytest.raises(ValueError):
            function(value)


def test_continuous_tent_integrals_and_isotropic_centroid_energy_tensor():
    row = c.tent_geometry_control()
    assert F(row['normalized_tent_mean']) == F(1, 3)
    assert F(row['normalized_tent_second_moment']) == F(1, 6)
    assert F(row['unit_side_tent_gradient_squared']) == 12
    assert F(row['native_limit_tent_energy_quotient_at_coarse_level_zero']) == F(3, 2)


@pytest.mark.parametrize('level', range(5))
def test_all_ancestral_tents_use_actual_dual_edges_and_vanishing_boundary_terms(level):
    row = c.tent_refinement_control(level)
    # Along each coarse boundary there are m edges, each joining tent
    # values1/m and0 for each distinct ancestral tent. The off-diagonal
    # rescaled energy is therefore exactly -1/(20m).
    assert F(row['max_dual_energy_off_diagonal']) == F(1, 20*2**level)
    assert row['source_cells'] == 20*4**level
    assert F(row['exact_normalized_patch_tent_second_moment']) >= F(1, 6)


def test_tent_readouts_approximate_continuous_moments_without_fitting_energy():
    coarse = c.tent_refinement_control(2)
    fine = c.tent_refinement_control(4)
    for key, target in (('exact_normalized_patch_tent_mean', F(1, 3)),
                        ('exact_normalized_patch_tent_second_moment', F(1, 6))):
        assert abs(F(fine[key])-target) < abs(F(coarse[key])-target)/10
    assert F(fine['max_dual_energy_diagonal_error_from_three_twentieths']) < F(1, 100)


def test_exact_full_port_lift_transports_contraction_and_cutoff_with_constant_mode():
    row = c.full_port_form_control()
    assert row['full_port_harmonic_energy_identity_checked']
    assert row['constant_zero_mode_adjoined_before_normal_contraction']
    assert 0 < F(row['internal_norm_squared_over_raw_energy']) < F(1, 2)
    assert 0 < F(row['contracted_energy_over_raw_energy']) < 1
    assert 0 < F(row['contraction_error_squared_over_internal_norm_squared']) < 1
    assert F(row['cutoff_energy']) <= F(row['cutoff_energy_upper'])


def test_changing_the_mean_elimination_kernel_breaks_the_full_native_lift(monkeypatch):
    graph, laplacian, r, integer_f, denominator = c.source.matrices(0)
    changed = integer_f+np.eye(len(integer_f), dtype=np.int64)
    monkeypatch.setattr(c.source, 'matrices', lambda level:
        (graph, laplacian, r, changed, denominator))
    with pytest.raises(ArithmeticError, match='harmonic energy lift'):
        c.full_port_form_control()


@pytest.mark.parametrize('level', (0, 1, 20, 100))
def test_full_port_internal_weight_vanishes_at_every_refinement(level):
    assert c.harmonic_lift_tail(level) == F(1, 2*4**level)
    assert c.harmonic_lift_tail(level+1) == c.harmonic_lift_tail(level)/4


@pytest.mark.parametrize('index', (1, 19, 20, 21, 79, 80, 81, 319, 320, 10**8, 10**100))
def test_global_eigenvalue_enclosures_follow_from_all_level_minmax_indices(index):
    lower, upper = c.continuum_eigenvalue_enclosure(index).values()
    lower_level = 0
    while 20*4**(lower_level+1) <= index:
        lower_level += 1
    actual_lower = (F(3, 4096) if index < 20
                    else 1/c.covariance_eigenvalue_tail(lower_level))
    upper_level = 0
    while c.tent_eigenvalue_bound(upper_level)['index'] < index:
        upper_level += 1
    actual_upper = c.tent_eigenvalue_bound(upper_level)['upper']
    assert lower <= actual_lower <= actual_upper <= upper


def test_regular_three_dimensional_eigen_growth_eventually_violates_native_lower_bound():
    # Exact exponent countercontrol for n^(2/3), with unit coefficient.
    # The theorem treats every fixed positive coefficient asymptotically;
    # this is not a fit to an actual laboratory spectrum.
    side = 49921
    index = side**3
    assert side**2 < c.continuum_eigenvalue_enclosure(index)['lower']


@pytest.mark.parametrize('index', (True, 0, -1, 1., '2'))
def test_continuum_eigenvalue_index_rejects_invalid_values(index):
    with pytest.raises(ValueError):
        c.continuum_eigenvalue_enclosure(index)


@pytest.mark.parametrize('level', (True, -1, 1., '2'))
def test_all_level_form_and_tent_bounds_reject_invalid_levels(level):
    for function in (c.harmonic_lift_tail, c.tent_eigenvalue_bound, c.tent_refinement_control):
        with pytest.raises(ValueError):
            function(level)


def test_exact_interpolation_integrals_and_both_source_cycle_certificates():
    row = c.interpolation_geometry_control()
    assert F(row['normalized_gradient_times_side_squared_per_dual_energy_upper']) == F(7, 30)
    assert F(row['normalized_mass_error_times_cell_count_per_dual_energy_upper']) == F(5, 18)
    for certificate, zero_count in zip(row['source_vertex_cycle_certificates'], (1, 3), strict=True):
        pivots = [F(x) for x in certificate['rational_nonnegative_pivots']]
        assert all(value >= 0 for value in pivots)
        assert pivots.count(0) == zero_count


@pytest.mark.parametrize('level', range(4))
def test_source_interpolation_uses_closed_native_stars_and_exact_energy_mass_bounds(level):
    row = c.source_interpolation_control(level)
    assert row['vertex_valence_counts'].get('5') == 12
    assert row['vertex_valence_counts'].get('6', 0) == 10*4**level-10
    assert F(row['normalized_gradient_times_side_squared']) <= F(row['normalized_gradient_times_side_squared_upper'])
    assert 0 < F(row['exact_normalized_interpolation_mass_error']) <= F(row['normalized_interpolation_mass_error_upper'])
    compare_receipt(json.loads(json.dumps(row)), row)


def test_missing_native_edge_fails_the_interpolation_dictionary_binding(monkeypatch):
    original = c.uniform.dual_edges
    monkeypatch.setattr(c.uniform, 'dual_edges', lambda faces: original(faces)[1:])
    with pytest.raises(ArithmeticError, match='dual-edge dictionaries differ'):
        c.source_interpolation_control(1)


def test_source_interpolation_does_not_silently_insert_a_boundary_condition(monkeypatch):
    vertices, faces = c.uniform.affine_prefix(1)
    monkeypatch.setattr(c.uniform, 'affine_prefix', lambda level: (vertices, faces[1:]))
    with pytest.raises(ArithmeticError, match='two cells at each primal edge'):
        c.source_interpolation_control(1)


def test_valence_seven_would_invalidate_the_specific_cycle_variance_constant():
    values = (10, 6, -2, -9, -9, -2, 6)
    assert sum(values) == 0
    variance = sum(x*x for x in values)
    energy = sum((values[i]-values[(i+1) % 7])**2 for i in range(7))
    assert variance > energy
    with pytest.raises(ValueError, match='valence five or six'):
        c.cycle_variance_control(7)


@pytest.mark.parametrize('valence', (True, 4, 5., '6'))
def test_cycle_variance_control_requires_an_actual_source_valence(valence):
    with pytest.raises(ValueError):
        c.cycle_variance_control(valence)


def test_interpolation_finite_control_size_is_bounded():
    with pytest.raises(ValueError, match='levels zero through three'):
        c.source_interpolation_control(4)


@pytest.mark.parametrize('degree', (1, 2, 3, 7))
def test_ball_harmonic_extension_energy_is_linear_in_degree_not_quadratic(degree):
    # Independent exact radial integral for a spherical harmonic of unit
    # angular L2 norm and angular energy j(j+1). Boundary area contributes R².
    r, radius = c.sp.symbols('r radius', positive=True)
    radial = (r/radius)**degree
    energy = c.sp.integrate(r*r*c.sp.diff(radial, r)**2
                           +degree*(degree+1)*radial**2, (r, 0, radius))
    assert c.sp.simplify(energy/radius**2-degree/radius) == 0


def test_positive_ball_steklov_index_retains_every_angular_multiplicity():
    import math

    spectrum = [degree for degree in range(1, 21) for _ in range(2*degree+1)]
    assert len(spectrum) == 20*22
    assert spectrum == [math.isqrt(index) for index in range(1, len(spectrum)+1)]


@pytest.mark.parametrize('scale,bounded_correction', [(F(1), F(0)), (F(7, 3), F(5)),
                                                     (F(1, 1000), F(999))])
def test_native_precision_eventually_exceeds_every_declared_ball_trace_upper_bound(scale, bounded_correction):
    # scale includes c*||T||²/R for a bounded form-congruence map. This is
    # an exact all-level bound witness, not a computed high-frequency spectrum.
    threshold = 49920*(scale+bounded_correction+1)
    degree = -((-threshold.numerator)//threshold.denominator)
    index = degree*degree
    native_lower = c.continuum_eigenvalue_enclosure(index)['lower']
    assert native_lower > scale*degree+bounded_correction
