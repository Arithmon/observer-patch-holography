"""Actual source-feedback family coefficient, Higgs identities and exclusion."""
import pytest
import sympy as s

from . import native_repair_mass as source
from . import native_repair_mass_attachment as m


@pytest.fixture(scope='module')
def data():
    return source.architecture()


@pytest.fixture(scope='module')
def attached(data):
    return m.report(data['Y'], projectors=data['projectors'])


def test_frozen_source_roots_are_positive_ordered_and_certified(attached):
    report = attached['source_spectrum']
    assert report['three_sign_changes_exhaust_the_degree_three_polynomial']
    assert report['empirical_inputs'] == []
    assert report['e2_exact'] == str(m.E2)
    intervals = [list(map(s.Rational, row)) for row in report['positive_eigenvalue_intervals_exact']]
    assert 0 < intervals[0][0] < intervals[0][1] < intervals[1][0] < intervals[1][1] < intervals[2][0] < intervals[2][1] < 1


def test_fixed_y_is_shared_across_every_gauge_invariant_channel(data, attached):
    rows = m.higgs_kernels(data['Y'])
    current = m.matter.representation()['matter_currents'][0]
    lifted = m.matter.representation()['lifted'][0]
    g = s.kronecker_product(current, s.eye(3))
    for row in rows:
        h = -lifted[3:, 3:].T if row['conjugate'] else lifted[3:, 3:]
        for k in (0, 1):
            c = row['components'][k]
            assert m._zero(g.T*c+c*g+sum((h[j, k]*row['components'][j] for j in (0, 1)), s.zeros(45)))
    report = attached['same_current_Higgs_attachment']
    assert report['complete_Higgs_component_identities_exact'] == 72
    assert report['neutral_K45_rank'] == 42 and report['neutral_K45_zero_dimension'] == 3
    assert report['equal_mass_degeneracy_across_up_down_and_charged_lepton_channels']
    assert report['quark_charged_current_family_mixing'].startswith('I3')
    assert not report['independent_species_strengths_or_empirical_yukawas_inserted']


def test_entire_geometric_memory_family_is_excluded_without_probability_tuning(attached):
    shape = attached['whole_probability_shape']
    comparison = attached['separate_heldout_comparison']
    assert s.simplify(s.sympify(shape['determinant_over_e2_cubed_lower_exact'])-(25-5*s.sqrt(5))/32) == 0
    assert comparison['fixed_source_ratio_intersects_each_window'] == [False, False]
    assert comparison['all_probabilities_and_both_fivefold_classes_excluded_for_joint_windows']
    assert not comparison['other_interacting_repair_laws_excluded']
    assert not comparison['windows_are_experimental_standard_uncertainties']


def test_neutral_family_copy_matches_the_exact_finite_family_band(attached):
    report = attached['native_family_band_theorem_connection']
    V = s.Matrix([[s.sympify(x) for x in row] for row in report['normalized_coefficient_isometry12_by3_exact']])
    assert V.shape == (12, 3) and m._zero(V.T*V-s.eye(3))
    assert report['actual_rotation_intertwining_count'] == 60
    assert report['actual_family_rotation_kernel_order'] == 1
    assert report['full_coefficient_band_dimensions_retained'] == [1, 3, 3, 5]
    assert not report['separate_neutral_family_copy_is_the_colored_current_triplet']
    assert not report['physical_family_number_forced_by_gauge_current_alone']


def test_comparison_data_are_not_used_for_source_or_higgs_construction(data, monkeypatch):
    original = m.strict_json.load
    def without_empirical(path, *args, **kwargs):
        if path == m.ROOT/m.BENCHMARK or str(path).endswith('calibrated_matter_parameters.json'):
            raise AssertionError('empirical packet read during source construction')
        return original(path, *args, **kwargs)
    monkeypatch.setattr(m.strict_json, 'load', without_empirical)
    m.spectral_certificate(data['Y'])
    m.gauge_attachment(data['Y'])
    m.whole_probability_shape_certificate(data['projectors'])


def test_cubic_certificate_rejects_a_changed_coefficient(data):
    with pytest.raises(ArithmeticError, match='polynomial'):
        m.spectral_certificate(data['Y']+s.eye(3)/1000)
    with pytest.raises(ValueError): m.spectral_certificate(data['Y'].evalf())
    for deficient in (s.zeros(3), s.diag(1, 0, 0)):
        with pytest.raises(ArithmeticError, match='polynomial'):
            m.gauge_attachment(deficient)


def test_shape_certificate_requires_all_actual_distinct_projectors(data):
    with pytest.raises(ValueError): m.whole_probability_shape_certificate(data['projectors'][:4])
    damaged = list(data['projectors']); damaged[1] = damaged[0]
    with pytest.raises(ArithmeticError, match='equiangular'):
        m.whole_probability_shape_certificate(damaged)
