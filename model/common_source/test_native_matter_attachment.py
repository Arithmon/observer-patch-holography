"""Same-current matter lifting, complete native records and custody controls."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import native_matter_attachment as m


@pytest.fixture(scope='module')
def built():
    return m.build()


def test_auxiliary_fock_car_and_all_sectors_are_retained():
    creation = m.fock_data()
    assert len(creation) == 5
    for i, a in enumerate(creation):
        for j, b in enumerate(creation):
            assert a*b+b*a == s.zeros(32)
            assert a.H*b+b*a.H == int(i == j)*s.eye(32)
    r = m.representation()
    assert len(r['selected']) == 15
    assert 0 not in r['selected'] and 31 not in r['selected']
    assert r['encoder'].T*r['encoder'] == s.eye(15)
    assert all(g.shape == (32, 32) for g in r['full_currents'])


def test_actual_native_center_and_spin_lift_are_not_an_inherited_so3_subrepresentation():
    r = m.representation()
    assert r['primitive_charge5'] == s.diag(-2, -2, -2, 3, 3)
    assert r['center5'] == s.diag(1, 1, 1, -s.Rational(3, 2), -s.Rational(3, 2))
    z = m.canonical_currents()[10]
    assert np.max(abs(expm(2*np.pi*np.asarray(z, complex))-np.eye(6))) < 1e-13
    assert np.max(abs(expm(2*np.pi*np.asarray(m.spin_lift(z), complex))-np.diag([1, 1, 1, -1, -1]))) < 1e-13


def test_same_port_matter_intertwiner_and_physical_car_frame_are_exact():
    p, vertices = m.action_intertwiner()
    r = m.representation(); old = m.existing_action._representation()
    assert p.T*p == p*p.T == s.eye(15)
    for port, current in enumerate(r['matter_currents']):
        assert m.zero(old[2][r['native']['native_to_fixture_port'][port]]*p-p*current)
    physical = s.kronecker_product(s.eye(2), p, s.eye(3))
    assert physical.shape == (90, 90) and physical.T*physical == s.eye(90)
    assert sorted(row['neutral_rank'] for row in vertices) == [2, 6, 6]


def test_same_kinetic_and_mass_action_is_transferred_without_replacing_its_spectrum(built):
    report = built[0]['same_kinetic_and_mass_action']
    assert report['same_existing_uncalibrated_action_kernel_exact']
    assert report['arbitrary_complex_family_action_intertwines_exactly']
    assert report['independent_complex_family_symbols'] == 27
    assert report['both_Higgs_components_and_all_offdiagonal_family_entries_checked']
    assert not report['calibrated_or_uncalibrated_parameter_producer_executed_for_this_identity']
    assert report['all_three_Weyl_principal_symbol_axes_intertwine_exactly']
    assert report['physical_CAR_frame_is_unitary_and_pairing_antisymmetric']
    assert not report['existing_action_hierarchy_failure_removed_by_frame_change']
    assert not report['mass_parameters_refitted_in_this_handoff']


def test_same_twelve_currents_allow_distinct_exact_positive_relative_spectra(built):
    report = built[0]['relative_spectrum_from_same_current_constraints']
    assert report['Higgs_invariant_line_count'] == 3
    assert report['complex_Yukawa_coefficient_dimension_before_family_basis_quotient'] == 27
    assert report['common_positive_Higgs_norm_cancels_all_mass_ratios']
    assert not report['physical_mass_scale_or_ratios_follow_from_twelve_currents_alone']
    assert not report['full_mass_actions_of_controls_are_identical']
    first, second = report['exact_positive_hierarchy_controls']
    for control in (first, second):
        assert len(control['complete_K45_mass_squared_multiset_exact']) == 45
        assert control['complete_K45_mass_squared_multiset_exact'].count('0') == 3
        assert control['same_verified_twelve_current_Higgs_tensor_factorization']
        for row in control['channels']:
            values = list(map(s.Rational, row['positive_Yukawa_singular_values_exact']))
            assert 0 < values[0] < values[1] < values[2]
    first_ratios = {row['channel']: row['relative_to_lightest_exact'] for row in first['channels']}
    second_ratios = {row['channel']: row['relative_to_lightest_exact'] for row in second['channels']}
    assert first_ratios['charged_lepton'] == ['1', '4/3', '5/3']
    assert second_ratios['charged_lepton'] == ['1', '8/3', '5']
    assert first_ratios != second_ratios


def test_internal_multiplets_are_decoded_from_joint_native_operator_records():
    rows = m.decoded_multiplets()
    assert [(r['q6Y'], r['dimension']) for r in rows] == [('-4', 3), ('-3', 2), ('1', 6), ('2', 3), ('6', 1)]
    assert all('species' not in r and 'internal_type' not in r for r in rows)
    assert sum(r['dimension'] for r in rows) == 15
    assert {r['weak_C2'] for r in rows} == {'0', '3/4'}
    assert {r['color_C2'] for r in rows} == {'0', '4/3'}


def test_higgs_lines_solve_the_new_full_current_equations():
    representation = m.representation()
    for row in m.invariant_higgs_lines():
        cs = row['components']
        assert all(c.T == c for c in cs)
        for g, lifted in zip(representation['matter_currents'], representation['lifted'], strict=True):
            h = -lifted[3:, 3:].T if row['conjugate'] else lifted[3:, 3:]
            for k in range(2):
                assert m.zero(g.T*cs[k]+cs[k]*g+sum((h[j, k]*cs[j] for j in range(2)), s.zeros(15)))


def test_new_native_primitive_controls_match_full_six_dimensional_actions():
    r = m.representation()
    for kind, coefficients in m.primitive_coefficients().items():
        generator = m.clean(sum((a*b for a, b in zip(coefficients, r['native']['generators'], strict=True)), s.zeros(6)))
        pulse = m.native_pulse(kind, 0.125)
        assert np.max(abs(pulse-expm(np.asarray(generator, complex)*0.125))) < 1e-14
        d = 4 if kind == 'J' else 2
        assert np.max(abs(pulse[d:, :d]), initial=0) < 1e-14


def test_complete_source_current_columns_and_spin_cover_holonomy(built):
    report, arrays = built
    assert len(report['source_current_executions']) == 14
    for index, current in enumerate(m.canonical_currents()):
        target = expm(0.125*np.asarray(m.exterior_derivative(m.spin_lift(current)), complex))
        actual = arrays[f'current_{index}_source']
        assert actual.shape == (32, 32)
        assert np.linalg.norm(actual-target, 2) < m.ADMISSION
        assert np.linalg.norm(actual.conj().T@actual-np.eye(32), 2) < m.ADMISSION
    sign = np.diag([(-1)**((mask >> 3).bit_count()) for mask in range(32)])
    assert np.linalg.norm(arrays['weak_2pi_source']-sign, 2) < m.ADMISSION
    assert np.linalg.norm(arrays['weak_4pi_source']-np.eye(32), 2) < m.ADMISSION


def test_new_current_executor_does_not_use_the_old_fixture_execution_backend(monkeypatch):
    def unavailable(*args, **kwargs): raise AssertionError('old fixture executor used')
    monkeypatch.setattr(m.compiler, '_native', unavailable)
    monkeypatch.setattr(m.compiler, 'execute', unavailable)
    rows, arrays = m.current_execution()
    assert len(rows) == 14 and arrays['current_0_source'].shape == (32, 32)


def test_complete_native_pointer_law_including_excluded_masks(built):
    report, arrays = built
    kraus = arrays['occupation_instrument']
    assert kraus.shape == (32, 32, 32)
    assert np.linalg.norm(sum((k.conj().T@k for k in kraus), np.zeros((32, 32), complex))-np.eye(32), 2) < m.ADMISSION
    for mask, k in enumerate(kraus):
        ideal = np.zeros((32, 32)); ideal[mask, mask] = 1
        assert np.linalg.norm(k-ideal, 2) < m.ADMISSION
    counts = arrays['occupation_counts']
    assert counts.dtype == np.dtype('int64') and counts.shape == (32,) and counts.sum() == m.TRIALS
    records = report['record_identification']
    assert len(records['excluded_outcomes_retained']) == 17
    assert [r['decoded_dimension'] for r in records['decoded_internal_multiplets']] == [3, 2, 6, 3, 1]
    assert not records['species_name_table_used_by_decoder']


def test_all_proper_code_complement_exports_are_retained(built):
    for dimension in (2, 4):
        kraus = built[1][f'proper_code_{dimension}_complete_export']
        assert kraus.shape == (7-dimension, dimension+1, 6)
        assert np.array_equal(sum((k.conj().T@k for k in kraus), np.zeros((6, 6), complex)), np.eye(6))
        assert all(np.count_nonzero(k[:, :dimension]) == 0 for k in kraus[1:])


def test_count_inverse_does_not_read_source_amplitudes_or_mass_parameters(built, monkeypatch):
    def unavailable(*args, **kwargs): raise AssertionError('forward source or mass action used by inverse')
    monkeypatch.setattr(m, 'execute_word', unavailable)
    monkeypatch.setattr(m.existing_action, 'mass_kernel', unavailable)
    decoded = m.decode_record_counts(built[1]['occupation_counts'])
    assert decoded == built[0]['record_identification']['decoded_internal_multiplets']
    with pytest.raises(ValueError): m.decode_record_counts(built[1]['occupation_counts'].astype(float))
    with pytest.raises(ValueError): m.decode_record_counts(built[1]['occupation_counts'][:15])


def test_fixed_draw_sampler_is_stable_under_source_probability_roundoff():
    with np.load(m.ROOT/'model/common_source/native_matter_attachment_results.npz', allow_pickle=False) as stored:
        probabilities = stored['occupation_complete_probabilities']
    counts, diagnostics = m.sample_complete_pointer_law(probabilities)
    assert counts.tolist() == [0, 0, 0, 17297, 0, 17313, 17587, 0, 0, 17530, 17493, 0,
                               17471, 0, 0, 17610, 0, 17388, 17604, 0, 17373, 0, 0,
                               17477, 17560, 0, 0, 17403, 0, 17607, 17431, 0]
    selected = np.flatnonzero(probabilities > 0.01)
    for sign in (-1, 1):
        for affected in (selected[:1], selected[::2]):
            perturbed = probabilities.copy()
            perturbed[affected] += sign*16*np.spacing(perturbed[affected])
            fresh, _ = m.sample_complete_pointer_law(perturbed)
            assert np.array_equal(fresh, counts)
    assert diagnostics['minimum_draw_to_source_CDF_distance'] > m.SAMPLER_BOUNDARY_GUARD
    assert diagnostics['grid_total_variation_upper'] < m.ADMISSION
    assert diagnostics['boundary_guard_failure_action'] == 'abort_without_resampling'
    assert not diagnostics['statistical_coverage_conditioned_on_boundary_guard_pass']
    assert np.random.PCG64(m.SEED).random_raw(4).tolist() == [18208654110976716130,
        17750534445809024300, 624749274061255245, 8174374305464244112]


def test_fixed_draw_sampler_rejects_a_numerically_fragile_boundary():
    boundary = int(np.random.PCG64(m.SEED).random_raw(1)[0])/2**64
    probabilities = np.zeros(32); probabilities[:2] = (boundary, 1-boundary)
    with pytest.raises(ArithmeticError, match='too close'):
        m.sample_complete_pointer_law(probabilities)


def test_fixed_draw_sampler_keeps_zero_outcomes_and_rejects_incomplete_laws():
    probabilities = np.zeros(32); probabilities[17] = 1
    counts, _ = m.sample_complete_pointer_law(probabilities)
    assert counts.dtype == np.dtype('int64') and counts.shape == (32,)
    assert counts[17] == counts.sum() == m.TRIALS
    for invalid in (probabilities[:15], probabilities.astype(np.float32),
                    np.full(32, np.nan), -probabilities):
        with pytest.raises(ValueError): m.sample_complete_pointer_law(invalid)
    with pytest.raises(ArithmeticError, match='mass differs'):
        m.sample_complete_pointer_law(2*probabilities)


def test_global_center_quotient_matches_the_same_lift(built):
    report = built[0]['global_cover_contract']
    assert report['compatible_group_quotient'] == '(SU3 x SU2 x U1_q)/Z6'
    assert len(report['central_kernel_triples_color_weak_qsixthroot']) == 6
    assert not report['weak_spin_lift_descends_to_original_SO3_without_a_cover']


def test_lift_carries_classical_record_weights_once_not_by_exterior_degree(built):
    report = built[0]['lifted_repair_contract']
    weights = list(map(s.Rational, report['complete_native_joint_record_weights_exact']))
    assert len(weights) == 128 and sum(weights) == 1
    assert report['chosen_lift_completeness_on_all32_auxiliary_inputs']
    assert report['lifted_control_path_or_cover_element_is_required']
    assert not report['six_by_six_endpoint_alone_determines_spin_cover_holonomy']
    assert not report['raw_Kraus_exterior_lift_is_trace_preserving']
    assert not report['literal_128_branch_auxiliary_native_repair_history_executed_here']


def test_chosen_architecture_and_unproved_physical_claims_are_explicit(built):
    report = built[0]
    assert not report['architecture']['native_port_carrier_changed']
    assert report['architecture']['spin_half_lift_and_central_trace_redistribution_are_chosen']
    assert report['architecture']['nonvacuum_even_exterior_selection_is_chosen']
    assert not report['architecture']['physical_Weyl_fermion_statistics_inferred_from_auxiliary_even_parity']
    assert report['physical_field_attachment']['physical_CAR_is_a_separate_architecture_functor']
    assert not report['physical_field_attachment']['physical_spatial_Spin_is_internal_weak_SU2']
    assert report['empirical_mass_or_coupling_inputs'] == []
    assert not report['scope']['measured_particle_mass_spectrum_predicted']
    assert not report['scope']['quantum_gauge_Higgs_loops_executed']


def test_noncommuting_onebody_exact_compilation_is_rejected():
    generator = s.I*(m.unit(5, 0, 0)+m.unit(5, 0, 1)+m.unit(5, 1, 0))
    with pytest.raises(ValueError, match='commuting'):
        m.onebody_word(generator, 0.1)


@pytest.mark.parametrize('value', [s.eye(3), 'not an operator'])
def test_exterior_derivative_requires_the_declared_five_mode_domain(value):
    with pytest.raises(ValueError): m.exterior_derivative(value)


def _receipt(path, report, arrays):
    artifact = path.with_suffix('.npz'); np.savez_compressed(artifact, **arrays)
    value = deepcopy(report); value['source_bindings'] = m.source_bindings()
    value['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    path.write_text(json.dumps(value, allow_nan=False))


@pytest.mark.parametrize('tamper', ['missing', 'extra', 'changed', 'nonfinite', 'dtype', 'shape', 'count'])
def test_full_array_custody_even_when_binary_is_rehashed(tmp_path, monkeypatch, tamper):
    report, arrays = {'count': 15}, {'matrix': np.eye(2), 'counts': np.array([5, 10], dtype=np.int64)}
    monkeypatch.setattr(m, 'build', lambda: (report, arrays))
    damaged = {key: value.copy() for key, value in arrays.items()}
    if tamper == 'missing': del damaged['matrix']
    elif tamper == 'extra': damaged['extra'] = np.eye(2)
    elif tamper == 'changed': damaged['matrix'][0, 0] += 1e-5
    elif tamper == 'nonfinite': damaged['matrix'][0, 0] = np.nan
    elif tamper == 'dtype': damaged['matrix'] = damaged['matrix'].astype(np.float32)
    elif tamper == 'shape': damaged['matrix'] = damaged['matrix'].reshape(4)
    else: damaged['counts'][0] += 1
    path = tmp_path/'receipt.json'; _receipt(path, report, damaged)
    with pytest.raises(ValueError): m.check(path)


def test_receipt_types_source_drift_and_exact_pass(tmp_path, monkeypatch):
    report, arrays = {'count': 1}, {'matrix': np.eye(2)}
    monkeypatch.setattr(m, 'build', lambda: (report, arrays))
    path = tmp_path/'receipt.json'; _receipt(path, report, arrays); m.check(path)
    value = json.loads(path.read_text()); value['count'] = True; path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='type'): m.check(path)
    _receipt(path, report, arrays)
    pins = m.source_bindings(); sequence = iter([pins, {**pins, 'changed': '1'}])
    monkeypatch.setattr(m, 'source_bindings', lambda: next(sequence))
    with pytest.raises(ValueError, match='changed during replay'): m.check(path)
