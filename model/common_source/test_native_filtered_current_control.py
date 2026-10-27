"""Exact logical current transport and the state-reset countercontrol."""
from fractions import Fraction

import numpy as np
import pytest
import sympy as s

from . import native_filtered_current_control as n
from . import native_sm_current as declared


@pytest.fixture(scope='module')
def architecture():
    return n.exact_architecture()


def test_native_source_derivatives_and_oriented_rotations_precede_remapping(architecture):
    d = architecture
    assert n._zero(d['P']*d['P']-d['P']) and s.trace(d['P']) == 3
    assert n._zero(d['class_sum']*d['P']-(2+2*s.sqrt(5))*d['P'])
    assert n._rank(n._columns(d['raw_derivatives'])) == 6
    assert n._rank(n._columns(d['filtered_hamiltonians'])) == 6
    for p, (K, H, F) in enumerate(zip(d['cayley_generators'], d['filtered_hamiltonians'], d['generators'])):
        index = d['oriented_rotation_indices'][p]
        row = d['rotation_rows'][index]
        assert row['order'] == 5 and row['port_action'][p] == p
        U = d['rotation_actions'][index]
        assert n._zero((s.eye(6)+K)*(s.eye(6)-K).inv()-U)
        assert n._zero(H-d['P']*d['raw_derivatives'][p]*d['P'])
        assert n._zero(F-K-s.I*H) and n._zero(K*H-H*K)


def test_bidirectional_exact_control_remap_is_not_primitive_equality(architecture):
    d = architecture
    assert n._zero(d['declared_from_filtered']*d['filtered_from_declared']-s.eye(12))
    old = declared.exact_architecture()['generators']
    for p in range(12):
        field = tuple(int(p == j) for j in range(12))
        mapped = n.to_declared_controls(field)
        assert n._zero(n.generator(field)-n._combination(mapped, old))
        assert n.from_declared_controls(mapped) == field
        assert not n._zero(d['generators'][p]-old[p])
    assert d['native_to_fixture_port'] == declared.exact_architecture()['native_to_fixture_port']


def test_all_brackets_covariance_and_observable_center_are_certified():
    report = n.exact_certificate()
    assert report['logical_current_dimension'] == report['logical_EndM6_response_rank'] == 12
    assert len(report['commutator_closure']) == 66
    assert report['native_port_covariance_identities'] == 720
    assert report['all_primitive_admission_and_reverse_remapping_checks'] == 24
    assert not report['primitive_port_drives_equal_declared_primitives']
    assert report['remapping_can_require_irrational_controls']


@pytest.mark.parametrize('values', [
    [0]*12, (0,)*11, (True,)+(0,)*11, (0.1,)+(0,)*11,
    (s.Float(1),)+(0,)*11, (s.Symbol('x', real=True),)+(0,)*11,
    (s.I,)+(0,)*11, (s.oo,)+(0,)*11, (s.pi,)+(0,)*11,
])
def test_controls_reject_coercion_and_values_outside_exact_domain(values):
    with pytest.raises(ValueError): n.generator(values)
    with pytest.raises(ValueError): n.to_declared_controls(values)
    with pytest.raises(ValueError): n.from_declared_controls(values)


def test_exact_algebraic_remapping_does_not_bypass_rational_source_domain():
    values = (1,)+(0,)*11
    with pytest.raises(ValueError, match='outside the declared rational source API'):
        n.declared_rational_controls(values)
    zero = (0,)*12
    assert n.declared_rational_controls(zero) == (Fraction(0),)*12
    mixed = (Fraction(1, 3), s.sqrt(5))+(0,)*10
    assert n.from_declared_controls(n.to_declared_controls(mixed)) == mixed
    target = (1,)+(0,)*11
    filtered = n.from_declared_controls(target)
    admitted = n.declared_rational_controls(filtered)
    assert admitted == (Fraction(1),)+(Fraction(0),)*11
    assert n._zero(declared.public_generator(admitted)-n.generator(filtered))


def test_operator_filter_is_not_actual_state_reset_or_its_adjoint(architecture):
    d = architecture; P, Q = d['P'], d['Q']
    assert n.state_cooling(Q) == P
    assert n.state_cooling_adjoint(P) == s.eye(6)
    assert n._zero(P*Q*P)
    for i in range(6):
        for j in range(6):
            unit = s.eye(6)[:, i]*s.eye(6)[j, :]
            cooled = n.state_cooling(unit)
            assert n._zero(Q*cooled) and n._zero(cooled*Q)
            assert s.trace(cooled) == s.trace(unit)
            assert n.state_cooling(cooled) == cooled
    # The state reset is trace preserving. Hamiltonian filtering sends Q to
    # zero, whereas cooling sends it to P. They cannot be interchangeable.
    report = n.cooling_countercontrol()
    assert report['continuing_work_state_response_rank'] == 8
    assert report['central_phase_invisible_in_each_classical_cooling_record']
    assert not report['full_recorded_instrument_rank_claimed']


def test_all_actual_classical_cooling_records_erase_relative_phase(architecture):
    from . import native_family_selection as family
    from . import native_code_transport as transport
    maps, labels = family.cooling_maps()
    center = np.zeros((64, 64), complex)
    T = np.asarray(architecture['geometric_frame'], complex)
    addresses = transport.product_addresses()
    center[np.ix_(addresses, addresses)] = 1j*T@np.asarray(architecture['P'], complex)@T.conj().T
    for K, label in zip(maps, labels, strict=True):
        # K exp(t center) differs from K only by a branchwise scalar phase.
        expected = 1j*K if label[0] == 0 else np.zeros_like(K)
        np.testing.assert_allclose(K@center, expected, atol=3e-14, rtol=0)
    assert len(labels) == 605


@pytest.mark.parametrize('bad', [np.eye(6), s.eye(5), s.eye(6)*s.Float(1), s.eye(6)*s.oo])
def test_state_countercontrol_rejects_other_operator_domains(bad):
    with pytest.raises(ValueError): n.state_cooling(bad)
    with pytest.raises(ValueError): n.state_cooling_adjoint(bad)


def test_conditional_spin_exterior_and_higgs_interfaces():
    report = n.theorem_interfaces()
    assert report['conditional_spin_lift_remapping_checks'] == 12
    assert report['conditional_exterior_lift_remapping_checks'] == 12
    assert report['finite_admissible_component_masks'] == [227, 796]
    assert report['selected_matter_dimension'] == 15
    assert report['one_doublet_invariant_lines'] == 3
    assert report['Higgs_tensor_substitution_checks'] == 72
    assert report['Higgs_kernel_dimensions_preserved_by_invertible_current_remapping']
    assert not report['physical_spin_statistics_or_family_count_derived']
    assert not report['Yukawa_coefficients_or_particle_masses_derived']


def test_report_never_claims_source_compilation_or_full_model_closure():
    report = n.run(include_theorem_interfaces=False)
    scope = report['scope']
    assert scope['Hamiltonian_operator_filter_is_an_added_control_choice']
    assert scope['complete_only_over_declared_logical_current_tangent']
    for name in ('source_compilation_executed', 'source_pulse_or_record_history_certificate_provided',
                 'unchanged_C_E_selects_this_current', 'all_physically_accessible_interventions_classified',
                 'full_model_A1_complete_accessible_response_discharged',
                 'physical_gauge_redundancy_or_SM_spectrum_derived',
                 'finite_noncommuting_pulse_products_identified_by_linear_remapping'):
        assert scope[name] is False
