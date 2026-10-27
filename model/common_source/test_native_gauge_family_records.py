"""Complete current records, invariant Higgs domains and ungated controls."""
import json
import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import native_gauge_family_records as m


@pytest.fixture(scope='module')
def instrument():
    return m.instrument()


def test_complete_current_instrument_preserves_trace_on_arbitrary_entangled_inputs(instrument):
    kraus = instrument['kraus']
    assert kraus.shape == (13, 15, 15)
    rng = np.random.default_rng(20261027)
    psi = rng.normal(size=(15, 4)) + 1j*rng.normal(size=(15, 4))
    psi /= np.linalg.norm(psi)
    outputs = np.stack([k @ psi for k in kraus])
    assert abs(np.vdot(outputs, outputs) - 1) < 1e-13
    assert instrument['completeness_error'] < 2e-15
    data = m.current_certificate()
    assert m._zero(sum((g.H*g for g in data['matter_frame']), s.zeros(15)) - data['casimir15'])


def test_all_adjoint_outcomes_are_needed_for_finite_gauge_covariance(instrument):
    errors = m.covariance_errors(instrument['kraus'])
    assert errors.shape == (12,) and np.max(errors) < 2e-14
    data = m.current_certificate()
    V = instrument['kraus'].reshape(195, 15)
    U = expm(float(m.PROTOCOL_ANGLE) * np.asarray(data['representation']['matter_currents'][0], complex))
    assert np.linalg.norm(np.kron(np.eye(13), U) @ V - V @ U, 2) > .1


def test_coarse_binary_channel_is_covariant_without_postselecting_jump(instrument):
    data = m.current_certificate()
    U = expm(float(m.PROTOCOL_ANGLE) * np.asarray(data['representation']['matter_currents'][3], complex))
    rng = np.random.default_rng(33)
    psi = rng.normal(size=15) + 1j*rng.normal(size=15)
    psi /= np.linalg.norm(psi)
    state = np.outer(psi, psi.conj())
    transformed = U @ state @ U.conj().T
    kraus = instrument['kraus']
    for group in (kraus[:1], kraus[1:]):
        first = sum((k @ transformed @ k.conj().T for k in group), np.zeros((15, 15), complex))
        second = U @ sum((k @ state @ k.conj().T for k in group), np.zeros((15, 15), complex)) @ U.conj().T
        np.testing.assert_allclose(first, second, atol=2e-15, rtol=0)


def test_higgs_supports_use_both_components_and_are_invariant_under_all_currents():
    data = m.current_certificate()
    supports = m.higgs_supports()
    assert len(supports) == 3
    assert sorted((r['left_dimension'], r['right_dimension']) for r in supports) == [(2, 1), (6, 3), (6, 3)]
    for row in supports:
        left, right = row['left_projector'], row['right_projector']
        assert left * right == s.zeros(15)
        for component in row['components']:
            assert m._zero((left + right) * component - component)
            assert m._zero(left * component * left) and m._zero(right * component * right)
        for g in data['representation']['matter_currents']:
            assert m._zero(left*g - g*left) and m._zero(right*g - g*right)
        # A neutral-component support alone fails the full weak covariance.
        neutral_rows = {i for (i, _), value in row['components'][1].todok().items() if value != 0}
        neutral = s.diag(*(int(i in neutral_rows) for i in range(15)))
        assert any(not m._zero(neutral*g-g*neutral) for g in data['representation']['matter_currents'])


def test_source_derived_scalar_effects_determine_gates_for_all_input_states(instrument):
    kraus = instrument['kraus']
    jump_effect = sum((k.conj().T @ k for k in kraus[1:]), np.zeros((15, 15), complex))
    p_native, q_native = m.native_event_rates()
    assert (p_native, q_native) == (s.Rational(1, 12), s.Rational(1, 6))
    casimirs = set()
    for row in m.sector_gates():
        casimirs.add((row['left_casimir'], row['right_casimir']))
        for leg in ('left', 'right'):
            P = np.asarray(row[leg + '_projector'], float)
            probability = float(row[leg + '_jump'])
            np.testing.assert_allclose(jump_effect @ P, probability * P, atol=1e-15, rtol=0)
        # This operator identity gives independence for arbitrary entangled L/R states.
        joint = np.kron(jump_effect, jump_effect)
        domain = np.kron(np.asarray(row['left_projector'], float), np.asarray(row['right_projector'], float))
        np.testing.assert_allclose(joint @ domain, float(row['left_jump']*row['right_jump']) * domain,
                                   atol=2e-15, rtol=0)
    assert casimirs == {(s.Rational(25, 8), s.Integer(4)),
                        (s.Rational(25, 8), s.Integer(3)),
                        (s.Rational(9, 8), s.Integer(3))}


def test_all_raw_native_and_binary_records_are_retained_and_factorize():
    for row in m.sector_gates():
        records = m.complete_binary_records(row['origin'])
        assert len(records) == 512
        assert len({r[:6] for r in records}) == 512
        assert any(r[-1] == 0 for r in records)
        assert s.simplify(sum(r[-1] for r in records) - 1) == 0
        p = sum(r[-1] for r in records if r[-3])
        q = sum(r[-1] for r in records if r[-2])
        joint = sum(r[-1] for r in records if r[-3] and r[-2])
        assert s.simplify(p - row['p']) == 0
        assert s.simplify(q - row['q']) == 0
        assert s.simplify(joint - row['p']*row['q']) == 0


def test_zero_angle_is_complete_idle_instrument_not_a_conditioned_channel():
    idle = m.instrument(0)['kraus']
    np.testing.assert_array_equal(idle[0], np.eye(15))
    np.testing.assert_array_equal(idle[1:], np.zeros((12, 15, 15)))
    assert all(r['p'] == r['q'] == 0 for r in m.sector_gates(0))


def test_declared_geometry_catalog_retains_every_nonidentity_order():
    catalog = m.protocol_catalog()
    assert tuple(catalog.values()) == (s.pi/5, s.pi/3, s.pi/2)
    for angle in catalog.values():
        assert m.instrument(angle)['completeness_error'] < 3e-15


def test_catalog_retains_exact_zero_family_directions_and_serializes_without_infinities():
    catalog = m.catalog_summaries()
    assert len(catalog) == 3 and all(len(row['Higgs_lines']) == 3 for row in catalog)
    zeros = [line for protocol in catalog for line in protocol['Higgs_lines']
             if line['family_rank_exact'] < 3]
    assert len(zeros) == 1
    assert zeros[0]['q'] == 0
    assert zeros[0]['family_eigenvalues'][0] == 0
    assert zeros[0]['family_eigenvalue_ratios'] is None
    json.dumps(catalog, allow_nan=False)


def test_named_field_handoff_uses_actual_higgs_intertwiner_without_changing_global_policy():
    channels = m.candidate_channels()
    assert set(channels) == {'up', 'down', 'charged_lepton'}
    assert len({(row['reset'], row['conjugate'], row['angle']) for row in channels.values()}) == 1
    for row in channels.values():
        assert row['Y'].shape == (3, 3)
        assert row['Higgs_components'].shape == (2, 15, 15)
        assert np.linalg.eigvalsh(row['Y']).min() > 0
        assert abs(np.trace(row['Y']) - 1) < 1e-14


def test_actual_higgs_tensors_have_unit_particle_modulus_and_sign_only_dictionary():
    rows = m.higgs_tensor_certificate()
    assert {row['action_channel']: row['neutral_rank'] for row in rows} == {
        'up': 6, 'down': 6, 'charged_lepton': 2}
    assert all(row['directed_nonzero_singular_value_exact'] == '1' for row in rows)
    assert all(row['dictionary_scale_is_sign_only'] for row in rows)
    assert {row['action_channel']: row['dictionary_scale_exact'] for row in rows} == {
        'up': '-1', 'down': '1', 'charged_lepton': '1'}


def test_unconditional_duty_traces_and_event_memory_correlation_are_exact():
    duties = m.duty_coefficients()
    assert tuple(duties) == ('always', 'C', 'E', 'CE')
    for index, model in enumerate(m.family_models()):
        p, q = model['gate']['p'], model['gate']['q']
        Y, P = model['family']['Y'], model['family']['projectors'][2]
        expected = {'always': 1, 'C': p, 'E': q, 'CE': p*q}
        for policy, trace in expected.items():
            row = duties[policy][index]
            assert row['raw_record_count'] == 512 and row['old_memory_input_count'] == 12
            assert s.cancel(row['trace']-trace) == 0
            assert s.cancel(s.trace(row['coefficient'])-trace) == 0
        correction = duties['C'][index]['coefficient'] - Y + (1-p)*P
        assert all(s.cancel(value) == 0 for value in correction)
        assert np.linalg.norm(np.asarray(duties['C'][index]['coefficient']-p*Y, float)) > .01


def test_event_only_duties_retain_exact_rank_two_geometric_support():
    duties = m.duty_coefficients()
    for policy in ('E', 'CE'):
        for row in duties[policy]:
            assert row['supported_memory_directions'] == (2, 10)
            assert row['family_rank_exact'] == 2
            assert row['family_eigenvalues'][0] == 0
            assert row['family_eigenvalue_ratios'] is None
    for policy in ('always', 'C'):
        assert all(row['family_rank_exact'] == 3 for row in duties[policy])
    catalog = m.duty_catalog()
    assert len(catalog) == 4
    assert all(row['actual_heaviest_span_diagnostic'] <= row['global_heaviest_ratio_bound_diagnostic']
               for row in catalog)
    assert all(not row['pulse_duty_renormalized'] for row in catalog)
    json.dumps(catalog, allow_nan=False)


@pytest.mark.parametrize('change', ['missing', 'nonfinite', 'missing_idle', 'missing_jump'])
def test_missing_current_outcomes_and_nonfinite_inputs_fail_closed(instrument, change):
    kraus = instrument['kraus'].copy()
    if change == 'missing':
        kraus = kraus[:-1]
    elif change == 'nonfinite':
        kraus[0, 0, 0] = np.nan
    elif change == 'missing_idle':
        kraus[0] = 0
    else:
        kraus[1] = 0
    with pytest.raises(ValueError):
        m.validate_kraus(kraus)


@pytest.mark.parametrize('angle', [s.oo, s.I, s.Symbol('t', real=True)])
def test_undefined_or_complex_protocol_angles_are_rejected(angle):
    with pytest.raises(ValueError):
        m.instrument(angle)
