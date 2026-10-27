"""Unknown-input custody, charge, inference and held-out source-record tests."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import inspect
import json
import math
from unittest.mock import patch

import numpy as np
import pytest
from scipy.linalg import expm

from . import native_charged_spectroscopy as h


@pytest.fixture(scope='module')
def execution():
    return h.execute()


def test_pulse_schedule_never_receives_a_spectral_pole():
    assert tuple(inspect.signature(h.cycle_word).parameters) == ('mass', 'exchange', 'step')
    assert not inspect.signature(h.preparation_word).parameters
    assert tuple(inspect.signature(h.detector_basis_word).parameters) == ('axis',)
    assert tuple(inspect.signature(h.infer_poles).parameters) == ('contrasts', 'sampling_interval', 'contrast_error')
    assert set(h.preparation_target().nonzero()[0]) == {2, 4}


def test_native_pulse_cycle_retains_interaction_on_all_fock_columns():
    j, masses = h.spectrum.native_parameters()
    for m in masses:
        for dt in h.STEPS:
            actual = h.compiler.execute(h.cycle_word(m, j, dt), 3)
            np.testing.assert_allclose(actual, h.ideal_cycle(m, j, dt), atol=1e-11, rtol=0)
            # Filled dimer energy is zero; a free one-particle lift fails this.
            np.testing.assert_allclose(actual[:, 3], np.eye(8)[:, 3], atol=1e-11, rtol=0)


def test_transfer_reads_preserve_charge_and_resolve_complex_correlations():
    n = np.diag([x.bit_count() for x in range(8)])
    j, masses = h.spectrum.native_parameters()
    seed = h.frame_target()@h.preparation_target()
    for axis in ('X', 'Y'):
        q = h.transfer_operator(axis)
        np.testing.assert_array_equal(q@n, n@q)
        for m in masses:
            u = expm(-1j*float(h.SPACING)*h.hamiltonian(m, j))
            actual = np.vdot(u@seed, q@(u@seed)).real
            expected = u[1, 1].real if axis == 'X' else u[1, 1].imag
            assert abs(actual-expected) < 1e-14


def test_strang_remainder_contains_independent_matrix_exponential():
    j, masses = h.spectrum.native_parameters()
    for m in masses:
        for dt in h.STEPS:
            proof = h.continuum_bound(m, j, dt)
            discrepancy = np.linalg.norm(h.ideal_cycle(m, j, dt)-expm(-1j*float(dt)*h.hamiltonian(m, j)), 2)
            assert discrepancy <= float(F(proof['single_cycle_operator_error_bound']))+3e-15


def test_analytic_cycle_matches_independent_onsite_exchange_exponentials():
    j, masses = h.spectrum.native_parameters()
    for m in masses:
        for dt in h.STEPS:
            onsite = h.hamiltonian(m, F(0))
            half = expm(-.5j*float(dt)*onsite)
            expected = half@expm(-1j*float(dt)*h.hamiltonian(F(0), j))@half
            np.testing.assert_allclose(h.ideal_cycle(m, j, dt), expected, atol=3e-16, rtol=0)


def test_compensated_cycle_power_preserves_accumulated_nonunitary_error(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('a backend-dependent matrix power reached replay')
    monkeypatch.setattr(np.linalg, 'matrix_power', forbidden)
    # This exact dyadic shear has a large, known nonunitary power. A polar
    # projection or normalization would destroy the measured source error.
    shear = np.array([[1, 1j/8], [0, 1]], complex)
    expected = np.array([[1, 100003j/8], [0, 1]], complex)
    np.testing.assert_array_equal(h._compensated_power(shear, 100003), expected)
    np.testing.assert_array_equal(h._compensated_power(shear, 0), np.eye(2))


def test_complete_native_records_inputs_and_numerical_source_limits(execution):
    report, arrays = execution
    r = report['execution']
    assert r['terminal_records_per_protocol'] == 158
    assert r['successful_prefix_records'] == 60
    assert r['halt_records'] == 38
    assert arrays['prefix_decoded'].shape == (60, 512, 7)
    assert arrays['prefix_halts'].shape == (38, 512, 7)
    assert r['maximum_complete_dilation_error_diagnostic'] < h.SOURCE_ADMISSION
    assert r['maximum_prepared_source_probability_defect_diagnostic'] < h.SOURCE_ADMISSION
    assert 0 < r['prepared_seed_retained_halt_probability_diagnostic'] < h.SOURCE_ADMISSION
    assert not report['scope']['measured_particle_spectrum_or_physical_SM_mass_attachment_derived']


def test_blind_record_poles_contain_continuous_charged_response(execution):
    report, _ = execution
    j, masses = h.spectrum.native_parameters()
    for family, m in enumerate(masses):
        exact_gap = float(m*m/(np.sqrt(float(j*j+m*m))+float(j)))
        for refinement, dt in enumerate(h.STEPS):
            inferred = report['execution']['record_only_pole_inference'][family][refinement]
            finite_error = float(F(report['continuum'][family][refinement]['continuous_energy_error_bound']))
            assert abs(inferred['charged_gap']-exact_gap) < inferred['record_only_energy_error_bound']+finite_error
            assert inferred['charged_gap'] > inferred['record_only_energy_error_bound']
            assert inferred['held_out_maximum_complex_residual'] < 1e-10
            assert inferred['zero_gap_control_held_out_maximum_residual'] > 1e-6
            assert inferred['no_source_J_mass_or_target_pole_supplied_to_inference']


def test_inference_is_independent_of_source_parameters(execution, monkeypatch):
    records = execution[0]['execution']['complex_contrasts'][2][2]
    values = np.array([complex(*pair) for pair in records])
    def forbidden():
        raise AssertionError('source parameters reached the observer estimator')
    monkeypatch.setattr(h.spectrum, 'native_parameters', forbidden)
    inferred = h.infer_poles(values, sampling_interval=8.)
    assert 4.3e-7 < inferred['charged_gap'] < 4.4e-7


def test_held_out_records_are_not_used_for_training_and_reject_wrong_law(execution):
    values = np.array([complex(*pair) for pair in execution[0]['execution']['complex_contrasts'][0][0]])
    values[5] += .001j
    with pytest.raises(ValueError, match='held-out'):
        h.infer_poles(values, sampling_interval=8.)


def test_two_positive_energies_are_not_a_positive_charged_gap():
    energies = np.array([.05, .2])
    records = np.array([.4*np.exp(-8j*k*energies[0])+.6*np.exp(-8j*k*energies[1]) for k in range(6)])
    with pytest.raises(ValueError, match='negative and one positive'):
        h.infer_poles(records, sampling_interval=8.)


@pytest.mark.parametrize('bad', ['masked', 'short', 'nonfinite', 'singular', 'uncertain'])
def test_record_inference_refuses_invalid_or_unidentifiable_inputs(execution, bad):
    values = np.array([complex(*pair) for pair in execution[0]['execution']['complex_contrasts'][0][0]])
    error = 0.
    if bad == 'masked':
        values = np.ma.array(values)
    elif bad == 'short':
        values = values[:-1]
    elif bad == 'nonfinite':
        values[0] = np.nan
    elif bad == 'singular':
        values[:] = 1.
    else:
        error = 1.
    with pytest.raises(ValueError):
        h.infer_poles(values, sampling_interval=8., contrast_error=error)


def test_source_record_law_keeps_every_unnormalized_abort(execution):
    report, arrays = execution
    changed = deepcopy(arrays)
    index = np.argmax(np.linalg.norm(changed['prefix_halts'].reshape(38, -1), axis=1))
    changed['prefix_halts'][index] /= np.linalg.norm(changed['prefix_halts'][index])
    with pytest.raises(ValueError, match='input/reference'):
        h.metrics(changed)
    assert not report['resources']['finite_shot_campaign_executed']


def test_full_input_phase_tampering_is_not_hidden_by_probabilities(execution):
    _, arrays = execution
    changed = deepcopy(arrays)
    changed['prefix_coefficients'][:, :, 0] *= 1j
    changed['prefix_decoded'][:, :, 0] *= 1j
    with pytest.raises(ValueError, match='input/reference'):
        h.metrics(changed)


@pytest.mark.parametrize('bad', ['record', 'padding', 'macro', 'masked'])
def test_binary_schema_and_operator_checks_cover_unprepared_rows(execution, bad):
    _, arrays = execution
    changed = deepcopy(arrays)
    if bad == 'record':
        changed['success_labels'][0, 0] = 2
    elif bad == 'padding':
        changed['prefix_decoded'][0, 511, 0] += .001
    elif bad == 'macro':
        changed['cycles'][0, 0, 7, 7] += .001
    else:
        changed['prefix_halts'] = np.ma.array(changed['prefix_halts'])
    with pytest.raises(ValueError):
        h.metrics(changed)


def test_factor_basis_gauge_is_not_a_new_physical_protocol(execution):
    report, arrays = execution
    changed = deepcopy(arrays)
    gauge = np.diag([1j, -1., -1j, 1., 1j, 1.])[:, [2, 0, 3, 1, 5, 4]]
    changed['prefix_basis'] = changed['prefix_basis']@gauge
    changed['prefix_coefficients'] = np.einsum('ij,njc->nic', gauge.conj().T, changed['prefix_coefficients'])
    own = deepcopy(report)
    own['execution'] = h.metrics(changed)
    h.compare_execution(own, changed, report, arrays)


def test_reported_diagnostics_are_bound_to_original_arrays(execution):
    report, arrays = execution
    altered = deepcopy(report)
    altered['execution']['maximum_complete_dilation_error_diagnostic'] = 0.
    with pytest.raises(ValueError, match='own binary'):
        h.compare_execution(altered, arrays, report, arrays)


def test_gram_diagnostic_preserves_small_terms_under_row_reordering():
    columns = np.asarray([[1., 1.], [1e-8, 1e-8], [1e-8, -1e-8],
                          [1e-8, 1e-8], [1e-8, -1e-8]], complex)
    expected = np.asarray([[math.fsum([1., 1e-16, 1e-16, 1e-16, 1e-16]), 1.],
                           [1., math.fsum([1., 1e-16, 1e-16, 1e-16, 1e-16])]], complex)
    assert np.array_equal(h._compensated_gram(columns), expected)
    assert np.array_equal(h._compensated_gram(columns[::-1]), expected)


def test_factored_gram_matches_complete_complex_record_instrument():
    rng = np.random.default_rng(414)
    basis = rng.normal(size=(73, 6))+1j*rng.normal(size=(73, 6))
    basis = np.linalg.qr(basis)[0]
    coefficients = (rng.normal(size=(11, 6, 7))+1j*rng.normal(size=(11, 6, 7)))/20
    complete = np.einsum('ar,nrc->nac', basis, coefficients, optimize=False)
    expected = h._compensated_gram(complete.reshape(-1, 7))
    gram = h._compensated_gram(coefficients.reshape(11, 42)).reshape(6, 7, 6, 7)
    actual = np.einsum('rs,rcsd->cd', h._compensated_gram(basis), gram, optimize=False)
    np.testing.assert_allclose(actual, expected, atol=8*np.finfo(float).eps, rtol=0)


def test_service_counts_include_every_reused_source_cycle(execution):
    r = execution[0]['resources']
    components = (r['preparation_native_pulses']+r['occupation_frame_native_pulses']
        +r['encoding_decoding_native_pulses']+sum(r['C_E_repair_recovery_native_pulses_upper'].values())
        +r['evolution_native_pulses_upper']+r['detector_native_pulses'])
    assert r['one_successful_protocol_native_pulses_upper'] == components
    assert r['required_proper_load_export_service_events_upper'] == 2*components
    assert not r['physical_load_export_protection_and_memory_services_executed']


def test_SM_spectator_boundary_is_explicit(execution):
    certificate = execution[0]['SM_spectator_interface']
    assert certificate['actual_imported_source_generators'] == 12
    assert certificate['all_retained_Kraus_maps_allow_this_algebraic_extension']
    assert not certificate['ninety_state_SM_source_execution_performed']
    assert not certificate['SM_charge_Higgs_Yukawa_or_species_attachment_derived_by_this_extension']
    assert not execution[0]['charge_interface']['occupation_charge_identified_with_SM_electric_charge']


def test_receipt_roundtrip_replays_actual_arrays(execution, tmp_path):
    report, arrays = execution
    target = tmp_path/'receipt.json'
    artifact = target.with_suffix('.npz')
    np.savez_compressed(artifact, **arrays)
    saved = deepcopy(report)
    saved['source_bindings'] = h.source_bindings()
    saved['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    target.write_text(json.dumps(saved))
    with patch.object(h, 'execute', return_value=execution):
        h.check(target)


@pytest.mark.parametrize('kind', ['missing', 'extra', 'nonfinite'])
def test_publication_fresh_array_schema_is_checked(execution, kind):
    report, arrays = execution
    fresh = deepcopy(arrays)
    if kind == 'missing':
        fresh.pop('prefix_decoded')
    elif kind == 'extra':
        fresh['undeclared_helper'] = np.zeros(1)
    else:
        fresh['prefix_decoded'][0, 0, 0] = np.nan
    with pytest.raises(ValueError):
        h.compare_execution(report, arrays, report, fresh)


def test_publication_replay_requires_unchanged_source_pins(execution, tmp_path):
    report, arrays = execution
    binary = tmp_path/'source-pins.npz'
    path = binary.with_suffix('.json')
    np.savez_compressed(binary, **arrays)
    pins = h.source_bindings()
    changed = {**pins, 'altered_source': '0'*64}
    saved = {**report, 'source_bindings': pins,
        'artifact': {'filename': binary.name, 'sha256': hashlib.sha256(binary.read_bytes()).hexdigest()}}
    path.write_text(json.dumps(saved))
    with patch.object(h, 'execute', return_value=execution), patch.object(h, 'source_bindings', side_effect=[pins, changed]):
        with pytest.raises(ValueError, match='changed during replay'):
            h.check(path)
