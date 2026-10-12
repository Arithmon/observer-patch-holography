"""Independent weighted-form and actual-state checks for internal sanity tests."""
import json
import math

import numpy as np
import pytest

from ..observer.scalar_response.field import ScalarField
from ..observer.scalar_response.local_evolution import LocalEvolution
from . import sanity as s


@pytest.fixture(scope='module')
def dispersion():
    return s.scalar_mode_dispersion()


@pytest.fixture(scope='module')
def clocks():
    return [s.clock_flow_comparison(d, steps=16) for d in (2, 4)]


@pytest.mark.parametrize('q', s.GRIDS)
def test_modes_come_from_independent_edge_energy_assembly(dispersion, q):
    field = ScalarField(q)
    incidence = np.zeros((q, q-1))
    for edge in range(q):
        if edge > 0:
            incidence[edge, edge-1] = -1
        if edge < q-1:
            incidence[edge, edge] = 1
    stiffness = incidence.T@np.diag(1/field.gaps)@incidence
    operator = stiffness/np.sqrt(np.outer(field.weights, field.weights))
    eigenvalues = np.linalg.eigvalsh(operator)
    level = next(row for row in dispersion['levels'] if row['q'] == q)
    for mode in level['modes']:
        expected = math.sqrt(1+sum(eigenvalues[n-1] for n in mode['standing_wave_label']))
        assert mode['actual_oscillator_frequency'] == pytest.approx(expected, abs=2e-11)
    assert level['nonuniform_grid']


def test_same_tensor_operator_is_used_by_actual_local_evolution(dispersion):
    local = LocalEvolution(13)
    mode = local.field.eigenvectors[:, 0]
    tensor = np.einsum('i,j,k->ijk', mode, mode, mode)
    acted = local.apply_action(tensor)
    omega = dispersion['levels'][0]['modes'][0]['actual_oscillator_frequency']
    np.testing.assert_allclose(acted, omega**2*tensor, atol=3e-12, rtol=0)


def test_resolved_modes_approach_massive_KG_and_do_not_assert_exact_linear_dispersion(dispersion):
    for index, label in enumerate(s.MODE_LABELS):
        rows = [level['modes'][index] for level in dispersion['levels']]
        errors = [row['relative_frequency_error'] for row in rows]
        assert errors[2] < errors[1] < errors[0]
        expected = math.sqrt(1+math.pi**2*sum(n*n for n in label))
        for row in rows:
            assert row['declared_continuum_KG_frequency'] == expected
            assert row['continuum_frequency_over_wavenumber'] > 1
            assert row['finite_action_identity_residual'] < 2e-12
    assert not dispersion['scope']['translation_momentum_operator_or_detector_constructed']
    assert not dispersion['scope']['independent_detector_frequency_measurement']
    assert not dispersion['scope']['Lorentz_invariance_derived']
    assert not dispersion['scope']['common_interacting_source_spatial_extension']


def test_actual_local_state_changes_but_its_modular_reference_is_stationary(clocks):
    for row in clocks:
        encoded = np.asarray(row['norm_one_witness'])
        witness = encoded[:, :, 0]+1j*encoded[:, :, 1]
        np.testing.assert_allclose(witness, witness.conj().T, atol=2e-15)
        assert np.linalg.norm(witness, 2) == pytest.approx(1, abs=2e-12)
        assert row['model_time_witness_expectation_derivative'] > 1e-3
        assert row['modular_state_derivative_frobenius'] < 2e-12
        assert row['central_difference_derivative_error'] < 2e-8
        assert not row['same_state_model_flow_equals_scalar_rescaled_modular_flow']
        assert not row['scope']['every_possible_observer_clock_rejected']
        assert not row['scope']['programmed_modular_clock_invalidated']


def test_report_pass_means_internal_checks_not_empirical_success(monkeypatch, dispersion, clocks):
    monkeypatch.setattr(s, 'scalar_mode_dispersion', lambda: dispersion)
    monkeypatch.setattr(s, 'clock_flow_comparison', lambda d, steps: clocks[(d//2)-1])
    report = s.report(16)
    json.dumps(report, allow_nan=False)
    assert report['all_pass']
    assert not report['scientific_verdict']['agreement_with_nature_established']
    assert report['scientific_verdict']['specified_model_time_modular_time_identification_rejected_on_tested_state']


@pytest.mark.parametrize('d,steps', [(True, 16), (8, 16), (2, True), (2, 0), (4, 1.5)])
def test_rejects_undeclared_execution_controls(d, steps):
    with pytest.raises(ValueError):
        s.clock_flow_comparison(d, steps)
