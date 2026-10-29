from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import action
from ..local_repair_family import local_bias as local, modular_clock


def counts():
    source = json.loads((action.ROOT/'model/local_repair_family/record_clock_results.json').read_text())
    return {tuple(b['seam']): (b['trials'], b['coin1_count']) for b in source['acquisition']['blocks']}


@pytest.mark.parametrize('parameters', [modular_clock.parameters(), local.Parameters(),
    local.Parameters(port_weights=(1, 2, 3, 4)*3, seam_weights=tuple(range(1, 16))*2,
                     c_weight=F(0), repair_probability=F(1, 4))])
def test_native_parent_and_ground_state(parameters):
    exact = action.native_binding(parameters)
    assert exact['nonzero_undirected_hops'] == 300
    weights = tuple(map(F, parameters.port_weights))
    parent = action.kinetic_parent(weights, seam_weights=parameters.seam_weights,
                                   repair_probability=parameters.repair_probability)
    # Independent matrix from source transitions, with no local-projector builder.
    kernel = np.array([[float(x) for x in row] for row in local.kernel(parameters)])
    pi = np.array(list(map(float, local.reference(parameters))))
    assert np.max(abs(parent-(np.eye(66)-kernel*np.sqrt(pi[:, None]/pi[None, :])))) < 2e-15
    assert np.linalg.norm(parent@np.sqrt(pi)) < 2e-15
    eigenvalues = np.linalg.eigvalsh(parent)
    assert abs(eigenvalues[0]) < 1e-14 and eigenvalues[1] > 0
    assert eigenvalues[-1] <= float(parameters.repair_probability)+1e-14


def test_local_projectors_and_interaction_are_actual_two_particle_blocks():
    p = modular_clock.parameters()
    for q in action.local_projectors(tuple(map(F, p.port_weights))):
        assert np.linalg.norm(q@q-q) < 2e-15
        assert abs(np.trace(q)-10) < 2e-15
    diagnostics = action.parent_diagnostics(p)
    assert len(diagnostics['exact_density_interactions']) == 30
    assert {r['density_interaction'] for r in diagnostics['exact_density_interactions']} == {F(-1, 30)}
    assert diagnostics['same_reference_gibbs_identification_excluded']
    # Uniform reference commutes, but its modular flow is scalar and cannot
    # produce the nonzero parent kinetics.
    uniform = action.parent_diagnostics(local.Parameters())
    assert not uniform['same_reference_gibbs_identification_excluded']
    assert uniform['affine_modular_identification_excluded']


def test_readonly_inference_parent_and_complete_output_bound():
    inferred = action.infer_parent(counts())
    assert inferred['status'] == 'accept'
    actual = action.kinetic_parent(tuple(map(F, modular_clock.parameters().port_weights)))
    learned = action.kinetic_parent(inferred['relative_weights'])
    bound = float(inferred['parent_operator_error_upper'])
    assert np.linalg.norm(actual-learned, 2) < bound < .03
    # Weyl's bound checked independently across the whole spectrum.
    assert max(abs(np.linalg.eigvalsh(actual)-np.linalg.eigvalsh(learned))) < bound
    # Entangled Choi input under both unitaries: complete quantum comparison,
    # without successful-output conditioning or a chosen detector branch.
    def unitary(h):
        values, vectors = np.linalg.eigh(h)
        return (vectors*np.exp(-1j*values))@vectors.conj().T
    u, v = unitary(actual), unitary(learned)
    overlap = np.trace(u.conj().T@v)/66
    choi_halftrace = np.sqrt(max(0, 1-abs(overlap)**2))
    assert choi_halftrace < bound


def test_projector_derivative_bound_checks_interval_corners():
    def q(p):
        return np.array([[1-p, -np.sqrt(p*(1-p))],[-np.sqrt(p*(1-p)), p]])
    for p in (.125, .2, .5, .8, .875):
        for r in (.125, .2, .5, .8, .875):
            assert np.linalg.norm(q(p)-q(r), 2) <= 2*abs(p-r)+1e-15


def test_stop_and_malformed_data_are_not_success():
    data = counts()
    data[next(iter(data))] = (2**20, 0)
    assert action.infer_parent(data)['status'] == 'stop'
    with pytest.raises(ValueError):
        action.infer_parent({})
    with pytest.raises(ValueError):
        action.kinetic_parent((1.0,)*12)
    with pytest.raises(ValueError):
        action.kinetic_parent((F(1),)*12, repair_probability=F(-1, 10**400))


def test_native_current_mutation_cannot_keep_energy_certificate(monkeypatch):
    monkeypatch.setattr(action.records, 'record_current', lambda *args: 0)
    with pytest.raises(ArithmeticError, match='current continuity'):
        action.native_binding(modular_clock.parameters())


def test_inference_does_not_read_the_simulator_weights(monkeypatch):
    data = counts()
    def forbidden(*args, **kwargs):
        raise AssertionError('private simulator weights must not enter inference')
    monkeypatch.setattr(action.modular_clock, 'parameters', forbidden)
    monkeypatch.setattr(action.local, 'reference', forbidden)
    assert action.infer_parent(data)['status'] == 'accept'
