"""Independent all-column checks of actual canonical scalar source words."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.relativistic_scalar import source
from model.relativistic_scalar.lattice import ScalarParameters


def oscillator(d, omega):
    # Assemble independently beyond all four intermediate ladder steps.
    a = np.zeros((d+6, d+6))
    for n in range(1, d+6):
        a[n-1, n] = math.sqrt(n)
    q = (a+a.T)/math.sqrt(2*omega)
    p = 1j*math.sqrt(omega/2)*(a.T-a)
    return q[:d, :d], (q@q)[:d, :d], (p@p)[:d, :d], np.linalg.matrix_power(q, 4)[:d, :d]


def check_ledger(metadata):
    resources = metadata["resources"]
    assert resources["native_pulses"] > 0
    assert resources["proper_code_load_events"] == resources["native_pulses"]
    assert resources["proper_code_export_events"] == resources["native_pulses"]
    assert resources["proper_code_service_events"] == 2*resources["native_pulses"]
    for key in source.COUNT_KEYS:
        assert resources[key] == sum(stage["resources"][key] for stage in metadata["source_stages"])
    assert metadata["full_phase_target_residual"] < 2e-10
    assert metadata["service_failure_outcomes_retained"]
    assert not metadata["physical_service_error_certified"]


@pytest.mark.parametrize("d", [2, 4, 8])
@pytest.mark.parametrize("coupling", [0., 1.])
def test_onsite_all_columns_use_physical_projected_powers(d, coupling):
    parameters = ScalarParameters(coupling=coupling)
    q, q2, p2, q4 = oscillator(d, parameters.omega)
    hamiltonian = p2/2+2*q2+parameters.g*q4
    matrix, metadata = source.onsite(d, 3, .073, parameters)
    np.testing.assert_allclose(matrix, expm(-.073j*hamiltonian), atol=3e-12, rtol=0)
    np.testing.assert_allclose(matrix.conj().T@matrix, np.eye(d), atol=3e-13, rtol=0)
    check_ledger(metadata)
    if coupling:
        wrong = p2/2+2*q2+parameters.g*np.linalg.matrix_power(q, 4)
        assert np.linalg.norm(matrix-expm(-.073j*wrong), 2) > 1e-3


def test_onsite_retains_absolute_phase_and_general_parameters():
    parameters = ScalarParameters(mass=1.25, speed=.7, spacing=.8, coupling=.3, omega=1.3)
    _, q2, p2, q4 = oscillator(4, parameters.omega)
    hamiltonian = p2/2+(parameters.mass**2+2*parameters.kappa)*q2/2+parameters.g*q4
    matrix, _ = source.onsite(4, 2, .12, parameters)
    np.testing.assert_allclose(matrix, expm(-.12j*hamiltonian), atol=4e-12, rtol=0)
    # A ground-energy subtraction preserves rays but violates this gate contract.
    shifted = expm(-.12j*(hamiltonian-np.linalg.eigvalsh(hamiltonian)[0]*np.eye(4)))
    assert np.linalg.norm(matrix-shifted, 2) > .05


@pytest.mark.parametrize("dimensions", [(2, 2), (8, 2), (8, 4), (4, 8)])
def test_bond_uses_separate_bases_and_little_endian_unequal_codes(dimensions):
    left, right = dimensions
    parameters = ScalarParameters(speed=.8, spacing=1.2, omega=1.6)
    qleft = oscillator(left, parameters.omega)[0]
    qright = oscillator(right, parameters.omega)[0]
    matrix, metadata = source.bond(left, right, .19, parameters)
    target = expm(.19j*parameters.kappa*np.kron(qright, qleft))
    np.testing.assert_allclose(matrix, target, atol=6e-12, rtol=0)
    check_ledger(metadata)
    backward, _ = source.bond(left, right, -.19, parameters)
    np.testing.assert_allclose(backward@matrix, np.eye(left*right), atol=8e-12, rtol=0)


@pytest.mark.parametrize("sign", [-1, 1])
def test_pointer_kraus_keep_both_outcomes_and_disturbed_reference_columns(sign):
    kraus, metadata = source.instrument(sign=sign)
    vector = np.zeros(8, complex)
    vector[0], vector[4] = 1/5, sign*1j*math.sqrt(24)/5
    projector = np.outer(vector, vector.conj())
    expected = (np.eye(8)-projector, projector)
    for actual, target in zip(kraus, expected, strict=True):
        np.testing.assert_allclose(actual, target, atol=4e-12, rtol=0)
    np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(8), atol=4e-13, rtol=0)
    rng = np.random.default_rng(137)
    reference = rng.normal(size=(8, 3))+1j*rng.normal(size=(8, 3))
    reference *= .37/np.linalg.norm(reference)
    for actual, target in zip(kraus, expected, strict=True):
        np.testing.assert_allclose(actual@reference, target@reference, atol=2e-12, rtol=0)
    assert abs(sum(np.linalg.norm(k@reference)**2 for k in kraus)-.37**2) < 1e-13
    # A projective read must keep the projected state, not just its probability.
    assert np.linalg.norm(kraus[0]@vector) < 4e-12
    np.testing.assert_allclose(kraus[1]@vector, vector, atol=4e-12, rtol=0)
    check_ledger(metadata)


def test_feedback_relative_phase_identity_control_and_mutation_isolation():
    matrix, metadata = source.feedback()
    expected = np.eye(8, dtype=complex)
    expected[4, 4] = 1j
    np.testing.assert_allclose(matrix, expected, atol=3e-12, rtol=0)
    check_ledger(metadata)
    baseline, baseline_metadata = source.feedback(phase=0)
    np.testing.assert_allclose(baseline, np.eye(8), atol=3e-12, rtol=0)
    assert baseline_metadata["resources"] == metadata["resources"]
    matrix[:] = 0
    metadata["resources"]["native_pulses"] = -1
    again, again_metadata = source.feedback()
    np.testing.assert_allclose(again, expected, atol=3e-12, rtol=0)
    assert again_metadata["resources"]["native_pulses"] > 0


def test_incorrect_source_phase_is_rejected_not_replaced_by_target(monkeypatch):
    original = source.c.execute

    def phase_fault(*args, **kwargs):
        return np.exp(.001j)*original(*args, **kwargs)

    monkeypatch.setattr(source.c, "execute", phase_fault)
    with pytest.raises(ArithmeticError, match="full-phase"):
        source.onsite(2, 3, .013579, ScalarParameters())


def test_declared_origin_code_and_input_validation():
    with pytest.raises(ValueError):
        source.instrument(4)
    with pytest.raises(ValueError):
        source.instrument(sign=True)
    with pytest.raises(ValueError):
        source.onsite(8, -1, .1)
    with pytest.raises(ValueError):
        source.bond(8, 4, float("nan"))
