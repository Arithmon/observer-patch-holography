"""Independent complete-Fock, source-export, reference and entropy controls."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import thermal as t
from ..full_fock_field import field
from ..observer.source_scalar import compiler
from ..reference_observation.observation import ReferenceObservation


def ladder(modes, at):
    # Independent Jordan--Wigner tensor assembly, highest bit first.
    annihilation = np.array([[0, 1], [0, 0]], complex)
    parity = np.diag([1, -1])
    factors = [annihilation if i == at else parity if i < at else np.eye(2)
               for i in reversed(range(modes))]
    out = factors[0]
    for factor in factors[1:]:
        out = np.kron(out, factor)
    return out


@pytest.mark.parametrize('axis', (0, 1, 2))
@pytest.mark.parametrize('interaction,beta', ((0., .7), (1., 1.3)))
def test_complete_Gibbs_reference_against_independent_CAR_exponential(axis, interaction, beta):
    sides = tuple(2 if i == axis else 1 for i in range(3))
    geometry = field.FieldGeometry(sides, speed=.25, interaction=interaction)
    case = t.ThermalField(geometry, beta)
    a = [ladder(4, i) for i in range(4)]
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))[axis]
    h = np.zeros((16, 16), complex)
    for i in range(2):
        for j in range(2):
            hop = 1j*.25/2*pauli[i, j]*(a[i].conj().T@a[j+2])
            h += hop+hop.conj().T
    number = [x.conj().T@x for x in a]
    h += interaction*(number[0]@number[1]+number[2]@number[3])
    density = expm(-beta*h); density /= np.trace(density)
    root = expm(-beta*h/2)/math.sqrt(float(np.trace(expm(-beta*h)).real))
    target = case.target()
    np.testing.assert_allclose(target['hamiltonian'], h, atol=2e-15)
    np.testing.assert_allclose(target['reference'], density, atol=2e-15)
    np.testing.assert_allclose(target['purification'], root.T.reshape(-1), atol=2e-15)
    prepared = case.prepare()
    np.testing.assert_allclose(prepared['state'], root.T.reshape(-1), atol=3e-11, rtol=0)
    np.testing.assert_allclose(prepared['success_field_density'], density, atol=3e-11, rtol=0)
    assert len(prepared['first_export_failures']) == 10249
    assert prepared['resources']['proper_code_service_events'] == 20498
    assert np.linalg.eigvalsh(density)[0] > 0
    assert all(density[mask, mask].real > 0 for mask in range(16))


def test_complete_export_vectors_retain_first_failure_and_untouched_reference(monkeypatch):
    # Deliberately replace the proper-code-preserving primitive by a leaking
    # unitary. Full Gram conservation must include the two abort blocks.
    def leaking(kind, angle):
        size = 4 if kind == 'J' else 2
        generator = np.zeros((6, 6), complex)
        generator[0, size] = generator[size, 0] = 1
        return expm(1j*angle*generator)
    monkeypatch.setattr(compiler, 'native_pulse', leaking)
    initial = np.array([[1, 1j], [2j, 3], [4, -2j], [3j, 5]], complex)/10
    word = [compiler.pulse('X', .2, [0]), compiler.pulse('J', -.4, [1, 0])]
    success, failures = t.execute_with_exports(word, 2, initial)
    first = initial.copy(); first[[0, 2]] *= math.cos(.2)
    expected = first.copy(); expected[0] *= math.cos(.4)
    np.testing.assert_allclose(success, expected, atol=1e-15)
    expected_first = np.zeros((8, 2), complex)
    expected_first[:2] = 1j*math.sin(.2)*initial[[0, 2]]
    expected_second = np.zeros((2, 2), complex)
    expected_second[0] = -1j*math.sin(.4)*first[0]
    np.testing.assert_allclose(failures[0], expected_first, atol=1e-15)
    np.testing.assert_allclose(failures[1], expected_second, atol=1e-15)
    total_gram = success.conj().T@success+sum(x.conj().T@x for x in failures)
    np.testing.assert_allclose(total_gram, initial.conj().T@initial, atol=1e-15)
    assert np.linalg.norm(success.conj().T@success-initial.conj().T@initial) > .005


def test_preparation_rejects_wrong_complete_source_amplitudes(monkeypatch):
    original = t.execute_with_exports
    def damaged(*args):
        state, failures = original(*args)
        state = state.copy(); state[0] *= 1j
        return state, failures
    monkeypatch.setattr(t, 'execute_with_exports', damaged)
    with pytest.raises(ArithmeticError):
        t.ThermalField().prepare()


def test_purifier_is_retained_and_conditioning_does_not_give_the_Gibbs_state():
    case = t.ThermalField()
    target = case.target()
    coefficients = target['purification'].reshape(16, 16).T
    branches = [np.outer(coefficients[:, j], coefficients[:, j].conj()) for j in range(16)]
    np.testing.assert_allclose(sum(branches), target['reference'], atol=1e-15)
    for branch in branches:
        weight = np.trace(branch).real
        assert weight > 0
        np.testing.assert_allclose(branch@branch, weight*branch, atol=1e-15)
        assert np.linalg.norm(branch/weight-target['reference']) > .1
    assert t.entropy(target['reference']) > 1
    # Coherent purification is pure on the joint algebra, although the field
    # functional is faithful. This distinction is not a numerical cutoff.
    state = target['purification']
    joint = np.outer(state, state.conj())
    np.testing.assert_allclose(joint@joint, joint, atol=2e-15)


@pytest.mark.parametrize('kept', ((0,), (3,), (1, 3), (0, 2, 3), (0, 1, 2, 3)))
def test_signed_restriction_preserves_the_actual_local_CAR_algebra(kept):
    kraus = t.car_restriction(4, kept)
    np.testing.assert_array_equal(sum(k.conj().T@k for k in kraus), np.eye(16))
    for index, mode in enumerate(kept):
        embedded = sum(k.conj().T@ladder(len(kept), index)@k for k in kraus)
        np.testing.assert_array_equal(embedded, ladder(4, mode))
    if kept == (1, 3):
        unsigned = tuple(abs(k) for k in kraus)
        wrong = sum(k.T@ladder(2, 0)@k for k in unsigned)
        assert np.linalg.norm(wrong-ladder(4, 1)) > 1


def test_full_CAR_state_and_reference_towers_on_all_units_and_complex_references():
    rho = t.ThermalField().target()['reference']
    first = t.car_restriction(4, (1, 3))
    second = t.car_restriction(2, (1,))
    direct = t.car_restriction(4, (3,))
    for i in range(16):
        for j in range(16):
            unit = np.zeros((16, 16), complex); unit[i, j] = 1
            np.testing.assert_array_equal(t.restrict(t.restrict(unit, first), second), t.restrict(unit, direct))
    a = ReferenceObservation(rho, first)
    b = ReferenceObservation(t.restrict(rho, first), second)
    ab = ReferenceObservation(rho, direct)
    random = np.random.default_rng(41)
    x = random.normal(size=(32, 32))+1j*random.normal(size=(32, 32))
    np.testing.assert_allclose(b.observe(a.observe(x, 2), 2), ab.observe(x, 2), atol=3e-14)
    np.testing.assert_allclose(ab.recover(ab.target_reference), rho, atol=2e-15)
    # Boundary hopping changes number-read outputs, so stationarity under H
    # does not make the reference invariant under every observer instrument.
    n = field.FullFock(4).number(0).toarray()
    measured = n@rho@n+(np.eye(16)-n)@rho@(np.eye(16)-n)
    assert np.linalg.norm(measured-rho) > .001


def test_same_action_noncommuting_first_law_and_complete_executed_read():
    result = t.ThermalField().report()
    check = result['checks']
    assert check['noncommuting_variation_norm'] > .001
    assert check['onsite_non_Gaussian_Wick_defect'] < -.05
    assert check['first_law_identity_residual'] < 1e-13
    assert 0 < check['quantum_relative_entropy'] < check['relative_entropy_remainder_upper']
    assert abs(check['complete_pointer_and_abort_probability_sum']-1) < 1e-10
    assert result['continuation_resources']['native_pulses'] == 1411
    assert not result['scope']['geometric_G1_G6_or_Einstein_premises_discharged']
    assert not result['scope']['joint_purification_faithful_on_full_algebra']
    assert not result['scope']['reference_observation_maps_source_compiled']


@pytest.mark.parametrize('beta', (0, -1, True, float('nan'), float('inf'), 1e5))
def test_bad_or_unresolved_faithful_temperature_fails_closed(beta):
    with pytest.raises(ValueError):
        t.ThermalField(beta=beta).target()


def test_bounded_domains_masks_and_incomplete_operations_fail_closed():
    with pytest.raises(ValueError):
        t.ThermalField(field.FieldGeometry((2, 2, 1)))
    for kept in ((), (1, 0), (1, 1), (True,), (4,), [1]):
        with pytest.raises(ValueError):
            t.car_restriction(4, kept)
    with pytest.raises(ValueError):
        t.restrict(np.eye(4)/4, [np.eye(4)/2])
    with pytest.raises(ValueError):
        t.restrict(np.eye(4)/4, [np.ma.array(np.eye(4), mask=True)])
    with pytest.raises(ValueError):
        t.entropy(np.diag([1., 0.]))
    with pytest.raises(ValueError):
        t.entropy(np.ma.array(np.eye(2)/2, mask=True))
    with pytest.raises(ValueError):
        t.execute_with_exports([], True, [1., 0.])
    with pytest.raises(ValueError):
        t.execute_with_exports([], 2, np.zeros((4, 0)))
    with pytest.raises(ValueError):
        t.execute_with_exports([compiler.pulse('X', 0., [0])]*100, 9, np.zeros((512, 512)))
