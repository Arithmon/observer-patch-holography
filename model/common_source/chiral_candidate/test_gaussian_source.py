"""Finite-Fock independent checks of the compressed source CAR evolution."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from ...observer.source_scalar import compiler
from . import gaussian_source as g


def ladders(modes):
    result = []
    for mode in range(modes):
        a = np.zeros((1 << modes, 1 << modes), complex)
        for mask in range(1 << modes):
            if mask & (1 << mode):
                a[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
        result.append(a)
    return result


def hamiltonian(h, delta):
    a = ladders(len(h))
    result = np.zeros((1 << len(h),) * 2, complex)
    for i in range(len(h)):
        for j in range(len(h)):
            result += h[i, j]*a[i].conj().T@a[j]
            result += (delta[i, j]*a[i].conj().T@a[j].conj().T
                       +delta[i, j].conjugate()*a[j]@a[i])/2
    return result


def covariance(state):
    a = ladders(int(math.log2(len(state))))
    w = [value for c in a for value in (c+c.conj().T, -1j*(c-c.conj().T))]
    return np.array([[(.5j*np.vdot(state, (left@right-right@left)@state)).real
                      for right in w] for left in w])


@pytest.mark.parametrize("kind,value", [("onsite", .37), ("hopping", .37),
                                        ("hopping", .31-.19j), ("pairing", .37),
                                        ("pairing", -.13+.41j)])
def test_actual_local_four_column_source_equals_independent_quadratic_Fock_action(kind, value):
    h, delta = np.zeros((2, 2), complex), np.zeros((2, 2), complex)
    if kind == "onsite":
        h[0, 0] = value
    elif kind == "hopping":
        h[0, 1], h[1, 0] = value, np.conj(value)
    else:
        delta[0, 1], delta[1, 0] = value, -value
    gate = g.local_gate(kind, value, .43)
    expected = expm(-.43j*hamiltonian(h, delta))
    np.testing.assert_allclose(gate.fock_unitary, expected, atol=3e-14, rtol=0)
    assert gate.fock_unitary.shape == (4, 4)
    assert gate.diagnostics["contracted_fock_columns"] == 4
    assert gate.diagnostics["majorana_linearity_residual"] < 1e-14
    assert gate.diagnostics["majorana_orthogonality_residual"] < 1e-13
    assert not gate.diagnostics["polar_or_norm_correction_applied"]
    # Every completed source map preserves all CAR and the even parity block.
    parity = np.diag([1., -1., -1., 1.])
    np.testing.assert_allclose(gate.fock_unitary@parity, parity@gate.fock_unitary, atol=2e-14)


def test_source_FSWAP_has_essential_double_occupation_minus_sign():
    swap = g.local_gate("fswap")
    desired = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, -1]])
    np.testing.assert_allclose(swap.fock_unitary, desired, atol=3e-14, rtol=0)
    np.testing.assert_allclose(swap.orthogonal, np.eye(4)[[2, 3, 0, 1]], atol=3e-14)
    assert abs(swap.fock_unitary[3, 3]+1) < 3e-14
    assert swap.resources["native_pulses"] == len(swap.word)


@pytest.mark.parametrize("kind", ["hopping", "pairing"])
def test_nonadjacent_source_route_matches_full_four_mode_CAR_and_charges_every_swap(kind):
    h, delta = np.zeros((4, 4), complex), np.zeros((4, 4), complex)
    coefficient, duration = .27+.13j, .31
    if kind == "hopping":
        h[0, 3], h[3, 0] = coefficient, coefficient.conjugate()
    else:
        delta[0, 3], delta[3, 0] = coefficient, -coefficient
    gate, swap = g.local_gate(kind, coefficient, duration/2), g.local_gate("fswap")
    route = [(2, 3), (1, 2)]
    half = [event for pair in route for event in compiler.remap(swap.word, list(pair))]
    half += compiler.remap(gate.word, [0, 1])
    half += [event for pair in reversed(route) for event in compiler.remap(swap.word, list(pair))]
    actual_fock = compiler.execute(half+half, 4)
    target_fock = expm(-1j*duration*hamiltonian(h, delta))
    np.testing.assert_allclose(actual_fock, target_fock, atol=2e-13, rtol=0)
    source = g.evolve(h, delta, duration)
    a = ladders(4)
    w = [value for c in a for value in (c+c.conj().T, -1j*(c-c.conj().T))]
    recovered = np.array([[np.trace(right@actual_fock.conj().T@left@actual_fock).real/16
                            for right in w] for left in w])
    np.testing.assert_allclose(source.orthogonal, recovered, atol=3e-13, rtol=0)
    assert source.resources["routing_fswap_calls"] == 8
    assert source.resources["completed_local_blocks"] == 10
    assert source.resources["native_pulses"] == 8*swap.resources["native_pulses"]+2*gate.resources["native_pulses"]
    assert source.resources["proper_code_service_events"] == 2*source.resources["native_pulses"]


def test_complex_noncommuting_source_Strang_converges_and_keeps_reference_cross_covariances():
    h = np.array([[.13, .21+.08j, -.11j], [.21-.08j, -.17, .19], [.11j, .19, .09]])
    delta = np.array([[0, .18j, .23], [-.18j, 0, -.07+.03j], [-.23, .07-.03j, 0]])
    # The fourth mode is an entangled Gaussian reference, not a traced reservoir.
    prep_h, prep_d = np.zeros((4, 4), complex), np.zeros((4, 4), complex)
    prep_d[0, 3], prep_d[3, 0] = 1., -1.
    psi = expm(-.51j*hamiltonian(prep_h, prep_d))[:, 0]
    initial = covariance(psi)
    full_h, full_d = np.zeros((4, 4), complex), np.zeros((4, 4), complex)
    full_h[:3, :3], full_d[:3, :3] = h, delta
    target = covariance(expm(-.4j*hamiltonian(full_h, full_d))@psi)
    errors = []
    for steps in (1, 2, 4):
        result = g.evolve(h, delta, .4, steps=steps, covariance=initial)
        errors.append(np.max(abs(result.covariance-target)))
        np.testing.assert_allclose(result.covariance[6:, 6:], initial[6:, 6:], atol=1e-14)
        assert np.linalg.norm(result.covariance[:6, 6:]) > .1
        assert result.resources["retained_reference_fermion_modes"] == 1
        assert not result.diagnostics["polar_or_covariance_normalization_applied"]
    assert errors[1] < errors[0]/3.9 and errors[2] < errors[1]/3.9


def test_exact_Gaussian_occupation_conditionals_match_all_dense_projector_outputs():
    h = np.array([[.1, .23j, 0], [-.23j, 0, .17], [0, .17, -.1]])
    delta = np.array([[0, .41, .12j], [-.41, 0, -.31], [-.12j, .31, 0]])
    psi = expm(-.57j*hamiltonian(h, delta))[:, 0]
    incoming = covariance(psi)
    read = g.occupation_read(incoming, 0, weight=2/7)
    for branch in read.branches:
        projected = psi.copy()
        projected[[mask for mask in range(8) if (mask & 1) != branch.occupation]] = 0
        probability = np.vdot(projected, projected).real
        assert branch.probability == pytest.approx(probability, abs=2e-15)
        assert branch.weight == pytest.approx(2*probability/7, abs=2e-15)
        np.testing.assert_allclose(branch.covariance, covariance(projected/math.sqrt(probability)), atol=3e-14)
    assert sum(branch.weight for branch in read.branches) == pytest.approx(2/7, abs=1e-15)
    np.testing.assert_allclose(read.source_kraus[0], np.diag([1., 0.]), atol=1e-14)
    np.testing.assert_allclose(read.source_kraus[1], np.diag([0., 1.]), atol=1e-14)
    assert read.resources["fresh_pointer_blanks"] == 1
    assert read.resources["source_mode_preparations"] == 0
    assert read.diagnostics["maximum_conditional_CAR_bound_excess"] < 1e-13


def test_zero_occupation_outcome_remains_and_invalid_probabilities_are_not_repaired():
    read = g.occupation_read(g.vacuum_covariance(2), 0, weight=.3)
    assert len(read.branches) == 2
    assert read.branches[0].weight == .3
    assert read.branches[1].weight == 0 and read.branches[1].covariance is None
    bad = g.vacuum_covariance(1)*(1+1e-11)
    with pytest.raises(ValueError, match="no clipping"):
        g.occupation_read(bad, 0)


def test_ninety_modes_execute_only_local_four_column_Fock_blocks(monkeypatch):
    original = compiler.execute
    calls = []
    def checked(word, qubits, *args, **kwargs):
        calls.append(qubits)
        assert qubits == 2
        return original(word, qubits, *args, **kwargs)
    monkeypatch.setattr(compiler, "execute", checked)
    h, delta = np.zeros((90, 90), complex), np.zeros((90, 90), complex)
    delta[0, 89], delta[89, 0] = .314159+.271828j, -.314159-.271828j
    result = g.evolve(h, delta, .1234567, covariance=g.vacuum_covariance(92))
    assert calls and set(calls) == {2}
    assert result.orthogonal.shape == (180, 180)
    assert result.covariance.shape == (184, 184)
    assert result.resources["routing_fswap_calls"] == 352
    assert result.resources["maximum_executed_local_Fock_dimension"] == 4
    assert result.resources["retained_reference_fermion_modes"] == 2
    assert result.diagnostics["full_Fock_dimension_if_allocated"] == str(2**90)
    assert not result.record()["full_Fock_state_allocated"]
    assert not result.record()["target_BdG_exponential_used_by_source"]
    expected = math.sin(abs(.314159+.271828j)*.1234567)**2
    assert g.occupations(result.covariance)[0] == pytest.approx(expected, abs=5e-12)
    assert g.occupations(result.covariance)[89] == pytest.approx(expected, abs=5e-12)


def test_cached_source_objects_and_input_domains_cannot_be_corrupted():
    gate = g.local_gate("pairing", .372, .43)
    gate.resources["native_pulses"] = 0
    gate.word[0]["angle"] = 123.
    again = g.local_gate("pairing", .372, .43)
    assert again.resources["native_pulses"] > 0
    assert again.word[0]["angle"] != 123.
    with pytest.raises(ValueError):
        again.fock_unitary.setflags(write=True)
    for kind, value, time in (("onsite", 1j, .2), ("fswap", 1., 0.),
                               ("pairing", math.nan, .2), ("hopping", .1, True)):
        with pytest.raises(ValueError):
            g.local_gate(kind, value, time)
    h, delta = np.zeros((2, 2), complex), np.zeros((2, 2), complex)
    with pytest.raises(ValueError):
        g.evolve(h, delta, .2, steps=True)
    with pytest.raises(ValueError):
        g.evolve(h+np.array([[0, 1], [0, 0]]), delta, .2)
    with pytest.raises(ValueError):
        g.evolve(h, delta+np.eye(2), .2)
    with pytest.raises(ValueError):
        g.occupation_read(g.vacuum_covariance(2), True)
