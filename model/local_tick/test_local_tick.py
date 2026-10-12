"""Independent stopping chains and a complete six-level quantum control."""
from fractions import Fraction as F
from itertools import combinations
import json

import numpy as np
from scipy.linalg import expm, null_space

from model.local_tick import experiment as subject


def classical_channels():
    states = tuple(combinations(range(12), 2))
    lookup = {state: index for index, state in enumerate(states)}
    channels = []
    for edge in subject.edges():
        swap = np.zeros((66, 66))
        a, b = edge
        for i, state in enumerate(states):
            outcome = tuple(sorted(b if v == a else a if v == b else v
                                   for v in state))
            swap[lookup[outcome], i] = 1
        channels.append((np.eye(66) + swap) / 2)
    return states, channels


def test_pinned_inputs_and_symbolic_identity():
    assert subject.verify_sources()
    assert subject.symbolic_shell_identity()


def test_all_observers_and_exact_closed_controls():
    for alpha in (F(0), F(1, 9), F(1, 2), F(1)):
        expected = subject.baseline(alpha)
        for observer in range(12):
            assert subject.mean_correlations(alpha, observer) == expected
        assert subject.corrected(*expected) == 0
    assert subject.baseline(1) == (F(5, 11), F(93, 407))


def test_remote_shell_matrix_from_actual_carrier():
    expected = np.array([[0, 0, 0, 0], [0, 2, -2, 0],
                         [0, -2, 3, -1], [0, 0, -5, 5]])
    seams = subject.edges()
    for observer in range(12):
        shells, seen = [[observer]], {observer}
        while len(seen) < 12:
            following = sorted({b if a in shells[-1] else a for a, b in seams
                                if (a in shells[-1] and b not in seen) or
                                   (b in shells[-1] and a not in seen)})
            assert following
            shells.append(following)
            seen.update(following)
        assert list(map(len, shells)) == [1, 5, 5, 1]
        inclusion = np.zeros((12, 4), dtype=int)
        for j, shell in enumerate(shells):
            inclusion[shell, j] = 1
        remote = np.zeros((12, 12), dtype=int)
        for a, b in seams:
            if observer not in (a, b):
                v = np.zeros(12, dtype=int)
                v[a], v[b] = 1, -1
                remote += np.outer(v, v)
        np.testing.assert_array_equal(remote @ inclusion, inclusion @ expected)


def test_full_66_state_endpoint_law_and_wrong_factor_order():
    states, channels = classical_channels()
    occupation = np.array([[int(p in state) for state in states]
                           for p in range(12)], dtype=float)
    for observer in (0, 3, 7):
        a = sum((channel / 30 for edge, channel in
                 zip(subject.edges(), channels) if observer not in edge),
                np.zeros((66, 66)))
        b = sum((channel / 30 for edge, channel in
                 zip(subject.edges(), channels) if observer in edge),
                np.zeros((66, 66)))
        resolvent = np.linalg.inv(np.eye(66) - a)
        stopped = b @ resolvent
        f = occupation[observer] - 1/6
        actual = tuple(float(f @ np.linalg.matrix_power(stopped, n) @ f /
                             (f @ f)) for n in (1, 2))
        np.testing.assert_allclose(actual, list(map(float, subject.baseline(1))),
                                   rtol=0, atol=2e-15)
        _, _, mean = subject.mean_maps(F(1), observer)
        np.testing.assert_allclose(occupation @ stopped,
                                   np.array(mean, dtype=float) @ occupation,
                                   rtol=0, atol=2e-15)
        # Correlations alone cannot catch a transposition of the full kernel.
        # Arbitrary preparations and all means do catch the reversed word.
        assert np.max(abs(occupation @ (resolvent @ b - stopped))) > 0.02
        np.testing.assert_allclose(stopped.sum(axis=0), 1, atol=2e-15)
        np.testing.assert_allclose(stopped.sum(axis=1), 1, atol=2e-15)


def test_finite_first_hit_series_keeps_the_failure_branch():
    _, channels = classical_channels()
    a = sum((channel / 30 for edge, channel in zip(subject.edges(), channels)
             if 0 not in edge), np.zeros((66, 66)))
    b = sum((channel / 30 for edge, channel in zip(subject.edges(), channels)
             if 0 in edge), np.zeros((66, 66)))
    full = b @ np.linalg.inv(np.eye(66) - a)
    for horizon in (1, 2, 7, 64):
        prefix, power = np.zeros((66, 66)), np.eye(66)
        for _ in range(horizon):
            prefix += b @ power
            power = a @ power
        remainder = full - prefix
        assert np.min(remainder) >= -2e-15
        np.testing.assert_allclose(remainder.sum(axis=0), (5/6)**horizon,
                                   rtol=0, atol=3e-15)


def small_quantum_channels():
    """Three one-particle ports; only middle port carries a private qubit.

    Both incident-pair algebras include that shared endpoint's private qubit.
    This dimension-six model checks the operator argument, not N4 numerics.
    Bare seam conditional expectations are constructed directly from Kraus
    matrices, then coherently dressed by a nonnormal exchange amplitude.
    """
    result = []
    for edge, amplitude in (((0, 1), np.array([[1, 2j], [0, 1]], complex)),
                            ((1, 2), np.array([[1j, 1], [0, 2]], complex))):
        a, b = edge
        generator = np.zeros((6, 6), complex)
        generator[2*a:2*a+2, 2*b:2*b+2] = amplitude
        generator[2*b:2*b+2, 2*a:2*a+2] = -amplitude.conj().T
        unitary = expm(0.37 * generator)
        spectator = next(v for v in range(3) if v not in edge)
        even = np.zeros((6, 6), complex)
        even[2*spectator:2*spectator+2, 2*spectator:2*spectator+2] = np.eye(2)
        kraus = [even]
        for left in edge:
            for right in edge:
                k = np.zeros((6, 6), complex)
                k[2*left:2*left+2, 2*right:2*right+2] = np.eye(2)/np.sqrt(2)
                kraus.append(unitary @ k @ unitary.conj().T)
        channel = sum(np.kron(k.conj(), k) for k in kraus)
        result.append(channel)
    return result


def test_complete_small_quantum_system_and_schur_identity():
    incident, remote = small_quantum_channels()
    for channel in (incident, remote):
        np.testing.assert_allclose(channel, channel.conj().T, atol=1e-15)
        np.testing.assert_allclose(channel @ channel, channel, atol=1e-15)
        assert np.linalg.eigvalsh(channel).min() > -1e-15
    fmat = np.diag([2/3, 2/3, -1/3, -1/3, -1/3, -1/3])
    second = np.diag([-1/3, -1/3, 2/3, 2/3, -1/3, -1/3])
    f = fmat.reshape(-1, order="F").astype(complex)
    basis = np.linalg.qr(np.column_stack((f, second.reshape(-1, order="F"))))[0]
    hidden = null_space(basis.conj().T)
    a, b, p = remote/2, incident/2, 0.5
    np.testing.assert_allclose(a @ f, (1-p)*f, atol=1e-15)
    dmat = np.eye(36) - a
    resolvent = np.linalg.inv(dmat)
    stopped = b @ resolvent  # Schrödinger convention, built from full Kraus maps.
    defect = (np.eye(36)-a-b) @ f
    norm = np.vdot(f, f).real
    r = np.vdot(f, defect).real/norm
    c1 = np.vdot(f, stopped @ f).real/norm
    c2 = np.vdot(f, stopped @ stopped @ f).real/norm
    assert abs(c1-(1-r/p)) < 1e-14
    expression = 1-2*r/p + np.vdot(defect, resolvent @ defect).real/(p*norm)
    assert abs(c2-expression) < 1e-14
    d00 = basis.conj().T @ dmat @ basis
    d10 = hidden.conj().T @ dmat @ basis
    d11 = hidden.conj().T @ dmat @ hidden
    h0, h1 = basis.conj().T @ defect, hidden.conj().T @ defect
    z = h1 - d10 @ np.linalg.solve(d00, h0)
    schur = d11 - d10 @ np.linalg.solve(d00, d10.conj().T)
    compressed = np.vdot(h0, np.linalg.solve(d00, h0)).real
    excess = (np.vdot(defect, resolvent @ defect).real - compressed)/(p*norm)
    independent = np.vdot(z, np.linalg.solve(schur, z)).real/(p*norm)
    assert excess > 0.0001
    assert abs(excess-independent) < 1e-14
    assert excess >= np.vdot(z, z).real/(p*norm) - 1e-14
    assert np.linalg.eigvalsh(schur).max() <= 1+1e-14
    # Independently sum all remote prefixes; this retains arbitrary coherences.
    prefix, power = np.zeros((36, 36), complex), np.eye(36)
    for _ in range(48):
        prefix += b @ power
        power = a @ power
    np.testing.assert_allclose(prefix, stopped, rtol=0, atol=1e-14)
    assert np.linalg.norm(stopped-stopped.conj().T) > 0.01


def test_full_result_replays_with_recomputed_interval_certificate():
    expected = json.loads((subject.HERE / "results.json").read_text())
    assert subject.calculate() == expected
    for row in expected["n4_certified_lower_bounds"]:
        assert F(row["certified_b_tick_lower_bound"]) > 0


def test_old_global_tick_correction_misclassifies_native_null():
    c1, c2 = subject.baseline(1)
    old_statistic = c2-c1*c1-(1-c1)**2/10
    assert old_statistic == F(-16, 2035)
    assert subject.corrected(c1, c2) == 0
