"""Independent controls for clock agreement, source geometry and native reads."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import comparison as c


@pytest.fixture(scope='module')
def source():
    kernel = c.transfer.native_counts()
    return kernel, c.clock_correctors(kernel)


def test_poisson_equations_hold_for_every_source_state_and_port(source):
    kernel, cert = source
    h = cert['correctors']
    for v in range(12):
        assert sum(p*h[i][v] for i, p in enumerate(cert['palm'])) == 0
        for i, row in enumerate(cert['transitions']):
            direct_drift = sum(p*(abs(kernel.states[j][v]-kernel.states[i][v])-F(1, 6))
                               for j, p in row)
            assert direct_drift == h[i][v]-sum(p*h[j][v] for j, p in row)
    assert cert['range'] == F(38, 25)
    assert cert['variance_ceiling'] == F(30341677, 47045000)
    assert cert['increment_bound'] == F(241, 150)


def test_clock_covariances_recover_all_native_distance_shells(source):
    kernel, cert = source
    rows = c.distance_shells(kernel, cert['covariance'])
    assert [row['pair_count'] for row in rows] == [30, 30, 6]
    assert [row['asymptotic_covariance_per_movement'] for row in rows] == [
        '11281/261900', '-98989/1309500', '-144091/1309500']
    assert [row['asymptotic_difference_variance_per_movement'] for row in rows] == [
        '50101/109125', '608/873', '287/375']
    assert {tuple(p) for p in rows[0]['unordered_pairs']} == set(kernel.seams)
    assert all(sum(row) == 0 for row in cert['covariance'])
    values = np.linalg.eigvalsh(np.array(cert['covariance'], dtype=float))
    assert min(values) > -1e-14
    assert values[1] > 0  # only the exact total-count direction has zero variance
    # A clock-symmetric model with no distance information fails this test.
    mutated = tuple(tuple(F(1) if i == j else F(-1, 11) for j in range(12)) for i in range(12))
    with pytest.raises(ArithmeticError, match='cannot distinguish'):
        c.distance_shells(kernel, mutated)


def _finite_count_moments(kernel, transitions, prior, n, port):
    """Different method: propagate source populations and count moments."""
    mass = list(prior)
    first = [F(0)]*66
    second = [F(0)]*66
    for _ in range(n):
        future = [[F(0)]*66 for _ in range(3)]
        for i, row in enumerate(transitions):
            for j, probability in row:
                tick = abs(kernel.states[j][port]-kernel.states[i][port])
                future[0][j] += probability*mass[i]
                future[1][j] += probability*(first[i]+tick*mass[i])
                future[2][j] += probability*(second[i]+2*tick*first[i]+tick*mass[i])
        mass, first, second = future
    mean = sum(first)
    return mean, sum(second)-mean**2


@pytest.mark.parametrize('initial', ['palm', 0, 10, 65])
def test_arbitrary_preparation_mean_bound_and_finite_martingale_bound(source, initial):
    kernel, cert = source
    prior = cert['palm'] if initial == 'palm' else tuple(F(int(i == initial)) for i in range(66))
    n = 16
    mean, variance = _finite_count_moments(kernel, cert['transitions'], prior, n, 2)
    assert abs(mean-F(n, 6)) <= cert['range']
    # E(D^2) <= 2 E(M^2)+2 R^2, with no assumption of stationarity.
    assert variance+(mean-F(n, 6))**2 <= 2*n*cert['variance_ceiling']+2*cert['range']**2
    if initial == 'palm':
        assert mean == F(n, 6)
    else:
        assert mean != F(n, 6)  # finite initial bias is retained, not silently reset


def test_asymptotic_covariance_agrees_with_independent_finite_count_moments(source):
    kernel, cert = source
    n = 128
    mean, variance = _finite_count_moments(kernel, cert['transitions'], cert['palm'], n, 2)
    assert mean == F(n, 6)
    assert abs(variance/n-cert['covariance'][2][2]) < F(1, 100)
    assert variance/n > F(5, 36)  # independent Bernoulli ticks have the wrong noise


def test_exact_finite_confidence_bound_decays_without_a_time_unit(source):
    _, cert = source
    probabilities = [c.confidence_bound(cert, n, F(1, 60)) for n in (10000, 100000, 1000000)]
    assert 0 < probabilities[2] < probabilities[1] < probabilities[0] < 1
    assert probabilities[-1] < F(1, 400)
    assert c.confidence_bound(cert, 1, F(1, 60)) == 1


@pytest.mark.parametrize('law', ['C', 'E'])
def test_native_joint_clock_readouts_execute_all_complete_records(source, law):
    kernel, _ = source
    actual = c.instrument(kernel, law)
    original = c.quantum.native_instrument(law)
    assert sorted(k.tobytes() for group in actual.kraus for k in group) == sorted(
        k.tobytes() for group in original.kraus for k in group)
    initial = np.zeros((66, 66), complex)
    initial[0, 0] = 1
    history = c.runtime.execute(c.runtime.Program('complete_count_control', 66, initial,
        tuple(c.runtime.Event(str(i), instruments=(actual,)) for i in range(3))))
    expected = c.exact_attempt_records(kernel, 3)
    assert len(history.branches) == len(expected) == 125
    for branch in history.branches:
        np.testing.assert_allclose(branch.density, np.diag(list(map(float, expected[branch.record]))), atol=2e-15, rtol=0)
        assert dict(branch.resources) == {'native_attempts': 3, 'joint_clock_reads': 3}
    assert sum(sum(row) for row in expected.values()) == 1
    assert sum(expected[(0, 0, 0)]) > 0


@pytest.mark.parametrize('law', ['C', 'E'])
def test_initial_quantum_coherences_do_not_change_this_classical_clock_family(source, law):
    kernel, _ = source
    read = c.instrument(kernel, law)
    vector = np.zeros(66, complex)
    vector[0], vector[10] = 1/np.sqrt(2), 1j/np.sqrt(2)
    coherent = np.outer(vector, vector.conj())
    diagonal = np.diag(np.diag(coherent))
    events = tuple(c.runtime.Event(str(i), instruments=(read,)) for i in range(2))
    histories = [c.runtime.execute(c.runtime.Program('coherence_control', 66, rho, events))
                 for rho in (coherent, diagonal)]
    for a, b in zip(histories[0].branches, histories[1].branches, strict=True):
        assert a.record == b.record
        assert abs(a.probability-b.probability) < 2e-15


def test_independent_global_idling_cancels_at_every_source_movement(source):
    kernel, cert = source
    # Exact full-kernel cancellation, including state-dependent escape rates.
    for activity in (F(1), F(1, 2), F(1, 100)):
        for i, row in enumerate(kernel.marked_counts[1]):
            escape = activity*F(sum(row), 60)
            assert tuple((j, activity*F(count, 60)/escape) for j, count in enumerate(row) if count) == cert['transitions'][i]
    # Local relative clocks remain scheduler-sensitive: double one seam.
    uniform_seam_flux = (F(1, 30),)*30
    biased_seam_flux = (F(2, 31),)+(F(1, 31),)*29
    first_port = kernel.seams[0][0]
    assert sum(w for seam, w in zip(kernel.seams, uniform_seam_flux) if first_port in seam) == F(1, 6)
    assert sum(w for seam, w in zip(kernel.seams, biased_seam_flux) if first_port in seam) == F(6, 31)


@pytest.mark.parametrize('n,eps', [(True, F(1)), (0, F(1)), (1, 0.1), (1, F(0)), (1, F(-1))])
def test_invalid_clock_comparisons_fail_closed(source, n, eps):
    with pytest.raises(ValueError):
        c.confidence_bound(source[1], n, eps)
