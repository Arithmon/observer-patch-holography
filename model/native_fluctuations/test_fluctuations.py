from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from . import source
from ..observer import native
from ..primitive_observer.transfer import native_counts


@pytest.fixture(scope='module')
def kernel():
    return native_counts()


@pytest.fixture(scope='module')
def moments(kernel):
    return source.source_moments(kernel)


def test_source_covariance_preserves_conservation_and_has_positive_innovations(moments):
    a, c, d = moments
    assert all(sum(row) == 1 for row in a)
    assert all(sum(row) == 0 for row in c+d)
    assert c[0][0] == F(5, 36)
    assert c[0][1] == -F(5, 396)
    eigenvalues = np.linalg.eigvalsh(np.asarray(d, float))
    assert abs(eigenvalues[0]) < 1e-14
    assert eigenvalues[1] > .013


def test_state_relabel_without_transporting_kernel_is_rejected(kernel):
    states = list(kernel.states)
    states[0], states[1] = states[1], states[0]
    with pytest.raises(ArithmeticError, match='conditional mean'):
        source.source_moments(replace(kernel, states=tuple(states)))


def test_two_record_law_against_direct_native_event_enumeration(kernel):
    index = {state: i for i, state in enumerate(kernel.states)}
    transitions = []
    for state in kernel.states:
        row = []
        for seam, coin in product(kernel.seams, (0, 1)):
            final, record = native.native(state, seam, coin)
            delta = -(int(seam[0] == 0)-int(seam[1] == 0))*record[1]
            row.append((index[final], delta))
        transitions.append(row)
    direct = Counter()
    for i in range(66):
        for j, first in transitions[i]:
            for k, second in transitions[j]:
                direct[first, second, k] += 1
    law = source.word_law(kernel, steps=2)
    for (first, second), terminal in law.items():
        for k, probability in enumerate(terminal):
            assert probability == F(direct[first, second, k], 66*60**2)
    assert sum(law[1, 1]) == 0
    variance = sum(word[0]**2*sum(row) for word, row in law.items())
    covariance = sum(word[0]*word[1]*sum(row) for word, row in law.items())
    assert variance == F(5, 198)
    assert covariance == -F(1, 792)
    assert covariance/variance == -F(1, 20)


def test_complete_runtime_output_and_checkpoint(kernel):
    report = source.execute_current_check(kernel)
    assert len(report['records']) == 27
    assert sum(F(row['probability']) for row in report['records']) == 1
    assert report['complete_terminal_density_residual'] < 1e-12
    assert report['checkpoint_distance'] < 1e-12
    assert any(F(row['probability']) == 0 for row in report['records'])


def test_stationary_initial_term_is_needed_for_the_history_precision(moments):
    a, c, d = moments
    for steps in (1, 3, 5):
        report = source.gaussian_action_check(a, c, d, steps)
        assert report['stationary_initial_term_included']
        assert report['precision_times_native_covariance_residual'] < 1e-11
    # Conditional increments alone leave the initial 11-dimensional state free.
    a = np.asarray(a, float)
    basis = np.linalg.qr(np.eye(12)[:, :11]-np.eye(12)[:, [11]])[0]
    ar = basis.T@a@basis
    transform = np.block([-ar, np.eye(11)])
    assert np.linalg.matrix_rank(transform.T@transform) == 11


def test_product_characteristic_agrees_with_complete_four_replica_convolution(kernel):
    law = source.score_distribution(source.word_law(kernel))
    aggregate = {0: F(1)}
    for _ in range(4):
        following = defaultdict(F)
        for x, p in aggregate.items():
            for y, q in law.items():
                following[x+y] += p*q
        aggregate = following
    assert sum(aggregate.values()) == 1
    mean = sum(x*p for x, p in aggregate.items())/2
    variance = sum(x*x*p for x, p in aggregate.items())/4
    native_variance, _, native_k4 = source.score_moments(law)
    assert mean == 0 and variance == native_variance
    assert sum(F(x**4, 16)*p for x, p in aggregate.items())-3*variance**2 == native_k4/4
    expected = sum(float(p)*np.exp(1j*x/2) for x, p in aggregate.items())
    actual = source.characteristic_comparison(law, 4)
    assert abs(expected-complex(actual['native_product_characteristic_real'],
                                actual['native_product_characteristic_imag'])) < 1e-14


@pytest.mark.parametrize('replicas', (1, 2, 16, 256))
@pytest.mark.parametrize('theta', (-2., -.1, 0., .5, 2.))
def test_analytic_characteristic_bound_also_covers_skew_nonGaussian_law(replicas, theta):
    result = source.characteristic_comparison({-1: F(2, 3), 2: F(1, 3)}, replicas, theta)
    assert result['difference_diagnostic'] <= result['analytic_characteristic_error_bound']+1e-14
    assert F(result['variance_exact']) == 2
    assert F(result['perfectly_correlated_replica_variance_control']) == 2*replicas


@pytest.mark.parametrize('invalid', (True, 0, -1, 1.5, 1000001))
def test_invalid_replication_is_rejected(invalid):
    with pytest.raises(ValueError):
        source.characteristic_comparison({-1: F(1, 2), 1: F(1, 2)}, invalid)


@pytest.mark.parametrize('invalid', (True, math.inf, math.nan, 11))
def test_invalid_Fourier_coordinate_is_rejected(invalid):
    with pytest.raises(ValueError):
        source.characteristic_comparison({-1: F(1, 2), 1: F(1, 2)}, 4, invalid)


@pytest.mark.parametrize('invalid', (True, -1, 12, .5))
def test_invalid_port_is_rejected(kernel, invalid):
    with pytest.raises(ValueError):
        source.current_instrument(kernel, invalid)


@pytest.mark.parametrize('invalid', (True, 0, 5, 1.5))
def test_invalid_record_horizon_is_rejected(kernel, invalid):
    with pytest.raises(ValueError):
        source.current_program(kernel, invalid)
    with pytest.raises(ValueError):
        source.word_law(kernel, invalid)


def test_noncentered_or_incomplete_score_is_rejected():
    with pytest.raises(ValueError, match='centered'):
        source.score_moments({0: F(1, 2), 1: F(1, 2)})
    with pytest.raises(ValueError, match='normalized'):
        source.score_moments({0: F(1, 2)})
