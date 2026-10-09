"""Original public-probe regressions, independent of its equilibrium solver."""

from fractions import Fraction as F

import mpmath
import pytest

import collar_matrix_realization_probe as probe


def test_slow_chain_equilibrium_is_not_a_small_residual_iterate():
    epsilon = F(1, 2**48)
    kernel = [[1-epsilon, epsilon], [2*epsilon, 1-2*epsilon]]
    # Balance the only cut: pi[0] * epsilon = pi[1] * 2*epsilon.
    assert probe.stationary_of(kernel) == [F(2, 3), F(1, 3)]


def test_public_audit_uses_the_actual_stationary_reference():
    epsilon = F(1, 2**48)
    kernel = [[1-epsilon, epsilon], [2*epsilon, 1-2*epsilon]]
    ctx = mpmath.mp.clone()
    ctx.dps = 70
    report = probe.audit_irreducible_chain(kernel, ctx)
    # D((1/2,1/2)||(2/3,1/3)) = log(9/8)/2, not numerical zero.
    expected = ctx.log(ctx.mpf(9)/8)/2
    assert abs(ctx.mpf(report['kl_to_stationary_initial'])-expected) < ctx.mpf('1e-12')


def test_loopless_period_one_chain_is_recognized():
    # Return walks 0->1->0 and 0->1->2->0 have coprime lengths.
    kernel = [[F(0), F(1), F(0)],
              [F(1, 2), F(0), F(1, 2)],
              [F(1), F(0), F(0)]]
    report = probe.audit_irreducible_chain(kernel, mpmath.mp.clone())
    assert report['aperiodic'] is True
    assert report['aperiodic_self_loop_witness'] is None


def test_count_aggregation_preserves_original_integer_masses():
    counts = [[2**54-1, 1], [2, 2**54-2]]
    labels = [[['record_family', 0]], [['record_family', 1]]]
    _, aggregated, kernel = probe.coarsen_counts(counts, labels, ('record_family',))
    assert aggregated == counts
    assert kernel == [[F(2**54-1, 2**54), F(1, 2**54)],
                      [F(2, 2**54), F(2**54-2, 2**54)]]


@pytest.mark.parametrize('exponent', [48, 400, 1000])
def test_public_entropy_loss_preserves_slow_positive_dissipation(exponent):
    epsilon = F(1, 2**exponent)
    kernel = [[1-epsilon, epsilon], [2*epsilon, 1-2*epsilon]]
    report = probe.audit_irreducible_chain(kernel, mpmath.mp.clone())
    # Independent closed-form trajectory; 500 digits resolve the difference
    # between the two full entropies even in the 2^-1000 case.
    ctx = mpmath.mp.clone()
    ctx.dps = 500
    eps = ctx.mpf(1)/2**exponent
    entropies = []
    for step in range(probe.KL_STEPS+1):
        p0 = ctx.mpf(2)/3-(1-3*eps)**step/6
        entropies.append(p0*ctx.log(3*p0/2)+(1-p0)*ctx.log(3*(1-p0)))
    expected = min(a-b for a, b in zip(entropies, entropies[1:]))
    actual = ctx.mpf(report['kl_min_stepwise_descent'])
    assert expected > 0
    assert abs(actual/expected-1) < ctx.mpf('1e-10')


def test_periodic_nonstationary_trajectory_can_have_exactly_zero_entropy_loss():
    kernel = [[F(0), F(1, 2), F(1, 2)],
              [F(1), F(0), F(0)], [F(1), F(0), F(0)]]
    # Uniform initial law oscillates with (2/3,1/6,1/6). Its density relative
    # to (1/2,1/4,1/4) is constant on each cyclic class, so no information is
    # lost by the transition, despite not being at equilibrium.
    ctx = mpmath.mp.clone()
    report = probe.audit_irreducible_chain(kernel, ctx)
    assert report['period'] == 2
    assert report['stationary_distribution_exact'] == ['1/2', '1/4', '1/4']
    assert ctx.mpf(report['kl_to_stationary_initial']) > 0
    assert report['kl_to_stationary_initial'] == report['kl_to_stationary_final']
    assert ctx.mpf(report['kl_min_stepwise_descent']) == 0


def test_nonreversible_entropy_loss_agrees_with_direct_high_precision_path():
    weights = [[8, 1, 0], [0, 4, 1], [1, 0, 3]]
    kernel = [[F(value, sum(row)) for value in row] for row in weights]
    # Equal circulation on three directed off-diagonal edges determines pi.
    pi = [F(9, 18), F(5, 18), F(4, 18)]
    ctx = mpmath.mp.clone()
    ctx.dps = 150
    probabilities = [ctx.mpf(1)/3]*3
    reference = [ctx.mpf(x.numerator)/x.denominator for x in pi]
    transition = [[ctx.mpf(value)/sum(row) for value in row] for row in weights]
    entropies = []
    for _ in range(probe.KL_STEPS+1):
        entropies.append(ctx.fsum(p*ctx.log(p/q) for p, q in zip(probabilities, reference)))
        probabilities = [ctx.fsum(probabilities[i]*transition[i][j] for i in range(3))
                         for j in range(3)]
    expected = min(a-b for a, b in zip(entropies, entropies[1:]))
    report = probe.audit_irreducible_chain(kernel, ctx)
    assert report['reversible'] is False
    assert report['stationary_distribution_exact'] == [str(x) for x in pi]
    assert expected > 0
    actual = ctx.mpf(report['kl_min_stepwise_descent'])
    assert abs(actual/expected-1) < ctx.mpf('1e-10')


def test_entropy_diagnostics_refuse_a_nonstationary_reference():
    kernel = [[F(3, 4), F(1, 4)], [F(1, 2), F(1, 2)]]
    with pytest.raises(probe.ThermoError, match='stationary reference'):
        probe.kl_diagnostics(kernel, [F(1, 2), F(1, 2)], mpmath.mp.clone())
