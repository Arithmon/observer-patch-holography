"""Original public-probe regressions, independent of its equilibrium solver."""

from fractions import Fraction as F

import mpmath

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
