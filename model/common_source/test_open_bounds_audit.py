"""Independent full-output and original-parameter audits of the noise bounds."""
from dataclasses import replace
from fractions import Fraction as F
import math

import mpmath as mp
import numpy as np
import pytest

from . import open_bounds as b, open_environment as e, open_protocol as runtime, protocol
from .test_open_protocol import direct_gaussian, entangled, short_program


@pytest.mark.parametrize('dimension', [2, 4])
def test_complete_noisy_product_history_obeys_reference_uniform_finite_budget(dimension):
    """Compare every quantum block with an independent continuous/Gaussian oracle.

    The oracle constructs full tensor maps and never invokes the production
    environment or history engine. The source channel in the actual run is
    the finite M=2 rule; its difference is not a mere detector-probability test.
    """
    program = replace(short_program(dimension), steps=128)
    events = e.Environment(tuple(e.PhaseNoise(after, site, F(1, 64), 2)
                                 for after, site in ((0, 1), (2, 0), (4, 1))))
    initial = entangled(dimension)
    actual = runtime.execute(program, events, backend='product', initial=initial)
    expected = direct_gaussian(program, events, initial)
    differences = [branch.density-expected[branch.record] for branch in actual.branches]
    full_distance = sum(np.linalg.svd(delta, compute_uv=False).sum() for delta in differences)/2
    classical_distance = sum(abs(np.trace(delta)) for delta in differences)/2
    report = b.budget(program, events)
    bound = F(report['ideal_finite_source_to_continuous_GAUSSIAN_history_error_upper'])
    assert 0 < full_distance <= float(bound)+3e-12 < 1
    assert full_distance > classical_distance+1e-7
    assert report['finite_comparison_includes_arbitrary_work_reference_inputs']
    assert not report['floating_execution_error_certified']


@pytest.mark.parametrize('site', [0, 1])
def test_noise_dilation_intertwines_cutoff_and_preserves_projected_error(site):
    """Check the algebra behind the unchanged moment/Duhamel induction.

    Use a larger rectangular occupation space with amplitudes both inside
    and outside the smaller one. No Hamiltonian solution or estimated tail
    is reused to establish the exact phase/projection intertwining.
    """
    small, large, order, variance = 2, 4, 4, F(1, 7)
    rows = [g+4*(n0+large*n1) for n1 in range(small)
            for n0 in range(small) for g in range(4)]
    projection = np.eye(4*large*large)[rows]
    channels = [e.source.channel(d, 1., float(variance), order) for d in (small, large)]
    maps = []
    for d, channel in zip((small, large), channels, strict=True):
        maps.append(tuple(np.kron(np.kron(u, np.eye(d)), np.eye(4)) if site == 1
                          else np.kron(np.kron(np.eye(d), u), np.eye(4))
                          for u in channel.unitaries))
    random = np.random.default_rng(318+site)
    infinite_like = random.normal(size=(4*large*large, 3))+1j*random.normal(size=(4*large*large, 3))
    finite = random.normal(size=(4*small*small, 3))+1j*random.normal(size=(4*small*small, 3))
    before = np.linalg.norm(finite-projection@infinite_like)**2
    after = 0.
    for weight, low, high in zip(channels[0].weights, *maps, strict=True):
        assert np.max(np.abs(low@projection-projection@high)) < 3e-15
        after += weight*np.linalg.norm(low@finite-projection@high@infinite_like)**2
    assert after == pytest.approx(before, rel=3e-15)


@pytest.mark.parametrize('dimension', [2, 3, 4, 8, 16, 32])
def test_refinement_rule_bounds_original_factorial_remainder(dimension):
    variances = (F(0), F(1, 64), F(1, 3))
    report = b.gaussian_refinement(dimension, variances)
    exact_remainder, geometric = F(0), F(0)
    for variance, order in zip(variances, report['orders'], strict=True):
        assert order >= dimension
        assert F(order) >= 2*variance*(dimension-1)**2
        exact_remainder += F(dimension*math.factorial(order), 2*math.factorial(2*order))*(
            variance*(dimension-1)**2)**order
        geometric += F(0) if variance == 0 else F(dimension, 2**(order+1))
    declared = F(report['ideal_exact_parameter_quadrature_error_upper'])
    assert min(F(1), exact_remainder) <= declared
    grid = F(report['outward_rounding_grid'])
    assert min(F(1), geometric) <= declared <= min(F(1), geometric+2*grid)
    assert report['parameter_rounding_control_required']


def test_large_refinement_orders_have_bounded_outward_reports_that_converge():
    previous = F(1)
    for dimension in (128, 256, 1024, 4096):
        report = b.gaussian_refinement(dimension, (F(1), F(0), F(10**100)))
        radius = F(report['ideal_exact_parameter_quadrature_error_upper'])
        grid = F(report['outward_rounding_grid'])
        assert radius == 2*grid
        assert 0 < radius < previous
        assert report['orders'][2] > 10**100
        assert len(report['ideal_exact_parameter_quadrature_error_upper']) < 100
        assert not report['source_execution_available']
        previous = radius


def test_order_one_identity_does_not_become_a_zero_gaussian_error():
    program = protocol.continuation_program(2, 128)
    schedule = e.Environment((e.PhaseNoise(4, 0, F(1, 4), 1),))
    report = b.budget(program, schedule, prepared_noiseless_error=F(1, 5))
    comparison = report['prepared_noiseless_reference_comparison']
    assert comparison['ideal_finite_noise_to_noiseless_finite_error_upper'] == '0'
    assert comparison['ideal_noisy_finite_to_NOISELESS_infinite_error_upper'] == '1/5'
    assert comparison['target_contains_phase_noise'] is False
    assert F(report['environment_budget']['ideal_finite_to_gaussian_half_diamond_upper']) > 0


def test_underflow_conversion_covers_original_positive_gaussian_effect():
    """A positive original-input error must survive a float-zero implementation."""
    context = mp.mp.clone()
    context.dps = 80
    variance = F(1, 2**1075)
    result = e.event_metadata(protocol.Program(), e.PhaseNoise(0, 0, variance, 2))
    assert result['represented_variance_exact'] == '0'
    declared = F(result['ideal_finite_to_gaussian_half_diamond_upper'])
    exact_variance = context.mpf(variance.numerator)/variance.denominator
    # Bell coherence on levels zero/one: its trace distance to the Gaussian
    # channel is (1-exp(-v/2))/2. expm1 preserves the nonzero tiny difference.
    bell_error = -context.expm1(-exact_variance/2)/2
    upper = context.mpf(declared.numerator)/declared.denominator
    assert 0 < bell_error <= upper
    neighbor = e.event_metadata(protocol.Program(), e.PhaseNoise(0, 0, 2*variance, 2))
    assert F(neighbor['represented_variance_exact']) == 2*variance
    assert F(neighbor['variance_conversion_half_diamond_upper']) == 0
    assert F(neighbor['ideal_quadrature_half_diamond_upper']) > 0
