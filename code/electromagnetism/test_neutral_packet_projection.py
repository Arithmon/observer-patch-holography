"""Original-input Gaussian controls for the neutral pointwise packet."""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_quantum_packet as packet


def _original_real(ctx, value):
    value = Fraction(value)
    return ctx.mpf(value.numerator)/value.denominator


@lru_cache(maxsize=None)
def _oscillatory_packet_oracle(momentum):
    """Original circle Gaussian integral, with positive overlap integration.

    For a zero center and one scalar momentum P, the pointwise circle
    integral at unit scalar radius is J0(P). For its squared norm change
    variables u=2P sin(theta/2) in the positive overlap integral. The tail
    beyond u=16 is negligible at the asserted tolerance for P>=10.
    Smaller P uses the untransformed positive overlap integral. Neither
    this check nor the producer's precision agreement is an interval bound.
    """
    ctx = mpmath.mp.clone()
    ctx.dps = 100+max(0, Fraction(momentum).numerator.bit_length()//3)
    p = _original_real(ctx, momentum)
    if p >= 10:
        density = lambda u: ctx.exp(-u*u/2)/ctx.sqrt(1-u*u/(4*p*p))
        norm_squared = ctx.quad(density, [0, 1, 4, 8, 16])/(ctx.pi*p)
    else:
        density = lambda theta: ctx.exp(-2*p*p*ctx.sin(theta/2)**2)
        norm_squared = ctx.quad(density, [0, ctx.pi/2, ctx.pi])/ctx.pi
    return (2*ctx.pi)**-14*ctx.exp(-ctx.mpf(1)/4)*ctx.besselj(0, p)/ctx.sqrt(norm_squared)


def _cancelling_phase_inputs(indices=(0, 4, 8)):
    point, center, momentum = [0]*56, [0]*56, [0]*56
    for index, value in zip(indices, (2**53, 1, -(2**53))):
        point[index], momentum[index] = 1, value
    assert sum(Fraction(p)*Fraction(x) for p, x in zip(momentum, point)) == 1
    return point, center, momentum


def test_seed_keeps_exact_phase_when_large_dot_products_cancel():
    point, center, momentum = _cancelling_phase_inputs()
    measured = packet.seed_log_half_density(point, center, momentum, 1)
    assert measured.imag == 1.0


def test_public_projection_keeps_the_same_cancelled_radiative_phase():
    point, center, momentum = _cancelling_phase_inputs()
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    expected = (2*ctx.pi)**-14*ctx.exp(-ctx.mpf(3)/4 + ctx.j)
    measured = packet.projected_half_density(point, center, momentum, 1)
    assert measured == pytest.approx(complex(expected), rel=2e-12, abs=0)


def test_original_displacement_is_retained_before_phase_cancellation():
    point, center, momentum = [0]*56, [0]*56, [0]*56
    point[0] = point[1] = 1
    center[0] = 2.**-54
    momentum[0], momentum[1] = 2**54, -(2**54)
    # Each supplied value is exactly binary64, but their first displacement
    # is not. Rounding it before the otherwise exact dot product loses a radian.
    assert point[0]-center[0] == 1
    delta = Fraction(point[0])-Fraction(center[0])
    phase = Fraction(momentum[0])*delta+Fraction(momentum[1])*point[1]
    assert phase == -1
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    expected_log = (-14*ctx.log(2*ctx.pi)
                    - (_original_real(ctx, delta)**2+1)/4
                    + ctx.j*_original_real(ctx, phase))
    measured = packet.seed_log_half_density(point, center, momentum, 1)
    assert measured.imag == -1
    assert measured.real == pytest.approx(float(expected_log.real), rel=2e-12, abs=0)
    assert packet.projected_half_density(point, center, momentum, 1) == pytest.approx(
        complex(ctx.exp(expected_log)), rel=2e-12, abs=0)


@pytest.mark.parametrize("momentum", [10., 256., 1e20, 1e100])
def test_public_projection_resolves_oscillatory_scalar_momentum(momentum):
    center, p, point = np.zeros(56), np.zeros(56), np.zeros(56)
    p[30], point[30] = momentum, 1.
    expected = complex(_oscillatory_packet_oracle(momentum))
    measured = packet.projected_half_density(point, center, p, 1, nodes=256)
    assert measured == pytest.approx(expected, rel=2e-12, abs=0)


def test_oscillatory_neutral_projection_is_rotation_invariant():
    center, p, point = np.zeros(56), np.zeros(56), np.zeros(56)
    p[30], point[30] = 256., 1.
    expected = complex(_oscillatory_packet_oracle(256.))
    angle = .373
    rotated = point.copy()
    rotated[30], rotated[43] = np.cos(angle), np.sin(angle)
    for evaluation_point in (point, rotated):
        measured = packet.projected_half_density(evaluation_point, center, p, 1, nodes=256)
        assert measured == pytest.approx(expected, rel=3e-12, abs=0)


def _general_inputs():
    point, center, momentum = [0]*56, [0]*56, [0]*56
    for index, y, q, p in ((0, -.5, .25, .75), (17, .25, -.125, .5),
            (30, .5, .75, -.25), (43, -.25, .125, .5),
            (33, -.375, .25, .75), (46, .625, -.5, -.125),
            (42, .25, -.125, .375), (55, -.5, .25, .125)):
        point[index], center[index], momentum[index] = y, q, p
    return point, center, momentum


def _original_circle_oracle(point, center, momentum, sigma, hbar):
    """Integrate the original rotated Gaussian and its Gaussian overlap.

    The pointwise numerator is an angular integral. The denominator uses
    the original pairwise Gaussian overlap, obtained by completing the
    square in all 56 coordinates. Neither integral uses Bessel invariants
    or any producer helper; these moderate data resolve by direct quadrature.
    """
    ctx = mpmath.mp.clone()
    ctx.dps = 70
    y, q, p = [[_original_real(ctx, v) for v in values]
               for values in (point, center, momentum)]
    s, h = _original_real(ctx, sigma), _original_real(ctx, hbar)
    pairs = [i for i in range(13)
             if any(values[30+i] or values[43+i] for values in (y, q, p))]
    active = [i for i in range(30) if y[i] or q[i] or p[i]]
    # Rotation creates the other coordinate even if its supplied value is
    # zero. Include the whole scalar plane in both Gaussian integrals.
    active += [index for i in pairs for index in (30+i, 43+i)]

    def rotated(values, cosine, sine):
        result = values.copy()
        for i in pairs:
            result[30+i] = cosine*values[30+i]-sine*values[43+i]
            result[43+i] = sine*values[30+i]+cosine*values[43+i]
        return result

    def numerator(theta):
        cosine, sine = ctx.cos(theta), ctx.sin(theta)
        qt, pt = rotated(q, cosine, sine), rotated(p, cosine, sine)
        displacement = {i: y[i]-qt[i] for i in active}
        decay = sum(v*v for v in displacement.values())/(4*s*s)
        phase = sum(pt[i]*displacement[i] for i in active)/h
        return ctx.exp(-decay+ctx.j*phase)

    def overlap(theta):
        cosine, sine = ctx.cos(theta), ctx.sin(theta)
        qt, pt = rotated(q, cosine, sine), rotated(p, cosine, sine)
        dq = {i: q[i]-qt[i] for i in active}
        dp = {i: p[i]-pt[i] for i in active}
        decay = (sum(v*v for v in dq.values())/(8*s*s)
                 + s*s*sum(v*v for v in dp.values())/(2*h*h))
        phase = sum((p[i]+pt[i])*dq[i] for i in active)/(2*h)
        return ctx.exp(-decay+ctx.j*phase)

    intervals = [0, ctx.pi/2, ctx.pi, 3*ctx.pi/2, 2*ctx.pi]
    norm_squared = ctx.quad(overlap, intervals)/(2*ctx.pi)
    assert norm_squared.real > 0
    assert abs(norm_squared.imag) < ctx.mpf('1e-60')*norm_squared.real
    prefactor = (2*ctx.pi*s*s)**-14
    return prefactor*ctx.quad(numerator, intervals)/(2*ctx.pi*ctx.sqrt(norm_squared.real))


@pytest.mark.parametrize("sigma,hbar", [(.5, .75), (2., .5)])
def test_general_complex_projection_matches_original_gaussian_integrals(sigma, hbar):
    point, center, momentum = _general_inputs()
    expected = complex(_original_circle_oracle(point, center, momentum, sigma, hbar))
    # Both components matter: a modulus-only computation loses coherent phase.
    assert abs(expected.real) > 0 and abs(expected.imag) > 0
    actual = packet.projected_half_density(point, center, momentum, sigma, hbar)
    assert actual.real == pytest.approx(expected.real, rel=2e-12, abs=0)
    assert actual.imag == pytest.approx(expected.imag, rel=2e-12, abs=0)
    # Permute radiative coordinates and scalar complex modes together.
    permutation = list(reversed(range(30)))
    scalar = list(reversed(range(13)))
    permutation += [30+i for i in scalar]+[43+i for i in scalar]
    permuted = [[values[i] for i in permutation] for values in (point, center, momentum)]
    assert packet.projected_half_density(*permuted, sigma, hbar) == pytest.approx(expected, rel=2e-12, abs=0)


@pytest.mark.parametrize("indices", [(0, 1, 2), (0, 8, 16), (0, 1, 29)])
def test_cancelled_phase_is_independent_of_radiative_coordinate_order(indices):
    point, center, momentum = _cancelling_phase_inputs(indices)
    assert packet.seed_log_half_density(point, center, momentum, 1).imag == 1
    ctx = mpmath.mp.clone()
    ctx.dps = 80
    expected = (2*ctx.pi)**-14*ctx.exp(-ctx.mpf(3)/4+ctx.j)
    assert packet.projected_half_density(point, center, momentum, 1) == pytest.approx(complex(expected), rel=2e-12, abs=0)


def test_original_integral_keeps_rotation_generated_coordinates():
    point, q, p = [0]*56, [0]*56, [0]*56
    point[0], point[30], q[30], p[30] = .25, .375, 1, .75
    expected = complex(_original_circle_oracle(point, q, p, .5, 1))
    actual = packet.projected_half_density(point, q, p, .5)
    assert actual.real == pytest.approx(expected.real, rel=2e-12, abs=0)
    assert actual.imag == pytest.approx(expected.imag, rel=2e-12, abs=0)


def test_projection_preserves_destructive_coherent_interference():
    point, center, momentum = np.zeros(56), np.zeros(56), np.zeros(56)
    point[30], momentum[30] = 1, 256
    oscillatory = packet.projected_half_density(point, center, momentum, 1)
    centered = packet.projected_half_density(point, center, np.zeros(56), 1)
    ctx = mpmath.mp.clone()
    ctx.dps = 80
    reference = (2*ctx.pi)**-14*ctx.exp(-ctx.mpf(1)/4)
    expected = complex(reference+_oscillatory_packet_oracle(256.))
    assert oscillatory.real < 0 < centered.real
    assert abs(centered+oscillatory)**2 == pytest.approx(abs(expected)**2, rel=2e-12, abs=0)


@pytest.mark.parametrize("center,sigma,hbar", [(100., 1., 1.), (1e100, 1., 1.), (1e100, .25, .5)])
def test_exact_chiral_normalization_compensates_unreportable_norm(center, sigma, hbar):
    q, p, point = [0]*56, [0]*56, [0]*56
    q[30], p[43] = center, hbar*center/(2*sigma*sigma)
    # p=(hbar/(2 sigma^2)) Jq exactly; P0 removes every nonconstant
    # holomorphic term, so the normalized result is a centered Gaussian.
    assert Fraction(p[43]) == Fraction(hbar)*Fraction(center)/(2*Fraction(sigma)**2)
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    s = _original_real(ctx, sigma)
    for radius in (0., .375):
        point[30], point[44] = radius, -radius
        expected = (2*ctx.pi*s*s)**-14*ctx.exp(-2*_original_real(ctx, radius)**2/(4*s*s))
        assert packet.projected_half_density(point, q, p, sigma, hbar) == pytest.approx(complex(expected), rel=2e-12, abs=0)


@pytest.mark.parametrize("exponent", [500, 1023])
@pytest.mark.parametrize("direction", [-1, 1])
def test_near_chiral_complex_projection_matches_original_fourier_series(exponent, direction):
    point, q, p = [0]*56, [0]*56, [0]*56
    center = 2.**exponent
    point[30], q[30], p[30], p[43] = 1, center, direction/center, center/2
    assert Fraction(p[30])*Fraction(center) == direction
    ctx = mpmath.mp.clone()
    ctx.dps = 180
    inverse = _original_real(ctx, p[30])

    def constant_fourier(product):
        # In exp(a*exp(i theta)+b*exp(-i theta)), precisely equal
        # Fourier powers survive averaging. Sum them without Bessel calls.
        term = total = ctx.mpc(1)
        for k in range(1, 1000):
            term *= product/(k*k)
            total += term
            if abs(term) < abs(total)*ctx.mpf('1e-150'):
                return total
        raise AssertionError('original Fourier series did not resolve')

    # Original one-plane Gaussian: q=(C,0), p=(direction/C,C/2), y=(1,0).
    # Its two Fourier coefficients multiply to (i*direction-C^-2)/4;
    # the Gaussian overlap coefficients multiply to (1+C^-4)/4.
    # These simplified products retain the near-chiral perturbation without
    # subtracting the huge center terms or reusing producer intermediates.
    numerator = constant_fourier((direction*ctx.j-inverse**2)/4)
    norm = constant_fourier((1+inverse**4)/4)
    expected = ((2*ctx.pi)**-14
                * ctx.exp(-ctx.mpf(1)/4+inverse**2/2-direction*ctx.j)
                * numerator/ctx.sqrt(norm))
    actual = packet.projected_half_density(point, q, p, 1)
    assert actual.real == pytest.approx(float(expected.real), rel=2e-12, abs=0)
    assert actual.imag == pytest.approx(float(expected.imag), rel=2e-12, abs=0)


@pytest.mark.parametrize("center", [1e100, 1e308])
def test_large_pure_position_projection_resolves_normalized_amplitude(center):
    q, p = [0]*56, [0]*56
    q[30] = center
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    # Laplace expansion of the two original positive angular integrals:
    # exp(-x) I0(x) ~ 1/sqrt(2 pi x). Relative corrections here are
    # smaller than 1e-190, far below the asserted binary64 tolerance.
    expected = (2*ctx.pi)**(-ctx.mpf(57)/4)/ctx.sqrt(_original_real(ctx, center))
    assert packet.projected_half_density(q, q, p, 1) == pytest.approx(complex(expected), rel=2e-12, abs=0)


def test_nearest_binary64_momenta_on_both_sides_of_a_bessel_zero():
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    nearest = float(ctx.besseljzero(0, 1))
    values = [np.nextafter(nearest, -np.inf), nearest, np.nextafter(nearest, np.inf)]
    actual = []
    for value in values:
        point, q, p = np.zeros(56), np.zeros(56), np.zeros(56)
        point[30], p[30] = 1, value
        expected = complex(_oscillatory_packet_oracle(float(value)))
        assert expected != 0
        result = packet.projected_half_density(point, q, p, 1)
        assert result == pytest.approx(expected, rel=2e-12, abs=0)
        actual.append(result.real)
    assert actual[0] > 0 > actual[-1]


def _large_radiative_phase_oracle():
    ctx = mpmath.mp.clone()
    ctx.dps = 430
    phase = _original_real(ctx, 1e300)/7
    return complex((2*ctx.pi)**-14*ctx.exp(-ctx.mpf(1)/4+ctx.j*phase))


def test_huge_radiative_phase_uses_original_ratio_before_phase_reduction():
    point, q, p = [0]*56, [0]*56, [0]*56
    point[0], p[0] = 1, 1e300
    expected = _large_radiative_phase_oracle()
    assert Fraction(float(Fraction(1e300)/7)) != Fraction(1e300)/7
    assert packet.projected_half_density(point, q, p, 1, 7) == pytest.approx(expected, rel=2e-12, abs=0)
    with pytest.raises(ValueError, match="phase.*(precision|range)"):
        packet.seed_log_half_density(point, q, p, 1, 7)
    # A huge but exactly representable unwrapped phase remains a valid log.
    exact_phase = Fraction(1e300)/3
    assert Fraction(float(exact_phase)) == exact_phase
    assert packet.seed_log_half_density(point, q, p, 1, 3).imag == float(exact_phase)


def test_maximal_momentum_over_minimal_hbar_preserves_reduced_phase():
    point, q, p = [0]*56, [0]*56, [0]*56
    point[0], p[0] = 1, float.fromhex('0x1.fffffffffffffp+1023')
    hbar = float.fromhex('0x0.0000000000001p-1022')
    phase = Fraction(p[0])/Fraction(hbar)
    assert phase > Fraction(p[0])
    ctx = mpmath.mp.clone()
    ctx.dps = 800  # Original phase has 632 decimal digits before reduction.
    expected = (2*ctx.pi)**-14*ctx.exp(-ctx.mpf(1)/4+ctx.j*_original_real(ctx, phase))
    with pytest.raises(ValueError, match='phase.*range'):
        packet.seed_log_half_density(point, q, p, 1, hbar)
    assert packet.projected_half_density(point, q, p, 1, hbar) == pytest.approx(
        complex(expected), rel=2e-12, abs=0)


def test_projection_is_independent_of_inherited_mpmath_precision():
    before = (mpmath.mp.dps, mpmath.mp.trap_complex)
    try:
        mpmath.mp.dps, mpmath.mp.trap_complex = 3, True
        point, q, p = [0]*56, [0]*56, [0]*56
        point[0], p[0] = 1, 1e300
        assert packet.projected_half_density(point, q, p, 1, 7) == pytest.approx(_large_radiative_phase_oracle(), rel=2e-12, abs=0)
        assert (mpmath.mp.dps, mpmath.mp.trap_complex) == (3, True)
    finally:
        mpmath.mp.dps, mpmath.mp.trap_complex = before


@pytest.mark.parametrize("sigma", [1e-8, .5, 1., 1e8])
def test_seed_prefactor_normalizes_all_56_original_gaussian_coordinates(sigma):
    zero = [0]*56
    value = packet.seed_log_half_density(zero, zero, zero, sigma)
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    s = _original_real(ctx, sigma)
    # Integral of |seed|^2 from 56 separate one-dimensional Gaussians.
    integrated_norm = ctx.exp(2*ctx.mpf(value.real))*(2*ctx.pi*s*s)**28
    assert float(integrated_norm) == pytest.approx(1., rel=2e-12, abs=0)


def test_seed_and_projection_keep_representable_result_when_width_squared_underflows():
    sigma = 2.**-550
    assert sigma*sigma == 0
    point, q, p = [0]*56, [0]*56, [0]*56
    point[0] = 206*sigma
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    expected_log = -14*ctx.log(2*ctx.pi*_original_real(ctx, sigma)**2)-ctx.mpf(206)**2/4
    assert packet.seed_log_half_density(point, q, p, sigma) == pytest.approx(complex(expected_log), rel=2e-12, abs=0)
    assert packet.projected_half_density(point, q, p, sigma) == pytest.approx(complex(ctx.exp(expected_log)), rel=2e-12, abs=0)


def test_seed_phase_reports_resolved_subnormal_and_refuses_unrepresentable_phase():
    point, q, p = [0]*56, [0]*56, [0]*56
    point[0], p[0] = 1, 1e-320
    assert packet.seed_log_half_density(point, q, p, 1).imag == 1e-320
    point[0], p[0] = .5, float.fromhex('0x0.0000000000001p-1022')
    with pytest.raises(ValueError, match="phase.*precision"):
        packet.seed_log_half_density(point, q, p, 1)


@pytest.mark.parametrize("radius,accepted", [(52.5, True), (53., False), (55., False)])
def test_pointwise_output_range_depends_on_final_amplitude_precision(radius, accepted):
    point, q, p = [0]*56, [0]*56, [0]*56
    point[0] = radius
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    expected = (2*ctx.pi)**-14*ctx.exp(-_original_real(ctx, radius)**2/4)
    loss = abs(ctx.mpf(float(expected))/expected-1)
    assert (loss <= ctx.mpf('1e-12')) == accepted
    if accepted:
        assert 0 < float(expected) < np.finfo(float).tiny
        assert packet.projected_half_density(point, q, p, 1) == pytest.approx(complex(expected), rel=1e-12, abs=0)
    else:
        with pytest.raises(ValueError, match="reporting range"):
            packet.projected_half_density(point, q, p, 1)


@pytest.mark.parametrize("nodes", [16, 32, 511, 4096])
def test_legacy_node_count_does_not_change_the_circle_integral(nodes):
    point, q, p = [0]*56, [0]*56, [0]*56
    point[30], p[30] = 1, 256
    assert packet.projected_half_density(point, q, p, 1, nodes=nodes) == pytest.approx(complex(_oscillatory_packet_oracle(256.)), rel=2e-12, abs=0)


@pytest.mark.parametrize("nodes", [True, 15, 16., "256"])
def test_invalid_legacy_node_count_is_not_silently_coerced(nodes):
    with pytest.raises(ValueError, match="integer circle nodes"):
        packet.projected_half_density([0]*56, [0]*56, [0]*56, 1, nodes=nodes)


@pytest.mark.parametrize("function", [packet.seed_log_half_density, packet.projected_half_density])
@pytest.mark.parametrize("position", [0, 1, 2, 3, 4])
@pytest.mark.parametrize("invalid", [Fraction(1, 3), 2**53+1])
def test_public_pointwise_calls_reject_silent_original_scalar_rounding(function, position, invalid):
    args = [[0]*56, [0]*56, [0]*56, 1, 1]
    if position < 3:
        args[position][0] = invalid
    else:
        args[position] = invalid
    with pytest.raises(ValueError, match="binary64"):
        function(*args)


@pytest.mark.parametrize("function", [packet.seed_log_half_density, packet.projected_half_density])
def test_exact_dyadic_rational_inputs_remain_valid(function):
    point, q, p = [0]*56, [0]*56, [0]*56
    point[0], q[0], p[0] = Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)
    result = function(point, q, p, Fraction(1, 2), Fraction(3, 4))
    expected = function([float(v) for v in point], [float(v) for v in q], [float(v) for v in p], .5, .75)
    assert result == expected


@pytest.mark.parametrize("invalid", [True, "1", 1+0j, np.nan, np.inf, np.ma.array(1., mask=False)])
@pytest.mark.parametrize("function", [packet.seed_log_half_density, packet.projected_half_density])
def test_pointwise_scalar_validation_precedes_array_coercion(function, invalid):
    point = [0]*56
    point[7] = invalid
    with pytest.raises(ValueError):
        function(point, [0]*56, [0]*56, 1)


@pytest.mark.parametrize("invalid", [[0]*55, [[0]]*56, np.ma.array(np.zeros(56), mask=False)])
def test_projection_rejects_invalid_original_vector_shape_or_mask(invalid):
    with pytest.raises(ValueError):
        packet.projected_half_density(invalid, [0]*56, [0]*56, 1)


@pytest.mark.parametrize("sigma,hbar", [(0, 1), (-1, 1), (1, 0), (1, -1)])
def test_projection_requires_positive_width_and_hbar(sigma, hbar):
    with pytest.raises(ValueError):
        packet.projected_half_density([0]*56, [0]*56, [0]*56, sigma, hbar)
