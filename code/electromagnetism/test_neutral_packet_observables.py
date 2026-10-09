"""Original-input controls for neutral packet norms and squared radii."""
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_quantum_packet as packet
import verify_whitney_quantum_packet as replay


def perpendicular(center, momentum):
    q, p = np.zeros(56), np.zeros(56)
    q[30], p[43] = center, momentum
    return q, p


def perpendicular_integral(center, momentum, sigma=1., hbar=1.):
    """Positive integral, no Bessel functions or subtracting large moments.

    For this family z=|X-Y|. Set u=2 sqrt(z) sin(theta/2)
    when z is large. The Gaussian tail beyond u=16 is negligible at the
    asserted 1e-12 precision; these checks are not interval certificates.
    """
    mp = mpmath.mp.clone(); mp.dps = 90
    x, p, s, h = map(mp.mpf, (center, momentum, sigma, hbar))
    xx, pp = x*x/(4*s*s), s*s*p*p/(h*h)
    z = abs(xx-pp)
    if z > 128:
        def density(u):
            return mp.exp(-u*u/2)/mp.sqrt(1-u*u/(4*z))
        mass = mp.quad(density, [0, 1, 4, 8, 16])
        deficit = mp.quad(lambda u: u*u/2*density(u), [0, 1, 4, 8, 16])/mass
        scaled_norm = mass/(mp.pi*mp.sqrt(z))
    else:
        def density(theta):
            return mp.exp(-2*z*mp.sin(theta/2)**2)
        mass = mp.quad(density, [0, mp.pi/2, mp.pi])
        deficit = mp.quad(lambda t: 2*z*mp.sin(t/2)**2*density(t), [0, mp.pi/2, mp.pi])/mass
        scaled_norm = mass/mp.pi
    return {
        'A': xx+pp, 'B': x*p/h,
        'norm_squared': mp.exp(-2*min(xx, pp))*scaled_norm,
        'radius': 2*s*s*(13+2*max(xx-pp, 0)-deficit),
    }


@pytest.mark.parametrize('center,momentum', [
    (0., 0.), (1., .5), (1., .5+2**-40), (1., .5-2**-40),
    (0., 1e3), (0., 1e6), (0., 1e8), (0., 1e10), (0., 1e75),
    (0., 1e100), (1., 1e8), (-1., 1e8), (1., -1e8),
    (2., 1e10), (20., 1e8), (2e8, .5), (1e-150, 1e-150),
])
def test_original_input_observables_against_positive_integrals(center, momentum):
    q, p = perpendicular(center, momentum)
    expected = perpendicular_integral(center, momentum)
    actual = packet.overlap_parameters(q, p, 1)
    for key in ('A', 'B', 'norm_squared'):
        assert actual[key] == pytest.approx(float(expected[key]), rel=1e-12, abs=0), key
    assert packet.scalar_radius_moment(q, p, 1) == pytest.approx(float(expected['radius']), rel=1e-12, abs=0)


@pytest.mark.parametrize('sigma', [1e-100, 1e-50, .25, 1., 1e50, 1e100])
def test_position_units_preserve_dimensionless_norm_and_radius(sigma):
    q, p = perpendicular(sigma, 1e8/sigma)
    expected = perpendicular_integral(q[30], p[43], sigma)
    assert packet.overlap_parameters(q, p, sigma)['norm_squared'] == pytest.approx(float(expected['norm_squared']), rel=1e-12, abs=0)
    assert packet.scalar_radius_moment(q, p, sigma) == pytest.approx(float(expected['radius']), rel=1e-12, abs=0)


def test_radius_survives_unrepresentable_projection_probability():
    # The normalized state exists and its radius is representable even when
    # its strictly positive projection probability is below binary64 range.
    q, p = perpendicular(100., 1e8)
    expected = perpendicular_integral(100., 1e8)
    with pytest.raises(ValueError, match='represent|range|precision'):
        packet.overlap_parameters(q, p, 1)
    assert packet.scalar_radius_moment(q, p, 1) == pytest.approx(float(expected['radius']), rel=1e-12, abs=0)


@pytest.mark.parametrize('bad', [True, '1', 1j, np.nan, np.inf, np.ma.masked])
def test_mixed_invalid_original_scalar_rejected(bad):
    q, p = perpendicular(0., 1.)
    q = q.tolist(); q[30] = bad
    with pytest.raises(ValueError):
        packet.overlap_parameters(q, p, 1)


def test_masked_container_is_not_silently_unmasked():
    q, p = perpendicular(1., 2.)
    q = np.ma.array(q, mask=False); q.mask[30] = True
    with pytest.raises(ValueError):
        packet.scalar_radius_moment(q, p, 1)


@pytest.mark.parametrize('seed', range(10))
@pytest.mark.parametrize('sigma,hbar', [(.25, .5), (1., 2.), (4., .125)])
def test_general_scalar_directions_against_independent_integrals(seed, sigma, hbar):
    rng = np.random.default_rng(seed)
    q, p = rng.normal(size=(2, 56))
    expected = replay.circle_observables(q, p, sigma, hbar)
    actual = packet.overlap_parameters(q, p, sigma, hbar)
    for key in actual:
        assert actual[key] == pytest.approx(expected[key], rel=1e-12, abs=0), key
    assert packet.scalar_radius_moment(q, p, sigma, hbar) == pytest.approx(expected['scalar_radius_numeric'], rel=1e-12, abs=0)


def test_original_mixed_large_integers_preserve_small_charge():
    q, p = [0.]*56, [0.]*56
    q[30], q[31], p[43], p[44] = 2**53+1, 2**53, 1, -1
    assert packet.overlap_parameters(q, p, 1)['B'] == 1


@pytest.mark.parametrize('index', range(26))
def test_every_scalar_direction_retained(index):
    q, p = perpendicular(0., 0.)
    q[30+index], p[30+index] = 2, 3
    expected = replay.circle_observables(q, p, .5)
    assert packet.overlap_parameters(q, p, .5)['norm_squared'] == pytest.approx(expected['norm_squared'], rel=1e-12, abs=0)
    assert packet.scalar_radius_moment(q, p, .5) == pytest.approx(expected['scalar_radius_numeric'], rel=1e-12, abs=0)


@pytest.mark.parametrize('center', [38., 40., 100.])
def test_saturated_norm_refuses_insufficient_float_precision(center):
    q, p = perpendicular(center, center/2)
    with pytest.raises(ValueError, match='reporting range'):
        packet.overlap_parameters(q, p, 1)
    assert packet.scalar_radius_moment(q, p, 1) == 26
    # The complete normalized chiral state is a centered Gaussian even
    # when its separate projection norm cannot be reported as binary64.
    assert packet.projected_half_density(np.zeros(56), q, p, 1) == pytest.approx(
        (2*np.pi)**-14, rel=1e-12, abs=0)


@pytest.mark.parametrize('momentum', [1e160, 1e200, 1e300])
def test_radius_does_not_require_representable_A(momentum):
    q, p = perpendicular(0., momentum)
    with pytest.raises(ValueError, match='reporting range'):
        packet.overlap_parameters(q, p, 1)
    assert packet.scalar_radius_moment(q, p, 1) == 25


@pytest.mark.parametrize('sigma', [1e-170, 1e170])
def test_unrepresentable_radius_is_not_zero_or_infinity(sigma):
    q, p = perpendicular(0., 0.)
    with pytest.raises(ValueError, match='reporting range'):
        packet.scalar_radius_moment(q, p, sigma)
    # A well-defined probability is independent of the radius reporting range.
    assert packet.overlap_parameters(q, p, sigma)['norm_squared'] == 1


def test_small_positive_observables_cannot_be_zeroed():
    for value in (1e-250, 1e-100, 1e-20):
        replay.close_observable(value, value, 'retained')
        with pytest.raises(ValueError):
            replay.close_observable(0., value, 'erased')
    replay.close_observable(0, 0., 'true zero')
    with pytest.raises(ValueError):
        replay.close_observable(1e-250, 0., 'invented')


def test_private_mpmath_context_and_uncharged_sign_symmetry():
    old = mpmath.mp.dps
    q, p = perpendicular(1., 1e10)
    first = packet.overlap_parameters(q, p, 1)
    second = packet.overlap_parameters(q, -p, 1)
    assert first['A'] == second['A'] and first['B'] == -second['B']
    assert first['norm_squared'] == second['norm_squared']
    assert packet.scalar_radius_moment(q, p, 1) == packet.scalar_radius_moment(q, -p, 1)
    assert mpmath.mp.dps == old


@pytest.mark.parametrize('bad', [True, np.ma.masked])
@pytest.mark.parametrize('reader', ['seed', 'projected', 'rotate', 'phase_input'])
def test_float_caller_rejects_original_invalid_entries(bad, reader):
    value = [0.]*56; value[0] = bad
    q, p = np.zeros(56), np.zeros(56)
    with pytest.raises(ValueError):
        if reader == 'seed':
            packet.seed_log_half_density(value, q, p, 1)
        elif reader == 'projected':
            packet.projected_half_density(value, q, p, 1, nodes=16)
        elif reader == 'rotate':
            packet.rotate(value, .5)
        else:
            packet.real_vector(value+[0.]*12, 68, 'phase configuration')


@pytest.mark.parametrize('reader', ['seed', 'projected', 'rotate', 'phase_input'])
def test_float_caller_rejects_masked_container(reader):
    value = np.ma.array(np.zeros(56), mask=False); value.mask[0] = True
    q, p = np.zeros(56), np.zeros(56)
    with pytest.raises(ValueError):
        if reader == 'seed':
            packet.seed_log_half_density(value, q, p, 1)
        elif reader == 'projected':
            packet.projected_half_density(value, q, p, 1, nodes=16)
        elif reader == 'rotate':
            packet.rotate(value, .5)
        else:
            packet.real_vector(value, 56, 'phase configuration')


def test_float_caller_refuses_lost_displacement_and_accepts_resolved_neighbor():
    q, p, point = [0.]*56, [0.]*56, [0.]*56
    q[0], point[0] = 2**53, 2**53+1
    with pytest.raises(ValueError, match='binary64'):
        packet.seed_log_half_density(point, q, p, 1)
    point[0] = 2**53+2
    assert packet.seed_log_half_density(point, q, p, 1)-packet.seed_log_half_density(q, q, p, 1) == -1


@pytest.mark.parametrize('bad', [True, np.bool_(False), np.ma.masked])
def test_float_caller_rejects_invalid_angle(bad):
    with pytest.raises(ValueError):
        packet.rotate(np.zeros(56), bad)


@pytest.mark.parametrize('parameter', ['width', 'hbar'])
def test_float_caller_refuses_changed_exact_parameter(parameter):
    from fractions import Fraction
    q, p = np.zeros(56), np.zeros(56)
    width, hbar = (Fraction(1, 3), 1) if parameter == 'width' else (1, Fraction(1, 3))
    # The exact scalar calculation continues to support the original value.
    assert packet.overlap_parameters(q, p, width, hbar)['norm_squared'] == 1
    with pytest.raises(ValueError, match='binary64'):
        packet.seed_log_half_density(q, q, p, width, hbar)


def power_series_reference(q, p, sigma, hbar):
    """Decimal positive series from the Gaussian overlap's Fourier powers.

    Independently sum both its normalization and k-weighted moment; no
    mpmath, Bessel functions, quadrature or producer intermediate values.
    """
    from decimal import Decimal, localcontext
    from fractions import Fraction as F
    x, p = list(map(F, q[30:])), list(map(F, p[30:]))
    s, h = F(sigma), F(hbar)
    xx, pp = sum(v*v for v in x), sum(v*v for v in p)
    jx = [-v for v in x[13:]]+x[:13]
    b = sum(u*v for u, v in zip(p, jx))/h
    a = xx/(4*s*s)+s*s*pp/(h*h)
    t = (a*a-b*b)/4
    assert 0 <= t <= 256**2/4  # bounded reference, no large-z approximation
    with localcontext() as ctx:
        ctx.prec = 350
        def d(v):
            return Decimal(v.numerator)/Decimal(v.denominator)
        total, weighted, term = Decimal(1), Decimal(0), Decimal(1)
        for k in range(1, 10000):
            term *= d(t)/Decimal(k*k)
            total += term
            weighted += k*term
            ratio = d(t)/Decimal((k+1)**2)
            if ratio < Decimal('.5') and (k+1)*term < total*Decimal('1e-120'):
                break
        else:
            raise AssertionError('reference did not converge')
        norm = (-d(a)).exp()*total
        radius = d(26*s*s+xx/2-2*s**4*pp/h**2)+4*d(s*s)*weighted/total
        return {'A': d(a), 'B': d(b), 'norm_squared': norm, 'radius': radius}


@pytest.mark.parametrize('case', range(24))
def test_independent_positive_power_series(case):
    from fractions import Fraction as F
    q, p = [0.]*56, [0.]*56
    sigma, hbar = 1, 1
    if case < 8:
        q[30], p[43] = 2**case, 2**(case-1)
    elif case < 12:
        # The exact discriminant is nonzero even when rounded A equals B.
        q[30], p[43] = F(2), F(1)+F((-1)**case, 2**(40+case))
    elif case < 16:
        q[30], p[43], p[44] = F(2), F(1), F(1, 2**(case+12))
        sigma = F(1, 2)**(case-13)
        hbar = F(1, 2)**(2*(case-13))
    else:
        rng = np.random.default_rng(case)
        q, p = [v/8 for v in rng.integers(-8,9,56)], [v/8 for v in rng.integers(-8,9,56)]
    expected = power_series_reference(q, p, sigma, hbar)
    probability = float(expected['norm_squared'])
    if probability == 0:
        with pytest.raises(ValueError, match='reporting range'):
            packet.overlap_parameters(q, p, sigma, hbar)
    else:
        result = packet.overlap_parameters(q, p, sigma, hbar)
        for key in result:
            assert result[key] == pytest.approx(float(expected[key]), rel=1e-12, abs=0), key
    assert packet.scalar_radius_moment(q, p, sigma, hbar) == pytest.approx(float(expected['radius']), rel=1e-12, abs=0)


@pytest.mark.parametrize('sigma', [1e-200, 1e200, 1e-20])
def test_pointwise_log_survives_unrepresentable_width_square(sigma):
    q = np.zeros(56)
    mp = mpmath.mp.clone(); mp.dps = 90
    expected = -14*mp.log(2*mp.pi*mp.mpf(sigma)**2)
    assert packet.seed_log_half_density(q, q, q, sigma) == pytest.approx(complex(float(expected)), rel=1e-14, abs=0)


@pytest.mark.parametrize('sigma', [1e-200, 1e200, 1e-20])
def test_unreportable_pointwise_amplitude_is_an_explicit_error(sigma):
    q = np.zeros(56)
    with pytest.raises(ValueError, match='range'):
        packet.projected_half_density(q, q, q, sigma, nodes=16)


def test_representable_pointwise_mean_does_not_overflow_its_sum():
    q = np.zeros(56); sigma = 10**-11.35
    mp = mpmath.mp.clone(); mp.dps = 90
    expected = float((2*mp.pi*mp.mpf(sigma)**2)**-14)
    assert packet.projected_half_density(q, q, q, sigma) == pytest.approx(expected, rel=2e-12, abs=0)


@pytest.mark.parametrize('value', [
    np.int8(3), np.uint8(255), np.int16(-300), np.uint16(65535),
    np.int32(2**30), np.uint32(2**32-1), np.int64(-(2**63)),
    np.int64(2**32), np.uint64(2**63+1),
])
def test_numpy_integer_coordinates_use_unbounded_rational_arithmetic(value):
    q, p = [0.]*56, [0.]*56; q[30] = value
    expected = perpendicular_integral(int(value), 0)
    for key, actual in packet.overlap_parameters(q, p, 1).items():
        assert actual == pytest.approx(float(expected[key]), rel=1e-12, abs=0)
    assert packet.scalar_radius_moment(q, p, 1) == pytest.approx(float(expected['radius']), rel=1e-12, abs=0)


@pytest.mark.parametrize('wrap_fraction', [False, True])
@pytest.mark.parametrize('field', ['center', 'momentum', 'width', 'hbar'])
def test_numpy_integer_parameters_and_fraction_components(field, wrap_fraction):
    from fractions import Fraction as F
    value = F(np.int64(3), np.int64(2)) if wrap_fraction else np.int64(3)
    native = F(3,2) if wrap_fraction else 3
    q, p, sigma, hbar = [0.]*56, [0.]*56, 1, 1
    q[30], p[43] = 2, 1
    if field == 'center': q[30] = value
    elif field == 'momentum': p[43] = value
    elif field == 'width': sigma = value
    else: hbar = value
    # Build the independent series from native Python integer components.
    qq, pp, ss, hh = q.copy(), p.copy(), sigma, hbar
    if field == 'center': qq[30] = native
    elif field == 'momentum': pp[43] = native
    elif field == 'width': ss = native
    else: hh = native
    expected = power_series_reference(qq, pp, ss, hh)
    for key, actual in packet.overlap_parameters(q, p, sigma, hbar).items():
        assert actual == pytest.approx(float(expected[key]), rel=1e-12, abs=0)
    assert packet.scalar_radius_moment(q, p, sigma, hbar) == pytest.approx(float(expected['radius']), rel=1e-12, abs=0)
