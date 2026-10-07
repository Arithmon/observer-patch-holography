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
    with pytest.raises(ValueError, match='reporting range'):
        packet.projected_half_density(q, q, p, 1)


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
