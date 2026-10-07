"""Original-input controls for neutral packet norms and squared radii."""
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_quantum_packet as packet


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
