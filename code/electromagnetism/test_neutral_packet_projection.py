"""Original-input Gaussian controls for the neutral pointwise packet."""
from fractions import Fraction
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_quantum_packet as packet


def _oscillatory_packet_oracle(momentum):
    """Original circle Gaussian integral, with positive overlap integration.

    For a zero center and one scalar momentum P, the pointwise circle
    integral at unit scalar radius is J0(P). For its squared norm change
    variables u=2P sin(theta/2) in the positive overlap integral. The tail
    beyond u=16 is negligible at the asserted tolerance for these P>=10
    controls; this independent check is not an interval certificate.
    """
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    p = ctx.mpf(momentum)
    density = lambda u: ctx.exp(-u*u/2)/ctx.sqrt(1-u*u/(4*p*p))
    norm_squared = ctx.quad(density, [0, 1, 4, 8, 16])/(ctx.pi*p)
    return (2*ctx.pi)**-14*ctx.exp(-ctx.mpf(1)/4)*ctx.besselj(0, p)/ctx.sqrt(norm_squared)


def _cancelling_phase_inputs():
    point, center, momentum = [0]*56, [0]*56, [0]*56
    for index, value in zip((0, 4, 8), (2**53, 1, -(2**53))):
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


@pytest.mark.parametrize("momentum", [10., 256.])
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
