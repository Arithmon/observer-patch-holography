"""Original-input potential and slice controls at large/small binary64 scales."""
from pathlib import Path
import sys
import warnings

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_interacting_quantum as quantum


@pytest.mark.parametrize('power', [0, 500, 600])
def test_large_pure_gradient_is_never_a_coulomb_configuration(power):
    mesh = quantum.geometry()
    gauge = np.r_[0., 1., -1., np.zeros(10)]
    scale = 2.**power
    # Positive edge mass implies gauge.T D.T M D gauge > 0: a nonzero
    # exact gradient cannot lie in ker(D.T M), independently of units.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        with pytest.raises(ValueError, match='Coulomb slice'):
            quantum.reduced_kinetic(scale*(mesh.d@gauge), np.full(13, scale),
                                    charge=1/scale, mesh=mesh, include_potential=False)


@pytest.mark.parametrize('power', [100, 500, 600])
def test_neutral_density_omits_unused_covariant_products(power):
    scale = 2.**power
    coordinates = np.r_[scale, np.zeros(29), np.full(13, scale), np.zeros(13)]
    mp = mpmath.mp.clone()
    mp.dps = 80
    # At zero charge gamma is independent of a and psi. This absolute
    # determinant is independently derived from exact cone simplex masses
    # in test_whitney_kinetic_precision.original_geometry.
    logdet = mp.mpf('-59.98181822720102208420414571293429952844')
    expected = float(-14*(mp.log(2*mp.pi)+2*power*mp.log(2))-mp.mpf(14)/4-logdet/4)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        actual = quantum.gaussian_state_log_amplitude(coordinates, scale, charge=0)
    assert actual == pytest.approx(expected, rel=0, abs=2e-10)


@pytest.mark.parametrize('amplitude,mass,quartic', [
    (2.**600, 2.**-800, 0.), (2.**300, 0., 2.**-800),
    (2.**-600, 2.**800, 0.), (2.**-300, 0., 2.**800),
])
@pytest.mark.parametrize('evaluate', [quantum.coefficients, quantum.reduced_kinetic])
@pytest.mark.parametrize('phase', [1, 1j, 1+1j])
def test_balanced_original_potential_products_remain_representable(amplitude, mass, quartic, evaluate, phase):
    mp = mpmath.mp.clone()
    mp.dps = 80
    volume = 10+10*mp.sqrt(5)/3
    squared = mp.mpf(amplitude)**2*(mp.re(phase)**2+mp.im(phase)**2)
    expected = float(volume*(mp.mpf(mass)*squared+mp.mpf(quartic)*squared**2/2))
    result = evaluate(np.zeros(42), np.full(13, phase*amplitude), charge=0,
                      mass_squared=mass, quartic=quartic)
    actual = result[1] if evaluate is quantum.coefficients else result.potential
    assert actual > 0
    assert actual == pytest.approx(expected, rel=2e-12, abs=0)


@pytest.mark.parametrize('amplitude', [2.**-600, 2.**600])
def test_truly_unreportable_positive_potential_is_refused(amplitude):
    with pytest.raises(ValueError, match='precision|range'):
        quantum.coefficients(np.zeros(42), np.full(13, amplitude),
                             charge=0, mass_squared=1, quartic=0)


@pytest.mark.parametrize('delta', [0., 1.])
@pytest.mark.parametrize('amplitude', [1., 1e15])
def test_original_face_circulations_preserve_magnetic_energy(amplitude, delta):
    mesh = quantum.geometry()
    gauge = np.r_[0., 1., -1., np.zeros(10)]
    edges = amplitude*(mesh.d@gauge)
    edges[12] += delta
    # Every exact gradient has zero oriented face circulation. A unit change
    # on one boundary edge meets two golden tetrahedra. Its magnetic energy
    # is vol*(5-2sqrt(5)) = (5-sqrt(5))/6, independently of the large gauge.
    _, actual = quantum.coefficients(edges, np.zeros(13), charge=0,
                                     mass_squared=0, quartic=0, mesh=mesh)
    expected = delta**2*(5-np.sqrt(5))/6
    if delta:
        assert actual == pytest.approx(expected, rel=2e-12, abs=0)
    else:
        assert actual == 0
