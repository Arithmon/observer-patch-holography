"""Independent source-validation and arbitrary-rank classical fidelity controls."""

from fractions import Fraction

import mpmath
import numpy as np
import pytest

from quantum_information.recovery import state_fidelity
from quantum_information.state_factors import source_state
from quantum_information.states import ATOL


def test_original_subnormal_source_is_independent_of_numpy_error_mode():
    unit = np.nextafter(0., 1.)
    rho = np.diag([1., unit, unit]).astype(complex)
    rho[1, 2] = 1j*unit
    sigma = np.zeros((3, 3), complex)
    sigma[1:, 1:] = [[.5, .5j], [-.5j, .5]]
    # Exact averaging retains coherence i*u/2. The pure sigma gives
    # F=Tr(rho sigma)=3u/2, which rounds ties-to-even to 2u.
    previous = np.geterr()
    with np.errstate(all="raise"):
        assert source_state(rho).rank == 3
        assert state_fidelity(rho, sigma) == 2*unit
        assert np.geterr() == dict.fromkeys(previous, "raise")
    assert np.geterr() == previous


def test_invalid_finite_huge_entries_are_rejected_without_floating_overflow():
    invalid = np.array([[.5, 1e300], [0., .5]])
    with np.errstate(all="raise"):
        with pytest.raises(ValueError, match="Hermitian"):
            source_state(invalid)


@pytest.mark.parametrize("drift,valid", ((2.**-42, True), (2.**-39, False)))
def test_original_trace_tolerance_is_preserved_without_normalization(drift, valid):
    rho = np.diag([.5+drift, .5])
    if valid:
        expected = float((Fraction(float(rho[0, 0]))+Fraction(.5))**2)
        assert state_fidelity(rho, rho) == expected
    else:
        with pytest.raises(ValueError, match="trace one"):
            state_fidelity(rho, rho)


@pytest.mark.parametrize("multiplier,valid", ((.25, True), (.5, False)))
def test_original_complex_hermiticity_tolerance(multiplier, valid):
    rho = np.diag([.5+1j*ATOL*multiplier, .5-1j*ATOL*multiplier])
    # ||rho-rho*||_F = 2 sqrt(2) ATOL multiplier, compared by exact
    # squared quantities; the Hermitian representative is exactly I/2.
    if valid:
        assert state_fidelity(rho, np.eye(2)/2) == 1
    else:
        with pytest.raises(ValueError, match="Hermitian"):
            state_fidelity(rho, np.eye(2)/2)


@pytest.mark.parametrize("bad", (np.zeros((0, 0)), np.ones((2, 3)), [1.],
                                  [[1., 0.], [0.]], [[True, 0.], [0., 0.]],
                                  np.ma.array(np.eye(2)/2, mask=True)))
def test_source_validation_still_rejects_malformed_or_missing_evidence(bad):
    with pytest.raises(ValueError):
        source_state(bad)


@pytest.mark.parametrize("multiple,rounded_units", ((1, 4), (3, 14)))
def test_rank_three_classical_fidelity_rounds_exact_midpoints(multiple, rounded_units):
    unit = np.nextafter(0., 1.)
    rho = np.diag([1., multiple*unit, 4*multiple*unit, 4*multiple*unit, 0.])
    sigma = np.diag([0., .5, .125, .125, .25])
    # All three common square-root products equal sqrt(m*u/2).
    # F=9m*u/2 exactly, with opposite ties-to-even directions for m=1,3.
    assert state_fidelity(rho, sigma) == rounded_units*unit
    assert state_fidelity(sigma, rho) == rounded_units*unit


@pytest.mark.parametrize("mass", (2.**-500, 2.**-1000, np.nextafter(0., 1.)))
def test_rank_three_nonsquare_classical_fidelity_keeps_its_output_scale(mass):
    rho = np.diag([.5, .25, .25, 0.])
    sigma = np.diag([mass, 4*mass, 9*mass, 1-14*mass])
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    expected = float(ctx.mpf(mass)*(ctx.sqrt(ctx.mpf('.5'))+ctx.mpf('2.5'))**2)
    assert state_fidelity(rho, sigma) == expected


@pytest.mark.parametrize("dense", (False, True))
def test_rank_eight_exact_cross_roots_and_complex_rotated_control(dense):
    rho = np.eye(8)/8
    if dense:
        v = np.array([1, 1j, -1, -1j, 1, 1j, -1, -1j])
        sigma = np.eye(8)/16 + np.outer(v, v.conj())/16
    else:
        sigma = np.diag([9/16]+[1/16]*7)
    # Spectrum sigma=(9,1,1,1,1,1,1,1)/16. The sqrt-products
    # (3,1,...,1)/sqrt(128) are irrational but their cross products
    # are rational: F=100/128. In the dense case, reading only the
    # diagonal would incorrectly return 1 instead of the spectral answer.
    assert state_fidelity(rho, sigma) == 100/128
    assert state_fidelity(sigma, rho) == 100/128
