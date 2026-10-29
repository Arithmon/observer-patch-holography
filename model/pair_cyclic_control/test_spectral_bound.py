"""Independent exact edge checks and complete-field spectral controls."""
from fractions import Fraction as F
from itertools import product
import numpy as np
import pytest
from .spectral_bound import enclosure, canonical_clock_bound
from ..bound_pair.fiber import PairFiber
from ..bound_pair.spectrum import complete_spectrum


def test_canonical_exact_inequalities_and_stronger_clock_coefficient():
    result=canonical_clock_bound()
    assert F(7,16)**2 > F(3,16)
    assert F(7,6)*(F(7,6)-1) > F(3,16)
    assert result == {'lower':-F(7,16),'upper':F(7,6),'center':F(35,96),
                      'centered_norm_upper':F(77,96),
                      'rate_error_halfdiamond_upper':F(77,2560)}
    assert result['rate_error_halfdiamond_upper'] < F(9,256)
    assert canonical_clock_bound(F(0))['rate_error_halfdiamond_upper']==0
    assert canonical_clock_bound(F(-1,2))['rate_error_halfdiamond_upper']==F(77,5120)


def test_full_physical_spectrum_keeps_all_momenta_and_scattering_states():
    result=canonical_clock_bound();count=0
    for momentum in product(range(4),repeat=3):
        fiber=PairFiber(list(np.pi*np.array(momentum)/2),side=4)
        energies=complete_spectrum(fiber)['energies']
        count+=len(energies)
        assert min(energies)>float(result['lower'])
        assert max(energies)<float(result['upper'])
        assert max(abs(energies-float(result['center'])))<float(result['centered_norm_upper'])
    assert count==8128


def test_generic_noncommuting_contact_projection_and_failed_stronger_edges():
    # A two-state example saturates the analytic upper-root formula; P and T
    # do not commute, so the proof must use compression rather than commute them.
    t=np.array([[0.,.25],[.25,0.]])
    p=np.diag([1.,0.]);h=t+p
    assert np.array_equal(p@t@p,np.zeros((2,2)))
    assert np.linalg.norm(p@t-t@p)>0
    assert np.linalg.eigvalsh(h)[-1]==pytest.approx((1+np.sqrt(1.25))/2)
    with pytest.raises(ValueError): enclosure(F(3,16),F(1),F(7,16),F(9,8))
    with pytest.raises(ValueError): enclosure(F(3,16),F(1),F(2,5),F(7,6))
    with pytest.raises(ValueError): canonical_clock_bound(True)
    with pytest.raises(ValueError): canonical_clock_bound(F(2))
    with pytest.raises(ValueError): canonical_clock_bound(beta_error=F(-1))
