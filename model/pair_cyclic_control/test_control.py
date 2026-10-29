"""Exact secular signs and an independent full physical-fiber intertwiner."""
from fractions import Fraction as F
import numpy as np
import pytest
from scipy.linalg import expm
from .control import contact_resolvent, pole_interval, matrices
from ..bound_pair.fiber import PairFiber
from ..bound_pair.spectrum import physical_basis, full_matrix


@pytest.mark.parametrize('speed,interaction',[(F(1,8),F(1)),(F(1,16),F(3,4)),(F(1,4),F(9,10))])
def test_exact_rational_pole_interval(speed,interaction):
    lo, hi = pole_interval(speed,interaction)
    assert hi-lo == 4*speed/2**64
    assert interaction*contact_resolvent(lo,speed) >= 1
    assert interaction*contact_resolvent(hi,speed) < 1
    _,_,h = matrices(speed,interaction)
    root = np.linalg.eigvalsh(h)[-1]
    assert abs(root-float((lo+hi)/2)) < 2e-14


def test_complete_cyclic_subspace_intertwines_full_physical_action_and_contact():
    f = PairFiber([0,0,0],side=4)
    basis = physical_basis(f)
    h = basis.conj().T@full_matrix(f)@basis
    d = basis.conj().T@f.d
    free, u, reduced = matrices()
    t = h-np.outer(d,d.conj())
    columns = []
    for eigenvalue, amplitude in zip(free,u,strict=True):
        # Independent spectral projector polynomial in the full free operator.
        vector = d.copy()
        for other in free:
            if other != eigenvalue:
                vector = (t@vector-other*vector)/(eigenvalue-other)
        columns.append(vector/amplitude)
    inclusion = np.column_stack(columns)
    assert h.shape == (120,120)
    np.testing.assert_allclose(inclusion.conj().T@inclusion,np.eye(7),atol=2e-13)
    np.testing.assert_allclose(h@inclusion,inclusion@reduced,atol=2e-13)
    np.testing.assert_allclose(inclusion@u,d,atol=2e-13)
    contact = np.outer(d,d.conj()); small_contact = np.outer(u,u)
    np.testing.assert_allclose(contact@inclusion,inclusion@small_contact,atol=2e-13)
    # Complete two-step contact/complement history, with a complex untouched
    # reference; this does not normalize or select surviving contact outcomes.
    rng = np.random.default_rng(1152)
    state = rng.normal(size=(7,2))+1j*rng.normal(size=(7,2))
    state /= np.linalg.norm(state)
    large = [inclusion@state]; small = [state]
    for time in (.37,-1.3):
        large = [p@expm(-1j*time*h)@v for v in large
                 for p in (contact,np.eye(120)-contact)]
        small = [p@expm(-1j*time*reduced)@v for v in small
                 for p in (small_contact,np.eye(7)-small_contact)]
    for actual, reduced_output in zip(large,small,strict=True):
        np.testing.assert_allclose(actual,inclusion@reduced_output,atol=3e-13)
    assert abs(sum(np.linalg.norm(v)**2 for v in small)-1) < 1e-13


def test_invalid_domains_do_not_inherit_the_exact_certificate():
    for value in (True,0.125,F(0),F(-1,8)):
        with pytest.raises(ValueError): pole_interval(value)
    with pytest.raises(ValueError): pole_interval(F(1),F(1))
    with pytest.raises(ValueError): pole_interval(iterations=True)
    with pytest.raises(ValueError): contact_resolvent(F(1,4))
