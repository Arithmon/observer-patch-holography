"""Independent dense, quadrature and finite-difference checks of the pair map."""
from fractions import Fraction as F
import numpy as np
import pytest
from .fiber import PairFiber, clock_subspace
from .mass import curvature_ratio_interval,curvature_quadrature


@pytest.mark.parametrize('side',[4,6,8])
def test_same_number_same_momentum_reference_and_pair(side):
    s=clock_subspace(side)
    assert max(s['residuals'].values())<3e-14
    assert np.allclose(s['restriction'],np.diag([0,1]))
    f=PairFiber([np.pi]*3,side=side)
    p=f.pole()
    assert abs(p['energy']-1)<1e-14
    assert abs(p['contact_residue']-1)<1e-14


@pytest.mark.parametrize('momentum',[[0,0,0],[np.pi/2,np.pi,np.pi],[np.pi]*3])
def test_resolvent_pole_against_complete_dense_fiber(momentum):
    f=PairFiber(momentum,side=4)
    dense=np.zeros((f.dimension,f.dimension),complex)
    for i,b in enumerate(f.blocks):
        dense[4*i:4*i+4,4*i:4*i+4]=b
    dense+=np.outer(f.d,f.d.conj())
    vals,vecs=np.linalg.eigh(dense)
    p=f.pole()
    assert abs(vals[-1]-p['energy'])<3e-13
    assert abs(abs(np.vdot(f.d,vecs[:,-1]))**2-p['contact_residue'])<3e-13
    v=f.eigenstate()
    assert np.linalg.norm(f.apply(v)-p['energy']*v)<3e-13
    assert abs(np.vdot(v,v)-1)<3e-13
    assert f.exchange_residual(v)<3e-13
    assert vals[-2]<=2*np.sqrt(3)/8+1e-13


def test_rational_curvature_bounds_and_independent_quadrature():
    lo,hi=curvature_ratio_interval()
    assert F(17,16)<lo<hi<F(14,13)
    assert hi-lo<F(1,10**12)
    q=curvature_quadrature()*64
    assert float(lo)-2e-14<q<float(hi)+2e-14
    assert 59.86<64/float(hi)<64/float(lo)<59.87
    eps=.002
    e=PairFiber([np.pi+eps,np.pi,np.pi],side=32).pole()['energy']
    assert abs(2*(e-1)/eps**2-curvature_quadrature())<2e-8


def test_deleted_interaction_and_invalid_clock_maps_rejected():
    with pytest.raises(ValueError): PairFiber([0,0,0],interaction=0.)
    with pytest.raises(ValueError): clock_subspace(5)
    with pytest.raises(ValueError): curvature_ratio_interval(speed=F(1),interaction=F(1))
    f=PairFiber([.1,0,0],side=4)
    assert not f.is_finite_torus_momentum
    with pytest.raises(ValueError): f.exchange_residual(f.d)


@pytest.mark.parametrize('k',[[np.pi]*3,[0,np.pi,np.pi],[0,0,0]])
def test_same_pair_interaction_energy_and_hellmann_feynman(k):
    f=PairFiber(k,side=8)
    v=f.eigenstate();p=f.pole();eps=1e-5
    derivative=(PairFiber(k,side=8,interaction=1+eps).pole()['energy']-
                PairFiber(k,side=8,interaction=1-eps).pole()['energy'])/(2*eps)
    assert abs(derivative-p['contact_residue'])<2e-10
    interaction=abs(np.vdot(f.d,v))**2
    kinetic=np.vdot(v,f.free_apply(v)).real
    assert abs(interaction+kinetic-p['energy'])<3e-13
    assert abs(interaction-p['contact_residue'])<3e-13
