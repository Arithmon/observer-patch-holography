"""Independent complete-matrix, exchange and continuum-response controls."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy.linalg import block_diag

from ..bound_pair.fiber import PairFiber, SINGLET
from ..bound_pair.spectrum import physical_basis
from ..interacting_field.lattice import symbol
from .scattering import ContactFiber, free_green, response, volume_enclosure


def independent_matrices(fiber):
    # Explicit Kronecker sum, independently assembled from the one-body law.
    blocks = [np.kron(symbol(list(k), fiber.speed), np.eye(2))+
              np.kron(np.eye(2), symbol(list(np.array(fiber.momentum)-k), fiber.speed))
              for k in fiber.k]
    h0 = block_diag(*blocks)
    d = np.concatenate([SINGLET/np.sqrt(fiber.volume) for _ in blocks])
    v = fiber.interaction*np.outer(d, d.conj())
    return h0, v, d


@pytest.mark.parametrize('momentum', [(0.,0.,0.), (np.pi,0.,np.pi), (np.pi/2,)*3])
@pytest.mark.parametrize('interaction', [0., .05, 1.])
def test_complete_resolvent_and_optical_operator_identity(momentum, interaction):
    f = ContactFiber(momentum, interaction=interaction)
    h0,v,d = independent_matrices(f)
    z = .17+.23j
    r0 = np.linalg.inv(z*np.eye(f.dimension)-h0)
    target = np.linalg.inv(z*np.eye(f.dimension)-h0-v)
    expected_t = v+v@target@v
    actual = f.resolvent(np.eye(f.dimension), z)
    actual_t = f.t_apply(np.eye(f.dimension), z)
    assert np.linalg.norm(actual-target,2) < 5e-13
    assert np.linalg.norm(actual_t-expected_t,2) < 5e-13
    assert np.linalg.norm(actual_t-actual_t.conj().T-
                          actual_t.conj().T@(r0-r0.conj().T)@actual_t,2) < 2e-12
    report = response(momentum, z, side=4, interaction=interaction)
    assert abs(report['free_green']-np.vdot(d,r0@d)) < 2e-14
    assert abs(report['contact_green']-np.vdot(d,target@d)) < 2e-14
    assert report['optical_residual'] < 2e-14
    assert report['finite_torus_momentum']
    if interaction == 0:
        assert np.count_nonzero(actual_t) == 0
        assert report['contact_green'] == report['free_green']
    else:
        assert np.linalg.norm(actual_t,2) > 0


def test_exchange_projector_and_every_physical_spin_channel_are_retained():
    f = ContactFiber((0.,np.pi,np.pi))
    old = PairFiber(list(f.momentum),side=4)
    basis = physical_basis(old)
    expected_projector = basis@basis.conj().T
    actual_projector = f.physical_project(np.eye(f.dimension))
    assert np.linalg.norm(actual_projector-expected_projector,2) < 1e-14
    assert abs(np.trace(actual_projector)-120) < 1e-12
    rng = np.random.default_rng(714)
    x = rng.normal(size=(f.dimension,3))+1j*rng.normal(size=(f.dimension,3))
    x = f.physical_project(x)
    out = f.resolvent(x,.2+.13j)
    assert np.linalg.norm(f.physical_project(out)-out) < 5e-13
    assert np.linalg.norm(f.physical_project(f.d)-f.d) < 1e-14
    # Channels orthogonal to the contact have zero T action; they are not deleted.
    dark = x-f.d[:,None]*(f.d.conj()@x)
    assert np.linalg.norm(dark)>1
    assert np.linalg.norm(f.t_apply(dark,.2+.13j)) < 1e-13
    # They can have nonzero free propagation, so T=0 does not mean zero state.
    assert np.linalg.norm(f.free_resolvent(dark,.2+.13j))>1


def test_same_bound_particle_is_a_pole_of_the_scattering_denominator():
    f = PairFiber([0,0,0],side=4)
    pole = f.pole(); e,z = pole['energy'],pole['contact_residue']
    eta = 1e-6
    r = response([0,0,0],e+1j*eta,side=4)
    assert abs((1j*eta)*r['tau']-f.interaction**2*z) < 2e-6
    assert abs((1j*eta)*r['contact_green']-z) < 2e-6
    # Direct full physical eigensystem gives the same singular residue.
    h0,v,d=independent_matrices(ContactFiber())
    b=physical_basis(f); values,vectors=np.linalg.eigh(b.conj().T@(h0+v)@b)
    assert abs(values[-1]-e)<1e-13
    assert abs(abs(np.vdot(b@vectors[:,-1],d))**2-z)<1e-13


def test_exact_rational_alias_bounds_and_volume_convergence():
    bound = volume_enclosure(128,F(1,4))
    assert type(bound['free_green_error']) is F
    assert bound['free_green_error'] < F(1295,10**10)
    assert bound['tau_error'] < F(3236,10**9)
    assert bound['contact_green_error'] == bound['tau_error']
    assert bound['on_shell_limit_claimed'] is False
    # Three independent Fourier quadratures: their difference must be below
    # the sum of the two certified errors to the common fixed-eta integral.
    total = (.31,.77,1.12); frequency = .17+1j
    rows = {n:response(total,frequency,side=n) for n in (16,32,64)}
    for a,b in ((16,32),(32,64),(16,64)):
        ea,eb=volume_enclosure(a,F(1)),volume_enclosure(b,F(1))
        for key,error in [('free_green','free_green_error'),('tau','tau_error'),
                          ('contact_green','contact_green_error')]:
            assert abs(rows[a][key]-rows[b][key]) < float(ea[error]+eb[error])+5e-14
    assert not rows[16]['finite_torus_momentum']
    # Atomic c=0 has no quadrature error, including a nonzero interaction.
    atom = volume_enclosure(4,F(1,8),speed=F(0),interaction=F(3))
    assert atom['free_green_error']==atom['tau_error']==atom['contact_green_error']==0
    r=response(total,frequency,side=4,speed=0,interaction=3)
    assert abs(r['contact_green']-1/(frequency-3))<1e-14


@pytest.mark.parametrize('z', [0,1j*0,-1j,complex(float('nan'),1),complex(0,float('inf')),True])
def test_nonretarded_or_nonfinite_frequencies_rejected(z):
    with pytest.raises(ValueError): response([0,0,0],z)


@pytest.mark.parametrize('bad', [0,-1,F(-1,3),.1,True])
def test_broadening_bound_requires_positive_exact_rational(bad):
    with pytest.raises(ValueError): volume_enclosure(16,bad)


def test_explicit_host_and_physical_domains():
    with pytest.raises(ValueError): free_green([0,0,0],1j,side=257)
    with pytest.raises(ValueError): ContactFiber(side=25)
    with pytest.raises(ValueError): ContactFiber((.1,0,0)).physical_project(np.eye(256))
    with pytest.raises(ValueError): ContactFiber().resolvent(np.zeros(3),1j)
    with pytest.raises(ValueError): response([0,0,0],1j,interaction=-1)
    with pytest.raises(ValueError): response(np.ma.array([0.,0.,0.],mask=[0,1,0]),1j)
    with pytest.raises(ValueError): ContactFiber().resolvent(np.full(256,np.nan),1j)
