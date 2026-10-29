"""Independent physical-fiber, detector-completeness and native-clock joins."""
import math

import numpy as np
import pytest
from scipy.linalg import block_diag

from ..interacting_field.lattice import symbol
from .fiber import PairFiber, SWAP
from . import spectrum as m


@pytest.mark.parametrize('momentum',[[0,0,0],[0,np.pi,np.pi],[np.pi]*3,[np.pi/2,0,0]])
def test_physical_fiber_is_complete_and_invariant(momentum):
    f=PairFiber(momentum,side=4)
    q=m.physical_basis(f)
    assert np.linalg.norm(q.conj().T@q-np.eye(q.shape[1]))<3e-14
    partner_index=[]
    for k in f.k:
        partner=np.rint((np.array(momentum)-k)*4/(2*np.pi)).astype(int)%4
        partner_index.append(int((partner[0]*4+partner[1])*4+partner[2]))
    # The independently assembled exchange projector fixes exactly this basis.
    exchange=np.zeros((f.dimension,f.dimension),complex)
    for i,j in enumerate(partner_index):
        exchange[4*i:4*i+4,4*j:4*j+4]=-SWAP
    projector=(np.eye(f.dimension)+exchange)/2
    assert np.linalg.norm(q@q.conj().T-projector)<3e-14
    # Build the free symbol independently, not from the cached fiber blocks.
    free=block_diag(*[np.kron(symbol(list(k)),np.eye(2))+
                      np.kron(np.eye(2),symbol(list(np.array(momentum)-k))) for k in f.k])
    h=free+np.outer(f.d,f.d.conj())
    assert np.linalg.norm(h-m.full_matrix(f))<2e-14
    assert np.linalg.norm(h@q-q@(q.conj().T@h@q))<3e-14
    assert np.linalg.norm(q@(q.conj().T@f.d)-f.d)<1e-14


@pytest.mark.parametrize('momentum',[[0,0,0],[0,np.pi,np.pi],[np.pi]*3])
def test_every_spectral_weight_and_dark_multiplicity_retained(momentum):
    f=PairFiber(momentum,side=4)
    s=m.complete_spectrum(f)
    atoms=m.spectral_atoms(s)
    assert sum(row['multiplicity'] for row in atoms)==len(s['energies'])==120
    assert sum(row['contact_weight'] for row in atoms)==pytest.approx(1,abs=2e-14)
    assert np.min(s['weights'])>=0
    # First and second moments independently probe the complete measure.
    hd=f.apply(f.d)
    assert np.dot(s['weights'],s['energies'])==pytest.approx(1,abs=2e-14)
    assert np.dot(s['weights'],s['energies']**2)==pytest.approx(np.vdot(hd,hd).real,abs=2e-14)
    p=f.pole()
    assert s['energies'][-1]==pytest.approx(p['energy'],abs=2e-14)
    assert s['weights'][-1]==pytest.approx(p['contact_residue'],abs=2e-14)
    assert s['energies'][-2]<=2*math.sqrt(3)/8+2e-14
    # At least one exact dark eigenspace survives in each retained fiber.
    assert any(row['contact_weight']<1e-25 for row in atoms)


def test_complete_contact_detector_matches_independent_evolution_and_rejects_pole_only():
    f=PairFiber([0,0,0],side=4)
    s=m.complete_spectrum(f)
    laws=[m.detector_law(f,t,s) for t in (0.,2.,7.,20.)]
    for law in laws:
        assert law['independent_state_residual']<2e-12
        assert law['complete_mass_residual']<5e-14
        assert law['orthogonal_outputs_residual']<2e-14
        assert sum(law['bound_filter_complete_probabilities'])==1
    assert laws[0]['probabilities'][0]==pytest.approx(1,abs=1e-14)
    assert laws[0]['single_pole_contact_probability_error']>.1
    assert max(x['probabilities'][1] for x in laws[1:])>.1
    assert max(x['probabilities'][0] for x in laws)-min(x['probabilities'][0] for x in laws)>.1
    # A normalized bound eigenstate is a distinct physical preparation.
    z=s['weights'][-1]
    assert abs(np.vdot(f.d,f.eigenstate()))**2==pytest.approx(z,abs=2e-14)
    assert abs(z-1)>.05


def test_eta_has_exact_contact_pole_but_energy_needs_relative_reference():
    f=PairFiber([np.pi]*3,side=4)
    s=m.complete_spectrum(f)
    assert np.sum(s['weights'][:-1])<1e-28
    for t in (0.,2.,7.):
        law=m.detector_law(f,t,s)
        assert law['probabilities'][0]==pytest.approx(1,abs=1e-14)
    # Contact return alone cannot read this global phase; use the N=2 r/d pair.
    a=m.clock_subspace()
    psi=(a['reference']+a['pair'])/math.sqrt(2)
    t=.7
    from scipy.sparse.linalg import expm_multiply
    evolved=expm_multiply(-1j*t*m.full_matrix(f),psi)
    assert abs(np.vdot(psi,evolved))**2==pytest.approx((1+math.cos(t))/2,abs=2e-14)


def test_common_eta_clock_uses_complete_unconditioned_native_records():
    result=m.native_clock_attachment()
    assert result['field_energy_gap']==1
    assert result['two_levels_have_the_same_particle_number']==2
    assert result['faithful_M2_Gibbs_residual']<1e-15
    for row in result['ideal_native_complete_readouts'].values():
        assert row['full_mass']==pytest.approx(1,abs=5e-14)
        assert row['energy_readout_formula_residual']<2e-14
        assert row['initial_outside']>.96
        assert row['continued_outside']>.004
        assert 0<row['energy_r']<.02
        assert 0<row['energy_d']<.02


def test_volume_poles_keep_actual_torus_momenta_and_open_continuum_scope():
    energies=[]
    for side in (4,6,8,12,16):
        f=PairFiber([0,0,0],side=side)
        p=f.pole()
        assert p['finite_torus_momentum']
        assert p['secular_residual']<1e-13
        energies.append(p['energy'])
    assert abs(energies[-1]-energies[-2])<1e-7
    assert abs(energies[0]-energies[-1])>1e-4
    with pytest.raises(ValueError,match='finite-torus'):
        m.complete_spectrum(PairFiber([.1,0,0],side=4))
    with pytest.raises(ValueError,match='L<=6'):
        m.complete_spectrum(PairFiber([0,0,0],side=8))
    with pytest.raises(ValueError):
        m.detector_law(PairFiber([0,0,0],side=4),float('nan'))


def test_interaction_deletion_removes_the_isolated_pole_on_same_physical_space():
    for momentum in ([0,0,0],[np.pi]*3):
        f=PairFiber(momentum,side=4)
        result=m.deleted_contact_control(f)
        assert result['eigenvalues_above_free_norm_bound']==0
        assert result['maximum_energy']<=result['free_norm_bound']+1e-12
        assert result['original_interacting_pole_energy']>=1-1e-14
        assert not result['original_bound_pole_survives']


def test_complete_contact_instrument_preserves_an_entangled_reference():
    f=PairFiber([0,np.pi,np.pi],side=4)
    s=m.complete_spectrum(f)
    v=s['basis']@s['vectors'][:,-1]
    orthogonal=f.d-v*np.vdot(v,f.d)
    orthogonal/=np.linalg.norm(orthogonal)
    # Complex reference correlations and unequal Schmidt coefficients.
    joint=np.column_stack((math.sqrt(.3)*v,1j*math.sqrt(.7)*orthogonal))
    t=3.
    from scipy.linalg import expm
    evolved=expm(-1j*t*m.full_matrix(f))@joint
    yes=np.outer(f.d,f.d.conj())@evolved
    no=evolved-yes
    recovered=yes.T@yes.conj()+no.T@no.conj()
    assert np.linalg.norm(recovered-joint.T@joint.conj())<3e-14
    assert abs(np.linalg.norm(yes)**2+np.linalg.norm(no)**2-1)<3e-14
