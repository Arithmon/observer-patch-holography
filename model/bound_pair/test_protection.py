"""Independent position-space and complete-outcome checks for one pair band."""
from fractions import Fraction as F
import math
import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply
from .fiber import PairFiber
from .protection import PeriodicPair,canonical_band,packet_coefficients,report,tail_bound


def test_all_fibers_are_eigenstates_of_independent_periodic_CAR_action():
    space,band,energy,residue,exchange=canonical_band()
    assert band.shape==(8128,64)
    assert exchange<3e-13
    assert np.linalg.norm(band.conj().T@band-np.eye(64))<3e-13
    assert np.linalg.norm(space.hamiltonian@band-band*energy)<3e-13
    assert np.all(energy>=1-2e-14)
    assert energy.max()-energy.min()>.08
    assert np.all(residue>0)
    assert np.all(residue<=1+2e-14)
    # A single hop cannot remain an onsite pair. This is the premise that
    # strengthens the elementary Weyl estimate to the gap U-b.
    contact=space.free_h[space.contact_indices][:,space.contact_indices]
    assert contact.nnz==0
    assert np.max(abs((space.single_h-space.single_h.conj().T).data),initial=0)==0


def test_complete_wavepacket_motion_and_unchanged_protected_band():
    space,band,energy,_,_=canonical_band()
    a=packet_coefficients()
    initial=band@a
    final=band@(a*np.exp(-16j*energy))
    independent=expm_multiply(-16j*space.hamiltonian,initial)
    assert np.linalg.norm(independent-final)<2e-12
    assert np.linalg.norm(final-band@(band.conj().T@final))<2e-13
    out=report()
    assert out['complete_contact_record_TV']>.01
    assert abs(out['initial_x_neighbor_asymmetry'])<2e-14
    assert abs(out['final_x_neighbor_asymmetry'])>.005
    assert out['flattened_band_density_change']<2e-14
    for record in out['records']:
        assert abs(record['trace']-1)<2e-13
        assert abs(sum(record['site_particle_number'])-2)<2e-13
        assert abs(sum(record['site_contact_probability'])+record['separated_probability']-1)<2e-13
        for row in record['relative_tails']:
            assert row['probability']<=float(F(row['rational_upper']))+2e-13


def test_band_gauge_changes_leave_the_same_physical_packet():
    _,band,energy,_,_=canonical_band()
    a=packet_coefficients()
    phase=np.exp(1j*np.arange(64)**2/19)
    state=band@(a*np.exp(-8j*energy))
    recharts=(band*phase)@(a/phase*np.exp(-8j*energy))
    assert np.linalg.norm(state-recharts)<2e-13


def test_tail_index_and_rational_outward_relaxation():
    assert tail_bound(0)==1
    assert tail_bound(1)==F(49,81)
    assert all(tail_bound(r+1)<tail_bound(r) for r in range(6))
    assert F(7,16)**2>F(3,16)
    with pytest.raises(ValueError): tail_bound(-1)
    with pytest.raises(ValueError): tail_bound(True)


@pytest.mark.parametrize('scale',[0.,.25,.5,1.])
def test_contact_section_and_gap_survive_contraction_to_onsite_pair(scale):
    f=PairFiber([np.pi/2,np.pi,np.pi],side=4,speed=scale/8)
    v=f.eigenstate(); p=f.pole()
    overlap=np.vdot(f.d,v)
    assert overlap.real>0 and abs(overlap.imag)<2e-14
    assert abs(overlap.real**2-p['contact_residue'])<3e-13
    assert p['energy']>=1-2e-14
    assert 1-2*math.sqrt(3)*scale/8>=1-2*math.sqrt(3)/8
    if scale==0:
        assert np.linalg.norm(v-f.d)<2e-14


def test_deleted_interaction_is_not_silently_kept_in_the_protection_class():
    with pytest.raises(ValueError): PairFiber([0,0,0],side=4,interaction=0.)
    free=PeriodicPair(interaction=0.)
    space,band,energy,_,_=canonical_band()
    assert np.linalg.norm(free.hamiltonian@band-band*energy)>.5
    # The periodic operator itself remains a valid free model; only the
    # particular isolated repulsive-band premise is lost.
    assert np.max(abs((free.hamiltonian-free.free_h).data),initial=0)==0


def test_unphysical_first_quantized_vector_cannot_be_compressed_as_a_pair():
    space=PeriodicPair()
    bad=np.zeros(4*space.volume,complex);bad[0]=1
    with pytest.raises(ValueError,match='exchange'):
        space.fiber_state([0,0,0],bad)
