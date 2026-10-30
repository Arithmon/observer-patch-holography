"""Independent complete-sector, spectral, particle-assignment and source controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from .spectrum import FieldSpectrum
from .dimer import certificate, setup, candidate_projectors, spin_generators, charge_conjugation
from . import source as readout
from ..full_fock_field.field import FieldGeometry
from ..observer.source_scalar import compiler


def jw(modes):
    minus = np.array([[0,1],[0,0]], complex)
    out = []
    for mode in range(modes):
        value = np.ones((1,1),complex)
        for j in reversed(range(modes)):
            value = np.kron(value, minus if j == mode else np.diag([1,-1]) if j < mode else np.eye(2))
        out.append(value)
    return out


def independent_h(geometry):
    a = jw(geometry.modes)
    single = geometry.single_h().toarray()
    h = sum(single[i,j]*(a[i].conj().T@a[j])
            for i in range(geometry.modes) for j in range(geometry.modes))
    h += sum(geometry.interaction*(a[2*x].conj().T@a[2*x])@(a[2*x+1].conj().T@a[2*x+1])
             for x in range(geometry.volume))
    return h, a


@pytest.mark.parametrize('sides',[(2,1,1),(1,2,1),(1,1,2),(2,2,1)])
def test_all_number_lehmann_response_agrees_with_independent_time_domain_and_sum_rules(sides):
    geometry = FieldGeometry(sides)
    s = FieldSpectrum(geometry, beta=1.3, chemical_potential=.2)
    h, a = independent_h(geometry)
    n = sum(x.conj().T@x for x in a)
    k = h-.2*n
    assert np.array_equal(k,s.generator)
    rho = expm(-1.3*k);rho /= np.trace(rho)
    assert np.linalg.norm(rho-s.reference,2)<2e-14
    for kind,index in [('single',0),('single',geometry.modes-1),('pair',0)]:
        op = a[index].conj().T if kind=='single' else a[0].conj().T@a[1].conj().T
        m = s.measure(op,1 if kind=='single' else 2)
        assert np.all(m.forward>=0) and np.all(m.reverse>=0)
        assert np.count_nonzero(m.forward[~m.allowed])==0
        for time in (-.71,0.,.13):
            u=expm(-1j*time*k);lower=op.conj().T
            direct_plus=np.trace(rho@(u.conj().T@lower@u)@op)
            direct_minus=np.trace(rho@op@(u.conj().T@lower@u))
            assert abs(m.correlation(time)-direct_plus)<3e-13
            assert abs(m.correlation(time,'reverse')-direct_minus)<3e-13
        r=s.certify(kind,index)
        assert r['maximum_moment_residual']<3e-13
        assert r['detailed_balance_residual']<3e-13
        assert abs(sum(row['forward_weight'] for row in m.lines())-m.moment(0))<3e-13
        if kind=='single':assert abs(m.moment(0)+m.moment(0,'reverse')-1)<3e-13


@pytest.mark.parametrize('beta,u',[(.25,0.),(1.,1.),(2.,.75)])
def test_atomic_hubbard_lines_and_pair_commutator_have_closed_thermal_weights(beta,u):
    s=FieldSpectrum(FieldGeometry((1,1,1),interaction=u),beta=beta)
    ratio=math.exp(-beta*u);z=3+ratio
    op,q=s.probe('single',0);m=s.measure(op,q)
    for t in (.2,.9):
        assert abs(m.correlation(t,'anticommutator')-(2+(1+ratio)*np.exp(-1j*u*t))/z)<1e-14
    op,q=s.probe('pair',0);m=s.measure(op,q)
    assert abs(m.moment(0)-1/z)<1e-14
    assert abs(m.moment(0,'reverse')-ratio/z)<1e-14
    assert abs(m.correlation(.3,'commutator')-(1-ratio)*np.exp(-.3j*u)/z)<1e-14


@pytest.mark.parametrize('axis',[0,1,2])
@pytest.mark.parametrize('speed,u',[(.125,1.),(.25,0.),(.4,.75)])
def test_two_particle_assignments_have_same_ground_charge_spin_parity_and_gap(axis,speed,u):
    sides=[1,1,1];sides[axis]=2;g=FieldGeometry(tuple(sides),speed,u)
    r=certificate(g);s=setup(g);p=candidate_projectors(s)
    assert r['maximum_matrix_and_residue_residual']<2e-12
    assert r['physical_identification'] is False
    assert r['absolute_physical_mass_predicted'] is False
    t=speed/2;delta=math.sqrt(u*u+16*t*t)
    assert abs(s.energies[0]+(u+delta)/2)<1e-13
    assert np.linalg.norm(s.number@s.reference-2*s.reference)<1e-13
    a=jw(4);n=sum(x.conj().T@x for x in a)
    for name,charge in [('charge_minus_one_hole',-1),('charge_plus_one_particle',1)]:
        assert np.linalg.norm((n-2*np.eye(16))@p[name]-charge*p[name])<1e-13
        assert abs(np.trace(p[name])-2)<1e-13
    # Removing the required half-filling subtraction breaks the particle/hole
    # gap match. It cannot silently count as the same vacuum/generator.
    if u:
        with pytest.raises(ValueError,match='half-filled'):
            candidate_projectors(FieldSpectrum(g,chemical_potential=.49*u,reference='ground'))


def test_particle_hole_unitary_and_internal_spin_are_actual_complete_operator_maps():
    s=setup();c=charge_conjugation(s);spin=spin_generators(s)
    assert np.array_equal(c.conj().T@c,np.eye(16))
    assert np.array_equal(c@s.generator@c.conj().T,s.generator)
    assert np.array_equal(c@s.number@c.conj().T,4*np.eye(16)-s.number)
    assert all(np.linalg.norm(x@s.generator-s.generator@x)==0 for x in spin)
    assert np.array_equal(spin[0]@spin[1]-spin[1]@spin[0],1j*spin[2])
    r=certificate();low=r['candidates'][:2]
    assert abs(low[0]['single_local_probe_residue']-(1+1/math.sqrt(17))/4)<1e-13
    assert abs(low[1]['single_local_probe_residue']-low[0]['single_local_probe_residue'])<1e-13
    assert r['neutral_triplet_gap']<low[0]['model_frequency_gap']
    assert abs(r['resolved_charge_one_frequency_splitting']-.125)<1e-14


def test_interaction_creates_a_resolved_second_charge_one_pole_on_the_same_ground_reference():
    free,interacting=setup(FieldGeometry((2,1,1),interaction=0.)),setup()
    for s,count in [(free,2),(interacting,4)]:
        op,q=s.probe('single',0);measure=s.measure(op,q)
        assert len([line for line in measure.lines() if line['visible']])==count
        assert abs(measure.correlation(0,'anticommutator')-1)<1e-14
    assert certificate(FieldGeometry((2,1,1),interaction=0.))['higher_single_local_probe_residue']==0
    assert certificate()['higher_single_local_probe_residue']>.18


@pytest.mark.parametrize('label',['XI','YX','ZZ','II'])
def test_controlled_pauli_source_word_keeps_scalar_phase_and_all_entangled_inputs(label):
    word=readout.controlled_pauli(label,2)
    actual=compiler.execute(word,3)
    expected=np.zeros((8,8),complex);expected[:4,:4]=np.eye(4);expected[4:,4:]=readout.pauli_matrix(label)
    assert np.linalg.norm(actual-expected,2)<3e-12
    assert np.linalg.norm(actual.conj().T@actual-np.eye(8),2)<3e-12


def test_ground_particle_and_hole_correlations_come_from_complete_source_records():
    report=readout.probe_readout(setup(),duration=F(1,4),steps=4)
    assert report['independent_preparations_required']==4
    assert report['preparation_source_pulses_per_run']==10249
    assert report['total_source_pulses_with_preparation']==52484
    assert report['maximum_correlation_residual']<7e-7
    assert len(report['records'])==4
    for row in report['records']:
        assert len(row['pointer_probabilities'])==2
        assert min(row['pointer_probabilities'])>=0
        assert abs(row['complete_probability']-1)<2e-9
        assert abs(sum(row['pointer_probabilities'])+row['preparation_abort_probability']
                   +row['continuation_abort_probability']-1)<2e-9
    plus=complex(*report['observed_forward']);minus=complex(*report['observed_reverse'])
    assert abs(plus-minus.conjugate())<2e-12


def test_pair_source_readout_reconstructs_hard_core_pair_correlator_without_conditioning():
    s=FieldSpectrum(FieldGeometry((1,1,1)),beta=1.)
    report=readout.probe_readout(s,'pair',0,F(1,4),1)
    assert report['independent_preparations_required']==16
    assert report['maximum_correlation_residual']<3e-11
    plus=complex(*report['observed_forward']);minus=complex(*report['observed_reverse'])
    assert abs(minus-math.exp(-1)*plus)<3e-11
    assert all(abs(row['complete_probability']-1)<2e-9 for row in report['records'])


def test_complete_source_readout_rejects_external_preparation_and_new_abort_amplitudes(monkeypatch):
    s=FieldSpectrum(FieldGeometry((1,1,1)))
    prepared=readout.prepare(s)
    with pytest.raises(TypeError,match='prepared'):
        readout.correlator(s,'XI','XI',prepared=prepared)
    with pytest.raises(ValueError):prepared.state[0]=0
    with pytest.raises(TypeError):prepared.resources['native_pulses']=0
    native=compiler.native_pulse
    def leaking(kind,angle):
        out=native(kind,angle).copy();mix=np.eye(6,dtype=complex)
        theta=.001;mix[0,0]=mix[5,5]=math.cos(theta);mix[0,5]=mix[5,0]=1j*math.sin(theta)
        return mix@out
    monkeypatch.setattr(compiler,'native_pulse',leaking)
    with pytest.raises(ArithmeticError,match='spectroscopy'):
        readout._correlator(s,'XI','XI',F(0),1,0.,prepared)


@pytest.mark.parametrize('bad',[True,float('nan'),float('inf'),0.,-1.])
def test_malformed_reference_parameters_fail_closed(bad):
    with pytest.raises(ValueError):FieldSpectrum(beta=bad)


def test_invalid_charge_reference_source_size_and_time_are_not_silently_changed():
    s=FieldSpectrum()
    op,q=s.probe('single',0)
    with pytest.raises(ValueError,match='charge'):s.measure(op,2)
    with pytest.raises(ValueError):s.measure(np.full((16,16),np.nan),1)
    with pytest.raises(ValueError):FieldSpectrum(FieldGeometry((2,1,1),speed=0),reference='ground')
    with pytest.raises(ValueError):readout.prepare(FieldSpectrum(FieldGeometry((2,2,1))))
    with pytest.raises(ValueError):readout.forward_word(s,.25,4)
    with pytest.raises(ValueError):readout.components('pair',True,4)
    with pytest.raises(ValueError):s.measure(op,True)
