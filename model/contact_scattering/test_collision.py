"""Full-source collision controls, with all exports and reference columns."""
from fractions import Fraction as F
import numpy as np
import pytest
from scipy.linalg import expm

from ..field_thermal.thermal import failure_probability
from ..full_fock_field.field import FullFock
from .collision import execute, exact_echo, geometry, preparation


def jw(mode, modes=4):
    # Independent entrywise CAR sign construction.
    out=np.zeros((1<<modes,1<<modes),complex)
    for column in range(1<<modes):
        if column&(1<<mode):
            parity=sum(bool(column&(1<<j)) for j in range(mode))
            out[column^(1<<mode),column]=(-1)**parity
    return out


def independent_h(axis, interaction):
    sigma=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))[axis]
    a=[jw(j) for j in range(4)];h=np.zeros((16,16),complex)
    for s in range(2):
        for t in range(2):
            term=1j/16*sigma[s,t]*a[s].conj().T@a[2+t]
            h+=term+term.conj().T
    for site in (0,1):
        p,q=a[2*site],a[2*site+1]
        h+=interaction*p.conj().T@p@q.conj().T@q
    return h


@pytest.mark.parametrize('axis', [0,1,2])
def test_interaction_echo_complete_source_reads_and_independent_hamiltonian(axis):
    value=execute(axis)
    report=value['report'];initial=preparation()
    target=expm(4j*independent_h(axis,0))@expm(-4j*independent_h(axis,1))@initial
    expected_empty=target.copy();expected_empty[1::2]=0
    expected_full=target.copy();expected_full[::2]=0
    assert np.linalg.norm(exact_echo(axis)-target)<1e-13
    error=np.sqrt(sum(np.linalg.norm(a-b)**2 for a,b in zip(
        value['outputs'],(expected_empty,expected_full),strict=True))+
        failure_probability(value['first_failure_events']))
    assert abs(report['complete_vector_error']-error)<1e-13
    assert error < float(F(report['ideal_word_half_diamond_bound']))
    assert report['pointer_probabilities'][0] > .046
    assert report['exact_empty_probability']-float(F(report['ideal_word_half_diamond_bound'])) > .036
    assert abs(sum(report['pointer_probabilities'])+report['failure_probability']-1)<2e-10
    assert report['reference_gram_residual']<2e-10
    assert report['source_resources']['native_pulses']==22531
    assert report['source_resources']['proper_code_service_events']==45062
    assert not report['asymptotic_scattering_or_physical_cross_section_claimed']


def test_free_control_keeps_complete_zero_collision_signal():
    value=execute(interaction=0)
    assert value['report']['pointer_probabilities'][0]<1e-20
    assert value['report']['exact_empty_probability']<1e-28
    assert np.linalg.norm(value['outputs'][1]-preparation())<1e-9
    assert value['report']['source_resources']['native_pulses']==22531


def test_complete_source_echo_accepts_all_number_coherences():
    rng=np.random.default_rng(919)
    initial=rng.normal(size=(16,2))+1j*rng.normal(size=(16,2))
    initial/=np.linalg.norm(initial)
    value=execute(duration=F(1,8),steps=2,initial=initial)
    expected=exact_echo(duration=.125,initial=initial)
    for pointer,out in enumerate(value['outputs']):
        target=expected.copy();target[(1-pointer)::2]=0
        assert np.linalg.norm(out-target)<float(F(value['report']['ideal_word_half_diamond_bound']))
    # Every particle-number sector is actually populated and retained.
    for n in range(5):
        indices=[i for i in range(16) if i.bit_count()==n]
        assert sum(np.linalg.norm(out[indices])**2 for out in value['outputs'])>0


@pytest.mark.parametrize('kwargs', [{'axis':3},{'steps':0},{'steps':65},{'duration':.5},
                                   {'duration':F(9)},{'initial':np.zeros((16,2))}])
def test_invalid_collision_domains_rejected(kwargs):
    with pytest.raises(ValueError): execute(**kwargs)


def test_rational_interaction_signal_enclosure_independent_of_float_execution():
    from .collision import collision_certificate
    value=collision_certificate()
    lo,hi=map(F,value['exact_empty_probability_interval'])
    target=expm(4j*independent_h(2,0))@expm(-4j*independent_h(2,1))@preparation()
    probability=np.linalg.norm(target[::2])**2
    assert float(lo)<probability<float(hi)
    assert F(value['ideal_word_interaction_signal_lower'])>F(36,1000)
    assert F(value['ideal_word_half_diamond_bound'])==F(11,1024)
