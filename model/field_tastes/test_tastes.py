"""Independent action, physical-sector, complete-outcome and thermal controls."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from copy import deepcopy
import json
import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from ..bound_pair.spectrum import physical_basis,full_matrix
from . import tastes,thermal


@pytest.mark.parametrize('side',[6,10,14])
def test_exact_pocket_partition_retains_every_node_and_exchange(side):
    labels=tastes.pocket_labels(side)
    assert np.array_equal(np.bincount(labels),np.full(8,side**3//8))
    indices=np.array(list(product(range(side),repeat=3)))
    partner=(-indices)%side
    mapping=(partner[:,0]*side+partner[:,1])*side+partner[:,2]
    assert np.array_equal(labels,labels[mapping])
    for i,bits in enumerate(product((0,1),repeat=3)):
        p=np.array(bits)*side//2
        index=(p[0]*side+p[1])*side+p[2]
        assert labels[index]==i


@pytest.mark.parametrize('side',[4,8,True,5,20])
def test_ambiguous_or_out_of_domain_pockets_rejected(side):
    with pytest.raises(ValueError): tastes.pocket_labels(side)


def test_complete_physical_action_commutator_singular_values():
    f,labels=tastes.setup()
    physical=physical_basis(f)
    h=physical.conj().T@full_matrix(f)@physical
    for pockets in ((0,),(0,2,7),tuple(range(8))):
        mask=tastes.mask_for(labels,pockets)
        q=physical.conj().T@(mask[:,None]*physical)
        comm=h@q-q@h
        value=np.linalg.svd(comm,compute_uv=False)[0]**2
        w=F(len(pockets),8)
        assert abs(value-float(w*(1-w)))<3e-13
        assert np.max(abs(q@q-q))<2e-14
    # The quarter-domain has both positive- and negative-chirality nodes.
    assert tastes.commutator_certificate(pockets=(0,1))['contact_weight']=='1/4'


def test_complete_time_outputs_against_independent_dense_physical_action():
    f,labels=tastes.setup()
    start,mask,w=tastes.initial(f,labels)
    physical=physical_basis(f)
    h=physical.conj().T@full_matrix(f)@physical
    duration=.5
    independent=physical@(expm(-1j*duration*h)@(physical.conj().T@start))
    changed=expm_multiply(-1j*duration*tastes.momentum_operator(f),start)
    assert np.linalg.norm(independent-changed)<5e-13
    outputs=[changed*tastes.mask_for(labels,(i,)) for i in range(8)]
    assert np.linalg.norm(sum(outputs)-changed)<1e-15
    probs=tastes.complete_probabilities(changed,labels)
    assert abs(sum(probs)-1)<2e-13
    assert min(probs)>0
    # A normalized retained pocket would replace the complete law by delta0.
    assert sum(probs[1:])>.02
    assert f.exchange_residual(changed)<2e-13


def test_free_control_phase_covariance_and_pocket_coarse_graining():
    f,labels=tastes.setup()
    start,mask,w=tastes.initial(f,labels)
    free=expm_multiply(-1j*tastes.momentum_operator(f,0.),start)
    assert sum(tastes.complete_probabilities(free,labels)[1:])==0
    changed=expm_multiply(-1j*tastes.momentum_operator(f),start)
    phase=expm_multiply(-1j*tastes.momentum_operator(f),np.exp(.37j)*start)
    assert np.max(abs(np.array(tastes.complete_probabilities(phase,labels))-
                         tastes.complete_probabilities(changed,labels)))<2e-13
    union=tastes.mask_for(labels,(0,1,2))
    measured=float(np.vdot(union*changed,union*changed).real)
    assert abs(measured-sum(tastes.complete_probabilities(changed,labels)[:3]))<2e-13


def test_complete_pocket_instrument_retains_an_untouched_complex_reference():
    f,labels=tastes.setup()
    first,_,_=tastes.initial(f,labels,(0,))
    second,_,_=tastes.initial(f,labels,(0,1))
    # Joint pure vector: columns label the reference, with nonzero complex coherence.
    joint=np.column_stack((first,1j*second))/np.sqrt(2)
    initial_reference=joint.conj().T@joint
    changed=expm_multiply(-.5j*tastes.momentum_operator(f),joint)
    outputs=[tastes.mask_for(labels,(i,))[:,None]*changed for i in range(8)]
    reduced=sum(branch.conj().T@branch for branch in outputs)
    assert abs(initial_reference[0,1])>.2
    assert np.linalg.norm(reduced-initial_reference)<2e-13
    assert abs(sum(float(np.vdot(branch,branch).real) for branch in outputs)-1)<2e-13
    assert np.linalg.norm(sum(outputs)-changed)<2e-13


def test_positive_exact_short_time_lower_bound_and_derivative():
    f,labels=tastes.setup()
    start,mask,w=tastes.initial(f,labels)
    bound=tastes.leakage_bound()
    duration=float(F(bound['model_time']))
    changed=expm_multiply(-1j*duration*tastes.momentum_operator(f),start)
    outside=float(np.linalg.norm((~mask)*changed)**2)
    assert 0<float(F(bound['unconditioned_outside_probability_lower']))<outside
    assert outside<float(F(bound['unconditioned_outside_probability_upper']))
    tiny=1e-4
    derivative=expm_multiply(-1j*tiny*tastes.momentum_operator(f),start)
    assert abs(np.linalg.norm((~mask)*derivative)**2/tiny**2-7/64)<1e-8


def test_positive_thermal_gram_and_nonnormal_KMS_pair():
    data=thermal.report()
    assert data['complete_physical_K0_two_fermion_dimension']==424
    assert data['smallest_Gibbs_eigenvalue']>0
    gram=np.array(data['reflection_Gram_real'])+1j*np.array(data['reflection_Gram_imag'])
    assert np.linalg.norm(gram-gram.conj().T)<2e-14
    assert np.min(np.linalg.eigvalsh(gram))>0
    # Independent complex test vector, retaining cross-operator imaginary parts.
    z=np.array([1+2j,-.4j,-3.])
    assert np.vdot(z,gram@z).real>0
    assert data['independent_trace_formula_residual']<2e-13
    assert data['KMS_residual']<2e-13
    assert data['wrong_tracial_reference_KMS_residual']>1e-4


def test_retained_source_custody_and_no_chiral_promotion():
    from .experiment import source_bindings
    path=Path(__file__).with_name('results.json')
    value=json.loads(path.read_text())
    assert value['source_bindings']==source_bindings()
    assert value['all_pass'] and all(value['checks'].values())
    assert value['integrated_gate_accepted'] is False
    assert value['physical_comparison_count_added']==0
    assert value['finite_thermal_positivity']['spatial_continuum_OS_reconstruction_or_chiral_QFT_claimed'] is False


def test_receipt_cannot_hide_lost_tastes_or_promote_a_physical_claim():
    from ..observer.scalar_response.experiment import compare_receipt
    retained=json.loads(Path(__file__).with_name('results.json').read_text())
    altered=deepcopy(retained)
    altered['taste_interaction']['complete_evolutions'][0]['all_eight_pocket_probabilities']=[1.]
    with pytest.raises(ValueError): compare_receipt(retained,altered)
    altered=deepcopy(retained)
    altered['finite_thermal_positivity']['spatial_continuum_OS_reconstruction_or_chiral_QFT_claimed']=True
    with pytest.raises(ValueError): compare_receipt(retained,altered)
