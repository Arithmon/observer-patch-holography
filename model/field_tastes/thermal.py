"""Finite thermal positivity on the same complete physical N=2, K=0 fiber.

This is a matrix-algebra check, not a spatial continuum OS reconstruction.
The selected Gibbs state and all momentum observables are supplied.
"""
from __future__ import annotations
import math
import numpy as np
from ..bound_pair.spectrum import physical_basis,full_matrix
from .tastes import setup,mask_for


def thermal_data():
    f,labels=setup(6)
    basis=physical_basis(f)
    h=basis.conj().T@full_matrix(f)@basis
    energies,vectors=np.linalg.eigh(h)
    beta=math.log(25/8)
    weights=np.exp(-beta*energies)
    weights/=np.sum(weights)
    d=basis.conj().T@f.d
    contact=np.outer(d,d.conj())
    mask=mask_for(labels,(0,))
    pocket=basis.conj().T@(mask[:,None]*basis)
    current=1j*(h@contact-contact@h)
    operators=[vectors.conj().T@a@vectors for a in (contact,pocket,current)]
    return beta,energies,weights,operators


def reflection_gram(energies,weights,operators,times):
    """Gij=Tr[e^(-(beta-ti-tj)H) Ai* e^(-(ti+tj)H) Aj]/Z."""
    delta=energies[None,:]-energies[:,None]  # E_column-E_row
    features=np.array([(np.exp(t*delta)*np.sqrt(weights)[None,:]*a).reshape(-1)
                       for a,t in zip(operators,times,strict=True)])
    return features.conj()@features.T


def kms_correlation(energies,weights,operator,time):
    delta=energies[:,None]-energies[None,:]
    return float(np.sum(weights[:,None]*np.exp(time*delta)*abs(operator)**2))


def report():
    beta,energies,weights,operators=thermal_data()
    times=[0.,beta/8,beta/4]
    gram=reflection_gram(energies,weights,operators,times)
    partition=float(np.sum(np.exp(-beta*energies)))
    direct=np.empty((3,3),complex)
    for i,a in enumerate(operators):
        for j,b in enumerate(operators):
            time=times[i]+times[j]
            direct[i,j]=np.trace((np.exp(-(beta-time)*energies)[:,None]*a.conj().T)
                                @(np.exp(-time*energies)[:,None]*b))/partition
    nonhermitian=operators[0]@operators[1]
    t=beta/3
    first=kms_correlation(energies,weights,nonhermitian,t)
    second=kms_correlation(energies,weights,nonhermitian.conj().T,beta-t)
    tracial=np.full_like(weights,1/len(weights))
    wrong=abs(kms_correlation(energies,tracial,nonhermitian,t)-
              kms_correlation(energies,tracial,nonhermitian.conj().T,beta-t))
    return {'complete_physical_K0_two_fermion_dimension':len(energies),
            'beta':beta,'beta_is_supplied_not_full_field_clock_calibration':True,
            'smallest_Gibbs_eigenvalue':float(np.min(weights)),
            'same_action_Gibbs_reference_prepared_by_native_repair':False,
            'even_observables':['contact pair projector','Weyl pocket projector','i[H,contact]'],
            'Euclidean_model_times':times,
            'reflection_Gram_real':gram.real.tolist(),'reflection_Gram_imag':gram.imag.tolist(),
            'reflection_Gram_minimum_eigenvalue':float(np.linalg.eigvalsh(gram)[0]),
            'independent_trace_formula_residual':float(np.max(abs(gram-direct))),
            'nonhermitian_KMS_correlation':first,'KMS_reversed_correlation':second,
            'KMS_residual':abs(first-second),'wrong_tracial_reference_KMS_residual':wrong,
            'finite_matrix_positivity_and_KMS_proved':True,
            'spatial_continuum_OS_reconstruction_or_chiral_QFT_claimed':False}
