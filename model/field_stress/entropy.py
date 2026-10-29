"""Full-sector Gibbs identities and actual CAR regional modular comparisons.

The reference is a mathematically constructed canonical Gibbs density of the
same interacting Hamiltonian. Its preparation and identification with the
native observer's reference are not asserted. Dense spectral calculations
below are numerical diagnostics; the finite identities are proved in README.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.linalg import eigh
from ..interacting_field.lattice import WeylLattice
from ..field_regions.regions import restrict_density, regional_basis


def entropy(values):
    values = np.asarray(values, float)
    if np.min(values) <= 0:
        raise ArithmeticError('faithful density required by this Gibbs calculation')
    return float(-np.dot(values, np.log(values)))


def regional_hamiltonian(lattice, kept_modes, labels):
    """Bare internal action on the retained 0/1/2-particle regional support.

Independent occupation-basis CAR assembly; no crossing hopping terms are
silently reclassified as internal regional energy.
"""
    modes = tuple(kept_modes)
    positions = {p:i for i,p in enumerate(modes)}
    masks = [sum(1 << positions[p] for p in label) for label in labels]
    index = {mask:i for i,mask in enumerate(masks)}
    single = lattice.single_h[np.ix_(modes,modes)].toarray()
    result = np.zeros((len(masks),len(masks)), complex)
    for column, mask in enumerate(masks):
        for p in range(len(modes)):
            for q in range(len(modes)):
                if not (mask >> q)&1:
                    continue
                middle = mask^(1 << q)
                if (middle >> p)&1:
                    continue
                final = middle|(1 << p)
                parity = ((mask & ((1 << q)-1)).bit_count()
                          +(middle & ((1 << p)-1)).bit_count())
                result[index[final],column] += (-1)**parity*single[p,q]
        for site in range(lattice.volume):
            first, second = 2*site, 2*site+1
            if first in positions and second in positions:
                if mask & (1 << positions[first]) and mask & (1 << positions[second]):
                    result[column,column] += lattice.interaction
    return result


def run(lattice=None):
    lattice = WeylLattice(1) if lattice is None else lattice
    beta = math.log(25/8)
    energies, vectors = eigh(lattice.hamiltonian.toarray(), check_finite=False)
    unnormalized = np.exp(-beta*energies)
    partition = float(sum(unnormalized))
    weights = unnormalized/partition
    reference = (vectors*weights)@vectors.conj().T
    modular = -(vectors*np.log(weights))@vectors.conj().T
    identity = np.eye(lattice.dimension)
    gibbs_residual = float(np.max(abs(modular-beta*lattice.hamiltonian.toarray()-math.log(partition)*identity)))
    energy_reference = float(weights@energies)
    initial = lattice.basis_state(lattice.origin,lattice.origin+1)
    target = np.outer(initial,initial.conj())
    direction = target-reference
    trace_direction = float(np.trace(direction).real)
    energy_direction = float(np.vdot(initial,lattice.hamiltonian@initial).real-energy_reference)
    entropy_derivative = float(np.einsum('ij,ji->',direction,modular-identity).real)
    modular_derivative = float(np.einsum('ij,ji->',direction,modular).real)
    epsilon = 1/64
    changed = reference+epsilon*direction
    changed_values = eigh(changed,eigvals_only=True,check_finite=False)
    entropy_difference = entropy(changed_values)-entropy(weights)
    # Tr rho' log rho' - Tr rho' log rho, without assuming commuting densities.
    relative_entropy = -entropy(changed_values)+float(np.einsum('ij,ji->',changed,modular).real)
    delta_energy = epsilon*energy_direction
    hs_direction_squared = float(np.vdot(direction,direction).real)
    remainder_upper = epsilon**2*hs_direction_squared/(2*(1-epsilon)*float(min(weights)))
    kept = (lattice.origin,lattice.origin+1,lattice.neighbor,lattice.neighbor+1)
    labels = regional_basis(kept)
    regional = restrict_density(lattice,reference,kept)
    region_values, region_vectors = eigh(regional)
    region_modular = -(region_vectors*np.log(region_values))@region_vectors.conj().T
    bare = regional_hamiltonian(lattice,kept,labels)
    number = np.diag([len(label) for label in labels])
    comparison = []
    for names, basis in (
        (['constant','bare_H','number'],[np.eye(len(labels)),bare,number]),
        (['constant','bare_H','number','number_squared'],[np.eye(len(labels)),bare,number,number@number]),
    ):
        design = np.stack([operator.ravel() for operator in basis],axis=1)
        coefficients = np.linalg.lstsq(design,region_modular.ravel(),rcond=None)[0]
        fitted = sum(coefficient*operator for coefficient,operator in zip(coefficients,basis,strict=True))
        comparison.append({'span':names,'coefficients_real':coefficients.real.tolist(),
                           'coefficient_imaginary_max':float(np.max(abs(coefficients.imag))),
                           'Hilbert_Schmidt_residual':float(np.linalg.norm(region_modular-fitted,'fro'))})
    commutator = region_modular@bare-bare@region_modular
    noncommuting = direction@reference-reference@direction
    return {'reference': 'exp(-beta H)/Z on the complete open R=1 N=2 sector',
            'reference_dimension': lattice.dimension,
            'beta': beta,'beta_choice': 'log(25/8), a declared comparison value; no full-field native thermal preparation',
            'reference_minimum_eigenvalue':float(min(weights)),
            'reference_trace_residual':abs(float(np.trace(reference).real)-1),
            'full_Gibbs_modular_max_entry_residual':gibbs_residual,
            'energy_minimum':float(energies[0]),'energy_maximum':float(energies[-1]),
            'reference_energy':energy_reference,
            'variation':'rho_epsilon=(1-epsilon)rho_beta+epsilon |origin contact pair><origin contact pair|',
            'variation_trace_residual':abs(trace_direction),
            'variation_reference_commutator_Hilbert_Schmidt':float(np.linalg.norm(noncommuting,'fro')),
            'energy_derivative':energy_direction,
            'entropy_derivative':entropy_derivative,
            'modular_derivative':modular_derivative,
            'first_law_derivative_residual':abs(entropy_derivative-beta*energy_direction),
            'epsilon_exact':'1/64','entropy_difference':entropy_difference,
            'energy_difference':delta_energy,'quantum_relative_entropy':relative_entropy,
            'exact_identity_numerical_residual':abs(entropy_difference+relative_entropy-beta*delta_energy),
            'analytic_relative_entropy_upper_evaluated_numerically':remainder_upper,
            'regional_modes':list(kept),'regional_dimension':len(labels),
            'regional_basis':[list(label) for label in labels],
            'regional_positive_support':'local particle number 0, 1, or 2; full four-mode Fock algebra is not faithful',
            'regional_minimum_eigenvalue':float(min(region_values)),
            'regional_trace_residual':abs(float(np.trace(regional).real)-1),
            'bare_regional_modular_comparisons':comparison,
            'regional_modular_bare_energy_commutator_Hilbert_Schmidt':float(np.linalg.norm(commutator,'fro')),
            'density_preparation_source_executed':False,
            'physical_clock_or_temperature_calibrated':False,
            'outward_certified_floating_error':False}
