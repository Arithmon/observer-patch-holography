"""Electric response of the same finite particle/hole candidates.

The external Peierls link is supplied. This checks its full-Fock Ward
identity and complete candidate dipole/current response, not a dynamical
photon theory. One open bond has no magnetic flux; its phase cannot test a
magnetic moment. The bare contact-kernel comparison excludes only its direct
identification with a Maxwell Coulomb kernel, not emergent electromagnetism.

All operators use Hamiltonian model coordinate, hbar_model=1. Under a
supplied physical dictionary, Q=e(N-2), D_phys=e*Lunit*D, and
J_phys=e*J/Tunit. No dimensionful small constant is used in residual tests.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction as F
import numpy as np

from model.full_fock_field.field import FieldGeometry, FullFock
from model.gauged_field.links import GaugedField
from model.interacting_field.lattice import real
from model.particle_response import dimer


REFERENCES = {
    "external_Peierls_coupling": "https://quantum.lassp.cornell.edu/node/335",
    "dynamical_QED_and_Coulomb_exchange": "https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S6.html",
}
TOLERANCE = 2e-10


def operators(geometry=FieldGeometry((2, 1, 1)), phase=0.):
    """Lift the existing external link law to every state of the dimer.

    The original GaugedField evaluator retains N=2. Its one-body symbol is
    used here with complete CAR matrices, so the N=1 and N=3 candidates are
    actually present. N=2 agreement is independently checked in tests.
    """
    if (not isinstance(geometry, FieldGeometry) or geometry.volume != 2
            or geometry.speed <= 0):
        raise ValueError('two-site open geometry with positive coordinate speed required')
    phase = real(phase, 'dimensionless external Peierls phase')
    links = GaugedField(geometry.sites, angles={(0, 1): phase},
                       speed=geometry.speed, interaction=geometry.interaction)
    space = FullFock(4)
    ladders = [space.annihilator(i).toarray() for i in range(4)]
    numbers = [a.conj().T@a for a in ladders]
    local = (numbers[0]+numbers[1], numbers[2]+numbers[3])
    h = links.single_h.toarray()
    kinetic = sum(h[i, j]*ladders[i].conj().T@ladders[j]
                  for i in range(4) for j in range(4))
    contact = geometry.interaction*(numbers[0]@numbers[1]+numbers[2]@numbers[3])
    forward = sum(h[i, j]*ladders[i].conj().T@ladders[j]
                  for i in range(2) for j in range(2, 4))
    number = local[0]+local[1]
    return {'hamiltonian': kinetic+contact,
            'generator': kinetic+contact-geometry.interaction*number/2,
            'kinetic': kinetic, 'contact': contact, 'local_numbers': local,
            'charge': number-2*np.eye(16),
            'dipole': (local[1]-local[0])/2,
            'current': 1j*(forward-forward.conj().T),
            'cycle_rank': len(links.axes)-links.volume+1}


def candidate_response(spectrum, projector, dipole, current, *, charge):
    """Complete spectral measure for a declared equally mixed spin doublet.

    The full sixteen-state final basis is retained, including dark and
    forbidden-number transitions. The initial rank-two state is declared;
    no detected response is conditioned or renormalized.
    """
    if type(charge) is not int or charge not in (-1, 1):
        raise ValueError('one of the two assigned charge-one candidates required')
    from model.particle_response.spectrum import matrix
    p = matrix(projector, spectrum.dimension)
    dipole, current = (matrix(op, spectrum.dimension) for op in (dipole, current))
    if any(np.linalg.norm(op-op.conj().T, 2) > 2e-12 for op in (dipole, current)):
        raise ValueError('Hermitian electric observables required')
    values, vectors = np.linalg.eigh(p)
    selected = values > .5
    if (np.count_nonzero(selected) != 2 or np.linalg.norm(p@p-p, 2) > 2e-12
            or np.linalg.norm(p-p.conj().T, 2) > 2e-12):
        raise ValueError('complete rank-two candidate projector required')
    basis = vectors[:, selected]
    rho = basis@basis.conj().T/2
    k = spectrum.generator
    energy = float(np.trace(rho@k).real)
    if (np.linalg.norm(k@basis-energy*basis, 2) > 2e-11
            or np.linalg.norm((spectrum.number-2*np.eye(16))@basis-charge*basis, 2) > 2e-11
            or abs(energy-min(spectrum.energies[spectrum.sectors == charge+2])) > 2e-11):
        raise ValueError('lowest same-energy doublet of the declared charge required')
    final = spectrum.vectors
    dipole_weights = np.sum(abs(final.conj().T@dipole@basis)**2, axis=1)/2
    current_weights = np.sum(abs(final.conj().T@current@basis)**2, axis=1)/2
    frequencies = spectrum.energies-energy
    moment = float(frequencies@dipole_weights)
    commutator = k@dipole-dipole@k
    double_commutator = dipole@commutator-commutator@dipole
    expected_moment = float(np.trace(rho@double_commutator).real/2)
    total_dipole = float(np.sum(dipole_weights))
    residuals = {
        'current_dipole_operator_identity': float(np.linalg.norm(current-1j*commutator, 2)),
        'dipole_total_weight': abs(total_dipole-float(np.trace(rho@dipole@dipole).real)),
        'current_total_weight': abs(float(np.sum(current_weights))-float(np.trace(rho@current@current).real)),
        'first_moment_double_commutator': abs(moment-expected_moment),
        'current_dipole_frequency_identity': float(np.max(abs(current_weights-frequencies**2*dipole_weights))),
        'forbidden_number_dipole_weight': float(np.sum(dipole_weights[spectrum.sectors != charge+2])),
        'forbidden_number_current_weight': float(np.sum(current_weights[spectrum.sectors != charge+2])),
    }
    if not all(math.isfinite(x) and x <= TOLERANCE for x in residuals.values()):
        raise ArithmeticError('complete candidate electric-response identities failed')
    return {'charge_in_supplied_units_e': charge,
            'initial_state': 'equally mixed complete lowest charged doublet',
            'final_states_retained': spectrum.dimension,
            'dipole_total_weight': total_dipole,
            'current_total_weight': float(np.sum(current_weights)),
            'dipole_first_frequency_moment': moment,
            'half_double_commutator_expectation': expected_moment,
            'dipole_expectation': float(np.trace(rho@dipole).real),
            'current_expectation': float(np.trace(rho@current).real),
            'transitions': [{'final_number': int(n), 'model_frequency': float(w),
                             'dipole_weight': float(d), 'current_weight': float(j)}
                            for n,w,d,j in zip(spectrum.sectors, frequencies,
                                             dipole_weights, current_weights, strict=True)],
            'residuals': residuals}


def bare_kernel_comparison(interaction=1., side=8):
    """Compare the actual contact shape to an external static Maxwell target.

    Opposite-spin contact V(r)=U delta_(r,0) has Fourier kernel U. The
    inverse nearest-neighbor lattice Laplacian has kernel 1/qhat^2 away
    from its zero mode. Ratios remove the arbitrary overall coefficient.
    These are bare kernels, not dressed quasiparticle scattering amplitudes.
    """
    u = real(interaction, 'supplied onsite interaction', nonnegative=True)
    if type(side) is not int or not 8 <= side <= 64:
        raise ValueError('declared periodic comparison side from8through64 required')
    labels = np.array([[1, 0, 0], [2, 0, 0], [1, 1, 0]], dtype=int)
    momenta = 2*np.pi*labels/side
    laplacian = 4*np.sum(np.sin(momenta/2)**2, axis=1)
    coulomb = 1/laplacian
    contact = np.full(3, u)
    shape_difference = float(np.max(abs(coulomb/coulomb[0]-1))) if u else None
    return {'comparison_side': side, 'nonzero_momentum_labels': labels.tolist(),
            'contact_fourier_kernel': contact.tolist(),
            'inverse_lattice_Laplacian_kernel': coulomb.tolist(),
            'contact_normalized_shape': [1., 1., 1.] if u else None,
            'Coulomb_normalized_shape': (coulomb/coulomb[0]).tolist(),
            'maximum_normalized_shape_difference': shape_difference,
            'bare_contact_equals_nonzero_Coulomb_kernel': bool(u > 0 and shape_difference <= 2e-12),
            'basis': 'actual onsite opposite-spin term versus supplied static lattice Maxwell kernel',
            'zero_mode_excluded': 'periodic Poisson zero mode requires neutral total source',
            'dressed_or_emergent_photon_no_go': False}


def benchmark(geometry=FieldGeometry((2, 1, 1))):
    spectrum = dimer.setup(geometry)
    ops = operators(geometry)
    k, d, j = ops['generator'], ops['dipole'], ops['current']
    n0, n1 = ops['local_numbers']
    conjugation = dimer.charge_conjugation(spectrum)
    norm = lambda value: float(np.linalg.norm(value, 2))
    residuals = {
        'same_full_Fock_action': norm(ops['hamiltonian']-spectrum.hamiltonian),
        'same_half_filled_generator': norm(k-spectrum.generator),
        'left_continuity': norm(1j*(k@n0-n0@k)+j),
        'right_continuity': norm(1j*(k@n1-n1@k)-j),
        'dipole_current_identity': norm(1j*(k@d-d@k)-j),
        'charge_conservation': norm(k@ops['charge']-ops['charge']@k),
        'charge_conjugate_dipole': norm(conjugation@d@conjugation.conj().T+d),
        'charge_conjugate_current': norm(conjugation@j@conjugation.conj().T+j),
    }
    for phase in (.271, -.983, math.pi/2):
        changed = operators(geometry, phase)
        gauge = np.diag(np.exp(1j*phase*np.diag(n0)))
        residuals[f'pure_gauge_action_{phase}'] = norm(changed['generator']-gauge@k@gauge.conj().T)
        residuals[f'pure_gauge_current_{phase}'] = norm(changed['current']-gauge@j@gauge.conj().T)
        residuals[f'pure_gauge_spectrum_{phase}'] = float(np.max(abs(
            np.linalg.eigvalsh(changed['generator'])-np.linalg.eigvalsh(k))))
    projectors = dimer.candidate_projectors(spectrum)
    candidates = [candidate_response(spectrum, projectors[name], d, j, charge=charge)
                  for name,charge in (('charge_minus_one_hole', -1), ('charge_plus_one_particle', 1))]
    for field in ('dipole_total_weight', 'current_total_weight', 'dipole_first_frequency_moment'):
        residuals['charge_conjugate_'+field] = abs(candidates[0][field]-candidates[1][field])
    if not all(math.isfinite(x) and x <= TOLERANCE for x in residuals.values()):
        raise ArithmeticError('same-candidate electric-response benchmark failed')
    if any(row['current_total_weight'] <= 1e-12 for row in candidates):
        raise ArithmeticError('candidate has no resolved nonzero electric-current response')
    return {'scope': 'SAME_FULL_FOCK_DIMER_EXTERNAL_ELECTRIC_RESPONSE',
            'carrier_dimension': spectrum.dimension, 'sides': list(geometry.sides),
            'chemical_potential': spectrum.chemical_potential,
            'same_reference': 'unique N=2 ground of K=H-(U/2)N; responses of its charged excitations',
            'physical_dictionary_is_supplied': True,
            'normalized_operator_units': {'charge': 'e', 'dipole': 'e*Lunit', 'current': 'e/Tunit',
                                          'frequency': '1/Tunit'},
            'candidate_responses': candidates, 'operator_residuals': residuals,
            'open_graph_cycle_rank': ops['cycle_rank'],
            'single_bond_phase_is_a_magnetic_field': False,
            'magnetic_moment_or_gyromagnetic_ratio': None,
            'magnetic_test_status': 'not_determined: no loop flux or magnetic-field/spin coupling in this benchmark',
            'bare_interaction_comparison': bare_kernel_comparison(geometry.interaction),
            'dynamical_gauge_field_status': 'external link parameters, no photon canonical variables in this finite carrier',
            'QED_identification': 'NOT_ESTABLISHED: external charge response does not supply Maxwell propagation, Coulomb exchange or annihilation photons',
            'physical_continuum_or_emergent_photon_excluded': False,
            'source_histories_reexecuted_here': False,
            'method': 'complete operator/eigensystem response; residual tolerance is not an interval or hardware certificate',
            'references': REFERENCES}


def build():
    """Deterministic retained result, after executing all complete-matrix guards.

    For the declared c=1/8 dimer, each N=1/N=3 doublet has an electric
    transition at 2t=c. D has off-diagonal magnitude1/2, so its total
    weight is1/4, the current weight is c²/4 and the first moment is c/4.
    Raw eigensolver residuals are admission diagnostics, not serialized
    evidence. Every final state contributes to a deterministic grouped row.
    """
    report = benchmark()
    exact = {'dipole_total_weight':F(1,4), 'current_total_weight':F(1,256),
             'dipole_first_frequency_moment':F(1,32)}
    candidates = []
    def rounded(value):
        result = round(float(value),10)
        return 0. if result == 0 else result
    for row in report['candidate_responses']:
        if any(abs(row[key]-float(value)) > 2e-12 for key,value in exact.items()):
            raise ArithmeticError('candidate electric response differs from closed dimer formula')
        grouped = {}
        for transition in row['transitions']:
            key = (transition['final_number'],rounded(transition['model_frequency']))
            group = grouped.setdefault(key, {'final_number':key[0], 'model_frequency':key[1],
                'final_states':0, 'dipole_weight':0., 'current_weight':0.})
            group['final_states'] += 1
            group['dipole_weight'] += transition['dipole_weight']
            group['current_weight'] += transition['current_weight']
        groups = [{**group, 'dipole_weight':rounded(group['dipole_weight']),
                   'current_weight':rounded(group['current_weight'])}
                  for _,group in sorted(grouped.items())]
        candidates.append({'charge_in_supplied_units_e':row['charge_in_supplied_units_e'],
            'initial_state':row['initial_state'], 'final_states_retained':row['final_states_retained'],
            'closed_formula_response':{key:str(value) for key,value in exact.items()},
            'nonzero_electric_transition_model_frequency':'1/8',
            'complete_grouped_transition_display':groups})
    bare = report['bare_interaction_comparison']
    if (bare['bare_contact_equals_nonzero_Coulomb_kernel']
            or abs(bare['maximum_normalized_shape_difference']-1/math.sqrt(2)) > 2e-12):
        raise ArithmeticError('bare contact versus Coulomb shape comparison changed')
    return {'schema':'oph.particle_comparison.electromagnetism.v1',
        'scope':report['scope'], 'carrier_dimension':16,
        'generator':'K=H-N/2; c=1/8,U=1; same half-filled dimer reference',
        'candidate_responses':candidates,
        'full_operator_continuity_and_charge_conjugation_checked':True,
        'all_final_number_sectors_retained':True,
        'numerical_admission_tolerance':'2e-10',
        'display_rounding_digits':10, 'numerical_display_is_interval_certificate':False,
        'normalized_operator_units':report['normalized_operator_units'],
        'physical_dictionary_is_supplied':True,
        'open_graph_cycle_rank':report['open_graph_cycle_rank'],
        'pure_gauge_link_spectral_invariance_checked':True,
        'magnetic_moment_or_gyromagnetic_ratio':None,
        'magnetic_test_status':report['magnetic_test_status'],
        'bare_interaction_comparison':{
            'nonzero_momentum_labels':bare['nonzero_momentum_labels'], 'comparison_side':8,
            'contact_normalized_shape':['1','1','1'],
            'Coulomb_normalized_shape':['1','1-1/sqrt(2)','1/2'],
            'maximum_normalized_shape_difference':'1/sqrt(2)',
            'bare_contact_equals_nonzero_Coulomb_kernel':bare['bare_contact_equals_nonzero_Coulomb_kernel'],
            'dressed_or_emergent_photon_no_go':False},
        'dynamical_gauge_field_status':report['dynamical_gauge_field_status'],
        'QED_identification':report['QED_identification'],
        'physical_continuum_or_emergent_photon_excluded':False,
        'source_histories_reexecuted_here':False, 'references':REFERENCES}


if __name__ == '__main__':
    print(json.dumps(build(), indent=2, allow_nan=False))
