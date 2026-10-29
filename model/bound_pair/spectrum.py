"""Complete finite contact spectra and unconditioned two-outcome readouts.

The input is the same Weyl/contact action used by ``fiber``.  Restrict to the
physical antisymmetric fiber before diagonalization.  Spectral weights of
all other eigenstates, including dark zero-overlap states, are retained.
The field-level projectors are supplied mathematical observables; this
module does not claim a source compiler for an entire spatial field.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
import math

import numpy as np
from scipy.linalg import block_diag
from scipy.sparse.linalg import expm_multiply

from ..interacting_field.lattice import real
from .fiber import PairFiber, SINGLET, SWAP, clock_subspace


def physical_basis(fiber):
    """An orthonormal basis of psi(K-k)=-SWAP psi(k), not a sector fit."""
    if not isinstance(fiber, PairFiber) or not fiber.is_finite_torus_momentum:
        raise ValueError('an untwisted finite-torus PairFiber is required')
    if fiber.side > 6:
        raise ValueError('complete dense spectral replay is retained for L<=6')
    side = fiber.side
    total = np.rint(fiber.momentum*side/(2*np.pi)).astype(int)
    indices = tuple(product(range(side), repeat=3))
    columns = []
    fixed = 0
    for i, index in enumerate(indices):
        partner = (total-np.array(index)) % side
        j = int((partner[0]*side+partner[1])*side+partner[2])
        if i > j:
            continue
        if i == j:
            column = np.zeros(fiber.dimension, complex)
            column[4*i:4*i+4] = SINGLET
            columns.append(column)
            fixed += 1
        else:
            for spin in range(4):
                column = np.zeros(fiber.dimension, complex)
                column[4*i+spin] = 1/math.sqrt(2)
                column[4*j:4*j+4] = -SWAP[:,spin]/math.sqrt(2)
                columns.append(column)
    basis = np.column_stack(columns)
    expected = 2*(fiber.volume-fixed)+fixed
    if basis.shape != (fiber.dimension, expected):
        raise ArithmeticError('physical fermion fiber dimension changed')
    return basis


def full_matrix(fiber):
    """Dense full spin/relative-momentum operator, before exchange reduction."""
    if not isinstance(fiber, PairFiber) or fiber.side > 6:
        raise ValueError('retained complete fiber side at most six required')
    return block_diag(*fiber.blocks)+fiber.interaction*np.outer(fiber.d,fiber.d.conj())


def complete_spectrum(fiber):
    """Complete physical spectrum and contact spectral measure.

    Weights are squared matrix elements, without discarding small or zero
    entries.  Degenerate eigenspaces are combined only for the portable
    report; the eigensystem returned here remains complete.
    """
    basis = physical_basis(fiber)
    h = basis.conj().T@full_matrix(fiber)@basis
    energies, vectors = np.linalg.eigh(h)
    contact = basis.conj().T@fiber.d
    coefficients = vectors.conj().T@contact
    weights = abs(coefficients)**2
    if (np.linalg.norm(h-h.conj().T) > 1e-12
            or abs(float(np.sum(weights))-1) > 1e-12):
        raise ArithmeticError('Hermitian physical spectrum or complete weight failed')
    return {'energies':energies,'weights':weights,'vectors':vectors,
            'basis':basis,'hamiltonian':h,'contact':contact,
            'coefficients':coefficients}


def spectral_atoms(spectrum):
    """Combine numerically degenerate energies, retaining every multiplicity.

    This display grouping is a floating diagnostic, not a certified equality
    test.  Time evolution always uses the original ungrouped eigensystem.
    """
    rows = []
    for energy, weight in zip(spectrum['energies'],spectrum['weights'],strict=True):
        if rows and abs(float(energy)-rows[-1]['energy']) < 1e-11:
            row = rows[-1]
            n = row['multiplicity']
            row['energy'] = (n*row['energy']+float(energy))/(n+1)
            row['multiplicity'] = n+1
            row['contact_weight'] += float(weight)
        else:
            rows.append({'energy':float(energy),'multiplicity':1,
                         'contact_weight':float(weight)})
    return rows


def detector_law(fiber, duration, spectrum=None):
    """Contact/complement instrument after the actual finite Hamiltonian.

    The initial state is the normalized contact vector d_K.  Both output
    vectors are unnormalized; no bound-state or detector postselection is
    performed.  The independent reference evolves the unreduced spin fiber.
    """
    duration = real(duration,'Hamiltonian model time')
    if spectrum is None:
        spectrum = complete_spectrum(fiber)
    phases = np.exp(-1j*duration*spectrum['energies'])
    state = spectrum['basis']@(spectrum['vectors']@(phases*spectrum['coefficients']))
    amplitude = np.sum(spectrum['weights']*phases)
    contact = amplitude*fiber.d
    complement = state-contact
    independent = expm_multiply(-1j*duration*full_matrix(fiber),fiber.d)
    probabilities = [float(np.vdot(contact,contact).real),
                     float(np.vdot(complement,complement).real)]
    pole_energy = spectrum['energies'][-1]
    pole_weight = spectrum['weights'][-1]
    pole_only_amplitude = pole_weight*np.exp(-1j*duration*pole_energy)
    return {'duration':duration,'probabilities':probabilities,
            'contact_amplitude':[float(amplitude.real),float(amplitude.imag)],
            'independent_state_residual':float(np.linalg.norm(state-independent)),
            'complete_mass_residual':abs(sum(probabilities)-1),
            'orthogonal_outputs_residual':float(abs(np.vdot(contact,complement))),
            'unrenormalized_single_pole_contact_probability':float(abs(pole_only_amplitude)**2),
            'single_pole_contact_probability_error':abs(probabilities[0]-float(abs(pole_only_amplitude)**2)),
            'postselected_single_pole_predicts_contact_probability':float(pole_weight),
            'postselected_state_is_a_different_preparation':True,
            'bound_filter_complete_probabilities':[float(pole_weight),float(1-pole_weight)]}


def native_clock_attachment():
    """Execute the existing ideal complete native map at its declared s=1.

    This is an independent join check, not a repeat of the expensive source
    word certification retained by particle_clock.  Every rejected branch
    stays in the native output and its original unconditioned mass.
    """
    from ..particle_clock import clock
    field = clock_subspace()
    h = field['restriction']
    beta = math.log(25/8)
    rho = np.diag([25/33,8/33])
    thermal = np.diag(np.exp(-beta*np.diag(h)))
    thermal /= np.trace(thermal)
    q = F(264,7505)
    source_reference = np.diag([float(x) for x in
                               clock.native.reference(clock.native_clock.parameters())]).astype(complex)
    weights = tuple(F(x) for x in clock.native_clock.parameters().port_weights)
    rows = {}
    for law, alpha in (('C',F(1519,1980)),('E',F(37,60))):
        outputs = clock.execute(source_reference,weights,law,raw=False)
        p_r = clock.record_probability(outputs,law,0)
        p_d = clock.record_probability(outputs,law,1)
        initial_outside = float(np.trace(outputs['initial_outside',]).real)
        continued_outside = sum(float(np.trace(v).real) for key,v in outputs.items()
                                if key[0]=='outside')
        abort = sum(float(np.trace(v).real) for key,v in outputs.items() if key[0]=='abort')
        expected = float(q)*(85/198+float(F(17,66)*alpha)*math.cos(beta))
        rows[law] = {'initial_outside':initial_outside,'continued_outside':continued_outside,
                     'energy_r':p_r,'energy_d':p_d,'abort':abort,
                     'full_mass':initial_outside+continued_outside+p_r+p_d+abort,
                     'energy_readout_formula_residual':abs(p_r-expected),
                     'live_and_abort_blocks':len(outputs)}
    return {'field_side':field['side'],'field_energy_gap':float(h[1,1]-h[0,0]),
            'modular_parameter':1,'Hamiltonian_model_time':-beta,
            'two_levels_have_the_same_particle_number':2,
            'two_levels_have_the_same_total_momentum':['pi','pi','pi'],
            'faithful_M2_Gibbs_residual':float(np.linalg.norm(rho-thermal)),
            'ideal_native_complete_readouts':rows,
            'source_word_certificate':'model/particle_clock/results.json',
            'field_isometry_and_macroscopic_detector_are_supplied':True,
            'full_field_Gibbs_preparation_or_clock_calibration_claimed':False}


def deleted_contact_control(fiber):
    """Remove only U P on the same physical fiber and keep its full spectrum."""
    basis=physical_basis(fiber)
    free=basis.conj().T@block_diag(*fiber.blocks)@basis
    energies=np.linalg.eigvalsh(free)
    b=2*math.sqrt(3)*fiber.speed
    return {'same_hopping_and_physical_fiber':True,
            'physical_dimension':len(energies),
            'minimum_energy':float(energies[0]),'maximum_energy':float(energies[-1]),
            'free_norm_bound':b,
            'eigenvalues_above_free_norm_bound':int(np.sum(energies>b+1e-12)),
            'original_interacting_pole_energy':fiber.pole()['energy'],
            'original_bound_pole_survives':False}


def report():
    selected = {'eta':[np.pi]*3,'one_zero':[0.,np.pi,np.pi],'origin':[0.,0.,0.]}
    fibers = []
    for name,momentum in selected.items():
        f = PairFiber(momentum,side=4)
        spectrum = complete_spectrum(f)
        pole = f.pole()
        actual_energy = float(spectrum['energies'][-1])
        actual_weight = float(spectrum['weights'][-1])
        fibers.append({'name':name,'side':4,'pi_coordinates':[float(k/np.pi) for k in momentum],
                       'physical_dimension':len(spectrum['energies']),
                       'all_contact_spectral_atoms':spectral_atoms(spectrum),
                       'all_spectral_multiplicities_retained':True,
                       'bound_pole_energy':actual_energy,'bound_contact_weight':actual_weight,
                       'scattering_contact_weight':float(np.sum(spectrum['weights'][:-1])),
                       'resolvent_energy_residual':abs(actual_energy-pole['energy']),
                       'resolvent_weight_residual':abs(actual_weight-pole['contact_residue']),
                       'energy_first_moment':float(np.dot(spectrum['energies'],spectrum['weights'])),
                       'interaction_deletion_control':deleted_contact_control(f),
                       'readouts':[detector_law(f,t,spectrum) for t in (0.,2.,7.,20.)]})
    volumes = []
    for side in (4,6,8,12,16):
        values = {}
        for name,momentum in selected.items():
            values[name] = PairFiber(momentum,side=side).pole()
        volumes.append({'side':side,'pole_only_volume_diagnostics':values})
    return {'finite_fibers':fibers,'fixed_spacing_volume_diagnostics':volumes,
            'native_clock_join':native_clock_attachment(),
            'scope':{'same_interacting_Weyl_contact_action':True,
                     'contact_and_complement_outputs_both_retained':True,
                     'all_physical_fiber_states_retained_at_L4':True,
                     'scattering_continuum_is_finite_volume_discrete_here':True,
                     'all_state_source_compiled_field_detector':False,
                     'spatial_continuum_or_physical_mass_calibration':False,
                     'Standard_Model_species_assignment':False}}
