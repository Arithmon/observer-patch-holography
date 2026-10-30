"""Particle/hole and pair assignments about one declared half-filled ground state.

These are lattice excitations with actual charge, parity, spin-algebra and
probe-residue maps. Their physical electron/positron interpretation requires
an additional electromagnetic/Lorentz/clock dictionary and mass calibration.
"""
from __future__ import annotations

import math
import numpy as np

from .spectrum import FieldSpectrum
from ..full_fock_field.field import FieldGeometry
from ..interacting_field.lattice import SIGMA


def setup(geometry=FieldGeometry((2, 1, 1))):
    if (not isinstance(geometry, FieldGeometry) or geometry.volume != 2
            or geometry.speed <= 0):
        raise ValueError('a two-site open dimer with positive hopping required')
    return FieldSpectrum(geometry, chemical_potential=geometry.interaction/2, reference='ground')


def spin_generators(spectrum):
    """The dimer's on-site untwisted internal SU(2), not Lorentz rotations."""
    if not isinstance(spectrum, FieldSpectrum) or spectrum.geometry.volume != 2:
        raise ValueError('two-site complete field spectrum required')
    axis = spectrum.geometry.sides.index(2)
    turn = 1j*SIGMA[axis]
    a = [spectrum.space.annihilator(j).toarray() for j in range(4)]
    generators = []
    for sigma in SIGMA:
        value = np.zeros((16, 16), complex)
        for site, local in enumerate((sigma, turn.conj().T@sigma@turn)):
            for i in range(2):
                for j in range(2):
                    value += local[i, j]*(a[2*site+i].conj().T@a[2*site+j])/2
        generators.append(value)
    return tuple(generators)


def charge_conjugation(spectrum):
    """An explicit unitary particle-hole involution up to scalar phase.

    In the untwisted spin frame d_0=c_0,d_1=i sigma_axis c_1, it sends
    d_x to (-1)^x d_x^*. This preserves K=H-U N/2 on the dimer.
    """
    if not isinstance(spectrum, FieldSpectrum) or spectrum.geometry.volume != 2:
        raise ValueError('two-site complete field spectrum required')
    axis = spectrum.geometry.sides.index(2)
    a = [spectrum.space.annihilator(j).toarray() for j in range(4)]
    product = np.eye(16, dtype=complex)
    for operator in a:
        product = product@(operator+operator.conj().T)
    # Product of four Majoranas sends every c_j to -c_j^*. The diagonal
    # phase supplies the required signs in the original (possibly complex) frame.
    signs = (1, 1, -1, -1) if axis == 1 else (1, 1, 1, 1)
    diagonal = np.array([math.prod(-signs[j] for j in range(4) if mask >> j & 1)
                         for mask in range(16)])
    return np.diag(diagonal)@product


def candidate_projectors(spectrum):
    if not isinstance(spectrum, FieldSpectrum) or spectrum.geometry.volume != 2:
        raise ValueError('two-site complete field spectrum required')
    u, t = spectrum.geometry.interaction, spectrum.geometry.speed/2
    if (t <= 0 or spectrum.reference_kind != 'ground'
            or spectrum.chemical_potential != u/2):
        raise ValueError('same unique half-filled ground reference required')
    delta = math.hypot(u, 4*t)
    out = {}
    for name, number, energy in (
            ('charge_minus_one_hole', 1, delta/2-t),
            ('charge_plus_one_particle', 3, delta/2-t),
            ('charge_minus_two_pair_hole', 0, (u+delta)/2),
            ('charge_plus_two_pair_particle', 4, (u+delta)/2)):
        keep = (spectrum.sectors == number)&(abs(spectrum.energies-spectrum.energies[0]-energy) < 1e-10)
        vectors = spectrum.vectors[:, keep]
        expected = 2 if number in (1, 3) else 1
        if vectors.shape[1] != expected:
            raise ArithmeticError('unresolved candidate spectral multiplicity')
        out[name] = vectors@vectors.conj().T
    return out


def certificate(geometry=FieldGeometry((2, 1, 1))):
    spectrum = setup(geometry)
    u, t = geometry.interaction, geometry.speed/2
    delta = math.hypot(u, 4*t)
    gap, pair_gap = delta/2-t, (u+delta)/2
    reference, number = spectrum.reference, spectrum.number
    s = spin_generators(spectrum); casimir = sum(x@x for x in s)
    conjugation = charge_conjugation(spectrum)
    candidates = candidate_projectors(spectrum)
    residuals = [np.linalg.norm(x@spectrum.generator-spectrum.generator@x, 2) for x in s]
    residuals += [np.linalg.norm(s[0]@s[1]-s[1]@s[0]-1j*s[2], 2),
                  np.linalg.norm(casimir@reference, 2),
                  np.linalg.norm(conjugation.conj().T@conjugation-np.eye(16), 2),
                  np.linalg.norm(conjugation@spectrum.generator@conjugation.conj().T-spectrum.generator, 2),
                  np.linalg.norm(conjugation@number@conjugation.conj().T-(4*np.eye(16)-number), 2)]
    parity = np.diag((-1.)**spectrum.numbers)
    rows = []
    for name, projector in candidates.items():
        charge = (-1 if 'minus' in name else 1)*(2 if 'two' in name else 1)
        spin2 = 0. if abs(charge) == 2 else .75
        energy = pair_gap if abs(charge) == 2 else gap
        residuals += [np.linalg.norm((number-2*np.eye(16))@projector-charge*projector, 2),
                      np.linalg.norm(casimir@projector-spin2*projector, 2),
                      np.linalg.norm(parity@projector-(-1)**abs(charge)*projector, 2),
                      np.linalg.norm((spectrum.generator-spectrum.energies[0]*np.eye(16))@projector-energy*projector, 2)]
        op, _ = spectrum.probe('pair' if abs(charge) == 2 else 'single', 0)
        if charge < 0: op = op.conj().T
        residue = float(np.trace(projector@op@reference@op.conj().T).real)
        expected = (1-u/delta)/4 if abs(charge) == 2 else (1+4*t/delta)/4
        residuals.append(abs(residue-expected))
        rows.append(dict(candidate=name,
                         conditional_particle_assignment=('electron-like hole' if charge == -1 else
                            'positron-like particle' if charge == 1 else 'unassigned charge-two lattice pair'),
                         excess_CAR_number=charge,
                         fermion_parity=int((-1)**abs(charge)),
                         internal_SU2_Casimir=spin2, multiplicity=int(round(np.trace(projector).real)),
                         model_frequency_gap=energy, single_local_probe_residue=residue,
                         supplied_speed_gap_over_c_squared=energy/geometry.speed**2,
                         physical_charge_spin_and_rest_mass_identified=False))
    residuals += [np.linalg.norm(conjugation@candidates['charge_minus_one_hole']@conjugation.conj().T
                                 -candidates['charge_plus_one_particle'], 2),
                  np.linalg.norm(conjugation@candidates['charge_minus_two_pair_hole']@conjugation.conj().T
                                 -candidates['charge_plus_two_pair_particle'], 2)]
    maximum = float(max(residuals))
    if not math.isfinite(maximum) or maximum > 2e-11:
        raise ArithmeticError('same-state particle assignment certificate failed')
    return dict(scope='DECLARED_HALF_FILLED_INTERACTING_DIMER_PARTICLE_RESPONSE',
                sides=list(geometry.sides), coordinate_speed=geometry.speed, onsite_interaction=u,
                chemical_potential=u/2, reference='unique N=2 pure ground of K=H-mu*N',
                ground_generator_energy=float(spectrum.energies[0]),
                reference_is_faithful_on_full_field_algebra=False,
                charged_gap_formula='sqrt(U^2+16t^2)/2-t; t=c/2',
                neutral_triplet_gap=(delta-u)/2,
                higher_charge_one_gap=delta/2+t,
                higher_single_local_probe_residue=(1-4*t/delta)/4,
                resolved_charge_one_frequency_splitting=2*t,
                candidates=rows, exact_particle_hole_symmetry_model_mass_ratio=1,
                maximum_matrix_and_residue_residual=maximum,
                conditional_dictionary='q_physical=e*(N-2), supplied internal-to-Lorentz spin map and common calibrated energy unit',
                physical_identification=False, absolute_physical_mass_predicted=False,
                topological_defect_assignment=False, protection_theorem=False,
                physical_clock_calibrated=False,
                Lorentz_dispersion_or_three_dimensional_spin_identified=False)
