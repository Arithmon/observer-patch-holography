"""Complete four-site periodic-field diagnostic with actual translation sectors.

This explicitly changes the open dimer's boundary and volume. Its coordinate
dispersion tests a supplied lattice-Dirac benchmark, not continuum Lorentz
invariance or a physical particle assignment.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import nnls

from ..full_fock_field.field import FieldGeometry, FullFock
from ..reference_observation.observation import frozen


ROOTS = (1., 1j, -1., -1j)
TOLERANCE = 2e-10


def translation():
    """T c†_(x,a) T* = c†_(x+1,a), including fermionic permutation signs."""
    result = np.zeros((256, 256), complex)
    for mask in range(256):
        moved = [2*((j//2+1) % 4)+j % 2 for j in range(8) if mask >> j & 1]
        inversions = sum(a > b for i, a in enumerate(moved) for b in moved[i+1:])
        result[sum(1 << j for j in moved), mask] = (-1)**inversions
    return frozen(result)


def _momentum_basis(addresses, shift, momentum):
    """One projected orbit vector per translation orbit, without an eigencut."""
    remaining = set(int(x) for x in addresses)
    columns = []
    while remaining:
        seed = min(remaining)
        vector = np.zeros(256, complex); current = np.zeros(256, complex)
        current[seed] = 1
        for power in range(4):
            vector += ROOTS[momentum]**(-power)*current/4
            remaining.difference_update(np.flatnonzero(current))
            current = shift@current
        if np.count_nonzero(vector):
            columns.append(vector/np.linalg.norm(vector))
    return np.column_stack(columns) if columns else np.zeros((256, 0), complex)


class PeriodicSpectrum:
    """All 256 states of the same z-Weyl/contact action on a four-site ring."""
    def __init__(self, speed=.125, interaction=1.):
        self.geometry = FieldGeometry((1, 1, 4), speed, interaction)
        if self.geometry.speed <= 0:
            raise ValueError('positive hopping required')
        self.space = FullFock(8)
        self.annihilators = tuple(self.space.annihilator(j).toarray() for j in range(8))
        self.shift = translation()
        up = np.array([sum(bool(mask >> j & 1) for j in range(0, 8, 2)) for mask in range(256)])
        down = np.array([mask.bit_count() for mask in range(256)])-up
        self.up, self.down = up, down
        self.number = np.diag(up+down)
        h = self.space.hamiltonian(self.geometry).toarray()
        single = self.geometry.single_h().toarray()
        # A distinct closing bond from site3 to site0; no duplicate short ring.
        for spin, sign in enumerate((1, -1)):
            coefficient = 1j*self.geometry.speed*sign/2
            left, right = self.annihilators[6+spin], self.annihilators[spin]
            h += coefficient*(left.conj().T@right)+coefficient.conjugate()*(right.conj().T@left)
            single[6+spin, spin] += coefficient
            single[spin, 6+spin] += coefficient.conjugate()
        self.hamiltonian, self.single_h = frozen(h), frozen(single)
        self.generator = frozen(h-self.geometry.interaction*self.number/2)
        sx = np.zeros((256, 256), complex); sy = sx.copy()
        for site in range(4):
            raising = self.annihilators[2*site].conj().T@self.annihilators[2*site+1]
            sx += (-1)**site*(raising+raising.conj().T)/2
            sy += (-1)**site*(raising-raising.conj().T)/(2j)
        sz = np.diag((up-down)/2)
        self.spin = tuple(frozen(s) for s in (sx, sy, sz))
        self.casimir = frozen(sum(s@s for s in self.spin))
        # For odd N=3,5 the only spins are 1/2,3/2. This is the exact
        # spin-half projector on those number sectors, not on the full algebra.
        self.odd_spin_half = frozen((3.75*np.eye(256)-self.casimir)/3)
        if any(not np.array_equal(s@self.generator, self.generator@s) for s in self.spin):
            raise ArithmeticError('periodic untwisted internal SU(2) was lost')
        rows = []
        basis_error = 0.
        for nu in range(5):
            for nd in range(5):
                addresses = np.flatnonzero((up == nu)&(down == nd))
                for momentum in range(4):
                    basis = _momentum_basis(addresses, self.shift, momentum)
                    if not basis.shape[1]:
                        continue
                    basis_error = max(basis_error,
                        float(np.linalg.norm(basis.conj().T@basis-np.eye(basis.shape[1]), 2)),
                        float(np.linalg.norm(self.shift@basis-ROOTS[momentum]*basis, 2)))
                    values, vectors = np.linalg.eigh(basis.conj().T@self.generator@basis)
                    for energy, state in zip(values, (basis@vectors).T, strict=True):
                        rows.append((float(energy), nu, nd, momentum, state))
        rows.sort(key=lambda x: x[:4])
        if len(rows) != 256:
            raise ArithmeticError('translation/number decomposition lost field states')
        self.energies = frozen(np.array([x[0] for x in rows]))
        self.labels = tuple(tuple(x[1:4]) for x in rows)
        self.vectors = frozen(np.column_stack([x[4] for x in rows]))
        self.ground_degeneracy = int(np.count_nonzero(abs(self.energies-self.energies[0]) < TOLERANCE))
        self.half_filled_ground_degeneracy = sum(
            nu+nd == 4 and abs(e-self.energies[0]) < TOLERANCE
            for e, (nu, nd, _) in zip(self.energies, self.labels, strict=True))
        self.numerical_residual = max(basis_error,
            float(np.linalg.norm(self.vectors.conj().T@self.vectors-np.eye(256), 2)),
            float(np.linalg.norm(self.generator@self.vectors-self.vectors*self.energies, 2)))
        if not math.isfinite(self.numerical_residual) or self.numerical_residual > TOLERANCE:
            raise ArithmeticError('complete periodic eigensystem failed admission')
        if not np.array_equal(self.shift@self.generator, self.generator@self.shift):
            raise ArithmeticError('periodic action lost exact translation symmetry')
        if not np.array_equal(self.number@self.generator, self.generator@self.number):
            raise ArithmeticError('periodic action lost exact CAR charge')

    def ground(self):
        if self.ground_degeneracy != 1 or self.labels[0][:2] != (2, 2):
            raise ValueError('unique half-filled singlet-sector ground required')
        if np.linalg.norm(self.casimir@self.vectors[:, 0]) > TOLERANCE:
            raise ArithmeticError('ground is not the declared internal singlet')
        return self.vectors[:, 0]

    def probe(self, momentum, spin, charge=1):
        if (type(momentum) is not int or not 0 <= momentum < 4
                or type(spin) is not int or spin not in (0, 1)
                or type(charge) is not int or charge not in (-1, 1)):
            raise ValueError('momentum index0..3, spin component0..1 and charge +/-1 required')
        # This convention gives T O†_k T*=exp(+ik)O†_k and h(k)=+c sigma_z sin(k).
        op = sum(ROOTS[momentum]**(-x)*self.annihilators[2*x+spin].conj().T for x in range(4))/2
        return frozen(op if charge == 1 else op.conj().T)

    def thermal_density(self, beta=1.):
        if type(beta) not in (int, float) or not math.isfinite(beta) or beta <= 0:
            raise ValueError('finite positive beta required')
        weights = np.exp(-beta*(self.energies-self.energies[0])); weights /= weights.sum()
        if weights.min() < 1e-12:
            raise ValueError('thermal diagnostic exceeds numerical faithfulness domain')
        return frozen((self.vectors*weights)@self.vectors.conj().T)

    def response(self, momentum, spin, charge):
        ground = self.ground(); op = self.probe(momentum, spin, charge)
        changed = op@ground
        if np.linalg.norm(self.casimir@changed-.75*changed) > TOLERANCE:
            raise ArithmeticError('single-fermion probe does not create the declared internal doublet')
        weights = abs(self.vectors.conj().T@changed)**2
        expected = (2+charge*(spin == 0), 2+charge*(spin == 1),
                    (self.labels[0][2]+charge*momentum) % 4)
        allowed = np.array([label == expected for label in self.labels])
        outside = float(weights[~allowed].sum())
        gaps = self.energies-self.energies[0]
        if outside > TOLERANCE:
            raise ArithmeticError('probe escaped its number/spin/translation sector')
        rows = []
        half_weights = np.real(np.einsum('ij,ij->j', self.vectors.conj(), self.odd_spin_half@self.vectors))
        for gap, weight, half in zip(gaps[allowed], weights[allowed], half_weights[allowed], strict=True):
            if not rows or gap-rows[-1]['minimum_gap'] > TOLERANCE:
                rows.append(dict(minimum_gap=float(gap), maximum_gap=float(gap), weight=0.,
                                 multiplicity=0, spin_half_multiplicity=0.))
            rows[-1]['maximum_gap'] = float(gap)
            rows[-1]['weight'] += float(weight)
            rows[-1]['multiplicity'] += 1
            rows[-1]['spin_half_multiplicity'] += float(half)
        norm = float(np.vdot(changed, changed).real)
        moment = float(np.vdot(changed, (self.generator-self.energies[0]*np.eye(256))@changed).real)
        residual = max(abs(float(weights.sum())-norm), abs(float(weights@gaps)-moment))
        spin_half_rows = [row for row in rows if row['spin_half_multiplicity'] > .5]
        if (residual > TOLERANCE or not spin_half_rows or spin_half_rows[0]['weight'] < 1e-8
                or spin_half_rows[0]['minimum_gap'] <= TOLERANCE):
            raise ArithmeticError('lowest momentum-sector pole is unresolved or dark')
        # Complete arrays survive; no small amplitudes or inaccessible groups are deleted.
        return dict(momentum=momentum, spin=spin, charge=charge,
                    excitation_momentum_transfer=(charge*momentum) % 4,
                    sector=expected, poles=rows, full_weights=frozen(weights),
                    lowest_spin_half_gap=spin_half_rows[0]['minimum_gap'],
                    lowest_probe_residue=spin_half_rows[0]['weight'],
                    gaps=frozen(gaps), outside_sector_weight=outside,
                    zeroth_moment=norm, first_moment=moment, moment_residual=residual)


def _rounded(value):
    """Stable retained diagnostics; tolerances are checked before rounding."""
    if isinstance(value, dict):
        return {key: _rounded(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_rounded(item) for item in value]
    if isinstance(value, (float, np.floating)):
        rounded = round(float(value), 10)
        return 0. if rounded == 0 else rounded
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def build():
    spectrum = PeriodicSpectrum(); ground = spectrum.ground()
    rows, maxima = [], [spectrum.numerical_residual]
    by_key = {}
    for momentum in range(4):
        for spin in range(2):
            pair = []
            for charge in (1, -1):
                response = spectrum.response(momentum, spin, charge)
                pair.append(response)
                maxima += [response['moment_residual'], response['outside_sector_weight']]
                row = {key: value for key, value in response.items()
                       if key not in ('full_weights', 'gaps', 'moment_residual', 'outside_sector_weight')}
                rows.append(row); by_key[momentum, spin, charge] = response
            maxima.append(abs(pair[0]['zeroth_moment']+pair[1]['zeroth_moment']-1))
    # Explicit charge conjugation, with all field states retained.
    conjugation = np.eye(256, dtype=complex)
    for operator in spectrum.annihilators:
        conjugation = conjugation@(operator+operator.conj().T)
    maxima += [float(np.linalg.norm(conjugation@spectrum.generator@conjugation.conj().T-spectrum.generator, 2)),
               float(np.linalg.norm(conjugation@spectrum.number@conjugation.conj().T-(8*np.eye(256)-spectrum.number), 2))]
    for momentum in range(4):
        for spin in range(2):
            plus = by_key[momentum, spin, 1]['poles']
            minus = by_key[(-momentum) % 4, spin, -1]['poles']
            if len(plus) != len(minus):
                raise ArithmeticError('particle-hole pole catalogue changed cardinality')
            for a, b in zip(plus, minus, strict=True):
                maxima += [abs(a['minimum_gap']-b['minimum_gap']), abs(a['weight']-b['weight'])]
    free = PeriodicSpectrum(interaction=0.)
    free_bands = []
    for momentum, root in enumerate(ROOTS):
        basis = np.zeros((8, 2), complex)
        for x in range(4):
            basis[2*x:2*x+2] = root**(-x)*np.eye(2)/2
        symbol = basis.conj().T@free.single_h@basis
        expected = free.geometry.speed*complex(root).imag*np.diag([1, -1])
        maxima.append(float(np.linalg.norm(symbol-expected, 2)))
        free_bands.append(dict(momentum_index=momentum, signed_energies=np.linalg.eigvalsh(symbol).tolist()))
    if not all(math.isfinite(x) and x < TOLERANCE for x in maxima):
        raise ArithmeticError('periodic spectral or symmetry check failed')
    c = spectrum.geometry.speed
    anchor = math.hypot(spectrum.geometry.interaction, 2*c)/2-c/2
    data = [(complex(ROOTS[momentum]).imag**2,
             by_key[momentum, spin, 1]['lowest_spin_half_gap']**2)
            for momentum in range(4) for spin in range(2)]
    x, y = np.array(data).T
    periodic_rest = math.sqrt(data[0][1])
    coefficient, _ = nnls(np.column_stack([np.ones(len(x)), x]), y)
    envelope = [min(by_key[momentum, spin, 1]['lowest_spin_half_gap'] for spin in range(2))
                for momentum in range(4)]
    envelope_coefficient, _ = nnls(np.column_stack([np.ones(4), [0., 1., 0., 1.]]), np.square(envelope))
    fixed_residual = float(np.max(abs(y-anchor**2-c*c*x)))
    same_volume_residual = float(np.max(abs(y-periodic_rest**2-c*c*x)))
    envelope_residual = float(np.max(abs(np.square(envelope)-periodic_rest**2-c*c*np.array([0., 1., 0., 1.]))))
    return _rounded(dict(
        scope='DECLARED_FOUR_SITE_PERIODIC_WEYL_CONTACT_COORDINATE_DIAGNOSTIC',
        supplied_boundary='periodic four-site z-axis ring; no twist',
        complete_Fock_dimension=256, coordinate_speed=c, onsite_interaction=1., chemical_potential=.5,
        ground_number=4, ground_momentum_index=spectrum.labels[0][2],
        ground_degeneracy=1, ground_energy=spectrum.energies[0],
        reference_is_open_dimer_reference=False,
        translation_convention='T c†_(x,a) T*=c†_(x+1,a); c†_k=sum exp(-ikx)c†_x/2; h(k)=+c sigma_z sin(k)',
        momentum_indices='k=2*pi*j/4; removal at probe index j transfers momentum -j',
        exact_complete_field_translation_and_number_maps=True,
        numerical_eigensystem_moments_symmetry_within_declared_tolerance=True,
        numerical_tolerance=TOLERANCE,
        poles=rows, particle_hole_relation='addition(j,spin)=removal(-j,spin), including all pole residues',
        free_control=dict(signed_bands=free_bands, massless_nodes=[0, 2], doublers_retained=True,
                          grand_canonical_ground_degeneracy=free.ground_degeneracy,
                          half_filled_ground_degeneracy=free.half_filled_ground_degeneracy,
                          reference='state-independent free CAR anticommutator; no unique free vacuum supplied'),
        benchmark=dict(formula='E(k)^2=m^2+c^2 sin(k)^2', dimer_mass_anchor=anchor,
            periodic_zero_momentum_gap=periodic_rest,
            periodic_to_dimer_rest_gap_ratio=periodic_rest/anchor,
            fixed_dimer_anchor_max_squared_residual=fixed_residual,
            same_volume_zero_anchor_max_squared_residual=same_volume_residual,
            lowest_envelope_same_volume_anchor_max_squared_residual=envelope_residual,
            all_spin_resolved_best_nonnegative_intercept=float(coefficient[0]),
            all_spin_resolved_best_nonnegative_slope=float(coefficient[1]),
            all_spin_resolved_best_fit_max_squared_residual=float(np.max(abs(y-coefficient[0]-coefficient[1]*x))),
            lowest_envelope=envelope,
            envelope_best_nonnegative_intercept=float(envelope_coefficient[0]),
            envelope_best_nonnegative_slope=float(envelope_coefficient[1]),
            envelope_fitted_velocity=float(np.sqrt(envelope_coefficient[1])),
            envelope_fitted_slope_over_supplied_c_squared=float(envelope_coefficient[1]/c**2),
            distinct_sin_squared_coordinates=2,
            envelope_fit_is_underdetermined_as_a_dispersion_test=True,
            effective_velocity_independently_recovered=False,
            emergent_Lorentz_with_renormalized_velocity_excluded=False,
            declared_dimer_anchor_and_speed_benchmark_satisfied=fixed_residual <= TOLERANCE,
            same_volume_anchor_and_declared_speed_benchmark_satisfied=same_volume_residual <= TOLERANCE,
            lowest_envelope_same_volume_benchmark_satisfied=envelope_residual <= TOLERANCE),
        source_preparation_or_probe_execution=False, physical_clock_recovered=False,
        physical_Lorentz_dispersion_established=False, whole_OPH_or_continuum_no_go=False))


if __name__ == '__main__':
    import json
    print(json.dumps(build(), indent=2, allow_nan=False))
