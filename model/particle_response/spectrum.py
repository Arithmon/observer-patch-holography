"""All-number Lehmann measures with explicit state, charge and generator.

A visible transition is a response pole of the declared finite lattice model.
No electromagnetic charge, Lorentz spin or vacuum-relative physical mass is
inferred from a CAR charge or a Hamiltonian model-frequency difference.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from ..full_fock_field.field import FieldGeometry, FullFock
from ..interacting_field.lattice import real
from ..reference_observation.observation import frozen


def matrix(value, dimension):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked operators are forbidden')
    raw = np.asarray(value)
    if (raw.shape != (dimension, dimension) or raw.dtype.kind not in 'iufc'
            or not np.isfinite(raw).all()):
        raise ValueError('finite operator on the full declared Fock algebra required')
    return np.asarray(raw, complex)


@dataclass(frozen=True)
class Measure:
    """Every transition, including zero-weight and forbidden-charge entries."""
    frequencies: np.ndarray
    forward: np.ndarray
    reverse: np.ndarray
    allowed: np.ndarray
    charge: int

    def correlation(self, duration, direction='forward'):
        duration = real(duration, 'Hamiltonian model duration')
        if direction not in ('forward', 'reverse', 'anticommutator', 'commutator'):
            raise ValueError('forward, reverse, anticommutator or commutator required')
        weights = {'forward': self.forward, 'reverse': self.reverse,
                   'anticommutator': self.forward+self.reverse,
                   'commutator': self.forward-self.reverse}[direction]
        return complex(np.sum(weights*np.exp(-1j*duration*self.frequencies)))

    def moment(self, order, direction='forward'):
        if type(order) is not int or not 0 <= order <= 4:
            raise ValueError('moment order in0..4 required')
        if direction not in ('forward', 'reverse'):
            raise ValueError('forward or reverse moment required')
        return float(np.sum(getattr(self, direction)*self.frequencies**order))

    def lines(self, tolerance=1e-10, visibility=1e-12):
        """Numerical display grouping; the complete transition arrays survive.

        Grouping is not an exact degeneracy proof. Visibility is a declared
        display threshold, not a detector sensitivity or a justified state cut.
        Zero-weight groups are included, and no weight is renormalized.
        """
        tolerance = real(tolerance, 'display frequency tolerance')
        visibility = real(visibility, 'display visibility')
        if not 0 < tolerance <= 1e-6 or visibility < 0:
            raise ValueError('positive bounded grouping tolerance and nonnegative visibility required')
        rows = sorted(zip(self.frequencies[self.allowed], self.forward[self.allowed],
                          self.reverse[self.allowed], strict=True))
        groups = []
        for frequency, forward, reverse in rows:
            if not groups or frequency-groups[-1]['minimum_frequency'] > tolerance:
                groups.append(dict(minimum_frequency=float(frequency), maximum_frequency=float(frequency),
                                   forward_weight=0., reverse_weight=0., transitions=0))
            row = groups[-1]
            row['maximum_frequency'] = float(frequency)
            row['forward_weight'] += float(forward)
            row['reverse_weight'] += float(reverse)
            row['transitions'] += 1
        for row in groups:
            row['frequency'] = (row['minimum_frequency']+row['maximum_frequency'])/2
            row['visible'] = row['forward_weight']+row['reverse_weight'] > visibility
        return groups


class FieldSpectrum:
    """Complete bounded spectrum, diagonalized in every conserved N sector.

    K=H-mu*N is the declared evolution/reference generator. Mu is a supplied
    chemical potential, not a fitted rest-energy subtraction. Gibbs beta>0
    or a unique pure ground reference is explicit in every instance.
    """
    def __init__(self, geometry=FieldGeometry((2, 1, 1)), beta=1., chemical_potential=0.,
                 reference='gibbs'):
        if not isinstance(geometry, FieldGeometry) or geometry.modes > 8:
            raise ValueError('complete spectral evaluator requires at most8field modes')
        self.geometry, self.space = geometry, FullFock(geometry.modes)
        self.dimension = self.space.dimension
        self.beta = real(beta, 'declared inverse temperature')
        self.chemical_potential = real(chemical_potential, 'declared chemical potential')
        if self.beta <= 0 or reference not in ('gibbs', 'ground'):
            raise ValueError('positive beta and Gibbs or unique-ground reference required')
        self.reference_kind = reference
        self.numbers = np.array([i.bit_count() for i in range(self.dimension)])
        self.number = np.diag(self.numbers)
        self.hamiltonian = frozen(self.space.hamiltonian(geometry).toarray())
        self.generator = frozen(self.hamiltonian-self.chemical_potential*self.number)
        if not np.isfinite(self.generator).all():
            raise ValueError('finite full-Fock generator required')
        if not np.array_equal(self.generator@self.number, self.number@self.generator):
            raise ArithmeticError('declared generator fails exact particle-number conservation')
        rows = []
        for number in range(geometry.modes+1):
            addresses = np.flatnonzero(self.numbers == number)
            energies, vectors = np.linalg.eigh(self.generator[np.ix_(addresses, addresses)])
            for energy, local in zip(energies, vectors.T, strict=True):
                full = np.zeros(self.dimension, complex); full[addresses] = local
                rows.append((float(energy), number, full))
        rows.sort(key=lambda row: (row[0], row[1]))
        self.energies = frozen(np.array([row[0] for row in rows]))
        self.sectors = frozen(np.array([row[1] for row in rows], dtype=int))
        self.vectors = frozen(np.column_stack([row[2] for row in rows]))
        if reference == 'gibbs':
            exponent = -self.beta*(self.energies-self.energies[0])
            if not np.isfinite(exponent).all() or exponent.min() < -600:
                raise ValueError('reference exceeds the faithful floating domain')
            weights = np.exp(exponent); weights /= weights.sum()
            if weights.min() < 1e-12:
                raise ValueError('reference lacks the declared numerical margin')
        else:
            if self.energies[1]-self.energies[0] <= 1e-10:
                raise ValueError('ground reference requires a unique resolved ground state')
            weights = np.zeros(self.dimension); weights[0] = 1.
        self.weights = frozen(weights)
        self.reference = frozen((self.vectors*weights)@self.vectors.conj().T)
        self.reference_root = frozen((self.vectors*np.sqrt(weights))@self.vectors.conj().T)
        residual = max(np.linalg.norm(self.vectors.conj().T@self.vectors-np.eye(self.dimension), 2),
                       np.linalg.norm(self.generator@self.vectors-self.vectors*self.energies, 2),
                       abs(np.trace(self.reference)-1))
        if not math.isfinite(float(residual)) or residual > 2e-11:
            raise ArithmeticError('complete sector eigensystem failed numerical admission')

    def probe(self, kind, index):
        if type(index) is not int or type(kind) is not str:
            raise ValueError('a named probe and integer index are required')
        if kind == 'single' and 0 <= index < self.geometry.modes:
            return frozen(self.space.annihilator(index).conj().T.toarray()), 1
        if kind == 'pair' and 0 <= index < self.geometry.volume:
            a, b = (self.space.annihilator(2*index+j) for j in (0, 1))
            return frozen((a.conj().T@b.conj().T).toarray()), 2
        raise ValueError('single-mode or onsite-pair probe index outside its domain')

    def measure(self, operator, charge):
        operator = matrix(operator, self.dimension)
        if type(charge) is not int or not 1 <= charge <= self.geometry.modes:
            raise ValueError('positive integer creation charge required')
        if np.linalg.norm(self.number@operator-operator@self.number-charge*operator, 2) > 1e-12:
            raise ValueError('probe does not have the claimed CAR number charge')
        amplitudes = self.vectors.conj().T@operator@self.vectors
        strength = np.abs(amplitudes)**2
        frequencies = self.energies[:, None]-self.energies[None, :]
        allowed = self.sectors[:, None]-self.sectors[None, :] == charge
        # Charge blocks are exact zeros because sector vectors have exact zero
        # support outside their own N sector. No amplitudes are projected away.
        if np.count_nonzero(amplitudes[~allowed]):
            raise ArithmeticError('forbidden charge transition survived exact sector support')
        forward, reverse = strength*self.weights[None, :], strength*self.weights[:, None]
        if not np.isfinite(forward+reverse).all():
            raise ArithmeticError('nonfinite transition weights')
        return Measure(frozen(frequencies), frozen(forward), frozen(reverse),
                       frozen(allowed), charge)

    def certify(self, kind='single', index=0):
        operator, charge = self.probe(kind, index)
        response = self.measure(operator, charge)
        commutator = self.generator@operator-operator@self.generator
        expected = [np.trace(self.reference@operator.conj().T@operator).real,
                    np.trace(self.reference@operator@operator.conj().T).real,
                    np.trace(self.reference@operator.conj().T@commutator).real,
                    np.trace(self.reference@commutator@operator.conj().T).real]
        actual = [response.moment(0), response.moment(0, 'reverse'),
                  response.moment(1), response.moment(1, 'reverse')]
        residual = max(abs(a-b) for a, b in zip(actual, expected, strict=True))
        balance = (float(np.max(abs(response.reverse-np.exp(-self.beta*response.frequencies)*response.forward)))
                   if self.reference_kind == 'gibbs' else None)
        if not math.isfinite(residual) or residual > 2e-11 or (balance is not None and balance > 2e-11):
            raise ArithmeticError('complete spectral sum rule or detailed balance failed')
        return dict(probe=kind, index=index, CAR_charge=charge,
                    full_Fock_dimension=self.dimension, reference=self.reference_kind,
                    chemical_potential=self.chemical_potential, beta=self.beta,
                    forward_weight=actual[0], reverse_weight=actual[1],
                    first_forward_moment=actual[2], first_reverse_moment=actual[3],
                    maximum_moment_residual=float(residual), detailed_balance_residual=balance,
                    lines=response.lines(), physical_particle_identification=False)
