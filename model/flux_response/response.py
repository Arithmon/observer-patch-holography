"""Finite loop-flux derivatives, complete spectral sums and gauge covariance.

Flux is a dimensionless supplied Peierls phase. Derivatives use Hamiltonian
model time/energy, not a measured electromagnetic flux or physical current.
The exact identities stated in README hold algebraically; their evaluations
here are floating diagnostics on every state of the 256-dimensional Fock space.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..full_fock_field.field import FullFock
from ..reference_observation.observation import frozen


DIMENSION = 256
TOLERANCE = 5e-11


def _real(value, name, *, nonnegative=False):
    if (type(value) not in (int, float) or not math.isfinite(value)
            or (nonnegative and value < 0)):
        raise ValueError(f'finite {"nonnegative " if nonnegative else ""}{name} required')
    return float(value)


class FluxRing:
    """Same four-site z-Weyl/contact action with one supplied loop flux.

    Bond x points from x to (x+1)%4 and has phase weights[x]*flux.
    The weights sum to one; changing their distribution changes the gauge.
    Zero hopping is admitted for the atomic control, not for a unique-ground
    response. No particle-number or low-energy projection is applied.
    """
    def __init__(self, flux=0., *, speed=.125, interaction=1., weights=(.25,)*4):
        self.flux = _real(flux, 'dimensionless flux')
        self.speed = _real(speed, 'coordinate speed', nonnegative=True)
        self.interaction = _real(interaction, 'onsite interaction', nonnegative=True)
        if type(weights) not in (tuple, list) or len(weights) != 4:
            raise ValueError('four finite flux-distribution weights required')
        self.weights = tuple(_real(w, 'flux-distribution weight') for w in weights)
        if not math.isclose(sum(self.weights), 1., rel_tol=0., abs_tol=1e-14):
            raise ValueError('flux-distribution weights must sum to one')
        if max(abs(self.flux), self.speed, self.interaction, *map(abs, self.weights)) > 100:
            raise ValueError('finite numerical domain bounds every input magnitude by 100')
        space = FullFock(8)
        a = tuple(space.annihilator(j) for j in range(8))
        self.occupations = tuple(frozen(np.diag(
            sum(((space.addresses >> (2*x+s)) & 1) for s in range(2)))) for x in range(4))
        self.number = frozen(sum(self.occupations))
        contact = sum(((space.addresses >> (2*x)) & 1)*
                      ((space.addresses >> (2*x+1)) & 1) for x in range(4))
        h = np.diag(self.interaction*contact).astype(complex)
        first = np.zeros_like(h); second = np.zeros_like(h); currents = []
        for x, weight in enumerate(self.weights):
            directed = sum((1j*self.speed*sign/2)*
                (a[2*x+s].conj().T@a[2*((x+1) % 4)+s])
                for s, sign in enumerate((1, -1))).toarray()
            directed *= np.exp(1j*weight*self.flux)
            bond = directed+directed.conj().T
            current = 1j*(directed-directed.conj().T)
            h += bond; first += weight*current; second -= weight**2*bond
            currents.append(frozen(current))
        self.hamiltonian = frozen(h)
        self.generator = frozen(h-self.interaction*self.number/2)
        self.current = frozen(first)
        self.diamagnetic = frozen(second)
        self.bond_currents = tuple(currents)
        values, vectors = np.linalg.eigh(self.generator)
        self.energies, self.vectors = frozen(values), frozen(vectors)
        self.scale = max(1., float(np.max(abs(values))))
        self.ground_degeneracy = int(np.count_nonzero(values-values[0] <= TOLERANCE*self.scale))
        self.residual = max(float(np.linalg.norm(self.generator@vectors-vectors*values, 2)),
                            float(np.linalg.norm(vectors.conj().T@vectors-np.eye(DIMENSION), 2)))
        if not math.isfinite(self.residual) or self.residual > TOLERANCE*self.scale:
            raise ArithmeticError('complete flux eigensystem failed numerical admission')

    def gauge_to(self, other):
        """G with H_other=G H_self G*: phase redistribution, not flux removal."""
        if (not isinstance(other, FluxRing) or self.flux != other.flux
                or self.speed != other.speed or self.interaction != other.interaction):
            raise ValueError('the same loop flux and local action are required')
        chi = [0.]
        for x in range(3):
            chi.append(chi[-1]+(self.weights[x]-other.weights[x])*self.flux)
        phases = sum(chi[x]*np.diag(self.occupations[x]) for x in range(4))
        return frozen(np.diag(np.exp(1j*phases)))

    def ground_response(self):
        """Complete nondegenerate perturbation response; never divides by a cutoff gap."""
        if self.ground_degeneracy != 1:
            raise ValueError('unique numerically resolved ground required; degeneracy is retained')
        ground = self.vectors[:, 0]
        gaps = self.energies-self.energies[0]
        amplitudes = self.vectors.conj().T@self.current@ground
        weights = abs(amplitudes)**2
        first = float(np.vdot(ground, self.current@ground).real)
        dia = float(np.vdot(ground, self.diamagnetic@ground).real)
        para = float(2*np.sum(weights[1:]/gaps[1:]))
        commutator = self.generator@self.current-self.current@self.generator
        double_commutator = self.current@commutator-commutator@self.current
        moment = float(.5*np.vdot(ground, double_commutator@ground).real)
        norm = float(np.vdot(self.current@ground, self.current@ground).real)
        residual = max(abs(float(weights.sum())-norm), abs(float(weights@gaps)-moment),
                       abs(float(weights[0])-first**2))
        if residual > TOLERANCE*self.scale**3:
            raise ArithmeticError('complete current spectral moment failed')
        return dict(ground_energy=float(self.energies[0]),
                    ground_number=float(np.vdot(ground, self.number@ground).real),
                    first_excitation_gap=float(gaps[1]),
                    flux_derivative=first, diamagnetic_expectation=dia,
                    paramagnetic_sum=para, flux_curvature=dia-para,
                    complete_gaps=frozen(gaps), complete_current_weights=frozen(weights),
                    current_zeroth_moment=norm, current_first_moment=moment,
                    stationary_bond_currents=tuple(float(np.vdot(ground, j@ground).real)
                                                  for j in self.bond_currents),
                    numerical_moment_residual=residual)


def _summary(ring):
    response = ring.ground_response()
    return dict(flux=ring.flux, flux_weights=ring.weights,
                **response,
                retained_spectral_entries=len(response['complete_gaps']))


def report():
    rows = [_summary(FluxRing(flux)) for flux in (0., .4, 1.2)]
    uniform = FluxRing(.4); seam = FluxRing(.4, weights=(0., 0., 0., 1.))
    change = uniform.gauge_to(seam)
    a, b = uniform.ground_response(), seam.ground_response()
    gauge = dict(hamiltonian_residual=float(np.linalg.norm(
        seam.generator-change@uniform.generator@change.conj().T, 2)),
        spectrum_residual=float(np.max(abs(uniform.energies-seam.energies))),
        flux_derivative_difference=abs(a['flux_derivative']-b['flux_derivative']),
        flux_curvature_difference=abs(a['flux_curvature']-b['flux_curvature']),
        uniform_diamagnetic=a['diamagnetic_expectation'],
        seam_diamagnetic=b['diamagnetic_expectation'],
        uniform_paramagnetic=a['paramagnetic_sum'], seam_paramagnetic=b['paramagnetic_sum'])
    free = FluxRing(.4, interaction=0.)
    free_response = free.ground_response()
    alpha = free.flux/4
    free_control = dict(flux=free.flux, ground_degeneracy=free.ground_degeneracy,
        computed_energy=free_response['ground_energy'],
        independent_Bloch_energy=-2*free.speed*(math.cos(alpha)+math.sin(alpha)),
        computed_curvature=free_response['flux_curvature'],
        independent_Bloch_curvature=free.speed*(math.cos(alpha)+math.sin(alpha))/8,
        zero_flux_ground_degeneracy=FluxRing(interaction=0.).ground_degeneracy)
    atom = FluxRing(.4, speed=0.)
    root = Path(__file__).resolve().parents[2]
    source_files = ('model/flux_response/response.py', 'model/flux_response/README.md',
                    'model/full_fock_field/field.py', 'model/reference_observation/observation.py')
    return dict(schema='oph.complete_finite_flux_response.v1',
        source_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                       for name in source_files},
        scope='SUPPLIED_FOUR_SITE_WEYL_CONTACT_LOOP_RESPONSE',
        complete_Fock_dimension=DIMENSION, coordinate_speed=.125, onsite_interaction=1.,
        chemical_potential=.5, flux_is_dimensionless_supplied_Peierls_phase=True,
        interacting_responses=rows, gauge_redistribution=gauge, free_control=free_control,
        atomic_control=dict(ground_degeneracy=atom.ground_degeneracy,
            current_norm=float(np.linalg.norm(atom.current, 2)),
            diamagnetic_norm=float(np.linalg.norm(atom.diamagnetic, 2))),
        reference='unique full-Fock ground of H(phi)-U*N/2 when numerically resolved',
        full_spectrum_and_zero_weight_entries_retained=True,
        exact_identity_numerical_evaluation=True,
        rigorous_numerical_error_bound=False,
        source_pulse_execution_added=False,
        dynamical_photon_or_physical_QED_identified=False,
        physical_clock_or_current_calibrated=False,
        thermodynamic_conductivity_or_Lorentz_limit_proved=False,
        integrated_completion_gate_accepted=False)


def _json_value(value):
    if isinstance(value, dict):
        return {k: _json_value(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if isinstance(value, float):
        return round(value, 11)+0.
    return value


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.write and args.check:
        parser.error('choose write or check')
    result = _json_value(report())
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if args.check:
        if json.loads(args.check.read_text()) != result:
            raise SystemExit('flux response receipt differs from fresh replay')
        print('complete finite flux response replay passed')
    elif args.write:
        with args.write.open('x') as output:
            output.write(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
