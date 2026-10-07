"""Pressure of the existing golden 64-site scalar action, under declared dilation.

Two distinct readouts: instantaneous mechanical stress from authenticated
classical records, and equilibrium thermodynamics of a supplied finite Bose
state. Neither supplies a native repair EoS or a physical calibration.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
RER = HERE.parents[2]
PARENT = 'code/source_scalar_execution/source_scalar_execution_receipt.json'
SCHEMA = 'oph.golden_scalar_metric_pressure.v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encoded(data):
    return (json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()


def parent_modules(rer):
    sys.path.insert(0, str(rer/'code/source_scalar_execution'))
    import scalar_execution_algebra as alg
    import source_scalar_execution as source
    import verify_source_scalar_execution as independent
    return alg, source, independent


def dense(rows, size=64):
    a = np.zeros((size, size))
    for i, row in enumerate(rows):
        for j, c in row:
            a[i, j] = c
    return a


def pressure_samples(K, G, U):
    # Hold the original field coordinates and momenta p=M*v fixed.
    # These are raw work samples, independent of the exact derivative below.
    eps = 1e-4
    scales = [1-eps, 1, 1+eps]
    energies = [K/s**3+G*s+U*s**3 for s in scales]
    p = -(energies[2]-energies[0])/(scales[2]**3-scales[0]**3)
    return {'scale': scales, 'energy_fixed_q_p': energies, 'pressure_secant': p}


def classical_readout(packet, source, alg):
    parse, bounds, Q = alg.parse, alg.rational_bounds, alg.Q
    masses = [parse(x) for x in packet['mass_Qphi']]
    rows = [[(j, parse(c)) for j, c in row] for row in packet['action_rows_Qphi']]
    reports = {}
    for name, trace in packet['traces'].items():
        layers = [[parse(x) for x in layer] for layer in trace['layers']]
        samples = []
        for n in range(22):
            previous, q = layers[n:n+2]
            force = source.action(rows, q)
            v = [(x-y)/source.TAU-source.TAU*f/2 for x, y, f in zip(q, previous, force)]
            K = source.inner(masses, v, v)/2
            U = source.inner(masses, q, q)/2
            G = source.inner(masses, q, [f-x for f, x in zip(force, q)])/2
            H = K+G+U
            p = K-G/3-U  # V0=1, pressure from fixed canonical q,p work
            modified = H-source.TAU**2*source.inner(masses, force, force)/8
            assert min(K.sign(), G.sign(), U.sign()) >= 0
            assert (H-p).sign() >= 0 and (H+p).sign() >= 0
            w = p/H if H.sign() > 0 else None
            def fl(x):
                lo, hi = bounds(x)
                return float((F(lo)+F(hi))/2)
            samples.append({
                'step': n, 'time_Qphi': (n*source.TAU).encode(),
                'kinetic_Qphi': K.encode(), 'gradient_Qphi': G.encode(),
                'mass_potential_Qphi': U.encode(), 'action_energy_Qphi': H.encode(),
                'modified_energy_Qphi': modified.encode(), 'pressure_Qphi': p.encode(),
                'energy_density_Qphi': H.encode(),
                'w_stress_Qphi': w.encode() if w is not None else None,
                'energy_interval': bounds(H), 'pressure_interval': bounds(p),
                'w_stress_interval': bounds(w) if w is not None else None,
                'work_samples': pressure_samples(fl(K), fl(G), fl(U)),
            })
        assert all(s['modified_energy_Qphi'] == samples[0]['modified_energy_Qphi'] for s in samples)
        reports[name] = samples
    return reports


def _thermal_input(value, name, ndim):
    """Accept finite real data exactly representable in the binary64 interface."""
    if np.ma.is_masked(value):
        raise ValueError(name+' must not contain masked values')
    if isinstance(value, (list, tuple)) and any(
            isinstance(x, (bool, np.bool_)) for x in np.asarray(value, dtype=object).flat):
        raise ValueError(name+' must not contain Booleans')
    raw = np.asarray(value)
    if raw.ndim != ndim or raw.size == 0 or raw.dtype.kind not in 'iuf':
        raise ValueError(name+' requires a nonempty real input of the declared dimension')
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        data = raw.astype(float)
    if not np.all(np.isfinite(data)):
        raise ValueError(name+' must be finite in binary64')
    changed = (any(int(x) != int(y) for x, y in zip(raw.flat, data.flat))
               if raw.dtype.kind in 'iu' else np.any(raw != data))
    if changed:
        raise ValueError(name+' loses information in binary64 conversion')
    return data


def thermal_point(spatial_eigenvalues, scale, mass, temperature):
    """Finite zero-point-subtracted Bose readout in supplied model units.

    Nonnegative spatial eigenvalues and mass, positive scale/temperature and
    strictly positive mode frequencies are required. Compute in a private
    90-digit context, then require every nonzero reported component to survive
    binary64 conversion with relative error <= 1e-12. This is a numerical
    precision check, not an interval certificate or a low-temperature cutoff.
    """
    spectrum = _thermal_input(spatial_eigenvalues, 'spatial eigenvalues', 1)
    scale, mass, temperature = [float(_thermal_input(value, name, 0)) for value, name in
                               ((scale, 'scale'), (mass, 'mass'), (temperature, 'temperature'))]
    if np.any(spectrum < 0) or scale <= 0 or mass < 0 or temperature <= 0:
        raise ValueError('nonnegative spectrum/mass and positive scale/temperature required')
    if mass == 0 and np.any(spectrum == 0):
        raise ValueError('a zero-frequency Bose mode has no finite thermal partition function')
    ctx = mp.mp.clone()
    ctx.dps = 90
    length, rest_mass, temp = map(ctx.mpf, (scale, mass, temperature))
    volume = length**3
    fields = {key: [] for key in ('frequencies', 'occupations', 'mode_thermal_energies',
                                  'mode_K', 'mode_G', 'mode_U')}
    free_modes = []
    for value in spectrum:
        gradient = ctx.mpf(float(value))/length**2
        omega_squared = rest_mass**2 + gradient
        omega = ctx.sqrt(omega_squared)
        x = omega/temp
        if x > 750:  # Even the occupation is then below binary64 range.
            raise ValueError('occupations cannot be resolved in binary64')
        occupation = 1/ctx.expm1(x)
        energy = omega*occupation
        # K=G+U in thermal equilibrium. Thus pV=K-G/3-U=2G/3:
        # keep the small gradient term instead of subtracting rest energies.
        values = (omega, occupation, energy, energy/2,
                  energy*gradient/(2*omega_squared), energy*rest_mass**2/(2*omega_squared))
        for key, result in zip(fields, values):
            fields[key].append(result)
        # Complementary charts preserve both small x and exponentially small F.
        log_partition_inverse = (ctx.log(-ctx.expm1(-x)) if x <= 1
                                 else ctx.log1p(-ctx.exp(-x)))
        free_modes.append(temp*log_partition_inverse)
    energy = ctx.fsum(fields['mode_thermal_energies'])
    pressure_volume = 2*ctx.fsum(fields['mode_G'])/3

    def output(value, label):
        result = float(value)
        if (not np.isfinite(result)
                or (value != 0 and abs(ctx.mpf(result)/value-1) > ctx.mpf('1e-12'))):
            raise ValueError(label+' cannot be resolved in binary64')
        return result

    totals = {'volume': volume, 'energy': energy, 'energy_density': energy/volume,
              'pressure': pressure_volume/volume, 'w': pressure_volume/energy,
              'helmholtz_free_energy': ctx.fsum(free_modes)}
    return {'scale': scale, 'mass': mass, 'temperature': temperature,
            **{key: [output(value, key) for value in values] for key, values in fields.items()},
            **{key: output(value, key) for key, value in totals.items()}}


def build(rer=RER):
    alg, source, independent = parent_modules(rer)
    packet = independent.load(rer/PARENT)
    replay = independent.verify(packet, root=rer)
    def numeric(pair):
        # Individual model coefficients are small; long histories below use
        # exact arithmetic and outward rational intervals instead.
        return float(F(pair[0])) + float(F(pair[1]))*(1+5**.5)/2
    M = np.array([numeric(x) for x in packet['mass_Qphi']])
    A = dense([[(j, numeric(c)) for j, c in row] for row in packet['action_rows_Qphi']])
    K0 = M[:, None]*(A-np.eye(64))
    lam, eigenvectors = eigh(K0, np.diag(M))
    assert lam.min() > 0 and np.max(abs(K0-K0.T)) < 1e-13
    thermal = [thermal_point(lam, L, m, T)
               for L in (1., 2., 8.) for m in (0., 1., 4.) for T in (.5, 2., 10., 50.)]
    histories = classical_readout(packet, source, alg)
    pins = [PARENT, 'code/source_scalar_execution/source_scalar_execution.py',
            'code/source_scalar_execution/scalar_execution_algebra.py',
            'code/source_scalar_execution/verify_source_scalar_execution.py']
    scope = {
        'native_repair_eos_established': False,
        'physical_energy_clock_calibrated': False,
        'thermalization_by_native_repair_demonstrated': False,
        'vacuum_energy_included': False,
        'quantum_state_preparation_demonstrated': False,
        'new_continuum_or_cutoff_removal_claim': False,
        'instantaneous_stress_is_equilibrium_eos': False,
        'new_observer_event_history_executed': False,
        'all_parent_histories_replayed': True,
    }
    return {
        'schema': SCHEMA, 'scope': scope,
        'assumptions': [
            'existing q=5 golden tensor Dirichlet scalar action and prepared observer records',
            'isotropic 3D box dilation with fixed topology: M(L)=L^3 M0 and Kgrad(L)=L K0',
            'pressure is minus Hamiltonian derivative with respect to V=L^3 at fixed original q,p',
            'model mass stays fixed during volume variation; hbar=kB=c=1 are model unit conventions',
            'thermal branch supplies a canonical Bose equilibrium state with zero chemical potential',
            'thermal branch subtracts zero point energy at every L; no cosmological vacuum prediction',
            '64 oscillators fixed when varying volume; no thermodynamic or continuum limit is asserted',
            'finite mode count does not bound bosonic occupation or Hilbert dimension; no finite observer-capacity realization supplied',
            'volume is the full Dirichlet box L^3, not the sum of interior dual masses',
        ],
        'formulas': {
            'Hamiltonian_at_fixed_q_p': 'H(L)=K/L^3+L*G+L^3*U at reference L=1',
            'instantaneous_pressure': 'pV=K-G/3-U; rho=(K+G+U)/V',
            'canonical_velocity': 'v_n=(q_n-q_(n-1))/tau-tau*A*q_n/2',
            'thermal_frequency': 'omega_j(L)^2=m^2+lambda_j/L^2',
            'thermal_free_energy': 'F=T sum_j log(1-exp(-omega_j/T)); zero-point energy subtracted',
            'thermal_pressure': 'p=-partial_V F at fixed T,m,64-mode regulator',
            'thermal_w_identity': 'w=(sum_j E_j lambda_j/(L^2 omega_j^2))/(3 sum_j E_j)',
            'boundary': 'Dirichlet boundary terms retained in K0; scalar rest mass is an action input',
        },
        'parent_replay': replay,
        'spatial_operator': {'mass_diagonal': M.tolist(), 'gradient_stiffness': K0.tolist(),
                             'generalized_eigenvalues': lam.tolist(),
                             'generalized_eigenvectors': eigenvectors.tolist(),
                             'minimum_eigenvalue': float(lam.min()), 'maximum_eigenvalue': float(lam.max())},
        'classical_record_readouts': histories,
        'thermal_cases': thermal,
        'provenance': {'source_sha256': {p: sha(rer/p) for p in pins},
                       'producer_sha256': sha(Path(__file__)), 'python': platform.python_version(),
                       'packages': {n: importlib.metadata.version(n) for n in ('numpy','scipy','mpmath')}},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rer-root', type=Path, default=RER)
    parser.add_argument('--out', type=Path, default=HERE/'scalar_eos_receipt.json')
    args = parser.parse_args()
    receipt = build(args.rer_root)
    args.out.write_bytes(encoded(receipt))
    print(json.dumps({'output': str(args.out), 'thermal_cases': len(receipt['thermal_cases']),
                      'classical_samples': sum(map(len,receipt['classical_record_readouts'].values())),
                      'native_eos': False, 'parent_replay': receipt['parent_replay']}, indent=2))


if __name__ == '__main__':
    main()
