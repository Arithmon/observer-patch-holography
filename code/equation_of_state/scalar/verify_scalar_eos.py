"""Independent replay of metric pressure on the declared golden scalar action.

Classical values are checked in the parent verifier's independent Q(sqrt(5))
arithmetic. The thermal check rebuilds a 4x4 one-dimensional operator at 70
decimal digits, takes all 64 tensor eigenvalue sums, and differentiates the
Helmholtz free energy. Thermal tolerances are numerical, not interval proofs.
No producer module is imported.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
RER = HERE.parents[2]
PARENT = 'code/source_scalar_execution/source_scalar_execution_receipt.json'
SCOPE = {
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
ASSUMPTIONS = [
    'existing q=5 golden tensor Dirichlet scalar action and prepared observer records',
    'isotropic 3D box dilation with fixed topology: M(L)=L^3 M0 and Kgrad(L)=L K0',
    'pressure is minus Hamiltonian derivative with respect to V=L^3 at fixed original q,p',
    'model mass stays fixed during volume variation; hbar=kB=c=1 are model unit conventions',
    'thermal branch supplies a canonical Bose equilibrium state with zero chemical potential',
    'thermal branch subtracts zero point energy at every L; no cosmological vacuum prediction',
    '64 oscillators fixed when varying volume; no thermodynamic or continuum limit is asserted',
    'finite mode count does not bound bosonic occupation or Hilbert dimension; no finite observer-capacity realization supplied',
    'volume is the full Dirichlet box L^3, not the sum of interior dual masses',
]
FORMULAS = {
    'Hamiltonian_at_fixed_q_p': 'H(L)=K/L^3+L*G+L^3*U at reference L=1',
    'instantaneous_pressure': 'pV=K-G/3-U; rho=(K+G+U)/V',
    'canonical_velocity': 'v_n=(q_n-q_(n-1))/tau-tau*A*q_n/2',
    'thermal_frequency': 'omega_j(L)^2=m^2+lambda_j/L^2',
    'thermal_free_energy': 'F=T sum_j log(1-exp(-omega_j/T)); zero-point energy subtracted',
    'thermal_pressure': 'p=-partial_V F at fixed T,m,64-mode regulator',
    'thermal_w_identity': 'w=(sum_j E_j lambda_j/(L^2 omega_j^2))/(3 sum_j E_j)',
    'boundary': 'Dirichlet boundary terms retained in K0; scalar rest mass is an action input',
}


def require(condition, label):
    if not condition:
        raise ValueError(label)


def equal(actual, expected, label):
    require(json.dumps(actual, sort_keys=True, allow_nan=False)
            == json.dumps(expected, sort_keys=True, allow_nan=False), label)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    def pairs(items):
        out = {}
        for k, v in items:
            require(k not in out, 'duplicate JSON key')
            out[k] = v
        return out
    def constant(value):
        raise ValueError('nonfinite JSON number: '+value)
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=constant)


@lru_cache(maxsize=2)
def parent_context(rer_string):
    """Cache an independently verified immutable parent, checking its hash later."""
    rer = Path(rer_string)
    name = 'eos_independent_parent_'+hashlib.sha256(str(rer).encode()).hexdigest()[:12]
    spec = importlib.util.spec_from_file_location(
        name, rer/'code/source_scalar_execution/verify_source_scalar_execution.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    packet = module.load(rer/PARENT)
    report = module.verify(packet, root=rer)
    # Cache only after actual full read/write replay, never a producer assertion.
    pins = {p: digest(rer/p) for p in [PARENT, *packet['source_pins']]}
    return module, packet, report, pins


def qfloat(value):
    with mp.workdps(70):
        a, b = value.a, value.b
        return float(mp.mpf(a.numerator)/a.denominator
                     + mp.sqrt(5)*mp.mpf(b.numerator)/b.denominator)


def numeric(actual, expected, label, atol=3e-10, rtol=3e-11):
    require(type(actual) in (int, float) and np.isfinite(actual), label+' finite number')
    require(abs(actual-float(expected)) <= atol+rtol*abs(float(expected)), label)


def enclosure(encoded, exact, algebra, label):
    require(type(encoded) is list and len(encoded) == 2, label+' interval shape')
    try:
        lo, hi = map(Fraction, encoded)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError(label+' interval rationals') from exc
    require(all(type(v) is str for v in encoded)
            and [str(lo), str(hi)] == encoded, label+' interval canonical')
    require(lo <= hi and hi-lo <= Fraction(1, 10**12), label+' interval width')
    require((exact-algebra.R(lo)).sign() >= 0
            and (algebra.R(hi)-exact).sign() >= 0, label+' exact enclosure')


def check_classical(receipt, parent, alg):
    mass = [alg.parse(x) for x in parent['mass_Qphi']]
    rows = [[(j, alg.parse(c)) for j, c in row] for row in parent['action_rows_Qphi']]
    histories = receipt['classical_record_readouts']
    equal(sorted(histories), sorted(parent['traces']), 'four classical histories')
    seen = 0
    max_secant = 0.
    for name, trace in parent['traces'].items():
        raw = histories[name]
        require(type(raw) is list and len(raw) == 22, 'complete 22-state history')
        layers = [[alg.parse(x) for x in row] for row in trace['layers']]
        first_modified = None
        for n, report in enumerate(raw):
            previous, q = layers[n:n+2]
            force = alg.apply(rows, q)
            v = [(x-y)/alg.TAU-alg.TAU*f/2 for x, y, f in zip(q, previous, force)]
            K = alg.inner(mass, v, v)/2
            U = alg.inner(mass, q, q)/2
            G = alg.inner(mass, q, force)/2-U
            H = K+G+U
            pressure = K-G/3-U
            modified = H-alg.TAU**2*alg.inner(mass, force, force)/8
            w = pressure/H if H.sign() else None
            equal(report['step'], n, 'classical step')
            for key, value in {
                'time_Qphi': n*alg.TAU, 'kinetic_Qphi': K, 'gradient_Qphi': G,
                'mass_potential_Qphi': U, 'action_energy_Qphi': H,
                'modified_energy_Qphi': modified, 'pressure_Qphi': pressure,
                'energy_density_Qphi': H,
            }.items():
                equal(report[key], value.encode(), 'exact classical '+key)
            equal(report['w_stress_Qphi'], None if w is None else w.encode(), 'exact stress ratio')
            enclosure(report['energy_interval'], H, alg, 'energy')
            enclosure(report['pressure_interval'], pressure, alg, 'pressure')
            if w is None:
                equal(report['w_stress_interval'], None, 'zero baseline ratio undefined')
            else:
                enclosure(report['w_stress_interval'], w, alg, 'stress ratio')
                require((H-pressure).sign() >= 0 and (H+pressure).sign() >= 0,
                        'instantaneous stress bound')
            if first_modified is None:
                first_modified = modified
            require(modified == first_modified, 'modified energy conservation')
            work = report['work_samples']
            equal(work['scale'], [0.9999, 1, 1.0001], 'fixed-q,p work protocol')
            k, g, u = map(qfloat, (K, G, U))
            energies = [k/s**3+g*s+u*s**3 for s in work['scale']]
            require(len(work['energy_fixed_q_p']) == 3, 'three work samples')
            for got, expected in zip(work['energy_fixed_q_p'], energies):
                numeric(got, expected, 'independent work energy', atol=2e-12, rtol=0)
            secant = -(energies[2]-energies[0])/(1.0001**3-0.9999**3)
            numeric(work['pressure_secant'], secant, 'work derivative', atol=4e-9, rtol=0)
            max_secant = max(max_secant, abs(secant-qfloat(pressure)))
            require(abs(secant-qfloat(pressure)) < 3e-8, 'secant agrees with exact stress')
            seen += 1
    equal(histories['ascending_baseline'], histories['descending_baseline'], 'baseline schedule control')
    equal(histories['ascending_intervention'], histories['descending_intervention'], 'intervention schedule control')
    return seen, max_secant


def tensor_spectrum():
    """Fresh golden orbit and edge assembly; all boundary springs retained."""
    phi = (1+mp.sqrt(5))/2
    points = sorted([b*phi-mp.floor(b*phi) for b in range(5)])+[mp.mpf(1)]
    gaps = [b-a for a, b in zip(points, points[1:])]
    weights = [(a+b)/2 for a, b in zip(gaps, gaps[1:])]
    B = mp.matrix(4)
    for i in range(4):
        B[i, i] = (1/gaps[i]+1/gaps[i+1])/weights[i]
        if i < 3:
            B[i, i+1] = B[i+1, i] = -1/(gaps[i+1]*mp.sqrt(weights[i]*weights[i+1]))
    one = list(mp.eigsy(B, eigvals_only=True))
    return sorted(a+b+c for a, b, c in product(one, repeat=3))


def check_operator(receipt, parent, alg, spectrum):
    op = receipt['spatial_operator']
    mass = np.array([qfloat(alg.parse(x)) for x in parent['mass_Qphi']])
    A = np.zeros((64, 64))
    for i, row in enumerate(parent['action_rows_Qphi']):
        for j, c in row:
            A[i, j] = qfloat(alg.parse(c))
    K = mass[:, None]*(A-np.eye(64))
    for key, expected in [('mass_diagonal', mass), ('gradient_stiffness', K)]:
        actual = np.asarray(op[key], dtype=float)
        require(actual.shape == expected.shape and np.all(np.isfinite(actual)), key+' shape')
        require(np.max(abs(actual-expected)) < 2e-10, key+' independent assembly')
    lam = np.asarray(op['generalized_eigenvalues'], dtype=float)
    require(lam.shape == (64,) and np.all(np.isfinite(lam)), 'full 64-mode spectrum')
    require(np.max(abs(lam-np.array([float(x) for x in spectrum]))) < 3e-9,
            'high precision tensor spectrum')
    numeric(op['minimum_eigenvalue'], spectrum[0], 'minimum eigenvalue')
    numeric(op['maximum_eigenvalue'], spectrum[-1], 'maximum eigenvalue')
    U = np.asarray(op['generalized_eigenvectors'], dtype=float)
    require(U.shape == (64, 64) and np.all(np.isfinite(U)), 'eigenvector matrix')
    residual = np.max(abs(K@U-(mass[:, None]*U)*lam[None, :]))
    orthogonal = np.max(abs(U.T@(mass[:, None]*U)-np.eye(64)))
    require(residual < 2e-8 and orthogonal < 2e-10, 'generalized eigenpair checks')
    return float(residual), float(orthogonal)


def free_energy(spectrum, volume, mass, temperature):
    length_squared = volume**(mp.mpf(2)/3)
    return temperature*mp.fsum(mp.log(-mp.expm1(-mp.sqrt(mass*mass+x/length_squared)/temperature))
                               for x in spectrum)


def check_thermal(receipt, spectrum):
    cases = receipt['thermal_cases']
    expected_cases = list(product((1., 2., 8.), (0., 1., 4.), (.5, 2., 10., 50.)))
    require(type(cases) is list and len(cases) == len(expected_cases), 'complete 36-case sweep')
    max_pressure = mp.mpf(0)
    min_w, max_w = mp.mpf(1), mp.mpf(0)
    derivative_samples = []
    for case, config in zip(cases, expected_cases):
        equal([case[k] for k in ('scale', 'mass', 'temperature')], list(config), 'frozen case order')
        scale, mass, temp = map(mp.mpf, config)
        volume = scale**3
        omega = [mp.sqrt(mass**2+x/scale**2) for x in spectrum]
        occupations = [1/mp.expm1(w/temp) for w in omega]
        energies = [w*n for w, n in zip(omega, occupations)]
        K = [e/2 for e in energies]
        G = [e*x/(2*scale**2*w**2) for e, x, w in zip(energies, spectrum, omega)]
        U = [e*mass**2/(2*w**2) for e, w in zip(energies, omega)]
        energy = mp.fsum(energies)
        free = free_energy(spectrum, volume, mass, temp)
        # Independent high precision differentiation of the partition sum.
        derivative = -mp.diff(lambda v: free_energy(spectrum, v, mass, temp), volume)
        pressure = mp.fsum(k-g/3-u for k, g, u in zip(K, G, U))/volume
        require(abs(derivative-pressure) < mp.mpf('1e-55')*(1+abs(pressure)),
                'mechanical/partition pressure identity')
        h = volume*mp.mpf('1e-6')
        minus = free_energy(spectrum, volume-h, mass, temp)
        plus = free_energy(spectrum, volume+h, mass, temp)
        secant = -(plus-minus)/(2*h)
        require(abs(secant-derivative) < mp.mpf('3e-10')*(1+abs(derivative)),
                'finite-difference pressure check')
        for key, values in [('frequencies', omega), ('occupations', occupations),
                            ('mode_thermal_energies', energies), ('mode_K', K), ('mode_G', G), ('mode_U', U)]:
            require(type(case[key]) is list and len(case[key]) == 64, key+' 64 modes')
            for actual, expected in zip(case[key], values):
                numeric(actual, expected, key, atol=3e-11, rtol=3e-10)
        for key, value in {'volume': volume, 'energy': energy, 'energy_density': energy/volume,
                           'helmholtz_free_energy': free, 'pressure': derivative,
                           'w': derivative*volume/energy}.items():
            numeric(case[key], value, 'thermal '+key, atol=3e-10, rtol=3e-10)
        w = derivative*volume/energy
        require(0 <= w <= mp.mpf(1)/3+mp.mpf('1e-60'), 'thermal pressure range')
        if mass == 0:
            require(abs(w-mp.mpf(1)/3) < mp.mpf('1e-60'), 'massless radiation scaling')
        max_pressure = max(max_pressure, abs(mp.mpf(case['pressure'])-derivative))
        min_w, max_w = min(min_w, w), max(max_w, w)
        derivative_samples.append({
            'scale': float(scale), 'mass': float(mass), 'temperature': float(temp),
            'volume_minus': str(volume-h), 'volume_plus': str(volume+h),
            'free_energy_minus': mp.nstr(minus, 30), 'free_energy_plus': mp.nstr(plus, 30),
            'pressure_secant': mp.nstr(secant, 30), 'pressure_derivative': mp.nstr(derivative, 30),
        })
    return float(max_pressure), float(min_w), float(max_w), derivative_samples


def verify(receipt, rer=RER):
    require(type(receipt) is dict, 'receipt object')
    equal(receipt['schema'], 'oph.golden_scalar_metric_pressure.v1', 'schema')
    equal(receipt['scope'], SCOPE, 'scientific scope')
    equal(receipt['assumptions'], ASSUMPTIONS, 'supplied assumptions')
    equal(receipt['formulas'], FORMULAS, 'metric conventions')
    alg, parent, report, pins = parent_context(str(Path(rer).resolve()))
    for path, pin in pins.items():
        require(digest(Path(rer)/path) == pin, 'cached parent byte drift')
    equal(receipt['parent_replay'], report, 'independent parent replay')
    expected_pins = {p: digest(Path(rer)/p) for p in (
        PARENT, 'code/source_scalar_execution/source_scalar_execution.py',
        'code/source_scalar_execution/scalar_execution_algebra.py',
        'code/source_scalar_execution/verify_source_scalar_execution.py')}
    equal(receipt['provenance']['source_sha256'], expected_pins, 'source provenance')
    equal(receipt['provenance']['producer_sha256'], digest(HERE/'scalar_eos.py'), 'producer provenance')
    samples, secant_error = check_classical(receipt, parent, alg)
    with mp.workdps(70):
        spectrum = tensor_spectrum()
        residual, orthogonal = check_operator(receipt, parent, alg, spectrum)
        thermal_error, min_w, max_w, differences = check_thermal(receipt, spectrum)
    return {
        'verdict': 'PASS', 'classical_exact_samples': samples,
        'parent_events_replayed': report['events_replayed'],
        'parent_dynamic_reads_replayed': report['dynamic_field_reads_replayed'],
        'classical_max_secant_pressure_error': secant_error,
        'thermal_cases': 36, 'thermal_modes_per_case': 64,
        'independent_spectrum_decimal_digits': 70,
        'thermal_pressure_max_absolute_discrepancy': thermal_error,
        'generalized_eigenpair_max_residual': residual,
        'generalized_orthonormality_max_error': orthogonal,
        'thermal_w_range': [min_w, max_w],
        'thermal_derivative_samples': differences,
        'thermal_numerics_are_rigorous_interval_certificate': False,
        'native_repair_eos_established': False,
        'instantaneous_stress_is_equilibrium_eos': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', type=Path, nargs='?', default=HERE/'scalar_eos_receipt.json')
    parser.add_argument('--rer-root', type=Path, default=RER)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    report = verify(load(args.path), args.rer_root)
    if args.out:
        args.out.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    summary = {k: v for k, v in report.items() if k != 'thermal_derivative_samples'}
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
