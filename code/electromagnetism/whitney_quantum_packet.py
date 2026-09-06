"""Full interacting Coulomb phase-space preparation and neutral Gaussian seeds.

Two classical samples supply preparation parameters, not a quantum history.
The kinetic reduction is the actual 68-coordinate Schur reduction. Every
56-dimensional seed direction is retained before exact circle projection.
Bessel evaluations and metric matrices are numerical diagnostics; the paper
proves exact initial cotangent/projection bounds and the operator domain.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import i0e, i1e

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'runtime/whitney_quantum_packet_receipt.json'
SCHEMA = 'oph.whitney_quantum_packet.v1'
SCOPE = 'FULL_56D_NEUTRAL_PACKET_PREPARATION__NO_QUANTUM_PROPAGATION_CERTIFICATE'
PARENT_PATH = 'code/electromagnetism/runtime/whitney_charged_dynamics_receipt.json'
PIN_PATHS = (
    'Lean/Screen/SeamCurrentEdge30Moment.lean',
    'Lean/Screen/SeamCurrentCarrierQuotient.lean',
    'Lean/ObserverPatchHolography/CoreAxioms.lean',
    'code/electromagnetism/verify_cone_whitney_bridge.py',
    'code/electromagnetism/whitney_interacting_quantum.py',
    'code/electromagnetism/verify_whitney_quantum_state.py',
    'code/electromagnetism/verify_whitney_charged_dynamics.py',
    'code/electromagnetism/whitney_quantum_packet.py',
    'code/electromagnetism/verify_whitney_quantum_packet.py',
    'code/electromagnetism/test_whitney_quantum_packet.py',
    'paper/tex_fragments/WHITNEY_INTERACTING_QUANTUM.tex',
    'paper/tex_fragments/WHITNEY_QUANTUM_PACKET.tex',
    PARENT_PATH,
)


def sibling(name, canonical_name=None):
    """Load a pinned sibling independently of the caller's working directory."""
    path = HERE/(name+'.py')
    module_name = canonical_name or 'packet_preparation_'+name
    existing = sys.modules.get(module_name)
    if existing is not None:
        if Path(existing.__file__).resolve() != path.resolve():
            raise ImportError('unexpected sibling module origin')
        return existing
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError('missing sibling '+name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


sibling('verify_cone_whitney_bridge', 'verify_cone_whitney_bridge')
quantum = sibling('whitney_interacting_quantum')
parent = sibling('verify_whitney_charged_dynamics')


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode('ascii')


def real_vector(value, size, name):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape != (size,):
        raise ValueError('finite real '+name+' required')
    result = np.asarray(raw, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError('finite real '+name+' required')
    return result


def positive(value, name):
    if isinstance(value, (bool, np.bool_, str, bytes)) or np.ndim(value) or not np.isrealobj(value):
        raise ValueError('positive finite '+name+' required')
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError('positive finite '+name+' required') from error
    if not np.isfinite(result) or result <= 0:
        raise ValueError('positive finite '+name+' required')
    return result


def phase_space(q68, velocity68, mesh=None, charge=.25):
    """True tangent rechart, then the full Schur cotangent at c=0.

    No Gauss condition is silently imposed. The difference between the
    transformed temporal potential and its Schur minimizer is returned.
    """
    mesh = quantum.geometry(4) if mesh is None else mesh
    q68 = real_vector(q68, 68, 'configuration')
    velocity68 = real_vector(velocity68, 68, 'velocity')
    charge = positive(charge, 'charge')
    a, psi = q68[:42], q68[42:55]+1j*q68[55:]
    av, pv = velocity68[:42], velocity68[42:55]+1j*velocity68[55:]
    b = mesh.mean_zero
    laplacian = b.T@mesh.d.T@mesh.mass@mesh.d@b
    xi = -b@np.linalg.solve(laplacian, b.T@mesh.d.T@mesh.mass@a)
    xidot = -b@np.linalg.solve(laplacian, b.T@mesh.d.T@mesh.mass@av)
    rotation = np.exp(1j*charge*xi)
    ac, pc = a+mesh.d@xi, rotation*psi
    acv, pcv = av+mesh.d@xidot, rotation*(pv+1j*charge*xidot*psi)
    full_q, full_v = np.r_[ac, pc.real, pc.imag], np.r_[acv, pcv.real, pcv.imag]
    q, velocity = mesh.slice.T@full_q, mesh.slice.T@full_v
    gamma, potential, eta_map, _ = quantum.reduced_coefficients(ac, pc, charge, .5, .25, mesh)
    momentum = gamma@velocity
    generator = np.r_[np.zeros(30), -pc.imag, pc.real]
    eta = mesh.mean_zero@eta_map@velocity
    sign, logdet = np.linalg.slogdet(gamma)
    if sign != 1:
        raise ValueError('nonpositive reduced metric')
    return {'q': q, 'velocity': velocity, 'momentum': momentum,
            'full_q': full_q, 'full_velocity': full_v, 'gauge_parameter': xi,
            'gauge_parameter_velocity': xidot, 'transformed_scalar_potential': -xidot,
            'schur_scalar_potential': eta, 'constant_moment_map': float(momentum@generator),
            'minimizer_defect': float(np.max(abs(eta+xidot))),
            'coulomb_defect': float(np.max(abs(mesh.d.T@mesh.mass@ac))),
            'cotangent_identity_defect': float(np.max(abs(momentum-gamma@velocity))),
            'log_rho': float(logdet/2), 'potential': potential,
            'gamma_min_eigenvalue': float(np.linalg.eigvalsh(gamma)[0])}


def overlap_parameters(center, momentum, sigma, hbar=1):
    center = real_vector(center, 56, 'packet center')
    momentum = real_vector(momentum, 56, 'packet momentum')
    sigma, hbar = positive(sigma, 'width'), positive(hbar, 'hbar')
    x, p = center[30:], momentum[30:]
    jx = np.r_[-x[13:], x[:13]]
    a = float(x@x/(4*sigma**2)+sigma**2*(p@p)/hbar**2)
    b = float(p@jx/hbar)
    discriminant = a*a-b*b
    if not np.isfinite([a, b, discriminant]).all() or discriminant < -1e-11*(1+a*a):
        raise ValueError('invalid overlap parameters')
    z = float(np.sqrt(max(0, discriminant)))
    norm_squared = float(np.exp(z-a)*i0e(z))
    return {'A': a, 'B': b, 'norm_squared': norm_squared}


def seed_log_half_density(point, center, momentum, sigma, hbar=1):
    point = real_vector(point, 56, 'packet evaluation point')
    center = real_vector(center, 56, 'packet center')
    momentum = real_vector(momentum, 56, 'packet momentum')
    sigma, hbar = positive(sigma, 'width'), positive(hbar, 'hbar')
    delta = point-center
    return complex(-14*np.log(2*np.pi*sigma**2)-delta@delta/(4*sigma**2), momentum@delta/hbar)


def rotate(vector, angle):
    vector = real_vector(vector, 56, 'circle vector')
    if not np.ndim(angle) == 0 or not np.isfinite(angle):
        raise ValueError('finite angle required')
    z = (vector[30:43]+1j*vector[43:])*np.exp(1j*angle)
    return np.r_[vector[:30], z.real, z.imag]


def projected_half_density(point, center, momentum, sigma, hbar=1, nodes=256):
    """Numerical circle quadrature; no certified pointwise quadrature error."""
    if type(nodes) is not int or nodes < 16:
        raise ValueError('at least sixteen integer circle nodes required')
    norm = overlap_parameters(center, momentum, sigma, hbar)['norm_squared']**.5
    if norm == 0:
        raise ValueError('projection norm underflows at these parameters')
    values = [np.exp(seed_log_half_density(point, rotate(center, a), rotate(momentum, a), sigma, hbar))
              for a in np.arange(nodes)*(2*np.pi/nodes)]
    return complex(sum(values)/nodes/norm)


def scalar_radius_moment(center, momentum, sigma, hbar=1):
    """Exact nodal-radius formula, evaluated numerically; not a spatial L2 norm."""
    row = overlap_parameters(center, momentum, sigma, hbar)
    sigma, hbar = positive(sigma, 'width'), positive(hbar, 'hbar')
    x, p = np.asarray(center)[30:], np.asarray(momentum)[30:]
    z = np.sqrt(max(0, row['A']**2-row['B']**2))
    ratio = float(i1e(z)/i0e(z))
    return float(26*sigma**2+(x@x)/2-2*sigma**4*(p@p)/hbar**2+2*sigma**2*z*ratio)


def exact_initial():
    momentum_norm = [Q(91, 6), Q(13, 2)]
    widths = []
    for sigma in (Q(1, 4), Q(1, 2), Q(1)):
        aa = [Q(13)/(4*sigma**2)+sigma**2*momentum_norm[0], sigma**2*momentum_norm[1]]
        # Vol<18, and 39 Vol^2/400 bounds the squared scalar momentum.
        upper = Q(13)/(4*sigma**2)+sigma**2*Q(39, 400)*18**2
        widths.append({'sigma': str(sigma), 'A_in_Qsqrt5': [str(x) for x in aa],
            'A_rational_upper': str(upper), 'B': '0', 'projection_norm_squared_lower': '1/64',
            'seed_position_variance': str(sigma**2), 'seed_momentum_variance_hbar1': str(1/(4*sigma**2))})
    return {'volume_in_Qsqrt5': ['10', '10/3'],
        'radial_velocity_in_Qsqrt5': ['7/80', '3/80'],
        'center_gauge_factor': '12/13', 'boundary_gauge_factor': '-1/13',
        'scalar_center_real': '1', 'scalar_boundary_real': '1',
        'momentum_center_imaginary_in_Qsqrt5': ['3', '1'],
        'momentum_boundary_imaginary_in_Qsqrt5': ['-1/4', '-1/12'],
        'scalar_momentum_norm_squared_in_Qsqrt5': [str(x) for x in momentum_norm],
        'radiative_center_and_momentum': '0', 'constant_moment_map': '0',
        'proof_labels': ['prop:whitney-packet-coulomb', 'prop:whitney-packet-neutral-projection'],
        'widths': widths}


def contract():
    return {'configuration_dimension': 56, 'radiative_dimension': 30, 'scalar_real_dimension': 26,
        'seed_covariance': 'sigma^2 I_56; full rank; no classical fixed-space restriction',
        'projected_covariance': 'not the seed covariance; circle averaging changes scalar moments',
        'seed_formula': '(2*pi*sigma^2)^(-14)*exp(-|x-q|^2/(4*sigma^2)+i*p.(x-q)/hbar)',
        'state_formula': 'F=rho^(-1/2)*P0(seed)/||P0(seed)||; rho=sqrt(det(gamma))',
        'quantization': 'same full interacting Laplace-Beltrami Hamiltonian and metric volume',
        'circle_projection': 'exact group integral; floating Bessel values are diagnostics',
        'physical_comparison': False, 'empirical_prediction': False, 'observer_preparation': False,
        'quantum_propagation_certified': False, 'covariance_evolution_computed': False,
        'full_56D_preparation': True, 'neutral_strong_domain_state': True,
        'classical_samples_are_quantum_means': False,
        'scalar_one_point_mean': '0 by residual U(1) invariance',
        'sample_role': 'two separate packet preparations from adjacent classical samples; not quantum propagation',
        'numeric_scope': 'float64 classical data, phase-space maps and exact-degree simplex integration; no interval trajectory enclosure',
        'requested_propagation_target': {'T': '1/40', 'hbar': '1', 'norm_error': '1/10',
            'status': 'NOT_CERTIFIED', 'residual_upper_bound': None, 'all_tail_bound': None}}


def build():
    raw = parent.load(ROOT/PARENT_PATH)
    parent.verify(raw)
    mesh = quantum.geometry(4)
    samples = []
    for index in (0, 1):
        source = raw['samples'][index]
        data = phase_space(source['q'], source['velocity'], mesh)
        row = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in data.items()}
        row.update({'parent_sample_index': index, 'model_time': source['t'], 'widths': []})
        for sigma in (.25, .5, 1.):
            values = overlap_parameters(data['q'], data['momentum'], sigma)
            values.update({'sigma': str(Q(sigma)), 'scalar_radius_numeric': scalar_radius_moment(data['q'], data['momentum'], sigma)})
            row['widths'].append(values)
        samples.append(row)
    return {'schema': SCHEMA, 'scope': SCOPE,
        'source_pins': {x: hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in PIN_PATHS},
        'parameters': {'e': '1/4', 'm_squared': '1/2', 'g': '1/4', 'hbar': '1'},
        'contract': contract(), 'initial_exact': exact_initial(),
        'orthonormal_coulomb_frame': mesh.slice[:42, :30].tolist(), 'samples': samples}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({'receipt': str(args.output), 'samples': 2, 'quantum_propagation_certified': False}))
