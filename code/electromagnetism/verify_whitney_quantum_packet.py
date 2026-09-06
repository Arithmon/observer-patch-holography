"""Independent Coulomb tangent, Schur cotangent and circle-projection replay.

Neither the packet producer nor interacting coefficient evaluator is imported.
The metric at the charged Coulomb representative a=0 is assembled from
simplex monomials; the gauge solve uses a grounded augmented Laplacian.
The parent trajectory is independently verified fresh after packet checks.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import combinations, product
import json
from math import factorial, prod
from pathlib import Path

import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'runtime/whitney_quantum_packet_receipt.json'
PARENT_PATH = 'code/electromagnetism/runtime/whitney_charged_dynamics_receipt.json'
PIN_PATHS = {
    'Lean/Screen/SeamCurrentEdge30Moment.lean', 'Lean/Screen/SeamCurrentCarrierQuotient.lean',
    'Lean/ObserverPatchHolography/CoreAxioms.lean',
    'code/electromagnetism/verify_cone_whitney_bridge.py',
    'code/electromagnetism/whitney_interacting_quantum.py',
    'code/electromagnetism/verify_whitney_quantum_state.py',
    'code/electromagnetism/verify_whitney_charged_dynamics.py',
    'code/electromagnetism/whitney_quantum_packet.py',
    'code/electromagnetism/verify_whitney_quantum_packet.py',
    'code/electromagnetism/test_whitney_quantum_packet.py',
    'paper/tex_fragments/WHITNEY_INTERACTING_QUANTUM.tex',
    'paper/tex_fragments/WHITNEY_QUANTUM_PACKET.tex', PARENT_PATH,
}


def sibling(name):
    spec = importlib.util.spec_from_file_location('packet_independent_'+name, HERE/(name+'.py'))
    if spec is None or spec.loader is None:
        raise ImportError('missing independent source verifier')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


source = sibling('verify_whitney_quantum_state')
parent = sibling('verify_whitney_charged_dynamics')
require = source.require


def load(path=OUTPUT):
    return source.load(path)


def keys(value, expected, name):
    require(type(value) is dict and set(value) == set(expected), name+' key census')


def equal(actual, expected, name):
    try:
        require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
                json.dumps(expected, sort_keys=True, allow_nan=False), name)
    except (TypeError, ValueError) as error:
        raise ValueError(name) from error


def numbers(value, shape, name):
    def valid(item):
        return all(valid(v) for v in item) if type(item) is list else type(item) in (int, float)
    require(valid(value), name+' numeric types')
    result = np.asarray(value, dtype=float)
    require(result.shape == shape and np.isfinite(result).all(), name+' shape/finiteness')
    return result


def close(value, expected, name, tolerance=3e-10):
    actual = numbers(value, np.shape(expected), name)
    require(np.max(abs(actual-expected), initial=0) <= tolerance*(1+np.max(abs(expected), initial=0)), name)
    return actual


def simplex(powers):
    return Q(6*prod(factorial(x) for x in powers), factorial(3+sum(powers)))


def metric_and_potential(psi, xyz, edges, tetrahedra):
    """a=0 monomial assembly independent of all element-quadrature code."""
    g = np.zeros((68, 68)); m = np.zeros((42, 42)); potential = 0.
    for tet in tetrahedra:
        points = xyz[list(tet)]
        gradients = np.linalg.inv(np.column_stack((np.ones(4), points)))[1:].T
        volume = abs(np.linalg.det(points[1:]-points[0]))/6
        columns, local_edges = [], []
        for edge, (u, v) in enumerate(edges):
            if u in tet and v in tet:
                i, j = tet.index(u), tet.index(v)
                local_edges.append((edge, i, j))
                columns.append((edge, .25j*(psi[u]-psi[v]), (i, j)))
        for i, node in enumerate(tet):
            columns.extend(((42+node, 1, (i,)), (55+node, 1j, (i,))))
        for e, i, j in local_edges:
            for f, k, l in local_edges:
                m[e, f] += volume/20*((1+int(i == k))*(gradients[j]@gradients[l])
                    -(1+int(i == l))*(gradients[j]@gradients[k])
                    -(1+int(j == k))*(gradients[i]@gradients[l])
                    +(1+int(j == l))*(gradients[i]@gradients[k]))
        for col, value, powers in columns:
            for other, val, powers2 in columns:
                degree = powers+powers2
                integral = volume*float(simplex(tuple(degree.count(i) for i in range(4))))
                g[col, other] += 2*np.real(np.conj(value)*val)*integral
        local = psi[list(tet)]
        gradient = local@gradients
        potential += volume*float(np.vdot(gradient, gradient).real)
        for i, j in product(range(4), repeat=2):
            potential += .5*volume*(local[i].conjugate()*local[j]).real*float(simplex(tuple((i,j).count(k) for k in range(4))))
        for i, j, k, l in product(range(4), repeat=4):
            power = (i,j,k,l)
            potential += .125*volume*(local[i].conjugate()*local[j]*local[k].conjugate()*local[l]).real*float(simplex(tuple(power.count(n) for n in range(4))))
    g[:42, :42] += m
    return g, m, potential


def replay_phase_space(q68, v68, xyz, edges, tetrahedra, frame):
    q68, v68 = np.asarray(q68), np.asarray(v68)
    psi = q68[42:55]+1j*q68[55:]
    d = np.zeros((42, 13))
    for row, (i, j) in enumerate(edges):
        d[row, i], d[row, j] = -1, 1
    _, mass, _ = metric_and_potential(np.ones(13, complex), xyz, edges, tetrahedra)
    laplacian = d.T@mass@d
    augmented = np.block([[laplacian, np.ones((13, 1))], [np.ones((1, 13)), np.zeros((1,1))]])
    xi = np.linalg.solve(augmented, np.r_[-d.T@mass@q68[:42], 0])[:13]
    xidot = np.linalg.solve(augmented, np.r_[-d.T@mass@v68[:42], 0])[:13]
    phase = np.exp(.25j*xi)
    scalar = phase*psi
    scalar_velocity = phase*(v68[42:55]+1j*v68[55:]+.25j*xidot*psi)
    ac, av = q68[:42]+d@xi, v68[:42]+d@xidot
    require(np.max(abs(ac)) < 1e-11, 'parent radial field has zero Coulomb representative')
    full_q, full_v = np.r_[ac,scalar.real,scalar.imag], np.r_[av,scalar_velocity.real,scalar_velocity.imag]
    section = np.zeros((68,56));section[:42,:30]=frame;section[42:,30:]=np.eye(26)
    close((frame.T@frame).tolist(), np.eye(30), 'orthonormal Coulomb frame', 2e-12)
    close((d.T@mass@frame).tolist(), np.zeros((13,30)), 'mass Coulomb frame', 2e-12)
    g, _, potential = metric_and_potential(scalar, xyz, edges, tetrahedra)
    # A different, nonorthonormal gauge basis verifies basis independence.
    basis = np.vstack((-np.ones((1,12)), np.eye(12)))
    vertical = np.vstack((d, -.25*np.diag(scalar.imag), .25*np.diag(scalar.real)))@basis
    inertia, coupling = vertical.T@g@vertical, vertical.T@g@section
    gamma = section.T@g@section-coupling.T@np.linalg.solve(inertia,coupling)
    q, velocity = section.T@full_q, section.T@full_v
    momentum = gamma@velocity
    eta = -basis@np.linalg.solve(inertia,coupling@velocity)
    _, logdet = np.linalg.slogdet(gamma)
    return {'q':q, 'velocity':velocity, 'momentum':momentum, 'full_q':full_q,
        'full_velocity':full_v, 'gauge_parameter':xi,'gauge_parameter_velocity':xidot,
        'transformed_scalar_potential':-xidot,'schur_scalar_potential':eta,
        'constant_moment_map':float(momentum@np.r_[np.zeros(30),-scalar.imag,scalar.real]),
        'minimizer_defect':float(np.max(abs(eta+xidot))),
        'coulomb_defect':float(np.max(abs(d.T@mass@ac))),
        'cotangent_identity_defect':float(np.max(abs(momentum-gamma@velocity))),
        'log_rho':float(logdet/2),'potential':float(potential),
        'gamma_min_eigenvalue':float(np.linalg.eigvalsh(gamma)[0])}


def circle_observables(q, p, sigma):
    """Direct one-dimensional overlap integration, not a Bessel implementation."""
    x, p = q[30:],p[30:]
    xx, pp = float(x@x),float(p@p)
    b = float(p@np.r_[-x[13:],x[:13]])
    a = xx/(4*sigma**2)+sigma**2*pp
    def overlap(theta):
        return np.exp(-a*(1-np.cos(theta))-1j*b*np.sin(theta))
    norm = quad(lambda theta: overlap(theta).real, -np.pi,np.pi,epsabs=2e-13,epsrel=2e-13)[0]/(2*np.pi)
    def radius(theta):
        multiplier = 26*sigma**2+xx/2*(1+np.cos(theta))-2*sigma**4*pp*(1-np.cos(theta))-2j*sigma**2*b*np.sin(theta)
        return (overlap(theta)*multiplier).real
    moment = quad(radius,-np.pi,np.pi,epsabs=2e-12,epsrel=2e-13)[0]/(2*np.pi*norm)
    return {'A':a,'B':b,'norm_squared':norm,'scalar_radius_numeric':moment,'sigma':str(Q(sigma))}


def initial_algebra(volume):
    # Compute the two scalar mass rows from simplex moments and node incidence.
    c = 2*(3*simplex((2,0,0,0))-3*simplex((1,1,0,0)))
    boundary = 2*Q(1,4)*(3*simplex((1,1,0,0))-(Q(1,4)-simplex((1,1,0,0))))
    require(c == Q(3,10) and boundary == -Q(1,40), 'initial simplex cotangent')
    center = [c*x for x in volume]; edge = [boundary*x for x in volume]
    momentum_norm = [(c*c+12*boundary**2)*(volume[0]**2+5*volume[1]**2),
                     (c*c+12*boundary**2)*2*volume[0]*volume[1]]
    widths=[]
    for sigma in (Q(1,4),Q(1,2),Q(1)):
        upper=Q(13)/(4*sigma**2)+(c*c+12*boundary**2)*sigma**2*18**2
        require(upper < 64,'projection bound range')
        widths.append({'sigma':str(sigma),'A_in_Qsqrt5':[str(Q(13)/(4*sigma**2)+sigma**2*momentum_norm[0]),str(sigma**2*momentum_norm[1])],
            'A_rational_upper':str(upper),'B':'0','projection_norm_squared_lower':'1/64',
            'seed_position_variance':str(sigma**2),'seed_momentum_variance_hbar1':str(Q(1)/(4*sigma**2))})
    return {'volume_in_Qsqrt5':[str(x) for x in volume], 'radial_velocity_in_Qsqrt5':['7/80','3/80'],
        'center_gauge_factor':'12/13','boundary_gauge_factor':'-1/13',
        'scalar_center_real':'1','scalar_boundary_real':'1',
        'momentum_center_imaginary_in_Qsqrt5':[str(x) for x in center],
        'momentum_boundary_imaginary_in_Qsqrt5':[str(x) for x in edge],
        'scalar_momentum_norm_squared_in_Qsqrt5':[str(x) for x in momentum_norm],
        'radiative_center_and_momentum':'0','constant_moment_map':'0',
        'proof_labels':['prop:whitney-packet-coulomb','prop:whitney-packet-neutral-projection'], 'widths':widths}


def expected_contract():
    return {'configuration_dimension':56,'radiative_dimension':30,'scalar_real_dimension':26,
        'seed_covariance':'sigma^2 I_56; full rank; no classical fixed-space restriction',
        'projected_covariance':'not the seed covariance; circle averaging changes scalar moments',
        'seed_formula':'(2*pi*sigma^2)^(-14)*exp(-|x-q|^2/(4*sigma^2)+i*p.(x-q)/hbar)',
        'state_formula':'F=rho^(-1/2)*P0(seed)/||P0(seed)||; rho=sqrt(det(gamma))',
        'quantization':'same full interacting Laplace-Beltrami Hamiltonian and metric volume',
        'circle_projection':'exact group integral; floating Bessel values are diagnostics',
        'physical_comparison':False,'empirical_prediction':False,'observer_preparation':False,
        'quantum_propagation_certified':False,'covariance_evolution_computed':False,
        'full_56D_preparation':True,'neutral_strong_domain_state':True,
        'classical_samples_are_quantum_means':False,'scalar_one_point_mean':'0 by residual U(1) invariance',
        'sample_role':'two separate packet preparations from adjacent classical samples; not quantum propagation',
        'numeric_scope':'float64 classical data, phase-space maps and exact-degree simplex integration; no interval trajectory enclosure',
        'requested_propagation_target':{'T':'1/40','hbar':'1','norm_error':'1/10','status':'NOT_CERTIFIED','residual_upper_bound':None,'all_tail_bound':None}}


def verify(packet):
    keys(packet, {'schema','scope','source_pins','parameters','contract','initial_exact','orthonormal_coulomb_frame','samples'},'packet')
    equal(packet['schema'],'oph.whitney_quantum_packet.v1','schema')
    equal(packet['scope'],'FULL_56D_NEUTRAL_PACKET_PREPARATION__NO_QUANTUM_PROPAGATION_CERTIFICATE','scope')
    keys(packet['source_pins'],PIN_PATHS,'source pins')
    for path in PIN_PATHS:
        equal(packet['source_pins'][path],hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),'source pin '+path)
    equal(packet['parameters'],{'e':'1/4','m_squared':'1/2','g':'1/4','hbar':'1'},'parameters')
    equal(packet['contract'],expected_contract(),'preparation and non-propagation contract')
    volume,xyz,edges,tets=source.exact_geometry()
    equal(packet['initial_exact'],initial_algebra(volume),'exact initial algebra and projection bound')
    proof=(ROOT/'paper/tex_fragments/WHITNEY_QUANTUM_PACKET.tex').read_text(encoding='utf-8')
    for label in initial_algebra(volume)['proof_labels']:
        require('\\label{'+label+'}' in proof,'proof reference')
    frame=numbers(packet['orthonormal_coulomb_frame'],(42,30),'Coulomb frame')
    require(type(packet['samples']) is list and len(packet['samples'])==2,'sample census')
    raw=parent.load(ROOT/PARENT_PATH)
    diagnostics=[]
    for index,row in enumerate(packet['samples']):
        require(type(row) is dict,'sample object')
        original=raw['samples'][index]
        expected=replay_phase_space(original['q'],original['velocity'],xyz,edges,tets,frame)
        keys(row,set(expected)|{'parent_sample_index','model_time','widths'},'sample')
        equal(row['parent_sample_index'],index,'parent index')
        equal(row['model_time'],index/40,'model sample time')
        for key,value in expected.items():
            close(row[key],value,'independent phase-space '+key)
        require(abs(expected['constant_moment_map'])<1e-9 and expected['minimizer_defect']<1e-9,'parent numerical Gauss compatibility')
        require(type(row['widths']) is list and len(row['widths'])==3,'width census')
        for width,sigma in zip(row['widths'],(.25,.5,1),strict=True):
            expected_width=circle_observables(expected['q'],expected['momentum'],sigma)
            keys(width,expected_width,'width')
            equal(width['sigma'],expected_width['sigma'],'width parameter')
            for key in ('A','B','norm_squared','scalar_radius_numeric'):
                close(width[key],expected_width[key],'independent circle '+key,2e-11)
        diagnostics.append({'model_time':row['model_time'],'gauss_abs':abs(expected['constant_moment_map']),
            'minimizer_defect':expected['minimizer_defect'],
            'projection_norm_squared_numeric':[x['norm_squared'] for x in row['widths']]})
    # Fresh full parent replay is deliberately last: malformed packets fail
    # before the longer numerical trajectory validation, never instead of it.
    parent_result=parent.verify(raw)
    require(parent_result['accepted'] is True,'independent parent trajectory')
    return {'accepted':True,'configuration_dimension':56,'radiative_dimension':30,
        'prepared_classical_samples':2,'widths':['1/4','1/2','1'],
        'initial_projection_norm_squared_lower':'1/64','full_56D_preparation':True,
        'neutral_strong_domain_state':True,'quantum_propagation_certified':False,
        'covariance_evolution_computed':False,'physical_comparison':False,
        'analytic_domain_proved_by_numeric_replay':False,'diagnostics':diagnostics}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,default=OUTPUT)
    args=parser.parse_args()
    print(json.dumps(verify(load(args.receipt)),sort_keys=True))
