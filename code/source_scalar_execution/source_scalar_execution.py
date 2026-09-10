"""Complete q=5 source-addressed scalar execution in exact Q(phi).

All substeps and authenticated reads are retained. The clock is a declared
model parameter; no old count-clock theorem is applied to this new ancestry.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import factorial
from pathlib import Path
import json

from scalar_execution_algebra import Q, rational_bounds

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'source_scalar_execution_receipt.json'
TAU = Q(F(-1, 35), F(2, 35))
STEPS = 21
PINS = (
    'code/source_scalar_execution/scalar_execution_algebra.py',
    'code/source_scalar_execution/source_scalar_execution.py',
    'code/source_scalar_execution/verify_source_scalar_execution.py',
    'code/source_scalar_execution/test_source_scalar_execution.py',
    'paper/tex_fragments/SOURCE_SCALAR_EXECUTION.tex',
    'paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex',
    'code/source_scalar_packet/source_common_scalar_receipt.json',
    'Lean/Screen/PrimitivePortFrameQuotient.lean',
    'Lean/Screen/SeamCurrentCarrierQuotient.lean',
    'Lean/Screen/PortFrameGram.lean',
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'),
                       allow_nan=False, ensure_ascii=True)+'\n').encode('ascii')


def digest(value):
    return sha256(canonical(value)).hexdigest()


def model():
    # Actual sorted orbit for b=0..4, followed by the source-coded boundary1.
    points = [Q(), Q(-3, 2), Q(-6, 4), Q(-1, 1), Q(-4, 3), Q(1)]
    gaps = [b-a for a, b in zip(points, points[1:])]
    w = [(a+b)/2 for a, b in zip(gaps, gaps[1:])]
    sites = list(product(range(4), repeat=3))
    def index(p):
        return (p[0]*4+p[1])*4+p[2]
    masses, rows, preparation, detector, addresses = [], [], [], [], []
    for site in sites:
        xyz = [points[i+1] for i in site]
        masses.append(w[site[0]]*w[site[1]]*w[site[2]])
        row = {index(site): Q(1)}
        for axis, i in enumerate(site):
            for direction, gap in ((-1, gaps[i]), (1, gaps[i+1])):
                coefficient = 1/(w[i]*gap)
                row[index(site)] += coefficient
                adjacent = list(site)
                adjacent[axis] += direction
                if 0 <= adjacent[axis] < 4:
                    row[index(adjacent)] = -coefficient
        rows.append(sorted(row.items()))
        polynomial = Q(64)
        for value in xyz:
            polynomial *= value*(1-value)
        left = (xyz[0]-F(9, 20)).sign() < 0
        right = (xyz[0]-F(11, 20)).sign() > 0
        preparation.append(4*polynomial if left else Q())
        detector.append(polynomial if right else Q())
        aa = [int(v.a) for v in xyz]
        bb = [int(v.b) for v in xyz]
        z = [bb[1]-aa[0], bb[1]+aa[0], bb[2]-aa[1], bb[2]+aa[1], bb[0]-aa[2], bb[0]+aa[2]]
        h = sum(z)//2
        currents = [-z[5], h-z[0], h-z[2]-z[3], h-z[1]-z[4]-z[5], h-z[2]-z[3]-z[4], z[3]]
        addresses.append({'site': list(site), 'coordinate_Qphi': [v.encode() for v in xyz],
                          'D6_control': z, 'signed_seam_currents': currents})
    return points, gaps, masses, rows, preparation, detector, addresses


def action(rows, values):
    return [sum((coefficient*values[j] for j, coefficient in row), Q()) for row in rows]


def inner(mass, u, v):
    return sum((m*a*b for m, a, b in zip(mass, u, v)), Q())


def energy(mass, rows, previous, current):
    force = action(rows, current)
    centered = [(q-p)/TAU-TAU*f/2 for p, q, f in zip(previous, current, force)]
    return inner(mass, centered, centered)/2+inner(mass, current, force)/2-TAU**2*inner(mass, force, force)/8


def execute(rows, velocity, schedule):
    if (type(schedule) is not list or any(type(i) is not int for i in schedule)
            or sorted(schedule) != list(range(64))):
        raise ValueError('complete exact site permutation required')
    if len(velocity) != 64 or any(type(v) is not Q for v in velocity):
        raise ValueError('64 exact algebraic velocity inputs required')
    previous, current = [-TAU*v for v in velocity], [Q() for _ in velocity]
    records, layers = [], [[v.encode() for v in previous], [v.encode() for v in current]]
    committed = {}
    chain = '0'*64
    def write(phase, site, reads, resource, version, value):
        nonlocal chain
        identity = [phase, site]
        material = {'id': identity, 'reads': reads,
                    'write': [resource, version, identity, value.encode()], 'ledger_parent': chain}
        chain = digest(material)
        material['hash'] = chain
        records.append(material)
        committed[(resource, version)] = (identity, value)
    for phase, values, parity in ((0, previous, 1), (1, current, 0)):
        for site in schedule:
            write(phase, site, [], parity*64+site, 1, values[site])
    for j in range(1, STEPS+1):
        following = [None]*64
        for site in schedule:
            reads = []
            old_resource = (j % 2)*64+site
            old_version = (j+1)//2
            writer, old = committed[(old_resource, old_version)]
            reads.append([old_resource, old_version, writer, old.encode()])
            for other, coefficient in rows[site]:
                resource = ((j-1) % 2)*64+other
                version = j//2+1
                writer, value = committed[(resource, version)]
                reads.append([resource, version, writer, value.encode()])
            force = sum((coefficient*current[other] for other, coefficient in rows[site]), Q())
            value = 2*current[site]-previous[site]-TAU**2*force
            write(j+1, site, reads, old_resource, old_version+1, value)
            following[site] = value
        previous, current = current, following
        layers.append([v.encode() for v in current])
    return {'schedule': schedule, 'events': records, 'layers': layers, 'final_hash': chain,
            'seed_writes': 128, 'update_writes': 1344, 'dynamic_field_reads': 8736}


def sqrt_bounds(value, denominator=10**12):
    if value.sign() < 0:
        raise ValueError('negative squared norm')
    lo, hi = 0, denominator
    while (Q(F(hi, denominator)**2)-value).sign() < 0:
        hi *= 2
    while hi-lo > 1:
        mid = (hi+lo)//2
        if (Q(F(mid, denominator)**2)-value).sign() <= 0:
            lo = mid
        else:
            hi = mid
    return F(lo, denominator), F(hi, denominator)


def widen(bounds, radius):
    scale = 10**12
    lo, hi = F(bounds[0])-radius, F(bounds[1])+radius
    return [str(F((lo*scale).numerator//(lo*scale).denominator, scale)),
            str(F(-((-hi*scale).numerator//(-hi*scale).denominator), scale))]


def diagnostics(mass, rows, velocity, detector, trace):
    from scalar_execution_algebra import parse
    layers = [[parse(v) for v in row] for row in trace['layers']]
    energy0 = inner(mass, velocity, velocity)/2
    energies = [energy(mass, rows, layers[j], layers[j+1]) for j in range(STEPS+1)]
    if any(v != energy0 for v in energies):
        raise ValueError('exact modified energy changed')
    # Whole64-mode continuous solution. The sine Taylor polynomial through
    # degree161 has an explicit uniform remainder on0<=t<=3/sqrt5<3/2.
    powers = [velocity]
    for _ in range(80):
        powers.append(action(rows, powers[-1]))
    tail = 4*F(3, 2)**163*614**81/factorial(163)
    if tail >= F(1, 10**30):
        raise ValueError('Taylor remainder budget exceeded')
    comparison = []
    for j in range(1, STEPS+1):
        polynomial = [Q() for _ in velocity]
        for k, vector in enumerate(powers):
            coefficient = TAU*F((-1)**k*j**(2*k+1), 245**k*factorial(2*k+1))
            for i, v in enumerate(vector):
                polynomial[i] += coefficient*v
        actual = layers[j+1]
        difference = [u-v for u, v in zip(actual, polynomial)]
        err = sqrt_bounds(inner(mass, difference, difference))
        step_detector = inner(mass, detector, actual)
        continuous_detector = inner(mass, detector, polynomial)
        comparison.append({'step': j, 'model_time_Qphi': (j*TAU).encode(),
                           'executed_detector_Qphi': step_detector.encode(),
                           'canonical_commutator_over_i_hbar_Qphi': (-step_detector/4).encode(),
                           'continuous_detector_interval': widen(rational_bounds(continuous_detector), tail),
                           'full_mass_norm_discretization_error': widen(err, tail)})
    ancestors = {}
    for event in trace['events']:
        parent_ids = [tuple(read[2]) for read in event['reads']]
        past = set(parent_ids)
        for parent in parent_ids:
            past.update(ancestors[parent])
        ancestors[tuple(event['id'])] = past
    lower = (1, 21)
    counts = []
    for step in (7, 14, 21):
        upper = (step+1, 21)
        interval = ancestors[upper] | {upper}
        count = sum(e == lower or lower in ancestors[e] for e in interval)
        counts.append({'lower': list(lower), 'upper': list(upper), 'inclusive_count': count})
    # Source-plane intervention can reach a detector at x=.854 in one old
    # count layer. This finite-step example concerns actual ancestry, not a
    # continuum or finite-detector no-go.
    target = 53
    outside_value = layers[8][target]
    if outside_value.sign() == 0:
        raise ValueError('the declared finite ancestry example is vacuous')
    return {'modified_energy_Qphi': energy0.encode(), 'modified_energy_interval': rational_bounds(energy0),
            'energy_equal_at_all_22_states': True,
            'taylor_remainder_upper': str(tail), 'comparisons': comparison,
            'actual_field_interval_counts': counts,
            'first_coarse_layer_outside_radius': {'target_site': target, 'step': 7,
                'source_plane_x_Qphi': Q(-3, 2).encode(), 'target_x_Qphi': Q(-4, 3).encode(),
                'old_radius_squared': '1/5', 'x_displacement_squared_Qphi': (Q(-1, 1)**2).encode(),
                'field_response_Qphi': outside_value.encode(), 'response_interval': rational_bounds(outside_value),
                'all_nonzero_preparation_sites_outside_old_radius': True}}


def build():
    points, gaps, masses, rows, velocity, detector, addresses = model()
    if TAU**2 != Q(F(1, 245)) or any((5*h*h-1).sign() >= 0 for h in gaps):
        raise ValueError('time or source edge certificate failed')
    if (Q(5, -3)-F(7, 50)).sign() <= 0:
        raise ValueError('minimum gap bound failed')
    if any((v-4).sign() > 0 or v.sign() < 0 for v in velocity):
        raise ValueError('prepared velocity bound failed')
    if any((v-1).sign() > 0 or v.sign() < 0 for v in detector):
        raise ValueError('detector bound failed')
    if (sum(masses, Q())-1).sign() > 0:
        raise ValueError('mass bound failed')
    traces = {}
    for name, schedule in [('ascending', list(range(64))), ('descending', list(reversed(range(64))))]:
        for preparation, values in [('baseline', [Q() for _ in velocity]), ('intervention', velocity)]:
            traces[name+'_'+preparation] = execute(rows, values, schedule)
    for preparation in ('baseline', 'intervention'):
        if traces['ascending_'+preparation]['layers'] != traces['descending_'+preparation]['layers']:
            raise ValueError('schedule changed field values')
    pins = {p: {'sha256': sha256((ROOT/p).read_bytes()).hexdigest(), 'bytes': (ROOT/p).stat().st_size}
            for p in PINS}
    return {'schema': 'oph.authenticated-golden-scalar-execution.v1',
        'model': {'q': 5, 'mutable_sites': 64, 'steps': STEPS, 'coarse_layers': 3, 'substeps_per_coarse_layer': 7,
                  'time_step_Qphi': TAU.encode(), 'time_step_squared': '1/245',
                  'source_length': 'L=2/sqrt(phi+2)', 'time_unit': 'c*t/L',
                  'dimensionless_mass': '1', 'boundary': 'Dirichlet cube',
                  'registers': 'two alternating immutable-version field registers per site',
                  'previous_seed': 'integrator memory -tau*v0; not a captured negative-time solution',
                  'preparation': 'v0=256*x*y*z*(1-x)*(1-y)*(1-z) on x<9/20; zero elsewhere',
                  'detector': '64*x*y*z*(1-x)*(1-y)*(1-z) on x>11/20; zero elsewhere',
                  'cfl_upper': '30049/12005', 'cfl_strict_margin': 'less than3, hence less than4'},
        'source_addresses': addresses,
        'mass_Qphi': [v.encode() for v in masses],
        'action_rows_Qphi': [[[i, v.encode()] for i, v in row] for row in rows],
        'initial_velocity_Qphi': [v.encode() for v in velocity], 'detector_Qphi': [v.encode() for v in detector],
        'traces': traces, 'diagnostics': diagnostics(masses, rows, velocity, detector, traces['ascending_intervention']),
        'scope': {'all_substeps_and_field_reads_stored': True, 'source_addresses_protected_preparation': True,
                  'all_arithmetic_exact_Qphi': True, 'action_boundary_model_time_supplied': True,
                  'source_action_family_preserved': True, 'same_q233_preparation_or_accuracy': False,
                  'read_values_writers_and_immutable_versions_replayed': True,
                  'serial_audit_parent_is_a_semantic_read_edge': False,
                  'old_count_clock_theorem_applied': False, 'native_repair_or_physical_clock_selected': False,
                  'quantum_scope': 'state-independent canonical commutator of the split unitary only',
                  'quantum_covariance_or_probability_queried': False,
                  'modified_hamiltonian_vacuum_substituted': False,
                  'continuous_action_comparison': 'all64modes at all21 executed times; no continuum error claim'},
        'source_pins': pins}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    raw = canonical(build())
    if args.check:
        if OUTPUT.read_bytes() != raw:
            raise SystemExit('source execution receipt drift')
    else:
        OUTPUT.write_bytes(raw)
    print('authenticated source scalar execution PASS')


if __name__ == '__main__':
    main()
