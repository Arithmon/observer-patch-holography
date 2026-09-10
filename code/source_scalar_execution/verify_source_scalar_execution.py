"""Independent exact replay: sqrt(5) arithmetic, edge energy assembly, Horner.

No producer or shared algebra code is imported. The proof is in the pinned
SOURCE_SCALAR_EXECUTION fragment; the replay checks its finite hypotheses and
all retained numerical and authenticated consequences.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import cmp_to_key
from hashlib import sha256
from itertools import product
from math import factorial
from pathlib import Path
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'source_scalar_execution_receipt.json'
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


def require(condition, message):
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class R:
    """Independent a+b*sqrt(5) field representation."""
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        require(type(self.a) in (int, F) and type(self.b) in (int, F), 'exact coefficient')
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def of(x):
        return x if type(x) is R else R(x)

    def __add__(self, other):
        y = R.of(other)
        return R(self.a+y.a, self.b+y.b)

    __radd__ = __add__

    def __neg__(self):
        return R(-self.a, -self.b)

    def __sub__(self, other):
        return self+-R.of(other)

    def __rsub__(self, other):
        return R.of(other)+-self

    def __mul__(self, other):
        y = R.of(other)
        return R(self.a*y.a+5*self.b*y.b, self.a*y.b+self.b*y.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        y = R.of(other)
        norm = y.a*y.a-5*y.b*y.b
        require(norm != 0, 'zero divisor')
        return self*R(y.a/norm, -y.b/norm)

    def __rtruediv__(self, other):
        return R.of(other)/self

    def __pow__(self, n):
        require(type(n) is int and n >= 0, 'exact nonnegative power')
        z, x = R(1), self
        while n:
            if n % 2:
                z *= x
            x *= x
            n //= 2
        return z

    def sign(self):
        if self.b == 0:
            return (self.a > 0)-(self.a < 0)
        if self.a == 0:
            return (self.b > 0)-(self.b < 0)
        if self.a*self.b > 0:
            return 1 if self.a > 0 else -1
        d = self.a*self.a-5*self.b*self.b
        return ((d > 0)-(d < 0))*(1 if self.a > 0 else -1)

    def encode(self):
        return [str(self.a-self.b), str(2*self.b)]


TAU = R(0, F(1, 35))


def rational(x):
    require(type(x) is str and len(x) <= 3000, 'rational string type/size')
    try:
        result = F(x)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError('invalid rational') from exc
    require(str(result) == x and result.numerator.bit_length() <= 10000
            and result.denominator.bit_length() <= 10000, 'canonical bounded rational')
    return result


def parse(value):
    require(type(value) is list and len(value) == 2, 'two algebraic coefficients')
    a, b = map(rational, value)
    return R(a+b/2, b/2)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'),
                       ensure_ascii=True, allow_nan=False)+'\n').encode('ascii')


def equal(actual, expected, label):
    # Canonical JSON distinguishes booleans, integers and strings recursively.
    require(canonical(actual) == canonical(expected), label)


def strict_pairs(pairs):
    out = {}
    for key, value in pairs:
        require(key not in out, 'duplicate JSON key')
        out[key] = value
    return out


def forbidden_number(token):
    raise ValueError('floating/nonfinite JSON token forbidden')


def load(path=OUTPUT):
    raw = Path(path).read_bytes()
    require(len(raw) <= 12_000_000, 'receipt size budget')
    return json.loads(raw.decode('ascii'), object_pairs_hook=strict_pairs,
                      parse_float=forbidden_number, parse_constant=forbidden_number)


def table(root, relative, name):
    text = (root/relative).read_text(encoding='utf-8')
    match = re.search(r'def '+name+r'\s*:[^=]*:=\s*!\[(.*?)\]', text, re.S)
    require(match is not None, 'missing source table '+name)
    return [int(v.strip()) for v in match[1].split(',')]


def primitive_basis(root):
    seam = 'Lean/Screen/SeamCurrentCarrierQuotient.lean'
    frame = 'Lean/Screen/PrimitivePortFrameQuotient.lean'
    left, right = table(root, seam, 'seamLeft'), table(root, seam, 'seamRight')
    relabel = table(root, frame, 'sourceToRERPort')
    positive = [relabel[i] for i in table(root, frame, 'positiveSourcePort')]
    axes = ((R(-2), R(1, 1), R()), (R(2), R(1, 1), R()),
            (R(), R(-2), R(1, 1)), (R(), R(2), R(1, 1)),
            (R(1, 1), R(), R(-2)), (R(1, 1), R(), R(2)))
    raw = (root/'Lean/Screen/PortFrameGram.lean').read_text(encoding='utf-8')
    neighbors = {int(i): [int(v) for v in row.split(',')]
                 for i, row in re.findall(r'\| (\d+) => \[([\d, ]+)\]', raw)}
    require(len(neighbors) == 12 and len(left) == len(right) == 30, 'primitive dimensions')
    for i, u in enumerate(positive):
        for j, v in enumerate(positive):
            gram5 = R(5) if u == v else (R(-5) if u+v == 11 else
                                        (R(0, 1) if v in neighbors[u] else R(0, -1)))
            require(sum((x*y for x, y in zip(axes[i], axes[j])), R())/4
                    == gram5*R(F(1, 2), F(1, 10)), 'Gram basis mismatch')
    return left, right, positive, axes


def model(root=ROOT):
    orbit = []
    for b in range(5):
        a = 0
        while (R(F(b, 2)-a-1, F(b, 2))).sign() >= 0:
            a += 1
        orbit.append(R(F(b, 2)-a, F(b, 2)))
    x = sorted(orbit, key=cmp_to_key(lambda a, b: (a-b).sign()))+[R(1)]
    gap = [x[i+1]-x[i] for i in range(5)]
    weights = [(x[i+2]-x[i])/2 for i in range(4)]
    sites = list(product(range(4), repeat=3))
    indices = {site: i for i, site in enumerate(sites)}
    mass = [weights[a]*weights[b]*weights[c] for a, b, c in sites]
    # Assemble the symmetric stiffness from undirected edge energies plus
    # the boundary and mass terms, then divide its rows by the dual masses.
    stiffness = [{i: m} for i, m in enumerate(mass)]
    for i, site in enumerate(sites):
        for axis in range(3):
            cross = R(1)
            for d in range(3):
                if d != axis:
                    cross *= weights[site[d]]
            k = site[axis]
            if k == 0:
                stiffness[i][i] += cross/gap[0]
            if k == 3:
                stiffness[i][i] += cross/gap[4]
            else:
                other = list(site)
                other[axis] += 1
                j = indices[tuple(other)]
                c = cross/gap[k+1]
                stiffness[i][i] += c
                stiffness[j][j] += c
                stiffness[i][j] = -c
                stiffness[j][i] = -c
    rows = [sorted((j, c/mass[i]) for j, c in row.items()) for i, row in enumerate(stiffness)]
    velocity, detector, addresses = [], [], []
    left, right, positive, axes = primitive_basis(root)
    for site in sites:
        xyz = [x[k+1] for k in site]
        product_value = R(1)
        for z in xyz:
            product_value *= z*(1-z)
        require((xyz[0]-F(9, 20)).sign() != 0 and (xyz[0]-F(11, 20)).sign() != 0,
                'ambiguous support boundary')
        velocity.append(256*product_value if (xyz[0]-F(9, 20)).sign() < 0 else R())
        detector.append(64*product_value if (xyz[0]-F(11, 20)).sign() > 0 else R())
        a = [int(z.a-z.b) for z in xyz]
        b = [int(2*z.b) for z in xyz]
        currents = [-a[2]-b[0], a[0]+b[0]+b[2], b[0]+b[1]-b[2],
                    -a[0]-b[0]+b[2], a[2]+b[1]-b[2], a[1]+b[2]]
        load = [0]*12
        for seam, value in zip((5, 8, 12, 13, 22, 24), currents):
            load[left[seam]] -= value
            load[right[seam]] += value
        require(sum(load) == 0, 'source conservation')
        z = [load[p]-load[11-p] for p in positive]
        require(sum(z) % 2 == 0, 'D6 parity')
        for axis in range(3):
            require(sum((c*ax[axis] for c, ax in zip(z, axes)), R()) == 4*xyz[axis],
                    'primitive source readback')
        addresses.append({'site': list(site), 'coordinate_Qphi': [v.encode() for v in xyz],
                          'D6_control': z, 'signed_seam_currents': currents})
    require(TAU**2 == R(F(1, 245)), 'step square')
    require(all((5*h*h-1).sign() < 0 for h in gap), 'source edges exceed old radius')
    require(all((h-F(7, 50)).sign() > 0 for h in gap), 'minimum gap')
    require(F(30049, 49) < 614 and F(30049, 12005) < 3, 'CFL bound')
    require((sum(mass, R())-1).sign() <= 0, 'mass budget')
    require(all(v.sign() >= 0 and (v-4).sign() <= 0 for v in velocity), 'velocity budget')
    require(all(v.sign() >= 0 and (v-1).sign() <= 0 for v in detector), 'detector budget')
    return mass, rows, velocity, detector, addresses


def apply(rows, vector):
    return [sum((c*vector[j] for j, c in row), R()) for row in rows]


def inner(mass, x, y):
    return sum((m*a*b for m, a, b in zip(mass, x, y)), R())


def verify_trace(trace, rows, velocity, schedule):
    require(type(trace) is dict and set(trace) == {'schedule', 'events', 'layers', 'final_hash',
            'seed_writes', 'update_writes', 'dynamic_field_reads'}, 'trace schema')
    equal(trace['schedule'], schedule, 'schedule census')
    require(type(trace['events']) is list and len(trace['events']) == 1472, 'event census')
    for key, value in [('seed_writes', 128), ('update_writes', 1344), ('dynamic_field_reads', 8736)]:
        equal(trace[key], value, 'execution costs')
    registers, ancestors, chain = {}, {}, '0'*64
    layers = [[R() for _ in range(64)] for _ in range(23)]
    reads_total = 0
    for position, event in enumerate(trace['events']):
        phase, offset = divmod(position, 64)
        site = schedule[offset]
        require(type(event) is dict and set(event) == {'id', 'reads', 'write', 'ledger_parent', 'hash'},
                'event schema')
        equal(event['id'], [phase, site], 'event identity')
        equal(event['ledger_parent'], chain, 'serial commitment parent')
        require(type(event['reads']) is list, 'read list')
        expected_reads = []
        if phase < 2:
            value = -TAU*velocity[site] if phase == 0 else R()
            resource, version = (1-phase)*64+site, 1
        else:
            step = phase-1
            resource = (step % 2)*64+site
            old_version = (step+1)//2
            keys = [(resource, old_version)]+[
                (((step-1) % 2)*64+j, step//2+1) for j, _ in rows[site]]
            read_values = []
            for key in keys:
                require(key in registers, 'read before actual writer')
                writer, rvalue = registers[key]
                expected_reads.append([*key, writer, rvalue.encode()])
                read_values.append(rvalue)
            # Compute from the authenticated read payload, including full row.
            old = read_values[0]
            local = next(v for (j, _), v in zip(rows[site], read_values[1:]) if j == site)
            force = sum((c*v for (_, c), v in zip(rows[site], read_values[1:])), R())
            value = 2*local-old-R(F(1, 245))*force
            version = old_version+1
        equal(event['reads'], expected_reads, 'complete exact read/writer/version/value set')
        equal(event['write'], [resource, version, [phase, site], value.encode()], 'exact field write')
        require((resource, version) not in registers, 'immutable version overwritten')
        registers[(resource, version)] = ([phase, site], value)
        layers[phase][site] = value
        parents = {tuple(read[2]) for read in expected_reads}
        past = set(parents)
        for parent in parents:
            require(parent in ancestors, 'missing semantic parent')
            past.update(ancestors[parent])
        ancestors[(phase, site)] = past
        material = {k: event[k] for k in ('id', 'reads', 'write', 'ledger_parent')}
        chain = sha256(canonical(material)).hexdigest()
        equal(event['hash'], chain, 'event byte commitment')
        reads_total += len(expected_reads)
    equal(trace['final_hash'], chain, 'final commitment')
    equal(trace['layers'], [[v.encode() for v in row] for row in layers], 'field readback layers')
    require(reads_total == 8736, 'complete read bill')
    return layers, ancestors


def value_bounds(x, scale=10**12):
    # Bracket algebraic values independently by monotone binary search.
    l, r = -scale, scale
    while (x-R(F(l, scale))).sign() < 0:
        l *= 2
    while (x-R(F(r, scale))).sign() > 0:
        r *= 2
    while r-l > 1:
        m = (l+r)//2
        if (x-R(F(m, scale))).sign() >= 0:
            l = m
        else:
            r = m
    return F(l, scale), F(r, scale)


def norm_bounds(square, scale=10**12):
    require(square.sign() >= 0, 'nonnegative norm')
    l, r = 0, scale
    while (R(F(r, scale)**2)-square).sign() < 0:
        r *= 2
    while r-l > 1:
        m = (l+r)//2
        if (R(F(m, scale)**2)-square).sign() <= 0:
            l = m
        else:
            r = m
    return F(l, scale), F(r, scale)


def widen(bounds, radius):
    scale = 10**12
    low, high = bounds[0]-radius, bounds[1]+radius
    return [str(F(low*scale//1, scale)), str(F(-(-high*scale//1), scale))]


def diagnostics(mass, rows, velocity, detector, layers, ancestors):
    e0 = inner(mass, velocity, velocity)/2
    for previous, current in zip(layers, layers[1:]):
        difference = [b-a for a, b in zip(previous, current)]
        energy = inner(mass, difference, difference)*F(245, 2)+inner(mass, current, apply(rows, previous))/2
        require(energy == e0, 'modified energy identity')
    powers = [velocity]
    for _ in range(80):
        powers.append(apply(rows, powers[-1]))
    tail = F(4)*F(3, 2)**163*614**81/factorial(163)
    require(tail < F(1, 10**30), 'full spectral remainder')
    comparisons = []
    coefficients = [F((-1)**k, factorial(2*k+1)) for k in range(81)]
    for j in range(1, 22):
        square_time = F(j*j, 245)
        polynomial = []
        for i in range(64):
            value = R()
            for k in range(80, -1, -1):
                value = value*square_time+coefficients[k]*powers[k][i]
            polynomial.append(value*(j*TAU))
        actual = layers[j+1]
        diff = [a-b for a, b in zip(actual, polynomial)]
        detector_value = inner(mass, detector, actual)
        comparisons.append({'step': j, 'model_time_Qphi': (j*TAU).encode(),
            'executed_detector_Qphi': detector_value.encode(),
            'canonical_commutator_over_i_hbar_Qphi': (-detector_value/4).encode(),
            'continuous_detector_interval': widen(value_bounds(inner(mass, detector, polynomial)), tail),
            'full_mass_norm_discretization_error': widen(norm_bounds(inner(mass, diff, diff)), tail)})
    counts = []
    lower = (1, 21)
    for j in (7, 14, 21):
        upper = (j+1, 21)
        count = sum(e == lower or lower in ancestors[e] for e in ancestors[upper] | {upper})
        counts.append({'lower': list(lower), 'upper': list(upper), 'inclusive_count': count})
    response = layers[8][53]
    require(response.sign() != 0 and (R(F(3, 2), F(-1, 2))-F(1, 5)).sign() > 0,
            'nonzero old-radius mismatch')
    return {'modified_energy_Qphi': e0.encode(), 'modified_energy_interval': [str(v) for v in value_bounds(e0)],
        'energy_equal_at_all_22_states': True, 'taylor_remainder_upper': str(tail), 'comparisons': comparisons,
        'actual_field_interval_counts': counts,
        'first_coarse_layer_outside_radius': {'target_site': 53, 'step': 7,
            'source_plane_x_Qphi': R(-2, 1).encode(), 'target_x_Qphi': R(F(-5, 2), F(3, 2)).encode(),
            'old_radius_squared': '1/5', 'x_displacement_squared_Qphi': R(F(3, 2), F(-1, 2)).encode(),
            'field_response_Qphi': response.encode(), 'response_interval': [str(v) for v in value_bounds(response)],
            'all_nonzero_preparation_sites_outside_old_radius': True}}


MODEL = {'q': 5, 'mutable_sites': 64, 'steps': 21, 'coarse_layers': 3, 'substeps_per_coarse_layer': 7,
    'time_step_Qphi': ['-1/35', '2/35'], 'time_step_squared': '1/245',
    'source_length': 'L=2/sqrt(phi+2)', 'time_unit': 'c*t/L', 'dimensionless_mass': '1',
    'boundary': 'Dirichlet cube', 'registers': 'two alternating immutable-version field registers per site',
    'previous_seed': 'integrator memory -tau*v0; not a captured negative-time solution',
    'preparation': 'v0=256*x*y*z*(1-x)*(1-y)*(1-z) on x<9/20; zero elsewhere',
    'detector': '64*x*y*z*(1-x)*(1-y)*(1-z) on x>11/20; zero elsewhere',
    'cfl_upper': '30049/12005', 'cfl_strict_margin': 'less than3, hence less than4'}
SCOPE = {'all_substeps_and_field_reads_stored': True, 'source_addresses_protected_preparation': True,
    'all_arithmetic_exact_Qphi': True, 'action_boundary_model_time_supplied': True,
    'source_action_family_preserved': True, 'same_q233_preparation_or_accuracy': False,
    'read_values_writers_and_immutable_versions_replayed': True,
    'serial_audit_parent_is_a_semantic_read_edge': False, 'old_count_clock_theorem_applied': False,
    'native_repair_or_physical_clock_selected': False,
    'quantum_scope': 'state-independent canonical commutator of the split unitary only',
    'quantum_covariance_or_probability_queried': False, 'modified_hamiltonian_vacuum_substituted': False,
    'continuous_action_comparison': 'all64modes at all21 executed times; no continuum error claim'}


def verify(packet, root=ROOT):
    require(type(packet) is dict and set(packet) == {'schema', 'model', 'source_addresses', 'mass_Qphi',
            'action_rows_Qphi', 'initial_velocity_Qphi', 'detector_Qphi', 'traces', 'diagnostics',
            'scope', 'source_pins'}, 'packet schema')
    equal(packet['schema'], 'oph.authenticated-golden-scalar-execution.v1', 'schema version')
    equal(packet['model'], MODEL, 'declared finite model')
    equal(packet['scope'], SCOPE, 'scientific scope')
    expected_pins = {p: {'sha256': sha256((root/p).read_bytes()).hexdigest(), 'bytes': (root/p).stat().st_size}
                     for p in PINS}
    equal(packet['source_pins'], expected_pins, 'source byte custody')
    mass, rows, velocity, detector, addresses = model(root)
    equal(packet['source_addresses'], addresses, 'primitive source addresses')
    equal(packet['mass_Qphi'], [v.encode() for v in mass], 'dual masses')
    equal(packet['action_rows_Qphi'], [[[j, c.encode()] for j, c in row] for row in rows], 'edge action')
    equal(packet['initial_velocity_Qphi'], [v.encode() for v in velocity], 'prepared velocity')
    equal(packet['detector_Qphi'], [v.encode() for v in detector], 'supplied detector')
    traces = packet['traces']
    require(type(traces) is dict and set(traces) == {'ascending_baseline', 'descending_baseline',
            'ascending_intervention', 'descending_intervention'}, 'four schedule/preparation controls')
    states = {}
    for name in sorted(traces):
        schedule = list(range(64)) if name.startswith('ascending') else list(reversed(range(64)))
        v = velocity if name.endswith('intervention') else [R() for _ in range(64)]
        states[name] = verify_trace(traces[name], rows, v, schedule)
    for preparation in ('baseline', 'intervention'):
        require(states['ascending_'+preparation][0] == states['descending_'+preparation][0],
                'schedule independence')
    values, ancestors = states['ascending_intervention']
    expected = diagnostics(mass, rows, velocity, detector, values, ancestors)
    equal(packet['diagnostics'], expected, 'independent full-spectrum/ancestry diagnostics')
    return {'verdict': 'PASS', 'mutable_sites': 64, 'executed_times': 21,
            'traces_replayed': 4, 'events_replayed': 5888, 'dynamic_field_reads_replayed': 34944,
            'full_mass_norm_error_final': expected['comparisons'][-1]['full_mass_norm_discretization_error'],
            'same_finite_operator_all_modes': True, 'old_count_clock_applied': False,
            'quantum_covariance_or_probability_queried': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.path)), sort_keys=True))


if __name__ == '__main__':
    main()
