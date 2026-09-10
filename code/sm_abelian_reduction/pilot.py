"""Local Hamiltonian operations on 64 prepared golden interior addresses.

Binary64 execution; exact preservation statements concern ideal subflows.
Real connection values are retained, never reconstructed from compact phases.
"""
import cmath
from copy import deepcopy
from hashlib import sha256
from itertools import product, combinations
import json
import math

H = 1/200
STEPS = 2
TARGET = 21


def canonical(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def digest(x):
    return sha256(canonical(x)).hexdigest()


def geometry():
    phi = (1+math.sqrt(5))/2
    order = sorted(range(5), key=lambda b: b*8 % 5)
    orbit = [(-((b+math.isqrt(5*b*b))//2), b) for b in order]
    points = [a+b*phi for a, b in orbit]+[1.0]
    gaps = [b-a for a, b in zip(points, points[1:])]
    weights = [(a+b)/2 for a, b in zip(gaps, gaps[1:])]
    coords = list(product(range(1, 5), repeat=3)); lookup = {x: i for i, x in enumerate(coords)}
    sites = []
    for xyz in coords:
        records = [orbit[k] for k in xyz]; a, b = zip(*records)
        z = [b[1]-a[0], b[1]+a[0], b[2]-a[1], b[2]+a[1], b[0]-a[2], b[0]+a[2]]
        half = sum(z)//2
        currents = [-z[5], half-z[0], half-z[2]-z[3], half-z[1]-z[4]-z[5], half-z[2]-z[3]-z[4], z[3]]
        volume = math.prod(weights[k-1] for k in xyz)
        boundary = float(sum(math.prod(weights[xyz[j]-1] for j in range(3) if j != axis)/gaps[0 if k == 1 else 4]
                       for axis, k in enumerate(xyz) if k in (1, 4)))
        sites.append({'index': list(xyz), 'address_Qphi': [list(v) for v in records],
                      'source_record': z, 'source_seam_currents': currents,
                      'mass': volume, 'dirichlet_stiffness': boundary})
    edges = []; edge_lookup = {}
    for i, xyz in enumerate(coords):
        for axis in range(3):
            if xyz[axis] == 4: continue
            other = list(xyz); other[axis] += 1
            c = math.prod(weights[xyz[j]-1] for j in range(3) if j != axis)/gaps[xyz[axis]]
            edge_lookup[(xyz, axis)] = len(edges)
            edges.append({'ends': [i, lookup[tuple(other)]], 'axis': axis,
                          'conductance': c, 'electric_mass': (8/3)*c})
    faces = []
    for xyz in coords:
        for a, b in combinations(range(3), 2):
            if xyz[a] == 4 or xyz[b] == 4: continue
            xa, xb = list(xyz), list(xyz); xa[a] += 1; xb[b] += 1
            boundary = [[edge_lookup[(xyz, a)], 1], [edge_lookup[(tuple(xa), b)], 1],
                        [edge_lookup[(tuple(xb), a)], -1], [edge_lookup[(xyz, b)], -1]]
            normal = 3-a-b
            faces.append({'boundary': boundary,
                          'magnetic_weight': (8/3)*weights[xyz[normal]-1]/(gaps[xyz[a]]*gaps[xyz[b]])})
    return {'q': 5, 'fibonacci_index': 5, 'orbit_order': order,
            'sites': sites, 'edges': edges, 'plaquettes': faces}


def encode(z):
    return [float(z.real), float(z.imag)] if isinstance(z, complex) else float(z)


def execute(geo, cohort='baseline', reverse_kick=False):
    if cohort not in ('baseline', 'phase_intervention', 'gauge_copy'):
        raise ValueError('unknown cohort')
    values = {}; writers = {}; versions = {}; events = []; previous = '0'*64
    def event(op, reads, writes, meta):
        nonlocal previous
        names = sorted(set(reads))
        record = {'id': len(events), 'operation': op, 'metadata': meta,
                  'reads': [{'port': k, 'version': versions[k], 'writer': writers[k],
                             'value': encode(values[k])} for k in names],
                  'parents': sorted({writers[k] for k in names}),
                  'writes': [], 'previous_hash': previous}
        for k, v in sorted(writes.items()):
            values[k] = v; versions[k] = versions.get(k, -1)+1; writers[k] = record['id']
            record['writes'].append({'port': k, 'version': versions[k], 'value': encode(v)})
        record['event_hash'] = digest(record); previous = record['event_hash']; events.append(record)
    for i in range(64):
        event('prepare_site', [], {f'psi/{i}': complex(1/5), f'pi/{i}': 0j}, {'site': i})
    for e in range(144):
        event('prepare_link', [], {f'a/{e}': 0., f'P/{e}': 0.}, {'edge': e})
    if cohort != 'baseline':
        names = [f'psi/{TARGET}', f'pi/{TARGET}']; writes = {k: values[k]*complex(3/5, 4/5) for k in names}
        if cohort == 'gauge_copy':
            angle = math.atan2(4, 3)
            for e, row in enumerate(geo['edges']):
                i, j = row['ends']
                if TARGET in (i, j):
                    key = f'a/{e}'; names.append(key); writes[key] = values[key]+angle*((i == TARGET)-(j == TARGET))
        event(cohort, names, writes, {'site': TARGET, 'unit_phase': ['3/5', '4/5'], 'gauge_parameter': 'atan2(4,3)' if cohort == 'gauge_copy' else None})
    checkpoints = []
    def checkpoint(step, phase):
        names = sorted(values)
        snapshot = {k: encode(values[k]) for k in names}
        event('readout', names, {}, {'step': step, 'phase': phase})
        checkpoints.append({'event': len(events)-1, 'step': step, 'phase': phase, 'values': snapshot})
    def kick(step, part):
        # These potential subflows commute: positions/connections stay fixed.
        order = list(range(144))
        if reverse_kick: order.reverse()
        for e in order:
            row = geo['edges'][e]; i, j = row['ends']; c = row['conductance']
            qi, qj, pi, pj, a, P = f'psi/{i}', f'psi/{j}', f'pi/{i}', f'pi/{j}', f'a/{e}', f'P/{e}'
            U = cmath.exp(1j*values[a]); delta = U*values[qj]-values[qi]
            J = -2*c*(values[qi].conjugate()*U*values[qj]).imag
            event('edge_kick', [qi, qj, pi, pj, a, P],
                  {pi: values[pi]+H*c*delta/2, pj: values[pj]-H*c*U.conjugate()*delta/2,
                   P: values[P]+H*J/2}, {'edge': e, 'step': step, 'part': part})
        for i, row in enumerate(geo['sites']):
            q, p = f'psi/{i}', f'pi/{i}'
            gradient = (row['dirichlet_stiffness']+row['mass']*(1+abs(values[q])**2/2))*values[q]
            event('onsite_kick', [q, p], {p: values[p]-H*gradient/2}, {'site': i, 'step': step, 'part': part})
        for f, row in enumerate(geo['plaquettes']):
            curl = sum(sign*values[f'a/{e}'] for e, sign in row['boundary'])
            names = [f'{p}/{e}' for e, _ in row['boundary'] for p in ('a', 'P')]
            event('magnetic_kick', names,
                  {f'P/{e}': values[f'P/{e}']-H*row['magnetic_weight']*sign*curl/2 for e, sign in row['boundary']},
                  {'plaquette': f, 'step': step, 'part': part})
    checkpoint(0, 'initial')
    for step in range(1, STEPS+1):
        kick(step, 'first')
        if step == 1: checkpoint(step, 'first_kick')
        for i, row in enumerate(geo['sites']):
            q, p = f'psi/{i}', f'pi/{i}'
            event('scalar_drift', [q, p], {q: values[q]+H*values[p]/row['mass']}, {'site': i, 'step': step})
        for e, row in enumerate(geo['edges']):
            a, P = f'a/{e}', f'P/{e}'
            event('electric_drift', [a, P], {a: values[a]+H*values[P]/row['electric_mass']}, {'edge': e, 'step': step})
        kick(step, 'last'); checkpoint(step, 'completed')
    return {'cohort': cohort, 'events': events, 'checkpoints': checkpoints, 'final_event_hash': previous}


def observables(geo, snapshot):
    z = lambda k: complex(*snapshot[k])
    charge = [2*(z(f'psi/{i}').conjugate()*z(f'pi/{i}')).imag for i in range(64)]
    gauss = charge.copy(); currents = []; energy = 0.
    for i, row in enumerate(geo['sites']):
        r = abs(z(f'psi/{i}'))**2
        energy += abs(z(f'pi/{i}'))**2/row['mass']+(row['mass']+row['dirichlet_stiffness'])*r+row['mass']*r*r/4
    for e, row in enumerate(geo['edges']):
        i, j = row['ends']; a, P = snapshot[f'a/{e}'], snapshot[f'P/{e}']; U = cmath.exp(1j*a)
        J = -2*row['conductance']*(z(f'psi/{i}').conjugate()*U*z(f'psi/{j}')).imag
        currents.append(J); gauss[i] += P; gauss[j] -= P
        energy += P*P/(2*row['electric_mass'])+row['conductance']*abs(U*z(f'psi/{j}')-z(f'psi/{i}'))**2
    curls = [sum(sign*snapshot[f'a/{e}'] for e, sign in row['boundary']) for row in geo['plaquettes']]
    energy += sum(row['magnetic_weight']*v*v/2 for row, v in zip(geo['plaquettes'], curls))
    return {'charge': charge, 'current': currents, 'gauss_residual': gauss,
            'scalar_intensity': [abs(z(f'psi/{i}'))**2 for i in range(64)],
            'energy': energy, 'max_weak_plaquette_angle': max(abs(4*v/3) for v in curls)}
