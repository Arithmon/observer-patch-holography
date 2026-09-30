"""Independent analytic/entrywise replay. No producer imports."""

from fractions import Fraction
import itertools
import math

import mpmath as mp
import numpy as np

from source_selection_model import verify_response as response
from source_selection_model.geometry import check_tower, oriented_key
from m1_operational_clocks.check import clock_interval_certificate


def need(condition, message):
    if not condition:
        raise ValueError(message)


def same(actual, expected, path='evidence'):
    need(type(actual) is type(expected), f'{path}: JSON type')
    if isinstance(expected, dict):
        need(actual.keys() == expected.keys(), f'{path}: keys')
        for k in expected:
            same(actual[k], expected[k], f'{path}.{k}')
    elif isinstance(expected, list):
        need(len(actual) == len(expected), f'{path}: length')
        for k, (a, e) in enumerate(zip(actual, expected)):
            same(a, e, f'{path}[{k}]')
    elif isinstance(expected, float):
        need(math.isfinite(actual) and math.isfinite(expected), f'{path}: nonfinite')
        # Only numerical matrix entries/eigenvalue zeros and residuals have
        # an absolute tolerance. Positive physical scales cannot become zero.
        absolute = any(x in path for x in ('eigenvalues', 'matrix_unit_images', 'diagonal',
            'three_step_probe', 'amplitudes', '.output', 'concurrence', 'error'))
        tol = 2e-10 if absolute else (2e-9*abs(expected) if expected else 1e-12)
        need(abs(actual-expected) <= tol, f'{path}: mismatch {actual} != {expected}')
    else:
        need(actual == expected, f'{path}: mismatch')


def encoded(z):
    return np.stack((np.asarray(z).real, np.asarray(z).imag), axis=-1).tolist()


def matrix(value, shape):
    need(type(value) is list, 'matrix list')
    a = np.asarray(value)
    need(a.shape == shape+(2,) and a.dtype.kind == 'f' and np.all(np.isfinite(a)), 'complex matrix schema')
    # Reject mixed integer/boolean cells even when NumPy would coerce them.
    def leaves(x):
        return all(leaves(v) for v in x) if type(x) is list else type(x) is float and math.isfinite(x)
    need(leaves(value), 'complex entries must be JSON floats')
    return a[..., 0]+1j*a[..., 1]


def verify_controls(row):
    need(type(row) is dict and row.keys() == {'gram_determinant', 'rows'}, 'control schema')
    gs = response.reconstructed_generators()
    targets = []
    for i in range(3):
        h = np.zeros((6, 6), complex)
        h[i, i] = 1j
        targets.append((f'phase_{i}', h))
    for i, j in itertools.combinations(range(3), 2):
        for kind in ('real', 'imaginary'):
            h = np.zeros((6, 6), complex)
            h[i, j], h[j, i] = (1, -1) if kind == 'real' else (1j, 1j)
            targets.append((f'{kind}_{i}{j}', h))
    targets.append(('relative_phase', np.diag([1j]*3+[0]*3)))
    need(type(row['rows']) is list and len(row['rows']) == len(targets), 'control coverage')
    for actual, (name, h) in zip(row['rows'], targets):
        need(type(actual) is dict and actual.keys() == {'name', 'coefficients'}, 'control row')
        same(actual['name'], name)
        need(type(actual['coefficients']) is list and len(actual['coefficients']) == 12, 'twelve controls')
        full = response.combine(gs, [response.real(x) for x in actual['coefficients']])
        for b, i, j in itertools.product(range(2), range(3), range(3)):
            z = h[3*b+i, 3*b+j]
            need(full[b][i][j] == response.C5(response.F5(int(z.real)), response.F5(int(z.imag))),
                 'exact closed-form source reconstruction')
    same(row['rows'][-1]['coefficients'], [['1', '0']]*12)
    gram = [[-sum((gs[p][b][i][j]*gs[q][b][j][i]
                   for b, i, j in itertools.product(range(2), range(3), range(3))), response.CZ).re
             for q in range(12)] for p in range(12)]
    det = response.O
    for p in range(12):
        pivot = gram[p][p]
        need(pivot != response.Z, 'Gram pivot')
        det *= pivot
        for i in range(p+1, 12):
            f = gram[i][p]/pivot
            for j in range(p+1, 12):
                gram[i][j] -= f*gram[p][j]
    need(response.real(row['gram_determinant']) == det, 'exact Gram determinant')


def expected_channels():
    rows = []
    for d in (2, 3, 4):
        for ez, ex in ((0., 0.), (0., (d-1)/d), (.02, .03), (0., .05)):
            p = [1-ez]+[ez/(d-1)]*(d-1)
            q = [1-ex]+[ex/(d-1)]*(d-1)
            values = [x*y for x in p for y in q]
            h = lambda r: -sum(x*math.log(x) for x in r if x)
            rows.append(dict(d=d, ez=ez, ex=ex, choi_eigenvalues=sorted(values),
                z_fidelity=1-ez, x_fidelity=1-ex, entanglement_fidelity=(1-ez)*(1-ex),
                entropy=h(p)+h(q), choi_half_trace_distance=ez+ex-ez*ex))
    return rows


def expected_instruments():
    rows = []
    for d in (2, 4):
        for coherent in (False, True):
            images = []
            for i, j in itertools.product(range(6), repeat=2):
                z = np.zeros((d+1, d+1), complex)
                if i < d and j < d and (coherent or i == j):
                    z[i, j] = 1
                elif i == j and i >= d:
                    z[d, d] = 1
                images.append(encoded(z.reshape(-1)))
            rows.append(dict(d=d, coherent=coherent, completeness_error=0., matrix_unit_images=images))
    return rows


def expected_entangler():
    return [dict(theta=t, code_diagonal=encoded(np.array([1, 1, 1, complex(math.cos(t), math.sin(t))])),
                 concurrence=abs(math.sin(t/2)), marginal_purity=1-math.sin(t/2)**2/2)
            for t in (0., math.pi/7, math.pi/2, math.pi)]


def coin():
    c = np.zeros((4, 4), complex)
    c[0, 1:] = c[1:, 0] = 1/math.sqrt(3)
    for i, j in ((1, 2), (2, 3), (3, 1)):
        c[i, j], c[j, i] = 1j/math.sqrt(3), -1j/math.sqrt(3)
    return c


def verify_coin(row):
    need(type(row) is dict and row.keys() == {'diagonal', 'stages', 'reconstruction_error'}, 'coin schema')
    d = matrix(row['diagonal'], (4,))
    need(np.max(abs(abs(d)-1)) < 2e-12, 'all four final phases')
    result = np.diag(d)
    need(type(row['stages']) is list and len(row['stages']) == 6, 'six rotations')
    for step, pair in zip(row['stages'], ([2, 3], [1, 2], [2, 3], [0, 1], [1, 2], [2, 3])):
        need(type(step) is dict and step.keys() == {'pair', 'matrix', 'native_unitarity_error'}, 'coin stage')
        same(step['pair'], pair)
        g = matrix(step['matrix'], (2, 2))
        need(np.linalg.norm(g.conj().T@g-np.eye(2)) < 2e-12 and abs(np.linalg.det(g)-1) < 2e-12,
             'native SU2 rotation')
        result[pair, :] = g@result[pair, :]
        same(step['native_unitarity_error'], 0., 'unitarity_error')
    need(np.linalg.norm(result-coin()) < 2e-12, 'entire coin including final phases')
    same(row['reconstruction_error'], 0., 'reconstruction_error')


def closed_spatial(q, a, mass):
    """Direct entries from the two possible flight paths, without QR/gates."""
    points = list(itertools.product(range(q), repeat=3))
    signs = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    c = coin()
    out = np.zeros((8*q**3, 8*q**3), complex)
    mu = mass*a/math.sqrt(3)
    for yi, y in enumerate(points):
        for incoming, p, r in itertools.product(range(8), range(4), range(2)):
            f, s = divmod(incoming, 4)
            delta = signs[p if f == 0 else s]
            x = tuple((y[j]+(1 if f == 0 else -1)*delta[j]) % q for j in range(3))
            xi = (x[0]*q+x[1])*q+x[2]
            out[8*xi+4*r+p, 8*yi+incoming] += c[p, s]*(math.cos(mu) if r == f else -1j*math.sin(mu))
    return out


def expected_spatial():
    rows = []
    for q, a, mass in itertools.product((2, 3), (.07, .013), (.7, 1.1)):
        u = closed_spatial(q, a, mass)
        size = len(u)
        probe = np.array([complex(math.cos(j*j/17), math.sin(j*j/17)) for j in range(size)])/math.sqrt(size)
        out = np.linalg.matrix_power(u, 3)@probe
        need(np.linalg.norm(u.conj().T@u-np.eye(size)) < 1e-11, 'complete spatial unitarity')
        rows.append(dict(q=q, a=a, mass=mass, modes=size, instructions=25, pulses_per_cell=24,
                         code_events_per_cell=64, flights_per_cell=8, column_norms=[1.]*size,
                         three_step_probe=encoded(out)))
    return rows


def expected_preparations():
    rows = []
    for n in (2, 4, 8):
        for holes in (False, True):
            target = np.zeros((n, n, n), complex)
            for i, j, k in itertools.product(range(n), repeat=3):
                if not holes or max(i, j, k) >= n//2:
                    r = sum((z-(n-1)/2)**2 for z in (i, j, k))
                    theta = (i+2*j-k)/3
                    target[i, j, k] = math.exp(-r/n)*complex(math.cos(theta), math.sin(theta))
            target /= np.linalg.norm(target)
            rows.append(dict(side=n, holes=holes, pulses=2*n**3-1, flights=8*(n**3-1)//7,
                             amplitudes=encoded(target)))
    return rows


def expected_noise():
    rows = []
    for modes in (4, 8, 16):
        psi = np.exp(1j*np.arange(modes)/3)/math.sqrt(modes)
        rho = np.outer(psi, psi.conj())
        group = np.equal.outer(np.arange(modes)//2, np.arange(modes)//2)
        for exposure in (.03, .4):
            s = math.exp(-exposure)
            out = s*s*rho+(s-s*s)*rho*group+(1-s)*np.diag(np.diag(rho))
            rows.append(dict(modes=modes, exposure=exposure, output=encoded(out),
                half_trace_distance=float(sum(abs(np.linalg.eigvalsh(out-rho)))/2),
                population_free_upper=-math.expm1(-2*exposure),
                vacuum_bitflip_failure=-math.expm1(modes*math.log1p(-.02))))
    return rows


def expected_resources():
    clock_interval_certificate()  # Parent's analytic signal, not a copied success flag.
    with mp.workdps(65):
        a, c, h, rate = map(mp.mpf, ('1e-9', '3', '.01', '1e-14'))
        tau = a/mp.sqrt(3)
        tick = (1+h)*tau
        ticks = int(mp.ceil(25*mp.pi/tau))
        side = 1
        while side*a < 2*(90+a):
            side *= 2
        depth = side.bit_length()-1
        omega, event = 96*mp.pi/(h*tau), h*tau/256
        layer = mp.pi/omega+2*event
        preparation = a*mp.sqrt(3)*(side-1)/6+(7*depth+31)*layer+(2*depth+2)*event
        read = 16*layer+event
        run = ticks*tick
        noise = 2*rate*(preparation+read+run)
        size = 1
        while size*a < 2048:
            size *= 2
        buffers = 32*size**3+16*(side**3-1)//7
        pulses, flights = 32*side**3-1, 8*(side**3-1)//7
        total = 2+buffers+pulses+flights+ticks*size**3
        bits = 0
        while 2**bits < total:
            bits += 1
        records = (1+8*bits)*(4*buffers+4*pulses+4*flights+(256*ticks+100)*size**3)
        result = dict(a=a, c=c, overhead_fraction=h, phase_rate=rate, flight_tick=tau, wall_tick=tick,
            ticks=ticks, preparation_side_cells=side, preparation_leaves=side**3, preparation_depth=depth,
            preparation_pulses=32*side**3-1, preparation_flights=8*(side**3-1)//7,
            workspace_side_cells=size, mode_buffers=32*size**3+16*(side**3-1)//7,
            active_processors=size**3+(side**3-1)//7,
            record_identifier_bits=bits, central_record_slots_upper=records,
            pulses_per_cell_tick=48, code_events_per_cell_tick=128, register_flights_per_cell_tick=16,
            drive_norm_bound=omega, code_event_time=event, preparation_time=preparation, run_time=run,
            read_time=read, matter_speed=1/(1+h), clock_velocity=mp.mpf('.6')/(1+h),
            beat_frequency=mp.mpf('.08')/(1+h), probability_noise_upper=noise,
            accounting_energy_error_upper=mp.pi/tau*noise, observable_swing_lower=mp.mpf('.6753')-2*noise,
            reporting_delay_upper=200*mp.sqrt(3))
        need(ticks*tau < 25*mp.pi+tau, 'parent comparison horizon')
        need(noise < mp.mpf('4e-12') and result['accounting_energy_error_upper'] < mp.mpf('.022'), 'finite budget')
        need(preparation < 105 and run < 80 and result['observable_swing_lower'] > mp.mpf('.6752'), 'useful full cycle')
        need(a*size/2 > 600+mp.mpf('.6')*run and a*size/2 > 90+c*run, 'finite domain contains detector and native cone')
        return {k: float(v) if isinstance(v, mp.mpf) else v for k, v in result.items()}


def verify_topology(rows):
    need(type(rows) is list and len(rows) == 4, 'support catalogue')
    check_tower([row['support'] for row in rows])
    old = None
    for level, row in enumerate(rows):
        need(type(row) is dict and row.keys() == {'support', 'fibres', 'carriers', 'coarsen', 'section', 'process_links'}, 'fibre schema')
        count = 10*4**level+2
        fibres = [2+level]+[1]*(count-1)
        carriers = [[v, j] for v in range(count) for j in range(fibres[v])]
        same(row['fibres'], fibres)
        same(row['carriers'], carriers)
        same(row['section'], [carriers.index([v, 0]) for v in range(count)])
        same(row['process_links'], [[i, j] for i, j in itertools.combinations(range(len(carriers)), 2)
                                   if carriers[i][0] == carriers[j][0]])
        if old is None:
            same(row['coarsen'], None)
        else:
            target = []
            for v, j in carriers:
                label = [v, j] if [v, j] in old['carriers'] else [row['support']['coarsen'][v], 0]
                target.append(old['carriers'].index(label))
            same(row['coarsen'], target)
            for i, parent in enumerate(target):
                need(old['carriers'][parent][0] == row['support']['coarsen'][carriers[i][0]], 'commuting bridge')
            for i, j in old['process_links']:
                a, b = old['carriers'][i], old['carriers'][j]
                need([carriers.index(a), carriers.index(b)] in row['process_links'], 'retained private link endpoints')
        old = row


def verify_evidence(packet):
    need(type(packet) is dict and packet.keys() == {'controls', 'channels', 'instruments', 'entangler', 'coin',
                                                  'spatial', 'preparations', 'noise', 'resources', 'topology'}, 'evidence schema')
    verify_controls(packet['controls'])
    same(packet['channels'], expected_channels(), 'channels')
    same(packet['instruments'], expected_instruments(), 'instruments')
    same(packet['entangler'], expected_entangler(), 'entangler')
    verify_coin(packet['coin'])
    same(packet['spatial'], expected_spatial(), 'spatial')
    same(packet['preparations'], expected_preparations(), 'preparations')
    same(packet['noise'], expected_noise(), 'noise')
    same(packet['resources'], expected_resources(), 'resources')
    verify_topology(packet['topology'])
