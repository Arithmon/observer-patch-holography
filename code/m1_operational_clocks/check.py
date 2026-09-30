"""Independent replay: entrywise gates, closed spectra and binary propagation.

Imports neither the producer nor its matrix/clock helpers. Numerical agreement
supports, but does not replace, the analytic proofs in the accompanying note.
"""

import itertools
import math

import mpmath as mp
import numpy as np


def need(condition, message):
    if not condition:
        raise ValueError(message)


def same(actual, expected, path='evidence'):
    need(type(actual) is type(expected), f'{path}: wrong JSON type')
    if isinstance(expected, dict):
        need(actual.keys() == expected.keys(), f'{path}: missing/extra keys')
        for key in expected:
            same(actual[key], expected[key], f'{path}.{key}')
    elif isinstance(expected, list):
        need(len(actual) == len(expected), f'{path}: wrong catalog length')
        for i, (a, b) in enumerate(zip(actual, expected)):
            same(a, b, f'{path}[{i}]')
    elif isinstance(expected, float):
        need(math.isfinite(actual) and math.isfinite(expected), f'{path}: nonfinite value')
        residual = any(tag in path for tag in ('coin_involution', 'carrier_mass', '.clifford[',
                       '.errors[', '.reverse_errors[', '.after_pre[', '.final_probabilities['))
        tolerance = 2e-10 if residual else (2e-9*abs(expected) if expected else 1e-12)
        need(abs(actual-expected) <= tolerance, f'{path}: independent replay differs ({actual}, {expected})')
    else:
        need(actual == expected, f'{path}: independent replay differs')


def frames():
    points = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    coin = np.zeros((4, 4), complex)
    for i in range(1, 4):
        coin[0, i] = coin[i, 0] = 1/math.sqrt(3)
    for i, j in ((1, 2), (2, 3), (3, 1)):
        coin[i, j], coin[j, i] = 1j/math.sqrt(3), -1j/math.sqrt(3)
    beta = np.zeros((8, 8), complex)
    carrier = np.zeros((8, 8), complex)
    velocity = [np.zeros((8, 8), complex) for _ in range(3)]
    for p in range(8):
        beta[p, (p+4) % 8] = 1
        for q in range(8):
            if p//4 == q//4:
                carrier[p, q] = coin[p % 4, q % 4]
                for axis in range(3):
                    x = (points[p % 4][axis] if p == q else 0)
                    x += sum(coin[p % 4, j]*points[j][axis]*coin[j, q % 4] for j in range(4))
                    velocity[axis][p, q] = (1 if p < 4 else -1)*math.sqrt(3)*x/2
    return points, coin, beta, carrier, velocity


def gates(k, a, mass, c):
    points, coin, beta, _, _ = frames()
    k = np.asarray(k)
    result = np.zeros(k.shape[:-1]+(8, 8), complex)
    mu = mass*math.sqrt(3)*a/c
    for row in range(8):
        for col in range(8):
            spinrow, spincol = row % 4, col % 4
            f = col//4
            if f == 0:
                flight_coin = np.exp(-1j*a*sum(k[..., d]*points[spinrow][d] for d in range(3)))*coin[spinrow, spincol]
            else:
                flight_coin = coin[spinrow, spincol]*np.exp(1j*a*sum(k[..., d]*points[spincol][d] for d in range(3)))
            result[..., row, col] = (math.cos(mu) if row//4 == f else -1j*math.sin(mu))*flight_coin
    return result


def binary_apply(matrix, vector, ticks):
    result = vector.copy()
    factor = matrix.copy()
    n = ticks
    while n:
        if n & 1:
            result = np.einsum('...ij,...j->...i', factor, result)
        n //= 2
        if n:
            factor = factor @ factor
    return result


def replay_algebra():
    _, coin, beta, carrier, alpha = frames()
    return dict(coin_involution=float(np.linalg.norm(coin@coin-np.eye(4))),
                carrier_mass=float(np.linalg.norm(carrier@beta-beta@carrier)),
                clifford=[float(np.linalg.norm(x@y+y@x-2*(i == j)*np.eye(8)))
                          for i, x in enumerate(alpha+[beta]) for j, y in enumerate(alpha+[beta])])


def replay_spectra():
    rows = []
    for q in (4, 6, 8):
        for mass in (.5, 1.25):
            tau = math.sqrt(3)/q
            energies = []
            for k in itertools.product(range(q), repeat=3):
                ss = sum(math.sin(2*math.pi*j/q)**2 for j in k)/3
                omega = math.acos(math.cos(tau*mass)*math.sqrt(max(0, 1-ss)))
                energies.append([omega/tau]*4+[(math.pi-omega)/tau]*4)
            need(len(energies) == q**3, 'complete Brillouin grid')
            need(min(min(row) for row in energies) >= mass-1e-12, 'mass gap')
            rows.append(dict(q=q, mass=mass, momenta=q**3, modes=8*q**3, energies=energies,
                             logz=[math.fsum(math.log1p(math.exp(-b*e)) for row in energies for e in row)
                                   for b in (1., 2.)]))
    return rows


def replay_dynamics():
    _, _, beta, carrier, alpha = frames()
    rows = []
    for a, mass, k, n in itertools.product((1/8, 1/32, 1/128), (1/4, 3/2),
                                         ((0., 0., 0.), (.2, -.3, .4), (1., 1., 1.), (-.8, .1, .2)),
                                         (0, 1, 2, 5, 16)):
        normk = math.sqrt(sum(z*z for z in k))
        energy = math.sqrt(mass**2+normk**2)
        h = mass*beta+sum(z*x for z, x in zip(k, alpha))
        t = n*math.sqrt(3)*a/3
        ideal = math.cos(t*energy)*np.eye(8)-1j*math.sin(t*energy)*h/energy
        if n % 2:
            ideal = carrier @ ideal
        errors = []
        for z in itertools.product(range(2), repeat=3):
            u = gates(np.array(k)+np.pi*np.array(z)/a, a, mass, 3.)
            actual = np.eye(8, dtype=complex)
            for _ in range(n):
                actual = u @ actual
            errors.append(float(np.linalg.norm(actual-(-1)**(n*sum(z))*ideal, 2)))
        b, mu = a*math.sqrt(3)*normk, mass*math.sqrt(3)*a/3
        bound = (n//2)*(b**2+4*mu*b)+(n % 2)*(b+mu*b)
        need(max(errors) <= bound+2e-12, 'full-channel analytic bound violated')
        rows.append(dict(a=a, mass=mass, k=list(k), ticks=n, errors=errors, bound=bound))
    return rows


def replay_velocities():
    rows = []
    for mu in (.1, .4, 1.):
        for p in ((.2, .3, .4), (.7, -.4, 1.), (1.4, 1.5, 1.6), (math.pi/2-1e-3,)*3):
            s = sum(math.sin(x)**2 for x in p)/3
            # |grad omega|^2 = cos(mu)^2 sum sin(2p)^2 /
            #                    [36 (1-s) (1-cos(mu)^2(1-s))].
            speed = math.sqrt(3*math.cos(mu)**2*sum(math.sin(2*x)**2 for x in p)/
                              (36*(1-s)*(1-math.cos(mu)**2*(1-s))))
            need(speed <= math.cos(mu)+2e-10, 'sharp all-band velocity bound')
            rows.append(dict(mu=mu, angular_k=list(p), speeds=[speed]*8, limit=math.cos(mu)))
    return rows


def replay_fast_records():
    points, coin, _, _, _ = frames()
    rows = []
    for f, channel in itertools.product(range(2), range(4)):
        after = [0.]*8
        after[4*f+channel] = 1.
        # Flight routes exactly this occupied mode; all subsequent mixing is
        # at its destination. Independent scalar probabilities suffice here.
        spin = np.eye(4)[:, channel] if f == 0 else coin[:, channel]
        prob = [(math.cos(.4)**2 if j//4 == f else math.sin(.4)**2)*float(abs(spin[j % 4])**2)
                for j in range(8)]
        need(abs(sum(prob)-1) < 1e-14, 'lossless native record')
        rows.append(dict(flavour=f, channel=channel, destination=[(1-2*f)*z for z in points[channel]],
                         time=math.sqrt(3), distance=math.sqrt(3), after_pre=after,
                         final_probabilities=prob, empty_probability=0.))
    return rows


def replay_analytic_clock():
    with mp.workdps(70):
        c, m1, m2, speed = map(mp.mpf, ('3', '100', '100.1', '.6'))
        sigma, radius, extent, a = map(mp.mpf, ('10', '100', '600', '1e-9'))
        gamma = 1/mp.sqrt(1-speed**2)
        gap = m2-m1
        cycle = 2*mp.pi*gamma/gap
        tau = mp.sqrt(3)*a/c
        t = cycle+tau
        k1, k2 = gamma*m1*speed, gamma*m2*speed
        cutoff = k2+1
        p = (1+(sigma/radius)**2)**(-mp.mpf('1.5'))
        s2 = 1/(1/sigma**2+1/radius**2)
        visibility = p*mp.exp(-s2*(k2-k1)**2/2)
        e = mp.sqrt(3)/(m1*sigma)+mp.sqrt(15)*t/(8*m1*sigma**2)
        d = a*t*(mp.sqrt(3)*c*cutoff**2/2+2*mp.sqrt(3)*m2*cutoff)+mp.sqrt(3)*a*cutoff+3*a*a*m2*cutoff/c
        moment = mp.sqrt(k2**4+5*k2**2/(2*sigma**2)+15/(16*sigma**4))
        alias = (a/mp.pi)**2*mp.sqrt(4*mp.pi**2+mp.pi**4/45)*moment
        tail = mp.sqrt(mp.erfc(mp.sqrt(2)*sigma)+2*mp.sqrt(2)*sigma*mp.exp(-2*sigma**2)/mp.sqrt(mp.pi))
        preparation = 2*mp.sqrt(3*mp.erfc((80-a)/(mp.sqrt(2)*sigma)))
        for z in (2*mp.pi**2*sigma**2/a**2, mp.pi**2*s2/(2*a**2)):
            # B(z)=6w+12w^2+8w^3, w=e^-z/(1-e^-3z).
            need(z > 1010, 'analytic alias estimate requires this lower bound')
        bound = 4*(e+d+4*alias+2*tail)+8*mp.exp(-1000)+mp.exp(-extent**2/(2*radius**2))+preparation
        need(bound < mp.mpf('.08') and visibility > mp.mpf('.74'), 'clock witness is not informative')
        need(a*(k2-k1) < mp.pi and cutoff < mp.pi/(2*a), 'clock domain')
        result = dict(c=c, masses=[float(m1), float(m2)], velocity=speed, sigma=sigma, radius=radius,
                      detector_extent=extent, a=a, gamma=gamma, cycle=cycle, horizon=t, cutoff=cutoff,
                      rigid_probability_offset=p/2, visibility=visibility,
                      continuum_error=e, dynamic_error=d, sampling_error=alias,
                      preparation_extent=80., preparation_error=preparation,
                      tail_norm=tail, probability_bound=bound, alias_upper_exponent=-1000)
        result = {k: float(v) if isinstance(v, mp.mpf) else v for k, v in result.items()}
        result['ticks_upper'] = int(mp.ceil(cycle/tau))
        result['interval_certificate'] = clock_interval_certificate()
        return result


def clock_interval_certificate():
    """Outward interval enclosure of the finite witness, not a float fit."""
    iv = mp.iv
    previous = iv.dps
    iv.dps = 50
    try:
        a, sigma, radius = iv.mpf('1e-9'), iv.mpf(10), iv.mpf(100)
        m1, m2, speed = iv.mpf(100), iv.mpf('100.1'), iv.mpf('.6')
        gamma = 1/iv.sqrt(1-speed**2)
        k1, k2 = gamma*m1*speed, gamma*m2*speed
        t = 2*iv.pi*gamma/(m2-m1)+a/iv.sqrt(3)
        cutoff = k2+1
        s2 = sigma**2*radius**2/(sigma**2+radius**2)
        e = iv.sqrt(3)/(m1*sigma)+iv.sqrt(15)*t/(8*m1*sigma**2)
        d = a*t*(iv.sqrt(3)*3*cutoff**2/2+2*iv.sqrt(3)*m2*cutoff)+iv.sqrt(3)*a*cutoff+a*a*m2*cutoff
        moment = iv.sqrt(k2**4+5*k2**2/(2*sigma**2)+15/(16*sigma**4))
        alias = (a/iv.pi)**2*iv.sqrt(4*iv.pi**2+iv.pi**4/45)*moment
        x = iv.sqrt(2)*sigma
        # Integral tail inequality erfc(x) <= exp(-x^2)/(x sqrt(pi)).
        tail = iv.sqrt(iv.exp(-x*x)/iv.sqrt(iv.pi)*(1/x+2*x))
        need((2*iv.pi**2*sigma**2/a**2).a > iv.mpf(1010), 'normalization alias domain')
        need((iv.pi**2*s2/(2*a**2)).a > iv.mpf(1010), 'readout alias domain')
        y = (80-a)/(iv.sqrt(2)*sigma)
        preparation = 2*iv.sqrt(3*iv.exp(-y*y)/(y*iv.sqrt(iv.pi)))
        bound = 4*(e+d+4*alias+2*tail)+8*iv.exp(-1000)+iv.exp(-18)+preparation
        visibility = (1+sigma**2/radius**2)**iv.mpf('-1.5')*iv.exp(-s2*(k2-k1)**2/2)
        need(bound.b < iv.mpf('0.035153'), 'outward clock-error enclosure')
        need(visibility.a > iv.mpf('0.7457'), 'outward reference fringe enclosure')
        phase_miss = (m2-m1)*a/(2*iv.sqrt(3)*gamma)
        swing = visibility*(1-phase_miss**2/4)-2*bound
        need(swing.a > iv.mpf('0.6753'), 'observable finite-clock contrast')
        return dict(probability_upper='0.035153', reference_fringe_lower='0.7457', observable_swing_lower='0.6753')
    finally:
        iv.dps = previous


def replay_packets():
    _, _, beta, carrier, alpha = frames()
    size, width, radius, extent, speed, c = 64., 3., 9., 20., .6, 3.
    masses = [4., 4.2]
    gamma = 1.25
    momenta = [gamma*m*speed for m in masses]
    h = masses[0]*beta+momenta[0]*alpha[0]
    spin = ((np.eye(8)+h/(gamma*masses[0]))@(np.eye(8)+carrier))[:, 0]
    spin /= np.linalg.norm(spin)
    rows = []
    for count in (512, 1024, 2048, 4096):
        a = size/count
        tau = a/math.sqrt(3)
        positions = np.array([(j-count//2)*a for j in range(count)])
        envelope = np.exp(-positions**2/(4*width**2))
        envelope /= math.sqrt(float(np.vdot(envelope, envelope).real))
        prepared = [envelope[:, None]*np.exp(1j*p*positions[:, None])*spin for p in momenta]
        frequencies = np.array([2*math.pi*(j if j < count//2 else j-count)/size for j in range(count)])
        vectors = np.stack((frequencies, np.zeros(count), np.zeros(count)), axis=-1)
        operators = [gates(vectors, a, m, c) for m in masses]
        initial = [np.fft.fft(p, axis=0, norm='ortho') for p in prepared]
        energy_matrices = []
        for m, op in zip(masses, operators):
            cosomega = math.cos(m*tau)*np.sqrt(1-np.sin(a*frequencies)**2/3)
            omega = np.arccos(cosomega)
            realpart = (op+op.conj().swapaxes(-1, -2))/2
            energy = math.pi/(2*tau)*np.eye(8)+(omega-math.pi/2)[:, None, None]/(tau*cosomega[:, None, None])*realpart
            energy_matrices.append(energy)
        initial_energies = [float(np.vdot(f, np.einsum('nij,nj->ni', e, f)).real)
                            for f, e in zip(initial, energy_matrices)]
        for m, f, energy in zip(masses, initial, initial_energies):
            weights = np.sum(abs(f)**2, axis=1)
            upper = m+math.pi*c/6*np.dot(weights, abs(frequencies))+4*math.sqrt(3)*math.pi*c*a*np.dot(weights, frequencies**2)
            need(energy <= upper+1e-10, 'full-spectrum energy bound')
        for target in (4., 10., 20.):
            n = round(target/tau)
            t = n*tau
            final = [binary_apply(u, f, n) for u, f in zip(operators, initial)]
            states = [np.fft.ifft(f, axis=0, norm='ortho') for f in final]
            norms = [float(np.linalg.norm(z)) for z in states]
            need(max(abs(z-1) for z in norms) < 2e-10, 'unitary packet evolution')
            displacement = (positions-speed*t+size/2) % size-size/2
            w = np.exp(-displacement**2/(2*radius**2))*(abs(displacement) <= extent)
            local_sum = (abs(states[0])**2+abs(states[1])**2).sum(axis=1)
            cross = (states[0].conj()*states[1]).sum(axis=1)
            probs = [float(np.dot(w, local_sum+2*np.real(np.exp(-1j*theta)*cross))/4)
                     for theta in (0., math.pi/2)]
            back = [binary_apply(u.conj().swapaxes(-1, -2), f, n) for u, f in zip(operators, final)]
            inverse_errors = [float(np.linalg.norm(p-q)) for p, q in zip(back, initial)]
            need(max(inverse_errors) < 2e-10, 'executed inverse failed')
            p = (1+width**2/radius**2)**(-.5)
            vis = p*math.exp(-width**2*radius**2/(width**2+radius**2)*(momenta[1]-momenta[0])**2/2)
            rigid = [float((p+vis*math.cos((masses[1]-masses[0])*t/gamma+theta))/2)
                     for theta in (0., math.pi/2)]
            rows.append(dict(sites=count, total_modes=count*8, a=a, ticks=n, time=t,
                             probabilities=probs, rigid=rigid, norms=norms, reverse_errors=inverse_errors,
                             initial_energies=initial_energies,
                             final_energies=[float(np.vdot(f, np.einsum('nij,nj->ni', e, f)).real)
                                             for f, e in zip(final, energy_matrices)]))
    return dict(length=size, transverse='constant normalized torus mode', sigma=width, radius=radius,
                extent=extent, speed=speed, c=c, masses=masses, rows=rows)


def reconstruct():
    return dict(algebra=replay_algebra(), spectra=replay_spectra(), velocities=replay_velocities(), dynamics=replay_dynamics(),
                fast_records=replay_fast_records(), analytic_clock=replay_analytic_clock(), packets=replay_packets())


def verify_evidence(candidate):
    same(candidate, reconstruct())
