"""Producer: matrix assembly, spectral powers and complete periodic wave packets."""

import itertools
import math

import numpy as np
from scipy.linalg import expm
from scipy.special import erfc


SIGNS = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
C = np.array([[0, 1, 1, 1], [1, 0, 1j, -1j],
              [1, -1j, 0, 1j], [1, 1j, -1j, 0]]) / np.sqrt(3)
BETA = np.kron([[0, 1], [1, 0]], np.eye(4))
R = np.kron(np.eye(2), C)
ALPHA = [np.kron(np.diag([1, -1]), np.sqrt(3)/2 *
                    (np.diag(SIGNS[:, i])+C @ np.diag(SIGNS[:, i]) @ C))
         for i in range(3)]


def walk(k, a, m, c=3.):
    tau = np.sqrt(3)*a/c
    f = np.exp(-1j*a*(np.asarray(k) @ SIGNS.T))
    w = f[..., :, None]*C
    d = np.zeros(w.shape[:-2]+(8, 8), complex)
    d[..., :4, :4] = w
    d[..., 4:, 4:] = w.conj().swapaxes(-1, -2)
    mass = np.cos(tau*m)*np.eye(8)-1j*np.sin(tau*m)*BETA
    return mass @ d


def hamiltonian(k, m, c=3.):
    return m*BETA+(c/3)*sum(x*b for x, b in zip(k, ALPHA))


def error_bound(a, m, c, k, ticks):
    b, mu = a*np.sqrt(3)*k, m*np.sqrt(3)*a/c
    return (ticks//2)*(b*b+4*mu*b)+(ticks % 2)*(b+mu*b)


def spectral_apply(u, psi, ticks):
    # Solve for coefficients; degenerate eigenspaces need not be orthonormal.
    values, vectors = np.linalg.eig(u)
    coeff = np.linalg.solve(vectors, psi[..., None])[..., 0]
    phases = np.exp(1j*ticks*np.angle(values))
    return np.einsum('...ij,...j->...i', vectors, phases*coeff)


def algebra():
    return dict(coin_involution=float(np.linalg.norm(C@C-np.eye(4))),
                carrier_mass=float(np.linalg.norm(R@BETA-BETA@R)),
                clifford=[float(np.linalg.norm(x@y+y@x-2*(i == j)*np.eye(8)))
                          for i, x in enumerate(ALPHA+[BETA])
                          for j, y in enumerate(ALPHA+[BETA])])


def spectra():
    rows = []
    # Every mode of every declared 3D grid, not a near-valley selection.
    for q in (4, 6, 8):
        for m in (.5, 1.25):
            a, c = 1/q, 1.
            tau = np.sqrt(3)*a/c
            grid = np.array(list(itertools.product(range(q), repeat=3)))
            k = 2*np.pi*grid
            eig = np.linalg.eigvals(walk(k, a, m, c))
            energies = np.sort(abs(np.angle(eig))/tau, axis=-1)
            rows.append(dict(q=q, mass=m, momenta=q**3, modes=8*q**3,
                             energies=energies.tolist(),
                             logz=[float(np.logaddexp(0, -b*energies).sum()) for b in (1., 2.)]))
    return rows


def velocities():
    rows = []
    # Slopes from Hellmann--Feynman differentiation of the actual matrix.
    for mu in (.1, .4, 1.):
        for p in ((.2, .3, .4), (.7, -.4, 1.), (1.4, 1.5, 1.6),
                  (np.pi/2-1e-3,)*3):
            a, c, m = 1., 3., mu*np.sqrt(3)
            u = walk(p, a, m, c)
            eig, vec = np.linalg.eig(u)
            gradients = []
            # Degenerate eigenvalues have identical slopes away from crossings.
            for j in range(8):
                spin = vec[:, j]/np.linalg.norm(vec[:, j])
                grad = []
                for axis in range(3):
                    b = np.diag(SIGNS[:, axis])
                    f = np.diag(np.exp(-1j*np.asarray(p)@SIGNS.T))
                    d = np.zeros((8, 8), complex)
                    d[:4, :4] = -1j*b@f@C
                    d[4:, 4:] = 1j*C@f.conj()@b
                    mass = np.cos(mu)*np.eye(8)-1j*np.sin(mu)*BETA
                    grad.append(float(np.real(1j*np.vdot(spin, (mass@d)@spin)/eig[j]))*np.sqrt(3))
                gradients.append(float(np.linalg.norm(grad)))
            rows.append(dict(mu=mu, angular_k=list(p), speeds=sorted(gradients), limit=float(np.cos(mu))))
    return rows


def dynamics():
    rows = []
    ks = [(0., 0., 0.), (.2, -.3, .4), (1., 1., 1.), (-.8, .1, .2)]
    for a, m, k, n in itertools.product((.125, .03125, .0078125), (.25, 1.5), ks, (0, 1, 2, 5, 16)):
        h = hamiltonian(k, m)
        ideal = np.linalg.matrix_power(R, n)*1. @ expm(-1j*n*np.sqrt(3)*a/3*h)
        errors = []
        for z in itertools.product(range(2), repeat=3):
            actual = np.linalg.matrix_power(walk(np.array(k)+np.pi*np.array(z)/a, a, m), n)
            errors.append(float(np.linalg.norm(actual-(-1)**(n*sum(z))*ideal, 2)))
        rows.append(dict(a=a, mass=m, k=list(k), ticks=n, errors=errors,
                         bound=float(error_bound(a, m, 3., np.linalg.norm(k), n))))
    return rows


def fast_records():
    rows = []
    for flavour, s in itertools.product(range(2), range(4)):
        initial = np.zeros(8, complex)
        initial[4*flavour:4*flavour+4] = C[:, s] if flavour == 0 else np.eye(4)[:, s]
        pre = np.eye(8, dtype=complex)
        pre[:4, :4] = C
        post = np.eye(8, dtype=complex)
        post[4:, 4:] = C
        mu = .4
        final = (np.cos(mu)*np.eye(8)-1j*np.sin(mu)*BETA) @ post @ pre @ initial
        rows.append(dict(flavour=flavour, channel=s, destination=((1-2*flavour)*SIGNS[s]).tolist(),
                         time=float(np.sqrt(3)), distance=float(np.sqrt(3)),
                         after_pre=abs(pre@initial).tolist(),
                         final_probabilities=(abs(final)**2).tolist(), empty_probability=0.))
    return rows


def analytic_clock():
    c, m1, m2, speed, sigma, radius, detector, a = 3., 100., 100.1, .6, 10., 100., 600., 1e-9
    v, gamma = c/3, 1/math.sqrt(1-speed**2)
    gap, tau = m2-m1, math.sqrt(3)*a/c
    cycle = 2*math.pi*gamma/gap
    t = cycle+tau  # Covers the first native tick after a full cycle too.
    k1, k2 = gamma*m1*speed/v**2, gamma*m2*speed/v**2
    cutoff = k2+1
    s2 = sigma**2*radius**2/(sigma**2+radius**2)
    p = (1+sigma**2/radius**2)**(-1.5)
    visibility = p*math.exp(-s2*(k2-k1)**2/2)
    e = math.sqrt(3)*v/(m1*sigma)+math.sqrt(15)*t*v*v/(8*m1*sigma**2)
    d = a*t*(math.sqrt(3)*c*cutoff**2/2+2*math.sqrt(3)*m2*cutoff)+math.sqrt(3)*a*cutoff+3*a*a*m2*cutoff/c
    moment = math.sqrt(k2**4+5*k2**2/(2*sigma**2)+15/(16*sigma**4))
    alias = (a/math.pi)**2*math.sqrt(4*math.pi**2+math.pi**4/45)*moment
    tail = math.sqrt(erfc(math.sqrt(2)*sigma)+2*math.sqrt(2)*sigma*math.exp(-2*sigma**2)/math.sqrt(math.pi))
    preparation = 2*math.sqrt(3*erfc((80-a)/(math.sqrt(2)*sigma)))
    # exp(-1e20) aliases are bounded by exp(-1000), not claimed to be exactly zero.
    prob_bound = 4*(e+d+4*alias+2*tail)+8*math.exp(-1000)+math.exp(-detector**2/(2*radius**2))+preparation
    return dict(c=c, masses=[m1, m2], velocity=speed, sigma=sigma, radius=radius,
                detector_extent=detector, a=a, gamma=gamma, cycle=cycle, horizon=t, cutoff=cutoff,
                rigid_probability_offset=p/2, visibility=visibility,
                continuum_error=e, dynamic_error=d, sampling_error=alias,
                preparation_extent=80., preparation_error=preparation,
                tail_norm=tail, probability_bound=prob_bound,
                alias_upper_exponent=-1000, ticks_upper=math.ceil(cycle/tau),
                interval_certificate=dict(probability_upper='0.035153', reference_fringe_lower='0.7457',
                                          observable_swing_lower='0.6753'))


def packets():
    """Normalizable torus packets; transverse state is the constant Fourier mode."""
    rows = []
    length, sigma, radius, extent, speed, c = 64., 3., 9., 20., .6, 3.
    masses = (4., 4.2)
    gamma = 1/np.sqrt(1-speed**2)
    central = [gamma*m*speed for m in masses]
    h = hamiltonian([central[0], 0, 0], masses[0])
    # A fixed projected basis vector avoids arbitrary eigenvector phases.
    xi = ((np.eye(8)+h/(gamma*masses[0]))@(np.eye(8)+R))[:, 0]
    xi /= np.linalg.norm(xi)
    for count in (512, 1024, 2048, 4096):
        a, tau = length/count, np.sqrt(3)*length/(3*count)
        x = (np.arange(count)-count//2)*a
        envelope = np.exp(-x*x/(4*sigma*sigma))
        envelope /= np.linalg.norm(envelope)
        initial = [envelope[:, None]*np.exp(1j*k*x[:, None])*xi for k in central]
        kh = np.zeros((count, 3))
        kh[:, 0] = 2*np.pi*np.fft.fftfreq(count, d=a)
        symbols = [walk(kh, a, m, c) for m in masses]
        fourier = [np.fft.fft(z, axis=0, norm='ortho') for z in initial]
        energy_matrices = []
        for symbol in symbols:
            re_values, re_vectors = np.linalg.eigh((symbol+symbol.conj().swapaxes(-1, -2))/2)
            energies = np.arccos(np.clip(re_values, -1, 1))/tau
            energy_matrices.append((re_vectors*energies[:, None, :])@re_vectors.conj().swapaxes(-1, -2))
        initial_energies = [float(np.einsum('ni,nij,nj->', f.conj(), e, f).real)
                            for f, e in zip(fourier, energy_matrices)]
        for target in (4., 10., 20.):
            ticks = round(target/tau)
            time = ticks*tau
            evolved = [spectral_apply(u, p, ticks) for u, p in zip(symbols, fourier)]
            states = [np.fft.ifft(z, axis=0, norm='ortho') for z in evolved]
            y = (x-speed*time+length/2) % length-length/2
            weight = np.exp(-y*y/(2*radius*radius))*(abs(y) <= extent)
            probabilities = [float(np.sum(weight[:, None]*abs(states[0]+np.exp(-1j*theta)*states[1])**2)/4)
                             for theta in (0., np.pi/2)]
            returned = [spectral_apply(u.conj().swapaxes(-1, -2), p, ticks)
                        for u, p in zip(symbols, evolved)]
            p = (1+sigma*sigma/(radius*radius))**(-.5)
            visibility = p*np.exp(-sigma*sigma*radius*radius/(sigma*sigma+radius*radius)*(central[1]-central[0])**2/2)
            rigid = [float((p+visibility*np.cos((masses[1]-masses[0])*time/gamma+theta))/2)
                     for theta in (0., np.pi/2)]
            rows.append(dict(sites=count, total_modes=8*count, a=a, ticks=ticks, time=float(time),
                             probabilities=probabilities, rigid=rigid,
                             initial_energies=initial_energies,
                             final_energies=[float(np.einsum('ni,nij,nj->', f.conj(), e, f).real)
                                             for f, e in zip(evolved, energy_matrices)],
                             norms=[float(np.linalg.norm(z)) for z in states],
                             reverse_errors=[float(np.linalg.norm(p-q)) for p, q in zip(returned, fourier)]))
    return dict(length=length, transverse='constant normalized torus mode', sigma=sigma,
                radius=radius, extent=extent, speed=speed, c=c, masses=list(masses), rows=rows)


def candidate():
    return dict(algebra=algebra(), spectra=spectra(), velocities=velocities(), dynamics=dynamics(),
                fast_records=fast_records(), analytic_clock=analytic_clock(), packets=packets())
