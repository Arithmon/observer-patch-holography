"""Independent packet propagation and outward finite-clock certificate."""

import math
import mpmath as mp
import numpy as np
from scipy.special import erfc
from m1_operational_clocks.check import frames, gates, binary_apply
from .spectrum_check import sea
from .format import need, keys, exact, close, real_vector


def verify_packets(row):
    keys(row, 'length sigma velocity masses transverse rows')
    exact({k: v for k, v in row.items() if k != 'rows'}, dict(length=64., sigma=3., velocity=.6,
          masses=[4., 4.2], transverse='constant normalized torus mode; finite control only'))
    rows = row['rows']
    need(type(rows) is list and len(rows) == 12, 'all lattice sizes and times')
    _, _, beta, r, alpha = frames()
    xi = ((np.eye(8)+(4*beta+3*alpha[0])/5)@(np.eye(8)+r))[:, 0]
    xi /= np.linalg.norm(xi)
    for index, count in enumerate((256, 512, 1024)):
        a, tau = 64/count, 64/(count*np.sqrt(3))
        x = np.arange(-count//2, count//2)*a
        envelope = np.exp(-(x/6)**2)
        envelope /= np.linalg.norm(envelope)
        us, hs, ps, initial = [], [], [], []
        for mass, center in ((4., 3.), (4.2, 3.15)):
            u = np.array([gates([k, 0., 0.], a, mass, 3.) for k in 2*np.pi*np.fft.fftfreq(count, d=a)])
            parts = [sea(matrix, tau) for matrix in u]
            us.append(u)
            hs.append(np.array([p[0] for p in parts]))
            ps.append(np.array([p[1] for p in parts]))
            initial.append(np.fft.fft(envelope[:, None]*np.exp(1j*center*x[:, None])*xi/np.sqrt(2), axis=0, norm='ortho'))
        projected = [np.einsum('nij,nj->ni', p, z) for p, z in zip(ps, initial)]
        blocked = float(sum(np.vdot(z, z).real for z in projected))
        positive = [(z-pz)/np.sqrt(1-blocked) for z, pz in zip(initial, projected)]
        energy = float(sum(np.einsum('ni,nij,nj->', z.conj(), h, z).real for z, h in zip(positive, hs)))
        need(0 < blocked < .1 and energy > 0, 'genuine sea background and positive added particle')
        for offset, target in enumerate((0., 4., 10., 20.)):
            item = rows[4*index+offset]
            keys(item, 'sites a ticks time blocked background probabilities energy norm rigid')
            ticks = round(target/tau)
            time = ticks*tau
            exact([item['sites'], item['a'], item['ticks'], item['time']], [count, a, ticks, float(time)])
            evolved = [binary_apply(u, z, ticks) for u, z in zip(us, positive)]
            y = (x-.6*time+32) % 64-32
            e = np.exp(-(y/6)**2)
            e /= np.linalg.norm(e)
            probes = [np.fft.fft(e[:, None]*np.exp(1j*k*x[:, None])*xi/np.sqrt(2), axis=0, norm='ortho') for k in (3., 3.15)]
            bg = float(sum(np.einsum('ni,nij,nj->', z.conj(), p, z).real for z, p in zip(probes, ps)))
            overlaps = [np.vdot(f, z) for f, z in zip(probes, evolved)]
            expected = [float(bg+abs(overlaps[0]+np.exp(-1j*t)*overlaps[1])**2) for t in (0., np.pi/2)]
            for key, value in [('blocked', blocked), ('background', bg), ('energy', energy), ('norm', 1.)]:
                exact(item[key], float(value))
            close(real_vector(item['probabilities'], 2), expected, 'complete sea-read probabilities')
            exact(item['rigid'], [float((1+np.cos(.2*time/1.25+t))/2) for t in (0., np.pi/2)])
            need(all(-1e-10 <= p <= 1+1e-10 for p in expected), 'physical bounded binary effect')
            if ticks == 0:
                close(expected[0], 1., 'deterministic initial click')
            # Full norm and energy conservation, not just fitting one fringe.
            close(sum(np.vdot(z, z).real for z in evolved), 1., 'complete packet normalization')
            close(sum(np.einsum('ni,nij,nj->', z.conj(), h, z).real for z, h in zip(evolved, hs)), energy, 'excitation conservation')


def interval_certificate():
    """Analytic upper bounds, evaluated with outward interval arithmetic."""
    iv, previous = mp.iv, mp.iv.dps
    iv.dps = 50
    try:
        a, m, high, sigma, speed = iv.mpf('1e-9'), iv.mpf(100), iv.mpf('100.1'), iv.mpf(10), iv.mpf('.6')
        tau, gamma = a/iv.sqrt(3), 1/iv.sqrt(1-speed**2)
        t = 2*iv.pi*gamma/(high-m)+tau
        center, boundary = gamma*high*speed, iv.mpf('1e-10')
        k = center+1
        e = iv.sqrt(3)/(m*sigma)+iv.sqrt(15)*t/(8*m*sigma**2)
        d = a*t*(iv.sqrt(3)*3*k*k/2+2*iv.sqrt(3)*high*k)+iv.sqrt(3)*a*k+a*a*high*k
        moment = iv.sqrt(center**4+5*center**2/(2*sigma**2)+15/(16*sigma**4))
        alias = (a/iv.pi)**2*iv.sqrt(4*iv.pi**2+iv.pi**4/45)*moment
        x = iv.sqrt(2)*sigma
        tail = iv.sqrt(iv.exp(-x*x)*(1/x+2*x)/iv.sqrt(iv.pi))
        y = (80-a)/(iv.sqrt(2)*sigma)
        cut = 4*iv.sqrt(3*iv.exp(-y*y)/(y*iv.sqrt(iv.pi)))
        need((2*iv.pi**2*sigma**2/a**2).a > iv.mpf(1010), 'uniform translated Gaussian alias domain')
        tiny = 8*iv.exp(-1000)  # never round the analytic alias allowance down to zero
        z = (high*tau+iv.sqrt(3)*a*k)**2/(2*tau)
        b = 3/(4*m*m*sigma*sigma)+iv.pi*z/(4*m)+2*tail+4*alias+4*cut+boundary+tiny
        error = b+iv.sqrt(b)+4*(e+d+4*alias+2*tail)+8*cut+boundary+tiny
        phase = (high-m)*tau/(2*gamma)
        swing = 1-phase**2/4-2*error
        need(error.b < iv.mpf('.038'), 'finite informative error enclosure')
        need(swing.a > iv.mpf('.923'), 'finite native-tick observable swing')
        return dict(probability_upper='0.038', observable_swing_lower='0.923')
    finally:
        iv.dps = previous


def verify_witness(row):
    keys(row, 'a masses sigma velocity c extent boundary_error cutoff horizon probability_error blocking_upper boundary_margin interval')
    exact({k: row[k] for k in ('a', 'masses', 'sigma', 'velocity', 'c', 'extent', 'boundary_error')},
          dict(a=1e-9, masses=[100., 100.1], sigma=10., velocity=.6, c=3., extent=80., boundary_error=1e-10))
    # Ordinary floats are a reproducible diagnostic; the certificate is independently outward.
    a, tau, k = 1e-9, 1e-9/math.sqrt(3), 76.075
    t = 2*math.pi*1.25/(100.1-100.)+tau
    e = math.sqrt(3)/1000+math.sqrt(15)*t/80000
    d = a*t*(math.sqrt(3)*3*k*k/2+2*math.sqrt(3)*100.1*k)+math.sqrt(3)*a*k+a*a*100.1*k
    moment = math.sqrt(75.075**4+5*75.075**2/200+15/160000)
    alias = (a/math.pi)**2*math.sqrt(4*math.pi**2+math.pi**4/45)*moment
    tail = math.sqrt(erfc(math.sqrt(2)*10)+20*math.sqrt(2)*math.exp(-200)/math.sqrt(math.pi))
    cut = 4*math.sqrt(3*erfc((80-a)/(math.sqrt(2)*10)))
    b = 3/4000000+math.pi*(100.1*tau+math.sqrt(3)*a*k)**2/(800*tau)+2*tail+4*alias+4*cut+1e-10
    bound = b+math.sqrt(b)+4*(e+d+4*alias+2*tail)+8*cut+1e-10
    exact([row['cutoff'], row['horizon'], row['probability_error'], row['blocking_upper']], [k, t, bound, b])
    margin = row['boundary_margin']
    need(type(margin) is int and margin > 0, 'integer finite boundary margin')
    # Check the ceiling with high precision; a float log1p is not the certificate.
    with mp.workdps(60):
        delta = mp.sin(mp.mpf(100)*mp.mpf('1e-9')/mp.sqrt(3))
        minimal = mp.log(4/(delta*mp.mpf('1e-10')))/mp.log1p(delta/2)
        need(margin-1 < minimal <= margin, 'certified integer collar')
        need(4/delta*mp.exp(-mp.log1p(delta/2)*margin) <= mp.mpf('1e-10'), 'boundary suppression')
    exact(row['interval'], interval_certificate())


def verify(row):
    keys(row, 'packets witness')
    verify_packets(row['packets'])
    verify_witness(row['witness'])
