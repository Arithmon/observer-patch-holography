"""Finite sea-clock controls, retaining every momentum and background mode."""

import math
import numpy as np
from scipy.special import erfc
from m1_operational_clocks import model as parent
from .spectrum import sea


def packets():
    rows = []
    length, sigma, speed = 64., 3., .6
    masses = (4., 4.2)
    gamma = 1/math.sqrt(1-speed**2)
    central = [gamma*m*speed for m in masses]
    h = parent.hamiltonian([central[0], 0, 0], masses[0])
    xi = ((np.eye(8)+h/(gamma*masses[0]))@(np.eye(8)+parent.R))[:, 0]
    xi /= np.linalg.norm(xi)
    for count in (256, 512, 1024):
        a, tau = length/count, length/(count*np.sqrt(3))
        x = (np.arange(count)-count//2)*a
        e = np.exp(-x*x/(4*sigma*sigma))
        e /= np.linalg.norm(e)
        initial = [e[:, None]*np.exp(1j*k*x[:, None])*xi/np.sqrt(2) for k in central]
        phi = [np.fft.fft(z, axis=0, norm='ortho') for z in initial]
        ks = np.zeros((count, 3))
        ks[:, 0] = 2*np.pi*np.fft.fftfreq(count, d=a)
        us = [parent.walk(ks, a, m) for m in masses]
        hs, ps = [], []
        for u in us:
            parts = [sea(matrix, tau) for matrix in u]
            hs.append(np.array([p[0] for p in parts]))
            ps.append(np.array([p[1] for p in parts]))
        blocked = float(sum(np.einsum('ni,nij,nj->', z.conj(), p, z).real for z, p in zip(phi, ps)))
        chi = [(z-np.einsum('nij,nj->ni', p, z))/np.sqrt(1-blocked) for p, z in zip(ps, phi)]
        energy = float(sum(np.einsum('ni,nij,nj->', z.conj(), h, z).real for z, h in zip(chi, hs)))
        for target in (0., 4., 10., 20.):
            ticks = round(target/tau)
            time = ticks*tau
            evolved = [parent.spectral_apply(u, z, ticks) for u, z in zip(us, chi)]
            y = (x-speed*time+length/2) % length-length/2
            env = np.exp(-y*y/(4*sigma*sigma))
            env /= np.linalg.norm(env)
            probes = [np.fft.fft(env[:, None]*np.exp(1j*k*x[:, None])*xi/np.sqrt(2), axis=0, norm='ortho') for k in central]
            background = float(sum(np.einsum('ni,nij,nj->', z.conj(), p, z).real for z, p in zip(probes, ps)))
            overlaps = [np.vdot(f, z) for f, z in zip(probes, evolved)]
            probs = [float(background+abs(overlaps[0]+np.exp(-1j*theta)*overlaps[1])**2) for theta in (0., np.pi/2)]
            rows.append(dict(sites=count, a=a, ticks=ticks, time=float(time), blocked=blocked,
                             background=background, probabilities=probs, energy=energy,
                             norm=float(sum(np.vdot(z, z).real for z in evolved)),
                             rigid=[float((1+np.cos((masses[1]-masses[0])*time/gamma+theta))/2) for theta in (0., np.pi/2)]))
    return dict(length=length, sigma=sigma, velocity=speed, masses=list(masses),
                transverse='constant normalized torus mode; finite control only', rows=rows)


def witness():
    a, m1, m2, sigma, speed = 1e-9, 100., 100.1, 10., .6
    gamma, tau = 1/math.sqrt(1-speed**2), a/math.sqrt(3)
    t = 2*math.pi*gamma/(m2-m1)+tau
    k = gamma*m2*speed+1
    e = math.sqrt(3)/(m1*sigma)+math.sqrt(15)*t/(8*m1*sigma*sigma)
    d = a*t*(math.sqrt(3)*3*k*k/2+2*math.sqrt(3)*m2*k)+math.sqrt(3)*a*k+a*a*m2*k
    center = gamma*m2*speed
    moment = math.sqrt(center**4+5*center**2/(2*sigma*sigma)+15/(16*sigma**4))
    alias = (a/math.pi)**2*math.sqrt(4*math.pi**2+math.pi**4/45)*moment
    tail = math.sqrt(erfc(math.sqrt(2)*sigma)+2*math.sqrt(2)*sigma*math.exp(-2*sigma*sigma)/math.sqrt(math.pi))
    cutoff = 4*math.sqrt(3*erfc((80-a)/(math.sqrt(2)*sigma)))
    boundary = 1e-10
    z = (m2*tau+math.sqrt(3)*a*k)**2/(2*tau)
    rho = math.pi*z/(4*m1)+2*tail+4*alias+4*cutoff+boundary
    b = 3/(4*m1*m1*sigma*sigma)+rho
    error = b+math.sqrt(b)+4*(e+d+4*alias+2*tail)+8*cutoff+boundary
    delta = math.sin(m1*tau)
    margin = math.ceil(math.log(4/(delta*boundary))/math.log1p(delta/2))
    return dict(a=a, masses=[m1, m2], sigma=sigma, velocity=speed, c=3., extent=80.,
                boundary_error=boundary, cutoff=k, horizon=t, probability_error=float(error),
                blocking_upper=float(b), boundary_margin=margin,
                interval=dict(probability_upper='0.038', observable_swing_lower='0.923'))


def candidate():
    return dict(packets=packets(), witness=witness())
