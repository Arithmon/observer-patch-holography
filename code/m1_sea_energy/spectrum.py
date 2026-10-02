"""Producer: principal Schur logarithms of the actual massive source walk."""

import itertools
import numpy as np
from scipy.linalg import schur
from m1_operational_clocks import model as walk
from m1_fermionic_source.walk import execute
from .format import pack


def sea(u, tau):
    t, z = schur(u, output='complex')
    energies = -np.angle(np.diag(t))/tau
    return (z*energies)@z.conj().T, z[:, energies < 0]@z[:, energies < 0].conj().T


def block_rows():
    result = []
    for a, mass, k in itertools.product((.2, .025), (.7, 1.1), ((0., 0., 0.), (.3, -.4, .7))):
        tau = np.sqrt(3)*a/3
        u = walk.walk(k, a, mass)
        h, p = sea(u, tau)
        ev = np.linalg.eigvalsh(h)
        f = np.exp(1j*np.arange(8)*.37)*np.arange(1, 9)
        f /= np.linalg.norm(f)
        v = np.eye(8)+(np.exp(.63j)-1)*np.outer(f, f.conj())
        hp = h@p
        energy = float(np.trace(h@(v@p@v.conj().T-p)).real)
        particles = float(np.trace((np.eye(8)-p)@v@p@v.conj().T).real)
        holes = float(np.trace(p@(np.eye(8)-v@p@v.conj().T)).real)
        result.append(dict(a=a, mass=mass, k=list(k), unitary=pack(u), h=pack(h), p=pack(p),
                           energies=ev.tolist(), vacuum_energy=float(np.trace(hp).real),
                           logz=[float(np.logaddexp(0, -b*abs(ev)).sum()) for b in (.3, 1.2)],
                           phase_energy=energy, excitations=particles+holes))
    return result


def boundary_rows():
    result = []
    # Full reflecting operators, including every internal mode and both masses.
    for side in (1, 2, 3):
        a = .13
        u, empty = execute(side, a)
        h, p = sea(u, np.sqrt(3)*a/3)
        center = ((side//2)*side**2+(side//2)*side+side//2)*16
        indices = list(range(center, center+16))
        # Store full columns, not only diagonal densities. Large matrices are recomputed.
        result.append(dict(side=side, a=a, modes=16*side**3,
                           central_columns=pack(p[:, indices]),
                           vacuum_energy=float(np.trace(h@p).real),
                           gap=float(np.min(abs(np.linalg.eigvalsh((u.conj().T-u)/(2j)))))))
    return result


def projection_rows():
    result = []
    for a, k, mass, valley in itertools.product((.08, .02, .005),
            ((.4, -.2, .1), (1.2, .3, -.6)), (.7, 1.1), (0, 1)):
        tau = np.sqrt(3)*a/3
        shifted = np.asarray(k)+np.array([valley*np.pi/a, 0., 0.])
        u = walk.walk(shifted, a, mass)
        _, p = sea(u, tau)
        h = walk.hamiltonian(k, mass)
        energy = np.sqrt(mass**2+np.dot(k, k))
        limiting = (np.eye(8)-(-1)**valley*walk.R@h/energy)/2
        wrong = (np.eye(8)-h/energy)/2
        z = (mass*tau+np.sqrt(3)*a*np.linalg.norm(k))**2/(2*tau)
        result.append(dict(a=a, mass=mass, k=list(k), valley=valley,
                           error=float(np.linalg.norm(p-limiting, 2)),
                           bound=float(min(1., np.pi*z/(4*mass))),
                           ordinary_dirac_error=float(np.linalg.norm(p-wrong, 2))))
    return result


def candidate():
    return dict(blocks=block_rows(), boundary=boundary_rows(), projection=projection_rows())
