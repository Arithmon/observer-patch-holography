"""Independent Hermitian-sign replay, direct occupations and reflecting flights.

No producer is imported. Re(U) gives |h| and sign Im(U*) gives its sign;
this is independent of the producer's complex Schur logarithm.
"""

import itertools
import numpy as np
from scipy.linalg import expm
from m1_operational_clocks.check import gates, frames
from m1_fermionic_source.check import annihilator, exterior
from m1_fermionic_source.geometry_check import closed_walk
from .format import need, keys, exact, close, unpack, real_vector


def sea(u, tau):
    s, w = np.linalg.eigh((u.conj().T-u)/(2j))
    sign = (w*np.sign(s))@w.conj().T
    r, v = np.linalg.eigh((u+u.conj().T)/2)
    absolute = (v*(np.arccos(np.clip(r, -1, 1))/tau))@v.conj().T
    return sign@absolute, (np.eye(len(u))-sign)/2


def occupations(h):
    n = len(h)
    # Direct occupation-bit action, avoiding a product of dense CAR matrices.
    out = np.zeros((1 << n, 1 << n), complex)
    for bits in range(1 << n):
        for j in range(n):
            if not (bits >> j & 1):
                continue
            rest = bits ^ (1 << j)
            sj = (-1)**((bits & ((1 << j)-1)).bit_count())
            for i in range(n):
                if not (rest >> i & 1):
                    si = (-1)**((rest & ((1 << i)-1)).bit_count())
                    out[rest | (1 << i), bits] += h[i, j]*sj*si
    return out


def fock_check(u, h, p, tau):
    ev = np.linalg.eigvalsh(h)
    vacuum = float(np.trace(h@p).real)
    many = occupations(h)-vacuum*np.eye(256)
    values = np.linalg.eigvalsh(many)
    expected = sorted(sum(abs(ev[j]) for j in range(8) if bits >> j & 1) for bits in range(256))
    close(values, expected, 'all Fock excitation energies', 8e-8)
    need(abs(values[0]) < 1e-8 and values[1] > .1, 'unique sea ground state')
    lifted = exterior(u)
    close(expm(-1j*tau*many), np.exp(1j*tau*vacuum)*lifted, 'full Fock source identity', 1e-8)
    # Replacing signed h by |h| on the old occupations fails even in one sector.
    sign = np.eye(8)-2*p
    need(np.linalg.norm(expm(-1j*tau*(sign@h))-u) > .01, 'wrong absolute-generator control')
    empty_branch = h+2*np.pi/tau*p
    same_sea_branch = h+2*np.pi/tau*sign
    close(expm(-1j*tau*empty_branch), u, 'empty-ground branch has same step')
    close(expm(-1j*tau*same_sea_branch), u, 'same-sea branch has same step')
    need(np.min(np.linalg.eigvalsh(empty_branch)) > 0, 'empty-ground branch positive')
    shifted, shifted_p = sea(expm(-1j*tau*same_sea_branch), tau)
    close(shifted_p, p, 'branch channel retains sea')
    close(shifted, h, 'principal reconstruction cannot detect integer branch')
    need(np.logaddexp(0, -abs(np.linalg.eigvalsh(same_sea_branch))).sum()
         < np.logaddexp(0, -abs(ev)).sum()/2, 'partition functions distinguish branch')


def verify_blocks(rows):
    catalog = list(itertools.product((.2, .025), (.7, 1.1), ((0., 0., 0.), (.3, -.4, .7))))
    need(type(rows) is list and len(rows) == len(catalog), 'complete block catalog')
    for number, (row, (a, mass, k)) in enumerate(zip(rows, catalog)):
        keys(row, 'a mass k unitary h p energies vacuum_energy logz phase_energy excitations')
        exact([row['a'], row['mass'], row['k']], [a, mass, list(k)])
        tau = a/np.sqrt(3)
        u = gates(k, a, mass, 3.)
        h, p = sea(u, tau)
        for name, target in [('unitary', u), ('h', h), ('p', p)]:
            close(unpack(row[name], (8, 8)), target, name)
        close(u.conj().T@u, np.eye(8), 'unitarity')
        close(p@p, p, 'sea projector')
        close(p@h, h@p, 'stationary vacuum')
        ev = np.linalg.eigvalsh(h)
        close(real_vector(row['energies'], 8), ev, 'complete signed spectrum')
        need(np.min(abs(ev)) >= mass-1e-9 and np.max(abs(ev)) <= np.pi/tau-mass+1e-9, 'both gaps')
        exact(row['vacuum_energy'], float(np.trace(h@p).real))
        exact(row['logz'], [float(np.logaddexp(0, -b*abs(ev)).sum()) for b in (.3, 1.2)])
        f = np.arange(1, 9)*np.exp(.37j*np.arange(8))
        f /= np.linalg.norm(f)
        pf, qf = p@f, (np.eye(8)-p)@f
        weight = float(np.vdot(pf, pf).real)
        positive = float(np.vdot(qf, h@qf).real)
        negative = float(-np.vdot(pf, h@pf).real)
        cost = 4*np.sin(.63/2)**2*(weight*positive+(1-weight)*negative)
        exact(row['phase_energy'], float(cost))
        exact(row['excitations'], float(8*np.sin(.63/2)**2*weight*(1-weight)))
        need(0 < cost <= 4*(positive+negative), 'phase energy positive and bounded')
        if number in (0, 3, 6):
            fock_check(u, h, p, tau)


def verify_boundary(rows):
    need(type(rows) is list and len(rows) == 3, 'reflecting catalog')
    _, _, beta, _, _ = frames()
    for row, side in zip(rows, (1, 2, 3)):
        keys(row, 'side a modes central_columns vacuum_energy gap')
        exact([row['side'], row['a'], row['modes']], [side, .13, 16*side**3])
        u = closed_walk(side, .13)
        tau = .13/np.sqrt(3)
        n = len(u)
        b = np.kron(np.eye(2*side**3), beta)
        m = np.kron(np.eye(side**3), np.diag(np.repeat([.7, 1.1], 8)))
        u0 = expm(1j*tau*m@b)@u
        close(b@u0@b, u0.conj().T, 'actual reflected chiral symmetry')
        close((u+u.conj().T)/2, np.cos(tau*m.diagonal())[:, None]*(u0+u0.conj().T)/2, 'both boundary gaps')
        h, p = sea(u, tau)
        center = 16*((side//2)*side**2+(side//2)*side+side//2)
        close(unpack(row['central_columns'], (n, 16)), p[:, center:center+16], 'all central covariance columns', 1e-8)
        if side == 1:
            # Every flight reflects: this particular cube has P=0 exactly.
            # Boundary filling is not assumed to equal the periodic half filling.
            exact(row['vacuum_energy'], 0.)
            close(p, np.zeros((n, n)), 'one-cell empty negative band')
        else:
            exact(row['vacuum_energy'], float(np.trace(h@p).real))
        gap = float(np.min(abs(np.linalg.eigvalsh((u.conj().T-u)/(2j)))))
        exact(row['gap'], gap)
        need(gap >= np.sin(.7*tau)-1e-10, 'reflected gap lower bound')


def verify_projection(rows):
    catalog = list(itertools.product((.08, .02, .005), ((.4, -.2, .1), (1.2, .3, -.6)), (.7, 1.1), (0, 1)))
    need(type(rows) is list and len(rows) == len(catalog), 'projection catalog')
    _, _, beta, carrier, alpha = frames()
    for row, (a, k, mass, valley) in zip(rows, catalog):
        keys(row, 'a mass k valley error bound ordinary_dirac_error')
        exact([row['a'], row['mass'], row['k'], row['valley']], [a, mass, list(k), valley])
        tau = a/np.sqrt(3)
        u = gates(np.asarray(k)+[valley*np.pi/a, 0., 0.], a, mass, 3.)
        _, p = sea(u, tau)
        h0 = mass*beta+sum(x*y for x, y in zip(k, alpha))
        e = np.sqrt(mass**2+np.dot(k, k))
        limiting = (np.eye(8)-(-1)**valley*carrier@h0/e)/2
        error = float(np.linalg.norm(p-limiting, 2))
        bound = float(min(1., np.pi*(mass*tau+np.sqrt(3)*a*np.linalg.norm(k))**2/(8*tau*mass)))
        wrong = float(np.linalg.norm(p-(np.eye(8)-h0/e)/2, 2))
        exact([row['error'], row['bound'], row['ordinary_dirac_error']], [error, bound, wrong])
        need(error <= bound and wrong > .9, 'convergent full projection, rejected wrong high band')


def verify(row):
    keys(row, 'blocks boundary projection')
    verify_blocks(row['blocks'])
    verify_boundary(row['boundary'])
    verify_projection(row['projection'])
