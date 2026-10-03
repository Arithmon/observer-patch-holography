"""Occupation minors, all preparation branches, both returned QND outcomes."""

import itertools
import math
import numpy as np
from m1_fermionic_source.check import annihilator, at, X, Z, word, exterior
from m1_fermionic_source.experiment_check import replay_native, replay_rotations
from .format import keys, need, exact, close, unpack, real_vector
from .spectrum_check import occupations


def replay(rows, parity, iso, basis, edges, loops):
    need(type(rows) is list and 0 < len(rows) <= 40, 'bounded nonempty source circuit')
    single = np.eye(4, dtype=complex)
    physical = np.eye(16, dtype=complex)
    for row in rows:
        keys(row, 'modes matrix scalar rotations')
        modes = row['modes']
        need(type(modes) is list and len(modes) in (1, 2) and len(set(modes)) == len(modes)
             and all(type(i) is int and 0 <= i < 4 for i in modes), 'one or two source modes')
        if len(modes) == 2:
            need(sorted(modes) in edges, 'only actual edges, including routing')
        g = unpack(row['matrix'], (len(modes), len(modes)))
        close(g.conj().T@g, np.eye(len(modes)), 'unitary source gate')
        counts = (1,) if len(modes) == 1 else (6, 8)
        allowed = sum(1 << e for e, pair in enumerate(edges) if any(v in pair for v in modes))
        u = replay_rotations({k: row[k] for k in ('scalar', 'rotations')}, 4, loops, counts, allowed)
        lifted = np.eye(4, dtype=complex)
        lifted[np.ix_(modes, modes)] = g
        close(u@iso, iso@exterior(lifted)[np.ix_(basis, basis)], 'source gate on entire parity sector')
        single = lifted@single
        physical = u@physical
    return single, physical


def verify_case(row, rank):
    keys(row, 'rank vertices edges parity orbitals covariance basis isometry gates branches kernel read qnd kraus probabilities mean_read_energy pulse_rotations')
    edges = [[0, 1], [0, 2], [0, 3], [1, 2]]
    parity = rank % 2
    basis = [s for s in range(16) if s.bit_count() % 2 == parity]
    exact([row['rank'], row['vertices'], row['edges'], row['parity'], row['basis']], [rank, 4, edges, parity, basis])
    identity = np.eye(16)
    bq = [identity.copy() for _ in range(4)]
    for vertex in range(4):
        for e, pair in enumerate(edges):
            if vertex in pair:
                bq[vertex] = bq[vertex]@at(Z, e, 4)
    aq = []
    for e, (i, j) in enumerate(edges):
        a = at(X, e, 4)
        for earlier in range(e):
            if i in edges[earlier] or j in edges[earlier]:
                a = a@at(Z, earlier, 4)
        aq.append(a)
    loop = -1j*aq[0]@aq[3]@(-aq[1])  # i^3 A01 A12 A20
    close(loop@loop, identity, 'fundamental loop sign')
    c = [annihilator(4, i) for i in range(4)]
    majorana = [a+a.conj().T for a in c]
    iso = unpack(row['isometry'], (16, 8))
    close(iso.conj().T@iso, np.eye(8), 'isometry')
    for i in range(4):
        fermion = identity-2*c[i].conj().T@c[i]
        close((-1 if parity and i == 0 else 1)*bq[i]@iso, iso@fermion[np.ix_(basis, basis)], 'occupation intertwiner')
    for edge, (i, j) in zip(aq, edges):
        fermion = -1j*majorana[i]@majorana[j]
        close(edge@iso, iso@fermion[np.ix_(basis, basis)], 'all Majorana intertwiners')
    close(loop@iso, iso, 'entire code image')
    z = unpack(row['orbitals'], (4, 4))
    close(z.conj().T@z, np.eye(4), 'orthonormal complete orbitals')
    reference = np.array([[np.exp(2j*np.pi*i*j/4+.17j*i*i)/2 for j in range(4)] for i in range(4)])
    p = reference[:, :rank]@reference[:, :rank].conj().T
    close(z[:, :rank]@z[:, :rank].conj().T, p, 'same prescribed occupied subspace')
    close(unpack(row['covariance'], (4, 4)), p, 'prepared target covariance')
    single, physical = replay(row['gates'], parity, iso, basis, edges, [loop])
    close(single, z, 'compiled full orbital unitary')
    expected = exterior(z)[:, (1 << rank)-1]
    target = iso@expected[basis]
    need(type(row['branches']) is list and len(row['branches']) == 2, 'both loop outcomes retained')
    for bit, branch in enumerate(row['branches']):
        keys(branch, 'outcomes probability density')
        exact(branch['outcomes'], [bit])
        exact(branch['probability'], .5)
        seed = (identity+(-1)**bit*loop)@identity[:, 0]/2
        close(np.vdot(seed, seed), .5, 'actual loop probability')
        seed *= math.sqrt(2)
        if bit:
            seed = at(Z, 3, 4)@seed
        # Star edges directly pair root with every occupied nonroot vertex.
        for vertex in range(1, rank):
            seed = aq[vertex-1]@seed
        state = physical@seed
        density = np.outer(state, state.conj())
        close(unpack(branch['density'], (16, 16)), density, 'every corrected source preparation')
        close(density, np.outer(target, target.conj()), 'full Slater ray, including occupation signs')
    f = np.array([1., .3j, -.4+.2j, .1], complex)
    f /= np.linalg.norm(f)
    close(unpack(row['kernel'], (4,)), f, 'declared complex read mode')
    need(type(row['read']) is list and len(row['read']) == 3, 'all tree eliminations')
    gathering, code_gather = replay(row['read'], parity, iso, basis, edges, [loop])
    close(gathering@f, np.eye(4)[:, 0], 'actual mode gather')
    need(type(row['qnd']) is list and 0 < len(row['qnd']) <= 29, 'charged bounded native QND')
    qnd = replay_native(row['qnd'], 5)
    root = (-1 if parity else 1)*bq[0]
    close(qnd[:16, :16], (identity+root)/2, 'native zero branch')
    close(qnd[16:, :16], (identity-root)/2, 'native one branch')
    a = sum(f[j].conjugate()*c[j] for j in range(4))
    number = a.conj().T@a
    need(type(row['kraus']) is list and len(row['kraus']) == 2, 'both QND outcomes')
    probabilities = []
    for bit, effect in enumerate((identity-number, number)):
        actual = code_gather.conj().T@qnd[16*bit:16*(bit+1), :16]@code_gather
        close(unpack(row['kraus'][bit], (16, 16)), actual, 'returned complete QND Kraus')
        close(actual@iso, iso@effect[np.ix_(basis, basis)], 'entire sector and all spectators')
        probabilities.append(float(np.vdot(actual@target, actual@target).real))
    close(real_vector(row['probabilities'], 2), probabilities, 'unsubtracted probabilities')
    close(probabilities[1], np.vdot(f, p@f), 'Slater covariance is actual click probability')
    close(sum(probabilities), 1., 'no postselection')
    need(.01 < probabilities[1] < .99, 'nontrivial filled-background read')
    # Native Kraus maps and a direct occupation Hamiltonian test the *read's*
    # excitation cost, not only the energy of the prepared input state.
    h = (reference*np.array([-2., -1., .9, 1.4]))@reference.conj().T
    energy = occupations(h)+3*identity
    change = number@energy@number+(identity-number)@energy@(identity-number)-energy
    need(np.linalg.norm(change, 2) <= 2*np.linalg.norm(h@f)+1e-10, 'full Fock QND energy-change bound')
    encoded_energy = iso@energy[np.ix_(basis, basis)]@iso.conj().T
    weighted = [float(np.vdot(unpack(k, (16, 16))@target, encoded_energy@unpack(k, (16, 16))@target).real)
                for k in row['kraus']]
    need(min(weighted) >= -1e-10, 'both nonnegative branch energies retained')
    exact(row['mean_read_energy'], float(sum(weighted)))
    need(sum(weighted) <= np.vdot(expected, energy@expected).real+2*np.linalg.norm(h@f)+1e-10,
         'complete mean energy, without postselection')
    count = sum(len(g['rotations']) for g in row['gates'])+2*sum(len(g['rotations']) for g in row['read'])
    exact(row['pulse_rotations'], count)
    # Old any-occupied-mode detector can be completely blinded by the sea.
    old_effect = identity-exterior(np.eye(4)-p)
    close(np.vdot(expected, old_effect@expected), 1., 'old vacuum-click saturation control')


def verify(rows):
    need(type(rows) is list and len(rows) == 2, 'even sea and odd added-particle preparations')
    for row, rank in zip(rows, (2, 3)):
        verify_case(row, rank)
