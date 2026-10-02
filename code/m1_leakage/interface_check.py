"""Matrix-unit channel oracle; imports neither source transfers nor producers."""

import itertools
import numpy as np
from .format import close, integer, keys, need, unpack


def expected(d, i, j, flag):
    out = np.zeros((d, d), complex)
    if flag == 0 and i < d and j < d:
        out[i, j] = 1
    if flag == 1 and i == j and i >= d:
        out[0, 0] = 1
    return out


def image(maps, i, j):
    return sum(np.outer(k[:, i], k[:, j].conj()) for k in maps)


def check_groups(groups, d):
    need(type(groups) is list and len(groups) == 2, 'both public outcomes')
    maps = []
    for group, count in zip(groups, (1, 6-d)):
        need(type(group) is list and len(group) == count, 'all private branches')
        maps.append([unpack(k, (d, 6)) for k in group])
    close(sum(k.conj().T@k for group in maps for k in group), np.eye(6), 'complete reduction TP')
    for i, j, flag in itertools.product(range(6), range(6), range(2)):
        close(image(maps[flag], i, j), expected(d, i, j, flag), 'complete flagged reduction')
    return maps


def oracle_gate(name):
    if name == 'cz':
        return np.diag([-1., -1., -1., 1.])
    if name == 'h_t':
        return np.kron(np.array([[1, 1], [1, -1]])/np.sqrt(2),
                       np.diag([1, (1+1j)/np.sqrt(2)]))
    raise ValueError('native gate catalog')


def check_pair(row, name):
    keys(row, 'name branches')
    need(row['name'] == name, 'native gate order')
    need(type(row['branches']) is list and len(row['branches']) == 4, 'all pair flags')
    u = oracle_gate(name)
    complete = np.zeros((36, 36), complex)
    for branch, (a, b) in zip(row['branches'], itertools.product(range(2), repeat=2)):
        keys(branch, 'flags maps')
        need(type(branch['flags']) is list and all(type(x) is int for x in branch['flags'])
             and branch['flags'] == [a, b], 'literal pair flags')
        count = (4 if a else 1)*(4 if b else 1)
        need(type(branch['maps']) is list and len(branch['maps']) == count, 'complete private pair branches')
        ks = [unpack(k, (4, 36)) for k in branch['maps']]
        complete += sum(k.conj().T@k for k in ks)
        for i in range(36):
            for j in range(36):
                ideal = np.kron(expected(2, i//6, j//6, a), expected(2, i % 6, j % 6, b))
                close(image(ks, i, j), u@ideal@u.conj().T, 'complete two-carrier gate instrument')
    close(complete, np.eye(36), 'two-carrier TP')


def verify(row):
    keys(row, 'reductions pair_gates')
    need(type(row['reductions']) is list and len(row['reductions']) == 2, 'native reduction catalog')
    for item, d in zip(row['reductions'], (2, 4)):
        keys(item, 'd groups')
        integer(item['d'], d, d)
        check_groups(item['groups'], d)
    need(type(row['pair_gates']) is list and len(row['pair_gates']) == 2, 'pair gate catalog')
    for item, name in zip(row['pair_gates'], ('cz', 'h_t')):
        check_pair(item, name)
