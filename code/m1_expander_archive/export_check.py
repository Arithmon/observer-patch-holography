"""Matrix-unit and analytic Choi checks, independent of the producer."""

import math
import numpy as np
from .format import keys, need, integer, number, close, unpack
from .graph_check import reference


def check_groups(groups, d):
    need(type(groups) is list and len(groups) == 2, 'success/failure public groups')
    need(all(type(g) is list for g in groups) and [len(g) for g in groups] == [1, 6-d],
         'complete private Kraus multiplicities')
    matrices = [[unpack(k, (d, 6)) for k in group] for group in groups]
    close(sum(k.conj().T @ k for g in matrices for k in g), np.eye(6), 'complete TP reduction')
    for i in range(6):
        for j in range(6):
            for outcome, group in enumerate(matrices):
                actual = sum(np.outer(k[:, i], k[:, j].conj()) for k in group)
                target = np.zeros((d, d), complex)
                if outcome == 0 and i < d and j < d:
                    target[i, j] = 1
                elif outcome == 1 and i == j and i >= d:
                    target[0, 0] = 1
                close(actual, target, 'native public instrument on every matrix unit')


def verify_reductions(rows):
    need(type(rows) is list and len(rows) == 2, 'complete native code catalog')
    for row, d in zip(rows, (2, 4)):
        keys(row, 'd groups')
        integer(row['d'], d, d)
        check_groups(row['groups'], d)


def verify_stages(rows):
    need(type(rows) is list and len(rows) == 4, 'complete staged-decoding norm controls')
    for row, k in zip(rows, (1, 2, 4, 8)):
        keys(row, 'levels strengths unitary difference_norm uniform_sum_bound raw_serial_flip')
        integer(row['levels'], k, k)
        need(type(row['strengths']) is list and len(row['strengths']) == k, 'all decoder stages')
        # Closed form, not the producer's recurrence or matrix products.
        strengths = [.05**(2**j)/2 for j in range(k)]
        for a, b in zip(row['strengths'], strengths):
            close(number(a), b, 'recursive stage strength')
        angle = sum(strengths)
        expected = np.array([[math.cos(angle), -1j*math.sin(angle)],
                             [-1j*math.sin(angle), math.cos(angle)]])
        close(unpack(row['unitary'], (2, 2)), expected, 'complete coherent stage composition')
        close(number(row['difference_norm']), 2*math.sin(angle/2), 'coherent operator distance')
        close(number(row['uniform_sum_bound']), 1/38, 'uniform geometric decoder bound')
        close(number(row['raw_serial_flip']), -math.expm1(3**k*math.log(.95))/2,
              'unprotected growing decoder exposure')
        need(angle <= row['uniform_sum_bound'], 'depth-independent bound')
    need(rows[-1]['raw_serial_flip'] > .499, 'raw decoder countercontrol must fail')


def verify_instruments(rows):
    need(type(rows) is list and len(rows) == 2, 'complete finite channel experiments')
    _, b = reference(3)
    for mask in range(512):
        votes = b @ np.array([mask >> j & 1 for j in range(9)], dtype=np.int64)
        need(np.all((2*votes > 8**12) == (mask.bit_count() >= 5)), 'finite powered majority map')
    for row, p in zip(rows, (.17, .31)):
        keys(row, 'p bits abort_label blocks wrong_label_probability shared_broadcast_probability')
        close(number(row['p']), p, 'fixed instrument noise control')
        integer(row['bits'], 9, 9)
        integer(row['abort_label'], 1, 1)
        need(type(row['blocks']) is list and len(row['blocks']) == 2
             and all(type(line) is list and len(line) == 2 for line in row['blocks']),
             'all true-abort/public-label branches')
        error = sum(math.comb(9, w)*p**w*(1-p)**(9-w) for w in range(5, 10))
        total = np.zeros((16, 16), complex)
        for true in (0, 1):
            # Explicit conditional signal/spectator map, including coherences.
            unitary = (np.array([[math.cos(.37), -math.sin(.37)], [math.sin(.37), math.cos(.37)]])
                       if true == 0 else np.array([[np.exp(-.61j), np.exp(-.61j)],
                                                   [np.exp(.61j), -np.exp(.61j)]])/math.sqrt(2))
            vec = np.zeros(16, complex)
            for source in range(2):
                for target in range(2):
                    i, o = 2*source+true, 2*target+true
                    vec[4*i+o] = unitary[target, source]
            ideal = np.outer(vec, vec.conj())
            for visible in (0, 1):
                actual = unpack(row['blocks'][true][visible], (16, 16))
                close(actual, (1-error if true == visible else error)*ideal,
                      'complete conditional quantum/public instrument')
                need(np.linalg.eigvalsh(actual).min() >= -2e-11, 'positive instrument branch')
                total += actual
        close(np.einsum('iaja->ij', total.reshape(4, 4, 4, 4)), np.eye(4), 'unconditional TP export')
        close(number(row['wrong_label_probability']), error, 'retained wrong-label mass')
        close(number(row['shared_broadcast_probability']), p, 'shared raw-bit floor')
        need(0 < error < p and np.trace(unpack(row['blocks'][1][0], (16, 16))).real > 1e-4,
             'wrong abort label retained, not selected away')


def verify(row):
    keys(row, 'reductions stages instruments fixed_noise')
    verify_reductions(row['reductions'])
    verify_stages(row['stages'])
    verify_instruments(row['instruments'])
    rows = row['fixed_noise']
    need(type(rows) is list and len(rows) == 2, 'nonzero full-carrier noise catalog')
    for item, p in zip(rows, (.01, .07)):
        keys(item, 'p kraus dilation_distance')
        close(number(item['p']), p, 'fixed noise strength control')
        need(type(item['kraus']) is list and len(item['kraus']) == 7, 'complete erasure environment')
        ks = [unpack(k, (6, 6)) for k in item['kraus']]
        close(sum(k.conj().T @ k for k in ks), np.eye(6), 'TP fixed-strength noise')
        for i in range(6):
            for j in range(6):
                target = np.zeros((6, 6))
                target[i, j] = 1-p
                if i == j:
                    target[5, 5] += p
                close(sum(np.outer(k[:, i], k[:, j].conj()) for k in ks), target,
                      'full-M6 nonzero leakage process')
        v = np.vstack(ks)
        ideal = np.zeros((42, 6))
        ideal[:6] = np.eye(6)
        distance = math.sqrt(2-2*math.sqrt(1-p))
        close(np.linalg.norm(v-ideal, 2), distance, 'actual dilation alignment')
        close(number(item['dilation_distance']), distance, 'fixed noise local amplitude')
        need(0 < distance <= math.sqrt(2*p), 'nonempty noise strength interval')
