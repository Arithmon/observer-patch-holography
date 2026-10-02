"""Complete native reductions, staged-error controls and public instruments."""

import math
import numpy as np
from m1_source_realization.model import code_transfer_kraus
from .graphs import walks
from .format import pack


def reductions():
    rows = []
    for d in (2, 4):
        continuation = np.zeros((d, d+1))
        continuation[:, :d] = np.eye(d)
        continuation[0, d] = 1
        maps = [continuation @ k for k in code_transfer_kraus(d)]
        rows.append(dict(d=d, groups=[[pack(maps[0])], [pack(k) for k in maps[1:]]]))
    return rows


def stage_controls():
    rows = []
    x = np.array([[0, 1], [1, 0]])
    for k in (1, 2, 4, 8):
        eta, stages, u = .025, [], np.eye(2, dtype=complex)
        for _ in range(k):
            stages.append(eta)
            u = (math.cos(eta)*np.eye(2)-1j*math.sin(eta)*x) @ u
            eta = 2*eta*eta
        rows.append(dict(levels=k, strengths=stages, unitary=pack(u),
                         difference_norm=float(np.linalg.norm(u-np.eye(2), 2)),
                         uniform_sum_bound=.025/(1-.05),
                         raw_serial_flip=(1-(1-.05)**(3**k))/2))
    return rows


def instruments():
    _, b = walks(3)
    d = 8**12
    h = np.array([[1, 1], [1, -1]])/math.sqrt(2)
    us = [np.array([[math.cos(.37), -math.sin(.37)], [math.sin(.37), math.cos(.37)]]),
          np.diag([np.exp(-.61j), np.exp(.61j)]) @ h]
    rows = []
    for p in (.17, .31):
        blocks = [[np.zeros((16, 16), complex) for _ in (0, 1)] for _ in (0, 1)]
        wrong = 0.
        for label in (0, 1):
            projector = np.diag([int(label == 0), int(label == 1)])
            k = np.kron(us[label], projector)
            vec = k.flatten(order='F')
            for mask in range(512):
                word = np.array([mask >> i & 1 for i in range(9)], dtype=np.int64)
                refreshed = (2*(b @ word) > d).astype(int)
                visible = int(2*int(refreshed.sum()) > 9)
                errors = (mask ^ (511*label)).bit_count()
                probability = p**errors*(1-p)**(9-errors)
                blocks[label][visible] += probability*np.outer(vec, vec.conj())
                wrong += probability*int(visible != label)/2
        rows.append(dict(p=p, bits=9, abort_label=1,
                         blocks=[[pack(a) for a in line] for line in blocks],
                         wrong_label_probability=wrong, shared_broadcast_probability=p))
    return rows


def fixed_noise():
    rows = []
    for p in (.01, .07):
        ks = [math.sqrt(1-p)*np.eye(6)]
        for j in range(6):
            k = np.zeros((6, 6))
            k[5, j] = math.sqrt(p)
            ks.append(k)
        v = np.concatenate(ks, axis=0)
        ideal = np.concatenate([np.eye(6), np.zeros((36, 6))], axis=0)
        rows.append(dict(p=p, kraus=[pack(k) for k in ks],
                         dilation_distance=float(np.linalg.norm(v-ideal, 2))))
    return rows


def candidate():
    return dict(reductions=reductions(), stages=stage_controls(), instruments=instruments(),
                fixed_noise=fixed_noise())
