"""Independent full-M6 matrix-unit oracle, without producer Kraus operators."""

import numpy as np
from .independent import keys, need, close, integer


def decode_matrix(row):
    keys(row, 'code_dimension environment real imag')
    integer(row['code_dimension'], 2, 4)
    d = row['code_dimension']
    need(d in (2, 4), 'native proper-code dimension')
    integer(row['environment'], 7-d, 7-d)
    shape = ((7-d)*2*(d+1), 6)
    arrays = []
    for part in ('real', 'imag'):
        value = row[part]
        need(type(value) is list and len(value) == shape[0]
             and all(type(line) is list and len(line) == shape[1]
                     and all(type(x) in (int, float) and np.isfinite(x) for x in line)
                     for line in value), 'finite rectangular matrix')
        arrays.append(np.array(value, dtype=float))
    return arrays[0]+1j*arrays[1]


def check_dilation(v, d):
    environment, public, out = 7-d, 2, d+1
    close(v.conj().T @ v, np.eye(6), 'complete native isometry')
    blocks = v.reshape(environment, public*out, 6)
    # Every matrix unit includes coherence between code and complement, not
    # just the protected input sector or the diagonal failure probabilities.
    for i in range(6):
        for j in range(6):
            actual = sum(np.outer(k[:, i], k[:, j].conj()) for k in blocks)
            expected = np.zeros((public*out, public*out), complex)
            if i < d and j < d:
                expected[i, j] = 1
            elif i == j and i >= d:
                expected[2*d+1, 2*d+1] = 1
            close(actual, expected, 'complete grouped native instrument')


def verify(rows):
    need(type(rows) is list and len(rows) == 2, 'complete native instrument catalog')
    for row, d in zip(rows, (2, 4)):
        v = decode_matrix(row)
        need(row['code_dimension'] == d, 'frozen native code catalog')
        check_dilation(v, d)
