"""Strict compact numeric serialization; no scientific model is shared here."""

import math
import numpy as np


def need(condition, message):
    if not condition:
        raise ValueError(message)


def keys(value, expected):
    need(type(value) is dict and set(value) == set(expected.split()), 'exact object keys')


def integer(value, low, high):
    need(type(value) is int and low <= value <= high, 'integer range')


def number(value):
    need(type(value) in (int, float), 'literal finite number')
    try:
        good = math.isfinite(value)
    except (OverflowError, TypeError):
        good = False
    need(good, 'literal finite number')
    return float(value)


def close(actual, expected, message, tolerance=3e-10):
    actual, expected = np.asarray(actual), np.asarray(expected)
    need(actual.shape == expected.shape and np.isfinite(actual).all()
         and np.isfinite(expected).all()
         and np.max(np.abs(actual-expected), initial=0) < tolerance, message)


def pack(matrix):
    matrix = np.asarray(matrix, dtype=complex)
    need(matrix.ndim == 2 and np.isfinite(matrix).all(), 'finite matrix')
    entries = []
    for i, j in np.ndindex(matrix.shape):
        real, imag = round(float(matrix[i, j].real), 13), round(float(matrix[i, j].imag), 13)
        if real or imag:
            entries.append([i, j, real, imag])
    return dict(shape=list(matrix.shape), entries=entries)


def unpack(row, shape):
    keys(row, 'shape entries')
    need(type(row['shape']) is list and row['shape'] == list(shape)
         and all(type(x) is int for x in row['shape']), 'matrix dimensions')
    need(type(row['entries']) is list and len(row['entries']) <= math.prod(shape), 'sparse matrix entries')
    out = np.zeros(shape, complex)
    previous = (-1, -1)
    for entry in row['entries']:
        need(type(entry) is list and len(entry) == 4, 'matrix entry')
        i, j, real, imag = entry
        integer(i, 0, shape[0]-1)
        integer(j, 0, shape[1]-1)
        need((i, j) > previous, 'unique sorted entries')
        real, imag = number(real), number(imag)
        need(real != 0 or imag != 0, 'nonzero sparse entry')
        out[i, j] = real+1j*imag
        previous = (i, j)
    return out
