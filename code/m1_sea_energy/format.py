"""Strict JSON primitives; structural checks survive optimized Python."""

import math
import numpy as np


def need(ok, message):
    if not ok:
        raise ValueError(message)


def keys(row, fields):
    need(type(row) is dict and set(row) == set(fields.split()), 'fields: '+fields)


def exact(actual, expected):
    need(type(actual) is type(expected), 'JSON type')
    if isinstance(expected, dict):
        need(actual.keys() == expected.keys(), 'dictionary coverage')
        for key in expected:
            exact(actual[key], expected[key])
    elif isinstance(expected, list):
        need(len(actual) == len(expected), 'list coverage')
        for a, b in zip(actual, expected):
            exact(a, b)
    elif type(expected) is float:
        need(math.isfinite(actual) and abs(actual-expected) <= 2e-9*abs(expected), 'relative scalar')
    else:
        need(actual == expected, 'exact value')


def pack(a):
    a = np.asarray(a, complex)
    return np.stack((a.real, a.imag), axis=-1).tolist()


def unpack(value, shape):
    def finite(x):
        return all(finite(a) for a in x) if type(x) is list else type(x) is float and math.isfinite(x)
    need(type(value) is list and finite(value), 'finite float tensor')
    a = np.asarray(value)
    need(a.shape == (*shape, 2), 'complex shape')
    return a[..., 0]+1j*a[..., 1]


def close(a, b, label, tol=2e-9):
    a, b = np.asarray(a), np.asarray(b)
    need(a.shape == b.shape and a.dtype.kind in 'iufc' and np.isfinite(a).all()
         and np.max(np.abs(a-b), initial=0) <= tol, label)


def real_vector(value, size):
    need(type(value) is list and len(value) == size
         and all(type(x) is float and math.isfinite(x) for x in value), 'real vector')
    return np.asarray(value)
