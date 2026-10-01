"""Strict finite receipt primitives (checks remain active under python -O)."""

import hashlib
import json
import math
import numpy as np


def need(condition, message):
    if not condition:
        raise ValueError(message)


def keys(value, fields):
    need(type(value) is dict and set(value) == set(fields.split()), 'exact fields: '+fields)


def integer(value, low, high):
    need(type(value) is int and low <= value <= high, 'bounded integer')
    return value


def number(value):
    need(type(value) in (int, float) and math.isfinite(value), 'finite real number')
    return value


def close(value, expected, name):
    actual = np.asarray(value)
    target = np.asarray(expected)
    need(actual.dtype.kind in 'iufc' and actual.shape == target.shape
         and np.isfinite(actual).all() and np.max(np.abs(actual-target), initial=0) <= 2e-11,
         name)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode('ascii')).hexdigest()


def pack(value):
    a = np.asarray(value, complex)
    return np.stack((a.real, a.imag), axis=-1).tolist()


def unpack(value, shape):
    def finite(item):
        return all(finite(x) for x in item) if type(item) is list else (
            type(item) in (int, float) and math.isfinite(item))
    need(type(value) is list and finite(value), 'numeric tensor')
    a = np.asarray(value)
    need(a.shape == (*shape, 2), 'complex tensor shape')
    return a[..., 0]+1j*a[..., 1]
