"""Transfer the frozen earlier preparation; separately invert zero-time data."""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
import zlib

import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]


def transferred_parameters(protocol):
    path = ROOT/protocol['preparation']['prior_receipt']
    raw = path.read_bytes()
    if path.suffix == '.gz':
        try:
            raw = gzip.decompress(raw)
        except (OSError, EOFError, zlib.error) as error:
            raise ValueError('Frozen earlier calibration receipt has invalid gzip data') from error
    if hashlib.sha256(raw).hexdigest() != protocol['preparation']['prior_receipt_sha256']:
        raise ValueError('Frozen earlier calibration receipt changed')
    return json.loads(raw)['calibration']['parameters']


def zero_time_displacement(position, thermal, h=.4):
    """Principal inverse of the measured finite-difference estimator, not <x>."""
    if any(isinstance(v, bool) or not math.isfinite(v) for v in (position, thermal, h)):
        raise ValueError('finite real calibration inputs required')
    if h <= 0 or thermal < 0:
        raise ValueError('positive probe and nonnegative thermal occupation required')
    with mp.workdps(80):
        x, n, b = (mp.mpf(float(v)) for v in (position, thermal, h))
        argument = x*mp.sqrt(2)*b*mp.exp((n+mp.mpf('.5'))*b*b)
        if abs(argument) > 1:
            raise ValueError('zero-time estimator is outside the declared thermal family')
        return float(mp.asin(argument)/(2*b))


def estimator(alpha, thermal, h=.4):
    return math.exp(-(thermal+.5)*h*h)*math.sin(2*alpha*h)/(math.sqrt(2)*h)
