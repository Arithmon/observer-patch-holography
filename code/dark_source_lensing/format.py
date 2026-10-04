"""Strict receipt parsing; mathematical inputs are bounded canonical strings."""
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path


def need(ok, message):
    if not ok:
        raise ValueError(message)


def keys(obj, names):
    need(type(obj) is dict and set(obj) == set(names.split()), 'exact object keys')


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('ascii')


def digest(obj):
    return hashlib.sha256(canonical(obj)).hexdigest()


def equal(actual, expected, label):
    need(canonical(actual) == canonical(expected), label)


def rational(raw):
    need(type(raw) is str and len(raw) <= 100 and
         re.fullmatch(r'-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?', raw) is not None,
         'canonical rational string')
    value = Fraction(raw)
    need(str(value) == raw, 'reduced rational required')
    return value


def decimal(ctx, raw):
    need(type(raw) is str and len(raw) <= 100 and
         re.fullmatch(r'-?(0|[1-9][0-9]*)(\.[0-9]+)?(e-?[0-9]+)?', raw) is not None,
         'finite decimal string')
    if 'e' in raw:
        need(abs(int(raw.split('e')[1])) <= 200, 'bounded decimal exponent')
    value = ctx.mpf(raw)
    need(ctx.isfinite(value), 'finite decimal')
    return value


def load(path):
    need(Path(path).stat().st_size <= 250_000, 'receipt size limit')
    raw = Path(path).read_bytes()
    need(len(raw) <= 250_000, 'receipt size limit')
    def pairs(items):
        result = {}
        for key, value in items:
            need(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def forbidden(value):
        raise ValueError('floating/nonfinite JSON token forbidden')
    def integer(value):
        need(len(value) <= 50, 'bounded integer')
        return int(value)
    return json.loads(raw.decode('ascii'), object_pairs_hook=pairs,
                      parse_float=forbidden, parse_constant=forbidden, parse_int=integer)
