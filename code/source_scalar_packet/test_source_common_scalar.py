"""Independent arithmetic, domain, full-leakage and forged-receipt controls."""
from copy import deepcopy
from fractions import Fraction as F
from math import factorial
import json

import pytest

import source_scalar_intervals as iv
import source_common_scalar as producer
import verify_source_common_scalar as verifier


def encloses(interval, value):
    return F(interval.lo, iv.SCALE) <= value <= F(interval.hi, iv.SCALE)


@pytest.mark.parametrize('a,b', [(F(1, 3), F(7, 11)), (F(-31, 13), F(5, 17)),
                                (F(-9, 7), F(-11, 3)), (F(0), F(3, 8))])
def test_directed_rational_arithmetic(a, b):
    x, y = iv.I.of(a), iv.I.of(b)
    assert encloses(x+y, a+b)
    assert encloses(x-y, a-b)
    assert encloses(x*y, a*b)
    assert encloses(x/y, a/b)


@pytest.mark.parametrize('a', [F(0), F(1, 3), F(2), F(10001, 7)])
def test_square_root_by_exact_squaring(a):
    z = iv.I.of(a).sqrt()
    assert F(z.lo, iv.SCALE)**2 <= a <= F(z.hi, iv.SCALE)**2


@pytest.mark.parametrize('bad', [True, False, 1.0, float('nan'), 1+2j, '1'])
def test_reject_inexact_inputs(bad):
    with pytest.raises(ValueError):
        iv.I.of(bad)


def test_zero_divisor_and_negative_domain():
    with pytest.raises(ValueError):
        iv.I.of(1)/iv.I(-1, 1)
    with pytest.raises(ValueError):
        iv.I.of(-1).sqrt()


def test_transcendentals_against_separate_exact_series():
    # A distinct alternating sine series and positive exponential series.
    for x in [F(1, 7), F(-2, 5), F(3, 2)]:
        a = sum(((-1)**k*x**(2*k+1)/factorial(2*k+1) for k in range(60)), F())
        rem = abs(x)**121/factorial(121)
        got = iv.sin(iv.I.of(x))
        assert F(got.lo, iv.SCALE) <= a-rem <= a+rem <= F(got.hi, iv.SCALE)
    x = F(1, 3)
    a = sum((x**k/factorial(k) for k in range(100)), F())
    rem = x**100/F(factorial(100))/(1-x/101)
    got = iv.exp(iv.I.of(x))
    assert F(got.lo, iv.SCALE) <= a <= a+rem <= F(got.hi, iv.SCALE)
    pi_lo = F('3.14159265358979323846264338327950288419716939937510')
    pi_hi = pi_lo+F(1, 10**50)
    assert iv.PI.lo/iv.SCALE > 3  # only a coarse noncertificate display check
    assert F(iv.PI.lo, iv.SCALE) <= pi_lo <= pi_hi <= F(iv.PI.hi, iv.SCALE)


@pytest.mark.parametrize('cut', [F(9, 20), F(11, 20)])
def test_ambiguous_slab_node_fails_closed(cut):
    bound = iv.I.of(cut)
    with pytest.raises(ValueError, match='slab'):
        producer.above(iv.I(bound.lo-1, bound.hi+1), cut)
    assert producer.above(iv.I.of(cut+F(1, 100)), cut)
    assert not producer.above(iv.I.of(cut-F(1, 100)), cut)


def test_three_by_three_resolvent_including_outside_band():
    # Positive A has a one-dimensional compressed eigenvalue2 but leaks to
    # the next node. A compressed residual would be zero, the full form isn't.
    diagonal = [iv.I.of(3), iv.I.of(4), iv.I.of(3)]
    off = [iv.I.of(-1), iv.I.of(-1)]
    rhs = [iv.I.of(0), iv.I.of(-1), iv.I.of(0)]
    a = producer.solve(diagonal, off, rhs)
    b = verifier.reverse_solve(diagonal, off, rhs)
    exact = [F(-1, 10), F(-3, 10), F(-1, 10)]
    assert all(encloses(x, y) for x, y in zip(a, exact))
    assert all(encloses(x, y) for x, y in zip(b, exact))
    assert encloses(iv.dot(rhs, a), F(3, 10))


@pytest.mark.parametrize('bad', [3, 8, 12, 15, True, 13.0])
def test_frozen_refinement_domain(bad):
    with pytest.raises(ValueError):
        producer.grid(bad)


def test_fresh_independent_receipt_replay():
    result = verifier.verify(verifier.load())
    assert result['verified']
    assert result['resolved_q'] == 233
    assert F(result['error_upper']) < F(result['resolved_signal_lower'])


@pytest.mark.parametrize('raw', [b'{"a":1,"a":2}', b'{"a":0.0}', b'{"a":NaN}',
                                b' '*500001], ids=['duplicate', 'float', 'nonfinite', 'oversize'])
def test_strict_receipt_loader(raw, tmp_path):
    path = tmp_path/'bad.json'
    path.write_bytes(raw)
    with pytest.raises((ValueError, UnicodeError)):
        verifier.load(path)


@pytest.mark.parametrize('kind', ['zero_leakage', 'zero_tail', 'forced_coarse_pass',
                                  'false_clock', 'wrong_coordinate', 'wrong_energy'])
def test_resealed_scientific_forgeries_rejected(kind):
    data = deepcopy(verifier.load())
    if kind == 'zero_leakage':
        data['levels'][-1]['full_resolvent_column_energies'][-1] = [['0', '0'], ['0', '0']]
    elif kind == 'zero_tail':
        data['levels'][-1]['sampled_detector_tail_mass'] = ['0', '0']
    elif kind == 'forced_coarse_pass':
        data['levels'][0]['error_smaller_than_resolved_signal'] = True
    elif kind == 'false_clock':
        data['scope']['count_clock_custody_for_field_execution_established'] = True
    elif kind == 'wrong_coordinate':
        data['levels'][-1]['ordered_source_coordinates_Qphi'][-1][0] += 1
    elif kind == 'wrong_energy':
        data['levels'][0]['compact_preparation_energy'] = ['8', '8']
    # JSON roundtrip removes any incidental in-memory type distinction; all
    # source pins remain valid, so rejection must concern scientific contents.
    data = json.loads(json.dumps(data))
    with pytest.raises(ValueError):
        verifier.verify(data)


def test_boolean_numeric_alias_rejected():
    data = verifier.load()
    data['preparation']['momentum_impulse'] = True
    with pytest.raises(ValueError):
        verifier.verify(data)
