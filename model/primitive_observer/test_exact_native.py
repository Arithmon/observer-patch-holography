"""Independent numerical witnesses and hostile controls for native certificates."""
from fractions import Fraction as F
from itertools import product

import mpmath
import numpy as np
import pytest

from ..common_source import quantum_native as source
from . import exact_native as exact


MODES = ('protected_qubit', 'passive_qubit')


def _mp(value, ctx):
    value = F(value)
    return ctx.mpf(value.numerator)/value.denominator


def _number(value, ctx):
    return ctx.mpc(_mp(float(value.real), ctx), _mp(float(value.imag), ctx))


@pytest.mark.parametrize('law,mode', tuple(product(('C', 'E'), MODES)))
def test_exact_adaptive_completeness_on_each_first_image(law, mode):
    result = exact.completeness(law, mode)
    assert result['current_weights'] == {-1: F(1, 4), 0: F(1, 2), 1: F(1, 4)}
    assert sum(result['first_record_weights'].values()) == 1
    for (coin, raw, selected), effect in result['conditional_final_effects'].items():
        weight = result['first_record_weights'][coin, raw]
        assert effect == {(0, 0): weight, (1, 1): weight, (0, 1): 0, (1, 0): 0}
        assert selected in (0, 1)
    assert result['arbitrary_finite_reference'] is True
    assert result['all_native_output_rows'] == 512


@pytest.mark.parametrize('law,mode', tuple(product(('C', 'E'), MODES)))
def test_complete_records_keep_zero_abort_targets_and_exact_integer_oracle(law, mode, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('rounded target entered exact certificate')
    monkeypatch.setattr(source, '_ideal_encoded_tree', forbidden)
    monkeypatch.setattr(source, '_logical_code', forbidden)
    exact._tree.cache_clear()
    tree = exact.ideal_tree(law, mode)
    assert all(type(weight) is F for row in tree.values() for weight in row.values())
    for record, abort in source._menu(law).items():
        row = exact.ideal_map(law, mode, record)
        if abort:
            assert row == {}
        elif record[2]:
            assert row == tree.get((record[0], record[1], record[2][2], *record[3:]), {})
    # Mutating a caller's sparse envelope cannot poison the cached target.
    key, row = next(iter(tree.items()))
    row.clear()
    assert exact.ideal_tree(law, mode)[key]


@pytest.mark.parametrize('mode', MODES)
def test_adaptive_completeness_rejects_deleted_or_changed_final_map(mode, monkeypatch):
    tree = exact.ideal_tree('C', mode)
    key = next(iter(tree))
    tree.pop(key)
    monkeypatch.setattr(exact, 'ideal_tree', lambda *args: tree)
    with pytest.raises(ArithmeticError, match='not complete'):
        exact.completeness('C', mode)


@pytest.mark.parametrize('law,mode', tuple(product(('C', 'E'), MODES)))
def test_original_executed_nonzero_and_abort_maps_receive_rational_bounds(law, mode):
    tree = source._native_tree(law, mode)
    first = next(row for row in tree if not row[3] and np.linalg.norm(row[4]) > .1)
    c1, raw1, _, _, _, endings = first
    for selected in (0, 1):
        c2, raw2, _, matrix = next(row for row in endings[selected] if np.linalg.norm(row[3]) > .1)
        record = (c1, raw1, (0, 0, selected, 0), c2, raw2)
        result = exact.compare_map(matrix, exact.ideal_map(law, mode, record))
        assert all(type(value) is F for value in result.values())
        assert 0 < result['distance_upper'] < F(1, 10**10)
        assert result['padding_frobenius_squared'] >= 0
    c1, raw1, _, _, matrix, _ = next(row for row in tree if row[3])
    result = exact.compare_map(matrix, exact.ideal_map(law, mode, (c1, raw1, (), None, None)))
    assert result['ideal_frobenius_squared'] == 0
    assert result['distance_upper'] < F(1, 10**20)


@pytest.mark.parametrize('hermitian', (False, True))
def test_high_precision_svd_on_entangled_reference_and_general_operator(hermitian):
    # The occupied output rows include a padding row. The independent SVD
    # computes the complete output on this invariant finite support.
    reported = np.zeros(exact.SHAPE, complex)
    reported[0] = (.25+.125j, -.5j)
    reported[1] = (.75j, .5-.25j)
    reported[511] = (.125, .25j)
    ideal = {(0, 0): F(1, 3), (1, 1): F(2, 7)}
    result = exact.compare_map(reported, ideal)
    ctx = mpmath.mp.clone(); ctx.dps = 100
    n = ctx.matrix([[_number(z, ctx) for z in reported[row]] for row in (0, 1, 511)])
    m = ctx.matrix([[ctx.sqrt(_mp(F(1, 3), ctx)), 0],
                    [0, ctx.sqrt(_mp(F(2, 7), ctx))], [0, 0]])
    eps = ctx.sqrt(sum(abs(z)**2 for z in n-m))
    assert _mp(result['frobenius_error_lower'], ctx) <= eps
    assert eps <= _mp(result['frobenius_error_upper'], ctx)
    assert _mp(result['frobenius_error_upper'], ctx)-eps < ctx.mpf('1e-50')
    # Expand by an untouched two-dimensional reference, including operators
    # that are neither positive nor Hermitian. This tests the cb trace norm
    # inequality rather than only a basis input or maximally mixed state.
    nr, mr = ctx.zeros(6, 4), ctx.zeros(6, 4)
    for i, j, r in product(range(3), range(2), range(2)):
        nr[2*i+r, 2*j+r], mr[2*i+r, 2*j+r] = n[i, j], m[i, j]
    if hermitian:
        v = ctx.matrix([1, ctx.j, 2, -ctx.j])/ctx.sqrt(7)
        x = v*v.H
    else:
        x = ctx.matrix([[1, ctx.j, 0, -2], [0, 1, 2*ctx.j, 0],
                        [3, 0, -ctx.j, 1], [0, 0, 1, 0]])
    input_norm = sum(ctx.svd(x, compute_uv=False))
    output = nr*x*nr.H-mr*x*mr.H
    measured = sum(ctx.svd(output, compute_uv=False))/2/input_norm
    assert measured <= _mp(result['distance_upper'], ctx)


def test_relative_phase_and_wrong_cptp_output_cannot_pass_by_matching_mass():
    ideal = {(0, 0): F(1), (1, 1): F(1)}
    correct = np.zeros(exact.SHAPE, complex)
    correct[0, 0] = correct[1, 1] = 1
    assert exact.compare_map(correct, ideal)['distance_upper'] == 0
    phase = correct.copy(); phase[1, 1] = -1
    flip = correct[:, ::-1].copy()
    for wrong in (phase, flip):
        # Both maps are isometries, so all-input completeness by itself is
        # insufficient to certify the intended channel.
        np.testing.assert_array_equal(wrong.conj().T@wrong, np.eye(2))
        result = exact.compare_map(wrong, ideal)
        assert result['reported_frobenius_squared'] == result['ideal_frobenius_squared'] == 2
        assert result['frobenius_error_lower'] >= 2
        assert result['distance_upper'] >= 1


def test_padding_and_abort_mass_never_projected_or_clipped():
    reported = np.zeros(exact.SHAPE, complex)
    reported[511, 0] = 4j
    assert 511 not in source.CODES
    result = exact.compare_map(reported, {})
    assert result['reported_frobenius_squared'] == result['padding_frobenius_squared'] == 16
    assert result['distance_upper'] == 8


def test_native_certificate_rejects_subclass_that_conceals_nonzero_entries():
    class Concealed(np.ndarray):
        def copy(self, *args, **kwargs):
            return np.zeros(self.shape, dtype=self.dtype)

    reported = np.zeros(exact.SHAPE).view(Concealed)
    reported[0, 0] = 1
    # Before the concrete-type guard, the overridden copy yielded a false
    # zero distance, although this map has half-diamond distance exactly 1/2.
    assert exact.compare_map(np.asarray(reported), {})['distance_upper'] == F(1, 2)
    with pytest.raises(ValueError, match='base ndarray'):
        exact.compare_map(reported, {})


def test_native_certificate_rejects_subclass_before_metadata_or_copy_dispatch():
    class Hostile(np.ndarray):
        def __getattribute__(self, name):
            if name in ('shape', 'dtype', 'ndim', 'copy'):
                raise AssertionError('subclass array dispatch reached')
            return super().__getattribute__(name)

    with pytest.raises(ValueError, match='base ndarray'):
        exact.compare_map(np.zeros(exact.SHAPE).view(Hostile), {})


def test_integer_components_remain_exact_above_binary64_precision():
    reported = np.zeros(exact.SHAPE, np.int64)
    reported[0, 0] = 2**53+1
    result = exact.compare_map(reported, {})
    assert result['reported_frobenius_squared'] == F((2**53+1)**2)
    assert result['distance_upper'] == F((2**53+1)**2, 2)


@pytest.mark.parametrize('dtype', (np.float16, np.float32, np.float64, np.complex64, np.complex128))
def test_binary_components_are_compared_at_original_precision(dtype):
    reported = np.zeros(exact.SHAPE, dtype=dtype)
    reported[0, 0] = .1
    result = exact.compare_map(reported, {})
    expected = F(float(reported[0, 0].real))**2
    assert result['reported_frobenius_squared'] == expected


@pytest.mark.parametrize('bad', [None, [], True, np.zeros((2, 512)),
                               np.zeros(exact.SHAPE, bool), np.zeros(exact.SHAPE, object),
                               np.ma.array(np.zeros(exact.SHAPE)),
                               np.full(exact.SHAPE, np.nan), np.full(exact.SHAPE, np.inf)])
def test_original_invalid_arrays_are_rejected(bad):
    with pytest.raises(ValueError):
        exact.compare_map(bad, {})


def test_extended_precision_is_not_silently_rounded():
    for dtype, maximum in ((np.longdouble, 8), (np.clongdouble, 16)):
        if np.dtype(dtype).itemsize > maximum:
            with pytest.raises(ValueError):
                exact.compare_map(np.zeros(exact.SHAPE, dtype), {})


@pytest.mark.parametrize('bad', [{(0, 0): 1}, {(0, 0): .5}, {(0, 0): True},
                               {(0, 0): F(-1)}, {(512, 0): F(1)}, {(0, 2): F(1)},
                               {(True, 0): F(1)}, {(0, np.int64(0)): F(1)}, [], None])
def test_invalid_squared_target_envelopes_are_rejected(bad):
    with pytest.raises(ValueError):
        exact.compare_map(np.zeros(exact.SHAPE), bad)


@pytest.mark.parametrize('warm', (False, True))
def test_labels_and_records_are_checked_before_cached_computation(warm):
    class Aliased(str):
        pass
    exact._tree.cache_clear(); exact._menu.cache_clear()
    if warm:
        exact.ideal_tree('C', MODES[0])
        exact.ideal_map('C', MODES[0], (0, 3, (), None, None))
    for law, mode in ((Aliased('C'), MODES[0]), ('C', Aliased(MODES[0])),
                      ('C', 'full_carrier'), (True, MODES[0])):
        with pytest.raises(ValueError):
            exact.ideal_tree(law, mode)
        with pytest.raises(ValueError):
            exact.completeness(law, mode)
    for record in ((True, 3, (), None, None), (0, 3., (), None, None),
                   (0, 0, (0, 0, True, 0), 0, 0), (0, 0, (0, 0, 0), 0, 0),
                   (0, 3, (0, 0, 0, 0), 0, 0), (0, 0, (), None, None),
                   (0, 0, [0, 0, 0, 0], 0, 0), (0, 0, (0, 0, 0, 0), False, 0)):
        with pytest.raises(ValueError):
            exact.ideal_map('C', MODES[0], record)
