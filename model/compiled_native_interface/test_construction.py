"""Independent source-word, complete-carrier and numerical-bound controls."""
from fractions import Fraction as F
import math
import warnings

import mpmath
import numpy as np
import pytest

from ..native_interface import construction as interface
from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from . import construction as c


def test_short_source_word_matches_generic_compiler_on_every_seam_column():
    independent = np.eye(8)
    independent[np.ix_((2, 1), (2, 1))] = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    generic = compiler.compile_real_orthogonal(independent)
    assert len(c.word()) == 34 < len(generic) == 1046
    np.testing.assert_allclose(c.local_matrix(), compiler.execute(generic, 3), atol=3e-13, rtol=0)
    np.testing.assert_allclose(c.local_matrix(), independent, atol=2e-14, rtol=0)
    assert np.any(c.local_matrix()[independent == 0] != 0)  # original defects remain


def test_exact_complete_unitarity_and_x_y_effects():
    result = c.exact_effect_audit()
    assert result['complete_unitarity_entries_checked'] == 192
    assert result['complete_pointer_effect_entries_checked'] == 192


def test_full_78_carrier_and_all_padding_columns_use_correct_spectator_wiring():
    full = c.full_matrix()
    local = c.local_matrix()
    np.testing.assert_allclose(full, compiler.execute(c.word(), 9), atol=2e-14, rtol=0)
    for spectator in range(64):
        np.testing.assert_array_equal(full[8*spectator:8*spectator+8, 8*spectator:8*spectator+8], local)
    codes = tuple(circuit.encode(state) for state in native.states())
    target = np.kron(np.eye(64), c.ideal_local_matrix())
    outside = sorted(set(range(512))-set(codes))
    assert not np.any(target[np.ix_(outside, codes)])
    # The padding extension is the same local operation, not falsely identity.
    assert target[8*63+1, 8*63+2] != 0
    embedding = c.embedding()
    np.testing.assert_allclose(target@embedding, embedding@interface.mixer(), atol=2e-16, rtol=0)
    assert set(embedding.sum(axis=0)) == {1}


@pytest.mark.parametrize('kind', c.KINDS)
def test_directed_binary_bound_encloses_independent_high_precision_frobenius(kind):
    ctx = mpmath.mp.clone(); ctx.dps = 95
    array = c.local_matrix(kind)
    targets = c.target(kind)
    def mpq(value):
        value = F(value)
        return ctx.mpf(value.numerator)/value.denominator
    squared = ctx.mpf(0)
    for i in range(8):
        for j in range(8):
            real, imag = targets.get((i, j), (c.ZERO, c.ZERO))
            wanted = ctx.mpc(mpq(real[0])+ctx.sqrt(2)*mpq(real[1]),
                             mpq(imag[0])+ctx.sqrt(2)*mpq(imag[1]))
            original = ctx.mpc(mpq(float(array[i, j].real)), mpq(float(array[i, j].imag)))
            squared += abs(original-wanted)**2
    bound = c.certify(kind=kind)
    assert mpq(bound['frobenius_error_lower']) <= ctx.sqrt(squared) <= mpq(bound['frobenius_error_upper'])
    assert bound['operator_error_upper'] < c.ERROR_GATE


def test_wrong_phase_or_swapped_orientation_cannot_pass_by_matching_unitarity():
    wrong_phase = -c.local_matrix()
    wrong_order = c.local_matrix().copy()
    wrong_order[[1, 2]] = wrong_order[[2, 1]]
    for wrong in (wrong_phase, wrong_order):
        np.testing.assert_allclose(wrong.conj().T@wrong, np.eye(8), atol=2e-14, rtol=0)
        assert c.certify(wrong)['operator_error_upper'] > 1


def test_off_support_defect_and_large_integer_entries_survive_exact_conversion():
    wrong = c.local_matrix().copy(); wrong[7, 0] = 2**-20
    assert c.certify(wrong)['operator_error_upper'] >= F(1, 2**20)
    large = np.zeros((8, 8), dtype='>u8'); large[0, 0] = 2**53+1
    assert c.certify(large)['operator_error_upper'] > F(2**53)


def test_two_calls_retain_binary_channel_norm_gain():
    # The exact algebra saturates for two amplified identity maps. Use the
    # scalar bound algebra independently of any native target gate.
    epsilon = F(1, 4)
    delta = epsilon*(2+epsilon)/2
    gain = (1+epsilon)**2
    assert delta*(1+gain) == ((1+epsilon)**4-1)/2
    assert delta*(1+gain) > 2*delta
    bound = c.certify()
    assert bound['two_call_half_diamond_upper'] == bound['half_diamond_upper']*(1+bound['channel_norm_upper'])


def test_source_mixer_preserves_reference_and_supports_complete_native_pointer():
    # Full 512-by-reference source application. No code compression is used.
    embedding = c.embedding()
    amplitudes = np.zeros((66, 2), complex)
    left, right = interface.phase_pairs()[0]
    amplitudes[left, 0], amplitudes[right, 1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    full_input = embedding@amplitudes
    actual = c.full_matrix()@full_input
    target = embedding@interface.mixer()@amplitudes
    assert np.linalg.norm(actual-target, 'fro') <= float(c.certify()['operator_error_upper'])+1e-16
    np.testing.assert_allclose(actual.conj().T@actual, np.eye(2)/2, atol=2e-14, rtol=0)
    for law in ('C', 'E'):
        _, pointer = interface.herald_and_pointer(law)
        # Preserve every pointer Kraus branch and a padding extension of the
        # complete instrument; projection failure is a separate retained row.
        columns = embedding.T@actual
        outputs = tuple(k@columns for group in pointer.kraus for k in group)
        padding = actual-embedding@columns
        assert abs(sum(np.vdot(v, v).real for v in outputs)+np.vdot(padding, padding).real-1) < 3e-14


def test_one_source_phase_pulse_distinguishes_conjugates_that_real_pointer_cannot():
    pointer = np.eye(8)/2; pointer[2, 2] = 1; pointer[1, 1] = 0
    real_effect = c.local_matrix().conj().T@pointer@c.local_matrix()
    complex_effect = c.local_matrix('y_analyzer').conj().T@pointer@c.local_matrix('y_analyzer')
    for sign in (-1, 1):
        state = np.zeros(8, complex); state[2] = 1/math.sqrt(2); state[1] = sign*1j/math.sqrt(2)
        assert abs(np.vdot(state, real_effect@state).real-.5) < 2e-14
        assert abs(np.vdot(state, complex_effect@state).real-(1+sign)/2) < 2e-14
    assert len(c.word('phase')) == 1


@pytest.mark.parametrize('bad', [[], True, np.zeros((66, 66)), np.zeros((8, 8), bool),
                               np.ones((8, 8), object), np.full((8, 8), np.nan)])
def test_invalid_original_arrays_rejected(bad):
    with pytest.raises(ValueError):
        c.certify(bad)


def test_subclasses_rejected_before_metadata_and_cache_cannot_be_poisoned():
    class Hostile(np.ndarray):
        def __getattribute__(self, name):
            if name in ('shape', 'dtype', 'ndim', 'copy'):
                raise AssertionError('array subclass dispatch')
            return super().__getattribute__(name)
    with pytest.raises(ValueError, match='base ndarray'):
        c.certify(np.eye(8).view(Hostile))
    first = c.local_matrix()
    with pytest.raises(ValueError):
        first[0, 0] = 10
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        first.shape = (64,)
    assert c.local_matrix().shape == (8, 8)
    tape = c.word(); tape[0]['angle'] = 0
    assert c.word()[0]['angle'] != 0


@pytest.mark.parametrize('bad', [True, 'unknown', None, 1])
def test_operation_labels_validated_before_cache_lookup(bad):
    for function in (c.word, c.local_matrix, c.full_matrix, c.target, c.resources):
        with pytest.raises(ValueError):
            function(bad)
