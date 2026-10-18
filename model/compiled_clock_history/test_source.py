"""Independent target, full-padding and hostile-input source controls."""
from fractions import Fraction as F
from itertools import product
import json
import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest

from ..quantum_observer_correspondence import construction as ideal
from . import source as s


def _number(value):
    return mp.mpf(value.numerator)/value.denominator


def _float_target(kind, gamma=1, reference=False):
    shape = ((16 if reference else 8,) if kind == 'preparation'
             else (16, 8) if kind == 'readout' else (8, 8))
    result = np.zeros(shape, complex)
    for key, (a, b) in s.target(kind, gamma=gamma, reference=reference).items():
        result[key] = float(a)+float(b)*np.sqrt(2)
    return result


def test_full_chart_is_fixed_injective_and_refuses_all446_padding_addresses():
    assert s.STATES == ideal.STATES
    assert len(set(s.codes())) == 66
    for state, code in zip(ideal.STATES, s.codes(), strict=True):
        assert s.decode(code) == state and s.encode(state) == code
        assert tuple((code >> bit) & 1 for bit in range(3)) == tuple(state[p] for p in (0, 1, 5))
    refused = 0
    for code in set(range(512))-set(s.codes()):
        with pytest.raises(ValueError):
            s.decode(code)
        refused += 1
    assert refused == 446
    assert s.decode(3) == ideal.STATES[ideal.address((0, 1))]
    assert s.decode(5) == ideal.STATES[ideal.address((0, 5))]


@pytest.mark.parametrize('gamma', (0, 1))
def test_every_exact_native_apparatus_column_matches_the_old66_state_target(gamma):
    target = s.target('apparatus', gamma=gamma)
    for i, rows in ideal.apparatus(gamma).items():
        code = s.codes()[i]
        observed = {(row+8*(code//8)): a for (row, column), (a, b) in target.items()
                    if column == code % 8 and a and b == 0}
        expected = {s.codes()[j]: value for j, value in rows}
        assert observed == expected
    # Exact complete local unitarity includes addresses with no native meaning.
    for i, j in product(range(8), repeat=2):
        gram = sum(target.get((k, i), (F(0), F(0)))[0]
                   *target.get((k, j), (F(0), F(0)))[0] for k in range(8))
        assert gram == int(i == j)


@pytest.mark.parametrize('kind,reference,gamma', (
    ('preparation', False, 1), ('preparation', True, 1),
    ('apparatus', False, 0), ('apparatus', False, 1), ('readout', False, 1)))
def test_directed_original_binary_certificate_encloses_independent_high_precision_error(kind, reference, gamma):
    raw = (s.prepared(reference) if kind == 'preparation' else s.readout_isometry()
           if kind == 'readout' else s.local_apparatus(gamma))
    result = s.certify(kind, raw, gamma=gamma, reference=reference)
    wanted = s.target(kind, gamma=gamma, reference=reference)
    with mp.workdps(100):
        squared = mp.mpf(0)
        for index in np.ndindex(raw.shape):
            a, b = wanted.get(index, (F(0), F(0)))
            delta = mp.mpc(float(raw[index].real), float(raw[index].imag))-_number(a)-_number(b)*mp.sqrt(2)
            squared += abs(delta)**2
        distance = mp.sqrt(squared)
        assert _number(result['frobenius_error_lower']) <= distance <= _number(result['frobenius_error_upper'])
    assert result['half_diamond_upper'] < F(1, 10**12)
    assert result['channel_norm_upper'] >= 1
    assert result['all_original_components_checked'] == raw.size
    assert result['fixed_blank_preparation_only'] == (kind == 'preparation')
    assert result['arbitrary_input_and_reference'] == (kind != 'preparation')


@pytest.mark.parametrize('reference', (False, True))
def test_actual_blank_preparation_and_reference_column_lift(reference):
    blank = np.zeros(16 if reference else 8, complex)
    blank[0] = 1
    direct = s.compiler.execute(s.word('preparation', reference=reference), 4 if reference else 3, blank)
    np.testing.assert_array_equal(direct, s.prepared(reference))
    lifted = s.preparation(reference)
    assert lifted.shape == (512, 2 if reference else 1)
    np.testing.assert_array_equal(lifted[:8], direct.reshape(-1, 8).T)
    assert np.count_nonzero(lifted[8:]) == 0
    np.testing.assert_allclose(direct, _float_target('preparation', reference=reference), atol=3e-15, rtol=0)


def test_actual_full512_lift_and_source_pointer_keep_every_padding_input():
    local, full = s.local_apparatus(), s.apparatus()
    for label in range(64):
        np.testing.assert_array_equal(full[8*label:8*label+8, 8*label:8*label+8], local)
    assert np.count_nonzero(full) == 64*np.count_nonzero(local)
    pointer = s.readout_isometry()
    np.testing.assert_array_equal(np.vstack(s.local_readout()), pointer)
    for column in range(8):
        expected = np.zeros(16)
        expected[column+8*((column >> 1) & 1)] = 1
        np.testing.assert_allclose(pointer[:, column], expected, atol=2e-15, rtol=0)


def test_wrong_unitary_phase_and_off_support_padding_entry_are_charged():
    wrong = s.local_apparatus().copy()
    wrong[7] *= -1  # The invalid local111 column is nevertheless in the certified map.
    assert np.linalg.norm(wrong.conj().T@wrong-np.eye(8), 2) < 1e-13
    assert s.certify('apparatus', wrong)['frobenius_error_lower'] > F(19, 10)
    wrong = s.local_apparatus().copy()
    wrong[0, 7] += .125j
    assert s.certify('apparatus', wrong)['frobenius_error_lower'] > F(1, 9)
    wrong = s.readout_isometry().copy()
    wrong[[0, 8]] = wrong[[8, 0]]  # Move a pointer label, preserve isometry.
    assert np.linalg.norm(wrong.conj().T@wrong-np.eye(8), 2) < 1e-13
    assert s.certify('readout', wrong)['frobenius_error_lower'] > 1


def test_public_arrays_words_targets_and_certificates_cannot_poison_cached_execution():
    raw = s.local_apparatus()
    expected = raw.tobytes()
    with pytest.raises(ValueError):
        raw.flags.writeable = True
    raw.shape = (64,)
    assert s.local_apparatus().shape == (8, 8) and s.local_apparatus().tobytes() == expected
    tape = s.word()
    tape[0]['angle'], tape[0]['targets'][0] = 0., 2
    assert s.word()[0]['angle'] != 0 and s.word()[0]['targets'][0] == 0
    target = s.target()
    target[0, 0][0]._numerator = 99
    assert s.target()[0, 0][0] == 1
    certificate = s.certify()
    saved = {key: str(value) for key, value in certificate.items()}
    for value in certificate.values():
        if type(value) is F:
            value._numerator = 0
    assert {key: str(value) for key, value in s.certify().items()} == saved


class ArrayAlias(np.ndarray):
    def copy(self, *args, **kwargs):
        raise AssertionError('hostile array method must not run')


@pytest.mark.parametrize('bad', (np.zeros((8, 8), bool), np.zeros((8, 8), object),
    np.full((8, 8), complex(np.inf, 0)), np.zeros((8, 7)), np.zeros((8, 8)).view(ArrayAlias),
    np.ma.array(np.zeros((8, 8))), [[0]*8]*8))
def test_certificate_refuses_original_bad_array_controls(bad):
    with pytest.raises(ValueError):
        s.certify('apparatus', bad)


@pytest.mark.parametrize('kwargs', ({'kind': True}, {'kind': 'unknown'}, {'gamma': True},
    {'gamma': 2}, {'gamma': 1.0}, {'reference': 1}))
def test_all_original_controls_are_checked_before_cache_lookup(kwargs):
    s.certify()
    for call in (s.word, s.resources, s.target, s.certify):
        with pytest.raises(ValueError):
            call(**kwargs)


@pytest.mark.parametrize('value', (True, -1, 512, 1.0))
def test_decode_rejects_original_invalid_addresses(value):
    with pytest.raises(ValueError):
        s.decode(value)


def test_resources_charge_zero_phase_requests_and_separate_ancillas():
    assert s.resources('preparation')['native_pulses'] == 24
    assert s.resources('preparation', reference=True)['native_pulses'] == 27
    assert s.resources('apparatus', gamma=0)['native_pulses'] == 24
    assert s.resources('apparatus', gamma=0)['zero_angle_requests'] == 2
    assert s.resources('apparatus')['proper_code_service_events'] == 48
    assert s.resources('readout')['native_pulses'] == 3
    assert s.resources('readout')['blank_pointer_qubits'] == 1
    assert s.resources('preparation', reference=True)['blank_reference_qubits'] == 1
    assert s.resources('readout')['pointer_reads'] == 1
    assert len(s.bindings()) == len(s.SOURCE_PATHS)


def test_static_source_bindings_cover_a_fresh_import_and_execution_profile():
    # This subprocess never imports pytest or the independent ideal oracle.
    # The production binding is static, so those different callers cannot
    # accidentally produce different receipt inventories.
    script = '''
from pathlib import Path
import json, sys
root=Path.cwd().resolve()
sys.path[:0]=[str(root),str(root/'code'),str(root/'tools')]
from model.compiled_clock_history import source as s
assert s.ROOT==root
for kind, ref, gamma in [('preparation',False,1),('preparation',True,1),
                         ('apparatus',False,0),('apparatus',False,1),('readout',False,1)]:
    s.certify(kind,reference=ref,gamma=gamma)
    s.resources(kind,reference=ref,gamma=gamma)
paths=set()
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        path=Path(name).resolve()
        if path.suffix=='.py' and path.is_relative_to(s.ROOT):
            relative=path.relative_to(s.ROOT)
            if relative.parts[0] in ('model','code','tools','evidence'):
                paths.add(relative.as_posix())
print(json.dumps(sorted(paths)))
'''
    run = subprocess.run([sys.executable, '-W', 'error', '-c', script],
                         cwd=s.ROOT, check=True, capture_output=True, text=True)
    assert set(json.loads(run.stdout)) <= set(s.SOURCE_PATHS)
