"""Independent effect, scheduler, composition and receipt rejection checks."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import subprocess
import sys

import pytest
from sympy.polys.domains import QQ_I

from ..observer import quantum_mixing
from . import history as h, receipt, source


def test_whole_producer_imports_have_static_source_custody():
    root = Path(__file__).resolve().parents[2]
    script = '''
import sys,json
from pathlib import Path
root=Path.cwd().resolve()
sys.path[:0]=[str(root),str(root/'code'),str(root/'tools')]
from model.compiled_clock_history import receipt,history
history.record_menu()
paths=[]
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if not name: continue
    try: path=Path(name).resolve().relative_to(root).as_posix()
    except ValueError: continue
    if path.endswith('.py') and path.startswith(('model/','code/','tools/')):
        paths.append(path)
print(json.dumps(sorted(set(paths))))
'''
    observed = json.loads(subprocess.check_output([sys.executable, '-c', script], cwd=root, text=True))
    retained = json.loads((root/'model/compiled_clock_history/results.json').read_text(encoding='utf8'))
    assert set(observed) <= set(retained['source_sha256'])


def _schedules(cap):
    """Enumerate geometric successes and retained timeouts independently."""
    single = [(F(2, 5)*F(3, 5)**n, n+1, True) for n in range(cap)]
    single.append((F(3, 5)**cap, cap, False))
    rows = []
    for weight, attempts, success in single:
        if not success:
            rows.append((weight, attempts, 0))
        else:
            rows.extend((weight*second, attempts+later, 1)
                        for second, later, _ in single)
    assert sum(row[0] for row in rows) == 1
    return tuple(rows)


def test_complete_native_effects_on_all4356_matrix_units_and_padding():
    # The success probability is derived from source incidence, not from a
    # scalar inserted by the receipt or from diagonal-only state checks.
    seams = quantum_mixing.seams()
    incident = sum(bool(set(seam) & {0, 1, 5}) for seam in seams)
    assert len(seams) == 30 and incident == 12
    remote, marks = h._counts()
    for i, j in product(source.codes(), repeat=2):
        remote_trace = sum(count for (a, b), count in remote.get((i, j), ()) if a == b)
        marked_trace = sum(count for table in marks.values()
                           for (a, b), count in table.get((i, j), ()) if a == b)
        assert remote_trace == 2*(len(seams)-incident)*int(i == j)
        assert marked_trace == 2*incident*int(i == j)
        for seam in ((0, 1), (0, 5)):
            feedback_trace = sum(count for table in h._counts(seam)[1].values()
                                 for (a, b), count in table.get((i, j), ()) if a == b)
            assert feedback_trace == 2*int(i == j)

    # Both cross-code orientations vanish; all446 padding populations have
    # effect one in the separately retained complement alternative.
    good = source.codes()[0]
    for bad in set(range(512))-set(source.codes()):
        state = {(bad, bad, 0, 0): QQ_I.one,
                 (good, bad, 0, 1): QQ_I(1, 2),
                 (bad, good, 1, 0): QQ_I(1, -2)}
        completed, waiting, aborts = h.interval(state, 1)
        assert not waiting and not any(completed.values())
        assert aborts == ({(bad, bad, 0, 0): QQ_I.one},)


@pytest.mark.parametrize('cap', (1, 2, 3))
@pytest.mark.parametrize('reference', (False, True))
def test_resource_expectations_are_full_scheduler_sums(cap, reference):
    rows = _schedules(cap)
    result = receipt.costs(cap=cap, reference=reference)
    preparation = 27 if reference else 24
    assert result['source_pulses_maximum'] == preparation+24+3
    assert result['proper_code_loads_maximum'] == preparation+24+3
    assert result['proper_code_exports_maximum'] == preparation+24+3
    assert result['scheduled_native_attempts_maximum'] == max(row[1] for row in rows)
    assert result['feedback_native_attempts_maximum'] == max(row[2] for row in rows)
    assert F(result['scheduled_native_attempts_ideal_expectation']) == sum(
        weight*attempts for weight, attempts, feedback in rows)
    assert F(result['feedback_native_attempts_ideal_expectation']) == sum(
        weight*feedback for weight, attempts, feedback in rows)
    assert F(result['source_pulses_ideal_expectation']) == sum(
        weight*(preparation+24*feedback+3) for weight, attempts, feedback in rows)
    assert result['peak_code_qubits_including_pointer_and_reference'] == 10+int(reference)
    assert result['complement_checks_before_native_attempts_maximum'] == 2*cap+1
    assert result['ideal_expectations_are_not_normalized_raw_masses'] is True


@pytest.mark.parametrize('cap', (1, 2, 3))
def test_composition_bound_is_saturated_by_exact_positive_scalar_errors(monkeypatch, cap):
    # Scale a normalized preparation, its successful-branch apparatus and
    # the final detector isometry. Every output block grows positively, so
    # half of the excess total mass IS the complete half trace-norm error.
    # This catches missing gains, including terms quadratic in tiny errors.
    scales = {kind: 1+F(n, 10**14) for n, kind in enumerate(
        ('preparation', 'apparatus', 'readout'), start=1)}

    def certificate(kind, **kwargs):
        scale = scales[kind]
        return {'operator_error_upper': scale-1,
                'half_diamond_upper': (scale*scale-1)/2,
                'channel_norm_upper': scale*scale}

    monkeypatch.setattr(source, 'certify', certificate)
    raw_mass = sum(weight*scales['preparation']**2*scales['readout']**2
                   *scales['apparatus']**(2*feedback)
                   for weight, attempts, feedback in _schedules(cap))
    assert receipt.bound(cap=cap) == (raw_mass-1)/2
    assert 0 < receipt.bound(cap=cap) < receipt.ERROR


@pytest.mark.parametrize('cap', (1, 2, 3))
def test_declared_record_buffer_encodes_every_word_injectively(cap):
    records = h.record_menu(cap)
    scheduled = sorted({mark for word, _, _ in records
                        for kind, mark in word if kind == 'scheduled'})
    statuses = sorted({status for _, status, _ in records})
    # The saved first mark selects the feedback seam and its three-outcome
    # menu. This field does not encode the five-label union of both menus.
    feedback = {mark: sorted({word[1][1] for word, _, _ in records
                             if len(word) > 1 and word[0][1] == mark})
                for mark in scheduled}
    assert len(scheduled) == 15 and max(map(len, feedback.values())) == 3
    status_bits = (len(statuses)-1).bit_length()
    encoded = set()
    for word, status, bit in records:
        first = scheduled.index(word[0][1]) if word else 15
        middle = feedback[word[0][1]].index(word[1][1]) if len(word) > 1 else 3
        last = scheduled.index(word[2][1]) if len(word) > 2 else 15
        encoded.add(first+16*middle+64*last+1024*statuses.index(status)
                    +(1024 << status_bits)*bit)
    assert len(encoded) == len(records)
    assert max(encoded) < 2**receipt.costs(cap=cap)['fixed_record_buffer_bits']


@pytest.fixture(scope='module')
def exact_reference_target():
    # A rational inherited execution is enough for verifier hostility tests;
    # this fixture does not duplicate an expensive emitted-array history.
    return h.ideal_history(gamma=1, cap=1, reference=True)


def _complete(target):
    return {record: dict(target.get(record, {})) for record in h.record_menu(1)}


def test_receipt_accepts_a_complete_exact_target(exact_reference_target):
    report = receipt.audit_output(_complete(exact_reference_target), exact_reference_target,
                                  cap=1, reference=True)
    assert report['complete_output_records'] == 1564
    assert report['inherited_nonabort_records'] == 1442
    assert report['explicit_code_complement_records'] == 122


class _IndexAlias(int):
    def __hash__(self):
        return hash(0)

    def __eq__(self, other):
        return other == 0


def test_transposed_scalar_cannot_repair_itself_and_rewrite_a_checked_block():
    menu = h.record_menu(1)
    first, second = sorted(record for record in menu if record[1] == 'complete')[:2]
    earlier = {(0, 0, 0, 0): QQ_I(F(1, 2))}
    later = {key: QQ_I(F(1, 4)) for key in (
        (0, 0, 0, 0), (0, 1, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0))}
    actual = {first: earlier, second: later}
    actual.update((record, {}) for record in menu if record not in actual)
    target = {record: dict(state) for record, state in actual.items()
              if not record[1].startswith('code_abort')}
    victim = later[1, 0, 0, 0]
    saved = victim.x
    # Keep the independent target's scalar separate from the crafted object.
    target[second][1, 0, 0, 0] = QQ_I(F(1, 4))

    class Component:
        calls = 0

        def __eq__(self, other):
            self.calls += 1
            victim.x = saved
            value = next(iter(earlier.values()))
            earlier.clear()
            earlier[_IndexAlias(511), _IndexAlias(511), 0, 0] = value
            return saved == other

    component = Component()
    victim.x = component
    with pytest.raises(ValueError):
        receipt.audit_output(actual, target, cap=1)
    assert component.calls == 0


def test_target_index_alias_cannot_make_orthogonal_states_compare_equal():
    record = next(record for record in h.record_menu(1) if record[1] == 'complete')
    actual = {key: {} for key in h.record_menu(1)}
    actual[record] = {(0, 0, 0, 0): QQ_I.one}
    target = {record: {(_IndexAlias(511), _IndexAlias(511), 0, 0): QQ_I.one}}
    for compare in (lambda: h.distance_upper(actual, target),
                    lambda: receipt.audit_output(actual, target, cap=1)):
        with pytest.raises(ValueError, match='original exact'):
            compare()


def test_original_record_types_are_checked_before_equality_dispatch():
    class Status(str):
        __hash__ = str.__hash__

        def __eq__(self, other):
            raise AssertionError('unvalidated status equality must never run')

    actual = {record: {} for record in h.record_menu(1)}
    del actual[(), 'timeout_first', 0]
    actual[(), Status('timeout_first'), 0] = {}
    with pytest.raises(ValueError, match='original tuple/string/integer'):
        h.validate_history(actual, cap=1)


def test_exact_components_are_type_checked_before_attribute_dispatch():
    class Component:
        def __getattribute__(self, name):
            raise AssertionError('unvalidated component attribute must never run')

    with pytest.raises(ValueError, match='canonical exact rational'):
        h.fraction(Component())


def test_validated_history_detaches_all_mutable_exact_scalar_slots():
    actual = {record: {} for record in h.record_menu(1)}
    record = ((), 'timeout_first', 0)
    actual[record] = {(0, 0, 0, 0): QQ_I(F(1, 2)),
                      (0, 1, 0, 0): QQ_I(F(1, 4), F(1, 8)),
                      (1, 0, 0, 0): QQ_I(F(1, 4), F(-1, 8)),
                      (1, 1, 0, 0): QQ_I(F(1, 2))}
    saved = h.validate_history(actual, cap=1)
    assert saved == actual and saved is not actual and saved[record] is not actual[record]
    for key, value in actual[record].items():
        copied = saved[record][key]
        assert copied is not value and copied.x is not value.x and copied.y is not value.y
    actual[record][0, 1, 0, 0].x.numerator = 9
    actual[record][0, 1, 0, 0].y.numerator = 9
    assert saved[record][0, 1, 0, 0] == QQ_I(F(1, 4), F(1, 8))


@pytest.mark.parametrize('mutation', ('extra', 'missing', 'abort'))
def test_target_menu_cannot_change_counts_by_empty_records(exact_reference_target, mutation):
    target = dict(exact_reference_target)
    if mutation == 'extra':
        target[(), 'unrelated', 0] = {}
    elif mutation == 'missing':
        del target[next(record for record, state in target.items() if not state)]
    else:
        target[(), 'code_abort_first_0', 0] = {}
    with pytest.raises(ValueError, match='complete.*record menu'):
        receipt.audit_output(_complete(exact_reference_target), target, cap=1, reference=True)


def test_receipt_refuses_boolean_record_alias(exact_reference_target):
    actual = _complete(exact_reference_target)
    record = next(record for record in actual if record[2] == 1)
    state = actual.pop(record)
    actual[record[:2]+(True,)] = state
    assert set(actual) == h.record_menu(1)  # Python equality alone is insufficient.
    with pytest.raises(ValueError, match='original tuple/string/integer'):
        receipt.audit_output(actual, exact_reference_target, cap=1, reference=True)


@pytest.mark.parametrize('abort', (False, True))
def test_receipt_rejects_a_missing_zero_mass_record(exact_reference_target, abort):
    actual = _complete(exact_reference_target)
    missing = next(record for record, state in actual.items()
                   if not state and record[1].startswith('code_abort') == abort)
    del actual[missing]
    # A distance-only checker misses deletion of an exactly zero block.
    assert h.distance_upper(actual, exact_reference_target) == 0
    with pytest.raises(ValueError, match='complete.*record menu'):
        receipt.audit_output(actual, exact_reference_target, cap=1, reference=True)


def test_receipt_rejects_reference_dephasing_with_unchanged_complete_record_law(
        exact_reference_target):
    actual = _complete(exact_reference_target)
    dephased = {record: {key: value for key, value in state.items() if key[2] == key[3]}
                for record, state in actual.items()}
    assert {record: h.trace(state) for record, state in actual.items()} == {
        record: h.trace(state) for record, state in dephased.items()}
    assert receipt.joint(actual) == receipt.joint(dephased)
    assert h.distance_upper(dephased, exact_reference_target) > F(1, 100)
    with pytest.raises(ArithmeticError, match='independent target'):
        receipt.audit_output(dephased, exact_reference_target, cap=1, reference=True)


def test_receipt_rejects_padding_error_with_unchanged_record_masses(exact_reference_target):
    actual = _complete(exact_reference_target)
    record = max(actual, key=lambda word: h.trace(actual[word]))
    # This is a valid local replacement channel: reset work to padding while
    # retaining the complete reduced reference and every classical mass.
    bad = next(code for code in range(512) if code not in h.CODE)
    reset = {}
    for (i, j, r, s), value in actual[record].items():
        if i == j:
            key = bad, bad, r, s
            reset[key] = reset.get(key, QQ_I.zero)+value
    assert h.trace(reset) == h.trace(actual[record])
    actual[record] = reset
    assert sum(map(h.trace, actual.values())) == 1
    with pytest.raises(ArithmeticError, match='independent target'):
        receipt.audit_output(actual, exact_reference_target, cap=1, reference=True)


@pytest.mark.parametrize('stage', ('preparation', 'apparatus', 'readout'))
def test_receipt_operation_gate_rejects_a_bad_source_component(monkeypatch, stage):
    certify = source.certify

    def changed(kind, **kwargs):
        result = certify(kind, **kwargs)
        if kind == stage:
            result = {**result, 'operator_error_upper': 2*receipt.COMPONENT_GATE}
        return result

    monkeypatch.setattr(source, 'certify', changed)
    with pytest.raises(ArithmeticError, match='operation budget'):
        receipt.bound()
