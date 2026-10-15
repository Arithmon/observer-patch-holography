"""Acceptance, portable arithmetic and receipt-boundary regressions."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import exact_output, output_certificate as c, protocol, smooth_history


def test_source_bindings_include_the_executed_compilers_data_inputs():
    from ..observer.source_scalar import compiler
    assert {compiler.RESPONSE, compiler.CONTROL_RECEIPT} <= set(c.SOURCE_PATHS)


@pytest.mark.parametrize('value', [F(0), F(1, 3), F(-1, 3), F(1, 10**400), F(-1, 10**400)])
def test_report_rounding_is_directed_from_exact_values(value):
    low, high = F(c.outward(value, upper=False)), F(c.outward(value))
    assert low <= value <= high
    assert high-low <= F(1, c.GRID)
    if value:
        assert high != 0 if value > 0 else low != 0


@pytest.mark.parametrize('value', [True, 0., float('nan'), float('inf'), -F(1), F(2)])
def test_acceptance_cannot_pass_nonfinite_wrong_type_or_wrong_bound(value):
    with pytest.raises(ArithmeticError):
        c._gate(value, F(1), 'hostile value')


def test_exact_event_mass_preserves_tiny_branches():
    vectors = {(0, 0, 0, 1): np.array([1e-200+0j]), (0, 0, 0, 0): np.array([1.+0j])}
    assert c._event_probability(vectors) == F(1e-200)**2 > 0


def test_locked_feedback_control_changes_only_both_enable_flags():
    original = protocol.continuation_program(4, 512)
    changed = c.feedback_control(original)
    assert changed.native_current == original.native_current
    assert sum(op.time for op in changed.operations) == F(5, 2)
    assert [op.record for op in changed.operations] == [op.record for op in original.operations]
    assert sum(op.kind in protocol.READS for op in changed.operations) == 4
    for before, after in zip(original.operations, changed.operations, strict=True):
        if before.kind == 'feedback':
            assert before.enabled and not after.enabled
        else:
            assert before == after


@pytest.mark.parametrize('wrong', [
    {'all_pass': 1, 'count': 16, 'bound': '1/10', 'extra': None},
    {'all_pass': True, 'count': 16., 'bound': '1/10', 'extra': None},
    {'all_pass': True, 'count': 16, 'bound': '0', 'extra': None},
    {'all_pass': True, 'count': 16, 'bound': '1/10'},
    {'all_pass': True, 'count': 16, 'bound': '1/10', 'extra': None, 'invented': True},
])
def test_receipt_cannot_replace_certificates_or_scalar_types(wrong):
    actual = {'all_pass': True, 'count': 16, 'bound': '1/10', 'extra': None}
    c.compare_receipt(actual, actual.copy())
    with pytest.raises(ArithmeticError):
        c.compare_receipt(wrong, actual)


def test_receipt_list_shape_and_elements_are_checked():
    for wrong in ([], [True], [1, 2], (1,)):
        with pytest.raises(ArithmeticError):
            c.compare_receipt(wrong, [1])


def test_certificate_rejects_wrong_but_complete_source_feedback(monkeypatch):
    program = protocol.Program(native_current=1, operations=(
        protocol.Operation('gauge_read'), protocol.Operation('feedback', record=0)))
    target = smooth_history.enclose(program, dimension=2)
    def compare():
        output = protocol.execute(program)
        return exact_output.compare_vectors(program,
            {branch.record: branch.state[:, 0] for branch in output.branches}, target)
    correct = compare()
    assert correct['distance_upper'] < F(1, 10**10)
    original = protocol._operation_maps
    def wrong(p, operation, bit, backend):
        maps, cost = original(p, operation, bit, backend)
        return ((np.eye(p.work_dimension),), cost) if operation.kind == 'feedback' else (maps, cost)
    monkeypatch.setattr(protocol, '_operation_maps', wrong)
    result = compare()
    assert result['distance_lower'] > F(1, 2)
    assert abs(result['source_mass_real']-1) < F(1, 10**10)
    with pytest.raises(ArithmeticError):
        c._summary(result, F(1, 1000))
