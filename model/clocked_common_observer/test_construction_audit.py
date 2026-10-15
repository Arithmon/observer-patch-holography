"""Independent record, continuation, reference and hostile-envelope controls."""
from dataclasses import replace

import numpy as np
import pytest

from ..common_source import protocol
from ..primitive_observer import common, runtime
from . import clock, construction as observer


@pytest.fixture(scope='module')
def work():
    return observer.execute_work(2, steps=2)


@pytest.fixture(scope='module')
def native():
    return clock.execute_native('C', 2)


def test_waiting_state_is_consumed_and_preparation_is_charged_once(work):
    checked = observer.validate_work(work)
    for n, current, history in checked.successes:
        assert all(dict(branch.resources)['work_preparations'] == 1 for branch in history.branches)
        assert dict(history.expanded_tree_resources)['work_preparations'] == 1
        reset = common.execute(protocol.continuation_program(2, 2, native_current=current))
        assert runtime.complete_distance(history, reset) > 1e-4
        assert len(history.program.events) == n+len(reset.program.events)
        np.testing.assert_array_equal(history.initial_density, checked.prefixes[0].initial_density)


def test_arbitrary_work_reference_survives_the_complete_continuation():
    # Nonmaximal complex work/reference entanglement, followed by a mixture.
    vector = np.zeros(16*3, complex)
    vector[[0, 4, 23, 47]] = [1, 2j, 1-1j, 3]
    vector /= np.linalg.norm(vector)
    initial = .75*np.outer(vector, vector.conj())
    initial[6, 6] += .25
    before = np.trace(initial.reshape(16, 3, 16, 3), axis1=0, axis2=2)
    supplied = observer.validate_work(observer.execute_work(1, steps=2, initial=initial))
    for _, _, history in supplied.successes:
        assert history.reference_dimension == 3
        assert history.preparation_mode == 'supplied'
        assert all(dict(branch.resources)['work_preparations'] == 0 for branch in history.branches)
        after = sum((np.trace(branch.density.reshape(16, 3, 16, 3), axis1=0, axis2=2)
                     for branch in history.branches), np.zeros((3, 3), complex))
        np.testing.assert_allclose(after, before, atol=2e-11, rtol=0)
        # The public replay boundary must reconstruct this complete prefix.
        restored = runtime.checkpoint_from_json(history.program, runtime.checkpoint_to_json(history))
        assert runtime.complete_distance(restored, history) == 0


def test_cap_extension_preserves_earlier_outputs_but_not_a_coarsened_timeout_state(work, native):
    small = tuple(observer.output_rows(work, clock.execute_native('C', 1), cap=1))
    large = tuple(observer.output_rows(work, native))
    assert len(small) == 33 and len(large) == 65
    for earlier, later in zip(small[:-1], large[:32], strict=True):
        assert earlier[0] == later[0]
        np.testing.assert_array_equal(earlier[1], later[1])
        np.testing.assert_array_equal(earlier[2], later[2])
    _, native_tail, work_tail = small[-1]
    old_timeout = np.trace(native_tail)*work_tail
    late = sum((np.trace(nr)*wr for (n, _, _), nr, wr in large if n == 2),
               np.zeros((16, 16), complex))
    # Classical survival coarsening holds, but it cannot undo later quantum
    # evolution, work measurement, or the newly triggered current kick.
    assert abs(np.trace(old_timeout)-np.trace(late)) < 2e-11
    assert np.linalg.norm(old_timeout-late, 'fro') > 1e-4


@pytest.mark.parametrize('change', ('missing', 'duplicate', 'boolean', 'prefix', 'history'))
def test_incomplete_or_relabelled_work_families_fail_closed(work, native, change):
    if change == 'missing':
        altered = replace(work, successes=work.successes[:-1])
    elif change == 'duplicate':
        altered = replace(work, successes=(work.successes[0], *work.successes[:-1]))
    elif change == 'boolean':
        _, j, history = work.successes[0]
        altered = replace(work, successes=((True, j, history), *work.successes[1:]))
    elif change == 'prefix':
        altered = replace(work, prefixes=work.prefixes[:-1])
    else:
        n, j, history = work.successes[0]
        history = replace(history, branches=history.branches[:-1])
        altered = replace(work, successes=((n, j, history), *work.successes[1:]))
    with pytest.raises(ValueError):
        observer.output_rows(altered, native)


def test_resealed_resource_claim_cannot_charge_a_second_preparation(work):
    n, j, history = work.successes[0]
    branch = history.branches[0]
    costs = dict(branch.resources)
    costs['work_preparations'] += 1
    branch = replace(branch, resources=tuple(sorted(costs.items())))
    changed = replace(history, branches=(branch, *history.branches[1:]), seal='')
    changed = replace(changed, seal=runtime._seal(changed))
    altered = replace(work, successes=((n, j, changed), *work.successes[1:]))
    with pytest.raises(ValueError, match='resources'):
        observer.validate_work(altered)


@pytest.mark.parametrize('change', ('missing', 'boolean', 'shape', 'nan', 'masked', 'mass'))
def test_native_output_menu_and_original_matrix_domain_are_checked(work, native, change):
    changed = dict(native)
    first = next(iter(changed))
    if change == 'missing':
        changed.pop(first)
    elif change == 'boolean':
        changed[True, first[1]] = changed.pop(first)
    elif change == 'shape':
        changed[first] = np.zeros((65, 65))
    elif change == 'nan':
        changed[first] = np.full((66, 66), np.nan)
    elif change == 'masked':
        changed[first] = np.ma.array(changed[first], mask=False)
    else:
        changed[first] = np.zeros((66, 66))
    with pytest.raises(ValueError):
        observer.output_rows(work, changed)


def test_lazy_output_detaches_mutable_native_arrays_and_work_envelopes(work, native):
    arrays = {record: value.copy() for record, value in native.items()}
    mutable = replace(work, successes=tuple((n, j, replace(history,
        branches=tuple(replace(branch) for branch in history.branches)))
        for n, j, history in work.successes))
    expected = next(observer.output_rows(work, native))
    stream = observer.output_rows(mutable, arrays)
    arrays[1, -1][:] = 0
    mutable.successes[0][2].branches[0].__dict__['density'] = np.zeros((16, 16))
    actual = next(stream)
    assert actual[0] == expected[0]
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(actual[2], expected[2])
    assert not actual[1].flags.writeable and not actual[2].flags.writeable
