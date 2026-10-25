"""The field protocol and endpoint observer continue on the same native carriers.

This checks the native interface. The field's effective cylinder dictionary
does not identify its hidden terminal state with the effective hidden state.
"""
from fractions import Fraction
import random

import pytest

from model.native_clock import passive_age
from model.observer_time_fields import preparation


@pytest.mark.parametrize('law', ('C', 'E'))
def test_field_feedback_then_current_only_age(law):
    generator = random.Random(112)
    work, auxiliary = (2,)+(0,)*11, (0, 2)+(0,)*10
    field = preparation.execute_history(work, auxiliary,
        iter(lambda: generator.getrandbits(1), None), law=law, feedback=True)
    assert field.flag == 'complete'
    assert field.terminal_work[1] == field.word[-1]
    assert field.read_clock == (1, 2+field.word[0], 3+field.word[0], 4+field.word[0])
    # No reconstruction from a field marginal and no reset between services.
    age = passive_age.execute_trials(1, trials=8, law=law, seed=119,
        initial_state=field.terminal_work, auxiliary_state=field.terminal_auxiliary)
    first = age['prefix'][0]
    start = passive_age.decode_current_tape(first['start_currents'], full=True)
    assert start.success
    assert start.initial == field.terminal_work
    assert age['initial_auxiliary_state'] == field.terminal_auxiliary
    counts = dict(age['records_and_counts'])
    inferred = passive_age.infer_counts(counts)
    assert inferred['all_trials_including_failures'] == 8
    # Eight trials do not establish the precise age; no silent precision upgrade.
    assert inferred['status'] == 'unresolved'
    assert not inferred['admission']['all_window_ages_uniquely_decodable_on_coverage']
    assert inferred['observed_mean'] == sum(
        (passive_age.score(record)*count for record, count in counts.items()), Fraction())/8
    used = field.native_attempts+age['main_native_attempts']+age['auxiliary_native_attempts']
    bound = (preparation.resources(feedback=True)['all_native_attempts_max']
             +passive_age.budget(8, 4)['total_native_attempts_upper'])
    assert 0 < used <= bound
