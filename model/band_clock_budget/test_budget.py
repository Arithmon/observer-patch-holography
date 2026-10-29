"""Independent schedule, resource and complete full-sector channel controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import budget
from ..band_clock import clock
from ..band_clock.experiment import retained_counts
from ..interacting_field import source


def test_smallest_split_count_under_the_previous_uniform_ceiling():
    target = F(clock.error_budget()['complete_average_halfdiamond_upper'])
    fixed = F(77, 2560)+F(clock.error_budget()['statistical_failure_upper'])
    coefficient = F(7, 4)**3*F(5, 4)**3/3
    assert budget.minimum_steps() == 21
    assert coefficient/21**2+fixed <= target < coefficient/20**2+fixed
    assert F(budget.error_budget()['complete_average_halfdiamond_upper']) == coefficient/21**2+fixed
    for parameter in (F(-1), F(-1, 2), F(0), F(1, 3), F(1)):
        n = budget.minimum_steps(parameter=parameter)
        assert F(budget.error_budget(parameter, n)['complete_average_halfdiamond_upper']) <= target
        if n > 1:
            assert F(budget.error_budget(parameter, n-1)['complete_average_halfdiamond_upper']) > target


def test_counted_source_savings_retain_both_pointer_services():
    space = clock.ClockField()
    duration = clock.rate_decision(retained_counts())['Hamiltonian_model_duration']
    old, new = (budget.resources(space, duration, n) for n in (32, 21))
    _, pointer = source.occupation()
    for key in source.COUNT_KEYS:
        fixed = 2*pointer['resources'][key]
        assert 32*(new[key]-fixed) == 21*(old[key]-fixed)
        assert type(new[key]) is int and 0 < new[key] < old[key]
    assert old['native_pulses'] == 26238982
    assert new['native_pulses'] == 17219334
    assert new['fresh_pointer_blanks'] == old['fresh_pointer_blanks'] == 2


@pytest.mark.parametrize('parameter', [F(-1), F(1, 2)])
def test_complete_field_histories_preserve_complex_reference(parameter):
    space = clock.ClockField()
    rng = np.random.default_rng(971)
    state = rng.normal(size=(space.dimension, 3))+1j*rng.normal(size=(space.dimension, 3))
    state /= np.linalg.norm(state)
    decision, outputs = budget.execute(retained_counts(), state, parameter, space)
    expected = clock.pointers(space, clock.exact_state(space, state, -float(parameter)*math.log(25/8)))
    assert decision['status'] == 'accept'
    assert set(outputs) == {('read', a, b) for a in (0, 1) for b in (0, 1)}
    assert clock.cq_distance(outputs.values(), expected) < float(F(
        budget.error_budget(parameter, budget.minimum_steps(parameter=parameter))[
            'complete_average_halfdiamond_upper']))
    np.testing.assert_allclose(sum(v.conj().T@v for v in outputs.values()),
                               state.conj().T@state, atol=2e-9)


def test_stops_keep_the_input_and_all_rate_failure_records():
    counts = {edge: (n, n//2) for edge, (n, _) in retained_counts().items()}
    space = clock.ClockField()
    state = np.zeros((space.dimension, 2), complex)
    state[0, 0], state[-1, 1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    decision, outputs = budget.execute(counts, state, space=space)
    assert decision['status'] == 'stop' and len(outputs) == 1
    assert np.array_equal(next(iter(outputs.values())), state)
    assert next(iter(outputs.values())) is not state


@pytest.mark.parametrize('target', [True, .04, F(0), F(1), F(1, 1000)])
def test_invalid_or_uncertifiable_targets_do_not_produce_a_schedule(target):
    with pytest.raises(ValueError):
        budget.minimum_steps(target)
