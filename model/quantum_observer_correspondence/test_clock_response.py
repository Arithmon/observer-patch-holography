"""Independent endpoint checks of the complete quantum clock/response law."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

import pytest

from . import clock_response as c
from . import construction as producer


def marginal(joint, positions):
    values = defaultdict(F)
    for key, value in joint.items():
        values[tuple(key[i] for i in positions)] += value
    return dict(values)


@pytest.mark.parametrize('gamma,reference,effective', product((0, 1), (False, True), (False, True)))
def test_independent_endpoints_match_every_output_matrix(gamma, reference, effective):
    expected = c.execute(gamma, reference=reference)
    actual = producer.execute(gamma, reference=reference, effective=effective)
    assert len(expected) == 1442
    assert expected == actual
    assert sum((c.trace(state) for state in expected.values()), F(0)) == 1
    assert any(not state for state in expected.values())


def test_clock_marginal_fixed_while_coherent_response_changes():
    laws = [c.clock_joint(c.execute(gamma)) for gamma in (0, 1)]
    assert marginal(laws[0], (0, 1, 2)) == marginal(laws[1], (0, 1, 2))
    assert c.total_variation(*laws) == F(1848, 15625)
    # The very same completed clock event has a different detector response.
    assert laws[0]['complete', 1, 3, 1]/F(16, 625) == F(23, 50)
    assert laws[1]['complete', 1, 3, 1]/F(16, 625) == F(7, 10)


@pytest.mark.parametrize('gamma,distance,covariance', [
    (0, F(168097, 6000000), F(3269, 983040)),
    (1, F(231073, 6000000), F(7877, 983040)),
])
def test_all_separate_status_conditioned_marginals_match_but_joint_fails(gamma, distance, covariance):
    law = c.clock_joint(c.execute(gamma))
    control = c.matched_marginal_control(law)
    assert sum(control.values(), F(0)) == 1
    for positions in ((0,), (0, 1, 2), (0, 3)):
        assert marginal(law, positions) == marginal(control, positions)
    assert c.total_variation(law, control) == distance
    assert c.clock_detector_covariance(law) == covariance
    assert c.clock_detector_covariance(control) == 0


def test_bijective_clock_relabeling_preserves_obstruction_but_coarsening_can_remove_it():
    law = c.clock_joint(c.execute(1))
    control = c.matched_marginal_control(law)
    def relabel(values):
        return {(status, 3*a+2, 5*b-7, bit): value
                for (status, a, b, bit), value in values.items()}
    assert c.total_variation(relabel(law), relabel(control)) == c.total_variation(law, control)
    assert c.clock_detector_covariance(relabel(law), ticks=(5, 8)) == c.clock_detector_covariance(law)
    # Forgetting both clocks is deliberately many-to-one. No invariance is
    # claimed for this loss of the observation that discriminates the laws.
    assert c.total_variation(marginal(law, (0, 3)), marginal(control, (0, 3))) == 0


def test_timeout_mass_retained_in_joint_response_and_marginal_control():
    expected = {('timeout_first',): F(9, 25),
                ('timeout_second',): F(144, 625), ('complete',): F(256, 625)}
    for gamma in (0, 1):
        law = c.clock_joint(c.execute(gamma))
        assert marginal(law, (0,)) == expected
        assert marginal(c.matched_marginal_control(law), (0,)) == expected
        assert sum(value for key, value in law.items() if key[0] != 'complete') == F(369, 625)


def test_diagonal_preparation_removes_full_phase_response():
    before = c.execute(0, dephased=True)
    after = c.execute(1, dephased=True)
    assert before == after
    assert before == producer.execute(0, initial=producer.preparation(dephased=True))
    assert c.execute(0) != before


def test_untouched_reference_trace_equals_diagonal_control_without_erasing_full_coherences():
    outputs = [c.execute(gamma, reference=True) for gamma in (0, 1)]
    for gamma, history in enumerate(outputs):
        reduced = {key: c.reference_trace(state) for key, state in history.items()}
        assert reduced == c.execute(gamma, dephased=True)
    assert c.clock_joint(outputs[0]) == c.clock_joint(outputs[1])
    assert outputs[0] != outputs[1]
    assert any(r != s for history in outputs for state in history.values() for i, j, r, s in state)


def test_timeout_reference_entanglement_survives_independent_endpoint_replay():
    state, expectation = c.retained_entanglement()
    expected = producer.execute(initial=producer.entangled_preparation())[(), 'timeout_first', 0]
    assert state == expected
    assert c.trace(state) == F(9, 25)
    assert expectation == -F(17, 200)
    assert expectation/c.trace(state) == -F(17, 72)


def test_receipt_preserves_observation_and_physical_scope():
    result = c.audit()
    assert result['independent_native_endpoint_oracle']['complete_source_or_effective_comparisons'] == 8
    assert result['full_joint_response_tv'] == '1848/15625'
    assert result['reference_full_matrix_entry_l1_response_not_trace_norm'] == '176/625'
    assert result['invariance_under_clock_coarsening_claimed'] is False
    assert result['physical_clock_or_field_identification_established'] is False
