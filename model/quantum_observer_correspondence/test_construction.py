from fractions import Fraction as F
from itertools import product

import pytest
import sympy as sp

from . import construction as q


@pytest.fixture(scope='module')
def receipt():
    return q.run()


def test_all_matrix_units_and_actual_local_service_commutators(receipt):
    certificate = receipt['all_input_structural_certificate']
    assert certificate['matrix_units_checked_against_analytic_native_reduction'] == 4356
    assert certificate['exact_remote_Kraus_local_service_commutators'] > 500
    assert certificate['occupation_phase_rectangle_ratio'] == '-1'


def test_resolvent_is_inverse_on_every_matrix_unit():
    blocks, _, _ = q.resolvent()
    assert sum(len(pairs) for pairs, _ in blocks) == 4356
    assert max(len(pairs) for pairs, _ in blocks) == 36
    for pairs, _ in blocks:
        for i, j in pairs:
            unit = {(i, j, 0, 1): F(1)}
            assert q.decode(q.lift(unit)) == unit


@pytest.mark.parametrize('reference,gamma', product((False, True), (0, 1)))
def test_full_outputs_preserve_coherence_reference_and_timeouts(reference, gamma):
    direct = q.execute(gamma, reference=reference)
    effective = q.execute(gamma, reference=reference, effective=True)
    checked = q.compare(direct, effective)
    assert checked['complete_output_words'] == 1442
    assert checked['timeout_probability'] == '369/625'


@pytest.mark.parametrize('cap', (1, 3))
def test_cap_changes_preserve_full_output_identity_and_exact_survival(cap):
    direct = q.execute(cap=cap)
    effective = q.execute(cap=cap, effective=True)
    checked = q.compare(direct, effective)
    assert F(checked['timeout_probability']) == 1-(1-F(3, 5)**cap)**2


def test_reference_entanglement_survives_in_unresolved_output(receipt):
    certificate = receipt['surviving_reference_entanglement']
    assert certificate['unnormalized_partial_transpose_expectation'] == '-17/200'
    assert certificate['conditional_partial_transpose_expectation'] == '-17/72'
    assert certificate['untouched_reference_marginal'] == 'I2/2'


def test_imaginary_bell_coherence_preserves_every_complex_reference_output():
    a, b = q.address((0, 1)), q.address((0, 5))
    initial = {(a, a, 0, 0): F(1, 2), (b, b, 1, 1): F(1, 2),
               (a, b, 0, 1): -sp.I/2, (b, a, 1, 0): sp.I/2}
    direct = q.execute(reference=True, initial=initial)
    effective = q.execute(reference=True, initial=initial, effective=True)
    checked = q.compare(direct, effective)
    assert checked['complete_output_words'] == 1442
    assert checked['timeout_probability'] == '369/625'
    assert any(sp.im(value) != 0 for state in direct.values() for value in state.values())


def test_coherent_interaction_changes_readout_but_not_dephased_control():
    histories = [q.execute(gamma) for gamma in (0, 1)]
    assert histories[0] != histories[1]
    detector = [sum((q.trace(state) for (_, _, bit), state in h.items() if bit), F(0))
                for h in histories]
    assert detector[0] != detector[1]
    assert q.execute(0, initial=q.preparation(dephased=True)) == q.execute(
        1, initial=q.preparation(dephased=True))


def test_quantum_encoder_required_for_nonstationary_preparation():
    i = q.address((2, 3))
    initial = {(i, i, 0, 0): F(1)}
    assert q.lift(initial) != initial
    assert q.compare(q.execute(initial=initial), q.execute(initial=initial, effective=True))[
        'total_probability'] == '1'
    _, native_timeout = q.interval(initial)
    _, wrong_effective_timeout = q.interval(initial, effective=True)
    assert q.decode(wrong_effective_timeout) != native_timeout


def test_cone_and_cross_boundary_controls(receipt):
    cone = receipt['transported_cone_control']
    assert F(cone['decoder_negative_population']) == -F(1, 24)
    assert F(cone['cap_one_negative_population']) < 0
    assert F(receipt['boundary_intervention_control']['exact_AU_minus_UA_entry']) == F(4, 375)


def test_complete_output_checker_rejects_missing_zero_and_reference_corruption():
    history = q.execute(reference=True)
    zero = next(key for key, state in history.items() if not state)
    missing = dict(history)
    del missing[zero]
    with pytest.raises(ArithmeticError, match='menu'):
        q.compare(history, missing)
    altered = {key: dict(state) for key, state in history.items()}
    key = next(key for key, state in altered.items() if state)
    # Preserve traces while corrupting a reference coherence of the output.
    altered[key][0, 0, 0, 1] = altered[key].get((0, 0, 0, 1), F(0))+F(1, 1000)
    with pytest.raises(ArithmeticError, match='source/reference'):
        q.compare(history, altered)


@pytest.mark.parametrize('gamma,cap', [(True, 2), (2, 2), (1, True), (1, 0), (1, 4)])
def test_protocol_controls_reject_invalid_parameters(gamma, cap):
    with pytest.raises(ValueError):
        q.execute(gamma, cap)
