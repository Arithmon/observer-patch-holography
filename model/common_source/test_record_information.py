"""Coherence, reference, stopping and independent native-image controls."""
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from . import record_information as info, record_code as code
from ..observer import native
from ..native_modular_reference import local_preparation
from ..native_equilibration import history_entropy as population


def _pair_density(left, right, *, coherent=True, reference=False):
    states = native.states(binary=True)
    indices = [states.index(tuple(int(i in pair) for i in range(12)))
               for pair in (left, right)]
    r = 2 if reference else 1
    vector = np.zeros(66*r, complex)
    vector[indices[0]*r] = 1/math.sqrt(2)
    vector[indices[1]*r+int(reference)] = 1j/math.sqrt(2)
    rho = np.outer(vector, vector.conj())
    if not coherent:
        rho = np.diag(np.diag(rho))
    return rho


def test_counting_envelope_extends_to_full_quantum_outputs_with_same_native_reference():
    rho = np.diag(list(map(float, local_preparation.reference())))
    ledger = info.native_history_ledger(rho)
    assert abs(ledger['input_relative_entropy']) < 2e-12
    assert abs(ledger['output_relative_entropy']) < 2e-12
    assert abs(ledger['quantum_counted_functional']-ledger['input_entropy']) < 2e-12
    assert abs(ledger['population_counted_functional']-ledger['input_entropy']) < 2e-12
    assert ledger['charge_resolved']
    assert 3 < ledger['mean_native_attempts'] < 5
    complete = info.native_history_ledger(rho, law='C')
    assert not complete['charge_resolved']
    assert complete['quantum_counted_functional'] == pytest.approx(ledger['input_entropy'])
    assert complete['population_counted_functional'] > complete['quantum_counted_functional']
    assert abs(complete['quantum_envelope_residual']) < 2e-12


@pytest.mark.parametrize('law', ['C', 'E'])
def test_original_dyadic_complex_mixture_checks_entropy_components_independently(law):
    # Two orthonormal complex vectors with original exact dyadic entries.
    # The mixed input spectrum is exactly (3/4,1/4), independently of the
    # producer's spectral calculations and its relative-entropy subtraction.
    first, second = np.zeros(132, complex), np.zeros(132, complex)
    first[[0, 1, 8, 9]] = [.5, .5j, .5, -.5j]
    second[[4, 5, 12, 13]] = [.5, -.5j, .5, .5j]
    rho = .75*np.outer(first, first.conj())+.25*np.outer(second, second.conj())
    entropy = -.75*math.log(.75)-.25*math.log(.25)
    sigma = np.repeat(np.array(list(map(float, local_preparation.reference())))/2, 2)
    cross = -sum(float(rho[i, i].real)*math.log(sigma[i]) for i in range(132))
    row = info.native_history_ledger(rho, law=law, reference_dimension=2)
    assert abs(row['input_entropy']-entropy) < 2e-13
    assert abs(row['input_relative_entropy']-(cross-entropy)) < 2e-13
    assert abs(row['quantum_counted_functional']-(entropy-math.log(2))
               -row['information_loss']) < 2e-13


def test_same_tape_can_hide_either_surviving_or_irreversibly_lost_coherence():
    for left, right, expected_loss in (((9, 10), (10, 11), 0.),
                                      ((0, 2), (1, 2), math.log(2))):
        pure = info.native_history_ledger(_pair_density(left, right))
        mixed = info.native_history_ledger(_pair_density(left, right, coherent=False))
        np.testing.assert_allclose(pure['record_probabilities'], mixed['record_probabilities'], atol=0)
        assert abs(pure['information_loss']-expected_loss) < 2e-12
        assert abs(mixed['information_loss']) < 2e-12
        assert abs(pure['population_counted_functional']-mixed['population_counted_functional']) < 2e-12
        assert abs(pure['input_relative_entropy']-pure['record_relative_entropy']
                   -pure['conditional_quantum_relative_entropy']-pure['information_loss']) < 2e-12
        assert abs(pure['quantum_envelope_residual']) < 2e-12


@pytest.mark.parametrize('law', ('C', 'E'))
def test_native_reference_entanglement_and_complete_timeouts_are_retained(law):
    rho = _pair_density((0, 2), (1, 2), reference=True)
    ledger = info.native_history_ledger(rho, law=law, reference_dimension=2)
    assert ledger['success_probability_exact'] == '1'
    assert ledger['reference_marginal_error'] < 2e-12
    assert ledger['information_loss'] > 0
    if law == 'C':
        assert not ledger['charge_resolved']
        flagged = info.native_history_ledger(rho, law=law, reference_dimension=2,
                                             success_only=True)
        assert flagged['success_probability_exact'] == '9/16'
        assert flagged['charge_resolved']
        assert abs(flagged['quantum_envelope_residual']) < 2e-12


@pytest.mark.parametrize('law', ('C', 'E'))
def test_stopped_corrected_code_retains_full_bell_information_while_tape_is_blind(law):
    bell = np.eye(6, dtype=complex).ravel()/math.sqrt(6)
    pure = info.protected_code_ledger(np.outer(bell, bell.conj()), law=law,
                                      reference_dimension=6)
    reset = info.protected_code_ledger(np.eye(36)/36, law=law, reference_dimension=6)
    assert pure['record_entropy'] == reset['record_entropy']
    assert pure['record_relative_entropy'] == 0
    assert pure['conditional_quantum_relative_entropy'] == pytest.approx(2*math.log(6))
    assert reset['conditional_quantum_relative_entropy'] == pytest.approx(0, abs=2e-12)
    assert pure['input_conditional_entropy'] == pytest.approx(-math.log(6))
    assert pure['output_conditional_entropy']-pure['record_entropy'] == pytest.approx(-math.log(6))
    assert pure['mean_attempts_exact'] == '91/36'
    assert pure['coin_only_subtraction_excess'] > 1
    assert pure['record_entropy'] == pytest.approx(float(F(91, 36))*pure['one_attempt_record_entropy'])


def _independent_passive_grams(law, word):
    """Original floor/ceiling primitive, explicit complex code columns."""
    columns = {}
    for logical, matching in enumerate(code.MATCHINGS):
        for pair in matching:
            initial = tuple(int(i in pair) for i in range(12))
            for coins in product((0, 1), repeat=len(word)):
                state, tape = initial, ()
                for (a, b), coin in zip(word, coins, strict=True):
                    low, odd = divmod(state[a]+state[b], 2)
                    following = list(state)
                    following[a], following[b] = low+coin*odd, low+(1-coin)*odd
                    label = state[a]-following[a] if law == 'C' else state[a]-state[b]
                    tape = (*tape, (a, b, coin, label))
                    state = tuple(following)
                rows = columns.setdefault(tape, {})
                vector = rows.setdefault(state, np.zeros(6, complex))
                vector[logical] += 1/math.sqrt(6*2**len(word))
    return {tape: np.stack(list(rows.values())).conj().T@np.stack(list(rows.values()))
            for tape, rows in columns.items()}


@pytest.mark.parametrize('law', ('C', 'E'))
def test_passive_two_event_loss_is_checked_against_full_quantum_branch_grams(law):
    word = ((0, 1), (2, 3))
    ledger = info.passive_code_information(law, word)
    grams = _independent_passive_grams(law, word)
    entropy = 0.
    for row in ledger['branches']:
        gram = grams[row['record']]
        expected = np.diag(row['logical_support_counts'])/(6*4)
        np.testing.assert_allclose(gram, expected, atol=2e-16)
        weight = np.trace(gram).real/6
        assert weight == pytest.approx(float(F(row['probability_exact'])))
        spectrum = np.linalg.eigvalsh(gram/np.trace(gram).real)
        entropy -= weight*sum(x*math.log(x) for x in spectrum if x > 0)
    assert ledger['information_loss'] == pytest.approx(math.log(6)-entropy)
    assert ledger['information_loss'] > .1
    assert ledger['positive_information_loss_exact']
    assert ledger['non_scalar_branch_grams'] > 0
    assert ledger['record_relative_entropy'] == 0
    if law == 'E':
        assert ledger['information_loss'] == pytest.approx(math.log(F(3, 2)))
        assert ledger['information_loss_prime_logs_exact'] == {'2': '-1', '3': '1'}
    else:
        assert ledger['information_loss_prime_logs_exact'] == {'2': '19/9', '3': '1', '5': '-55/36'}
    # A one-event uncorrected branch is still an isometry on the code. The
    # second incompatible observation, not the scalar coin, loses coherence.
    single = info.passive_code_information(law, word[:1])
    assert single['non_scalar_branch_grams'] == 0
    assert abs(single['information_loss']) < 2e-12


def test_small_positive_populations_are_retained_and_bad_original_inputs_rejected():
    for bad in (np.eye(6), np.zeros((6, 6)), np.diag([-.1]+[.22]*5),
                np.ones((6, 6), bool), [[True]+[0.]*5]+[[0.]*6]*5,
                np.full((6, 6), np.nan), np.ma.array(np.eye(6)/6, mask=False)):
        with pytest.raises(ValueError):
            info.protected_code_ledger(bad)
    for tiny in (1e-10, 1e-14, 1e-20):
        rho = np.diag([1-tiny, tiny, 0, 0, 0, 0])
        actual = info.protected_code_ledger(rho)['input_entropy']
        expected = -(1-tiny)*math.log(1-tiny)-tiny*math.log(tiny)
        assert actual > 0
        assert actual == pytest.approx(expected, rel=1e-14, abs=0)


@pytest.mark.parametrize('tiny', [1e-13, 1e-20, np.nextafter(0., 1.)])
def test_negative_original_populations_are_not_eigensolver_roundoff(tiny):
    rho = np.diag([1+tiny, -tiny, 0, 0, 0, 0])
    with pytest.raises(ValueError, match='original diagonal'):
        info.protected_code_ledger(rho)
    valid = np.diag([1-tiny, tiny, 0, 0, 0, 0])
    entropy = info.protected_code_ledger(valid)['input_entropy']
    assert entropy > 0
    expected = -(1-tiny)*math.log(1-tiny)-tiny*math.log(tiny)
    assert entropy == pytest.approx(expected, rel=1e-14, abs=0)


@pytest.mark.parametrize('option,value', [
    ('prefix', 3.0), ('prefix', F(3)), ('prefix', True), ('prefix', np.int64(3)),
    ('deadline', True), ('deadline', 3.0), ('deadline', F(3)), ('deadline', None),
    ('stop_on_current', 1), ('stop_on_current', 'false'), ('stop_on_current', np.bool_(True)),
    ('law', np.str_('E')), ('success_only', 0), ('reference_dimension', True),
])
def test_public_history_options_refused_with_cold_and_warm_population_cache(option, value):
    rho = np.eye(66)/66
    kwargs = {option: value}
    try:
        population.histories.cache_clear()
        with pytest.raises(ValueError):
            info.native_history_ledger(rho, **kwargs)
        for prefix in (1, 3):
            population.histories('E', prefix, 3, True)
        with pytest.raises(ValueError):
            info.native_history_ledger(rho, **kwargs)
    finally:
        population.histories.cache_clear()


@pytest.mark.parametrize('law', [None, True, ['C'], np.str_('C')])
def test_information_laws_preserve_original_types(law):
    for call in (lambda: info.code_stopped_weights(law),
                 lambda: info.protected_code_ledger(np.eye(6)/6, law=law),
                 lambda: info.passive_code_information(law),
                 lambda: info.field_reference_information(law)):
        with pytest.raises(ValueError):
            call()


@pytest.mark.parametrize('word', [((False, True),), ((0., 1.),), ([0, 1],), ((0, F(1)),)])
def test_passive_information_word_requires_original_native_integer_ports(word):
    with pytest.raises(ValueError):
        info.passive_code_information(word=word)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_same_field_protocol_separates_native_randomness_and_destructive_read_loss(law):
    from . import record_field as field
    actual = info.field_reference_information(law)
    _, arrays = field.detector_experiment()
    conditional = 0.
    for row in actual['complete_quantum_terminal_blocks']:
        label = row['terminal']
        key = ('initial' if label == 'first_selector_timeout' else
               f"after_first_{label[-1]}" if label.startswith('second_selector') else label)
        state = arrays[key]
        norm = np.vdot(state, state).real
        if norm < 1e-24:
            assert row['probability'] < 1e-24
            continue
        # Evaluate the conditional relative entropy with explicit joint
        # matrices in a thin orthonormal work basis. This uses the actual
        # full source amplitudes, including every reference column.
        _, reduced = np.linalg.qr(state, mode='reduced')
        vector = reduced.ravel()/math.sqrt(norm)
        reference = np.kron(reduced@reduced.conj().T/norm, np.eye(2)/2)
        values, vectors = np.linalg.eigh(reference)
        populations = abs(vectors.conj().T@vector)**2
        assert populations[values < 1e-12].sum() < 1e-12
        divergence = -sum(float(p)*math.log(float(v)) for p, v in
                          zip(populations, values, strict=True) if v > 1e-12)
        assert divergence == pytest.approx(row['conditional_joint_relative_entropy'])
        conditional += row['probability']*divergence
    assert conditional == pytest.approx(actual['conditional_quantum_relative_entropy'])
    assert actual['information_loss'] == pytest.approx(F(15, 16)*math.log(2))
    assert actual['information_loss'] < math.log(2)  # Keeping the first timeout matters.
    assert actual['native_recovery_information_loss'] == 0
    assert actual['second_read_additional_information_loss'] == 0
    assert actual['record_relative_entropy'] == 0
    assert actual['complete_reference_marginal_error'] < 2e-12
    native_entropy = (math.log(12)-5/6*math.log(5) if law == 'C'
                      else math.log(3)+2/3*math.log(2))
    tape = F(31, 16)*math.log(32)+6*F(465, 256)*native_entropy+F(15, 16)*math.log(2)
    assert actual['complete_joint_output_entropy'] == pytest.approx(tape)
    assert actual['native_randomness_subtracted_quantum_entropy'] == pytest.approx(actual['information_loss'])
