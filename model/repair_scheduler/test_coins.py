"""Independent source-map and adversarial controls for native randomness."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from model.observer import native, quantum_mixing
from model.predictions import observer_instrument as existing
from . import coins


SEAMS = quantum_mixing.seams()


def _oracle_maps(law, seam, binary):
    """Relabel the pre-existing complete instrument, not the new constructor."""
    ports = (*seam, *(p for p in range(12) if p not in seam))
    states = native.states(binary=binary)
    index = {state: i for i, state in enumerate(states)}

    def relabel(state):
        output = [0]*12
        for old, new in enumerate(ports):
            output[new] = state[old]
        return tuple(output)

    result = {}
    for record, support in existing.kraus_support(law).items():
        mapped = []
        for source, target in support:
            left, right = relabel(existing.STATES[source]), relabel(existing.STATES[target])
            if left in index:
                mapped.append((index[left], index[right]))
        if mapped:
            result[record] = tuple(sorted(mapped))
    return result


def _matrices(law, seam, binary):
    d = 66 if binary else 78
    result = {}
    for record, support in _oracle_maps(law, seam, binary).items():
        matrix = np.zeros((d, d), dtype=np.int64)
        for source, target in support:
            matrix[target, source] = 1
        result[record] = matrix
    return result


def _marginal(operator, dimension):
    ref = len(operator)//dimension
    return np.trace(operator.reshape(dimension, ref, dimension, ref), axis1=0, axis2=2)


def _opposite_pair(binary, seam=SEAMS[0]):
    states = native.states(binary=binary)
    a, b = seam
    spectator = next(p for p in range(12) if p not in seam)
    x, y = [0]*12, [0]*12
    x[a] = x[spectator] = y[b] = y[spectator] = 1
    return states.index(tuple(x)), states.index(tuple(y))


def _entangled(binary, reference=2):
    d = 66 if binary else 78
    x, y = _opposite_pair(binary)
    vector = np.zeros(d*reference, complex)
    vector[x*reference], vector[y*reference+1] = .6, .8j
    return np.outer(vector, vector.conj())


@pytest.mark.parametrize('binary', [False, True])
@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('seam', SEAMS)
def test_all_source_seams_match_existing_instrument_and_exact_coin_effect(law, seam, binary):
    actual = {(coin, label): tuple(sorted(mapping))
              for coin, label, mapping in coins.coin_maps(law, seam, binary)}
    assert actual == _oracle_maps(law, seam, binary)
    d = 66 if binary else 78
    for coin in (0, 1):
        gram = np.zeros((d, d), dtype=np.int64)
        for (c, _), matrix in _matrices(law, seam, binary).items():
            if c == coin:
                gram += matrix.T@matrix
        # Actual Kraus matrices are M/sqrt(2); this exact integer equality
        # checks every diagonal and off-diagonal effect, including collisions.
        np.testing.assert_array_equal(gram, np.eye(d, dtype=np.int64))


def test_effect_receipt_counts_actual_complete_source_inventory():
    report = coins.effect_certificate()
    assert report['source_dimensions'] == [78, 66]
    assert report['partial_injection_branches_checked'] == 780
    assert report['input_coin_addresses_checked'] == 17280
    assert report['coarse_coin_effect'] == 'I/2'


@pytest.mark.parametrize('binary', [False, True])
@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('coin', [0, 1])
@pytest.mark.parametrize('reference', [1, 2, 3])
def test_array_action_matches_independent_kraus_arrays_on_arbitrary_operators(
        law, coin, binary, reference):
    d = 66 if binary else 78
    rng = np.random.default_rng(3301)
    operator = (rng.normal(size=(d*reference, d*reference))
                + 1j*rng.normal(size=(d*reference, d*reference)))/32
    original = operator.copy()
    expected = np.zeros_like(operator)
    seam = SEAMS[7]
    for (c, _), matrix in _matrices(law, seam, binary).items():
        if c == coin:
            full = np.kron(matrix, np.eye(reference))
            expected += full@operator@full.T/2
    actual = coins.apply_coin(operator, law, coin, seam=seam, binary=binary)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-15)
    np.testing.assert_array_equal(operator, original)
    np.testing.assert_allclose(_marginal(actual, d), _marginal(original, d)/2,
                               rtol=0, atol=2e-15)


@pytest.mark.parametrize('binary', [False, True])
def test_reused_correlated_auxiliary_keeps_every_coin_word_fair_with_active_reference(binary):
    d = 66 if binary else 78
    initial = _entangled(binary)
    marginal = _marginal(initial, d)
    x = np.array([[0, 1], [1, 0]], complex)
    rotations = (3*np.eye(2)/5+4j*x/5, np.diag([1, (3+4j)/5]))
    assert np.linalg.norm(rotations[0]@rotations[1]-rotations[1]@rotations[0]) > .5
    for word in product((0, 1), repeat=3):
        state, expected = initial.copy(), marginal.copy()
        prefix = 0
        for step, coin in enumerate(word):
            law = 'C' if (step+prefix) % 2 == 0 else 'E'
            state = coins.apply_coin(state, law, coin, seam=SEAMS[prefix % 30], binary=binary)
            unitary = rotations[(step+coin) % 2]
            full = np.kron(np.eye(d), unitary)
            state = full@state@full.conj().T
            expected = unitary@expected@unitary.conj().T/2
            prefix = 2*prefix+coin+1
            np.testing.assert_allclose(_marginal(state, d), expected, rtol=0, atol=2e-16)
        assert np.trace(state).real == pytest.approx(1/8, abs=2e-16)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_actual_selector_reuses_auxiliary_and_retains_rejection_timeout(law):
    report = coins.auxiliary_execution(law, rounds=3)
    assert report['auxiliary_replaced_between_attempts'] is False
    assert report['initial_auxiliary_reference_entangled'] is True
    records = report['all_accepted_round_address_records']
    assert {(row['round'], row['address']) for row in records} == set(product(range(1, 4), range(30)))
    for row in records:
        assert row['probability'] == pytest.approx(float(F(1, 32)*F(1, 16)**(row['round']-1)), abs=1e-16)
    assert report['timeout_probability_diagnostic'] == pytest.approx(1/4096, abs=1e-16)
    assert report['maximum_subnormalized_reference_distance_diagnostic'] < 2e-15


@pytest.mark.parametrize('size, rounds', [(2, 1), (3, 4), (30, 1), (30, 7), (255, 32), (256, 3)])
def test_selector_exact_weights_agree_with_independent_word_enumeration(size, rounds):
    bits = (size-1).bit_length()
    counts = Counter(coins.decode_coin_word(word, size) for word in product((0, 1), repeat=bits))
    assert all(counts[address] == 1 for address in range(size))
    reject = F(counts[None], 2**bits)
    report = coins.selector_certificate(size, rounds)
    total = F(report['timeout_probability'])
    assert total == reject**rounds
    for r, row in enumerate(report['round_records']):
        assert row['round'] == r+1
        assert F(row['probability_per_address']) == reject**r/F(2**bits)
        total += size*F(row['probability_per_address'])
    assert total == 1
    assert F(report['conditional_address_probability']) == F(1, size)
    assert report['maximum_auxiliary_native_attempts'] == bits*rounds
    assert F(report['expected_auxiliary_native_attempts']) == bits*sum(
        (reject**r for r in range(rounds)), F())


def test_modulo_and_shared_hidden_bit_are_not_valid_uniform_seam_sources():
    modulo = Counter(value % 30 for value in range(32))
    distance = sum((abs(F(modulo[value], 32)-F(1, 30)) for value in range(30)), F())/2
    assert distance == F(7, 120)
    # Every position of this shared-bit process is marginally fair, but the
    # all-repeated five-bit word is not uniform and rejects with mass 1/2.
    words = {(0,)*5: F(1, 2), (1,)*5: F(1, 2)}
    assert all(sum((weight*word[k] for word, weight in words.items()), F()) == F(1, 2)
               for k in range(5))
    timeout = sum((weight for word, weight in words.items()
                   if coins.decode_coin_word(word) is None), F())
    assert timeout == F(1, 2) != F(coins.selector_certificate(rounds=1)['timeout_probability'])


@pytest.mark.parametrize('law', ['C', 'E'])
def test_fine_auxiliary_record_steers_reference_and_same_carrier_is_disturbed(law):
    d = 66
    x, y = _opposite_pair(True)
    correlated = np.zeros((2*d, 2*d), complex)
    correlated[2*x, 2*x] = correlated[2*y+1, 2*y+1] = .5
    outputs = []
    for (coin, _), matrix in _matrices(law, SEAMS[0], True).items():
        if coin == 0:
            full = np.kron(matrix, np.eye(2))
            outputs.append(_marginal(full@correlated@full.T/2, d))
    np.testing.assert_array_equal(sum(outputs), np.eye(2)/4)
    assert any(np.array_equal(value, np.diag([.25, 0])) for value in outputs)
    # Reading only the coin is fair, but its action on that same carrier is
    # not an identity channel: the selected occupied endpoint moves.
    carrier = np.zeros((d, d), complex)
    carrier[x, x] = 1
    conditional = 2*coins.apply_coin(carrier, law, 0)
    assert np.trace(conditional).real == 1
    assert np.linalg.norm(conditional-carrier, 'nuc')/2 == 1


@pytest.mark.parametrize('bad', [True, False, 1., np.int64(30), F(30), None, '30', 0, 1, 257])
def test_selector_rejects_nonoriginal_or_out_of_range_sizes(bad):
    with pytest.raises(ValueError):
        coins.selector_certificate(bad, 2)


@pytest.mark.parametrize('bad', [True, False, 1., np.int64(2), F(2), None, '2', 0, -1, 33])
def test_selector_rejects_nonoriginal_or_out_of_range_rounds(bad):
    with pytest.raises(ValueError):
        coins.selector_certificate(30, bad)


@pytest.mark.parametrize('bad', [[0]*5, (0,)*4, (0,)*6, (0, 0, 0, 0, True),
                                  (0, 0, 0, 0, 1.), (0, 0, 0, 0, np.int64(1)),
                                  (0, 0, 0, 0, 2), None, '00000'])
def test_decoder_rejects_malformed_words(bad):
    with pytest.raises(ValueError):
        coins.decode_coin_word(bad)


@pytest.mark.parametrize('bad', [1, 0, np.bool_(True), None, 'binary'])
def test_coin_map_domain_checks_precede_cache_hits(bad):
    coins.coin_maps('C', SEAMS[0], True)
    coins.coin_maps('C', SEAMS[0], False)
    with pytest.raises(ValueError):
        coins.coin_maps('C', SEAMS[0], bad)


@pytest.mark.parametrize('seam', [list(SEAMS[0]), (), (0,), (0, 0), (0, 12),
                                   tuple(map(float, SEAMS[0])), tuple(map(F, SEAMS[0])),
                                   tuple(map(np.int64, SEAMS[0])), (False, True), None])
def test_coin_map_rejects_nonoriginal_source_seams_after_valid_cache(seam):
    coins.coin_maps('C', SEAMS[0], True)
    with pytest.raises(ValueError):
        coins.coin_maps('C', seam)


@pytest.mark.parametrize('law', [None, True, 1, 'C0', 'c', [], np.str_('C')])
def test_coin_map_rejects_invalid_laws_even_when_equal_to_cached_text(law):
    coins.coin_maps('C', SEAMS[0], True)
    with pytest.raises(ValueError):
        coins.coin_maps(law, SEAMS[0])


@pytest.mark.parametrize('coin', [True, False, 0., np.int64(1), F(1), -1, 2, None])
def test_array_action_requires_original_coin(coin):
    with pytest.raises(ValueError):
        coins.apply_coin(np.eye(66), 'C', coin)


@pytest.mark.parametrize('density', [
    [], np.zeros(66), np.zeros((0, 0)), np.zeros((65, 65)), np.zeros((66, 67)),
    np.zeros((264, 264)), np.eye(66, dtype=int), np.eye(66, dtype=bool),
    np.eye(66, dtype=object), np.full((66, 66), np.nan), np.full((66, 66), np.inf),
    np.ma.array(np.eye(66), mask=False),
])
def test_array_action_rejects_incomplete_nonfinite_or_ambiguous_operators(density):
    with pytest.raises(ValueError):
        coins.apply_coin(density, 'C', 0)


def test_array_action_rejects_nonfinite_computed_output_without_warning_dependency():
    density = np.full((78, 78), np.finfo(float).max)
    with np.errstate(over='ignore', invalid='ignore'):
        with pytest.raises(ArithmeticError, match='nonfinite'):
            coins.apply_coin(density, 'C', 0, binary=False)


def test_changed_source_coin_record_is_not_certified(monkeypatch):
    source = native.native

    def wrong_record(state, seam, coin):
        output, (_, current) = source(state, seam, coin)
        return output, (1-coin, current)

    monkeypatch.setattr(native, 'native', wrong_record)
    # Bypass only the memoization so the changed primitive is actually read.
    with pytest.raises(ArithmeticError, match='coin record'):
        coins._coin_maps.__wrapped__('C', SEAMS[0], False)
