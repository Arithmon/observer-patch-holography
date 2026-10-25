"""Native passive-age identity, complete controllers and blind confidence."""
from copy import deepcopy
from fractions import Fraction as F
import math
import random

import numpy as np
import pytest

from . import passive_age as p


def test_source_eigenfield_uses_every_native_column_and_operator_unit():
    proof = p.exact_identity()
    assert p.pairs() == ((0, 3), (1, 2), (4, 7), (5, 6), (8, 11), (9, 10))
    assert len(proof['rows']) == 12
    assert all(row['diagonal_columns'] == 66 and row['complete_operator_matrix_units'] == 4356
               for row in proof['rows'])
    assert proof['untouched_reference_extension']
    assert proof['native_auxiliary_coin_effect']['coarse_coin_effect'] == 'I/2'
    assert proof['native_auxiliary_coin_effect']['independent_auxiliary_initial_state_required'] is False


def test_exact_identity_rejects_a_changed_actual_native_output(monkeypatch):
    original = p.coins.coin_maps
    def altered(law, seam, binary=True):
        branches = list(original(law, seam, binary))
        coin, label, mapping = branches[0]
        before, after = mapping[0]
        branches[0] = coin, label, ((before, (after+1) % 66),)+mapping[1:]
        return tuple(branches)
    monkeypatch.setattr(p.coins, 'coin_maps', altered)
    with pytest.raises(ArithmeticError):
        p.exact_identity()


@pytest.mark.parametrize('law', ['C', 'E'])
def test_coarse_coin_preserves_an_arbitrary_entangled_reference(law):
    rng = np.random.default_rng(17)
    v = rng.normal(size=132)+1j*rng.normal(size=132)
    v /= np.linalg.norm(v)
    rho = np.outer(v, v.conj())
    marginal = np.einsum('aras->rs', rho.reshape(66, 2, 66, 2))
    for coin in (0, 1):
        out = p.coins.apply_coin(rho, law, coin)
        np.testing.assert_allclose(np.einsum('aras->rs', out.reshape(66, 2, 66, 2)), marginal/2, atol=1e-14)


def test_current_only_decoder_reconstructs_discovery_from_full_collision_input():
    rng = random.Random(71)
    execution = p.discovery.execute((2,)+(0,)*11, p._coin_stream(rng), full=True, law='E')
    decoded = p.decode_current_tape(execution.outcome.currents, full=True)
    assert decoded == execution.outcome
    assert decoded.success and decoded.initial == (2,)+(0,)*11
    assert decoded.final == execution.output
    assert p.discovery.word_support(execution, full=True) == {decoded.initial: decoded.final}
    pair, initial_y = p.choose_pair(decoded)
    assert initial_y in (1, 2)
    assert sum(decoded.final[i] for i in p.pairs()[pair]) == initial_y
    with pytest.raises(TypeError):
        p.decode_current_tape(execution.outcome.currents, age=2)
    with pytest.raises(ValueError):
        p.decode_current_tape(execution.outcome.currents[:-1], full=True)
    with pytest.raises(ValueError):
        p.decode_current_tape(execution.outcome.currents+execution.outcome.currents[:1], full=True)


class Bits:
    def __init__(self, values):
        self.values = iter(values)
    def getrandbits(self, count):
        assert count == 1
        return next(self.values)


def test_scheduler_executes_native_coins_and_retains_rejection_and_timeout(monkeypatch):
    import hashlib
    actual = p.discovery.native.native
    calls = []
    def counted(state, seam, coin):
        calls.append((state, seam, coin))
        return actual(state, seam, coin)
    monkeypatch.setattr(p.discovery.native, 'native', counted)
    initial = p.discovery.native.states(binary=True)[0]
    address, final, attempts, rejected = p._selector(initial, Bits([1]*5+[0]*5), 'C', 2, hashlib.sha256())
    assert address == 0 and attempts == 10 and rejected == 1 and len(calls) == 10
    state = initial
    for _, seam, coin in calls:
        state, _ = actual(state, seam, coin)
    assert final == state
    address, _, attempts, rejected = p._selector(final, Bits([1]*40), 'E', 8, hashlib.sha256())
    assert address is None and attempts == 40 and rejected == 8


@pytest.mark.parametrize('law', ['C', 'E'])
def test_actual_small_history_has_no_reset_and_replays_public_current_prefix(law):
    initial = p.discovery.native.states(binary=True)[23]
    auxiliary = (0,)*4+(2,)+(0,)*7
    sample = p.execute_trials(3, trials=12, law=law, seed=733,
                              initial_state=initial, auxiliary_state=auxiliary)
    assert sample['initial_state'] == initial
    assert sample['initial_auxiliary_state'] == auxiliary
    assert sample['same_main_and_auxiliary_carriers_reused']
    assert sample['microscopic_native_calls_executed']
    assert sum(value for _, value in sample['records_and_counts']) == 12
    assert sum(count for _, _, count in sample['terminal_state_counts']) == 12
    assert sample['main_native_attempts'] > sample['passive_main_attempts'] == 36
    assert sample['auxiliary_native_attempts'] >= 5*36
    for entry in sample['prefix']:
        start = p.decode_current_tape(entry['start_currents'], full=entry['start_full'])
        end = p.decode_current_tape(entry['end_currents'])
        assert p.endpoint_record(start, end) == entry['record']
        assert end.final == entry['terminal_state']
    assert sample == p.execute_trials(3, trials=12, law=law, seed=733,
                                       initial_state=initial, auxiliary_state=auxiliary)


def test_complete_finite_discovery_failures_remain_in_actual_history():
    sample = p.execute_trials(2, trials=32, seed_scans=1, retries=1, seed=119)
    counts = dict(sample['records_and_counts'])
    assert set(counts) == set(p.alphabet())
    assert sum(n for record, n in counts.items() if record[0] != 'read') > 0
    result = p.infer_counts(counts, seed_scans=1, retries=1)
    assert result['all_trials_including_failures'] == 32
    assert result['status'] == 'unresolved'
    assert not result['successful_records_conditioned_on']


def _model_counts(age):
    # Algebraic input fixture only, not an execution receipt.
    theta = F(9, 10)**age
    upper_probability = (2*theta+1)/3
    denominator = upper_probability.denominator
    n = denominator*math.ceil(p.TRIALS/denominator)
    counts = dict.fromkeys(p.alphabet(), 0)
    counts['read', 0, 1, 1] = int(n*upper_probability)
    counts['read', 0, 1, 0] = n-counts['read', 0, 1, 1]
    return counts


@pytest.mark.parametrize('age', range(5))
def test_blind_inverse_identifies_each_admitted_age_from_complete_counts(age):
    result = p.infer_counts(_model_counts(age))
    assert result['candidate_attempt_counts'] == (age,)
    assert result['status'] == 'identified'
    assert result['observed_mean'] == F(9, 10)**age
    assert not result['age_or_hidden_configuration_input']
    assert result['admission']['coverage_failure_upper'] == F(1, 4096)


def test_inverse_has_no_age_hidden_state_or_target_weight_access(monkeypatch):
    counts = _model_counts(3)
    def forbidden(*args, **kwargs):
        raise AssertionError('inverse accessed an execution target')
    monkeypatch.setattr(p, 'execute_trials', forbidden)
    monkeypatch.setattr(p.discovery.native, 'native', forbidden)
    assert p.infer_counts(counts)['candidate_attempt_counts'] == (3,)
    with pytest.raises(TypeError):
        p.infer_counts(counts, age=3)


@pytest.mark.parametrize('mutation', ['missing', 'extra', 'bool', 'float', 'negative', 'key_bool', 'zero'])
def test_blind_inverse_rejects_incomplete_or_coerced_counts(mutation):
    counts = _model_counts(2)
    if mutation == 'missing': counts.pop(('start_failure', 'seed_timeout'))
    elif mutation == 'extra': counts['unregistered',] = 0
    elif mutation == 'bool': counts['start_failure', 'seed_timeout'] = False
    elif mutation == 'float': counts['start_failure', 'seed_timeout'] = 0.0
    elif mutation == 'negative': counts['start_failure', 'seed_timeout'] = -1
    elif mutation == 'zero': counts = dict.fromkeys(counts, 0)
    else:
        value = counts.pop(('read', 0, 1, 1))
        counts['read', False, 1, 1] = value
    with pytest.raises(ValueError):
        p.infer_counts(counts)


def test_failure_counts_are_not_discarded_or_success_renormalized():
    counts = _model_counts(2)
    initial = p.infer_counts(counts)
    counts['start_failure', 'seed_timeout'] += sum(counts.values())
    changed = p.infer_counts(counts)
    assert changed['observed_mean'] == initial['observed_mean']/2
    assert changed['status'] == 'incompatible'


def test_exact_admission_counts_every_native_service_and_has_nontrivial_window():
    result = p.admission()
    budget = result['budget']
    assert result['all_window_ages_uniquely_decodable_on_coverage']
    assert 2*result['total_mean_radius_upper'] < F(1, 10)*F(9, 10)**3
    assert budget['read_failure_upper'] == F(66, 2**32)
    assert budget['selector_timeout_per_attempt'] == F(1, 2**32)
    assert budget['conditional_score_bias_upper'] == F(7, 2)*F(66, 2**32)+10*F(1, 2**32)
    assert budget['total_native_attempts_upper'] == p.TRIALS*(5040+4+160)+30
    assert budget['coherent_source_pulses'] == budget['auxiliary_reset_calls'] == 0
    larger = p.admission(1 << 18, 12)
    assert larger['all_window_ages_uniquely_decodable_on_coverage']
    assert not p.admission(p.TRIALS, 64)['all_window_ages_uniquely_decodable_on_coverage']
    r = result['statistical_radius_upper']
    assert 2*p.TRIALS*r*r >= 9*(p.FAILURE_BITS+1)


def test_random_window_can_have_the_same_mean_and_different_mean_duration():
    witness = p.random_window_obstruction()
    assert witness['mean_score'] == F(9, 10)
    assert witness['mean_attempt_count'] != witness['same_mean_as_fixed_age']
    assert not witness['equality_of_complete_record_distributions_claimed']
    assert not witness['unknown_random_or_adaptively_stopped_windows_identified']


def test_low_statistics_return_unresolved_instead_of_a_guessed_integer():
    counts = dict.fromkeys(p.alphabet(), 0)
    counts['read', 0, 1, 1] = 2
    assert p.infer_counts(counts)['status'] == 'unresolved'


@pytest.fixture(scope='module')
def receipt():
    return p.strict_json.load(p.OUTPUT)


def test_returned_receipt_does_not_expose_the_scope_contract(receipt, monkeypatch):
    by_law = {r['law']: r['execution'] for r in receipt['alternatives']}
    def execution(age, *, law):
        sample = deepcopy(by_law[law])
        sample['records_and_counts'] = [(tuple(r), n) for r, n in sample['records_and_counts']]
        return sample
    monkeypatch.setattr(p, 'execute_trials', execution)
    expected = dict(p.SCOPE)
    returned = p.run()
    returned['scope']['decoder_receives_interval_attempt_count'] = True
    assert p.SCOPE == expected
    assert p.run()['scope'] == expected


@pytest.fixture
def replayed_receipt(receipt, monkeypatch):
    # Isolate receipt comparison from the large source execution. CI separately
    # replays every default trial through the real CLI, without this fixture.
    monkeypatch.setattr(p, 'run', lambda: deepcopy(receipt))
    return receipt


def test_receipt_verifies_against_owned_source_replay(replayed_receipt, monkeypatch):
    calls = []
    def replay():
        calls.append(True)
        return deepcopy(replayed_receipt)
    monkeypatch.setattr(p, 'run', replay)
    assert p.compare_receipt(replayed_receipt, replayed_receipt)
    assert calls == [True]


@pytest.mark.parametrize('mutation', ['scope', 'proof', 'hash', 'missing_failure', 'count', 'states',
                                     'prefix', 'native_cost', 'inference', 'auxiliary', 'unknown_field'])
def test_receipt_rejects_hostile_contracts_even_against_itself(replayed_receipt, mutation):
    bad = deepcopy(replayed_receipt)
    row = bad['alternatives'][0]
    sample = row['execution']
    if mutation == 'scope': bad['scope']['decoder_receives_interval_attempt_count'] = True
    elif mutation == 'proof': bad['source_proof']['eigenvalue'] = '1'
    elif mutation == 'hash': bad['source_sha256']['model/native_clock/passive_age.py'] = '0'*64
    elif mutation == 'missing_failure': sample['records_and_counts'].pop(0)
    elif mutation == 'count': sample['records_and_counts'][0][1] = True
    elif mutation == 'states': sample['terminal_state_counts'][0][2] += 1
    elif mutation == 'prefix': sample['prefix'][0]['start_currents'].pop()
    elif mutation == 'native_cost': sample['auxiliary_native_attempts'] = -1
    elif mutation == 'inference': row['inference']['candidate_attempt_counts'] = [0]
    elif mutation == 'auxiliary': sample['final_auxiliary_state'] = [0]*12
    else: bad['physical_seconds'] = True
    with pytest.raises(ValueError):
        p.compare_receipt(bad, bad)


@pytest.mark.parametrize('mutation', ['impossible_cell', 'auxiliary_identity', 'main_identity',
    'discovery_identity', 'passive_identity', 'prefix_initial', 'prefix_continuation',
    'prefix_endpoint', 'execution_field', 'alternative_field', 'prefix_field'])
def test_receipt_rejects_fabricated_source_support_and_resource_identities(replayed_receipt, mutation):
    bad = deepcopy(replayed_receipt)
    row = bad['alternatives'][1 if mutation == 'prefix_endpoint' else 0]
    sample = row['execution']
    if mutation == 'impossible_cell':
        support = p._terminal_support(row['validation_only_interval_attempts'])
        for cell in sample['terminal_state_counts']:
            if len(support[cell[0]]) < 66:
                cell[1] = next(j for j in range(66) if j not in support[cell[0]])
                break
        else:
            raise AssertionError('source support did not exclude any impossible terminal cell')
    elif mutation == 'auxiliary_identity': sample['auxiliary_native_attempts'] += 1
    elif mutation == 'main_identity': sample['main_native_attempts'] += 1
    elif mutation == 'discovery_identity': sample['discovery_calls'] -= 1
    elif mutation == 'passive_identity':
        sample['passive_main_attempts'] -= 1
        sample['main_native_attempts'] -= 1
        sample['auxiliary_native_attempts'] -= 5
    elif mutation == 'prefix_initial':
        # This is a valid full-sector discovery from a different carrier.
        other = p.discovery.execute((0, 2)+(0,)*10, p._coin_stream(random.Random(9)), full=True)
        sample['prefix'][0]['start_currents'] = p._json(other.outcome.currents)
    elif mutation == 'prefix_continuation':
        initial = tuple(sample['prefix'][0]['terminal_state'])
        different = next(x for x in p.discovery.native.states(binary=True) if x != initial)
        other = p.discovery.execute(different, p._coin_stream(random.Random(9)))
        sample['prefix'][1]['start_currents'] = p._json(other.outcome.currents)
    elif mutation == 'prefix_endpoint':
        # A complete valid read is insufficient if no passive source path can
        # connect it to this trial's actual start configuration.
        first = sample['prefix'][0]
        start = p.decode_current_tape(p._tape(first['start_currents']), full=True)
        reachable = p._reachable(start.final, row['validation_only_interval_attempts'])
        different = next(x for x in p.discovery.native.states(binary=True) if x not in reachable)
        other = p.discovery.execute(different, p._coin_stream(random.Random(9)))
        first['end_currents'] = p._json(other.outcome.currents)
        first['record'] = p._json(p.endpoint_record(start, other.outcome))
        first['terminal_state'] = list(other.output)
    elif mutation == 'execution_field': sample['known_age'] = 3
    elif mutation == 'alternative_field': row['known_age'] = 3
    else: sample['prefix'][0]['known_age'] = 3
    with pytest.raises(ValueError):
        p.compare_receipt(bad, bad)


@pytest.mark.parametrize('mutation', ('stream_hash', 'auxiliary_endpoint', 'discovery_cost'))
def test_matching_fabrications_are_not_evidence_of_source_execution(replayed_receipt, mutation):
    bad = deepcopy(replayed_receipt)
    sample = bad['alternatives'][0]['execution']
    if mutation == 'stream_hash':
        sample['source_record_sha256'] = '0'*64
    elif mutation == 'auxiliary_endpoint':
        sample['final_auxiliary_state'] = [0]*10+[1, 1]
    else:
        sample['discovery_native_attempts'] = sample['discovery_calls']+30
        sample['main_native_attempts'] = sample['discovery_native_attempts']+sample['passive_main_attempts']
    # These internally consistent records satisfy the structural checks; only
    # the module-owned execution distinguishes them from actual observations.
    p._validate(bad)
    with pytest.raises(ValueError):
        p.compare_receipt(bad, bad)
    with pytest.raises(ValueError):
        p.compare_receipt(replayed_receipt, bad)
