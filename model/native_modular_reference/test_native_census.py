"""Native census observation, full quantum failures and empirical inference."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import native_census as c


def _counts(n=1 << 20):
    values = (64000, 88000, 7100)
    return {label: (values[i] if i < 3 else n-sum(values) if i == 3 else 0)
            for i, label in enumerate(c.LABELS)}


def test_source_supports_and_record_only_decoder():
    audit = c.exact_identification()
    assert audit['successful_original_states_checked'] == 66*2048
    assert audit['charge_effect_factor'] == (1-F(1, 2**31))**11
    assert set(audit['laws']) == {'C', 'E'}
    assert all(row['matrix_unit_effect_checks'] == 11*4*4356 for row in audit['laws'].values())
    assert not audit['controller_reads_hidden_occupations']


@pytest.mark.parametrize('state', [0, 1, 17, 65])
def test_decoder_reconstructs_changed_state_from_currents(state):
    initial = c.native.states(binary=True)[state]
    records, final = c.execute_record(initial, tuple(i % 2 for i in range(11)))
    assert c.decode(records) == (initial, final)
    assert sum(final) == 2
    bad = list(records)
    bad[0] = (1, 1, 1)
    with pytest.raises(ValueError):
        c.decode(tuple(bad))


@pytest.mark.parametrize('deadline', [2, 3, 32, 64])
def test_finite_population_mass_and_charge_effects(deadline):
    maps = c.population_maps(deadline)
    mass = c.success_probability(deadline)
    charges = [sum(x[p] for p in c.local.CORE) for x in c.native.states(binary=True)]
    for q in range(3):
        np.testing.assert_allclose(maps[f'success_q{q}'].sum(axis=0),
                                   [float(mass) if charge == q else 0 for charge in charges], atol=2e-15)
    np.testing.assert_allclose(maps['census_timeout'].sum(axis=0), float(1-mass), atol=2e-15)
    np.testing.assert_allclose(sum(maps.values()).sum(axis=0), 1, atol=2e-15)
    if deadline == 2:
        assert np.count_nonzero(maps['census_timeout']-np.eye(66)*float(1-mass)) > 0


def test_quantum_instrument_preserves_entangled_reference_and_all_failures():
    states = c.native.states(binary=True)
    selected = [next(i for i, x in enumerate(states) if sum(x[p] for p in c.local.CORE) == q)
                for q in (0, 1)]
    vector = np.zeros(132, complex)
    vector[2*selected[0]] = 1/np.sqrt(2)
    vector[2*selected[1]+1] = 1j/np.sqrt(2)
    rho = np.outer(vector, vector.conj())
    for law in ('C', 'E'):
        outputs = c.execute(rho, law=law, deadline=2)
        assert set(outputs) == {'success_q0', 'success_q1', 'success_q2', 'census_timeout'}
        assert sum(np.trace(m).real for m in outputs.values()) == pytest.approx(1)
        assert np.trace(outputs['census_timeout']).real == pytest.approx(float(1-c.success_probability(2)))
        reduced = sum(np.einsum('aras->rs', m.reshape(66, 2, 66, 2)) for m in outputs.values())
        np.testing.assert_allclose(reduced, np.eye(2)/2, atol=1e-14)
        np.testing.assert_allclose(np.einsum('aras->rs', outputs['success_q0'].reshape(66, 2, 66, 2)),
                                   np.diag([float(c.success_probability(2))/2, 0]), atol=1e-14)
        assert all(np.linalg.eigvalsh(m).min() > -1e-14 for m in outputs.values())


def test_quantum_diagonals_agree_with_finite_population_maps():
    initial = np.arange(1, 67, dtype=float)
    initial /= initial.sum()
    maps = c.population_maps(3)
    for law in ('C', 'E'):
        out = c.execute(np.diag(initial).astype(complex), law=law, deadline=3)
        for label in out:
            np.testing.assert_allclose(out[label], np.diag(maps[label]@initial), atol=2e-15)


def test_compressed_edge_equals_all_microscopic_short_deadline_records():
    """Independent fine-record replay checks compression on coherent inputs."""
    rng = np.random.default_rng(93)
    vector = rng.normal(size=66)+1j*rng.normal(size=66)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    seam, deadline = c.TREE[0], 4
    alpha = 1-F(1, 2**(deadline-1))
    def apply(mapping, rho):
        before, after = zip(*mapping, strict=True)
        output = np.zeros_like(rho)
        output[np.ix_(after, after)] = rho[np.ix_(before, before)]
        return output
    for law in ('C', 'E'):
        native_branches = c.interface.fixed_maps(law, seam)
        fine = {False: np.zeros_like(density), True: np.zeros_like(density)}
        for word in product((0, 1), repeat=deadline):
            stop = next((i for i in range(1, deadline) if word[i] != word[0]), None)
            prefix = word if stop is None else word[:stop+1]
            # Summing all full words gives each stopped prefix exactly its
            # probability, including the unexecuted suffix multiplicity.
            branches = [density]
            for coin in prefix:
                branches = [apply(branch.mapping, rho) for rho in branches
                            for branch in native_branches if branch.coin == coin]
            fine[stop is not None] += sum(branches)/2**deadline
        for complete in (False, True):
            compressed = sum((apply(mapping, density)
                for coin in (0, 1)
                for _, _, mapping in c._edge_branches(law, seam, coin, complete=complete)),
                np.zeros_like(density))
            compressed *= float(alpha if complete else 1-alpha)/2
            np.testing.assert_allclose(compressed, fine[complete], atol=2e-15)


@pytest.mark.parametrize('mutation', ['phase', 'dephasing', 'output'])
def test_native_kraus_audit_rejects_hidden_quantum_and_output_changes(monkeypatch, mutation):
    original = c.interface.fixed_maps
    def changed(law, seam):
        branches = list(original(law, seam))
        branch = next(b for b in branches if len(b.mapping) > 2)
        index = branches.index(branch)
        if mutation == 'phase':
            matrix = branch.matrix.copy()
            i, j = branch.mapping[0]
            matrix[j, i] *= -1
            branches[index] = replace(branch, matrix=matrix)
        elif mutation == 'dephasing':
            branches.pop(index)
            for mapping in (branch.mapping[:1], branch.mapping[1:]):
                matrix = np.zeros_like(branch.matrix)
                for i, j in mapping:
                    matrix[j, i] = 1/np.sqrt(2)
                branches.append(replace(branch, mapping=mapping, matrix=matrix))
        else:
            mapping = list(branch.mapping)
            (i, a), (j, b) = mapping[:2]
            mapping[:2] = ((i, b), (j, a))
            matrix = np.zeros_like(branch.matrix)
            for i, j in mapping:
                matrix[j, i] = 1/np.sqrt(2)
            branches[index] = replace(branch, mapping=tuple(mapping), matrix=matrix)
        return tuple(branches)
    monkeypatch.setattr(c.interface, 'fixed_maps', changed)
    with pytest.raises(ArithmeticError):
        c._source_edge_audit('C', c.TREE[0])
    with pytest.raises(ArithmeticError):
        c.exact_identification()


def test_empirical_radius_proves_positive_root_bound_and_is_target_blind(monkeypatch):
    expected = c.local.budget()
    monkeypatch.setattr(c.local, 'reference', lambda: (_ for _ in ()).throw(AssertionError('target lookup')))
    row = c.infer_counts(_counts())
    assert row['status'] == 'faithful'
    n, k = row['trials'], row['confidence_bits']
    epsilon = row['systematic_halftrace_upper']
    coefficient, intercept = F(2*k, n), epsilon+F(2*k, 3*n)
    for q, radius in enumerate(row['charge_probability_radii']):
        frequency = F(_counts()[f'success_q{q}'], n)
        assert radius >= intercept+coefficient/2
        assert (radius-intercept)**2 >= coefficient*(frequency+epsilon+radius)
        assert row['central_subnormalized_weights'][q] == frequency/c.DIMENSIONS[q]
    assert row['source_native_attempts_per_trial_upper'] == expected['total_native_attempts_upper']+353
    assert row['source_expected_native_attempts_per_trial_upper'] == expected['expected_total_native_attempts_upper']+34
    assert row['confidence_failure_upper'] == F(3, 512)
    assert sum(row['census_controller_register_bits'].values()) == row['census_controller_register_bits_subtotal'] == 139
    assert row['training_terminal_counter_bits'] == 168
    assert row['training_trial_counter_bits'] == 21
    assert row['memory_counts_are_logical_layout_not_measured_python_memory']
    assert not row['host_receipt_terminal_state_counts_are_observer_inputs']


@pytest.mark.parametrize('mutate', [
    lambda d: d.pop('census_timeout'), lambda d: d.update(extra=0),
    lambda d: d.update(success_q0=True), lambda d: d.update(success_q0=1.0),
    lambda d: d.update(success_q0=-1), lambda d: d.update({k: 0 for k in d}),
])
def test_incomplete_or_invalid_counts_rejected(mutate):
    counts = _counts()
    mutate(counts)
    with pytest.raises(ValueError):
        c.infer_counts(counts)


@pytest.mark.parametrize('bias', [True, -1, 2, 0.0, float('nan'), '0'])
def test_invalid_bias_is_not_a_certificate(bias):
    with pytest.raises(ValueError):
        c.infer_counts(_counts(), sampling_bias=bias)


def test_missing_rare_sector_cannot_be_clipped_to_a_clock():
    counts = _counts()
    counts['herald_failure'] += counts['success_q2']
    counts['success_q2'] = 0
    row = c.infer_counts(counts)
    assert row['status'] == 'unresolved'
    assert row['subnormalized_weight_intervals'][2][0] < 0
    assert 'central_weights' not in row


def test_stateful_residual_sampling_is_chunk_invariant():
    kernel = np.zeros((528, 66), dtype=np.uint64)
    for state in range(66):
        kernel[66+3, state] = c.GRID//2
        kernel[2*66+(state+1) % 66, state] += c.GRID//2
    first = c.sample_kernel(kernel, 1000, seed=77, initial=9, chunk_size=1)
    later = c.sample_kernel(kernel, 1000, seed=77, initial=9, chunk_size=137)
    assert first == later
    assert first['common_subkernel_probability'] == F(1, 2)
    assert 400 < first['residual_draws'] < 600
    deterministic = np.zeros((528, 66), dtype=np.uint64)
    for state in range(66):
        deterministic[6*66+(state+1) % 66, state] = c.GRID
    row = c.sample_kernel(deterministic, 137, initial=17, chunk_size=13)
    assert row['final_state'] == (17+137) % 66
    assert row['counts']['census_timeout'] == 137


@pytest.fixture(scope='module')
def execution():
    return c.execute_trials()


@pytest.fixture(scope='module')
def receipt():
    return c.old._json(c.run())


def test_executed_continuing_training_has_useful_precision(execution):
    sample, inferred = execution['execution'], execution['inference']
    assert sum(sample['counts'].values()) == c.TRAINING_TRIALS
    assert sum(map(sum, sample['terminal_state_counts'])) == c.TRAINING_TRIALS
    assert sample['same_carrier_reused'] and not sample['microscopic_attempt_tapes_executed']
    assert inferred['status'] == 'faithful'
    weights, intervals = inferred['central_subnormalized_weights'], inferred['subnormalized_weight_intervals']
    lower = [c.old.common.logarithm_interval(w/hi)[0] for w, (_, hi) in zip(weights, intervals, strict=True)]
    upper = [c.old.common.logarithm_interval(w/lo)[1] for w, (lo, _) in zip(weights, intervals, strict=True)]
    assert max(upper[i]-lower[j] for i in range(3) for j in range(3) if i != j) < F(8, 100)
    assert inferred['source_native_attempts_all_trials_upper'] < 7*10**14
    assert inferred['source_expected_native_attempts_all_trials_upper'] < 111*10**9


def test_both_native_laws_and_exact_sampler_columns():
    first, diagnostic = c.quantized_kernel('C')
    second, _ = c.quantized_kernel('E')
    np.testing.assert_array_equal(first, second)
    audit = c.audit_kernel(first)
    assert audit['exact_dyadic_columns_checked'] == 66
    assert audit['complete_total_variation_upper'] < c.SAMPLER_GATE
    assert all(x <= c.SAMPLER_GATE for x in audit['charge_mean_bias_upper'])
    assert 0 <= diagnostic['binary_population_mass_residual_diagnostic'] <= 5e-10


def test_receipt_proof_fields_cannot_be_rounded_or_self_approved(receipt):
    assert c.compare_receipt(receipt, receipt)
    for field in ('complete_total_variation_upper', 'charge_mean_bias_upper'):
        changed = deepcopy(receipt)
        changed['continuing_execution']['sampler_certificate'][field] = '0'
        assert not c.compare_receipt(changed, receipt)
        assert not c.compare_receipt(changed, changed)
    changed = deepcopy(receipt)
    changed['continuing_execution']['numerical_kernel']['binary_population_mass_residual_diagnostic'] = -1e-20
    assert not c.compare_receipt(changed, changed)
    for field in ('confidence_failure_upper', 'source_native_attempts_all_trials_upper'):
        changed = deepcopy(receipt)
        changed['continuing_execution']['inference'][field] = '0'
        assert not c.compare_receipt(changed, changed)
    changed = deepcopy(receipt)
    changed['continuing_execution']['execution']['terminal_state_counts'][0][0] += 1
    assert not c.compare_receipt(changed, changed)


@pytest.mark.parametrize('mutation', ['effect', 'matrix_units', 'scope', 'source', 'seed',
    'schema', 'hash', 'boolean', 'prefix', 'grid', 'validation'])
def test_full_census_self_verifier_rejects_false_proof_and_metadata(receipt, mutation):
    changed = deepcopy(receipt)
    row = changed['continuing_execution']
    if mutation == 'effect':
        changed['identification']['charge_effect_factor'] = '0'
    elif mutation == 'matrix_units':
        changed['identification']['laws']['C']['matrix_unit_effect_checks'] = 0
    elif mutation == 'scope':
        changed['scope']['finite_physical_service_accuracy_established'] = True
    elif mutation == 'source':
        changed['source_sha256'] = {}
    elif mutation == 'seed':
        row['execution']['seed'] = -1
    elif mutation == 'schema':
        changed['schema'] = 'false'
    elif mutation == 'hash':
        row['execution']['record_and_terminal_state_sha256'] = 'z'*64
    elif mutation == 'boolean':
        row['execution']['microscopic_attempt_tapes_executed'] = 0
    elif mutation == 'prefix':
        row['execution']['prefix'][0][1] = True
    elif mutation == 'grid':
        row['numerical_kernel']['grid'] -= 1
    else:
        row['validation_only_known_target_in_confidence_box'] = False
    assert not c.compare_receipt(changed, changed)


def test_cached_identification_is_not_mutable_by_its_callers():
    evidence = c.exact_identification()
    evidence['laws']['C']['matrix_unit_effect_checks'] = 0
    evidence['charge_effect_factor'] = F()
    fresh = c.exact_identification()
    assert fresh['laws']['C']['matrix_unit_effect_checks'] == 11*4*4356
    assert fresh['charge_effect_factor'] == c.success_probability()
