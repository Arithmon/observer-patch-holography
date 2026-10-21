"""Independent native maps, sharp witnesses and rational preparation floors."""
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import coherence_budget as budget


def distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def endpoint_maps(law, seam, states):
    """Direct endpoint averaging, independent of native.native/fixed_maps."""
    index = {state: i for i, state in enumerate(states)}
    a, b = seam
    groups = defaultdict(dict)
    for i, state in enumerate(states):
        total = state[a]+state[b]
        for coin in (0, 1):
            output = list(state)
            output[a] = total//2+coin*(total % 2)
            output[b] = total-output[a]
            current = state[a]-output[a]
            label = current if law == 'C' else state[a]-state[b]
            groups[coin, label][i] = index[tuple(output)]
    return groups


@pytest.fixture(scope='module')
def reference():
    return budget.prepared_reference_floor()


def test_actual_native_partial_injections_match_independent_endpoint_census():
    native = budget.native_interface
    states = native.quantum.STATES
    for law in ('C', 'E'):
        for seam in native.quantum.mixing.seams():
            independent = endpoint_maps(law, seam, states)
            actual = {(branch.coin, branch.label): dict(branch.mapping)
                      for branch in native.fixed_maps(law, seam)}
            assert actual == independent
            assert all(len(set(mapping.values())) == len(mapping)
                       for mapping in independent.values())


def test_adaptive_branch_products_retain_the_dephasing_identity():
    native = budget.native_interface
    states = native.quantum.STATES
    seams = native.quantum.mixing.seams()
    rng = np.random.default_rng(20261021)
    vector = rng.normal(size=len(states))+1j*rng.normal(size=len(states))
    vector /= np.linalg.norm(vector)
    rho = np.outer(vector, vector.conj())
    diagonal = np.diag(np.diag(rho))
    for law in ('C', 'E'):
        for (coin, label), first in endpoint_maps(law, seams[0], states).items():
            # The adaptive seam uses only committed labels. Some branches
            # stop early; others are composed without branch normalization.
            following = ({i: i for i in range(len(states))},) if label < 0 else tuple(
                endpoint_maps(law, seams[1+coin], states).values())
            for second in following:
                pairs = [(i, second[j]) for i, j in first.items() if j in second]
                matrix = np.zeros((len(states), len(states)), complex)
                for i, j in pairs:
                    matrix[j, i] = 1/np.sqrt(2) if label < 0 else .5
                output = matrix@rho@matrix.conj().T
                other = matrix@diagonal@matrix.conj().T
                np.testing.assert_allclose(np.diag(np.diag(output)), other,
                                           atol=1e-16, rtol=0)
                assert np.trace(output).real == pytest.approx(np.trace(other).real, abs=1e-15)


def test_native_sign_targets_have_sharp_half_distance_and_complete_five_over_66():
    row = budget.native_phase_floor()
    native = budget.native_interface
    pairs = native.phase_pairs()
    size = len(native.quantum.STATES)
    diagonal = np.zeros((size, size))
    targets = {}
    for sign in (-1, 1):
        target = np.zeros((size, size))
        witness = np.zeros((size, size))
        for a, b in pairs:
            target[a, a] = target[b, b] = 1/(2*len(pairs))
            target[a, b] = target[b, a] = sign/(2*len(pairs))
            witness[a, b] = witness[b, a] = sign
            diagonal[a, a] = diagonal[b, b] = 1/(2*len(pairs))
        targets[sign] = target
        assert distance(target, diagonal) == pytest.approx(.5, abs=1e-15)
        assert np.count_nonzero(np.diag(witness)) == 0
        assert np.max(abs(np.linalg.eigvalsh(witness))) == 1
        # This zero-diagonal witness proves the lower bound against every
        # diagonal output, even with a different (unnormalized) success mass.
        assert np.trace(witness@target) == pytest.approx(1, abs=1e-15)
    assert row['conditional_distance_to_diagonal'] == F(1, 2)
    assert row['complete_herald_work_distance_floor'] == F(5, 66)

    for law in ('C', 'E'):
        branches = endpoint_maps(law, native.quantum.mixing.seams()[0], native.quantum.STATES)
        original = np.eye(size)/size
        first = {sign: np.zeros((size, size)) for sign in (-1, 0, 1)}
        for (coin, label), mapping in branches.items():
            current = label if law == 'C' else (label-(2*coin-1)*(abs(label) % 2))//2
            for i, j in mapping.items():
                first[current][j, j] += original[i, i]/2
        total = 0.
        for sign in (-1, 1):
            actual = np.zeros((size, size))
            for mapping in branches.values():
                for i, j in mapping.items():
                    actual[j, j] += first[sign][i, i]/2
            target = float(row['per_sign_success_probability'])*targets[sign]
            total += distance(target, actual)
        # Failed heralds are unchanged, and the extra repair is unread.
        assert total == pytest.approx(float(F(5, 66)), abs=2e-15)


def test_reference_enclosure_matches_independent_matrix_eigenvalues_and_twirl(reference):
    laplacian = np.zeros((12, 12))
    for a, b in budget.native_interface.quantum.mixing.seams():
        vector = np.zeros(12); vector[a], vector[b] = 1, -1
        laplacian += np.outer(vector, vector)
    lazy = np.eye(12)-laplacian/120
    matrix = sum(float(weight)*np.linalg.matrix_power(lazy, 2*n)
                 for n, weight in enumerate(reference['weights']))
    success = np.trace(matrix)/12
    rho = matrix/np.trace(matrix)
    uniform = np.eye(12)/12
    actual = distance(rho, uniform)
    lower, upper = reference['conditional_distance_to_diagonal_interval']
    assert float(lower)-2e-15 <= actual <= float(upper)+2e-15
    assert success == pytest.approx(float(reference['success_probability']), abs=2e-15)
    residual_eigenvalues = np.linalg.eigvalsh(rho-uniform)
    assert np.count_nonzero(residual_eigenvalues > 1e-12) == 4
    assert np.count_nonzero(residual_eigenvalues < -1e-12) == 8

    rng = np.random.default_rng(714)
    for _ in range(30):
        diagonal = np.diag(rng.dirichlet(np.ones(12)))
        assert distance(rho, diagonal) >= actual-2e-15
    average = np.zeros((12, 12))
    for permutation in budget.source().rotations():
        transform = np.eye(12)[:, permutation]
        average += transform@diagonal@transform.T/60
        np.testing.assert_allclose(transform@rho@transform.T, rho, atol=2e-16, rtol=0)
    np.testing.assert_allclose(average, uniform, atol=2e-16, rtol=0)
    assert distance(rho, average) <= distance(rho, diagonal)+2e-15


def test_reference_spectral_signs_and_sqrt_bounds_are_exact(reference):
    lower, upper = reference['sqrt5_interval']
    assert lower**2 < 5 < upper**2
    p = reference['success_probability']
    masses = reference['spectral_mass_intervals']
    assert masses['0'][0] > p and masses['5-sqrt5'][0] > p
    assert masses['6'][1] < p and masses['5+sqrt5'][1] < p
    lo, hi = reference['complete_success_failure_distance_floor_interval']
    assert F(4494803, 10**9) < lo <= hi < F(4494804, 10**9)
    assert (lo/p, hi/p) == reference['conditional_distance_to_diagonal_interval']
    wider = budget.prepared_reference_floor(3)
    wlo, whi = wider['complete_success_failure_distance_floor_interval']
    assert wlo <= lo <= hi <= whi


def test_changed_success_probabilities_do_not_evade_complete_flag_floor(reference):
    # Independent two-level example with target |+>, arbitrary diagonal
    # conditional outputs, and every failed-trial probability retained.
    rho = np.ones((2, 2))/2
    for p in (F(1), F(1, 3), F(1, 1000000)):
        for q in (0., float(p)/2, float(p), .9):
            for population in (0., .3, .5, 1.):
                delta = np.diag([population, 1-population])
                complete = distance(float(p)*rho, q*delta)+abs(float(p)-q)/2
                assert complete >= float(budget.herald_floor(p, F(1, 2)))-1e-15
        assert budget.herald_floor(p, F(1, 2))/p == F(1, 2)
        assert budget.herald_floor(p, F(1, 2), F(1, 10)) == p*F(2, 5)
    floor = reference['complete_success_failure_distance_floor_interval'][0]
    assert budget.robust_floor(floor, target_error=floor/7, input_error=floor/7,
                               per_step_error=floor/14, attempt_cap=2) == 4*floor/7
    assert budget.robust_floor(floor, per_step_error=F(1), attempt_cap=2) == 0


def test_source_guard_arithmetic_matches_existing_receipt_without_claiming_new_guard_check(reference):
    allowance = budget.prepared_source_allowance(reference['weights'], reference['success_probability'])
    receipt = json.loads((budget.ROOT/'model/predictions/prepared_clock_results.json').read_text())
    assert allowance == F(receipt['error_budget']['complete_half_trace_upper_under_guards'])
    result = budget.run()
    comparison = result['existing_prepared_source_comparison']
    assert not comparison['source_word_numerical_guards_reexecuted_here']
    assert not comparison['physical_fidelity_or_directed_rounding_certificate']
    floor = reference['complete_success_failure_distance_floor_interval'][0]
    assert F(comparison['distance_from_native_only_outputs_lower_under_those_guards']) == floor-allowance
    assert not result['scope']['preparation_floor_applied_to_terminal_measured_detector_output']
    assert not result['scope']['uncommitted_coherent_record_dilations_are_incoherent']
    assert not result['scope']['diagonal_marginals_of_correlated_accessible_registers_suffice']
    assert all(len(value) == 64 for value in result['source_sha256'].values())


@pytest.mark.parametrize('bad', (True, 0, .5, '1/2', None, float('inf'), float('nan'), F(-1), F(2)))
def test_error_allowances_reject_nonoriginal_or_out_of_range_values(bad):
    with pytest.raises(ValueError):
        budget.herald_floor(bad, F(1, 2))
    with pytest.raises(ValueError):
        budget.robust_floor(F(1, 2), target_error=bad)


@pytest.mark.parametrize('bad', (True, 1., -1, 10**9+1))
def test_attempt_caps_are_original_bounded_integers(bad):
    with pytest.raises(ValueError):
        budget.robust_floor(F(1, 2), attempt_cap=bad)


@pytest.mark.parametrize('bad', (True, 9., 0, 31))
def test_enclosure_precision_is_bounded_and_exact(bad):
    with pytest.raises(ValueError):
        budget.sqrt5_interval(bad)


def test_native_premise_failure_is_not_replaced_by_an_assumed_incoherent_model(monkeypatch):
    original = budget.native_interface.exact_native_audit

    def changed(law):
        result = original(law)
        result['native_outcomes_commute_with_work_dephasing'] = False
        return result

    monkeypatch.setattr(budget.native_interface, 'exact_native_audit', changed)
    with pytest.raises(ArithmeticError, match='partial-injection premise'):
        budget.native_phase_floor()


def test_uncommitted_native_isometry_is_not_diagonal_preserving():
    # One occupied endpoint and one spectator. The two fair coin alternatives
    # are different joint work/record basis states. Keeping their amplitudes
    # gives a coherent pure state, unlike committing the record.
    states = budget.native_interface.quantum.STATES
    a, b = budget.native_interface.quantum.mixing.seams()[0]
    spectator = next(i for i in range(12) if i not in (a, b))
    initial = tuple(int(i in (a, spectator)) for i in range(12))
    maps = endpoint_maps('C', (a, b), states)
    source_index = states.index(initial)
    images = [(label, mapping[source_index]) for label, mapping in maps.items() if source_index in mapping]
    assert len(images) == len(set(images)) == 2
    coherent_joint = np.ones((2, 2))/2
    assert coherent_joint[0, 1] != 0
    assert distance(coherent_joint, np.diag(np.diag(coherent_joint))) == .5


def test_diagonal_ancillary_marginals_do_not_replace_joint_accessible_diagonality():
    # A basis permutation can turn preexisting joint coherence into local
    # coherence. This is outside the premise despite both initial marginals
    # being diagonal; the test does not claim CNOT is an available native map.
    bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
    before = np.outer(bell, bell)
    tensor = before.reshape(2, 2, 2, 2)
    np.testing.assert_allclose(np.trace(tensor, axis1=1, axis2=3), np.eye(2)/2)
    np.testing.assert_allclose(np.trace(tensor, axis1=0, axis2=2), np.eye(2)/2)
    permutation = np.eye(4)[:, [0, 1, 3, 2]]
    after = (permutation@before@permutation.T).reshape(2, 2, 2, 2)
    work = np.trace(after, axis1=1, axis2=3)
    assert work[0, 1] == pytest.approx(.5)


@pytest.fixture(scope='module')
def receipt():
    return budget.run()


@pytest.mark.parametrize('mutation', (
    'missing_key', 'extra_key', 'source_hash', 'exact_bound', 'fraction_spelling',
    'count_type', 'scope_value', 'scope_type', 'omitted_band', 'nan', 'infinity',
))
def test_cli_rejects_mutated_complete_receipts(tmp_path, receipt, mutation):
    changed = deepcopy(receipt)
    prepared = changed['prepared_reference_target']
    if mutation == 'missing_key':
        del changed['robustness']
    elif mutation == 'extra_key':
        changed['physical_fidelity_certified'] = True
    elif mutation == 'source_hash':
        changed['source_sha256'][next(iter(changed['source_sha256']))] = '0'*64
    elif mutation == 'exact_bound':
        key = 'complete_success_failure_distance_floor_interval'
        prepared[key][0] = str(F(prepared[key][0])+F(1, 10**100))
    elif mutation == 'fraction_spelling':
        prepared['success_probability'] = '0'+prepared['success_probability']
    elif mutation == 'count_type':
        prepared['positive_multiplicity'] = 4.0
    elif mutation in ('scope_value', 'scope_type'):
        changed['scope']['classical_adaptation_and_heralding_included'] = (
            False if mutation == 'scope_value' else 1)
    elif mutation == 'omitted_band':
        prepared['spectral_mass_intervals'].pop('6')
    else:
        prepared['success_probability'] = float('nan' if mutation == 'nan' else 'inf')
    path = tmp_path/'mutated.json'
    path.write_text(json.dumps(changed), encoding='utf-8')
    with pytest.raises(ValueError, match='Receipt'):
        budget.main(['--check', str(path)])


def test_cli_rejects_duplicate_keys_before_comparison(tmp_path, receipt):
    text = budget.canonical(receipt)
    text = text.replace('"positive_multiplicity": 4',
                        '"positive_multiplicity": 3, "positive_multiplicity": 4')
    path = tmp_path/'duplicate.json'
    path.write_text(text, encoding='utf-8')
    with pytest.raises(ValueError, match='duplicate JSON key'):
        budget.main(['--check', str(path)])


def test_cli_reproduces_saved_receipt_and_refuses_overwrite(tmp_path, receipt, capsys):
    saved = budget.ROOT/'model/native_interface/coherence_results.json'
    budget.main(['--check', str(saved)])
    assert capsys.readouterr().out == 'Native coherence-budget receipt reproduced\n'
    path = tmp_path/'fresh.json'
    budget.main(['--write', str(path)])
    assert path.read_text(encoding='utf-8') == budget.canonical(receipt)
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        budget.main(['--write', str(path)])
    assert path.read_bytes() == before
