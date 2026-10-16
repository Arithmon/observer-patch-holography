from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from ..observer import native
from ..native_quantum_history import construction as quantum
from ..observer_reconstruction import native_geometry as ng
from ..primitive_observer import transfer
from . import binary_dynamics as d


@pytest.fixture(scope='module')
def source():
    kernel = transfer.native_counts()
    calibration = d.execute_calibration(kernel.seams)
    recovered = d.recover_response(calibration)
    response = [[F(x) for x in row] for row in recovered['response_exact']]
    return kernel, calibration, recovered, response


def test_current_only_binary_calibration_recovers_the_existing_geometry(source):
    kernel, calibration, recovered, _ = source
    assert calibration['success']
    assert len(calibration['current_records']) == 578
    assert all(len(row) == 4 for row in calibration['current_records'])
    assert all(d._loads(pair) in kernel.states for pair in d.PREPARATIONS)
    expected = ng.recover(ng.execute_calibration())
    assert recovered['gram_qsqrt5'] == expected['gram_qsqrt5']
    assert recovered['laplacian_exact'] == expected['laplacian_exact']
    assert not recovered['probability_or_coin_records_used']
    assert recovered['declared_fair_coin_law_used']


def test_current_frequencies_are_not_silently_used_to_fit_fairness(source):
    _, calibration, recovered, _ = source
    repeated = deepcopy(calibration)
    repeated['current_records'] = [row for row in calibration['current_records']
                                  for _ in range(1 if row[-1] == 0 else 7)]
    assert d.recover_response(repeated)['response_exact'] == recovered['response_exact']


def test_timeout_is_explicit_and_refuses_geometry(source):
    kernel, _, _, _ = source
    failed = d.execute_calibration(kernel.seams, cap=3, coin_source=lambda: 0)
    assert not failed['success']
    assert len(failed['timeout_settings']) == 110
    assert len(failed['current_records']) == 250+110*3
    with pytest.raises(ValueError, match='timed-out'):
        d.recover_response(failed)


@pytest.mark.parametrize('cap', [1, 2, 3, 4, 5])
def test_timeout_formula_against_all_fair_words_on_every_crossing_setting(source, cap):
    kernel, _, _, _ = source
    crossing, fixed = 0, 0
    for pair in d.PREPARATIONS:
        state = d._loads(pair)
        for seam in kernel.seams:
            delta = state[seam[0]]-state[seam[1]]
            if delta == 0:
                fixed += 1
                continue
            crossing += 1
            failures = 0
            for word in product((0, 1), repeat=cap):
                seen = {native.native(state, seam, coin)[1][1] for coin in word}
                failures += seen != {0, delta}
            assert F(failures, 2**cap) == F(2, 2**cap)
    certificate = d.stopping_certificate(kernel.seams, cap)
    assert certificate['crossing_settings'] == crossing == 110
    assert certificate['constant_current_settings'] == fixed == 250


@pytest.mark.parametrize('mutation', [
    lambda rows: rows[0].append(0),
    lambda rows: rows[0].__setitem__(0, True),
    lambda rows: rows[0].__setitem__(3, 7),
    lambda rows: rows.__setitem__(slice(None), [row for row in rows if row[0] != 0]),
])
def test_nonrecord_fields_bad_currents_and_missing_preparations_fail(source, mutation):
    calibration = deepcopy(source[1])
    mutation(calibration['current_records'])
    with pytest.raises(ValueError):
        d.recover_response(calibration)


def test_missing_one_current_support_point_cannot_be_replaced_by_empirical_frequency(source):
    calibration = deepcopy(source[1])
    key = next(tuple(row[:3]) for row in calibration['current_records'] if row[-1] != 0)
    calibration['current_records'] = [row for row in calibration['current_records']
                                      if tuple(row[:3]) != key or row[-1] == 0]
    with pytest.raises(ValueError, match='both possible current outcomes'):
        d.recover_response(calibration)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_exact_occupation_closure_matches_actual_quantum_kraus_compiler_on_coherent_state(source, law):
    kernel, _, _, response = source
    rng = np.random.default_rng(2113)
    psi = rng.normal(size=66)+1j*rng.normal(size=66)
    psi /= np.linalg.norm(psi)
    instrument = quantum.native_instrument(law)
    output_populations = sum((np.abs(kraus@psi)**2 for group in instrument.kraus for kraus in group),
                             np.zeros(66))
    occupations = np.asarray(kernel.states, float)
    actual = output_populations@occupations
    expected = np.asarray(response, float)@(np.abs(psi)**2@occupations)
    assert np.max(np.abs(actual-expected)) < 2e-14
    assert np.linalg.norm(np.outer(psi, psi.conj())-np.diag(np.abs(psi)**2)) > .1


def test_exact_GNS_covariance_and_strict_native_contraction(source):
    kernel, _, _, response = source
    result = d.heisenberg_and_fluctuation(kernel, response)
    assert result['constant_free_GNS_dimension'] == 11
    assert result['slow_response_GNS_subspace_dimension'] == 3
    assert [row['native_attempt_lag'] for row in result['lag_checks']] == [0, 1, 2, 4, 8]
    witness = result['explicit_occupation_contraction']
    assert F(witness['initial_GNS_norm_squared']) == F(10, 33)
    assert F(witness['one_attempt_GNS_norm_squared']) == F(1459, 5940)
    assert not result['global_no_wave_or_no_quantum_field_theorem']


def test_causal_response_front_matches_full_native_configuration_propagation(source):
    kernel, _, _, response = source
    witness = d.propagation_certificate(response, kernel.seams)['binary_sector_preparation_contrast']
    observed = []
    for name in ('preparation_A', 'preparation_B'):
        initial = d._loads(witness[name])
        law = [F(state == initial) for state in kernel.states]
        means = []
        for _ in range(4):
            means.append(sum(p*state[witness['receiver_port']] for p, state in zip(law, kernel.states)))
            law = [sum(law[i]*F(kernel.counts[i][j], 60) for i in range(66)) for j in range(66)]
        observed.append(means)
    difference = [a-b for a, b in zip(*observed)]
    assert difference == list(map(F, witness['occupation_mean_difference_by_attempt_0_to_3']))
    assert difference == [F(0), F(0), F(1, 1800), F(41, 27000)]


def test_complete_binary_dynamics_receipt():
    result = d.run()
    assert result['all_pass']
    assert result['scope']['main_calibration_and_dynamic_states_in_binary_66_state_sector']
    assert result['scope']['integer_calibration_used_only_as_external_comparison_control']
    assert result['calibration_stopping_certificate']['probability_of_any_timeout_upper_bound'] == '55/4194304'
    assert not result['scope']['modular_parameter_identified_with_native_attempts']
