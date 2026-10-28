"""Composition, actual preparation and hostile receipt controls."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest
from scipy.sparse import coo_matrix

from . import observer as o


@pytest.fixture(scope='module')
def replay():
    return o.report()


def test_complete_receipt_replays(replay):
    o._compare(o.OUTPUT, replay)


def test_full_channel_preserves_trace_and_untouched_quantum_reference():
    p = replace(o.mixing.canonical_parameters(), c_weight=F(1, 2))
    channel = o.source_channel(p)
    # A complex Bell pair, so a classical-only implementation cannot pass.
    v = np.zeros((66, 2), complex)
    v[0, 0], v[5, 1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    rho = np.einsum('ia,jb->ijab', v, v.conj())
    output = (channel@rho.reshape(66**2, 4)).reshape(66, 66, 2, 2)
    assert np.allclose(np.einsum('iiab->ab', output), np.eye(2)/2, atol=1e-13, rtol=0)
    assert np.linalg.norm(output[:, :, 0, 1]) > .01
    # Independent complete branch execution on every reference block.
    expected = np.zeros_like(rho)
    for branch in o.continuing.branches(p):
        for i, j in branch.pairs:
            for k, l in branch.pairs:
                expected[j, l] += float(branch.weight)*rho[i, k]
    assert np.allclose(output, expected, atol=1e-14, rtol=0)


def test_all_matrix_units_against_endpoint_oracle_and_wrong_cptp_control():
    p = replace(o.mixing.canonical_parameters(), c_weight=F(3, 7),
                repair_probability=F(1, 16), seam_weights=(1, 16)*15)
    states = o.local.basis(); index = {state: i for i, state in enumerate(states)}
    rows = list(range(66**2)); columns = rows.copy()
    values = [float(1-p.repair_probability)]*66**2
    # Independent floor/ceiling endpoint rule and original C-current/E-contrast
    # partition: no continuing branches, raw_maps or producer channel kernel.
    for name, mixture in (('C', p.c_weight), ('E', 1-p.c_weight)):
        for (a, b), selection in zip(o.local.SEAMS, p.seam_weights, strict=True):
            for coin in (0, 1):
                groups = {}
                for i, state in enumerate(states):
                    target = list(state); total = state[a]+state[b]
                    target[a], target[b] = (total+coin)//2, (total+1-coin)//2
                    label = state[a]-target[a] if name == 'C' else state[a]-state[b]
                    groups.setdefault(label, []).append((i, index[tuple(target)]))
                weight = (p.repair_probability*mixture*F(selection, sum(p.seam_weights))
                          *F(p.port_weights[a if coin else b],
                             p.port_weights[a]+p.port_weights[b]))
                for pairs in groups.values():
                    for i, j in pairs:
                        for k, l in pairs:
                            rows.append(66*j+l); columns.append(66*i+k)
                            values.append(float(weight))
    expected = coo_matrix((values, (rows, columns)), shape=(66**2, 66**2)).tocsr()
    actual = o.source_channel(p)
    difference = actual-expected
    assert max(abs(difference.data), default=0.) < 3e-15
    # Pre-dephasing is CPTP and has exactly the same population transition,
    # stationary reference and eventual replacement target. It still loses
    # native finite coherences and must fail the complete operator comparison.
    population_columns = np.eye(66).reshape(-1)
    wrong = actual.multiply(population_columns[None, :])
    assert np.count_nonzero((wrong-actual).multiply(population_columns[None, :]).data) == 0
    assert max(abs((wrong-expected).data)) > .9


def test_noncommuting_prepared_reference_log_and_unitary_bounds():
    budget = o.composition_budget({'complete_instrument_halfdiamond_upper': '0',
                                   'channel_norm_upper': '1'})
    p = o.mixing.canonical_parameters()
    pi = o.local.reference(p)
    # Unequal eigenvalues are essential: an off-diagonal perturbation within
    # one degenerate population sector would not challenge noncommutativity.
    high, low = float(max(pi)), float(min(pi))
    delta = float(F(budget['preparation_half_trace_upper']))
    target = np.diag([high, low]).astype(complex)
    actual = target+np.array([[0., 1j*delta], [-1j*delta, 0.]])
    assert np.linalg.norm(actual@target-target@actual, 2) > 1e-8
    assert abs(np.linalg.norm(actual-target, 'nuc')/2-delta) < 1e-20
    epsilon = float(F(budget['relative_operator_sandwich_epsilon']))
    assert np.linalg.eigvalsh(actual-(1-epsilon)*target).min() >= -1e-17
    assert np.linalg.eigvalsh((1+epsilon)*target-actual).min() >= -1e-17
    eigenvalues, frame = np.linalg.eigh(actual)
    actual_log = (frame*np.log(eigenvalues))@frame.conj().T
    log_bound = float(F(budget['modular_unitary_error_per_abs_parameter_upper']))
    assert 1e-6 < np.linalg.norm(actual_log-np.diag(np.log([high, low])), 2) <= log_bound
    for parameter in (-1., -.5, .5, 1.):
        actual_unitary = (frame*eigenvalues.astype(complex)**(-1j*parameter))@frame.conj().T
        target_unitary = np.diag(np.array([high, low], complex)**(-1j*parameter))
        assert 1e-6 < np.linalg.norm(actual_unitary-target_unitary, 2) <= abs(parameter)*log_bound


def test_actual_long_preparation_and_resource_accounting(replay):
    run = replay['preparation_execution']
    assert run['repair_attempts']+run['FIFO_filling_attempts'] == 1601
    assert run['complete_raw_slots_per_attempt'] == 721
    resources = replay['record_resources_through_two_repairs_and_detector']
    assert resources['exported_records'] == 1603
    assert resources['retained_records'] == 1
    assert resources['complete_archive_bits_upper'] == 16040
    assert resources['fifo_record_bits_upper'] == 10


def test_positive_reference_stability_and_budget_composition(replay):
    b = replay['composition_budget']
    d = F(b['preparation_half_trace_upper'])
    minimum = F(b['reference_minimum_eigenvalue'])
    eps = F(b['relative_operator_sandwich_epsilon'])
    assert minimum == F(64, 7505) and eps == d/minimum < F(1, 10000)
    assert F(b['modular_unitary_error_per_abs_parameter_upper']) == eps/(1-eps)
    assert d < F(1, 10**6)
    assert F(b['complete_nominal_clock_readout_half_trace_upper']) < F(1, 10000)
    assert F(b['complete_actual_source_reference_clock_comparison_at_s1_upper']) < F(1, 1000)
    assert not b['physical_service_error_included']


def test_sandwich_survives_rare_record_without_postselection_normalization():
    p = o.mixing.canonical_parameters()
    pi = np.array(list(map(float, o.local.reference(p))))
    perturbation = np.zeros((66, 66), complex)
    perturbation[0, 1], perturbation[1, 0] = 1j*1e-7, -1j*1e-7
    epsilon = 1e-7/min(pi)
    rho = np.diag(pi)+perturbation
    for branch in o.continuing.branches(p):
        if not branch.weight or not branch.pairs:
            continue
        indices = [i for i, _ in branch.pairs]
        actual = float(branch.weight)*rho[indices][:, indices]
        target = float(branch.weight)*np.diag(pi[indices])
        assert np.linalg.eigvalsh(actual-(1-epsilon)*target).min() >= -1e-15
        assert np.linalg.eigvalsh((1+epsilon)*target-actual).min() >= -1e-15


@pytest.mark.parametrize('mutation', ['empty', 'false_is_zero', 'uncounted_exports', 'changed_error', 'missing_source', 'extra_claim', 'nan'])
def test_receipt_rejects_trash_and_scope_promotions(replay, tmp_path, mutation):
    value = deepcopy(replay)
    if mutation == 'empty': value = {}
    elif mutation == 'false_is_zero': value['scope']['weights_source_selected'] = 0
    elif mutation == 'uncounted_exports': value['record_resources_through_preparation']['exported_records'] = 0
    elif mutation == 'changed_error': value['composition_budget']['preparation_half_trace_upper'] = '0'
    elif mutation == 'missing_source': value['source_sha256'].pop(next(iter(value['source_sha256'])))
    elif mutation == 'extra_claim': value['scope']['universal_physical_theory'] = True
    elif mutation == 'nan': value['composition_budget']['preparation_half_trace_upper'] = float('nan')
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(value), encoding='utf-8')
    with pytest.raises(ValueError): o._compare(path, replay)


def test_duplicate_keys_and_changed_frozen_contract_fail(tmp_path):
    path = tmp_path/'receipt.json'
    path.write_text('{"scope":false,"scope":false}', encoding='utf-8')
    with pytest.raises(ValueError): o._compare(path, {})
    path.write_text(o.CONTRACT.read_text().replace('1600', '16'), encoding='utf-8')
    with pytest.raises(ValueError): o.contract(path)


@pytest.mark.parametrize('error,gain', [('-1', '1'), ('0', '0'), ('NaN', '1'),
                                      (True, '1'), ('0/2', '1'), ('0', '1.0')])
def test_composition_refuses_invalid_scientific_bounds(error, gain):
    with pytest.raises(ValueError):
        o.composition_budget({'complete_instrument_halfdiamond_upper': error,
                              'channel_norm_upper': gain})
