"""Composition, actual preparation and hostile receipt controls."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest

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
