"""Actual-source, quantum-domain and failure controls for local bit inference."""
from fractions import Fraction as F
import json
import random

import numpy as np
import pytest
import sympy as sp

from . import bit_readout as b
from . import clock_bridge


def test_only_signed_records_enter_the_decoder_and_zeros_do_not_guess():
    assert b.decode([0]*20)['initial_bit'] is None
    result = b.decode([0, 1, 0, -1, 1])
    assert result == {'initial_bit': 0, 'current_bit': 1,
                      'resolved_at_local_tick': 2, 'local_ticks': 5}
    assert b.decode([-1, 0])['initial_bit'] == 1
    for word in ([1, 1], [-1, -1], [True], [0.0], [2], '1'):
        with pytest.raises(ValueError):
            b.decode(word)


def test_actual_native_records_reconstruct_initial_and_final_local_bits():
    kernel = b.transfer.native_counts()
    rng = random.Random(44771)
    successes, unresolved = 0, 0
    for initial in kernel.states:
        state = initial
        words = [[] for _ in range(12)]
        for attempt in range(80):
            seam = kernel.seams[rng.randrange(30)]
            output, (_, current) = b.transfer.native.native(state, seam, rng.randrange(2))
            # Observer-facing record uses orientation, not output occupations.
            words[seam[0]].append(-current)
            words[seam[1]].append(current)
            state = output
            if attempt in (0, 9, 79):
                for port, word in enumerate(words):
                    decoded = b.decode(word)
                    if decoded['current_bit'] is None:
                        unresolved += 1
                        assert all(j == 0 for j in word)
                    else:
                        successes += 1
                        assert decoded['initial_bit'] == initial[port]
                        assert decoded['current_bit'] == state[port]
    assert successes > 0 and unresolved > 0


def test_unsigned_ticks_do_not_identify_the_local_field():
    kernel = b.transfer.native_counts()
    seam = next(edge for edge in kernel.seams if 0 in edge)
    witnesses = {}
    for state in kernel.states:
        for coin in (0, 1):
            out, (_, current) = b.transfer.native.native(state, seam, coin)
            if current:
                witnesses.setdefault(state[0], (abs(current), out[0]))
    assert witnesses[0][0] == witnesses[1][0] == 1
    assert witnesses[0][1] != witnesses[1][1]


def test_complete_stopped_law_matches_independent_66_state_resolvent():
    kernel, blocks, _, small = b.source(0)
    full = clock_bridge.marked_source(kernel, 0)
    for block in blocks:
        for i in block:
            for dest_index, dest in enumerate(blocks):
                assert sum(full['Q0'][i, j] for j in dest) == small[0][blocks.index(block), dest_index]
                for sign in (-1, 1):
                    value = sum(full['Q1'][i, j] for j in dest
                                if kernel.states[j][0]-kernel.states[i][0] == sign)
                    assert value == small[sign][blocks.index(block), dest_index]


def test_resolution_keeps_failure_mass_and_conditional_sampling_changes_the_prior():
    kernel, blocks, _, q = b.source(0)
    one = sp.ones(len(blocks), 1)
    unresolved = one
    resolved = sp.zeros(len(blocks), 1)
    for n in range(16):
        resolved += q[0]**n*(q[-1]+q[1])*one
        unresolved = q[0]*unresolved
        assert resolved+unresolved == one
        assert min(unresolved) > 0
    assert max(unresolved) <= F(3, 8)
    prior = sp.Matrix([[sp.Rational(len(block), 66) for block in blocks]])
    bit = sp.Matrix([kernel.states[block[0]][0] for block in blocks])
    success = (q[-1]+q[1])*one
    # The initial bit among immediately resolved histories is not its prior 1/6.
    conditioned_initial = (prior*sp.diag(*list(bit))*success)[0]/(prior*success)[0]
    assert (prior*bit)[0] == F(1, 6)
    assert conditioned_initial == F(1, 2)


def test_both_quantum_laws_preserve_the_projector_domain_with_entangled_reference():
    kernel = b.transfer.native_counts()
    index = {state: i for i, state in enumerate(kernel.states)}
    rng = np.random.default_rng(119)
    psi = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    psi /= np.linalg.norm(psi)
    bit = np.array([state[0] for state in kernel.states])
    for law in ('C', 'E'):
        total = 0.0
        for edge in kernel.seams:
            for coin in (0, 1):
                groups = {}
                for i, state in enumerate(kernel.states):
                    output, (_, current) = b.transfer.native.native(state, edge, coin)
                    label = current if law == 'C' else state[edge[0]]-state[edge[1]]
                    groups.setdefault(label, []).append((i, index[output], current))
                for entries in groups.values():
                    out = np.zeros_like(psi)
                    for i, j, _ in entries:
                        out[j] = psi[i]/np.sqrt(60)
                    total += float(np.vdot(out, out).real)
                    current = entries[0][2]
                    local = -current if edge[0] == 0 else current if edge[1] == 0 else 0
                    if local:
                        expected = int(local == 1)
                        assert np.linalg.norm(out[bit != expected]) == 0
                        assert all(bit[i] == 1-expected for i, _, _ in entries)
        assert abs(total-1) < 2e-14
    assert b.quantum_support_audit(0)['C']['source_columns'] == 3960


def test_complete_record_interval_contains_the_original_field_without_postselection():
    kernel, blocks, _, q = b.source(0)
    one = sp.ones(len(blocks), 1)
    unknown = q[0]**16*one
    for block, failure in zip(blocks, unknown):
        initial = kernel.states[block[0]][0]
        resolved = F(1-failure)
        lo, hi = b.occupation_interval(resolved*(1-initial), resolved*initial, F(failure))
        assert lo <= initial <= hi
        assert hi-lo == failure
    for bad in ((F(1, 2), F(1, 2), F(1, 2)), (0, -1, 2), (0.0, 1, 0)):
        with pytest.raises(ValueError):
            b.occupation_interval(*bad)


def test_waiting_bounds_use_local_incidents_and_no_stationary_input():
    result = b.certificate()
    assert F(result['first_incident_success_min']) == F(2, 111)
    assert F(result['expected_local_ticks_max']) == F(129960, 8261)
    assert b.timeout_bound(128) == F(6561, 16777216) < F(1, 1000)
    assert 12*b.timeout_bound(160) < F(1, 1000)
    assert b.timeout_bound(15) == 1  # No invented fractional-block interpolation.
    for invalid in (-1, True, 1.0):
        with pytest.raises(ValueError):
            b.timeout_bound(invalid)


def test_complete_local_readout_intervals_feed_the_native_tensor_decoder():
    from ..common_source import null_tomography as tensor

    c = tensor.source()
    for pair in ((0, 3), (0, 1)):
        # All windows retrodict this common checkpoint. Each port has its own
        # sixteen-incident window; their finishing times are not synchronized.
        initial = tuple(int(port in pair) for port in range(12))
        law, means = tensor.native_record_law(initial, 8)
        intervals = []
        for port in range(12):
            kernel, blocks, _, q = b.source(port)
            unknown = q[0]**16*sp.ones(len(blocks), 1)
            masses = [F(0), F(0), F(0)]
            for block, failure in zip(blocks, unknown):
                weight = sum(law[i] for i in block)
                bit = kernel.states[block[0]][port]
                masses[bit] += weight*F(1-failure)
                masses[2] += weight*F(failure)
            interval = b.occupation_interval(*masses)
            assert interval[0] <= means[port] <= interval[1]
            intervals.append(interval)
        result = tensor.interval_tensor_estimate(intervals, native_attempts=8)
        assert not result['native_bound_excluded_by_intervals']
        assert not result['scope']['unknown_mean_exact_tensor_compatibility_asserted']
        coefficient_error = tuple(map(F, result['canonical_component_error_upper_Qsqrt5']))
        midpoint_tensor = [tuple(map(F, row)) for row in result['canonical_tensor_estimate_Qsqrt5']]
        true_tensor = [c.q5(0)]+[tensor._dot(row, [c.q5(x) for x in means])
                                for row in tensor._null_tomography_decoder()]
        for actual, center in zip(true_tensor, midpoint_tensor):
            error = c.q5_sub(actual, center)
            if c.q5_sign(error) < 0:
                error = c.q5_sub(c.q5(0), error)
            assert c.q5_sign(c.q5_sub(error, coefficient_error)) <= 0


def test_retained_receipt_and_all_port_covariance():
    result = b.run()
    assert result == json.loads(b.OUTPUT.read_text())
    assert result['ports_checked'] == 12
    assert result['contract']['initial_bit_disclosure_required'] is False
    assert result['contract']['physical_clock_or_field_identified'] is False
