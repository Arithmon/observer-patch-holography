"""Hostile readout controls and independent exact refinement comparisons."""
from collections import Counter
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import refinement as r
from .source import binary_step, federation


PAIRS = ((0, 1), (1, 2), (0, 2))


def fraction(value):
    return F(value['numerator'], value['denominator'])


def laplacian(level):
    """Independent edge-incidence construction, not the production rows."""
    graph = federation(level)
    rows, columns, data = [], [], []
    for edge, (a, b) in enumerate(graph.seams):
        rows.extend((edge, edge)); columns.extend((a, b)); data.extend((1, -1))
    incidence = sparse.coo_matrix((data, (rows, columns)),
                                 shape=(len(graph.seams), graph.ports)).tocsr()
    return incidence.T @ incidence


@pytest.mark.parametrize('coarse,fine', PAIRS)
def test_static_conditional_expectation_and_state(coarse, fine):
    nc, nf = federation(coarse).ports, federation(fine).ports
    x = tuple(F((7*i) % 13, 11) for i in range(nc))
    y = tuple(F((17*i) % 19, 7) for i in range(nf))
    jx = r.lift(x, coarse, fine)
    assert r.restrict(jx, coarse, fine) == x
    ey = r.lift(r.restrict(y, coarse, fine), coarse, fine)
    assert r.lift(r.restrict(ey, coarse, fine), coarse, fine) == ey
    assert sum(jx)/nf == sum(x)/nc
    assert sum(ey) == sum(y)
    assert r.restrict(tuple(a*b for a, b in zip(jx, y)), coarse, fine) == tuple(
        a*b for a, b in zip(x, r.restrict(y, coarse, fine)))
    assert r.static_certificate(coarse, fine)['full_fixed_particle_configuration_state_claimed'] is False


def test_parent_lift_and_restriction_compose():
    assert r.parent_map(0, 2) == tuple(r.parent_map(0, 1)[v] for v in r.parent_map(1, 2))
    x = tuple(F(i, 3) for i in range(federation(0).ports))
    y = tuple(F(i % 17, 5) for i in range(federation(2).ports))
    assert r.lift(r.lift(x, 0, 1), 1, 2) == r.lift(x, 0, 2)
    assert r.restrict(r.restrict(y, 1, 2), 0, 1) == r.restrict(y, 0, 2)
    assert r.lift(x, 0, 0) == x
    assert r.restrict(x, 0, 0) == x


def test_adjoint_preservation_and_conditional_variance():
    x = tuple(complex(i % 3, i % 5) for i in range(federation(0).ports))
    assert r.lift(tuple(z.conjugate() for z in x), 0, 1) == tuple(
        z.conjugate() for z in r.lift(x, 0, 1))
    y = tuple(F(i % 7) for i in range(federation(1).ports))
    ry = r.restrict(y, 0, 1)
    squares = r.restrict(tuple(v*v for v in y), 0, 1)
    assert all(a >= b*b for a, b in zip(squares, ry))


def test_static_site_state_does_not_identify_full_fixed_two_configuration_laws():
    # Count unordered particle pairs in each descendant fibre, independently.
    for coarse, fine in PAIRS:
        sizes = Counter(r.parent_map(coarse, fine))
        nf = federation(fine).ports
        pairs = sum(b*(b-1)//2 for b in sizes.values())
        probability = F(pairs, nf*(nf-1)//2)
        control = r.static_certificate(coarse, fine)['fixed_two_static_projection_control']
        assert probability == fraction(control['fine_probability_two_particles_share_one_parent_port'])
        assert probability > fraction(control['coarse_binary_probability_two_particles_share_one_port'])


@pytest.mark.parametrize('coarse,fine', PAIRS)
def test_actual_parent_graph_not_conflated_with_cell_graph(coarse, fine):
    result = r.graph_parent_certificate(coarse, fine)
    assert result['cell_parent_graph_morphism']
    assert not result['port_parent_graph_morphism']
    witness = result['first_port_edge_failure']
    assert tuple(witness['fine_seam']) in federation(fine).seams
    assert tuple(sorted(witness['mapped_ports'])) not in federation(coarse).seams
    assert sum(result['port_edge_categories'].values()) == len(federation(fine).seams)


@pytest.mark.parametrize('coarse,fine', PAIRS)
def test_full_sparse_generator_defect_independently(coarse, fine):
    nc, nf = federation(coarse).ports, federation(fine).ports
    j = sparse.coo_matrix((np.ones(nf, dtype=np.int64),
                           (np.arange(nf), np.array(r.parent_map(coarse, fine)))),
                          shape=(nf, nc)).tocsr()
    for rate in (1, 4):
        actual = (laplacian(fine) @ j - rate * j @ laplacian(coarse)).tocsr()
        expected = r.defect_rows(coarse, fine, F(rate))
        for i, row in enumerate(expected):
            begin, end = actual.indptr[i:i+2]
            found = {int(k): F(int(v), 2)
                     for k, v in zip(actual.indices[begin:end], actual.data[begin:end]) if v}
            assert found == row
        assert r.defect_norm(coarse, fine, F(rate)) == (1 if rate == 1 else 19)
    assert r.defect_norm(coarse, coarse) == 0


@pytest.mark.parametrize('coarse,fine', ((0, 1), (1, 2)))
def test_nonlumpability_with_complete_actual_fixed_two_source_attempt(coarse, fine):
    witness = r.nonlumpability_witness(coarse, fine)
    first, second = map(tuple, witness['initial_fixed_two_states'])
    assert r.block_counts(first, coarse, fine) == r.block_counts(second, coarse, fine)
    target = witness['target_parent_port']
    pmap = r.parent_map(coarse, fine)
    totals = []
    for state in (first, second):
        total = 0
        for seam in federation(fine).seams:
            for coin in (0, 1):
                final, _ = binary_step(state, seam, coin)
                total += sum(pmap[v] == target for v in final)
        totals.append(total)
    difference = F(totals[0]-totals[1], 2*len(federation(fine).seams))
    assert difference == fraction(witness['one_source_attempt_block_mean_difference'])
    assert difference != 0
    assert difference * len(federation(fine).seams) == F(1, 2)


@pytest.mark.parametrize('coarse,fine', PAIRS)
def test_exact_common_step_comparison_and_continuous_duhamel_diagnostic(coarse, fine):
    result = r.finite_horizon_certificate(coarse, fine)
    assert fraction(result['finite_horizons'][0]['max_indicator_error']) == F(1, 12)
    for row in result['finite_horizons']:
        assert fraction(row['max_indicator_error']) <= fraction(row['telescope_indicator_bound'])
    target = result['target_block']
    x = np.zeros(federation(coarse).ports); x[target] = 1
    jx = np.array(r.lift(x, coarse, fine))
    for time in (.1, .4):
        fine_value = expm_multiply(-time*laplacian(fine)/2, jx)
        coarse_value = expm_multiply(-time*laplacian(coarse)/2, x)
        error = np.max(np.abs(fine_value-np.array(r.lift(coarse_value, coarse, fine))))
        assert 0 < error <= time/2 + 1e-13


def test_every_endpoint_pattern_current_law_including_zero_records():
    graph = federation(1)
    pmap = r.parent_map(0, 1)
    seams = (graph.seams[0], graph.inter_seams[0], graph.seams[-1])
    for seam in seams:
        for bits in product((0, 1), repeat=2):
            occupied = tuple(v for v, bit in zip(seam, bits) if bit)
            initial = r.block_counts(occupied, 0, 1)
            for coin in (0, 1):
                final, current = binary_step(occupied, seam, coin)
                record = (1, *seam, coin, current)
                assert r.continue_counts(initial, (record,), 0, 1) == r.block_counts(final, 0, 1)
                assert sum(r.continue_counts(initial, (record,), 0, 1)) == len(occupied)
                if pmap[seam[0]] == pmap[seam[1]]:
                    assert r.continue_counts(initial, (record,), 0, 1) == initial


def test_all_history_prefixes_checkpoint_and_omission_control():
    receipt = r.current_certificate(0, 1)
    assert receipt['native_local_cases_checked'] == 8*len(federation(1).seams)
    assert receipt['complete_coin_histories'] == 64
    assert receipt['prefix_block_readout_equalities'] == 384
    assert receipt['checkpoint_continuations_equal'] == 64
    assert receipt['zero_current_events_retained'] > 0
    assert receipt['nonzero_current_events'] > 0
    assert receipt['omitted_boundary_record_control_detects_mismatch']
    assert receipt['markov_closure_from_boundary_records_claimed'] is False


def test_checkpoint_retains_event_cursor_and_rejects_replay():
    state = (0, 1)
    seam = (1, 38)
    moved, current = binary_step(state, seam, 0)
    first = (1, *seam, 0, current)
    initial = r.block_counts(state, 0, 1)
    checkpoint = r.continue_counts(initial, (first,), 0, 1)
    assert checkpoint == r.block_counts(moved, 0, 1)
    with pytest.raises(ValueError):
        r.continue_counts(checkpoint, (first,), 0, 1, after_event=1)
    assert r.continue_counts(checkpoint, (), 0, 1, after_event=1) == checkpoint


@pytest.mark.parametrize('record', ((1,1,38,1,1), (1,1,38,0,-1),
                                    (0,1,38,0,1), (1,38,1,0,1),
                                    (1,1,2,0,1), (True,1,38,0,1)))
def test_invalid_current_records_rejected(record):
    with pytest.raises(ValueError):
        r.continue_counts(r.block_counts((0,1), 0, 1), (record,), 0, 1)


@pytest.mark.parametrize('coarse,fine', ((-1,1), (1,0), (0,3), (False,1), (0,True)))
def test_invalid_level_arguments_rejected(coarse, fine):
    with pytest.raises(ValueError):
        r.parent_map(coarse, fine)


def test_invalid_inputs_and_capacity_rejected():
    with pytest.raises(ValueError):
        r.parent(-1, 0, 1)
    with pytest.raises(ValueError):
        r.lift((1,), 0, 1)
    with pytest.raises(ValueError):
        r.restrict((1,), 0, 1)
    with pytest.raises(ValueError):
        r.defect_rows(0, 1, F(-1))
    with pytest.raises(ValueError):
        r.finite_horizon_certificate(0, 1, (2, 1))
    with pytest.raises(ValueError):
        r.block_counts((1, 1), 0, 1)
    with pytest.raises(ValueError):
        r.continue_counts((0,)*240, ((1,1,38,0,1),), 0, 1)
    with pytest.raises(ValueError):
        r.continue_counts((0,)*240, (), 0, 1, after_event=True)
