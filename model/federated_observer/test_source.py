"""Independent finite checks of the elementary-to-federated source map."""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations, product
import json

import numpy as np
import pytest

from . import source as s


def test_pinned_federation_graphs_retain_the_elementary_carriers():
    s.verify_source_pins()
    assert len(s.SOURCE_PINS) == 49
    for level, cells in enumerate((20, 80, 320)):
        graph = s.federation(level)
        assert (graph.cells, graph.ports, len(graph.seams)) == (cells, 12*cells, 63*cells//2)
        assert len(graph.inter_seams) == 3*cells//2
        assert graph.max_degree == 6
        assert len({p for seam in graph.inter_seams for p in seam}) == 3*cells
        assert all(a//12 != b//12 for a, b in graph.inter_seams)
        assert sum(map(len, graph.neighbours)) == 2*len(graph.seams)
        with pytest.raises(FrozenInstanceError):
            graph.cells = 1


def test_actual_cell_ancestry_and_composition():
    first, second, both = s.cell_parents(0, 1), s.cell_parents(1, 2), s.cell_parents(0, 2)
    assert first == tuple(i//4 for i in range(80))
    assert second == tuple(i//4 for i in range(320))
    assert both == tuple(first[i] for i in second)
    assert s.cell_parents(1, 1) == tuple(range(80))


def test_all_endpoint_patterns_against_independent_fair_nearest_agreement():
    for x, y, coin in product((0, 1), repeat=3):
        occupied = tuple(p for p, v in ((2, x), (7, y), (20, 1)) if v)
        target, current = s.binary_step(occupied, (2, 7), coin)
        first = (x+y)//2+coin*((x+y) % 2)
        second = x+y-first
        assert int(2 in target) == first and int(7 in target) == second
        assert current == x-first
        assert 20 in target and len(target) == len(occupied)


def test_integer_local_kraus_certificate_preserves_coherence_and_references():
    matrices = tuple(np.asarray(m, dtype=int) for m in s.local_integer_kraus())
    np.testing.assert_array_equal(sum(m.T@m for m in matrices), 2*np.eye(4, dtype=int))
    np.testing.assert_array_equal(sum(m@m.T for m in matrices), 2*np.eye(4, dtype=int))
    # These two sectors share the same committed no-current outcome. A hidden
    # configuration measurement would destroy their reference coherence.
    vector = np.zeros(8, complex); vector[0] = vector[7] = 1/np.sqrt(2)
    density = np.outer(vector, vector.conj())
    output = sum(np.kron(m, np.eye(2))@density@np.kron(m.T, np.eye(2))/2 for m in matrices)
    np.testing.assert_allclose(output, density, atol=1e-15, rtol=0)
    assert abs(output[0, 7]) > .49


def test_fixed_total_exchange_generator_and_stationarity_on_all_active_inputs():
    graph, vertices, tape, basis, _ = s.tape_data()
    seams = tuple(sorted(set(tape)))
    size = len(basis)
    counts = np.zeros((size, size), dtype=int)
    generator = np.zeros((size, size), dtype=int)
    for a, b in seams:
        swap = np.zeros((size, size), dtype=int)
        for column, occupied in enumerate(basis):
            changed = tuple(sorted((b if p == a else a if p == b else p) for p in occupied))
            swap[basis.index(changed), column] = 1
            for coin in (0, 1):
                target, _ = s.binary_step(occupied, (a, b), coin)
                counts[basis.index(target), column] += 1
        generator += np.eye(size, dtype=int)-swap
    np.testing.assert_array_equal(counts, 2*len(seams)*np.eye(size, dtype=int)-generator)
    np.testing.assert_array_equal(counts, counts.T)
    np.testing.assert_array_equal(counts.sum(axis=1), np.full(size, 2*len(seams)))
    assert graph.ports == 240 and size == 6
    assert all(len(occupied) == 2 for occupied in basis)
    assert {sum(p//12 == tape[0][0]//12 for p in occupied) for occupied in basis} == {0, 1, 2}


def test_Bernoulli_stationarity_does_not_require_local_fixed_total_two():
    vertices = (0, 1, 2, 3)
    configurations = tuple(tuple(p for p, bit in zip(vertices, bits, strict=True) if bit)
                           for bits in product((0, 1), repeat=4))
    rho = F(1, 3)
    prior = {state: rho**len(state)*(1-rho)**(4-len(state)) for state in configurations}
    following = {state: F() for state in configurations}
    for state in configurations:
        for seam in ((0, 1), (1, 2), (2, 3)):
            for coin in (0, 1):
                target, _ = s.binary_step(state, seam, coin)
                following[target] += prior[state]/6
    assert following == prior
    assert sum(prior.values()) == 1


def test_complete_federated_readout_and_disconnected_control():
    result = s.tape_execution()
    assert result['complete_record_words'] == 81
    assert result['positive_words'] == 8
    assert result['second_carrier_final_occupation_law'] == {'0': '5/16', '1': '5/8', '2': '1/16'}
    assert result['mean_second_carrier_occupation'] == '3/4'
    assert result['disconnected_control_second_carrier_occupation'] == '0'
    assert sum(F(row['exact_probability']) for row in result['complete_record_law']) == 1
    assert result['full_global_fixed_total_dimension'] == 28680
    assert result['executed_exact_invariant_support_dimension'] == 6
    assert result['full_embedded_terminal_density_residual_diagnostic'] < 1e-13
    assert result['serialized_checkpoint_complete_distance_diagnostic'] == 0


def test_visible_current_instrument_preserves_global_coherences_inside_support():
    _, vertices, tape, basis, _ = s.tape_data()
    a, b = tape[0]
    other = [p for p in vertices if p not in (a, b)]
    first, second = (tuple(sorted((a, p))) for p in other)
    initial = np.zeros(12, complex)
    initial[2*basis.index(first)] = 1/np.sqrt(2)
    initial[2*basis.index(second)+1] = 1/np.sqrt(2)
    density = np.outer(initial, initial.conj())
    current_plus = s.tape_program().events[0].instruments[0].kraus[2]
    output = sum(np.kron(k, np.eye(2))@density@np.kron(k.conj().T, np.eye(2)) for k in current_plus)
    targets = [basis.index(tuple(sorted((b, p)))) for p in other]
    assert np.trace(output).real == pytest.approx(.5, abs=1e-15)
    assert abs(output[2*targets[0], 2*targets[1]+1]) == pytest.approx(.25, abs=1e-15)


def test_source_subsection_has_no_tuple_json_boundary():
    value = s.local_audit()
    assert json.loads(json.dumps(value)) == value


def test_changed_retained_source_hash_is_rejected(monkeypatch):
    monkeypatch.setattr(s, 'SOURCE_PINS', {'oph_exact/federation.py': '0'*64})
    with pytest.raises(ValueError, match='pinned federation source changed'):
        s.verify_source_pins()


@pytest.mark.parametrize('bad', [True, -1, 3, 0., '0', None])
def test_federation_levels_fail_closed_even_after_valid_cache_entries(bad):
    s.federation(0)
    with pytest.raises(ValueError):
        s.federation(bad)


@pytest.mark.parametrize('occupied,seam,coin', [
    ((True,), (0, 1), 0), ((0, 0), (0, 1), 0), ((2, 1), (0, 1), 0),
    ((-1,), (0, 1), 0), ([0], (0, 1), 0), ((0,), [0, 1], 0),
    ((0,), (0, True), 0), ((0,), (1, 0), 0), ((0,), (0, 0), 0),
    ((0,), (0, 1), True), ((0,), (0, 1), .0), ((0,), (0, 1), 2),
])
def test_binary_inputs_fail_closed(occupied, seam, coin):
    with pytest.raises(ValueError):
        s.binary_step(occupied, seam, coin)


@pytest.mark.parametrize('coarse,fine', [(True, 1), (0, 1.), (2, 1), (-1, 1), (0, 3)])
def test_parent_domain_is_declared(coarse, fine):
    with pytest.raises(ValueError):
        s.cell_parents(coarse, fine)
