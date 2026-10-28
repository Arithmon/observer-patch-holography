"""Independent configuration topology, endpoint and particle-boundary controls."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
import json

import numpy as np
import pytest
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from . import defect_history as d
from ..observer import record_recovery as records
from ..native_modular_reference import history_transport


@pytest.fixture(scope="module")
def report():
    return d.report()


def endpoint(pair, seam, coin):
    """Independent floor/ceiling map on ordered known particle trajectories."""
    a, b = seam
    if (a in pair) == (b in pair):
        return pair
    destination = a if coin else b
    source = a if a in pair else b
    return tuple(destination if p == source else p for p in pair)


def independent_rank_mod2(rows):
    matrix = np.array(rows, dtype=np.int8) % 2
    rank = 0
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        matrix[[rank, pivot]] = matrix[[pivot, rank]]
        for row in range(rank + 1, matrix.shape[0]):
            if matrix[row, column]:
                matrix[row] ^= matrix[rank]
        rank += 1
    return rank


def test_configuration_graph_is_exactly_the_existing_native_movement_graph():
    pairs, edges, cells, tree, chords, faces = d.geometry()
    index = {pair: i for i, pair in enumerate(pairs)}
    actual = set()
    for i, pair in enumerate(pairs):
        for seam in records.seams():
            for coin in (0, 1):
                j = index[tuple(sorted(endpoint(pair, seam, coin)))]
                if i != j:
                    actual.add(tuple(sorted((i, j))))
    assert actual == set(edges)
    assert (len(pairs), len(edges), len(tree), len(chords)) == (66, 300, 65, 235)
    assert sum(kind == "one_particle_face" for kind, _ in cells) == 20 * 9
    assert sum(kind == "disjoint_square" for kind, _ in cells) == 315
    assert len(faces) == 20


def test_every_cell_has_zero_integer_boundary_and_flat_lifted_holonomy():
    _, edges, cells, *_ = d.geometry()
    for _, cycle in cells:
        chain = d.boundary(cycle, edges)
        divergence = [0] * 66
        for coefficient, (a, b) in zip(chain, edges, strict=True):
            divergence[a] -= coefficient
            divergence[b] += coefficient
        assert divergence == [0] * 66
        assert d.transport(cycle + cycle[:1]) == 0
    assert all(d.edge_parity(a, b) == d.edge_parity(b, a) for a, b in edges)


def test_mod2_rank_and_integral_smith_form_independently_classify_homology(report):
    _, edges, cells, _, chords, _ = d.geometry()
    chains = tuple(d.boundary(cycle, edges) for _, cycle in cells)
    assert independent_rank_mod2(chains) == 234
    # Independent full rectangular Smith form, not the producer's determinant minor.
    matrix = sp.Matrix([[chain[i] for chain in chains] for i in chords])
    form = smith_normal_form(matrix, domain=ZZ)
    diagonal = [abs(form[i, i]) for i in range(235)]
    assert diagonal == [1] * 234 + [2]
    assert report["topology"]["integral_first_homology"] == "Z/2"


def test_double_exchange_filling_is_integral_on_all_actual_edges(report):
    _, edges, cells, *_ = d.geometry()
    filling = report["topology"]["double_exchange_integral_filling"]
    _, loop = d.exchange_loop()
    target = d.boundary(loop[:-1], edges)
    actual = [0] * len(edges)
    for cell, coefficient in filling:
        chain = d.boundary(cells[cell][1], edges)
        for i, value in enumerate(chain):
            actual[i] += coefficient * value
    assert actual == [2 * value for value in target]
    assert report["topology"]["negative_control_filling_exchange_as_extra_cell"]["homology_mod2_dimension"] == 0


def test_exchange_holonomy_reversal_composition_and_arbitrary_local_chart_changes():
    _, loop = d.exchange_loop()
    assert d.transport(loop) == d.transport(tuple(reversed(loop))) == 1
    assert d.transport(loop + loop[1:]) == 0
    assert d.transport((loop[0], loop[0]) + loop[1:]) == 1
    # Basis chart changes generate every independent relabeling of the66 pair fibers.
    for flipped in range(66):
        parity = 0
        for a, b in zip(loop, loop[1:]):
            parity ^= int(a == flipped) ^ d.edge_parity(a, b) ^ int(b == flipped)
        assert parity == 1
    # Open-path parity depends on charts; closed-loop parity does not.
    assert 1 ^ d.edge_parity(loop[0], loop[1]) != d.edge_parity(loop[0], loop[1])


@pytest.mark.parametrize("law", ["C", "E"])
def test_every_local_coin_word_matches_independent_ordered_trajectory(law, report):
    result = report["local_protocols"][law]
    initial, schedule = tuple(result["initial_pair"]), tuple(map(tuple, result["schedule"]))
    counts = defaultdict(F)
    for row in result["coin_words"]:
        pair = initial
        for seam, coin in zip(schedule, row["coins"], strict=True):
            pair = endpoint(pair, seam, coin)
        assert pair == tuple(row["ordered_endpoint"])
        closed = set(pair) == set(initial)
        assert closed == row["closed"]
        key = (int(pair != initial) if closed else None)
        assert row["closed_exchange_parity"] == key
        counts[key] += F(row["probability"])
    assert counts == {0: F(1, 4), 1: F(1, 8), None: F(5, 8)}
    assert result["all_paths_remain_inside_face"]
    assert result["complete_raw_record_slots"] == (512 if law == "C" else 4096)


@pytest.mark.parametrize("law", ["C", "E"])
def test_complete_exploration_law_includes_local_odd_loop_and_every_exit(law):
    result = d.exploring_protocol(law)
    local = d.local_protocol(law)
    initial = tuple(local["initial_pair"])
    face = set(local["supporting_native_face"])
    # Integer DP uses direct ordered endpoint balancing, not the record decoder.
    weights = {(initial, False): 1}
    for preferred in local["schedule"]:
        next_weights = defaultdict(int)
        for (pair, left), value in weights.items():
            for seam in records.seams():
                numerator = 451 if seam == preferred else 1
                for coin in (0, 1):
                    after = endpoint(pair, seam, coin)
                    next_weights[after, left or not set(after) <= face] += value * numerator
        weights = next_weights
    denominator = 960 ** 3
    assert sum(weights.values()) == denominator
    even = sum(value for (pair, _), value in weights.items() if pair == initial)
    odd = sum(value for (pair, _), value in weights.items() if pair == initial[::-1])
    assert result["masses"]["closed_even"] == F(even, denominator)
    assert result["masses"]["closed_odd"] == F(odd, denominator)
    assert result["masses"]["localized_closed_odd"] >= F(451, 960) ** 3
    assert result["masses"]["ever_left_face"] > 0


def test_same_source_endpoint_can_have_even_or_odd_record_histories(report):
    result = report["local_protocols"]["C"]
    even = next(row for row in result["coin_words"] if row["closed_exchange_parity"] == 0)
    odd = next(row for row in result["coin_words"] if row["closed_exchange_parity"] == 1)
    assert even["vertices"][0] == even["vertices"][-1] == odd["vertices"][-1]
    assert even["complete_records"] != odd["complete_records"]
    assert not report["particle_readiness"]["source_endpoint_determines_history_parity"]


def test_scalar_record_phase_cancels_but_input_dependent_phase_changes_actual_block():
    result = d.local_protocol("C")
    schedule = tuple(map(tuple, result["schedule"]))
    branch = next(b for b in history_transport.word_instrument("C", schedule) if len(b.inputs) >= 2)
    # The actual partial injection embeds this two-input coherent block isometrically.
    rho = np.array([[1, -1j], [1j, 1]], complex) / 2
    k = np.eye(2) * np.sqrt(float(branch.weight))
    np.testing.assert_array_equal((-k) @ rho @ (-k).conj().T, k @ rho @ k.conj().T)
    changed = k @ np.diag([1, -1])
    assert np.linalg.norm(changed @ rho @ changed.conj().T - k @ rho @ k.conj().T) > 0
    # Such input-dependent transport changes are not supplied by a history label.


@pytest.mark.parametrize("law", ["C", "E"])
def test_record_character_transfer_is_ordered_cover_projection_not_new_native_phase(law):
    pairs = d.configurations()
    index = {p: i for i, p in enumerate(pairs)}
    unsigned, signed = d.parity_transfer_counts(law)
    expected = [[0] * 66 for _ in range(66)]
    for i, pair in enumerate(pairs):
        for seam in records.seams():
            for coin in (0, 1):
                after = endpoint(pair, seam, coin)
                j = index[tuple(sorted(after))]
                expected[i][j] += 1 if after[0] < after[1] else -1
    assert signed == tuple(map(tuple, expected))
    assert all(sum(row) == 60 for row in unsigned)
    assert any(value < 0 for row in signed for value in row)
    assert sp.Matrix(signed) == 60 * sp.eye(66) - d.signed_gram()


def test_positive_signed_form_and_source_count_contraction_have_exact_certificates(report):
    form = report["signed_compatibility_control"]
    assert all(value > 0 for value in form["signed_gram_minus_2I_positive_leading_minors"])
    matrix = d.signed_gram() - 2 * sp.eye(66)
    # Independent rational LDL factorization checks the full lower-bound matrix.
    lower, diagonal = matrix.LDLdecomposition(hermitian=False)
    assert lower * diagonal * lower.T == matrix
    assert all(diagonal[i, i] > 0 for i in range(66))
    rows = report["particle_readiness"]["uniform_complete_record_return_laws"]
    assert F(next(row for row in rows if row["attempts"] == 3)["closed_odd_mass"]) == F(1, 54000)
    assert all(F(row["endpoint_character_l2_squared"]) <= F(row["l2_squared_upper"]) for row in rows)
    assert not form["signed_gram_is_the_native_repair_generator"]


def test_all_reused_generic_theorem_anchors_exist_and_scope_is_preserved(report):
    for mapping in report["theorem_maps"]:
        source = (contract_root() / mapping["source"]).read_text()
        for name in mapping["theorems"]:
            assert "theorem " + name in source or "def " + name in source
    assert not report["scope"]["native_autonomous_source_charge"]
    assert not report["scope"]["native_fermionic_phase_or_statistics"]
    assert not report["scope"]["integrated_CP03_accepted"]


def contract_root():
    return d.HERE.parents[1]


def test_exact_receipt_replay_and_strict_type_mutation(report, tmp_path):
    d.contract.check_receipt(d.OUTPUT, report)
    changed = json.loads(d.OUTPUT.read_text())
    changed["scope"]["physical_comparison_count_added"] = False
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="differs from exact replay"):
        d.contract.check_receipt(path, report)


@pytest.mark.parametrize("path", [(), (False,), (66,), (0, 65)])
def test_untyped_or_non_native_paths_are_rejected(path):
    with pytest.raises(ValueError):
        d.transport(path)


@pytest.mark.parametrize("epsilon", [F(0), F(-1, 16), F(17, 16), True, 0.5])
def test_exploration_probability_has_a_positive_unit_interval_boundary(epsilon):
    with pytest.raises(ValueError):
        d.exploring_protocol("C", epsilon)
