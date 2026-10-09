"""Independent finite-chain controls for the supplied weighted-count domain.

The stationary oracle enumerates directed spanning trees. It neither solves
the stationary linear system nor iterates the producer's transition kernel.
"""

from fractions import Fraction as F
import copy
import itertools
import math

import numpy as np
import pytest

import exact_count_chain as chain
import collar_matrix_realization_probe as probe


def tree_stationary(kernel):
    """Markov-chain tree theorem, with every nonroot edge directed to root."""
    n = len(kernel)
    root_weights = []
    for root in range(n):
        sources = [i for i in range(n) if i != root]
        choices = [
            [j for j in range(n) if j != i and kernel[i][j] > 0]
            for i in sources
        ]
        total = F(0)
        for targets in itertools.product(*choices):
            parent = dict(zip(sources, targets))
            reaches_root = True
            for start in sources:
                seen = set()
                at = start
                while at != root and at not in seen:
                    seen.add(at)
                    at = parent[at]
                if at != root:
                    reaches_root = False
                    break
            if reaches_root:
                weight = F(1)
                for source, target in parent.items():
                    weight *= kernel[source][target]
                total += weight
        root_weights.append(total)
    normalizer = sum(root_weights)
    assert normalizer > 0, "oracle requires a unique recurrent class"
    return [weight / normalizer for weight in root_weights]


def assert_stationary_control(kernel, expected=None):
    result = chain.stationary_distribution(kernel)
    assert list(result) == tree_stationary(kernel)
    if expected is not None:
        assert list(result) == expected
    assert all(isinstance(value, F) for value in result)
    assert sum(result) == 1
    assert all(
        sum(result[i] * kernel[i][j] for i in range(len(kernel))) == result[j]
        for j in range(len(kernel))
    )


@pytest.mark.parametrize("a,b", list(itertools.product(range(1, 5), repeat=2)))
def test_exhaustive_two_state_quarter_grid_against_tree_oracle(a, b):
    kernel = [[1 - F(a, 4), F(a, 4)], [F(b, 4), 1 - F(b, 4)]]
    assert_stationary_control(kernel, [F(b, a + b), F(a, a + b)])
    assert chain.detailed_balance_defect(kernel, chain.stationary_distribution(kernel)) == 0


@pytest.mark.parametrize(
    "weights",
    [
        [[0, 1, 0], [1, 0, 1], [1, 0, 0]],
        [[9, 2, 0, 1], [0, 3, 5, 0], [1, 0, 1, 4], [2, 3, 0, 7]],
        [[2**42, 1, 0, 0, 0], [0, 3, 2, 0, 0], [1, 0, 2, 4, 0],
         [0, 1, 0, 0, 5], [2, 0, 3, 0, 7]],
    ],
)
def test_three_to_five_state_weighted_chains_against_tree_oracle(weights):
    # This normalization is independently derived from the original scalars.
    kernel = [[F(value, sum(row)) for value in row] for row in weights]
    assert chain.kernel_from_weights(weights) == kernel
    assert_stationary_control(kernel)


def test_slow_chain_has_resolved_stationary_population_despite_tiny_residual():
    denominator = 2**42
    weights = [[denominator - 1, 1], [2, denominator - 2]]
    kernel = chain.kernel_from_weights(weights)
    uniform_residual = F(1, 2 * denominator)
    assert uniform_residual < F(1, 10**12)
    assert_stationary_control(kernel, [F(2, 3), F(1, 3)])


def test_periodic_chain_stationary_law_exists_without_power_convergence():
    kernel = [[F(0), F(1, 2), F(1, 2)], [F(1), F(0), F(0)],
              [F(1), F(0), F(0)]]
    assert chain.irreducible_period(kernel) == 2
    assert_stationary_control(kernel, [F(1, 2), F(1, 4), F(1, 4)])


def test_loopless_coprime_cycles_are_aperiodic():
    # Return paths 0->1->0 and 0->1->2->0 have lengths 2 and 3.
    kernel = [[F(0), F(1), F(0)], [F(1, 2), F(0), F(1, 2)],
              [F(1), F(0), F(0)]]
    assert all(kernel[i][i] == 0 for i in range(3))
    assert chain.irreducible_period(kernel) == 1
    assert_stationary_control(kernel, [F(2, 5), F(2, 5), F(1, 5)])


def test_exhaustive_three_state_support_periods_by_simple_cycle_enumeration():
    """All 144 strongly connected labelled three-state support graphs.

    Connectivity uses direct/simple two-edge paths; the period oracle lists
    every simple directed cycle, independently of the producer's BFS levels.
    """
    period_counts = {1: 0, 2: 0, 3: 0}
    for bits in itertools.product((0, 1), repeat=9):
        support = [list(bits[3 * i:3 * (i + 1)]) for i in range(3)]
        connected = all(
            i == j or support[i][j] or any(
                support[i][k] and support[k][j] for k in range(3)
            )
            for i in range(3) for j in range(3)
        )
        if not connected:
            continue
        period = 0
        for length in range(1, 4):
            for cycle in itertools.permutations(range(3), length):
                if all(support[cycle[i]][cycle[(i + 1) % length]] for i in range(length)):
                    period = math.gcd(period, length)
        assert period > 0
        kernel = [[F(value, sum(row)) for value in row] for row in support]
        assert chain.irreducible_period(kernel) == period, support
        assert_stationary_control(kernel)
        period_counts[period] += 1
    assert period_counts == {1: 139, 2: 3, 3: 2}


def test_unique_absorbing_class_and_nonunique_stationarity_are_distinct():
    kernel = [[F(1, 2), F(1, 2), F(0)], [F(0), F(1, 2), F(1, 2)],
              [F(0), F(0), F(1)]]
    assert_stationary_control(kernel, [F(0), F(0), F(1)])
    with pytest.raises(ValueError):
        chain.irreducible_period(kernel)
    with pytest.raises(ValueError):
        chain.stationary_distribution([[F(1), F(0)], [F(0), F(1)]])


def test_exact_reversibility_and_tiny_irreversible_circulation():
    epsilon = F(1, 2**48)
    symmetric = [[F(1, 2), F(1, 4), F(1, 4)],
                 [F(1, 4), F(1, 2), F(1, 4)],
                 [F(1, 4), F(1, 4), F(1, 2)]]
    circulating = [row[:] for row in symmetric]
    for i in range(3):
        circulating[i][(i + 1) % 3] += epsilon
        circulating[i][(i - 1) % 3] -= epsilon
    uniform = [F(1, 3)] * 3
    assert_stationary_control(symmetric, uniform)
    assert_stationary_control(circulating, uniform)
    assert chain.detailed_balance_defect(symmetric, uniform) == 0
    assert chain.detailed_balance_defect(circulating, uniform) == 2 * epsilon / 3


@pytest.mark.parametrize("epsilon", [F(0), F(1, 2**48), F(1, 2**1100)])
def test_public_reversibility_preserves_exact_flux_decision(epsilon):
    import mpmath

    # Uniform stationarity follows directly from the unit column sums. Each
    # directed cycle edge has excess stationary flux 2*epsilon/3.
    kernel = [[F(1, 2), F(1, 4), F(1, 4)],
              [F(1, 4), F(1, 2), F(1, 4)],
              [F(1, 4), F(1, 4), F(1, 2)]]
    for i in range(3):
        kernel[i][(i + 1) % 3] += epsilon
        kernel[i][(i - 1) % 3] -= epsilon
    report = probe.audit_irreducible_chain(kernel, mpmath.mp.clone())
    expected_defect = 2 * epsilon / 3
    assert report["stationary_distribution_exact"] == ["1/3"] * 3
    assert F(report["detailed_balance_defect_exact"]) == expected_defect
    assert report["reversible"] is (epsilon == 0)
    assert report["detailed_balance_max_err"] == float(expected_defect)
    # Both an exact zero and a positive defect below binary64 have a zero
    # display, while the public scientific classifications remain different.
    if epsilon == F(1, 2**1100):
        assert report["detailed_balance_max_err"] == 0.0
        assert report["reversible"] is False


def test_exact_lumpability_and_arbitrarily_small_escape_defect():
    epsilon = F(1, 2**48)
    kernel = [[F(1, 2), F(1, 4), F(1, 4)],
              [F(1, 4), F(1, 2), F(1, 4)],
              [F(1, 4), F(1, 4), F(1, 2)]]
    blocks = [[0, 1], [2]]
    assert chain.lumpability_defect(kernel, blocks) == 0
    perturbed = [row[:] for row in kernel]
    perturbed[1][1] -= epsilon
    perturbed[1][2] += epsilon
    assert chain.lumpability_defect(perturbed, blocks) == epsilon


@pytest.mark.parametrize("defect", [F(0), F(1, 2**48), F(1, 2**1100)])
def test_lumpability_report_keeps_exact_and_tolerance_decisions_distinct(monkeypatch, defect):
    """Challenge the report boundary, not the physics of the pinned source.

    The preceding original-kernel control establishes the helper's exact
    zero/nonzero distinction. Inject those outcomes here to test the receipt
    consumer independently of the retained table's large lumpability defect.
    """
    original = chain.lumpability_defect

    def controlled_defect(kernel, blocks):
        if len(blocks) == 8:
            return defect
        return original(kernel, blocks)

    monkeypatch.setattr(chain, "lumpability_defect", controlled_defect)
    receipt = probe.build_probe()
    report = receipt["raw_coarsening_audit"]["selected_raw_equilibrium_probe"]
    assert F(report["fine_chain_strong_lumpability_defect_exact"]) == defect
    assert report["fine_chain_strongly_lumpable"] is (defect == 0)
    assert report["fine_chain_strongly_lumpable_at_tolerance"] is True
    assert report["fine_chain_strong_lumpability_max_err"] == float(defect)
    if defect == F(1, 2**1100):
        assert report["fine_chain_strong_lumpability_max_err"] == 0.0
        assert report["fine_chain_strongly_lumpable"] is False


def test_weight_scaling_preserves_kernel_and_original_scalar_precision():
    weights = [[2**53 + 1, 1.0], [F(1, 3), F(2, 3)]]
    expected = [[F(2**53 + 1, 2**53 + 2), F(1, 2**53 + 2)],
                [F(1, 3), F(2, 3)]]
    assert chain.kernel_from_weights(weights) == expected
    exact = [[F(value) for value in row] for row in weights]
    for scale in [F(1, 2**1100), F(2**1100)]:
        assert chain.kernel_from_weights([[scale * value for value in row] for row in exact]) == expected


def test_binary64_weights_preserve_smallest_population_and_avoid_sum_overflow():
    tiny = float.fromhex("0x0.0000000000001p-1022")
    kernel = chain.kernel_from_weights([[1.0, tiny], [1e308, 1e308]])
    tiny_fraction = F(1, 2**1074)
    assert kernel[0] == [1 / (1 + tiny_fraction), tiny_fraction / (1 + tiny_fraction)]
    assert kernel[1] == [F(1, 2), F(1, 2)]
    assert kernel[0][1] > 0


def test_empty_count_row_requires_an_explicit_absorbing_convention():
    weights = [[0, 0], [2, 1]]
    with pytest.raises(ValueError):
        chain.kernel_from_weights(weights)
    assert chain.kernel_from_weights(weights, empty_rows="absorbing") == [
        [F(1), F(0)], [F(2, 3), F(1, 3)]
    ]


@pytest.mark.parametrize("value", [True, np.bool_(False), -1, float("nan"), float("inf")])
def test_invalid_original_weight_scalars_are_not_silently_coerced(value):
    with pytest.raises((TypeError, ValueError)):
        chain.kernel_from_weights([[1, value], [1, 1]])


def test_masked_count_evidence_is_not_coerced_to_observed_values():
    weights = np.ma.array([[1, 2], [3, 4]], mask=[[False, True], [False, False]])
    with pytest.raises((TypeError, ValueError)):
        chain.kernel_from_weights(weights)


@pytest.fixture(scope="module")
def pinned_source():
    return probe.load_inputs()[0]


def test_pinned_weighted_source_binding_and_symmetric_degree_stationarity(pinned_source):
    supplied = [[F(value) for value in row] for row in pinned_source["counts"]]
    counts, raw, reversible = probe.source_kernels(pinned_source)
    assert counts == supplied
    n = len(supplied)
    assert raw == [[value / sum(row) for value in row] for row in supplied]
    symmetric = [[supplied[i][j] + supplied[j][i] for j in range(n)] for i in range(n)]
    degree = [sum(row) for row in symmetric]
    expected_pi = [value / sum(degree) for value in degree]
    assert reversible == [[value / sum(row) for value in row] for row in symmetric]
    assert probe.stationary_of(reversible) == expected_pi
    assert chain.detailed_balance_defect(reversible, expected_pi) == 0


@pytest.mark.parametrize("name", ["raw", "rev"])
def test_binding_rejects_smallest_added_edge_at_an_exact_source_zero(pinned_source, name):
    changed = copy.deepcopy(pinned_source)
    i, j = next((i, j) for i, row in enumerate(changed[name])
                for j, value in enumerate(row) if value == 0)
    changed[name][i][j] = float.fromhex("0x0.0000000000001p-1022")
    with pytest.raises(probe.ThermoError, match="recorded weighted-count export"):
        probe.source_kernels(changed)


def test_binding_rejects_original_weight_change_even_if_roundoff_sized(pinned_source):
    changed = copy.deepcopy(pinned_source)
    # A positive count changes by one ulp; the stored transitions stay fixed.
    i, j = next((i, j) for i, row in enumerate(changed["counts"])
                for j, value in enumerate(row) if value > 0 and sum(x > 0 for x in row) > 1)
    changed["counts"][i][j] = float(np.nextafter(changed["counts"][i][j], np.inf))
    with pytest.raises(probe.ThermoError, match="recorded weighted-count export"):
        probe.source_kernels(changed)


@pytest.mark.parametrize("large_diagonal, lost_export", [(1.0, "rev"), (1e200, "raw")])
def test_binding_refuses_export_that_lost_an_original_positive_edge(large_diagonal, lost_export):
    tiny = float.fromhex("0x0.0000000000001p-1022")
    counts = np.array([[large_diagonal, tiny], [0.0, 1.0]])
    symmetric = (counts + counts.T) / 2
    with np.errstate(under="ignore"):
        raw = counts / counts.sum(axis=1)[:, None]
        reversible = symmetric / symmetric.sum(axis=1)[:, None]
    data = {"counts": counts.tolist(), "raw": raw.tolist(), "rev": reversible.tolist()}
    assert data[lost_export][0][1] == 0
    if lost_export == "rev":
        assert raw[0, 1] == tiny  # Symmetrization's half operation loses this edge.
    # Original count normalization can resolve it, but the archived exported
    # dynamics has different support and must not inherit that classification.
    assert chain.kernel_from_weights(counts)[0][1] > 0
    with pytest.raises(probe.ThermoError, match="source export lost count-kernel support"):
        probe.source_kernels(data)


def test_tiny_nonzero_total_variation_is_not_rounded_before_reporting():
    tiny = F(1, 2**1100)
    kernel = [[F(1), F(0)], [1 - tiny, tiny]]
    assert probe.pairwise_row_tv_max(kernel) == tiny
    assert isinstance(probe.pairwise_row_tv_max(kernel), F)


def test_coarsening_preserves_distinct_typed_categorical_values():
    labels = [[["category", 1]], [["category", True]]]
    weights = [[2, 1], [1, 3]]
    assert probe.coordinate_partition_blocks(labels, ("category",)) == [[0], [1]]
    coarse_labels, coarse_weights, kernel = probe.coarsen_counts(weights, labels, ("category",))
    assert type(coarse_labels[0][0][1]) is int
    assert type(coarse_labels[1][0][1]) is bool
    assert coarse_weights == [[F(2), F(1)], [F(1), F(3)]]
    assert kernel == [[F(2, 3), F(1, 3)], [F(1, 4), F(3, 4)]]


def determinant_by_permutations(rows):
    """Independent Leibniz determinant over integers (at most 7x7 here)."""
    result = 0
    for columns in itertools.permutations(range(len(rows))):
        term = 1
        for i, j in enumerate(columns):
            term *= rows[i][j]
            if not term:
                break
        if term:
            inversions = sum(columns[i] > columns[j]
                             for i in range(len(rows)) for j in range(i + 1, len(rows)))
            result += -term if inversions % 2 else term
    return result


def weighted_tree_stationary(weights):
    """Tree cofactors of integer masses, including row-normalization factors.

    For a root r, each nonroot row contributes one 1/row_mass factor to a
    transition-tree weight. Thus stationary mass is proportional to
    row_mass[r] times the weighted-graph Laplacian's r-th cofactor.
    """
    scale = math.lcm(*(value.denominator for row in weights for value in row))
    integers = [[int(value * scale) for value in row] for row in weights]
    n = len(integers)
    masses = [sum(row) for row in integers]
    laplacian = [[(masses[i] if i == j else 0) - integers[i][j]
                  for j in range(n)] for i in range(n)]
    cofactors = []
    for root in range(n):
        minor = [[laplacian[i][j] for j in range(n) if j != root]
                 for i in range(n) if i != root]
        cofactors.append(masses[root] * determinant_by_permutations(minor))
    return [F(value, sum(cofactors)) for value in cofactors]


def test_pinned_eight_state_coarsening_has_independent_tree_stationarity(pinned_source):
    fine_weights = [[F(value) for value in row] for row in pinned_source["counts"]]
    buckets = [dict(label)["repair_load_bucket"] for label in pinned_source["labels"]]
    values = sorted(set(buckets))
    expected_weights = [[sum(fine_weights[i][j]
                             for i in range(len(buckets)) if buckets[i] == left
                             for j in range(len(buckets)) if buckets[j] == right)
                         for right in values] for left in values]
    labels, weights, kernel = probe.coarsen_counts(
        pinned_source["counts"], pinned_source["labels"], ("repair_load_bucket",))
    assert labels == [[["repair_load_bucket", value]] for value in values]
    assert weights == expected_weights
    expected_pi = weighted_tree_stationary(expected_weights)
    assert probe.stationary_of(kernel) == expected_pi
    assert all(value > 0 for value in expected_pi)


def test_kl_near_boundary_preserves_finite_divergence_without_rounding_u_to_minus_one():
    import mpmath

    ctx = mpmath.mp.clone()
    ctx.dps = 90
    tiny = F(1, 2**1000)
    actual = probe.kl([tiny, 1 - tiny], [F(1, 2), F(1, 2)], ctx)
    # Independent binary-entropy expression uses original p, not p/q - 1.
    t = ctx.mpf(1) / ctx.power(2, 1000)
    expected = ctx.log(2) + t * ctx.log(t) + (1 - t) * ctx.log1p(-t)
    assert ctx.isfinite(actual)
    assert ctx.almosteq(actual, expected)
