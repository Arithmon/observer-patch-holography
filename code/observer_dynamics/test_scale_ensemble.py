"""Independent controls for the retained scale and angular estimator population.

These checks start from the pinned integer mass or original small vectors.
They do not claim to reconstruct the unavailable large terminal arrays.
"""
from fractions import Fraction
from collections import Counter
from itertools import combinations
import json
from math import comb
from statistics import mean

import numpy as np
import pytest

from .archive_adapter import SIM
from . import retained_measurements as retained


@pytest.fixture
def population():
    return json.loads((SIM / "data/codex_audit_20260925/scale_ensemble.json").read_text())


def reject(population):
    with pytest.raises(ValueError):
        retained.scale_ensemble(population)


def angular_first(population):
    return next(iter(population["levels"][0]["angular"].values()))["rows"][0]


def test_complete_frozen_population_and_independent_parent_nulls(population):
    result = retained.scale_ensemble(population)
    assert result["terminal_digest_crosspins"] == 48
    assert result["scale_rows"] == 51
    assert result["raw_terminal_arrays_rehashed"] is False
    checked = 0
    for level in population["levels"]:
        parent = json.loads((SIM / (
            "data/codex_audit_20260925/inputs/data/receipts/"
            f"L{level['level']}/receipt.json"
        )).read_text())
        ports = parent["ports"]
        raised = parent["expected"]["component_total"][0] % ports
        for scale in level["scales"]:
            if not scale["name"].startswith("ancestor_"):
                continue
            group_ports = 12 * int(scale["name"].split("_")[1])
            # Scalar hypergeometric variance of an exhaustive group's mean.
            # The population-wide mean is fixed by the original integer mass.
            expected = Fraction(
                raised * (ports - raised) * (ports - group_ports),
                ports * ports * group_ports * (ports - 1),
            )
            assert scale["equilibrium_variance"] == pytest.approx(float(expected), rel=1e-13, abs=0)
            checked += 1
    assert checked == 27


def test_coherent_null_rescaling_cannot_reverse_enhancement(population):
    scale = population["levels"][0]["scales"][0]
    assert mean(row["variance_over_exchangeable_equilibrium"] for row in scale["per_schedule"]) > 1
    scale["equilibrium_variance"] *= 2
    for row in scale["per_schedule"]:
        row["variance_over_exchangeable_equilibrium"] /= 2
    assert mean(row["variance_over_exchangeable_equilibrium"] for row in scale["per_schedule"]) < 1
    reject(population)


def test_differing_responses_require_positive_schedule_noise(population):
    scale = population["levels"][0]["scales"][0]
    responses = [Fraction.from_float(row["response"]) for row in scale["per_schedule"]]
    center = sum(responses) / len(responses)
    initial = Fraction.from_float(scale["per_schedule"][0]["initial_variance"])
    # Projection onto the fixed initial vector cannot increase variance.
    minimum = initial * sum((value - center) ** 2 for value in responses) / (len(responses) - 1)
    assert float(minimum) == pytest.approx(7.028510514184177e-6, rel=1e-13)
    scale["conditional_schedule_variance_unbiased"] = 0.0
    scale["cross_schedule_common_power_unbiased"] = scale["total_terminal_variance"]
    scale["schedule_noise_fraction"] = 0.0
    reject(population)


@pytest.mark.parametrize("field", ["initial_variance_over_iid_null", "groups", "cells_min", "cells_max"])
def test_ancestral_geometry_and_iid_normalization_are_bound_to_declared_groups(population, field):
    population["levels"][0]["scales"][0][field] *= 2
    reject(population)


def test_reported_coherence_must_agree_with_response_and_power(population):
    population["levels"][0]["scales"][0]["per_schedule"][0]["coherence_squared"] *= 0.5
    reject(population)


def test_schedules_share_one_initial_field(population):
    row = population["levels"][0]["scales"][0]["per_schedule"][0]
    row["initial_variance"] *= 1.1
    row["power_ratio"] = row["terminal_variance"] / row["initial_variance"]
    row["coherence_squared"] = row["response"] ** 2 / row["power_ratio"]
    row["residual_power_ratio"] = row["power_ratio"] - row["response"] ** 2
    reject(population)


@pytest.mark.parametrize("field", ["ports", "carriers", "balanced_fraction", "initial_load_seed"])
def test_original_parent_population_metadata_cannot_change(population, field):
    population["levels"][0][field] *= 2
    reject(population)


def test_geometry_digest_matches_parent_receipt(population):
    population["levels"][0]["geometry_file_sha256"] = "0" * 64
    reject(population)


@pytest.mark.parametrize("change", ["remove", "duplicate", "empty"])
def test_specification_requires_the_complete_scale_population(population, change):
    scales = population["levels"][0]["scales"]
    if change == "remove":
        scales.pop(0)  # Keep cube pairs intact: the census must catch an absent ancestor.
    elif change == "duplicate":
        scales[-1] = scales[0]
    else:
        scales.clear()
    reject(population)


def test_changed_embedded_specification_cannot_excuse_missing_scale(population):
    population["spec"]["primary_ancestral_group_cells"].remove(65536)
    for level in population["levels"]:
        level["scales"] = [scale for scale in level["scales"] if scale["name"] != "ancestor_65536"]
    reject(population)


@pytest.mark.parametrize("change", ["grid", "row", "duplicate_ell", "empty"])
def test_specification_requires_all_angular_grids_and_multipoles(population, change):
    angular = population["levels"][0]["angular"]
    first = next(iter(angular))
    if change == "grid":
        del angular[first]
    elif change == "row":
        angular[first]["rows"].pop()
    elif change == "duplicate_ell":
        angular[first]["rows"][1]["ell"] = 1
    else:
        angular.clear()
    reject(population)


@pytest.mark.parametrize("field", [
    "C_terminal_mean", "cross_response_mean", "cross_response_schedule_sd",
    "total_power_ratio_mean", "coherence_squared_mean", "orthogonal_residual_power_ratio_mean",
])
def test_angular_observables_agree_with_retained_schedule_population(population, field):
    angular_first(population)[field] *= 2
    reject(population)


def test_angular_differing_responses_cannot_have_zero_schedule_noise(population):
    row = angular_first(population)
    assert len(set(row["per_schedule_response"])) > 1
    row["conditional_schedule_noise_power"] = 0.0
    row["cross_schedule_common_power"] = row["C_terminal_mean"]
    reject(population)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), True])
def test_nonphysical_coherence_scalars_are_rejected(population, bad):
    population["levels"][0]["scales"][0]["per_schedule"][0]["coherence_squared"] = bad
    reject(population)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), True])
def test_nonphysical_angular_scalars_are_rejected(population, bad):
    angular_first(population)["coherence_squared_mean"] = bad
    reject(population)


def test_reordered_schedule_entries_cannot_break_digest_alignment(population):
    rows = population["levels"][0]["scales"][0]["per_schedule"]
    rows[0], rows[1] = rows[1], rows[0]
    reject(population)


def test_reordered_terminal_inputs_cannot_break_parent_alignment(population):
    rows = population["levels"][0]["inputs"]
    rows[0], rows[1] = rows[1], rows[0]
    reject(population)


def original_vector_statistics(coefficients):
    """Compute input observables from four-port, mean-zero exact vectors."""
    initial = (1, -1, 1, -1)
    residual_u = (1, 1, -1, -1)
    residual_v = (1, -1, -1, 1)
    fields = [tuple(Fraction(a * x + b * u + c * v)
                    for x, u, v in zip(initial, residual_u, residual_v))
              for a, b, c in coefficients]
    count = len(fields)
    average = tuple(sum(field[index] for field in fields) / count for index in range(4))
    initial_power = sum(Fraction(value * value, 4) for value in initial)
    responses = [sum(x * y for x, y in zip(initial, field)) / (4 * initial_power)
                 for field in fields]
    powers = [sum(value * value for value in field) / (4 * initial_power)
              for field in fields]
    # Start from the original coordinate-wise sample variance, independently
    # of the residual-polygon reduction used by the verifier.
    noise = sum((field[index] - average[index]) ** 2
                for field in fields for index in range(4)) / (4 * (count - 1) * initial_power)
    return [float(value) for value in responses], [float(value) for value in powers], float(noise)


@pytest.mark.parametrize("coefficients", [
    [(1, 0, 0), (-1, 0, 0)],  # Negative unbiased common power is valid.
    [(2, 1, 0), (-1, 2, 0), (0, 3, 0)],  # Maximum residual-sum length.
    [(1, 4, 0), (2, -1, 0), (-1, -1, 0)],  # Positive minimum sum length.
    [(1, 1, 0), (-1, 0, 1), (0, -1, -1)],  # Zero minimum sum length.
    [(1, 1, 0), (1, 1, 0)],  # Identical schedules: genuine zero noise.
])
def test_ensemble_bounds_accept_original_vector_populations_and_boundaries(coefficients):
    responses, powers, noise = original_vector_statistics(coefficients)
    retained.ensemble_bounds(responses, powers, noise)


def test_negative_unbiased_common_power_is_not_clamped():
    responses, powers, noise = original_vector_statistics([(1, 0, 0), (-1, 0, 0)])
    assert mean(powers) - noise == -1
    retained.ensemble_bounds(responses, powers, noise)


@pytest.mark.parametrize("coefficients,noise_change", [
    ([(2, 1, 0), (-1, 2, 0), (0, 3, 0)], -0.001),
    ([(1, 4, 0), (2, -1, 0), (-1, -1, 0)], 0.001),
    ([(1, 1, 0), (-1, 0, 1), (0, -1, -1)], 0.001),
])
def test_ensemble_bounds_reject_both_sides_outside_original_vector_boundaries(coefficients, noise_change):
    responses, powers, noise = original_vector_statistics(coefficients)
    with pytest.raises(ValueError):
        retained.ensemble_bounds(responses, powers, noise + noise_change)


@pytest.mark.parametrize("responses,powers,noise", [
    ([1.0, -1.0], [0.5, 1.0], 1.0),  # Individual Cauchy violation.
    ([1.0, -1.0], [1.0, 1.0], 0.0),  # Differing projections cannot be noiseless.
    ([0.0, 0.0], [1.0, 1.0], -0.1),
    ([0.0, 0.0], [1.0, float("nan")], 1.0),
    ([0.0, 0.0], [1.0, float("inf")], 1.0),
    ([0.0, 0.0], [1.0, True], 1.0),
    ([False, 0.0], [1.0, 1.0], 1.0),
    ([0.0, 0.0], [1.0, 1.0], True),
    ([0.0, 0.0], [1.0], 1.0),
    ([0.0], [1.0], 1.0),
])
def test_ensemble_bounds_reject_invalid_covariance_inputs(responses, powers, noise):
    with pytest.raises(ValueError):
        retained.ensemble_bounds(responses, powers, noise)


def rescaled_statistics(responses, powers, noise, gain):
    """Changing every terminal vector's units rescales its actual moments."""
    return ([gain * response for response in responses],
            [gain**2 * power for power in powers], gain**2 * noise)


@pytest.mark.parametrize("delta", [5e-9, 8e-9, 1e-8])
@pytest.mark.parametrize("gain", [1e-120, 1.0, 1e120])
def test_original_nearly_parallel_vectors_keep_their_residual_roundoff_budget(delta, gain):
    initial = np.array([1.0, 1.0, -1.0, -1.0])
    transverse = np.array([1.0, -1.0, 1.0, -1.0])
    fields = np.array([transverse] * 8 + [initial + delta * transverse] * 8)
    responses = np.mean(fields * initial, axis=1).tolist()
    powers = np.mean(fields**2, axis=1).tolist()
    # This calculation starts from coordinate-wise fluctuations. It neither
    # subtracts squared responses nor uses a residual-resultant identity.
    noise = float(np.mean(np.var(fields, axis=0, ddof=1)))
    exact_last = [Fraction(float(entry)) for entry in fields[-1]]
    exact_response = sum(Fraction(int(x)) * y for x, y in zip(initial, exact_last)) / 4
    exact_power = sum(y * y for y in exact_last) / 4
    assert exact_power - exact_response**2 > 0
    # Squaring and averaging the original binary64 data loses this positive
    # residual. A final tolerance on the resultant alone then falsely rejects
    # these vectors: sqrt near zero amplifies the permitted input roundoff.
    assert powers[-1] - responses[-1]**2 == 0.0
    retained.ensemble_bounds(*rescaled_statistics(responses, powers, noise, gain))


@pytest.mark.parametrize("gain", [1e-120, 1.0, 1e120])
def test_identical_original_vectors_allow_only_roundoff_negative_noise(gain):
    initial = np.array([1.0, 1.0, -1.0, -1.0])
    field = np.array([0.5, 0.1, -0.3, -0.9])
    field = field - field.mean()
    fields = np.tile(field, (16, 1))
    response = float(np.mean(initial * field))
    power = float(np.mean(field**2))
    # The retained estimator subtracts these two independently accumulated
    # powers. Exact original vectors give zero schedule variance because
    # every field equals its population mean coordinate by coordinate.
    exact_field = tuple(Fraction(float(value)) for value in field)
    exact_mean = tuple(sum([value] * 16) / 16 for value in exact_field)
    assert sum((value - average)**2 for value, average in zip(exact_field, exact_mean)) == 0
    summed_power = sum([power] * 16)
    power_of_sum = float(np.mean(np.sum(fields, axis=0)**2))
    noise = (summed_power - power_of_sum / 16) / 15
    assert noise == pytest.approx(-5.921189464667501e-17, rel=1e-14, abs=0)
    responses, powers = [response] * 16, [power] * 16
    retained.ensemble_bounds(*rescaled_statistics(responses, powers, noise, gain))
    with pytest.raises(ValueError):
        retained.ensemble_bounds(*rescaled_statistics(responses, powers, -1e-4 * power, gain))


@pytest.mark.parametrize("weighted", [False, True])
def test_geometry_factor_from_original_small_port_populations(weighted):
    cells = [4, 4, 4, 8]
    group_ports = [12 * count for count in cells]
    # Twelve ports outside these disjoint groups exercise the partial-coverage
    # case. Group sizes satisfy the archived cube estimator's minimum of four.
    port_groups = [group for group, size in enumerate(group_ports) for _ in range(size)] + [None] * 12
    ports = len(port_groups)
    weights = cells if weighted else [1] * len(cells)

    def centered_variance(values):
        center = sum(w * value for w, value in zip(weights, values)) / sum(weights)
        return sum(w * (value - center)**2 for w, value in zip(weights, values)) / sum(weights)

    coefficients = []
    for raised in (1, 2):
        populations = Counter()
        for occupied in combinations(range(ports), raised):
            counts = [0] * len(cells)
            for port in occupied:
                group = port_groups[port]
                if group is not None:
                    counts[group] += 1
            populations[tuple(counts)] += 1
        expected = sum(multiplicity * centered_variance(
            [Fraction(count, size) for count, size in zip(counts, group_ports)])
            for counts, multiplicity in populations.items()) / comb(ports, raised)
        fraction = Fraction(raised, ports)
        coefficients.append(expected / (fraction * (1 - fraction) * Fraction(ports, ports - 1)))
    assert coefficients[0] == coefficients[1] == (Fraction(1, 80) if weighted else Fraction(7, 512))

    # For independent uniform {0,...,5} port loads, centered innovations have
    # zero cross terms. Sum each original port's six quadratic contributions.
    iid_variance = Fraction(0)
    for group in port_groups:
        for value in range(6):
            means = [Fraction(0)] * len(cells)
            if group is not None:
                means[group] = (Fraction(value) - Fraction(5, 2)) / group_ports[group]
            iid_variance += centered_variance(means) / 6
    assert iid_variance == Fraction(35, 12) * coefficients[0]
    geometry = {
        "name": "cube_2_weighted" if weighted else "cube_2_unweighted",
        "groups": len(cells), "cells_min": min(cells), "cells_max": max(cells),
        "cells_arithmetic_mean": sum(cells) / len(cells),
        "cells_harmonic_mean": float(Fraction(len(cells)) / sum(Fraction(1, count) for count in cells)),
    }
    assert float(retained.geometry_factor(geometry, ports // 12)) == pytest.approx(
        float(coefficients[0]), rel=1e-14, abs=0)


def test_scale_consumer_preserves_negative_common_power_from_discrete_ports():
    initial_groups = ([2] * 3 + [3] * 9, [2] * 9 + [3] * 3) * 2
    terminal_groups = ([0] * 3 + [1] * 9, [0] * 9 + [1] * 3) * 2
    opposite_groups = tuple([1 - port for port in group] for group in terminal_groups)

    def centered_means(groups):
        values = [Fraction(sum(group), len(group)) for group in groups]
        center = sum(values) / len(values)
        return [value - center for value in values]

    def inner(left, right):
        return sum(x * y for x, y in zip(left, right)) / len(left)

    x = centered_means(initial_groups)
    fields = [centered_means(groups) for groups in (terminal_groups, opposite_groups)]
    assert sum(map(sum, terminal_groups)) == sum(map(sum, opposite_groups)) == 24
    initial = inner(x, x)
    powers = [inner(field, field) for field in fields]
    average = [(left + right) / 2 for left, right in zip(*fields)]
    noise = sum(inner([value - center for value, center in zip(field, average)],
                      [value - center for value, center in zip(field, average)])
                for field in fields)
    common = inner(*fields)
    assert initial == Fraction(1, 16) and noise == Fraction(1, 8) and common == -initial
    coefficient, factor = Fraction(12, 47), Fraction(1, 16)
    equilibrium = coefficient * factor
    rows = []
    for seed, (field, power) in enumerate(zip(fields, powers)):
        response = inner(x, field) / initial
        rows.append({
            "seed": seed, "initial_variance": float(initial), "terminal_variance": float(power),
            "response": float(response), "power_ratio": float(power / initial),
            "coherence_squared": float(inner(x, field)**2 / (initial * power)),
            "residual_power_ratio": 0.0,
            "variance_over_exchangeable_equilibrium": float(power / equilibrium),
            "legacy_unweighted_global_centered_response": float(response),
        })
    scale = {
        "name": "ancestor_1", "groups": 4, "cells_min": 1, "cells_max": 1,
        "cells_arithmetic_mean": 1.0, "cells_harmonic_mean": 1.0,
        "initial_variance_over_iid_null": float(initial / (Fraction(35, 12) * factor)),
        "equilibrium_variance": float(equilibrium), "per_schedule": rows,
        "response_mean": 0.0, "response_schedule_sd": 2**0.5, "response_schedule_se": 1.0,
        "total_terminal_variance": float(sum(powers) / 2),
        "conditional_schedule_variance_unbiased": float(noise),
        "cross_schedule_common_power_unbiased": float(common), "schedule_noise_fraction": 2.0,
    }
    retained.check_scale(scale, [0, 1], 4, coefficient)


def test_angular_consumer_preserves_negative_common_power_from_opposite_harmonics():
    # Original real ell=1 coefficient vectors ordered by m=-1,0,1.
    initial = (0, 1, 0)
    fields = (initial, tuple(-value for value in initial))

    def inner(left, right):
        return sum(Fraction(x * y, 3) for x, y in zip(left, right))

    power = inner(initial, initial)
    common = inner(*fields)
    average = tuple(Fraction(left + right, 2) for left, right in zip(*fields))
    noise = sum(inner(tuple(value - center for value, center in zip(field, average)),
                      tuple(value - center for value, center in zip(field, average)))
                for field in fields)
    assert power == Fraction(1, 3) and noise == Fraction(2, 3) and common == -power
    row = {
        "ell": 1, "C_initial": float(power),
        "per_schedule_response": [float(inner(initial, field) / power) for field in fields],
        "per_schedule_power_ratio": [float(inner(field, field) / power) for field in fields],
        "cross_response_mean": 0.0, "cross_response_schedule_sd": 2**0.5,
        "total_power_ratio_mean": 1.0, "C_terminal_mean": float(power),
        "coherence_squared_mean": 1.0, "orthogonal_residual_power_ratio_mean": 0.0,
        "conditional_schedule_noise_power": float(noise), "cross_schedule_common_power": float(common),
    }
    retained.check_angular(row, 2)


@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("cells", [
    [4, 5, 5, 8],  # Smallest reciprocal sum after reserving both extrema.
    [4, 4, 6, 8],  # Largest reciprocal sum with one interior count.
    [4, 4, 8, 8],  # Largest reciprocal sum entirely at the extrema.
    [4, 4, 4, 4],  # Coincident extrema.
])
def test_integer_cube_bounds_keep_attainable_extrema_and_partial_coverage(cells, weighted):
    count = len(cells)
    arithmetic = Fraction(sum(cells), count)
    harmonic = Fraction(count) / sum(Fraction(1, value) for value in cells)
    weights = cells if weighted else [1] * count
    total_weight = sum(weights)
    # Start from the actual original groups and their normalized weights.
    expected = sum((Fraction(weight, total_weight) - Fraction(weight, total_weight)**2)
                   / (12 * size) for weight, size in zip(weights, cells))
    geometry = {
        "name": "cube_2_weighted" if weighted else "cube_2_unweighted",
        "groups": count, "cells_min": min(cells), "cells_max": max(cells),
        "cells_arithmetic_mean": float(arithmetic), "cells_harmonic_mean": float(harmonic),
    }
    assert float(retained.geometry_factor(geometry, sum(cells) + 1)) == pytest.approx(
        float(expected), rel=1e-14, abs=0)


@pytest.mark.parametrize("arithmetic,harmonic", [
    (4.0, 4.0),  # The declared maximum eight must occur at least once.
    (8.0, 8.0),  # The declared minimum four must occur at least once.
    (5.1, 5.0),  # Four integer counts cannot sum to 20.4.
    (5.5, 4.9),  # Reciprocal sum above every possible four-count population.
    (5.5, 5.3),  # Reciprocal sum below every possible four-count population.
])
def test_integer_cube_bounds_reject_impossible_summary_populations(arithmetic, harmonic):
    geometry = {
        "name": "cube_2_unweighted", "groups": 4, "cells_min": 4, "cells_max": 8,
        "cells_arithmetic_mean": arithmetic, "cells_harmonic_mean": harmonic,
    }
    with pytest.raises(ValueError, match="integer cube"):
        retained.geometry_factor(geometry, 40)


def test_unused_carriers_do_not_expand_cube_mean_roundoff_budget():
    geometry = {
        "name": "cube_2_unweighted", "groups": 4, "cells_min": 4, "cells_max": 8,
        "cells_arithmetic_mean": 5.00001, "cells_harmonic_mean": float(Fraction(32, 7)),
    }
    # The arithmetic mean fixes a noninteger total 20.00004 even though a
    # much larger parent population lies outside the four retained groups.
    with pytest.raises(ValueError, match="integer cube population"):
        retained.geometry_factor(geometry, 10**9)


@pytest.mark.parametrize("direction", [-float("inf"), float("inf")])
def test_integer_cube_means_allow_last_bit_roundoff(direction):
    from math import nextafter

    cells = [4, 5, 5, 8]
    harmonic = Fraction(len(cells)) / sum(Fraction(1, value) for value in cells)
    geometry = {
        "name": "cube_2_unweighted", "groups": len(cells),
        "cells_min": min(cells), "cells_max": max(cells),
        "cells_arithmetic_mean": nextafter(sum(cells) / len(cells), direction),
        "cells_harmonic_mean": nextafter(float(harmonic), direction),
    }
    retained.geometry_factor(geometry, sum(cells) + 1)


def test_coherent_cube_renormalization_requires_possible_integer_geometry(population):
    level = population["levels"][0]
    parent = json.loads((SIM / (
        "data/codex_audit_20260925/inputs/data/receipts/L8/receipt.json"
    )).read_text())
    count, total, low, high = 632, parent["carriers"], 2073, 2074
    high_count = total - count * low
    low_count = count - high_count
    assert (low_count, high_count) == (48, 584)
    arithmetic = Fraction(total, count)
    harmonic = Fraction(count) / (Fraction(low_count, low) + Fraction(high_count, high))
    # The total integer count is uniquely determined even with the specified
    # arithmetic tolerance. The implied harmonic discrepancy exceeds it.
    assert retained.RTOL * total < Fraction(1, 2)
    assert (arithmetic - harmonic) / harmonic > 10 * retained.RTOL
    ports = parent["ports"]
    raised = parent["expected"]["component_total"][0] % ports
    coefficient = Fraction(raised * (ports - raised), ports * (ports - 1))
    for scale in level["scales"]:
        if not scale["name"].startswith("cube_12_"):
            continue
        scale.update(cells_min=low, cells_max=high,
                     cells_arithmetic_mean=float(arithmetic), cells_harmonic_mean=float(arithmetic))
        factor = Fraction(count - 1, 12 * count) / arithmetic
        scale["equilibrium_variance"] = float(coefficient * factor)
        scale["initial_variance_over_iid_null"] = (
            scale["per_schedule"][0]["initial_variance"] / float(Fraction(35, 12) * factor))
        for row in scale["per_schedule"]:
            row["variance_over_exchangeable_equilibrium"] = (
                row["terminal_variance"] / scale["equilibrium_variance"])
    with pytest.raises(ValueError, match="integer cube harmonic"):
        retained.scale_ensemble(population)


@pytest.mark.parametrize("attained", ["minimum", "maximum"])
def test_angular_extrema_must_fit_the_full_pixel_population(population, attained):
    level = population["levels"][0]
    grid = level["angular"]["z128_phi256"]
    pixels, carriers = 128 * 256, level["carriers"]
    assert carriers % pixels == 0
    if attained == "maximum":
        grid["counts_min"], grid["counts_max"] = 1, carriers
        assert grid["counts_max"] + (pixels - 1) * grid["counts_min"] > carriers
    else:
        grid["counts_min"], grid["counts_max"] = 1, carriers // pixels
        assert grid["counts_min"] + (pixels - 1) * grid["counts_max"] < carriers
    with pytest.raises(ValueError, match="angular pixel population"):
        retained.scale_ensemble(population)


@pytest.mark.parametrize("attained", ["minimum", "maximum", "both"])
def test_angular_extrema_accept_original_integer_pixel_counts(population, attained):
    level = population["levels"][0]
    grid = level["angular"]["z128_phi256"]
    pixels, carriers = 128 * 256, level["carriers"]
    average = carriers // pixels
    if attained == "minimum":
        counts = [average - 1] * (pixels - 1) + [average + pixels - 1]
    elif attained == "maximum":
        counts = [average - 1, average + 1] + [average] * (pixels - 2)
    else:
        counts = [average] * pixels
    assert sum(counts) == carriers and min(counts) > 0
    grid["counts_min"], grid["counts_max"] = min(counts), max(counts)
    assert retained.scale_ensemble(population)["angular_rows"] == 240
