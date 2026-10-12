"""Small dimensionless readout benchmarks, not calibration against nature.

Clock samples are synthetic Bernoulli trials of the analytic modular clock law,
checked separately against the actual source probabilities. Geometry compares
two representations of the same native response carrier. Energy uses independent
finite oscillator matrices, including larger
cutoffs that have not been executed through the source compiler.
"""
import json
import math
import random

import numpy as np

from . import experiment
from . import modular_clock
from . import slow_geometry as geometry


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"positive integer {name} required")
    return value


def binomial_parameter_interval(count, shots, clock, settings=5, failure_probability=.05):
    """Invert a simultaneous Hoeffding interval on the declared first window.

    The input probability is normalized, whereas source receipts retain total
    weights. Intersect with the model window before inversion; an empty
    intersection is reported explicitly, never forced into an interval.
    """
    shots = _positive_integer(shots, "shot count")
    settings = _positive_integer(settings, "setting count")
    if type(count) is not int or not 0 <= count <= shots:
        raise ValueError("integer count between zero and shots required")
    if type(failure_probability) not in (int, float) or not 0 < failure_probability < 1:
        raise ValueError("failure probability strictly between zero and one required")
    contrast = clock["contrast"] / clock["mass"]
    frequency = clock["modular_frequency"]
    if not 0 < contrast <= 1 or not math.isfinite(frequency) or not frequency > 0:
        raise ValueError("nontrivial normalized modular clock required")
    window = clock.get("monotone_parameter_window")
    if (not isinstance(window, (tuple, list)) or len(window) != 2
            or any(type(value) not in (int, float) or not math.isfinite(value) for value in window)
            or not 0 <= window[0] < window[1]
            or frequency * window[1] > math.pi / 2 + 2e-15):
        raise ValueError("declared increasing first-branch modular parameter window required")
    observed = count / shots
    radius = math.sqrt(math.log(2 * settings / failure_probability) / (2 * shots))
    probability_interval = [max(0., observed - radius), min(1., observed + radius)]
    allowed = [(1 + contrast * math.sin(frequency * endpoint)) / 2 for endpoint in window]
    lower, upper = max(allowed[0], probability_interval[0]), min(allowed[1], probability_interval[1])

    def inverse(probability):
        if probability <= allowed[0]:
            return window[0]
        if probability >= allowed[1]:
            return window[1]
        value = min(1., max(0., (2 * probability - 1) / contrast))
        # Clamp endpoint roundoff to the actual declared parameter interval,
        # which can be narrower than the full first monotone sine branch.
        return min(window[1], max(window[0], math.asin(value) / frequency))

    estimate_probability = min(allowed[1], max(allowed[0], observed))
    return {
        "shots": shots, "number_one_count": count, "observed_frequency": observed,
        "hoeffding_probability_radius": radius,
        "normalized_probability_interval": probability_interval,
        "declared_parameter_window": list(window),
        "declared_normalized_probability_window": allowed,
        "inferred_modular_parameter": inverse(estimate_probability),
        "point_estimate_projected_to_declared_window": estimate_probability != observed,
        "modular_parameter_interval": [inverse(lower), inverse(upper)] if lower <= upper else None,
        "interval_intersects_declared_model_range": lower <= upper,
    }


def clock_benchmark(steps=512, shots=100000, seed=20261012):
    steps = _positive_integer(steps, "split count")
    shots = _positive_integer(shots, "shot count")
    if type(seed) is not int or seed < 0:
        raise ValueError("nonnegative integer sampling seed required")
    history = experiment.history(2, steps)
    emitted = modular_clock.receipt(history)
    clock = emitted["clock_parameters"]
    rng = random.Random(seed)
    settings = []
    for row in emitted["settings"]:
        mass = sum(row["clock_bit_marginal"])
        probability = row["clock_bit_marginal"][1] / mass
        sampling_probability = (1 + clock["contrast"] / clock["mass"] * math.sin(row["phase"])) / 2
        # Sample the checked analytic law with explicit uniform thresholds.
        # A binomial implementation can reflect its seeded draw at p > 1/2;
        # one-ulp source roundoff at the zero-phase setting then changes the
        # retained count macroscopically. Keep the actual source probability
        # separate, without rounding it or changing the replay tolerance.
        count = sum(rng.random() < sampling_probability for _ in range(shots))
        inferred = binomial_parameter_interval(count, shots, clock, len(emitted["settings"]))
        interval = inferred["modular_parameter_interval"]
        truth = row["modular_parameter"]
        settings.append({
            "programmed_modular_parameter": truth, "programmed_modular_phase": row["phase"],
            "normalized_source_probability": probability,
            "synthetic_sampling_probability": sampling_probability,
            "source_retained_total_weight": mass,
            "source_law_probability_error": row["analytic_probability_error"],
            **inferred,
            "absolute_parameter_estimation_error": abs(inferred["inferred_modular_parameter"] - truth),
            "programmed_parameter_inside_interval": interval is not None and interval[0] <= truth <= interval[1],
        })
    return {
        "dimension": 2, "steps_per_segment": steps, "shots_per_setting": shots,
        "seed": seed, "settings": settings,
        "simultaneous_confidence_lower_bound": .95,
        "confidence_method": "two-sided Hoeffding with union bound over all five independent settings",
        "monotone_modular_parameter_window": clock["monotone_parameter_window"],
        "state_derived_modular_frequency": clock["modular_frequency"],
        "state_derived_normalized_contrast": clock["contrast"] / clock["mass"],
        "all_programmed_parameters_inside_reported_intervals": all(row["programmed_parameter_inside_interval"] for row in settings),
        "all_source_instrument_checks_pass": emitted["all_pass"],
        "sampling_is_synthetic": True,
        "sampling_method": "Python Random.random uniform-threshold Bernoulli trials of the checked analytic law",
        "fresh_complete_history_per_shot_and_setting": True,
        "sequential_use_of_one_system": False,
        "confidence_covers_hardware_or_state_estimation_error": False,
        "physical_seconds_calibrated": False,
    }


def _q5_float(value):
    return float(value[0]) + float(value[1]) * math.sqrt(5)


def distance_benchmark():
    source = geometry.source()
    # Run the existing exact native/response binding before reading distances.
    existing = geometry.calibration()
    projector = source.slow_band_projector_exact()
    gram = source.intrinsic_gram_exact()
    frame = geometry.raw_frame()
    radius_squared = source.q5_add(source.GOLDEN, source.q5(2))
    inverse_radius = source.q5_inv(radius_squared)
    seams = set(source.seams())
    pairs, largest_error = [], 0.
    for i in range(12):
        for j in range(12):
            if source.q5_mul(source.q5(4), projector[i][j]) != gram[i][j]:
                raise AssertionError("response projector and normalized Gram differ")
        for j in range(i + 1, 12):
            distance_squared = source.q5_mul(source.q5(4), source.q5_sub(
                source.q5_add(projector[i][i], projector[j][j]), source.q5_mul(source.q5(2), projector[i][j])))
            displacement = geometry._subtract(frame[j], frame[i])
            coordinate_squared = source.q5_mul(geometry._dot(displacement, displacement), inverse_radius)
            if distance_squared != coordinate_squared:
                raise AssertionError("exact projector and coordinate distances differ")
            coordinate_length = float(np.linalg.norm([_q5_float(x) for x in displacement])) / math.sqrt(_q5_float(radius_squared))
            response_length = math.sqrt(_q5_float(distance_squared))
            largest_error = max(largest_error, abs(coordinate_length - response_length))
            pairs.append({"ports": [i, j], "native_edge": (i, j) in seams,
                          "squared_distance_Qsqrt5": source.q5_str(distance_squared),
                          "response_distance": response_length,
                          "coordinate_distance": coordinate_length})
    edge_squared = source.q5_mul(source.q5(4), inverse_radius)
    checked = 0
    for seam in source.seams():
        for current in (-1, 0, 1):
            displacement = geometry.current_displacement(seam, current)
            squared = source.q5_mul(geometry._dot(displacement, displacement), inverse_radius)
            if squared != source.q5_mul(source.q5(current**2), edge_squared):
                raise AssertionError("actual signed current displacement has the wrong edge length")
            checked += 1
    return {
        "pairs": pairs, "exact_pair_comparisons": len(pairs),
        "native_edge_count": len(seams), "native_current_displacements_checked_exactly": checked,
        "expected_raw_icosahedron_edge_length": 2.,
        "expected_unit_radius_edge_squared_Qsqrt5": source.q5_str(edge_squared),
        "expected_unit_radius_edge_length": math.sqrt(_q5_float(edge_squared)),
        "maximum_numerical_coordinate_distance_error": largest_error,
        "exact_projector_gram_and_coordinate_agreement": True,
        "existing_native_calibration_checks_pass": existing["all_pass"],
        "distance_unit": "dimensionless unit-radius response carrier",
        "physical_length_calibrated": False,
        "three_dimensional_spatial_continuum_established": False,
    }


def energy_benchmark(dimensions=(2, 4, 8, 16, 32)):
    rows = []
    for d in dimensions:
        if type(d) is not int or d < 2 or d % 2:
            raise ValueError("even oscillator cutoff at least two required")
        q = np.diag(np.sqrt(np.arange(1, d) / 2), 1)
        q += q.T
        eigenvalues, basis = np.linalg.eigh(q)
        plus = (basis * np.exp(-.5j * eigenvalues)) @ basis.T
        minus = plus.conj().T
        maps = ((plus - 1j * minus) / 2, (plus + 1j * minus) / 2)
        vacuum = np.eye(d)[:, 0]
        energies = np.arange(d) + .5
        vectors = [operator @ vacuum for operator in maps]
        probabilities = [float(np.vdot(vector, vector).real) for vector in vectors]
        final_energy = float(sum(np.dot(abs(vector)**2, energies) for vector in vectors))
        increment = final_energy - .5
        rows.append({"dimension": d, "initial_energy": .5,
                     "unconditional_final_energy": final_energy,
                     "unconditional_energy_increment": increment,
                     "infinite_canonical_increment_exact": "1/8",
                     "absolute_error_from_infinite_increment": abs(increment - .125),
                     "outcome_probabilities": probabilities,
                     "all_outcomes_retained": True,
                     "execution_kind": "independent finite oscillator matrix calculation",
                     "source_compiler_execution": False})
    return {"initial_state": "unit-frequency oscillator vacuum",
            "energy_operator": "compressed N+1/2", "rows": rows,
            "detector": "same smooth Ramsey Kraus functions of compressed Q",
            "physical_energy_calibrated": False,
            "common_interacting_vacuum_used_in_this_component_benchmark": False}


def report(steps=512, shots=100000, seed=20261012):
    clock = clock_benchmark(steps, shots, seed)
    distances = distance_benchmark()
    energy = energy_benchmark()
    checks = {"emitted_clock_law_passes": clock["all_source_instrument_checks_pass"],
              "native_geometry_exact_comparison_passes": distances["existing_native_calibration_checks_pass"],
              "energy_keeps_both_normalized_outcomes": all(abs(sum(row["outcome_probabilities"]) - 1) < 2e-12 for row in energy["rows"]),
              "large_cutoff_energy_matches_infinite_value_numerically": energy["rows"][-1]["absolute_error_from_infinite_increment"] < 2e-12}
    if not all(checks.values()):
        raise AssertionError(checks)
    return {"schema": "oph.common_source.dimensionless_calibration_benchmarks.v1",
            "clock_sampling": clock, "response_distances": distances,
            "smooth_energy_cutoff": energy, "checks": checks, "all_pass": True,
            "scope": "synthetic dimensionless readout benchmarks; no experimental data, physical units, autonomous clock or agreement with nature"}


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, sort_keys=True))
