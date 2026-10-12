"""Complete occupation bound for the explicitly separate smooth instrument.

The prior interacting history is unchanged. Its final detector is replaced by
the Ramsey instrument with effects (I +/- sin(Q0))/2. The entire unnormalized
record/quantum output now has a finite-cutoff error bound, including the
pointer interaction. The old sharp-sign protocol and modular-clock attachment
are not covered. See SMOOTH_REFINEMENT.md for the proof and energy identity.
"""
from fractions import Fraction as F
import json

from ..observer.source_scalar.bounds import root_upper
from ..observer.source_scalar.cutoff import exponential_upper
from . import refinement as r
from . import scalar_coupling


POINTER_LINEAR_COEFFICIENT_UPPER = F(1, 2)
POINTER_MOMENT_RATE_UPPER = F(2)
POINTER_AUXILIARY_DURATION = F(1)
INFINITE_ENERGY_INJECTION = F(1, 8)


def smooth_occupation_bound(d, time=F(3, 2)):
    """Full ideal smooth instrument versus its infinite two-oscillator limit.

    The scalar preparation is the existing normalized total-number truncation.
    Both gauge measurement branches, conditional control, both Ramsey outcomes
    and their quantum states remain unnormalized. The bound is uniform over
    the three native-current gauge preparations and on/off gauge feedback.
    """
    d, t = r._dimension(d), r._time(time)
    scale = 10**30 * (d + 2)**6
    before = r.occupation_evolution_bound(d, t)
    projected_before = F(before["projected_history_vector_error_upper"])
    moment = F(before["final_fourth_number_moment_upper"])
    amplitude = root_upper(moment, scale=scale)
    # The auxiliary pointer parameter s runs from 0 to 1, independently of
    # the Hamiltonian model-time intervals. With C_pointer<=2, sqrt(M) grows
    # by at most exp(s); its integral is bounded by exp(1)-1.
    growth = exponential_upper(F(1))
    inverse_root = root_upper(F(1, d + 2), scale=scale)
    pointer_residual = amplitude * (growth - 1) * POINTER_LINEAR_COEFFICIENT_UPPER * inverse_root / (d + 2)
    projected = projected_before + pointer_residual
    outside = amplitude * growth / (d + 2)**2
    vector = root_upper(projected**2 + outside**2, scale=scale)
    distance = min(F(1), vector)
    return {
        "dimension": d,
        "interacting_model_time": str(t),
        "pointer_auxiliary_duration": str(POINTER_AUXILIARY_DURATION),
        "pointer_generator": "Z_ancilla Q0/2",
        "pointer_linear_coefficient_upper": str(POINTER_LINEAR_COEFFICIENT_UPPER),
        "pointer_fourth_moment_growth_rate_upper": str(POINTER_MOMENT_RATE_UPPER),
        "initial_fourth_number_moment_exact": str(r.initial_fourth_moment()),
        "pre_pointer_projected_history_error_upper": str(projected_before),
        "pointer_projected_duhamel_residual_integral_upper": str(pointer_residual),
        "final_projected_history_vector_error_upper": str(projected),
        "final_outside_vector_norm_upper": str(outside),
        "final_fourth_number_moment_upper": str(moment * growth**2),
        "complete_stinespring_vector_error_upper": str(vector),
        "complete_instrument_trace_distance_upper": str(distance),
        "trace_distance_convention": "one half of the trace norm of the state difference",
        "complete_classical_history_tv_error_upper": str(distance),
        "detector_probability_error_upper": str(distance),
        "includes_final_smooth_pointer_and_all_unnormalized_branches": True,
        "uniform_over_declared_native_currents_and_gauge_feedback_switch": True,
        "native_record_conditioned_continuation_preserves_bound": True,
        "bound_is_an_arbitrary_input_diamond_norm_claim": False,
        "includes_source_splitting_or_floating_execution": False,
        "includes_spatial_or_physical_service_error": False,
        "applies_to_old_sharp_sign_instrument_or_modular_clock": False,
        "scope": "complete ideal smooth instrument for the declared initial scalar state on fixed two-site interacting dynamics",
    }


def complete_error_budget(d, steps=512):
    """Add ideal finite temporal error at the actually implemented dimensions.

    The occupation term compares finite continuous evolution with the infinite
    history. The split term compares finite source product evolution with that
    same finite continuous history. Both comparisons use the smooth detector.
    Floating execution and physical-service errors are still excluded.
    """
    if type(d) is not int or d not in (2, 4):
        raise ValueError("implemented source dimensions 2 or 4 required")
    if type(steps) is not int or steps < 1:
        raise ValueError("positive integer split count required")
    occupation = smooth_occupation_bound(d)
    segments = [scalar_coupling.time_error_bound(d, steps, time=t, coupling=r.COUPLING)
                for t in (F(1), F(1, 2))]
    temporal = sum((F(item["ideal_channel_half_diamond_error_upper"]) for item in segments), F(0))
    cutoff = F(occupation["complete_instrument_trace_distance_upper"])
    combined = min(F(1), temporal + cutoff)
    return {
        "dimension": d, "steps_per_segment": steps,
        "occupation_certificate": occupation,
        "ideal_temporal_segment_certificates": segments,
        "ideal_temporal_instrument_error_upper": str(min(F(1), temporal)),
        "complete_ideal_product_to_infinite_instrument_error_upper": str(combined),
        "complete_classical_history_tv_error_upper": str(combined),
        "floating_execution_or_physical_service_error_included": False,
        "nonvacuous_combined_bound": combined < 1,
    }


def status():
    return {
        "protocol": "separate smooth Ramsey instrument; old sharp-sign history preserved",
        "full_instrument_occupation_error_rate": "O(1/d) at fixed total model time",
        "same_history_ideal_time_and_occupation_limit": "n(d)=d^2 steps per model-time segment suffices",
        "infinite_unconditional_energy_injection_exact": str(INFINITE_ENERGY_INJECTION),
        "energy_identity_scope": "normalized input in the infinite Hamiltonian quadratic-form domain; retain both outcomes",
        "finite_cutoff_canonical_energy_identity_claimed": False,
        "finite_post_measurement_fourth_moment_preserved_as_finite": True,
        "conditional_branch_energy_increment_of_one_eighth_claimed": False,
        "compiled_dimensions": [2, 4],
        "larger_occupation_source_execution_claimed": False,
        "old_sharp_readout_gap_closed_by_this_result": False,
        "spatial_continuum_or_empirical_completion": False,
    }


def report(dimensions=(4, 16, 64, 256), steps=512):
    return {"status": status(),
            "occupation_bounds": [smooth_occupation_bound(d) for d in dimensions],
            "executed_dimension_ideal_budgets": [complete_error_budget(d, steps) for d in (2, 4)]}


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, sort_keys=True))
