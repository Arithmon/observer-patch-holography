"""Exact rational bounds for a supplied spatial Bose--Hubbard observer family.

Certificates concern ideal unitary splitting and fixed retained instruments.
They exclude floating execution, physical services, continuum removal and
physical clock calibration.  See DERIVATION.md for the hypotheses and proofs.
"""
from fractions import Fraction as F
from math import factorial


def _rational(value, name, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{qualifier} exact rational {name} required")
    return F(value)


def _integer(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError(f"integer {name} >= {minimum} required")
    return value


def _number(receipt, key, value):
    """Exact strings carry the certificate; floats are display values only."""
    value = F(value)
    receipt[f"{key}_exact"] = str(value)
    try:
        receipt[f"{key}_float"] = float(value)
    except OverflowError:
        receipt[f"{key}_float"] = None


def _scope(kind):
    return {
        "schema": f"oph.spatial_interaction.{kind}.v1",
        "particle_sector": "exactly two; invariant under the declared operations",
        "mass_common_phase_omitted": True,
        "floating_execution_error_included": False,
        "physical_service_error_included": False,
        "spatial_continuum_claimed": False,
        "physical_clock_calibration_claimed": False,
        "time_parameter": "supplied Hamiltonian model coordinate",
    }


def generator_norm_bound(hopping, interaction):
    """S=B=24J+U bounds H-2mI and the six-matching split norm sum."""
    j = _rational(hopping, "hopping")
    u = _rational(interaction, "repulsive interaction")
    return 24 * j + u


def strang_bound(steps, model_time, hopping, interaction):
    """All-input operator bound for a six-matching/onsite palindrome.

    Every matching kinetic generator has norm <=4J in the two-particle
    sector.  The interaction generator has norm <=U.  Integral Taylor
    remainders and unitary telescoping give S^3 t^3/(3 steps^2).
    """
    steps = _integer(steps, "split count")
    time = _rational(model_time, "model time")
    budget = generator_norm_bound(hopping, interaction)
    error = budget**3 * time**3 / (3 * steps**2)
    receipt = _scope("strang_bound")
    receipt.update({"steps": steps, "bound_is_uniform_in_volume": True,
                    "method": "degree-two agreement, integral Taylor remainder, unitary telescoping"})
    for key, value in (("model_time", time), ("generator_norm_sum_upper", budget),
                       ("unclipped_operator_error_upper", error),
                       ("operator_error_upper", min(F(2), error)),
                       ("channel_half_diamond_error_upper", min(F(1), error))):
        _number(receipt, key, value)
    return receipt


def history_bound(segments, hopping, interaction):
    """Bound the entire unpostselected record/output law, not each branch.

    ``segments`` is a finite sequence of (exact model time, split count).
    Preparations, retained instruments and feedback must be identical in
    both compared protocols.  Each segment's unitary bound is uniform in
    the incoming state, including an entangled retained record/reference.
    """
    generator_norm_bound(hopping, interaction)
    if not isinstance(segments, (tuple, list)):
        raise ValueError("a finite list or tuple of (model time, split count) pairs required")
    rows, total_time, total_error = [], F(0), F(0)
    for pair in segments:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError("each Hamiltonian segment requires (model time, split count)")
        time, steps = pair
        row = strang_bound(steps, time, hopping, interaction)
        rows.append(row)
        total_time += F(row["model_time_exact"])
        total_error += F(row["operator_error_upper_exact"])
    receipt = _scope("history_bound")
    receipt.update({"segments": rows,
                    "all_outcome_branches_retained": True,
                    "branch_count_multiplier": 1,
                    "trace_distance_convention": "one half of trace norm",
                    "bound_is_uniform_in_volume": True})
    for key, value in (("total_model_time", total_time),
                       ("stinespring_vector_error_upper", min(F(2), total_error)),
                       ("complete_instrument_trace_distance_upper", min(F(1), total_error)),
                       ("complete_history_tv_upper", min(F(1), total_error))):
        _number(receipt, key, value)
    return receipt


def boundary_bound(radius, total_model_time, hopping, interaction):
    """Centered [-R,R]^3 box versus larger boxes or infinite Z^3.

    All initially occupied sites lie in A={origin, origin+e_x}.  Every
    fixed instrument/intervention preserves particles outside A and maps
    the two-particle sector to itself.  It may mix configurations inside
    A, including both doublon/split Ramsey reads.  No postselection.
    The first potentially unequal Taylor degree is R, not R+1.
    """
    radius = _integer(radius, "box radius")
    time = _rational(total_model_time, "total model time")
    budget = generator_norm_bound(hopping, interaction)
    error = 2 * (budget * time)**radius / factorial(radius)
    receipt = _scope("boundary_bound")
    receipt.update({"radius": radius, "box_side": 2 * radius + 1,
                    "apparatus_sites": [[0, 0, 0], [1, 0, 0]],
                    "distance_from_apparatus_to_boundary": radius - 1,
                    "first_potentially_unequal_taylor_order": radius,
                    "comparison": "same fixed-spacing experiment on larger boxes or infinite cubic lattice",
                    "all_outcome_branches_retained": True,
                    "includes_arbitrary_nonlocal_interventions": False,
                    "temporal_splitting_error_included": False,
                    "method": "joint Stinespring Taylor agreement below boundary-reaching order"})
    for key, value in (("total_model_time", time), ("hamiltonian_norm_upper", budget),
                       ("unclipped_stinespring_vector_error_upper", error),
                       ("stinespring_vector_error_upper", min(F(2), error)),
                       ("complete_instrument_trace_distance_upper", min(F(1), error)),
                       ("complete_history_tv_upper", min(F(1), error))):
        _number(receipt, key, value)
    return receipt


def interaction_witness_bound(hopping, interaction):
    """Positive interaction-minus-free probability at t=1/(32(24J+U)).

    Prepare (|2_origin>+|1_origin,1_neighbor>)/sqrt(2).  Read the projector
    onto (|2_origin>+i|1_origin,1_neighbor>)/sqrt(2), with no intermediate
    occupation measurement.  The U-deletion control changes U only.
    This is a separate short protocol, not the measured-feedback history.
    """
    j = _rational(hopping, "hopping")
    u = _rational(interaction, "repulsive interaction", positive=True)
    budget = generator_norm_bound(j, u)
    time = F(1, 32) / budget
    first_order = u * time / 2
    remainder = 2 * u * time / 15
    lower = first_order - remainder
    receipt = _scope("interaction_witness_bound")
    receipt.update({"readout": "projector onto (doublon + i split)/sqrt(2)",
                    "control": "same preparation, hopping and readout with U=0",
                    "intermediate_measurement_included": False,
                    "bound_is_uniform_in_volume": True,
                    "finite_box_and_infinite_lattice_covered": True,
                    "temporal_splitting_error_included": False,
                    "witness_is_strictly_positive": True})
    for key, value in (("hopping", j), ("interaction", u),
                       ("hamiltonian_norm_upper", budget), ("model_time", time),
                       ("first_order_probability_difference", first_order),
                       ("remainder_absolute_upper", remainder),
                       ("interaction_minus_free_probability_lower", lower)):
        _number(receipt, key, value)
    return receipt


def gap_bound(mass, hopping, interaction):
    """Supplied full-Fock vacuum spectral gap, not a two-particle readout."""
    mass = _rational(mass, "vacuum gap parameter", positive=True)
    _rational(hopping, "hopping")
    _rational(interaction, "repulsive interaction")
    receipt = _scope("vacuum_gap")
    receipt.update({"particle_sector": "full bosonic Fock space",
                    "mass_common_phase_omitted": False,
                    "vacuum_is_unique": True,
                    "finite_open_box_gap_is_exact": True,
                    "infinite_cubic_lattice_spectral_threshold_is_exact": True,
                    "gap_generated_by_native_repair": False,
                    "gap_measured_by_fixed_two_particle_protocol": False,
                    "relativistic_mass_identification_claimed": False,
                    "method": "H >= mN, constant finite one-particle mode, infinite-volume wave packets"})
    _number(receipt, "vacuum_spectral_gap", mass)
    return receipt
