"""Transfer the sharp-code classification without changing the native law.

The proof is in README.md. Positive endpoint coins rescale retained effects;
they do not change exact flagged correctability. Every admitted sharp-q=1
qubit therefore has core-occupation compression I2/2. This derives the core
kernel independently of any particular encoder or recovery implementation.
"""
from collections import Counter
from fractions import Fraction as F

from ..local_repair_family import modular_clock
from ..moving_defect import transport
from ..biased_reset import reset


def kernel(parameters=None):
    """All60 legal rows, all66 destinations, including forbidden arrivals.

    A hold is one attempted repair with no move. The seam scheduler and
    endpoint weights are fixed and source independent. The formula permits
    history-dependent exact code choices and retained C/E choices.
    """
    p = modular_clock._parameters(transport.parameters() if parameters is None else parameters)
    total = sum(p.seam_weights)
    answer = []
    for core in transport.LIVE:
        row = Counter()
        for seam, seam_weight in zip(transport.SEAMS, p.seam_weights, strict=True):
            occupied_core = set(core) & set(seam)
            if len(occupied_core) != 1:
                continue
            a = next(iter(occupied_core))
            b = next(v for v in seam if v != a)
            destination = tuple(sorted((set(core) - {a}) | {b}))
            row[destination] += (p.repair_probability * F(seam_weight, total)
                                 * F(1, 2) * F(p.port_weights[b],
                                              p.port_weights[a] + p.port_weights[b]))
        row[core] = 1 - sum(row.values())
        answer.append(tuple(row[c] for c in transport.CORES))
    return tuple(answer)


def reset_scope_certificate():
    """Attach the existing exact potential to the universal biased kernel.

    This strengthens the code-choice quantifier only. The unchanged source
    law, stopping set, reset cost and ideal recovery premises are essential.
    The biased all-finite-horizon optimum remains unproved.
    """
    p, full, _, _, h, maximum, best = reset.exact_chain()
    derived = kernel(p)
    if derived != full:
        raise ArithmeticError("all-code kernel differs from the retained biased process")
    if derived != transport.kernel(p, "E"):
        raise ArithmeticError("retained E law has a different all-code kernel")
    for i, row in enumerate(derived):
        drift = sum(row[transport.CORES.index(c)] * h[j]
                    for j, c in enumerate(transport.LIVE))
        if drift != h[i] - 1:
            raise ArithmeticError("universal code kernel does not satisfy the exact potential")
    return {
        "legal_cores": len(transport.LIVE),
        "excluded_antipodal_cores": len(transport.STOP),
        "all_complete_C_E_core_kernel_entries_checked": 2 * 60 * 66,
        "maximum_mean_excursion_native_attempts": str(maximum),
        "maximizing_reset_cores": [list(c) for c in best],
        "optimal_asymptotic_resets_per_native_attempt": str(1 / maximum),
        "adaptive_arbitrary_exact_source_qubit_codes_allowed": True,
        "arbitrary_untouched_reference_allowed": True,
        "positive_endpoint_coin_weights_required": True,
        "all_thirty_seams_exactly_correctable_at_every_ready_boundary": True,
        "core_resets_counted_including_preemptive_and_randomized_choices": True,
        "persistent_auxiliary_logical_storage_allowed": False,
        "all_finite_horizon_biased_optimality_proved": False,
        "physical_clock_energy_mass_or_particle_identification": False,
    }
