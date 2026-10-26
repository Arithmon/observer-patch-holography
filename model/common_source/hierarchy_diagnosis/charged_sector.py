"""Native charged response on the adjacent binary load sectors.

The total-two source has no internal occupation-changing CAR operator.
This diagnostic extends the *existing endpoint rule*, without changing it,
to binary loads zero through three. That is enough for the two insertions
of a charged field about the total-two reference; it is not the full Fock
space and is not claimed to select a physical fermion representation.

Exact current records preserve the signed difference of a matrix unit.
Consequently enlarging the load domain alone cannot give this channel a
propagating native-port charged response. Coherent recombination of distinct
current records and encoded/collective fields are outside this conclusion.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb

from ...native_modular_reference import local_preparation
from ...observer import native, quantum_mixing
from ...repair_scheduler import coins


PORTS = 12
SOURCE_LOAD = 2


def occupation(mask):
    return tuple((mask >> port) & 1 for port in range(PORTS))


def mask_of(loads):
    if len(loads) != PORTS or any(value not in (0, 1) for value in loads):
        raise ValueError("twelve binary native loads required")
    return sum(value << port for port, value in enumerate(loads))


def masks(load):
    if type(load) is not int or not 0 <= load <= 3:
        raise ValueError("a diagnostic load sector zero through three required")
    return tuple(sum(1 << port for port in occupied)
                 for occupied in combinations(range(PORTS), load))


def creation(mask, port):
    """The declared ordered CAR completion, including its spectator sign."""
    if type(port) is not int or not 0 <= port < PORTS:
        raise ValueError("native port zero through eleven required")
    if mask & (1 << port):
        return None
    sign = (-1) ** ((mask & ((1 << port)-1)).bit_count())
    return mask | (1 << port), sign


def annihilation(mask, port):
    if type(port) is not int or not 0 <= port < PORTS:
        raise ValueError("native port zero through eleven required")
    if not mask & (1 << port):
        return None
    sign = (-1) ** ((mask & ((1 << port)-1)).bit_count())
    return mask ^ (1 << port), sign


@lru_cache(maxsize=128)
def event_map(law, seam, coin):
    """Actual native outputs and retained C/E labels on the adjacent sectors."""
    if law not in ("C", "E") or type(coin) is not int or coin not in (0, 1):
        raise ValueError("native C/E law and binary placement coin required")
    if seam not in quantum_mixing.seams():
        raise ValueError("an existing native seam required")
    result = {}
    for load in range(4):
        for mask in masks(load):
            state = occupation(mask)
            output, record = native.native(state, seam, coin)
            target = mask_of(output)
            if target.bit_count() != load or record[0] != coin:
                raise ArithmeticError("native primitive changed load or coin record")
            label = record[1] if law == "C" else state[seam[0]]-state[seam[1]]
            result[mask] = target, label
    return result


def source_reference():
    return {mask_of(state): probability for state, probability in
            zip(native.states(binary=True), local_preparation.reference(), strict=True)}


def field_insertion(port):
    """Sparse {c_port†,sigma}; only N=1,2,3 are visited by the correlator."""
    result = defaultdict(F)
    for mask, probability in source_reference().items():
        added = creation(mask, port)
        if added is not None:
            target, sign = added
            result[target, mask] += sign*probability
        removed = annihilation(mask, port)
        if removed is not None:
            target, sign = removed
            result[mask, target] += sign*probability
    return dict(result)


def apply_event(operator, law, seam, coin):
    """One unweighted recorded-current sum; the coin probability is separate."""
    mapping = event_map(law, seam, coin)
    result = defaultdict(F)
    for (left, right), value in operator.items():
        a, label_a = mapping[left]
        b, label_b = mapping[right]
        if label_a == label_b:
            before = tuple(x-y for x, y in zip(occupation(left), occupation(right), strict=True))
            after = tuple(x-y for x, y in zip(occupation(a), occupation(b), strict=True))
            if before != after:
                raise ArithmeticError("a recorded native branch changed charge support")
            result[a, b] += value
    return dict(result)


def read_field(operator, port):
    """Exact trace(c_port operator), with all spectator CAR signs retained."""
    value = F()
    for (left, right), coefficient in operator.items():
        removed = annihilation(left, port)
        if removed is not None and removed[0] == right:
            value += removed[1]*coefficient
    return value


def response(law, *, selector_rounds=1):
    """One actual completion proposal, with accepted and terminal outputs apart.

    This is a regression correlator of the specified CP continuation. It is
    not an executed charged-field detector or a physical rest-mass claim.
    In particular a terminal timeout is never restarted.
    """
    if law not in ("C", "E"):
        raise ValueError("native C/E law required")
    selection = local_preparation.selectors(selector_rounds)
    weights = selection["completion_accepted_seam_weights"]
    timeout = selection["completion_timeout_probability"]
    accepted = [[F() for _ in range(PORTS)] for _ in range(PORTS)]
    terminal = [[F() for _ in range(PORTS)] for _ in range(PORTS)]
    for source_port in range(PORTS):
        inserted = field_insertion(source_port)
        output = defaultdict(F)
        for seam, weight in sorted(weights.items()):
            for coin in (0, 1):
                for key, value in apply_event(inserted, law, seam, coin).items():
                    output[key] += weight*value/2
        for detector_port in range(PORTS):
            accepted[detector_port][source_port] = read_field(output, detector_port)
            terminal[detector_port][source_port] = timeout*read_field(inserted, detector_port)
    return accepted, terminal


def primitive_certificate():
    """Map the extension back to every retained total-two native branch."""
    states = native.states(binary=True)
    input_addresses = 0
    for law in ("C", "E"):
        for seam in quantum_mixing.seams():
            maps = coins.coin_maps(law, seam)
            expected = {(coin, mask_of(states[before])): (mask_of(states[after]), label)
                        for coin, label, mapping in maps for before, after in mapping}
            for coin in (0, 1):
                actual = event_map(law, seam, coin)
                for state in states:
                    mask = mask_of(state)
                    if actual[mask] != expected[coin, mask]:
                        raise ArithmeticError("adjacent extension changed a retained source branch")
                    input_addresses += 1
    return input_addresses


def reference_extension(eta):
    """Faithful full-binary references with the same conditioned source state.

    Per-configuration weight is w_q eta**(N-2)**2, q the protected-core
    occupation. This is an explicit witness family, not a selected ensemble.
    Exact binomial counting avoids materializing the 4096-state density.
    """
    eta = F(eta)
    if eta <= 0:
        raise ValueError("a strictly positive fugacity factor required")
    actual = source_reference()
    core_weights = {}
    for mask, probability in actual.items():
        core = sum(occupation(mask)[port] for port in local_preparation.CORE)
        if core in core_weights and core_weights[core] != probability:
            raise ArithmeticError("source weight depends on more than protected-core occupation")
        core_weights[core] = probability
    weights = tuple(core_weights[core] for core in range(3))
    sector_weights = []
    for load in range(PORTS+1):
        weight = sum((comb(2, core)*comb(10, load-core)*weights[core]
                      for core in range(3) if 0 <= load-core <= 10), F())
        sector_weights.append(weight*eta**((load-SOURCE_LOAD)**2))
    partition = sum(sector_weights, F())
    if sector_weights[SOURCE_LOAD] != 1:
        raise ArithmeticError("extended reference failed to preserve source normalization")
    for mask, probability in actual.items():
        core = sum(occupation(mask)[port] for port in local_preparation.CORE)
        if probability != weights[core]:
            raise ArithmeticError("extended reference changed a source conditional weight")
    return {
        "eta_exact": str(eta),
        "full_binary_dimension": 2**PORTS,
        "partition_function_exact": str(partition),
        "number_sector_probabilities_exact": [str(weight/partition) for weight in sector_weights],
        "conditional_total_two_reference_identical_exact": True,
        "per_configuration_weight": "w_q eta^((N-2)^2), q=protected-core occupation",
        "modular_generator_difference_up_to_scalar": "-log(eta)*(N-2)^2",
        "two_copy_pair_frequency_shift": "-2*log(eta)",
    }


def report():
    certificates = primitive_certificate()
    rows = {}
    for law in ("C", "E"):
        accepted, terminal = response(law)
        if any(accepted[i][j] != 0 or terminal[i][j] != 0
               for i in range(PORTS) for j in range(PORTS) if i != j):
            raise ArithmeticError("native charged response unexpectedly propagated")
        rows[law] = {
            "accepted_response_diagonal_exact": [str(accepted[i][i]) for i in range(PORTS)],
            "terminal_response_diagonal_exact": [str(terminal[i][i]) for i in range(PORTS)],
            "off_port_entries_zero_exact": True,
        }
    return {
        "schema": "oph.hierarchy_diagnosis.charged_sector.v1",
        "native_rule": "model.observer.native.native, unchanged on binary loads zero through three",
        "sector_dimensions": [len(masks(load)) for load in range(4)],
        "retained_total_two_branch_addresses_checked": certificates,
        "field_representation": "ordered twelve-mode CAR restricted to insertions about total-two; no physical identification",
        "one_proposal": {"selector_rounds": 1, "outcomes": rows,
                         "terminal_timeout_restarted": False,
                         "accepted_outcomes_conditioned_on_success": False},
        "exact_support_theorem": {
            "identity": "x'=x-j(e_a-e_b), y'=y-j(e_a-e_b) implies x'-y'=x-y",
            "domain": "all nonnegative integer loads for C/E records with the same placement coin",
            "composition": "every finite recorded seam history and record-controlled schedule",
            "diagonal_reference_cross_port_charged_two_point_response_zero": True,
            "coherent_recombination_of_distinct_current_records_covered": False,
            "arbitrary_encoded_collective_fields_excluded": False,
        },
        "sector_phase_ambiguity": {
            "conditional_witness": "H_kappa=H_0+kappa(N-2)^2 with [H_0,N]=0",
            "faithful_reference_witnesses": [reference_extension(F(1)), reference_extension(F(1, 2))],
            "stationarity": "completion seams preserve N and protected-core q; each local population swap joins equal weights",
            "all_total_two_histories_and_clocks_unchanged": True,
            "addition_and_removal_energy_shift": "kappa each",
            "two_copy_charge_conserving_pair_frequency_shift": "2*kappa",
            "common_chemical_potential_contribution_to_pair_frequency": "0",
            "admissibility_as_local_Lorentzian_OPH_dynamics_proved": False,
            "purpose": "same-shell data cannot select an unmeasured coherent sector extension",
        },
        "native_port_charged_propagation_derived": False,
        "Standard_Model_mass_prediction": False,
        "empirical_inputs": [],
    }
