"""Exact native current noise, a record-preparation response and its Ward law.

The current is the actual transition record j in x'=x-B_e*j. Its POVM
second moment is not the square of its mean-effect operator. All rates in
this module use admitted native proposals; modular s is named separately.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

import sympy as s

from ...observer import native
from ...native_modular_reference import local_preparation
from . import quantum_response


def _matrix_json(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


@lru_cache(maxsize=4)
def source(selector_rounds=1):
    """Build every native branch before using the reduced fluctuation law."""
    weights, timeout, survival = quantum_response.selector(selector_rounds)
    seams = tuple(sorted(weights))
    probabilities = tuple(s.Rational(weights[e]/survival) for e in seams)
    states = native.states(binary=True)
    index = {x: i for i, x in enumerate(states)}
    sigma = s.Matrix(list(map(s.Rational, local_preparation.reference())))
    x = s.Matrix(12, 66, lambda port, column: states[column][port])
    b = s.zeros(12, len(seams))
    t = s.zeros(66)
    mean_current = s.zeros(len(seams), 66)
    covariance0 = s.zeros(len(seams))
    after_current = s.zeros(66, len(seams))
    core_port = local_preparation.CORE[0]
    core_occupation = s.diag(*(state[core_port] for state in states))
    collective_x = x*core_occupation
    collective_x[core_port, :] = s.zeros(1, 66)
    collective_covariance0 = s.zeros(len(seams))
    collective_after_current = s.zeros(66, len(seams))
    first_effects, second_effects = [], []
    branch_count = 0
    for edge, ((a, z), probability) in enumerate(zip(seams, probabilities, strict=True)):
        b[a, edge], b[z, edge] = 1, -1
        first, second = [], []
        for i, state in enumerate(states):
            m1 = m2 = s.Rational(0)
            for coin in (0, 1):
                output, (_, current) = native.native(state, (a, z), coin)
                target = index[output]
                if any(output[p] != state[p]-b[p, edge]*current for p in range(12)):
                    raise ArithmeticError("native pathwise charge conservation failed")
                t[target, i] += probability/2
                mean_current[edge, i] += probability*current/2
                covariance0[edge, edge] += sigma[i]*probability*current**2/2
                after_current[target, edge] += sigma[i]*probability*current/2
                collective_covariance0[edge, edge] += sigma[i]*probability*state[core_port]*current**2/2
                collective_after_current[target, edge] += sigma[i]*probability*state[core_port]*current/2
                m1 += s.Rational(current, 2)
                m2 += s.Rational(current**2, 2)
                branch_count += 1
            first.append(m1); second.append(m2)
        first_effects.append(first); second_effects.append(second)
    d = s.diag(*probabilities)
    laplacian = b*d*b.T
    m = s.eye(12)-laplacian/2
    mean = x*sigma
    covariance = x*s.diag(*sigma)*x.T-mean*mean.T
    scale = covariance0[0, 0]/probabilities[0]
    collective_mean = collective_x*sigma
    collective_covariance = collective_x*s.diag(*sigma)*collective_x.T-collective_mean*collective_mean.T
    collective_scale = collective_covariance0[0, 0]/probabilities[0]
    collective_current = mean_current*core_occupation
    if (t*sigma != sigma or t != t.T or x*t != m*x
            or covariance0 != scale*d or covariance*b != scale*b
            or mean_current != d*b.T*x/2
            or x*after_current != -scale*b*d/2):
        raise ArithmeticError("native fluctuation/mean identities failed")
    if (collective_x*t != m*collective_x
            or s.ones(1, 12)*collective_x != s.ones(1, 66)*core_occupation
            or collective_covariance*b != collective_scale*b
            or collective_covariance0 != collective_scale*d
            or collective_current != d*b.T*collective_x/2
            or collective_x*collective_after_current != -collective_scale*b*d/2):
        raise ArithmeticError("collective local-charge current identities failed")
    return {"seams": seams, "probabilities": probabilities,
            "selector_survival": s.Rational(survival), "selector_timeout": s.Rational(timeout),
            "states": states, "reference": sigma, "occupations": x,
            "incidence": b, "weights": d, "laplacian": laplacian,
            "mean_map": m, "transition": t, "density_covariance": covariance,
            "current_mean_map": mean_current, "current_covariance0": covariance0,
            "after_current": after_current, "current_scale": scale,
            "collective_core_port": core_port,
            "collective_occupations": collective_x,
            "collective_current_mean_map": collective_current,
            "collective_current_covariance0": collective_covariance0,
            "collective_after_current": collective_after_current,
            "collective_current_scale": collective_scale,
            "first_effects": first_effects, "second_effects": second_effects,
            "native_branches_checked": branch_count}


def current_covariance(lag, selector_rounds=1, *, collective=False):
    """E[j_(n+lag) j_n^T] per surviving proposal, with zero equilibrium mean."""
    if type(lag) is not int or lag < 0 or type(collective) is not bool:
        raise ValueError("a nonnegative integer lag is required")
    data = source(selector_rounds)
    prefix = "collective_" if collective else ""
    if lag == 0:
        return data[prefix+"current_covariance0"]
    b, d, m = (data[key] for key in ("incidence", "weights", "mean_map"))
    return -data[prefix+"current_scale"]*d*b.T*m**(lag-1)*b*d/4


def full_state_covariance(lag, selector_rounds=1, *, collective=False):
    """Independent 66-state branch contraction for a positive lag."""
    if type(lag) is not int or not 1 <= lag <= 4 or type(collective) is not bool:
        raise ValueError("one through four admitted proposals required")
    data = source(selector_rounds)
    prefix = "collective_" if collective else ""
    return data[prefix+"current_mean_map"]*data["transition"]**(lag-1)*data[prefix+"after_current"]


def cycle_noise(selector_rounds=1):
    """Exact long-record covariance, including its finite-sum remainder."""
    data = source(selector_rounds)
    b, d, laplacian, m = (data[key] for key in ("incidence", "weights", "laplacian", "mean_map"))
    active = [i for i in range(12) if any(b[i, :])]
    n = len(active)
    local = laplacian.extract(active, active)
    constants = s.ones(n)/n
    inverse = (local+constants).inv()-constants
    if local*inverse != s.eye(n)-constants:
        raise ArithmeticError("the active complement graph is not connected")
    green = s.zeros(12)
    for i, a in enumerate(active):
        for j, z in enumerate(active):
            green[a, z] = inverse[i, j]
    c = d-d*b.T*green*b*d
    if c != c.T or b*c != s.zeros(12, d.rows) or c*d.inv()*c != c:
        raise ArithmeticError("weighted cycle covariance identity failed")
    scale = data["current_scale"]
    # Exact finite Green--Kubo remainder, not an assertion that a finite
    # correlation sum equals its infinite limit.
    partial = current_covariance(0, selector_rounds)
    remainders = []
    for lag in range(1, 4):
        gamma = current_covariance(lag, selector_rounds)
        partial += gamma+gamma.T
        remainder = scale*d*b.T*m**lag*green*b*d
        if partial-scale*c != remainder:
            raise ArithmeticError("finite current-noise remainder identity failed")
        remainders.append(str(s.trace(remainder)))
    return {"matrix": scale*c, "weighted_projector": c,
            "active_ports": active, "incidence_rank": b.rank(),
            "cycle_rank": c.rank(), "finite_remainder_traces": remainders,
            "positive_form_identity": "v^T C v=(C v)^T D^-1(C v), C=D-D B^T L^+ B D"}


def fixed_seam_protocol(sign=1, *, auxiliary_index=0):
    """Two continuing work events, with two native coin reads between them.

    The middle flag accepts with (1+sign*j0/2)/2. Both flag values and all
    current-zero records are retained. The same work state enters event two.
    This is an addressed preparation experiment, not an imposed force law.
    """
    if type(sign) is not int or sign not in (-1, 0, 1):
        raise ValueError("sign must be -1, 0 or 1")
    auxiliary_states = native.states()
    if type(auxiliary_index) is not int or not 0 <= auxiliary_index < len(auxiliary_states):
        raise ValueError("an auxiliary native-state index is required")
    seam = source()["seams"][0]
    auxiliary_seam = (0, 1)
    joint = defaultdict(F)
    terminal = {False: [F() for _ in range(66)], True: [F() for _ in range(66)]}
    states = native.states(binary=True)
    index = {state: i for i, state in enumerate(states)}
    branches = 0
    for probability, state in zip(local_preparation.reference(), states, strict=True):
        for first_coin in (0, 1):
            first, (_, j0) = native.native(state, seam, first_coin)
            for bit0 in (0, 1):
                aux1, record0 = native.native(auxiliary_states[auxiliary_index], auxiliary_seam, bit0)
                for bit1 in (0, 1):
                    _, record1 = native.native(aux1, auxiliary_seam, bit1)
                    accepted = 2*record0[0]+record1[0] < 2+sign*j0
                    for last_coin in (0, 1):
                        last, (_, j1) = native.native(first, seam, last_coin)
                        weight = probability/16
                        joint[j0, accepted, j1] += weight
                        terminal[accepted][index[last]] += weight
                        branches += 1
    masses = {flag: sum(state, F()) for flag, state in terminal.items()}
    means = {flag: sum((weight*j1 for (_, f, j1), weight in joint.items() if f == flag), F())/masses[flag]
             for flag in (False, True)}
    covariance = sum((weight*j0*j1 for (j0, _, j1), weight in joint.items()), F())
    if (sum(joint.values(), F()) != 1 or any(value != F(1, 2) for value in masses.values())
            or means[True] != sign*covariance/2 or means[False] != -sign*covariance/2):
        raise ArithmeticError("complete native record-preparation response failed")
    for i, probability in enumerate(local_preparation.reference()):
        if terminal[False][i]+terminal[True][i] != probability:
            raise ArithmeticError("forgetting the flag changed the stationary reference")
    return {"seam": list(seam), "bias": str(F(sign, 2)),
            "native_histories": branches, "work_events_per_history": 2,
            "auxiliary_native_events_per_history": 2,
            "auxiliary_fine_records_coarse_grained_after_flag": True,
            "complete_joint_current_flag_law": [
                {"first_current": j0, "accepted": flag, "second_current": j1, "probability": str(weight)}
                for (j0, flag, j1), weight in sorted(joint.items())],
            "acceptance_probability": str(masses[True]),
            "accepted_next_current_mean": str(means[True]),
            "rejected_next_current_mean": str(means[False]),
            "current_lag_one_covariance": str(covariance),
            "unconditional_signed_flag_current_mean": str((means[True]-means[False])/2),
            "reset_of_work_between_events": False,
            "conditioning_is_a_force_perturbation_or_conductivity_measurement": False}


def modular_dictionary():
    """Distinguish the transition-current instrument from a modular operator."""
    data = source()
    for first, second in zip(data["first_effects"], data["second_effects"], strict=True):
        if any(b != 2*a*a for a, b in zip(first, second, strict=True)):
            raise ArithmeticError("binary native-current moment relation changed")
    q = data["first_effects"][0]
    operator_variance = sum((w*a*a for w, a in zip(data["reference"], q, strict=True)), s.Rational(0))
    instrument_variance = 2*operator_variance
    return {"native_mean_effect": "Q_e=(n_a-n_b)/2",
            "native_second_moment_effect": "M2_e=2 Q_e^2",
            "operator_variance": str(operator_variance),
            "record_variance": str(instrument_variance),
            "instrument_noise_above_operator_variance": str(operator_variance),
            "modular_convention": "sigma_s(A)=sigma^(i s) A sigma^(-i s); s is the modular parameter",
            "exact_modular_charge_current": "i[log(sigma),n_a]=0 for every native port",
            "retarded_mean_current_commutator": "-i theta(s) Tr(sigma[Q_e(s),Q_f])=0 for all seams and s",
            "repair_proposal_current_is_modular_charge_derivative": False,
            "electromagnetic_conclusion": "Transition-record noise cannot be inserted as a photon propagator or Thomson pole residue without a different derived charge/current and spacetime dictionary."}


def report():
    data = source()
    noise = cycle_noise()
    for lag in (1, 2, 3):
        if full_state_covariance(lag) != current_covariance(lag):
            raise ArithmeticError("reduced current response differs from full native state histories")
        if full_state_covariance(lag, collective=True) != current_covariance(lag, collective=True):
            raise ArithmeticError("collective current response differs from full native state histories")
    budget = local_preparation.budget()
    error = budget["conditional_reference_uniform_halftrace_upper"]
    # Conservative dyadic enclosure, keeping the exact preparation receipt
    # as the authority without serializing enormous rational intermediates.
    scale = 2**60
    error_upper = F((error.numerator*scale+error.denominator-1)//error.denominator, scale)
    protocol = fixed_seam_protocol()
    signal = abs(F(protocol["unconditional_signed_flag_current_mean"]))
    full_error = budget["complete_flagged_reference_halftrace_upper"]
    full_upper = F((full_error.numerator*scale+full_error.denominator-1)//full_error.denominator, scale)
    herald = budget["uniform_ideal_herald_probability"]
    attempted_signal = herald*signal
    margin = attempted_signal/2-2*full_upper
    trials = 2_000_000
    exponent = trials*margin**2/2
    if not margin > 0 or not exponent > 5 or not sum((F(5**k, factorial(k)) for k in range(7)), F()) > 100:
        raise ArithmeticError("finite attempted-trial response guarantee failed")
    return {
        "schema": "oph.full_repair_response.native_current_response.v1",
        "reference": "same ideal faithful local_preparation.reference() with its finite preparation error contract",
        "native_branches_checked": data["native_branches_checked"],
        "seams": [list(edge) for edge in data["seams"]],
        "admitted_seam_probabilities": list(map(str, data["probabilities"])),
        "selector_rounds": 1,
        "selector_survival": str(data["selector_survival"]),
        "selector_terminal_timeout": str(data["selector_timeout"]),
        "current_noise_scale": str(data["current_scale"]),
        "identities": {
            "pathwise_continuity": "x_(n+1)=x_n-B j_n",
            "conditional_mean": "M=I-B D B^T/2",
            "lag_zero": "Gamma_0=d D",
            "positive_lag": "Gamma_k=-(d/4) D B^T M^(k-1) B D, k>=1",
            "actual_stopped_covariance": "E[J_k J_0^T]=q^(k+1) Gamma_k, with zero readout after terminal timeout",
            "dc_noise": "S(0)=d[D-D B^T(B D B^T)^+ B D]",
            "parameter": "admitted native work-proposal count; not modular s or seconds",
        },
        "current_covariance_lag_zero": _matrix_json(current_covariance(0)),
        "current_covariance_lag_one": _matrix_json(current_covariance(1)),
        "collective_local_charge_current": {
            "protected_core_port": data["collective_core_port"],
            "density": "Q_j=n_0*n_j for j!=0, Q_0=0; sum_j Q_j=n_0 on the binary N=2 carrier",
            "current": "J_e=n_0*j_e on complement seams",
            "pathwise_ward_identity": "Q_(n+1)=Q_n-B J_n, since complement repair preserves n_0",
            "current_noise_scale": str(data["collective_current_scale"]),
            "same_covariance_and_cycle_form_with_replaced_scale": True,
            "full_66_state_lags_checked": [1, 2, 3],
            "modular_charge_derivative": "i[log(sigma),Q_j]=0",
            "protected_bit_is_part_of_raw_current_record": False,
            "additional_readout_requirement": "Reading J_e requires the protected n_0 bit as well as the native seam-current record; the raw-current flag protocol does not silently supply that bit.",
            "physical_electromagnetic_identification": False},
        "cycle_noise": {"active_ports": noise["active_ports"],
                        "incidence_rank": noise["incidence_rank"], "cycle_rank": noise["cycle_rank"],
                        "exact_psd_identity": noise["positive_form_identity"],
                        "finite_sum_remainder_traces": noise["finite_remainder_traces"],
                        "dc_trace": str(s.trace(noise["matrix"])),
                        "infinite_covariance_belongs_to_survivor_normalized_diagnostic": True},
        "record_preparation_protocol": protocol,
        "unbiased_flag_control": fixed_seam_protocol(0),
        "reversed_flag_control": fixed_seam_protocol(-1),
        "finite_preparation_error": {
            "conditioned_on_preparation_success_only": True,
            "preparation_failure_and_timeout_flags_keep_their_own_original_masses": True,
            "conditional_reference_halftrace_upper": str(error_upper),
            "subsequent_complete_law_TV_upper": str(error_upper),
            "bounded_signed_current_mean_error_upper": str(2*error_upper),
            "ideal_signed_current_mean_magnitude": str(signal),
            "finite_preparation_bound_below_signal": 2*error_upper < signal,
            "new_finite_preparation_execution_claimed": False,
            "physical_service_or_sampling_error_included": False},
        "modular_operator_dictionary": modular_dictionary(),
        "finite_attempted_trial_contract": {
            "attempted_preparations": trials,
            "worst_case_native_attempts_including_preparation": trials*(budget["total_native_attempts_upper"]+4),
            "readout": "Y=(2*accept-1)*j1 on prepared success, Y=0 on every preparation failure/timeout",
            "range": [-1, 1], "ideal_biased_mean": str(-attempted_signal),
            "ideal_unbiased_control_mean": "0", "decision_threshold": str(-attempted_signal/2),
            "complete_preparation_halftrace_upper": str(full_upper),
            "conditional_mean_margin_lower": str(margin),
            "Hoeffding_exponent_lower": str(exponent),
            "each_direction_error_strictly_below": "1/100",
            "proof": "Each history-conditioned mean has margin m and Y lies in [-1,1]; the bounded-increment exponential bound is exp(-N*m^2/2). N*m^2/2>5 and sum(k=0..6)5^k/k!>100.",
            "exact_reset_or_actual_trial_independence_assumed": False,
            "preparation_failures_count_as_trials": True,
            "complete_sampled_campaign_executed": False,
            "native_attempt_count_is_physical_duration": False,
            "physical_service_errors_included": False,
            "C_and_E_native_laws_distinguished": False},
        "physical_scope": {
            "integer_charge_and_current_relative_normalization_from_native_continuity": True,
            "same_history_current_preparation_response_computed": True,
            "additional_fitted_mass_or_coupling": False,
            "empirical_inputs": [],
            "Ward_identity_selects_Maxwell_kinetic_coefficient": False,
            "Chain_scalar_identified_with_current_noise_or_photon_residue": False,
            "missing_electromagnetic_map": "A propagating conserved Lorentzian current coupled to a massless transverse field with derived kinetic normalization and a physical-clock pole/residue limit.",
            "normalization_freedom": "Charge conservation fixes relative source current units; multiplying a candidate Maxwell quadratic coefficient by any positive constant preserves that conservation law and changes its canonically normalized coupling.",
        },
    }
