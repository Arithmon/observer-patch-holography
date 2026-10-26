"""Propagating collective charge transfer in the actual recorded native law.

X_j=|a,j><b,j| moves its common occupation j while its signed endpoints a,b
stay fixed. This is an even matrix-unit field, not a multi-mode CAR field.
The complete native records, first incident clock, and terminal outcomes are
kept. All coefficients come from the existing completion selector/reference.
No coherent hopping Hamiltonian or physical spacetime identification is added.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math

import numpy as np
import sympy as s

from ...native_modular_reference import local_preparation as preparation
from ...observer import native
from ...repair_scheduler import coins


STATES = native.states(binary=True)
INDEX = {state: index for index, state in enumerate(STATES)}


def _pair(pair):
    if (type(pair) is not tuple or len(pair) != 2
            or any(type(port) is not int or not 0 <= port < 12 for port in pair)
            or pair[0] == pair[1]):
        raise ValueError("two distinct ordered native ports required")
    return pair


def address(a, b):
    if a == b:
        raise ValueError("two different occupied ports required")
    return INDEX[tuple(int(port in (a, b)) for port in range(12))]


def fields(pair=(0, 5)):
    a, b = _pair(pair)
    return tuple((j, address(a, j), address(b, j)) for j in range(12) if j not in pair)


@lru_cache(maxsize=8)
def events(law, selector_rounds=1):
    if law not in ("C", "E"):
        raise ValueError("native C or E instrument required")
    packet = preparation.selectors(selector_rounds)
    rows = []
    for seam, weight in sorted(packet["completion_accepted_seam_weights"].items()):
        for coin, label, mapping in coins.coin_maps(law, seam):
            rows.append(((*seam, coin, label), F(weight)/2, dict(mapping)))
    return tuple(rows)


def _matrix(rows):
    return [[str(value) for value in row] for row in rows.tolist()]


@lru_cache(maxsize=8)
def transfer(law, pair=(0, 5), selector_rounds=1):
    """Exact accepted transfer and its incident/remote decomposition."""
    _pair(pair)
    modes = fields(pair)
    where = {(x, y): index for index, (_, x, y) in enumerate(modes)}
    incident, remote = s.zeros(10), s.zeros(10)
    cap = set(preparation.CAPS[0])
    for record, probability, mapping in events(law, selector_rounds):
        target = incident if cap.intersection(record[:2]) else remote
        for column, (_, x, y) in enumerate(modes):
            if x in mapping and y in mapping:
                output = mapping[x], mapping[y]
                if output not in where:
                    raise ArithmeticError("collective signed-difference block is not closed")
                target[where[output], column] += s.Rational(probability)
    accepted = incident+remote
    weights = preparation.selectors(selector_rounds)["completion_accepted_seam_weights"]
    q = s.Rational(sum(weights.values(), F()))
    tau = 1-q
    a, b = pair
    delta = tuple(int(port == a)-int(port == b) for port in range(12))
    killing = sum((s.Rational(weight)*(
        s.Rational(abs(delta[u]-delta[v]), 2) if law == "C" else int(delta[u] != delta[v]))
        for (u, v), weight in weights.items()), s.Rational(0))
    common = [j for j, _, _ in modes]
    laplacian = s.zeros(10)
    for (u, v), weight in weights.items():
        if u in common and v in common:
            i, j = common.index(u), common.index(v)
            laplacian[i, i] += s.Rational(weight)
            laplacian[j, j] += s.Rational(weight)
            laplacian[i, j] -= s.Rational(weight)
            laplacian[j, i] -= s.Rational(weight)
    if accepted != (q-killing)*s.eye(10)-laplacian/2:
        raise ArithmeticError("native branch transfer disagrees with the independent graph formula")
    reference = preparation.reference()
    ratios = tuple(reference[x]/reference[y] for _, x, y in modes)
    for i in range(10):
        for j in range(10):
            if accepted[i, j] and ratios[i] != ratios[j]:
                raise ArithmeticError("collective transport changed modular frequency")
    return {"fields": modes, "accepted": s.ImmutableMatrix(accepted),
            "incident": s.ImmutableMatrix(incident), "remote": s.ImmutableMatrix(remote),
            "laplacian": s.ImmutableMatrix(laplacian), "survival": q,
            "timeout": tau, "killing": killing, "modular_weight_ratios": ratios}


def algebra(pair=(0, 5)):
    modes = fields(pair)
    core = set(preparation.CORE)
    charges = [sum(STATES[x][p]-STATES[y][p] for p in core) for _, x, y in modes]
    if len(set(charges)) != 1:
        raise ArithmeticError("collective field lacks a fixed protected-core charge")
    # Pairwise-disjoint two-state supports establish the displayed matrix-unit
    # products exactly. No numerical rank or assumed statistics enters.
    supports = [{x, y} for _, x, y in modes]
    if len(set.union(*supports)) != 20:
        raise ArithmeticError("collective two-level supports overlap")
    if any(STATES[x][p] != STATES[y][p] for _, x, y in modes for p in range(12) if p not in pair):
        raise ArithmeticError("field insertion changed its common occupation")
    return {
        "field": "X_j=|a,j><b,j|",
        "pair": list(pair), "common_occupation_sites": [j for j, _, _ in modes],
        "conserved_core_number_charge": charges[0],
        "total_native_number_charge": 0,
        "products": "X_j X_k^dagger=delta_jk |a,j><a,j|; X_j^dagger X_k=delta_jk |b,j><b,j|; X_j X_k=0",
        "same_site_anticommutator": "{X_j,X_j^dagger}=P_j, a rank-two projector, not the identity",
        "different_site_commutator": "[X_j,X_k^dagger]=0 for j!=k",
        "generated_field_algebra": "direct sum of ten M2(C) blocks, plus scalar identity on the 46-dimensional complement",
        "global_sum_on_20_dimensional_support": "sum_j X_j is one CAR mode with tenfold coordinate multiplicity",
        "ten_independent_CAR_modes_realized": False,
        "odd_fermion_or_boson_statistics_derived": False,
        "physical_electric_charge_or_spin_identified": False,
    }


def charge_continuity():
    """Actual correlation density/current, with local field charge identities."""
    protected = preparation.CORE[0]
    sites = tuple(j for j in range(12) if j != protected)
    weights = preparation.selectors(1)["completion_accepted_seam_weights"]
    checked = 0
    for state in STATES:
        density = {j: state[protected]*state[j] for j in sites}
        if sum(density.values()) != state[protected]:
            raise ArithmeticError("correlation density does not sum to protected occupation")
        for u, v in weights:
            for coin in (0, 1):
                output, (_, current) = native.native(state, (u, v), coin)
                for j in sites:
                    change = output[protected]*output[j]-density[j]
                    expected = state[protected]*current*(int(j == v)-int(j == u))
                    if change != expected:
                        raise ArithmeticError("source correlation current fails branchwise continuity")
                checked += 1
    for j, x, y in fields():
        for location in sites:
            difference = (STATES[x][protected]*STATES[x][location]
                          -STATES[y][protected]*STATES[y][location])
            if difference != int(location == j):
                raise ArithmeticError("collective field is not locally charged under correlation density")
    return {
        "protected_port": protected,
        "density": "Q_j=n_0 n_j, j!=0",
        "total_density_identity": "sum_{j!=0} Q_j=n_0 on binary total-two states",
        "current_on_recorded_seam": "I_uv=n_0*j_native",
        "branchwise_continuity": "Q'_u-Q_u=-I_uv, Q'_v-Q_v=I_uv; other Q_j unchanged",
        "local_charge_identity": "[Q_l,X_j]=delta_lj X_j for X_j=|0,j><5,j|",
        "source_branch_input_addresses_checked": checked,
        "charge_density_and_Ramsey_quadrature_joint_classical_values_assumed": False,
        "physical_electric_charge_identified": False,
    }


def incident_kernel(law, pair=(0, 5), *, selector_rounds=1, deadline=2):
    if type(deadline) is not int or not 1 <= deadline <= 8:
        raise ValueError("finite work-proposal deadline one through eight required")
    data = transfer(law, pair, selector_rounds)
    waiting = sum((data["remote"]**k for k in range(deadline)), s.zeros(10))
    return {"success": data["incident"]*waiting,
            "selector_timeout": data["timeout"]*waiting,
            "deadline_timeout": data["remote"]**deadline}


def _branch(rho, mapping, probability):
    out = defaultdict(F)
    for (x, y), value in rho.items():
        if x in mapping and y in mapping:
            out[mapping[x], mapping[y]] += probability*value
    return dict(out)


def detector(rho, pair=(0, 5)):
    """All twenty Ramsey/location outputs and the complete complementary bin."""
    outputs = []
    supported = F()
    for j, x, y in fields(pair):
        diagonal = rho.get((x, x), F())+rho.get((y, y), F())
        interference = rho.get((x, y), F())+rho.get((y, x), F())
        plus, minus = (diagonal+interference)/2, (diagonal-interference)/2
        if min(plus, minus) < 0:
            raise ArithmeticError("negative complete detector probability")
        outputs.extend(((j, "+", plus), (j, "-", minus)))
        supported += diagonal
    trace = sum((value for (x, y), value in rho.items() if x == y), F())
    complement = trace-supported
    if complement < 0:
        raise ArithmeticError("detector support exceeds the complete work state")
    outputs.append((None, "complement", complement))
    if sum((row[2] for row in outputs), F()) != trace:
        raise ArithmeticError("the detector omitted an output")
    return tuple(outputs)


def local_effects(region, pair=(0, 5), target_site=7):
    """Three local number-preserving Lüders effects in a declared access region."""
    if (type(region) is not tuple or len(set(region)) != len(region)
            or any(type(port) is not int or not 0 <= port < 12 for port in region)
            or not set((*pair, target_site)).issubset(region)):
        raise ValueError("a native access region containing the detector's three ports required")
    size = 2**len(region)
    local_address = lambda occupied: sum(1 << k for k, port in enumerate(region) if port in occupied)
    x, y = local_address((pair[0], target_site)), local_address((pair[1], target_site))
    plus = s.SparseMatrix(size, size, {(x, x): s.Rational(1, 2), (y, y): s.Rational(1, 2),
                                     (x, y): s.Rational(1, 2), (y, x): s.Rational(1, 2)})
    minus = s.SparseMatrix(size, size, {(x, x): s.Rational(1, 2), (y, y): s.Rational(1, 2),
                                      (x, y): -s.Rational(1, 2), (y, x): -s.Rational(1, 2)})
    return plus, minus, s.SparseMatrix(s.eye(size))-plus-minus


def embed_local(effect, region):
    """Tensor with identity outside the region, then restrict to actual N=2."""
    outside = tuple(port for port in range(12) if port not in region)
    indices = [sum(state[port] << k for k, port in enumerate(region)) for state in STATES]
    entries = {}
    for i, left in enumerate(STATES):
        for j, right in enumerate(STATES):
            if all(left[p] == right[p] for p in outside):
                value = effect[indices[i], indices[j]]
                if value:
                    entries[i, j] = value
    return s.ImmutableSparseMatrix(66, 66, entries)


def detector_consensus():
    """Exact agreement of two overlapping local instruments at the same event."""
    regions = (tuple(preparation.CAPS[0]), (0, 5, 7, 10))
    globals_ = []
    for region in regions:
        effects = local_effects(region)
        number = s.diag(*[mask.bit_count() for mask in range(2**len(region))])
        if any(effect*number != number*effect for effect in effects):
            raise ArithmeticError("local detector changes total native occupation")
        embedded = tuple(embed_local(effect, region) for effect in effects)
        if sum(embedded, s.zeros(66)) != s.eye(66):
            raise ArithmeticError("local detector lost the complementary outcome")
        for i, left in enumerate(embedded):
            for j, right in enumerate(embedded):
                if left*right != (left if i == j else s.zeros(66)):
                    raise ArithmeticError("same-event Lüders readout is not repeatable")
        globals_.append(embedded)
    if globals_[0] != globals_[1]:
        raise ArithmeticError("the two observer restrictions assign different shared effects")
    return {
        "detector_support": [0, 5, 7], "observer_regions": [list(region) for region in regions],
        "shared_region": sorted(set(regions[0]).intersection(regions[1])),
        "effects": "E_±=1/2 (|0,7>±|5,7>)(<0,7|±<5,7|); E_rest=I-E_+-E_-",
        "local_embeddings_identical_on_actual_source_exact": True,
        "number_preserving_local_instrument_exact": True,
        "state_update": "rho -> E_o rho E_o / Tr(E_o rho), when the selected outcome has nonzero probability",
        "repeatability": "E_v E_u=delta_uv E_u; two same-event reads give p(u,v)=delta_uv p(u)",
        "all_terminal_flags_and_rest_outcomes_retained": True,
        "later_native_evolution_of_work_state_allowed": True,
        "archive_agreement_requires_protected_classical_records": True,
        "same_quantity_unchanged_after_arbitrary_intervening_repairs": False,
        "second_region_and_coherent_apparatus_are_declared": True,
        "source_pulse_execution_or_physical_detector_calibration_claimed": False,
    }


def histories(law, pair=(0, 5), *, source_site=2, selector_rounds=1, deadline=2, phase_sign=1):
    """Complete declared work-record paths through the first cap incident tick.

    The existing selector subinstrument is marginalized to accepted seam or
    timeout at this boundary; its auxiliary coin transcripts are not expanded.

    Source preparation and Ramsey readout are declared coherent apparatus.
    The actual C/E Kraus law transports them; these routines do not compile
    another native pulse implementation or identify a physical time unit.
    """
    if type(deadline) is not int or not 1 <= deadline <= 3:
        raise ValueError("complete path enumeration uses deadline one through three")
    if type(phase_sign) is not int or phase_sign not in (-1, 1):
        raise ValueError("opposite real Ramsey input phases required")
    if source_site not in [j for j, _, _ in fields(pair)]:
        raise ValueError("a common-occupation source site required")
    a, b = pair
    x, y = address(a, source_site), address(b, source_site)
    initial = {(x, x): F(1, 2), (y, y): F(1, 2),
               (x, y): F(phase_sign, 2), (y, x): F(phase_sign, 2)}
    cap = set(preparation.CAPS[0])
    timeout = preparation.selectors(selector_rounds)["completion_timeout_probability"]
    live = [((), initial)]
    for attempt in range(1, deadline+1):
        next_live = []
        for prefix, rho in live:
            stopped = {key: timeout*value for key, value in rho.items()}
            yield (*prefix, ("selector_timeout",)), "selector_timeout", attempt, 0, detector(stopped, pair)
            for record, probability, mapping in events(law, selector_rounds):
                output = _branch(rho, mapping, probability)
                trace = (*prefix, record)
                if cap.intersection(record[:2]):
                    yield trace, "success", attempt, 1, detector(output, pair)
                elif attempt == deadline:
                    yield trace, "deadline_timeout", attempt, 0, detector(output, pair)
                else:
                    next_live.append((trace, output))
        live = next_live


def clock_experiment(law, pair=(0, 5), *, phase_sign=1):
    sums = defaultdict(F)
    contrast = defaultdict(F)
    shared_detector = defaultdict(F)
    count = zero = 0
    digest = hashlib.sha256()
    for record, flag, attempts, ticks, outputs in histories(law, pair, phase_sign=phase_sign):
        count += 1
        total = sum((row[2] for row in outputs), F())
        zero += total == 0
        sums[flag] += total
        for site, sign, probability in outputs:
            if site is not None:
                contrast[flag, site] += (1 if sign == "+" else -1)*probability
            shared_detector[flag, sign if site == 7 else "rest"] += probability
        digest.update(json.dumps([record, flag, attempts, ticks,
                                 [[site, sign, str(p)] for site, sign, p in outputs]],
                                separators=(",", ":")).encode())
    if sum(sums.values(), F()) != 1:
        raise ArithmeticError("the full stopped history lost probability")
    kernels = incident_kernel(law, pair)
    positions = [j for j, _, _ in fields(pair)]
    source = positions.index(2)
    for flag, kernel in kernels.items():
        for row, site in enumerate(positions):
            if contrast[flag, site] != phase_sign*F(kernel[row, source]):
                raise ArithmeticError("complete record/clock detector differs from the collective response")
    return {
        "source_site": 2, "input_relative_phase_sign": phase_sign,
        "clock_cap": list(preparation.CAPS[0]), "selector_rounds": 1, "work_proposal_deadline": 2,
        "clock": "first cap-incident work event; timeout flags have zero incident ticks",
        "enumerated_terminal_paths": count, "zero_probability_paths_retained": zero,
        "complete_record_and_detector_sha256": digest.hexdigest(),
        "record_boundary": "complete work coin/current-or-contrast, accepted seam, stopping flag and detector law; auxiliary selector transcripts marginalized to seam/timeout",
        "terminal_probabilities_exact": {key: str(value) for key, value in sorted(sums.items())},
        "terminal_signed_Ramsey_response_exact": {
            flag: {str(j): str(contrast[flag, j]) for j in positions} for flag in kernels},
        "transported_site_7_signed_success_response_exact": str(contrast["success", 7]),
        "shared_local_detector_joint_law_exact": {
            flag: {outcome: str(shared_detector[flag, outcome]) for outcome in ("+", "-", "rest")}
            for flag in kernels},
        "joint_law_interpretation": "each listed probability is on the equal-outcome diagonal for the two same-event local observers; all unequal pairs have probability zero",
        "all_21_detector_outcomes_per_path_retained": True,
        "source_selected_preparation_and_readout": False,
        "native_pulse_compilation_reexecuted": False,
        "physical_or_modular_time_calibrated": False,
    }


def spectral_report(law, pair=(0, 5)):
    data = transfer(law, pair)
    t = data["accepted"]
    q = data["survival"]
    eig = np.linalg.eigvalsh(np.array(t/q, float))
    if min(eig) <= 0 or max(eig) > 1+1e-13:
        raise ArithmeticError("collective conditional transfer lost its positive contraction")
    lower = min(t[i, i]-sum(abs(t[i, j]) for j in range(10) if j != i) for i in range(10))
    return {
        "accepted_transfer_exact": _matrix(t),
        "incident_transfer_exact": _matrix(data["incident"]),
        "remote_transfer_exact": _matrix(data["remote"]),
        "accepted_probability_exact": str(q),
        "selector_timeout_exact": str(data["timeout"]),
        "coherence_killing_per_proposal_exact": str(data["killing"]),
        "graph_formula": "T_accepted=(q-killing)I-L_common/2",
        "accepted_transfer_positive_floor_exact": str(lower),
        "conditional_transfer_eigenvalues_diagnostic": eig.tolist(),
        "conditional_log_decay_rates_diagnostic": [-math.log(float(x)) for x in eig if x < 1-1e-12],
        "generating_resolvent": "sum_n z^n T_accepted^n=(I-z T_accepted)^(-1), |z| ||T_accepted||<1",
        "local_return_measurement": "average the signed Ramsey return over ten separate coordinate preparations and local detectors, retaining every timeout",
        "averaged_return_resolvent": "trace((I-z*T_accepted)^(-1))/10=sum_r (1/10)/(1-z*lambda_r), counting multiplicity",
        "positive_return_spectral_weight_per_mode_exact": "1/10",
        "return_residue_in_z_at_inverse_eigenvalue": "-multiplicity/(10*lambda); every accepted transfer eigenvalue has nonzero return residue",
        "characteristic_polynomial_exact": [str(x) for x in t.charpoly().all_coeffs()],
        "modular_weight_ratios_by_common_site": {
            str(j): str(ratio) for (j, _, _), ratio in zip(data["fields"], data["modular_weight_ratios"], strict=True)},
        "phase_under_separately_implemented_modular_flow": "exp(i*s*log(weight_ratio)); s is the modular parameter",
        "native_repair_transfer_oscillation_phase": 0,
        "Lorentzian_dispersion_or_rest_mass_derived": False,
    }


def linear_color_attachment_test(law):
    """Exact degeneracy test for color on the ten native collective modes.

    If an unbroken continuous SU(3) action acts unitarily on this *linear*
    response space and commutes with its accepted transfer, each eigenspace
    must carry an SU(3) representation.  No nontrivial such representation
    fits in dimension one or two: a nonzero Lie-algebra homomorphism from
    the simple eight-dimensional su(3) into u(2), dimension four, would
    have to be injective.  This does not test composite or enlarged fields.
    """
    data = transfer(law)
    z = s.Symbol("z")
    polynomial = s.Poly(data["accepted"].charpoly(z).as_expr(), z)
    repeated = polynomial.gcd(polynomial.diff()).monic()
    double = data["survival"] - data["killing"]
    if (polynomial.degree() != 10 or repeated != s.Poly(z-double, z)
            or polynomial.eval(double) != 0 or polynomial.diff().eval(double) != 0
            or polynomial.diff().diff().eval(double) == 0):
        raise ArithmeticError("collective eigenspace multiplicity certificate changed")
    if data["accepted"] != data["accepted"].T:
        raise ArithmeticError("collective transfer is no longer diagonalizable by symmetry")
    return {
        "linear_response_dimension": 10,
        "one_double_eigenvalue_exact": str(double),
        "characteristic_polynomial_gcd_with_derivative_exact": str(repeated.as_expr()),
        "other_simple_eigenvalues": 8,
        "maximum_eigenspace_dimension": 2,
        "unbroken_commuting_continuous_SU3_can_act_as_a_triplet": False,
        "scope": "The ten X_j modes with the named accepted C/E transfer, assuming an unbroken unitary SU(3) action that commutes with it; no claim about composite, enlarged, gauge-dressed or alternate coherent fields.",
    }


def conditional_dispersion_test(law):
    """Test a stated sharp-pole dictionary, without assigning physical momentum.

    In the connected outer response block all modular frequencies coincide,
    whereas the graph Laplacian has a zero mode and nonzero modes. Identifying
    those frequencies with one field's energies and that Laplacian with one
    positive multiple of momentum squared cannot give a nonzero-speed free
    Klein-Gordon/Weyl dispersion. Interacting/thermal descriptions and other
    field, momentum or clock dictionaries are outside this necessary test.
    """
    data = transfer(law)
    outer = [i for i, (j, _, _) in enumerate(data["fields"]) if j not in preparation.CORE]
    laplacian = data["laplacian"].extract(outer, outer)
    if len(outer) != 9 or laplacian*s.ones(9, 1) != s.zeros(9, 1) or s.trace(laplacian) <= 0:
        raise ArithmeticError("outer graph lacks the required zero/nonzero mode distinction")
    if any(data["modular_weight_ratios"][i] != F(25, 8) for i in outer):
        raise ArithmeticError("outer collective modular frequency is not constant")
    # Exact positive weights and connectivity establish one zero eigenvalue,
    # with eight positive eigenvalues, without a floating rank decision.
    reached = {0}
    while True:
        more = reached | {j for i in reached for j in range(9) if i != j and laplacian[i, j] < 0}
        if more == reached:
            break
        reached = more
    if len(reached) != 9:
        raise ArithmeticError("the outer collective response graph is disconnected")
    return {
        "clock_parameters": "native proposal count or finite incident-event stopping clock; separate modular parameter s",
        "outer_response_modes": 9,
        "graph_zero_modes_exact": 1,
        "graph_positive_modes_exact": 8,
        "graph_laplacian_trace_exact": str(s.trace(laplacian)),
        "modular_frequency_for_every_outer_mode": "log(25/8)",
        "native_transfer_phases_for_every_outer_mode": "0",
        "tested_extra_dictionary": "one common calibrated modular frequency as E; positive common multiple of graph Laplacian eigenvalue as p^2; one sharp free field with nonzero propagation speed",
        "contradiction": "equal E^2 at ell=0 and ell>0 cannot satisfy E^2=m^2+c^2*ell with c!=0",
        "tested_sharp_free_dispersion_attachment_passes": False,
        "graph_eigenvalue_identified_as_continuum_momentum": False,
        "native_counter_identified_as_proper_or_Euclidean_time": False,
        "interacting_thermal_or_alternative_continuum_fields_excluded": False,
    }


def hierarchy_bounds(law):
    """Source-only rational bounds for specified transfer-log assignments.

    These are attenuation spectra in named proposal/stopping conventions;
    no mass identification, empirical ratio or free energy calibration enters.
    """
    data = transfer(law)
    q, killing = data["survival"], data["killing"]
    t, j, r = data["accepted"], data["incident"], data["remote"]
    floor = min(t[i, i]-sum(abs(t[i, k]) for k in range(10) if k != i) for i in range(10))
    j_floor = min(j.eigenvals())
    r_floor = min(r[i, i]-sum(abs(r[i, k]) for k in range(10) if k != i) for i in range(10))
    r_top = sum(r[0, k] for k in range(10))
    j_top = sum(j[0, k] for k in range(10))
    if (j_floor <= 0 or r_floor < 0 or j != j.T or r != r.T
            or any(sum(r[i, k] for k in range(10)) != r_top
                   or sum(j[i, k] for k in range(10)) != j_top for i in range(10))):
        raise ArithmeticError("finite incident-clock spectral certificate changed")
    # J(I+R) is similar to (I+R)^(1/2) J (I+R)^(1/2).
    clock_floor, clock_top = j_floor, j_top*(1+r_top)
    # Selection depends only on the declared cap/seam schedule, not on
    # the probe state or record law. Include the terminal deadline.
    weights = preparation.selectors(1)["completion_accepted_seam_weights"]
    cap = set(preparation.CAPS[0])
    incident_weight = sum((s.Rational(w) for edge, w in weights.items() if cap.intersection(edge)), s.Rational(0))
    remote_weight = q-incident_weight
    success_probability = incident_weight*(1+remote_weight)
    intervals = {
        "accepted_proposal": (floor, q-killing),
        "surviving_proposal": (floor/q, (q-killing)/q),
        "finite_first_incident_success": (clock_floor, clock_top),
        "conditional_finite_first_incident_success": (clock_floor/success_probability,
                                                     clock_top/success_probability),
    }
    rows = {}
    for convention, (lower, upper) in intervals.items():
        if not 0 < lower <= upper < 1:
            raise ArithmeticError("log-rate interval lacks strict finite positivity")
        # 1-u <= -log(lambda) <= (1-l)/l on [l,u].
        bound = s.factor((1-lower)/(lower*(1-upper)))
        rows[convention] = {"transfer_lower_exact": str(lower), "transfer_upper_exact": str(upper),
                            "positive_log_rate_spread_upper_exact": str(bound)}
    return {
        "exact_log_inequalities": "1-u <= -log(lambda) <= (1-l)/l for 0<l<=lambda<=u<1",
        "clock_spectral_argument": "J(I+R) similar to (I+R)^(1/2)J(I+R)^(1/2); J>=j_min I, R>=0",
        "conventions": rows,
        "field_coordinate_spread_or_copy_normalization_changes_poles": False,
        "outer_modular_frequency_count": 1,
        "full_collective_positive_modular_frequency_count": 2,
        "full_collective_modular_frequencies": ["log(8/5)", "log(25/8)"],
        "physical_mass_assignment_established": False,
    }


def fluctuation_certificate(copies=64):
    """Exact finite-copy moments on the declared independent source reference."""
    if type(copies) is not int or copies < 1:
        raise ValueError("a positive integer source-copy count required")
    reference = preparation.reference()
    plus, minus = reference[address(0, 2)], reference[address(5, 2)]
    delta = plus-minus
    if delta <= 0:
        raise ArithmeticError("source reference has no positive fluctuation polarization")
    occupancy = minus/delta
    variance_coefficient = (plus+minus-delta**2)/delta**2
    return {
        "independent_source_copies": copies,
        "reference_contract": "tensor-independent ideal sigma copies; joint-history independence is an additional dynamical premise",
        "active_weights_exact": [str(plus), str(minus)],
        "polarization_exact": str(delta),
        "fluctuation_field": "F_j=(M*delta)^(-1/2) sum_r X_j^(r), nine outer common sites",
        "mean_commutator_exact": "< [F_j,F_k^dagger] >= delta_jk",
        "off_site_commutators_exact": "[F_j,F_k^dagger]=0 for j!=k",
        "mean_occupation_exact": str(occupancy),
        "mean_creation_annihilation_reverse_exact": str(occupancy+1),
        "commutator_mean_square_error_exact": str(variance_coefficient/copies),
        "commutator_error_coefficient_over_M_exact": str(variance_coefficient),
        "normal_ordered_fourth_moment_exact": str(2*F(copies-1, copies)*occupancy**2),
        "thermal_Gaussian_fourth_moment_limit_exact": str(2*occupancy**2),
        "KMS_weight_ratio_exact": str(plus/minus),
        "modular_frequency": "log(25/8)",
        "ordered_Weyl_limit_under_IID_contract": "exp(-(2*n_B+1)*|sum v|^2/4-i*sum_{r<s}symplectic(v_r,v_s)/2)",
        "limit_justification": "same bounded IID ordered-Weyl product/Taylor argument as model/observer_time_fields/REFINEMENT.md, applied to this nonuniform native reference",
        "finite_operator_norm_CCR_or_Lorentzian_limit_claimed": False,
        "growing_copy_preparation_or_detector_executed": False,
        "cost_contract": "M separate reference preparations and records; their independence, joint control and preparation error need accounted services",
        "independent_copy_scaling_changes_two_point_transfer_poles": False,
        "two_point_evolution_contract": "identical separate-copy native evolution at one declared per-copy parameter; cross-copy means vanish, leaving the same T poles",
    }


def ordered_weyl(vectors, copies):
    """Exact single-copy matrix formula raised to M; floating diagnostic only."""
    cert = fluctuation_certificate(copies)
    wa, wb = map(lambda value: float(F(value)), cert["active_weights_exact"])
    delta = float(F(cert["polarization_exact"]))
    sx = np.array([[0., 1.], [1., 0.]], complex)
    sy = np.array([[0., -1j], [1j, 0.]], complex)
    product = np.eye(2, dtype=complex)
    for u, v in vectors:
        radius = math.hypot(float(u), float(v))
        angle = radius/math.sqrt(2*copies*delta)
        local = np.eye(2, dtype=complex) if radius == 0 else (
            math.cos(angle)*np.eye(2)+1j*math.sin(angle)*(u*sx+v*sy)/radius)
        product = product@local
    single = 1-wa-wb+wa*product[0, 0]+wb*product[1, 1]
    total = np.sum(np.array(vectors, float), axis=0)
    phase = sum(a*d-b*c for index, (a, b) in enumerate(vectors) for c, d in vectors[index+1:])
    n = float(F(cert["mean_occupation_exact"]))
    limit = np.exp(-(2*n+1)*np.dot(total, total)/4-1j*phase/2)
    return complex(single**copies), complex(limit)


def report():
    return {
        "schema": "oph.hierarchy_diagnosis.collective_response.v1",
        "source": "unchanged total-two C/E instruments, mixed-cap selector and protected-core reference",
        "modular_reference": "ideal prepared reference; finite preparation error is not propagated here",
        "protected_pair_control": {"algebra": algebra((0, 1)),
                                   "C_spectrum": spectral_report("C", (0, 1)),
                                   "E_spectrum": spectral_report("E", (0, 1))},
        "core_charged_collective_field": {"algebra": algebra((0, 5)),
                                         "C_spectrum": spectral_report("C"),
                                         "E_spectrum": spectral_report("E")},
        "collective_charge_continuity": charge_continuity(),
        "complete_incident_clock_experiments": {law: clock_experiment(law) for law in ("C", "E")},
        "overlapping_observer_detector": detector_consensus(),
        "conditional_linear_color_attachment": {law: linear_color_attachment_test(law)
                                                for law in ("C", "E")},
        "conditional_sharp_pole_attachment": {law: conditional_dispersion_test(law) for law in ("C", "E")},
        "conditional_hierarchy_bounds": {law: hierarchy_bounds(law) for law in ("C", "E")},
        "finite_copy_bosonic_fluctuations": fluctuation_certificate(),
        "positive_result": "a core-charged even collective coherence propagates in its common-occupation coordinate through complete native recorded histories",
        "native_raw_port_charge_support_moves": False,
        "collective_field_has_nonzero_cross_coordinate_response": True,
        "multiple_fermion_modes_or_3p1_Lorentzian_particles_identified": False,
        "empirical_mass_or_coupling_inputs": [],
    }
