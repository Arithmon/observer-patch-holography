"""An equivariant moving sharp-q=1 memory with optimal boundary resets.

The source is the existing binary total-two carrier. A complete native raw
flag determines a transported classical core. Immediate coherent recovery
encodes the unknown qubit at that core, unless it is antipodal; then the
attempted seam becomes the new core. The quantum recovery is a separately
supplied service, implemented by mobile_service. No particle or energy is
identified by this controlled motion.

Universality of the core walk does not assume our equal-amplitude code.
The prior all-seam flagged-code theorem forces nonedge support and scalar
compressions V* n_p V for every exact qubit. At sharp q=1 only the safe
exterior ports survive. For each such j, the two rows at {a,j},{b,j} form
a square B_j with B_j* B_j=w_j I_2. Therefore B_j B_j*=w_j I_2, also for
w_j=0. Each core row has squared norm w_j and sum_j w_j=1. Scalarity of
V* n_a V then forces its trace to be one, hence V* n_a V=I_2/2, likewise b.
Every external neighbor move of either core endpoint thus has probability
1/(30*2*2)=1/120 under uniform seam selection and the fair native coin.

This proves the same walk for arbitrary exact codes, unknown logical states,
untouched references, and history-selected code choices. The sole impossibility
boundary is the six antipodal cores, where there is no safe exterior port.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache

from . import screening
from ..common_source import slow_geometry
from ..exchange_sector import sector
from ..observer import record_recovery as records
from ..observer.source_discrimination import circuit


SEAMS = tuple(records.seams())
CORE = (0, 1)
CLASSES = ("edge", "distance_two", "antipodal")
KILLED = ((F(29, 30), F(1, 30)), (F(1, 30), F(19, 20)))
CONTINUING = ((F(29, 30), F(1, 30)), (F(1, 20), F(19, 20)))
POTENTIAL = {"edge": 150, "distance_two": 120, "antipodal": 0}


def core_class(center):
    center = sector.core(center)
    return {4: "edge", 2: "distance_two", 0: "antipodal"}[
        len(screening.safe_exterior(center))]


def legal_cores():
    return tuple(c for c in screening.CORES if core_class(c) != "antipodal")


def _legal(center):
    center = sector.core(center)
    if core_class(center) == "antipodal":
        raise ValueError("an edge or distance-two core supporting an exact qubit required")
    return center


@lru_cache(maxsize=1)
def rotations():
    result = tuple(sorted(tuple(map(int, p)) for p in slow_geometry.source().rotations()))
    if len(result) != 60 or len(set(result)) != 60:
        raise ArithmeticError("native proper rotation group changed")
    return result


def prototype(center):
    return CORE if core_class(_legal(center)) == "edge" else (0, 2)


def chart(center):
    center = _legal(center)
    base = prototype(center)
    return next(p for p in rotations() if tuple(sorted(p[i] for i in base)) == center)


def _rotate_supports(columns, p):
    return tuple(tuple(sorted(tuple(sorted(p[i] for i in pair)) for pair in row))
                 for row in columns)


def supports(center):
    """Actual port rotations of the two existing sharp-code prototypes.

    Their two-element stabilizers fix each logical column. Consequently the
    result is independent of the chosen chart and has no hidden frame register.
    """
    center = _legal(center)
    return _rotate_supports(screening.sharp_code(prototype(center)), chart(center))


def _label(center, law, seam, coin, raw):
    center, law = _legal(center), records.candidate(law)
    if type(seam) is not tuple or seam not in SEAMS or any(type(p) is not int for p in seam):
        raise ValueError("a sorted seam from the complete native menu required")
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("an integer native binary coin required")
    if type(raw) is not int or not 0 <= raw < (1 << circuit.RECORD_BITS[law]):
        raise ValueError("a raw integer slot from the complete C/E menu required")
    return center, law, seam, coin, raw


def branch_images(center, law, seam, coin, raw):
    """Two disjoint source-image address lists, including zero/abort slots."""
    return _images(*_label(center, law, seam, coin, raw))


@lru_cache(maxsize=None)
def _images(center, law, seam, coin, raw):
    rows = []
    label = circuit.record_values(law, coin).get(raw)
    for column in supports(center):
        row = []
        for pair in column:
            state = tuple(int(p in pair) for p in range(12))
            output, record = records.recorded_step(law, state, seam, coin)
            if label == record[-1]:
                row.append(circuit.encode(output))
        rows.append(tuple(row))
    if len(rows[0]) != len(rows[1]) or len(set(rows[0]+rows[1])) != sum(map(len, rows)):
        raise ArithmeticError("moving code lost scalar injective native branch")
    return tuple(rows)


def branch_weight(center, law, seam, coin, raw):
    """Raw branch probability includes the fair coin, but not seam selection."""
    rows = branch_images(center, law, seam, coin, raw)
    return F(len(rows[0]), 2*len(supports(center)[0]))


def transition(center, law, seam, coin, raw):
    """Complete record-only controller, including off-code zero and abort flags.

    Invalid raw slots halt. Valid zero-code-weight slots select identity and
    retain their flag; their zero ideal probability is not postselection.
    On a positive flag the native core uses only its decoded signed current.
    The sole reset changes an antipodal core to the attempted native seam.
    """
    center, law, seam, coin, raw = _label(center, law, seam, coin, raw)
    admitted = raw in circuit.record_values(law, coin)
    weight = branch_weight(center, law, seam, coin, raw)
    moved = target = None
    reset = False
    if admitted:
        moved = center
        if weight:
            value = circuit.record_values(law, coin)[raw]
            moved = sector.moved_core(center, law, (*seam, coin, value))
        reset = core_class(moved) == "antipodal"
        target = seam if reset else moved
    return {"core": center, "law": law, "seam": seam, "coin": coin, "raw": raw,
            "admitted": admitted, "positive": bool(weight), "weight": weight,
            "native_core": moved, "target_core": target, "reset_performed": reset,
            "halt": not admitted, "identity_zero_branch": admitted and not weight}


def branches(center, law, seam):
    _label(center, law, seam, 0, 0)
    return tuple(transition(center, law, seam, coin, raw)
                 for coin in (0, 1) for raw in range(1 << circuit.RECORD_BITS[law]))


def replay(history, *, initial=CORE):
    """Replay a complete selected-law history under uniform seam selection.

    A law choice is conditioned upon, not assigned a hidden random weight.
    Every raw flag is retained. An abort halts; any supplied suffix after a
    halt is rejected. Valid zero-probability histories remain probability zero.
    """
    current = _legal(initial)
    if type(history) not in (tuple, list):
        raise ValueError("a complete history sequence required")
    weight, resets, archive, halted = F(1), 0, [], False
    for event in history:
        if halted:
            raise ValueError("a halted raw history cannot perform another repair")
        if type(event) not in (tuple, list) or len(event) != 4:
            raise ValueError("each event must be (law, seam, coin, raw)")
        row = transition(current, *event)
        archive.append((row["law"], row["seam"], row["coin"], row["raw"]))
        weight *= row["weight"]/30
        resets += row["reset_performed"]
        halted = row["halt"]
        if not halted:
            current = row["target_core"]
    return {"initial_core": initial, "terminal_core": current, "probability": weight,
            "reset_count": resets, "halted": halted, "archive": tuple(archive)}


def native_kernel(center):
    """Universal uncorrected core kernel for ANY exact sharp-q=1 qubit."""
    center = _legal(center)
    result = Counter({center: F(1)})
    for seam in SEAMS:
        if len(set(center) & set(seam)) == 1:
            a, b = seam
            target = tuple(sorted(b if p == a else a if p == b else p for p in center))
            result[target] += F(1, 120)
            result[center] -= F(1, 120)
    return dict(sorted(result.items()))


def continuing_kernel(center):
    """All60-core transition row, with the specified boundary-only reset."""
    center = _legal(center)
    result = Counter({center: F(1)})
    for seam in SEAMS:
        if len(set(center) & set(seam)) == 1:
            a, b = seam
            moved = tuple(sorted(b if p == a else a if p == b else p for p in center))
            target = seam if core_class(moved) == "antipodal" else moved
            result[target] += F(1, 120)
            result[center] -= F(1, 120)
    return dict(sorted(result.items()))


def _attempts(value):
    if type(value) is not int or value < 0:
        raise ValueError("a nonnegative integer attempt count required")
    return value


def survival(attempts, *, initial=CORE):
    """Probability that the natural transported core has not become impossible."""
    attempts, initial = _attempts(attempts), _legal(initial)
    row = (F(1), F()) if core_class(initial) == "edge" else (F(), F(1))
    for _ in range(attempts):
        row = tuple(sum(row[i]*KILLED[i][j] for i in range(2)) for j in range(2))
    return sum(row)


def expected_resets(attempts, *, initial=CORE):
    attempts, initial = _attempts(attempts), _legal(initial)
    start_d = F(core_class(initial) == "distance_two")
    return F(attempts, 150)+(start_d-F(2, 5))*F(1, 5)*(1-F(11, 12)**attempts)


def reset_lower_bound(attempts, *, initial=CORE):
    """Universal bound, allowing earlier, adaptive and randomized resets."""
    attempts, initial = _attempts(attempts), _legal(initial)
    return max(F(), F(attempts-POTENTIAL[core_class(initial)], 150))


def finite_horizon_certificate():
    """Exact Bellman proof for every fixed finite horizon, not a sampled optimum.

    Let v(n)=n*a+(1-r**n)*b with a=(1/150,1/150),
    b=(-2/25,3/25), r=11/12. The checked coefficient identities establish
    v(0)=0 and v(n+1)=g+P*v(n) for the boundary-reset policy. Furthermore
    v_D(n)-v_E(n)=(1-r**n)/5 lies in [0,1/5), so a legal-class reset costing
    one is strictly worse than retaining the native core. At an antipodal
    core a reset is mandatory and the edge class minimizes continuation.
    Induction therefore proves Bellman optimality, even for record-adaptive
    randomized policies and arbitrary exact code choices. Initial voluntary
    resets and repeated resets at one boundary cannot improve the optimum.
    """
    a, b, r = (F(1, 150),)*2, (F(-2, 25), F(3, 25)), F(11, 12)
    g = (F(), F(1, 60))
    product = lambda v: tuple(sum(CONTINUING[i][j]*v[j] for j in range(2))
                             for i in range(2))
    if product(a) != a or product(b) != tuple(r*x for x in b):
        raise ArithmeticError("finite-horizon coefficient eigen-identities changed")
    if tuple(a[i]+(1-r)*b[i] for i in range(2)) != g:
        raise ArithmeticError("finite-horizon immediate-reset coefficients changed")
    if b[1]-b[0] != F(1, 5) or not 0 < r < 1:
        raise ArithmeticError("Bellman voluntary-reset margin changed")
    return {"fixed_horizon_global_minimum": True,
            "minimum_expected_resets_from_edge": "n/150-(2/25)*(1-(11/12)^n)",
            "minimum_expected_resets_from_distance_two": "n/150+(3/25)*(1-(11/12)^n)",
            "continuation_value_difference": "(1-(11/12)^n)/5",
            "voluntary_reset_cost": 1,
            "maximal_future_value_improvement_strict_upper": "1/5",
            "coefficient_a": [str(x) for x in a], "coefficient_b": [str(x) for x in b],
            "ratio": str(r), "immediate_reset_vector": [str(x) for x in g],
            "verified_identities": ["P*a=a", "P*b=(11/12)*b", "a+(1/12)*b=g"],
            "horizon_scope": "A nonnegative finite number of attempted repairs fixed before the random history; every branch must restore a legal exact sharp-q=1 code; no postselection or halting to evade attempts.",
            "optimizes": "number of classical core-label interventions, not source pulses, physical work, or record storage"}


def geometry_certificate():
    checks, stabilizers = 0, {}
    for base in (CORE, (0, 2)):
        proto = _rotate_supports(screening.sharp_code(base), tuple(range(12)))
        fixed = tuple(p for p in rotations() if tuple(sorted(p[i] for i in base)) == base)
        if len(fixed) != 2 or any(_rotate_supports(proto, p) != proto for p in fixed):
            raise ArithmeticError("prototype stabilizer changes logical columns")
        stabilizers[str(base)] = len(fixed)
    for center in legal_cores():
        columns = supports(center)
        for p in rotations():
            target = tuple(sorted(p[i] for i in center))
            if _rotate_supports(columns, p) != supports(target):
                raise ArithmeticError("moving code family is not equivariant")
            checks += 1
        for endpoint in center:
            if any(F(sum(endpoint in pair for pair in column), len(column)) != F(1, 2)
                   for column in columns):
                raise ArithmeticError("prototype core occupation is not one half")
    return {"legal_cores": 60, "impossible_antipodal_cores": 6,
            "proper_rotations": 60, "equivariance_checks": checks,
            "prototype_stabilizers": stabilizers,
            "logical_frame_register_required": False,
            "core_endpoint_occupation_compression": "I2/2",
            "universal_occupation_proof": "For each safe exterior j, square B_j obeys B_j*B_j=w_j I2, hence B_j B_j*=w_j I2; sum_j w_j=1 and scalar V*n_aV imply V*n_aV=I2/2."}


def exact_certificate():
    """All raw-slot census, native/reset kernels, stationarity and optimality."""
    geometry = geometry_certificate()
    census = Counter()
    for center in legal_cores():
        natural, controlled = native_kernel(center), continuing_kernel(center)
        for law in ("C", "E"):
            raw_natural, raw_controlled = Counter(), Counter()
            for seam in SEAMS:
                rows = branches(center, law, seam)
                if sum(row["weight"] for row in rows) != 1:
                    raise ArithmeticError("complete source flag weights do not sum to one")
                for row in rows:
                    census[(law, "raw")] += 1
                    census[(law, "positive")] += row["positive"]
                    census[(law, "abort")] += not row["admitted"]
                    if row["positive"]:
                        raw_natural[row["native_core"]] += row["weight"]/30
                        raw_controlled[row["target_core"]] += row["weight"]/30
                        # Both logical branch images really retain q=1 in the
                        # transported core before the conditional source service.
                        for image in branch_images(center, law, seam, row["coin"], row["raw"]):
                            if any(sector.charge(circuit.decode(i), row["native_core"]) != 1 for i in image):
                                raise ArithmeticError("native branch lost sharp relational charge")
            if dict(raw_natural) != natural or dict(raw_controlled) != controlled:
                raise ArithmeticError("raw instrument differs from universal core walk")
        kinds = Counter()
        for target, weight in natural.items():
            kinds[core_class(target)] += weight
        expected = ((F(29, 30), F(1, 30), F()) if core_class(center) == "edge"
                    else (F(1, 30), F(19, 20), F(1, 60)))
        if tuple(kinds[k] for k in CLASSES) != expected:
            raise ArithmeticError("native class-chain row changed")
        if sum(w*POTENTIAL[core_class(c)] for c, w in natural.items()) != POTENTIAL[core_class(center)]-1:
            raise ArithmeticError("first-hit potential drift changed")
    legal = legal_cores()
    stationary = {c: F(1, 50) if core_class(c) == "edge" else F(1, 75) for c in legal}
    incoming = Counter()
    kernels = {c: continuing_kernel(c) for c in legal}
    for center in legal:
        for target, weight in kernels[center].items():
            incoming[target] += stationary[center]*weight
    if sum(stationary.values()) != 1 or dict(incoming) != stationary:
        raise ArithmeticError("full moving-core stationary distribution changed")
    one_way = [(c, d) for c in legal for d, w in kernels[c].items()
               if w and not kernels[d].get(c, F())]
    if len(one_way) != 60:
        raise ArithmeticError("reset-driven directed links changed")
    witness = one_way[0]
    return {"geometry": geometry,
            "raw_census": {law: {kind: census[law, kind] for kind in ("raw", "positive", "abort")}
                           for law in ("C", "E")},
            "uniform_seam_probability": "1/30", "native_coin_probability": "1/2",
            "per_external_core_neighbor_move": "1/120",
            "killed_class_order": ["edge", "distance_two"],
            "killed_class_kernel": [[str(v) for v in row] for row in KILLED],
            "first_impossible_core_mean_attempts": {"edge": 150, "distance_two": 120},
            "survival_eigenvalues": ["(115+sqrt(17))/120", "(115-sqrt(17))/120"],
            "survival_from_edge": {str(n): str(survival(n)) for n in (0, 1, 2, 10, 100)},
            "reset_policy": "Preserve the recorded native core unless antipodal; then use the attempted seam.",
            "continuing_class_kernel": [[str(v) for v in row] for row in CONTINUING],
            "stationary_class_probabilities": ["3/5", "2/5"],
            "stationary_probability_per_core": {"edge": "1/50", "distance_two": "1/75"},
            "stationary_reset_probability_per_attempt": "1/150",
            "expected_resets_start_edge": "(n-12+12*(11/12)^n)/150",
            "universal_reset_lower_bound": "E[N_n] >= max(0,(n-h(initial_core))/150)",
            "finite_horizon_optimality": finite_horizon_certificate(),
            "optimality_domain": "Every complete uniform-seam/fair-coin attempt starts in an exact all-seam sharp-q=1 qubit; native current transports the core; arbitrary adaptive/randomized legal-core resets are counted, including preemptive resets; immediate accurate quantum recovery is supplied.",
            "optimality_proof": "h(edge)=150,h(distance_two)=120,h(antipodal)=0 drops by one in conditional expectation per native attempt. Every counted reset raises h by at most150. Telescope to obtain the universal bound. Boundary-only antipodal-to-edge resets raise h by150 and attain asymptotic rate1/150.",
            "full60_detailed_balance": False, "one_way_links": len(one_way),
            "one_way_witness": {"from": list(witness[0]), "to": list(witness[1]),
                                "forward": str(kernels[witness[0]][witness[1]]), "reverse": "0"},
            "record_account": {"core_controller_bits": 7,
                               "per_attempt_fixed_law_complete_flag_bits": {"C": 8, "E": 9},
                               "tagged_law_choice_extra_bits": 1,
                               "reset_decision_and_target_replayable_from_complete_tape": True,
                               "independent_reset_archive_required": False,
                               "growing_archive_retained": True},
            "physical_scope": {"autonomous_particle": False, "localized_screened_defect": False,
                               "coherent_spatial_position": False, "energy_or_mass_derived": False,
                               "physical_duration_calibrated": False,
                               "source_pulse_cost_is_reset_count": False}}
