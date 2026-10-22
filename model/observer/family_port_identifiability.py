"""Exact W5 port-information audit of the fixed two-event native recorder.

This asks whether its complete committed coin/current tape determines even
the norm of the icosahedral W5 component of an initial binary port load. It
does not supply a charged-family operator, energy generator, or mass map.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from model.common_source import slow_geometry
from model.observer import native


HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = HERE / "results" / "family_port_identifiability.json"
WITNESS_PORTS = ((4, 7), (4, 8))
PROBE_SEAM = (7, 8)


def antipodes() -> tuple[int, ...]:
    """Derive the six opposite port pairs in the checked native source chart."""
    frame = slow_geometry.raw_frame()
    arithmetic = slow_geometry.source()
    zero = arithmetic.q5(0)
    partners = []
    for i, point in enumerate(frame):
        matches = [j for j, other in enumerate(frame)
                   if all(arithmetic.q5_add(a, b) == zero
                          for a, b in zip(point, other, strict=True))]
        if len(matches) != 1 or matches[0] == i:
            raise ArithmeticError("native port chart has no unique antipode")
        partners.append(matches[0])
    result = tuple(partners)
    if any(result[result[i]] != i for i in range(12)):
        raise ArithmeticError("native port antipodes are not an involution")
    return result


def binary_load(ports: tuple[int, int]) -> tuple[int, ...]:
    if (len(ports) != 2 or ports[0] == ports[1]
            or any(type(p) is not int or not 0 <= p < 12 for p in ports)):
        raise ValueError("two distinct native ports are required")
    return tuple(int(i in ports) for i in range(12))


def tape_law(loads: tuple[int, ...]) -> dict[tuple[tuple[int, int], ...], Q]:
    """Complete committed two-event law, with zero-current outcomes retained."""
    law: dict[tuple[tuple[int, int], ...], Q] = defaultdict(Q)
    for first_weight, first in native.tick(native.State(loads)):
        for second_weight, last in native.tick(first):
            law[last.prefix] += first_weight * second_weight
    if sum(law.values(), Q(0)) != 1:
        raise ArithmeticError("two-event tape law is not normalized")
    return dict(law)


def w5_norm_squared(loads: tuple[int, ...], partner: tuple[int, ...]) -> Q:
    """||[(I+A)/2-11^T/12] loads||², without a floating chart or fit."""
    if len(loads) != 12 or sum(loads) != 2:
        raise ValueError("native total-two load required")
    mean = Q(sum(loads), 12)
    return sum((Q(loads[i] + loads[partner[i]], 2) - mean) ** 2
               for i in range(12))


def w5_projector(partner: tuple[int, ...]) -> tuple[tuple[Q, ...], ...]:
    """Exact port quadrupole projector in the native source ordering."""
    if len(partner) != 12 or any(partner[partner[i]] != i or partner[i] == i
                                 for i in range(12)):
        raise ValueError("twelve paired native ports required")
    return tuple(tuple(Q(int(i == j) + int(partner[i] == j), 2) - Q(1, 12)
                       for j in range(12)) for i in range(12))


def _matmul(left, right):
    return tuple(tuple(sum((a*b for a, b in zip(row, column, strict=True)), Q(0))
                       for column in zip(*right, strict=True)) for row in left)


def native_w5_response() -> dict:
    """Source-selected port response, plus its stationary selection boundary.

    Repeated fair-coin, uniform-seam attempts are a separately declared
    protocol. Its linear response is reconstructed from the actual native
    output of all twelve doubled-port probes; every total-two input then
    verifies the result. The fixed feedback recorder is a different protocol.
    """
    seams = tuple(slow_geometry.source().seams())
    if len(seams) != 30:
        raise ArithmeticError("native seam set changed")
    partner = antipodes()
    projector = w5_projector(partner)
    probes = []
    for port in range(12):
        loads = tuple(2*int(i == port) for i in range(12))
        outputs = [native.native(loads, seam, coin)[0]
                   for seam in seams for coin in (0, 1)]
        probes.append(tuple(Q(sum(output[i] for output in outputs), 120)
                            for i in range(12)))
    # Column j is the response to one unit initially at port j. The actual
    # doubled-port probes have two units, hence the denominator 120 above.
    transition = tuple(tuple(probes[j][i] for j in range(12)) for i in range(12))
    for loads in native.states():
        outputs = [native.native(loads, seam, coin)[0]
                   for seam in seams for coin in (0, 1)]
        for i in range(12):
            observed = Q(sum(output[i] for output in outputs), 60)
            expected = sum((transition[i][j]*loads[j] for j in range(12)), Q(0))
            if observed != expected:
                raise ArithmeticError("native port response is not linear on total two")
    source_adjacency = tuple(tuple(int(i != j and tuple(sorted((i, j))) in seams)
                                   for j in range(12)) for i in range(12))
    expected_transition = tuple(tuple(Q(int(i == j)) -
                                      Q(5*int(i == j)-source_adjacency[i][j], 60)
                                      for j in range(12)) for i in range(12))
    if transition != expected_transition:
        raise ArithmeticError("native response differs from icosahedral Laplacian")
    eigencheck = _matmul(transition, projector)
    if eigencheck != tuple(tuple(Q(9, 10)*value for value in row) for row in projector):
        raise ArithmeticError("W5 is not a 9/10 native port-response eigenspace")
    if _matmul(projector, projector) != projector or sum(projector[i][i] for i in range(12)) != 5:
        raise ArithmeticError("native W5 projector is not rank five")

    binary = native.states(binary=True)
    if len(binary) != 66:
        raise ArithmeticError("native binary total-two carrier changed")
    index = {loads: i for i, loads in enumerate(binary)}
    counts = [[0]*len(binary) for _ in binary]
    for i, loads in enumerate(binary):
        for seam in seams:
            for coin in (0, 1):
                counts[i][index[native.native(loads, seam, coin)[0]]] += 1
    if any(sum(row) != 60 for row in counts) or any(
        counts[i][j] != counts[j][i] for i in range(66) for j in range(66)
    ):
        raise ArithmeticError("native binary transfer is not symmetric stochastic")
    reached, frontier = {0}, [0]
    while frontier:
        i = frontier.pop()
        for j, count in enumerate(counts[i]):
            if count and j not in reached:
                reached.add(j)
                frontier.append(j)
    if len(reached) != 66:
        raise ArithmeticError("native binary transfer is disconnected")
    # Symmetric, stochastic, connected finite transfer has the unique uniform
    # stationary law. Its port covariance is computed from that exact law.
    mean = tuple(sum((Q(loads[i], 66) for loads in binary), Q(0)) for i in range(12))
    covariance = tuple(tuple(sum(((Q(loads[i])-mean[i])*(Q(loads[j])-mean[j])/66
                                  for loads in binary), Q(0)) for j in range(12))
                       for i in range(12))
    projected = _matmul(_matmul(projector, covariance), projector)
    if any(mean[i] != Q(1, 6) for i in range(12)) or projected != tuple(
        tuple(Q(5, 33)*value for value in row) for row in projector
    ):
        raise ArithmeticError("stationary W5 mean or covariance changed")
    return {
        "protocol": "fair-coin independent uniform choice among all 30 native seams",
        "linear_port_response": "T=I-L_icosahedron/60, reconstructed from 12 doubled-port native probes and verified on all 78 total-two loads",
        "W5_projector_rank": 5,
        "W5_one_attempt_response_eigenvalue": "9/10",
        "W5_n_attempt_response_eigenvalue": "(9/10)^n for n nonnegative uniform-seam native attempts",
        "stationary_binary_law": "uniform on 66 total-two binary configurations; unique by symmetry and connectedness",
        "stationary_W5_mean": "0",
        "stationary_W5_covariance": "(5/33)P5",
        "stationary_expected_W5_norm_squared": "25/33",
        "selection_boundary": (
            "Under the declared uniform-seam protocol, the native rule fixes a "
            "dimensionless W5 response rate and isotropic stationary fluctuations, "
            "not a nonzero family-shape direction, "
            "a charged-family action or a chiral mass coupling."
        ),
    }


def probe_law(loads: tuple[int, ...], seam: tuple[int, int]) -> dict[tuple[int, int], Q]:
    """One separately admitted, labelled native seam attempt."""
    law: dict[tuple[int, int], Q] = defaultdict(Q)
    for coin in (0, 1):
        _, record = native.native(loads, seam, coin)
        law[record] += Q(1, 2)
    return dict(law)


def all_seam_packet(loads: tuple[int, ...]) -> dict[tuple[int, int], tuple[tuple[int, int], ...]]:
    """Executor: ideal exact branch laws from fresh identical preparations."""
    return {seam: tuple(native.native(loads, seam, coin)[1] for coin in (0, 1))
            for seam in slow_geometry.source().seams()}


def recover_loads_from_packet(
    packet: dict[tuple[int, int], tuple[tuple[int, int], ...]], *, total: int = 2,
) -> tuple[int, ...]:
    """Estimator: the two exact coin-conditioned currents give each contrast.

    For oriented seam (a,b), j_0+j_1=x_a-x_b. Connected seam incidence
    determines all loads up to one constant; conserved total fixes it. The
    recovered loads must also reproduce every branch under the native rule.
    This packet is an ideal ensemble law, not a finite continuing tape.
    """
    seams = tuple(slow_geometry.source().seams())
    if type(total) is not int or total < 0 or set(packet) != set(seams):
        raise ValueError("complete labelled native seam packet and total required")
    neighbors: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for a, b in seams:
        records = packet[a, b]
        if (type(records) not in (tuple, list) or len(records) != 2
                or any(type(record) not in (tuple, list) or len(record) != 2
                       for record in records)
                or set(record[0] for record in records) != {0, 1}
                or any(type(record[0]) is not int or type(record[1]) is not int
                       for record in records)):
            raise ValueError("each seam needs both exact coin/current records")
        contrast = sum(record[1] for record in records)
        neighbors[a].append((b, -contrast))
        neighbors[b].append((a, contrast))
    relative = {0: 0}
    frontier = [0]
    while frontier:
        a = frontier.pop()
        for b, increment in neighbors[a]:
            candidate = relative[a] + increment
            if b in relative:
                if relative[b] != candidate:
                    raise ValueError("inconsistent native edge contrasts")
            else:
                relative[b] = candidate
                frontier.append(b)
    if len(relative) != 12:
        raise ValueError("native seam packet does not connect all ports")
    shift = Q(total - sum(relative.values()), 12)
    recovered = tuple(relative[i] + shift for i in range(12))
    if any(value.denominator != 1 or value < 0 for value in recovered):
        raise ValueError("native packet has no nonnegative integer load solution")
    loads = tuple(int(value) for value in recovered)
    for seam in seams:
        expected = {coin: native.native(loads, seam, coin)[1]
                    for coin in (0, 1)}
        observed = {coin: (coin, current) for coin, current in packet[seam]}
        if observed != expected:
            raise ValueError("coin-conditioned record violates native repair")
    return loads


def _encoded_law(law: dict[tuple[tuple[int, int], ...], Q]) -> list[dict]:
    return [{"records": [list(r) for r in key], "probability": str(weight)}
            for key, weight in sorted(law.items())]


def _encoded_probe(law: dict[tuple[int, int], Q]) -> list[dict]:
    return [{"coin_current": list(key), "probability": str(weight)}
            for key, weight in sorted(law.items())]


def run() -> dict:
    partner = antipodes()
    response = native_w5_response()
    if len(native.states(binary=True)) != 66:
        raise ArithmeticError("binary source carrier changed")
    classes: dict[tuple, list[tuple[int, ...]]] = defaultdict(list)
    for loads in native.states(binary=True):
        law = tape_law(loads)
        classes[tuple(sorted(law.items()))].append(loads)
    class_sizes = sorted((len(group) for group in classes.values()), reverse=True)

    left, right = (binary_load(ports) for ports in WITNESS_PORTS)
    law_left, law_right = tape_law(left), tape_law(right)
    norm_left = w5_norm_squared(left, partner)
    norm_right = w5_norm_squared(right, partner)
    probe_left = probe_law(left, PROBE_SEAM)
    probe_right = probe_law(right, PROBE_SEAM)
    reconstructed = 0
    for loads in native.states():
        if recover_loads_from_packet(all_seam_packet(loads)) != loads:
            raise AssertionError("all-seam source readout did not reconstruct loads")
        reconstructed += 1
    checks = {
        "antipodes_derived_from_checked_native_port_chart":
            len({tuple(sorted((i, j))) for i, j in enumerate(partner)}) == 6,
        "all_binary_preparations_partition_by_complete_tape_law":
            sum(class_sizes) == 66 and class_sizes == [46, 10, 9, 1],
        "witness_has_identical_complete_native_tape_law":
            law_left == law_right and len(law_left) == 4,
        "witness_has_distinct_exact_W5_norms":
            (norm_left, norm_right) == (Q(5, 3), Q(2, 3)),
        "labelled_native_seam_probe_separates_witness":
            probe_left != probe_right and PROBE_SEAM in slow_geometry.source().seams(),
        "fresh_prepared_all_seam_laws_recover_all_78_native_loads": reconstructed == 78,
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return {
        "schema": "oph.native_family_port_identifiability.v1",
        "status": "two_event_tape_does_not_determine_W5_port_norm",
        "all_pass": True,
        "checks": checks,
        "source": {
            "native_update": "model/observer/native.py: native, tick, selected_seam",
            "port_chart": "model/common_source/slow_geometry.py: raw_frame",
            "binary_preparations": 66,
            "record_law_class_sizes": class_sizes,
            "antipodal_pairs": [list((i, j)) for i, j in enumerate(partner) if i < j],
            "mass_reference_data_consumed": False,
        },
        "port_resolved_readout": {
            "protocol": "ideal exact fair-coin branch law at each of 30 labelled seams, each on fresh identical preparations",
            "edge_contrast_identity": "j_coin0+j_coin1=x_a-x_b",
            "branch_validation": "both coin-conditioned currents replayed through native repair after load recovery",
            "source_loads_reconstructed_exactly": reconstructed,
            "readout": "complete load vector and hence its W5 projector",
            "premises": ["identical repeatable preparation", "labelled all-seam access",
                         "ideal exact ensemble branch laws, not finite tape samples",
                         "known conserved total two"],
        },
        "uniform_seam_W5_response": response,
        "witness": {
            "initial_occupied_ports": [list(WITNESS_PORTS[0]), list(WITNESS_PORTS[1])],
            "complete_coin_current_tape_law": _encoded_law(law_left),
            "W5_projector": "(I+A_antipode)/2-11^T/12",
            "W5_norm_squared": [str(norm_left), str(norm_right)],
            "separating_native_probe_seam": list(PROBE_SEAM),
            "probe_laws": [_encoded_probe(probe_left), _encoded_probe(probe_right)],
        },
        "scope": (
            "The fixed two-event controller reads only its committed coin/current tape. "
            "Its record law does not determine even the initial W5 port norm. "
            "Fresh, identically prepared all-seam records recover the load W5 projection exactly. "
            "The charged defect W5 record and any physical family/mass interpretation "
            "of this load observable remain separate attachments."
        ),
        "remaining_interface": (
            "A candidate charged-family route must emit a refinement-stable port-resolved "
            "W5 record, select a family action and source coupling from the same native "
            "dynamics, attach a chiral mass-response operator and its observable spectrum, "
            "then provide a clock/energy and physical mass-scheme map."
        ),
    }


def canonical(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", type=Path, default=DEFAULT_RECEIPT)
    group.add_argument("--write", type=Path)
    args = parser.parse_args()
    encoded = canonical(run())
    if args.write is not None:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    elif args.check.read_text(encoding="utf-8") != encoded:
        raise SystemExit("Family-port identifiability receipt mismatch")
    print(json.dumps({"all_pass": True, "sha256": hashlib.sha256(encoded.encode()).hexdigest()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
