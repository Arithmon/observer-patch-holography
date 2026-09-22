"""Independent labeled-tensor replay; imports no producer mathematics."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ZERO, ONE = (Q(0), Q(0)), (Q(1), Q(0))


def plus(x, y):
    return x[0]+y[0], x[1]+y[1]


def times(x, y):
    return x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]


def conjugate(x):
    return x[0], -x[1]


def total(xs):
    out = ZERO
    for x in xs:
        out = plus(out, x)
    return out


def magnitude(x):
    return x[0]**2+x[1]**2


def encoded(x):
    return [str(x[0]), str(x[1])]


def check(condition, reason):
    if not condition:
        raise ValueError(reason)


def equal(actual, expected, reason):
    # JSON comparison distinguishes booleans from 0/1 and canonical fractions.
    check(json.dumps(actual, sort_keys=True) == json.dumps(expected, sort_keys=True), reason)


def no_duplicate_pairs(pairs):
    out = {}
    for key, value in pairs:
        check(key not in out, "duplicate JSON key")
        out[key] = value
    return out


def load(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=no_duplicate_pairs)


def mm(a, b):
    return [[total(times(a[i][k], b[k][j]) for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def tensor(a):
    return [[times(a[i//2][j//2], a[i % 2][j % 2]) for j in range(4)]
            for i in range(4)]


def tensor_result(u, sign):
    # Labeled Hilbert space: LL, LR, RL, RR. Avoid sqrt(2) by dividing
    # expectations by norm 2 for LR +/- RL, and norm 1 for labeled LR.
    v = [ZERO, ONE, (Q(sign), Q(0)), ZERO]
    w = [total(times(row[j], v[j]) for j in range(4)) for row in tensor(u)]
    n = Q(1+sign*sign)
    check(sum(magnitude(x) for x in w) == n, "tensor state norm")
    p = {"20": str(magnitude(w[0])/n),
         "11": str((magnitude(w[1])+magnitude(w[2]))/n),
         "02": str(magnitude(w[3])/n)}
    rho = []
    for a in range(2):
        row = []
        for b in range(2):
            z = total(plus(times(conjugate(w[2*a+j]), w[2*b+j]),
                           times(conjugate(w[2*j+a]), w[2*j+b])) for j in range(2))
            row.append(encoded((z[0]/n, z[1]/n)))
        rho.append(row)
    return {"probabilities": p, "one_body_density": rho, "norm_squared": "1"}


def hidden_label_result(u, sign, eta):
    """Partial-traced labeled density, with offdiagonal sign*eta/2.

    eta is the squared overlap of two normalized unobserved internal states;
    no irrational choice of their amplitudes enters the reduced density.
    """
    t = tensor(u)
    density = []
    for i in range(4):
        row = []
        for j in range(4):
            direct = plus(times(t[i][1], conjugate(t[j][1])),
                          times(t[i][2], conjugate(t[j][2])))
            cross = plus(times(t[i][1], conjugate(t[j][2])),
                         times(t[i][2], conjugate(t[j][1])))
            row.append(times((Q(1, 2), Q(0)), plus(direct, times((sign*eta, Q(0)), cross))))
        density.append(row)
    check(total(density[i][i] for i in range(4)) == ONE, "hidden-label density trace")
    p = {"20": str(density[0][0][0]),
         "11": str(density[1][1][0]+density[2][2][0]),
         "02": str(density[3][3][0])}
    rho = [[encoded(total(plus(density[2*b+j][2*a+j], density[2*j+b][2*j+a])
                          for j in range(2))) for b in range(2)] for a in range(2)]
    return {"probabilities": p, "one_body_density": rho, "norm_squared": "1"}


def verify(receipt):
    equal(set_to_list(receipt), sorted([
        "schema", "source_pins", "binding", "splitter", "preparation", "results",
        "one_particle_bilinears", "spin_rotation", "gap", "error_contract", "scope",
        "partial_distinguishability",
    ]), "receipt fields")
    equal(receipt["schema"], "oph.spin_exchange.v1", "schema")
    paths = ["code/sm_fermion_current/quantum_link_receipt.json",
             "code/sm_fermion_current/quantum_link.py",
             "code/sm_fermion_current/current.py",
             "paper/tex_fragments/FERMION_SOURCE_CURRENT.tex",
             "code/spin_exchange/build.py", "code/spin_exchange/verify.py",
             "code/spin_exchange/test_spin_exchange.py"]
    equal(receipt["source_pins"], {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                                   for p in paths}, "parent source pins")
    parent = load(ROOT/paths[0])
    rows = [x for x in parent["channels"] if x["multiplet"] == "e_c"]
    check(len(rows) == 1, "unique e_c channel")
    row = rows[0]
    equal(receipt["binding"], {
        "parent_channel_id": row["id"], "multiplet": "e_c", "step": row["step"],
        "kappa": row["kappa"],
        "reuse": "Cayley coefficients only; no parent quantum-link Gauss or state transfer"
    }, "parent binding")
    u = [[tuple(Q(x) for x in z) for z in line] for line in receipt["splitter"]]
    check(len(u) == 2 and all(len(line) == 2 for line in u) and
          all(len(z) == 2 for line in u for z in line), "matrix dimensions")
    x, y = tuple(map(Q, row["left_amplitude"])), tuple(map(Q, row["right_amplitude"]))
    equal(receipt["splitter"], [[encoded(x), encoded(y)], [encoded(y), encoded(x)]],
          "parent hopping amplitudes")
    identity = [[ONE, ZERO], [ZERO, ONE]]
    adjoint = [[conjugate(u[j][i]) for j in range(2)] for i in range(2)]
    check(mm(adjoint, u) == identity, "one-particle unitary")
    half = Q(row["step"])*Q(row["kappa"])/2
    lhs = [[ONE, (Q(0), half)], [(Q(0), half), ONE]]
    rhs = [[ONE, (Q(0), -half)], [(Q(0), -half), ONE]]
    check(mm(lhs, u) == rhs, "Cayley defining equation")
    equal(receipt["preparation"], {
        "spatial_occupation": [1, 1], "identical_internal_and_spin_channel": True,
        "spin_frame": "same declared channel in local port frames; raw spin vectors may differ",
        "two_particles_in_same_channel": True,
        "distinguishable_control": "orthogonal unobserved internal labels",
        "normalization": "normalized symmetric/antisymmetric LR tensor; labeled LR for control",
    }, "same-channel preparation contract")
    results = {"fermion": tensor_result(u, -1), "boson": tensor_result(u, 1),
               "distinguishable": tensor_result(u, 0)}
    equal(receipt["results"], results, "labeled tensor probabilities and one-body density")
    for result in results.values():
        check(sum(map(Q, result["probabilities"].values())) == 1, "probability sum")
        equal(result["one_body_density"], [[encoded(z) for z in line] for line in identity],
              "all one-body readouts agree")
    # On either statistics' one-particle sector, a_i^dagger a_j = |i><j|.
    bilinears = [[[[encoded(ONE if r == i and c == j else ZERO)
                    for c in range(2)] for r in range(2)] for j in range(2)] for i in range(2)]
    equal(receipt["one_particle_bilinears"], bilinears, "one-particle matrix-unit restriction")
    swap = [[ONE if j == (2*(i % 2)+i//2) else ZERO for j in range(4)] for i in range(4)]
    check(mm(swap, tensor(u)) == mm(tensor(u), swap), "swap covariance")
    minus_identity = [[(Q(-1), Q(0)), ZERO], [ZERO, (Q(-1), Q(0))]]
    check(tensor(minus_identity) == [[ONE if i == j else ZERO for j in range(4)]
                                     for i in range(4)], "central spin sign on pair")
    equal(receipt["spin_rotation"], {
        "supplied_single_particle_central_sign": -1,
        "two_particle_sign_in_both_statistics": 1,
        "exchange_commutes_with_tensor_square": True,
    }, "spin rotation contract")
    gap = Q(results["fermion"]["probabilities"]["11"])-Q(results["boson"]["probabilities"]["11"])
    check(gap > 0, "strict exchange signal")
    partial = []
    for eta in [Q(0), Q(1, 4), Q(1)]:
        cases = {"fermion": hidden_label_result(u, -1, eta),
                 "boson": hidden_label_result(u, 1, eta)}
        difference = Q(cases["fermion"]["probabilities"]["11"])-Q(cases["boson"]["probabilities"]["11"])
        partial.append({"squared_hidden_state_overlap": str(eta), "results": cases,
                        "coincidence_gap": str(difference)})
    equal(receipt["partial_distinguishability"], partial, "partial distinguishability density replay")
    equal(receipt["gap"], {"fermion_minus_boson_coincidence": str(gap),
                           "strict_equal_error_threshold": str(gap/2)}, "gap and threshold")
    error = Q(1, 50)
    intervals = {k: [str(max(Q(0), Q(v["probabilities"]["11"])-error)),
                     str(min(Q(1), Q(v["probabilities"]["11"])+error))]
                 for k, v in results.items()}
    check(Q(intervals["fermion"][0]) > Q(intervals["boson"][1]), "disjoint model intervals")
    equal(receipt["error_contract"], {
        "preparation_trace_distance_at_most": "1/100",
        "readout_effect_operator_norm_error_at_most": "1/100",
        "total_probability_error_at_most": "1/50",
        "coincidence_intervals": intervals, "fermion_boson_intervals_disjoint": True,
        "errors_are_declared_bounds_not_measured": True, "transport_error_included": False,
    }, "declared error contract")
    equal(receipt["scope"], {
        "finite_declared_statistics_models": True, "exact_theoretical_probabilities": True,
        "one_particle_current_blindness_on_conserved_one_excitation_channels": True,
        "source_selected_preparation_or_statistics": False,
        "physical_spin_statistics_theorem": False, "A1_A3_no_go": False,
        "native_or_laboratory_detector_outcomes": False,
        "parent_quantum_link_operator_Gauss_inherited": False,
        "continuum_or_interacting_QFT_result": False,
    }, "scientific scope")
    return {"verified": True, "coincidence_gap": str(gap), "arithmetic": "exact rational"}


def set_to_list(value):
    check(isinstance(value, dict), "receipt must be an object")
    return sorted(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", nargs="?", type=Path, default=HERE/"receipt.json")
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True))


if __name__ == "__main__":
    main()
