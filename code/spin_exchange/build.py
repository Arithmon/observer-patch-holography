"""Exact creation-polynomial statistics comparison on a supplied two-port factor.

This computes theoretical probabilities. It does not simulate detector outcomes
or select a statistics law from observer axioms.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = "code/sm_fermion_current/quantum_link_receipt.json"


@dataclass(frozen=True)
class C:
    re: F = F(0)
    im: F = F(0)

    def __add__(self, other):
        return C(self.re + other.re, self.im + other.im)

    def __mul__(self, other):
        return C(self.re*other.re-self.im*other.im,
                 self.re*other.im+self.im*other.re)

    def conj(self):
        return C(self.re, -self.im)

    def norm(self):
        return self.re*self.re+self.im*self.im

    def encode(self):
        return [str(self.re), str(self.im)]


def create(state, mode, fermion):
    result = {}
    for occ, coeff in state.items():
        if fermion and occ[mode]:
            continue
        sign = -1 if fermion and sum(occ[:mode]) % 2 else 1
        out = list(occ)
        out[mode] += 1
        key = tuple(out)
        result[key] = result.get(key, C()) + C(F(sign))*coeff
    return result


def annihilate(state, mode, fermion):
    result = {}
    for occ, coeff in state.items():
        if not occ[mode]:
            continue
        factor = ((-1)**sum(occ[:mode])) if fermion else occ[mode]
        out = list(occ)
        out[mode] -= 1
        result[tuple(out)] = C(F(factor))*coeff
    return result


def inner(left, right, fermion):
    answer = C()
    for occ, x in left.items():
        weight = 1 if fermion else factorial(occ[0])*factorial(occ[1])
        answer = answer + C(F(weight))*x.conj()*right.get(occ, C())
    return answer


def linear_create(state, column, fermion):
    answer = {}
    for mode, coefficient in enumerate(column):
        for occ, value in create(state, mode, fermion).items():
            answer[occ] = answer.get(occ, C()) + coefficient*value
    return answer


def fock_result(matrix, fermion):
    vacuum = {(0, 0): C(F(1))}
    # c_L^dagger c_R^dagger |0>: rightmost operator acts first.
    state = linear_create(vacuum, [row[1] for row in matrix], fermion)
    state = linear_create(state, [row[0] for row in matrix], fermion)
    norm = inner(state, state, fermion)
    assert norm == C(F(1))
    probabilities = {}
    for occ in [(2, 0), (1, 1), (0, 2)]:
        weight = 1 if fermion else factorial(occ[0])*factorial(occ[1])
        probabilities[str(occ[0])+str(occ[1])] = str(weight*state.get(occ, C()).norm())
    rho = [[inner(state, create(annihilate(state, j, fermion), i, fermion),
                  fermion).encode() for j in range(2)] for i in range(2)]
    return {"probabilities": probabilities, "one_body_density": rho,
            "norm_squared": str(norm.re)}


def one_particle_bilinears(fermion):
    basis = [{(1, 0): C(F(1))}, {(0, 1): C(F(1))}]
    return [[[ [inner(basis[row], create(annihilate(basis[col], j, fermion),
                    i, fermion), fermion).encode() for col in range(2)]
                    for row in range(2)] for j in range(2)] for i in range(2)]


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def build():
    parent = json.loads((ROOT/PARENT).read_text())
    channel = next(row for row in parent["channels"] if row["multiplet"] == "e_c")
    step, kappa = F(channel["step"]), F(channel["kappa"])
    r = step*kappa/2
    a, b = (1-r*r)/(1+r*r), 2*r/(1+r*r)
    assert channel["left_amplitude"] == [str(a), "0"]
    assert channel["right_amplitude"] == ["0", str(-b)]
    matrix = [[C(a), C(F(0), -b)], [C(F(0), -b), C(a)]]
    fermion, boson = fock_result(matrix, True), fock_result(matrix, False)
    distinguishable = {
        "probabilities": {"20": str(a*a*b*b), "11": str(a**4+b**4),
                          "02": str(a*a*b*b)},
        "one_body_density": [[["1", "0"], ["0", "0"]],
                             [["0", "0"], ["1", "0"]]],
        "norm_squared": "1",
    }
    one_boson, one_fermion = one_particle_bilinears(False), one_particle_bilinears(True)
    assert one_boson == one_fermion
    gap = F(fermion["probabilities"]["11"])-F(boson["probabilities"]["11"])
    preparation, detector = F(1, 100), F(1, 100)
    error = preparation+detector
    intervals = {}
    for name, result in [("fermion", fermion), ("boson", boson),
                         ("distinguishable", distinguishable)]:
        p = F(result["probabilities"]["11"])
        intervals[name] = [str(max(F(0), p-error)), str(min(F(1), p+error))]
    distinguishability = []
    for eta in [F(0), F(1, 4), F(1)]:
        cases = {}
        for name, ideal in [("fermion", fermion), ("boson", boson)]:
            cases[name] = {
                "probabilities": {k: str(eta*F(v)+(1-eta)*F(distinguishable["probabilities"][k]))
                                  for k, v in ideal["probabilities"].items()},
                "one_body_density": ideal["one_body_density"], "norm_squared": "1",
            }
        distinguishability.append({"squared_hidden_state_overlap": str(eta),
                                   "results": cases,
                                   "coincidence_gap": str(eta*gap)})
    return {
        "schema": "oph.spin_exchange.v1",
        "source_pins": {path: sha(path) for path in [PARENT,
            "code/sm_fermion_current/quantum_link.py",
            "code/sm_fermion_current/current.py",
            "paper/tex_fragments/FERMION_SOURCE_CURRENT.tex",
            "code/spin_exchange/build.py", "code/spin_exchange/verify.py",
            "code/spin_exchange/test_spin_exchange.py"]},
        "binding": {"parent_channel_id": channel["id"], "multiplet": "e_c",
                    "step": str(step), "kappa": str(kappa),
                    "reuse": "Cayley coefficients only; no parent quantum-link Gauss or state transfer"},
        "splitter": [[z.encode() for z in row] for row in matrix],
        "preparation": {
            "spatial_occupation": [1, 1],
            "identical_internal_and_spin_channel": True,
            "spin_frame": "same declared channel in local port frames; raw spin vectors may differ",
            "two_particles_in_same_channel": True,
            "distinguishable_control": "orthogonal unobserved internal labels",
            "normalization": "normalized symmetric/antisymmetric LR tensor; labeled LR for control",
        },
        "results": {"fermion": fermion, "boson": boson,
                    "distinguishable": distinguishable},
        "partial_distinguishability": distinguishability,
        "one_particle_bilinears": one_boson,
        "spin_rotation": {
            "supplied_single_particle_central_sign": -1,
            "two_particle_sign_in_both_statistics": 1,
            "exchange_commutes_with_tensor_square": True,
        },
        "gap": {"fermion_minus_boson_coincidence": str(gap),
                "strict_equal_error_threshold": str(gap/2)},
        "error_contract": {
            "preparation_trace_distance_at_most": str(preparation),
            "readout_effect_operator_norm_error_at_most": str(detector),
            "total_probability_error_at_most": str(error),
            "coincidence_intervals": intervals,
            "fermion_boson_intervals_disjoint": intervals["boson"][1] != intervals["fermion"][0]
                and F(intervals["boson"][1]) < F(intervals["fermion"][0]),
            "errors_are_declared_bounds_not_measured": True,
            "transport_error_included": False,
        },
        "scope": {
            "finite_declared_statistics_models": True,
            "exact_theoretical_probabilities": True,
            "one_particle_current_blindness_on_conserved_one_excitation_channels": True,
            "source_selected_preparation_or_statistics": False,
            "physical_spin_statistics_theorem": False,
            "A1_A3_no_go": False,
            "native_or_laboratory_detector_outcomes": False,
            "parent_quantum_link_operator_Gauss_inherited": False,
            "continuum_or_interacting_QFT_result": False,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build(), indent=2, sort_keys=True)+"\n"
    if args.write:
        (HERE/"receipt.json").write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
