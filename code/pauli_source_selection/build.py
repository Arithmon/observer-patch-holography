"""Positive finite source-law controls; no source-selected statistics claim."""
from __future__ import annotations
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = "code/sm_fermion_current/current_receipt.json"
HISTORY = "Lean/InformationProjection/SourceHistoryPacket.lean"
PINS = [PARENT, "code/sm_fermion_current/current.py", HISTORY,
        "paper/tex_fragments/PAULI_SOURCE_SELECTION_BOUNDARY.tex",
        "code/pauli_source_selection/build.py", "code/pauli_source_selection/verify.py",
        "code/pauli_source_selection/test_pauli_source_selection.py"]


@dataclass(frozen=True)
class C:
    r: F = F(0)
    i: F = F(0)

    def __add__(self, other):
        return C(self.r+other.r, self.i+other.i)

    def __mul__(self, other):
        return C(self.r*other.r-self.i*other.i, self.r*other.i+self.i*other.r)

    def conj(self):
        return C(self.r, -self.i)

    def encode(self):
        return [str(self.r), str(self.i)]


Z, ONE = C(), C(F(1))


def zero(n, m=None):
    return [[Z for _ in range(n if m is None else m)] for _ in range(n)]


def product(a, b):
    out = zero(len(a), len(b[0]))
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x == Z:
                continue
            for j, y in enumerate(b[k]):
                if y != Z:
                    out[i][j] = out[i][j]+x*y
    return out


def combine(terms):
    size = len(terms[0][1])
    return [[sum((s*a[i][j] for s, a in terms), Z) for j in range(size)] for i in range(size)]


def encoded(a):
    return [[z.encode() for z in row] for row in a]


def basis():
    out = [(0,)*4]
    out += [tuple(int(i == j) for i in range(4)) for j in range(4)]
    out += [tuple(int(k == i)+int(k == j) for k in range(4))
            for i in range(4) for j in range(i, 4)]
    return out


def creators_annihilators(occupations):
    index = {x: i for i, x in enumerate(occupations)}
    creators, annihilators = [], []
    for mode in range(4):
        c, a = zero(15), zero(15)
        for j, occ in enumerate(occupations):
            if sum(occ) < 2:
                target = list(occ)
                target[mode] += 1
                c[index[tuple(target)]][j] = ONE
            if occ[mode]:
                target = list(occ)
                target[mode] -= 1
                a[index[tuple(target)]][j] = C(F(occ[mode]))
        creators.append(c)
        annihilators.append(a)
    return creators, annihilators


def lift_unitary(u, occupations):
    index = {occ: i for i, occ in enumerate(occupations)}
    result = zero(15)
    for j, occ in enumerate(occupations):
        terms = {(0,)*4: ONE}
        for mode, count in enumerate(occ):
            for _ in range(count):
                new = {}
                for power, coefficient in terms.items():
                    for i in range(4):
                        target = list(power)
                        target[i] += 1
                        key = tuple(target)
                        new[key] = new.get(key, Z)+coefficient*u[i][mode]
                terms = new
        for target, coefficient in terms.items():
            result[index[target]][j] = coefficient
    return result


def source_module():
    spec = importlib.util.spec_from_file_location("pauli_source_parent", ROOT/"code/sm_fermion_current/current.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def basis_only_control():
    models = {}
    for name, signed in [("commuting_bits", False), ("signed_CAR", True)]:
        creators = []
        for mode in range(2):
            c = zero(4)
            for mask in range(4):
                if not (mask & (1 << mode)):
                    sign = -1 if signed and mode == 1 and (mask & 1) else 1
                    c[mask | (1 << mode)][mask] = C(F(sign))
            creators.append(c)
        mixed = combine([(C(F(3, 5)), creators[0]), (C(F(4, 5)), creators[1])])
        square = product(mixed, mixed)
        adjoint = [[mixed[j][i].conj() for j in range(4)] for i in range(4)]
        n = product(mixed, adjoint)
        defect = combine([(ONE, product(n, mixed)), (C(F(-1)), product(mixed, n)), (C(F(-1)), mixed)])
        models[name] = {"creation": list(map(encoded, creators)),
                        "mixed_double_creation_vacuum_amplitude": square[3][0].encode(),
                        "mixed_phase_defect_on_mask1_amplitude": defect[3][1].encode(),
                        "sequential_basis_creation_occupation_probability": "1"}
    return {"basis_order": ["00", "10", "01", "11"], "real_mixing_coefficients": ["3/5", "4/5"],
            "models": models, "individual_basis_phase_and_contraction_pass_both": True,
            "signed_CAR_basis_and_real_imaginary_pair_acceptance": True,
            "occupation_record_probabilities_select_relative_sign": False}


def build():
    parent = json.loads((ROOT/PARENT).read_text())
    native = source_module()
    edge = parent["geometry"]["edges"][0]
    kappa, step = F(parent["law"]["hopping"]), F(parent["law"]["step"])
    assert native.HOPPING == kappa and native.STEP == step
    h = zero(4)
    for j in range(2):
        v = [native.C(F(i == j)) for i in range(2)]
        upper = native.spin_action(edge["axis"], native.C(1), v)
        lower = native.spin_action(edge["axis"], native.C(1), v, adjoint=True)
        for i in range(2):
            h[i][j+2] = C(upper[i].real*kappa, upper[i].imag*kappa)
            h[i+2][j] = C(lower[i].real*kappa, lower[i].imag*kappa)
    r = step*kappa/2
    diagonal, off = (1-r*r)/(1+r*r), -2*r/(kappa*(1+r*r))
    u = [[C(diagonal if i == j else F(0))+C(F(0), off)*h[i][j]
          for j in range(4)] for i in range(4)]
    occupations = basis()
    creators, annihilators = creators_annihilators(occupations)
    gram = [F(factorial(occ[0])*factorial(occ[1])*factorial(occ[2])*factorial(occ[3]))
            for occ in occupations]
    H = combine([(h[i][j], product(creators[i], annihilators[j])) for i in range(4) for j in range(4)])
    gamma = lift_unitary(u, occupations)
    vectors = [[ONE, Z, Z, Z], [C(F(3, 5)), C(F(4, 5)), Z, Z],
               [C(F(1, 2)), C(F(0), F(1, 2)), C(F(1, 2)), C(F(0), F(1, 2))],
               [C(F(1, 2)), C(F(1, 2)), C(F(0), F(1, 2)), C(F(0), F(1, 2))]]
    controls = []
    for v in vectors:
        c = combine([(v[i], creators[i]) for i in range(4)])
        a = combine([(v[i].conj(), annihilators[i]) for i in range(4)])
        n = product(c, a)
        double = [row[0] for row in product(c, c)]
        norm = sum((gram[i]*(z.r*z.r+z.i*z.i) for i, z in enumerate(double)), F(0))
        controls.append({"mode": [z.encode() for z in v], "double_creation_norm_squared": str(norm),
                         "number_phase_identity": product(n, c) == combine([(ONE, product(c, n)), (ONE, c)])})
    sum_effect = combine([(ONE, product(creators[0], annihilators[0])),
                          (ONE, product(annihilators[0], creators[0]))])
    idle = [1-sum_effect[i][i].r/4 for i in range(15)]
    code = (ROOT/HISTORY).read_text()
    actions = list(map(int, re.search(r"def sourceAction : Fin 8 → ℕ := !\[([^]]+)\]", code).group(1).split(",")))
    return {
        "schema": "oph.pauli_source_selection.v1",
        "source_pins": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},
        "basis": [list(occ) for occ in occupations], "gram_diagonal": list(map(str, gram)),
        "creation": list(map(encoded, creators)), "annihilation": list(map(encoded, annihilators)),
        "source_edge": {"index": 0, "ends": edge["ends"], "axis": edge["axis"],
                        "hopping": str(kappa), "step": str(step), "unit_link_phase": True,
                        "hamiltonian": encoded(h), "cayley": encoded(u)},
        "quadratic_hamiltonian": encoded(H), "lifted_cayley": encoded(gamma),
        "normalized_mode_controls": controls,
        "occupation_sectors": {"dimensions": [1, 4, 10], "maximum_total_number": 2,
                               "full_dimension": 15,
                               "energy_spectrum_multiplicities": [
                                   {"energy": str(-2*kappa), "multiplicity": 3},
                                   {"energy": str(-kappa), "multiplicity": 2},
                                   {"energy": "0", "multiplicity": 5},
                                   {"energy": str(kappa), "multiplicity": 2},
                                   {"energy": str(2*kappa), "multiplicity": 3}]},
        "scalar_q_failure": {"vacuum_requires": "no restriction", "single_v_occupation_requires_q": "1",
                             "double_v_occupation_requires_q": "-1/2", "one_common_q_exists": False},
        "three_outcome_instrument": {
            "selected_coordinate_mode": 0, "addition_and_removal_amplitude_scale": "1/2",
            "idle_effect_diagonal": list(map(str, idle)),
            "idle_Kraus": "positive diagonal square root of idle_effect_diagonal",
            "minimum_idle_effect": str(min(idle)), "effect_sum_is_identity": True,
            "exhaustive_two_adjoint_outcomes": False,
            "unscaled_field_number_law_inherited_by_scaled_branch": False},
        "binary_instrument_route": {
            "if_exactly_a_and_a_dagger_are_complete_Kraus_then_diagonal_CAR": True,
            "exact_field_branch_identification_and_two_outcome_completeness_supplied": True,
            "source_derived": False},
        "contractive_phase_route": {
            "two_dimensional_creator": [[0, 0], [1, 0]],
            "number_diagonal": [0, 1], "creator_squared_zero": True,
            "paired_adjoint_effects_sum_to_identity": True,
            "capped_unscaled_creator_norm_squared": "2",
            "capped_half_scaled_creator_norm_squared": "1/2",
            "capped_half_scaled_number_phase_coefficient": "1/4",
            "same_physical_creation_is_a_contraction_and_unit_phase_raiser_supplied": True,
            "source_derivation_claimed": False},
        "basis_only_control": basis_only_control(),
        "source_history_number_obstruction": {
            "sourceAction": actions, "eigenspace_dimensions": [actions.count(i) for i in range(3)],
            "hypothesis": "repair-count matrix equals a_dagger*a and [N,a_dagger]=a_dagger",
            "lowering_from_energy_one_must_be_isometric": True,
            "required_domain_dimension": actions.count(1), "available_target_dimension": actions.count(0),
            "hypothesis_inconsistent_on_this_eight_dimensional_carrier": True,
            "all_source_Hamiltonians_or_enlarged_carriers_excluded": False},
        "scope": {
            "finite_positive_capped_symmetric_model": True,
            "linear_creation_with_nonzero_double_occupation": True,
            "same_supplied_one_particle_source_hopping": True,
            "finite_controls_for_analytic_all_mode_covariance_and_phase": True,
            "all_modes_and_all_unitaries_exhaustively_numerically_tested": False,
            "scalar_q_law_derived_from_source": False,
            "physical_spin_statistics_theorem": False,
            "full_A1_A3_countermodel": False,
            "native_physical_matter_or_instrument": False,
            "local_quantum_Gauss_or_relativistic_field_locality": False,
            "same_parent_fifteen_particle_state": False,
            "positive_binary_instrument_route_is_conditional": True},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(build(), indent=2, sort_keys=True)+"\n"
    if args.write:
        (HERE/"receipt.json").write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
