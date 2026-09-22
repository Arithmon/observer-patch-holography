"""Finite exact witnesses for the conditional scalar-q Pauli selector.

The all-occupation theorem is analytic. This certificate checks its short-word
algebra, its supplied negative source-edge mode, and its decisive escape cases.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = "code/sm_fermion_current/current_receipt.json"
PINS = [PARENT, "code/sm_fermion_current/current.py",
        "paper/tex_fragments/FERMION_SOURCE_CURRENT.tex",
        "paper/tex_fragments/PAULI_STABILITY_SELECTION.tex",
        "code/pauli_stability/build.py", "code/pauli_stability/verify.py",
        "code/pauli_stability/test_pauli_stability.py"]


def poly_add(a, b):
    out = [0]*max(len(a), len(b))
    for p in (a, b):
        for i, x in enumerate(p):
            out[i] += x
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


@lru_cache(None)
def normal(word):
    """Normal order by aa† = 1 + q a†a; c denotes creation."""
    at = word.find("ac")
    if at < 0:
        return {word: (1,)}
    result = dict(normal(word[:at]+word[at+2:]))
    for key, p in normal(word[:at]+"ca"+word[at+2:]).items():
        result[key] = poly_add(result.get(key, (0,)), (0,)+p)
    return result


def expectation(word):
    return normal(word).get("", (0,))


def value(poly, q):
    return sum(F(c)*q**i for i, c in enumerate(poly))


def enc(z):
    return [str(z.real), str(z.imag)]


def source_module():
    spec = importlib.util.spec_from_file_location("pauli_parent_current", ROOT/"code/sm_fermion_current/current.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def build():
    source = source_module()
    parent = json.loads((ROOT/PARENT).read_text())
    kappa = F(parent["law"]["hopping"])
    assert kappa == source.HOPPING and kappa > 0
    edge = parent["geometry"]["edges"][0]
    axis = edge["axis"]
    basis = [[source.C(F(i == j)) for i in range(2)] for j in range(2)]
    K = [source.spin_action(axis, source.C(1), v) for v in basis]
    Kadj = [source.spin_action(axis, source.C(1), v, adjoint=True) for v in basis]
    h = [[source.C() for _ in range(4)] for _ in range(4)]
    for i in range(2):
        for j in range(2):
            h[i][j+2] = source.C(kappa)*K[j][i]
            h[i+2][j] = source.C(kappa)*Kadj[j][i]
    modes = []
    for spin_sign in [1, -1]:
        chi = [source.C(F(1, 2)), source.C(F(spin_sign, 2))]
        tail = source.spin_action(axis, source.C(1), chi, adjoint=True)
        for energy_sign in [-1, 1]:
            v = chi+[source.C(F(energy_sign))*x for x in tail]
            modes.append({"energy": str(energy_sign*kappa), "vector": list(map(enc, v))})
    norm2, n_expectation = expectation("aacc"), expectation("aacacc")
    residual = poly_add(n_expectation, tuple(-2*x for x in norm2))
    generator_defect = dict(normal("cac"))
    for word in ["cca", "c"]:
        for key, p in normal(word).items():
            generator_defect[key] = poly_add(generator_defect.get(key, (0,)), tuple(-x for x in p))
    generator_defect = {key: list(p) for key, p in generator_defect.items() if any(p)}
    fields = [{key: row[key] for key in ["name", "multiplicity", "integer_charge"]}
              for row in parent["multiplets"]]
    channels = sum(row["multiplicity"] for row in fields)
    total_charge = sum(row["multiplicity"]*row["integer_charge"] for row in fields)
    assert total_charge == 0
    samples = [{"occupation_per_channel": n, "total_charge": n*total_charge,
                "energy": str(-channels*kappa*n)} for n in [0, 1, 2, 4, 16]]
    states = []
    for mask in range(4):
        occupied_minus, occupied_plus = mask & 1, (mask >> 1) & 1
        energy = kappa*(occupied_plus-occupied_minus)
        states.append({"occupation_minus_plus": [occupied_minus, occupied_plus],
                       "energy": str(energy), "sea_shifted_energy": str(energy+kappa),
                       "particle_hole_count": occupied_plus+1-occupied_minus})
    return {
        "schema": "oph.pauli_stability.v1",
        "source_pins": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},
        "word_certificate": {
            "alphabet": "a annihilation; c creation; coefficients ascending powers of real q",
            "rule": "ac = 1 + q ca",
            "vacuum_norm": 1,
            "double_creation_norm_squared": list(norm2),
            "double_creation_number_expectation": list(n_expectation),
            "exact_generator_residual_expectation": list(residual),
            "generator_defect_normal_order": generator_defect,
            "admissible_real_q_roots": [-1, 1],
            "q_controls": [{"q": str(q), "double_creation_norm_squared": str(value(norm2, q)),
                             "generator_residual_expectation": str(value(residual, q))}
                            for q in [F(-2), F(-1), F(0), F(1, 2), F(1)]],
        },
        "source_edge": {"edge_index": 0, "ends": edge["ends"], "axis": axis,
                        "hopping": str(kappa), "declared_unit_link_phase": ["1", "0"],
                        "basis": "L spin0, L spin1, R spin0, R spin1",
                        "hamiltonian": [[enc(z) for z in row] for row in h],
                        "orthonormal_eigenmodes": modes,
                        "source_selection_claimed": False,
                        "edge_is_full_federation_generator": False},
        "neutral_bosonic_descent": {"multiplets": fields, "channels": channels,
            "charge_sum": total_charge, "energy_slope_per_common_occupation": str(-channels*kappa),
            "samples": samples, "arbitrary_n_proof": "analytic Fock ladder; samples are witnesses only",
            "charge_chemical_potential_changes_energy": False,
            "operator_Gauss_or_gauge_field_energy_included": False},
        "CAR_two_mode_sea": {"single_particle_energies": [str(-kappa), str(kappa)],
                             "states": states, "ground_energy": str(-kappa),
                             "energy_shift": str(kappa)},
        "CAR_full_edge_sea": {"negative_modes_per_internal_channel": 2,
                              "occupied_negative_modes": 2*channels,
                              "ground_energy": str(-2*channels*kappa),
                              "total_charge": 2*total_charge},
        "finite_dimension_alternative": {
            "argument": "q=1 requires trace([a,a_dagger])=trace(I); finite matrix traces give 0=dimension",
            "dimension_two_CAR_witness": {"annihilation": [[0, 1], [0, 0]],
                                          "creation": [[0, 0], [1, 0]]},
            "finite_trace_controls": [{"dimension": d, "commutator_trace": 0,
                                       "identity_trace": d} for d in [1, 2, 3, 4]],
            "alternative_to_energy_stability_not_an_additional_requirement": True,
            "full_matter_Hilbert_finite_dimension_supplied": True,
            "A1_finite_response_implies_this_matter_realization": False,
        },
        "truncated_boson_control": {
            "occupation_basis": [0, 1, 2], "gram_diagonal": ["1", "1", "2"],
            "creation": [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            "annihilation": [[0, 1, 0], [0, 0, 2], [0, 0, 0]],
            "number_diagonal": [0, 1, 2], "negative_energy_diagonal": ["0", str(-kappa), str(-2*kappa)],
            "double_creation_vacuum_norm_squared": "2",
            "exact_number_generator": True, "no_single_real_q_satisfies_scalar_relation": True,
        },
        "escapes": {
            "fixed_total_number_N_energy_lower_bound": "-kappa*N",
            "uniform_number_shift": {"coefficient": str(kappa),
                                     "shifted_one_particle_energies": ["0", str(2*kappa)],
                                     "preserves_fixed_N_continuous_number_conserving_readouts": True,
                                     "same_step_Cayley_map_inherited": False},
            "finite_occupation_caps_allow_nonzero_double_creation": True,
            "boundedness_alone_selects_CAR": False,
        },
        "premises": {
            "one_common_real_scalar_q_for_all_modes": True,
            "positive_Hilbert_metric_and_normalized_annihilation_vacuum": True,
            "common_invariant_algebraic_operator_domain": True,
            "number_is_a_dagger_a_with_exact_unit_phase_generator": True,
            "linear_creation_map_and_normalized_mode_relation": True,
            "unrestricted_ladder_domain_and_fixed_signed_generator": True,
            "reference_vacuum_is_generator_eigenvector": True,
            "negative_mode_creation_has_fixed_energy_commutator": True,
            "global_energy_lower_bound_on_that_domain": True,
            "these_premises_are_source_derived": False,
        },
        "scope": {
            "conditional_scalar_q_class_selection": True,
            "same_mode_Pauli_and_polarized_CAR_are_analytic_consequences": True,
            "finite_word_and_edge_witnesses_exact": True,
            "spin_blind_energy_selector": True,
            "physical_spin_statistics_theorem": False,
            "A1_A3_derive_CAR": False,
            "all_statistics_or_all_observer_models_excluded": False,
            "native_preparation_or_physical_energy_identified": False,
            "finite_samples_prove_unboundedness": False,
            "full_gauge_constrained_source_Hamiltonian": False,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(build(), indent=2, sort_keys=True)+"\n"
    if args.write:
        (HERE/"receipt.json").write_text(data)
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
