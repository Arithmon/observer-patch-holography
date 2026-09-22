"""Independent ladder-polynomial and matrix verification of the finite witnesses."""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ZERO, ONE, I = (Q(0), Q(0)), (Q(1), Q(0)), (Q(0), Q(1))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def equal(a, b, message):
    require(json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True), message)


def pairs(rows):
    result = {}
    for key, value in rows:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def load(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def cj(a):
    return a[0], -a[1]


def sumz(zs):
    result = ZERO
    for z in zs:
        result = add(result, z)
    return result


def matmul(a, b):
    return [[sumz(mul(a[i][k], b[k][j]) for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def adjoint(a):
    return [[cj(a[j][i]) for j in range(len(a))] for i in range(len(a[0]))]


def ident(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def decode(a):
    return [[tuple(map(Q, z)) for z in row] for row in a]


def encode(a):
    return [[[str(z[0]), str(z[1])] for z in row] for row in a]


def integer_matrix(a):
    return [[(Q(x), Q(0)) for x in row] for row in a]


def matrix_add(a, b, sign=1):
    return [[add(a[i][j], mul((Q(sign), Q(0)), b[i][j]))
             for j in range(len(a[0]))] for i in range(len(a))]


def polymul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def polyvalue(p, x):
    out = Q(0)
    for coefficient in reversed(p):
        out = out*x+coefficient
    return out


def verify(data):
    equal(sorted(data), sorted([
        "schema", "source_pins", "word_certificate", "source_edge", "neutral_bosonic_descent",
        "CAR_two_mode_sea", "CAR_full_edge_sea", "finite_dimension_alternative",
        "truncated_boson_control", "escapes", "premises", "scope"
    ]), "schema fields")
    equal(data["schema"], "oph.pauli_stability.v1", "schema version")
    paths = ["code/sm_fermion_current/current_receipt.json", "code/sm_fermion_current/current.py",
             "paper/tex_fragments/FERMION_SOURCE_CURRENT.tex", "paper/tex_fragments/PAULI_STABILITY_SELECTION.tex",
             "code/pauli_stability/build.py",
             "code/pauli_stability/verify.py", "code/pauli_stability/test_pauli_stability.py"]
    equal(data["source_pins"], {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
          "source pins")
    parent = load(ROOT/paths[0])
    # Independently use the q-ladder [n+1]_q=1+q[n]_q. At level two,
    # norm=(1+q), <N>_unnormalized=(1+q)^2, so exact unit generation gives q²−1.
    norm = [1, 1]
    number = polymul(norm, norm)
    residual = [number[i]-(2*norm[i] if i < len(norm) else 0) for i in range(len(number))]
    equal(data["word_certificate"], {
        "alphabet": "a annihilation; c creation; coefficients ascending powers of real q",
        "rule": "ac = 1 + q ca", "vacuum_norm": 1,
        "double_creation_norm_squared": norm,
        "double_creation_number_expectation": number,
        "exact_generator_residual_expectation": residual,
        "generator_defect_normal_order": {"cca": [-1, 1]},
        "admissible_real_q_roots": [-1, 1],
        "q_controls": [{"q": str(q), "double_creation_norm_squared": str(polyvalue(norm, q)),
                        "generator_residual_expectation": str(polyvalue(residual, q))}
                       for q in [Q(-2), Q(-1), Q(0), Q(1, 2), Q(1)]],
    }, "scalar-q word and generator certificate")
    edge, kappa = parent["geometry"]["edges"][0], Q(parent["law"]["hopping"])
    require(kappa > 0, "strictly positive supplied hopping")
    paulis = [integer_matrix([[0, 1], [1, 0]]),
              [[ZERO, (Q(0), Q(-1))], [I, ZERO]], integer_matrix([[1, 0], [0, -1]])]
    K = [[mul(I, z) for z in row] for row in paulis[edge["axis"]]]
    h = [[ZERO for _ in range(4)] for _ in range(4)]
    for i in range(2):
        for j in range(2):
            h[i][j+2] = mul((kappa, Q(0)), K[i][j])
            h[i+2][j] = mul((kappa, Q(0)), cj(K[j][i]))
    modes = []
    for sign in [1, -1]:
        chi = [(Q(1, 2), Q(0)), (Q(sign, 2), Q(0))]
        bottom = [sumz(mul(cj(K[j][i]), chi[j]) for j in range(2)) for i in range(2)]
        for energy_sign in [-1, 1]:
            v = chi+[mul((Q(energy_sign), Q(0)), z) for z in bottom]
            modes.append({"energy": str(energy_sign*kappa),
                          "vector": [[str(z[0]), str(z[1])] for z in v]})
    equal(data["source_edge"], {
        "edge_index": 0, "ends": edge["ends"], "axis": edge["axis"], "hopping": str(kappa),
        "declared_unit_link_phase": ["1", "0"], "basis": "L spin0, L spin1, R spin0, R spin1",
        "hamiltonian": encode(h), "orthonormal_eigenmodes": modes,
        "source_selection_claimed": False, "edge_is_full_federation_generator": False,
    }, "supplied source edge and eigenmodes")
    require(h == adjoint(h), "Hermitian hopping")
    require(matmul(h, h) == [[mul((kappa*kappa, Q(0)), z) for z in row] for row in ident(4)],
            "h squared is kappa squared identity")
    vectors = [[tuple(map(Q, z)) for z in row["vector"]] for row in modes]
    for i, v in enumerate(vectors):
        require([sumz(mul(row[j], v[j]) for j in range(4)) for row in h] ==
                [mul((Q(modes[i]["energy"]), Q(0)), z) for z in v], "eigen equation")
        for j, w in enumerate(vectors):
            require(sumz(mul(cj(x), y) for x, y in zip(v, w)) == (ONE if i == j else ZERO),
                    "orthonormal source modes")
    fields = [{key: row[key] for key in ["name", "multiplicity", "integer_charge"]}
              for row in parent["multiplets"]]
    channels = sum(row["multiplicity"] for row in fields)
    charge = sum(row["multiplicity"]*row["integer_charge"] for row in fields)
    require(charge == 0 and channels == 15, "complete neutral source multiplets")
    equal(data["neutral_bosonic_descent"], {
        "multiplets": fields, "channels": channels, "charge_sum": charge,
        "energy_slope_per_common_occupation": str(-channels*kappa),
        "samples": [{"occupation_per_channel": n, "total_charge": n*charge,
                     "energy": str(-channels*kappa*n)} for n in [0, 1, 2, 4, 16]],
        "arbitrary_n_proof": "analytic Fock ladder; samples are witnesses only",
        "charge_chemical_potential_changes_energy": False,
        "operator_Gauss_or_gauge_field_energy_included": False,
    }, "neutral ladder and arbitrary occupation boundary")
    # Reconstruct actual finite CAR matrices on the four occupation bitmasks.
    creators = []
    for mode in range(2):
        c = [[ZERO for _ in range(4)] for _ in range(4)]
        for mask in range(4):
            if not (mask >> mode) & 1:
                sign = -1 if (mask & ((1 << mode)-1)).bit_count() % 2 else 1
                c[mask | (1 << mode)][mask] = (Q(sign), Q(0))
        creators.append(c)
    ns = [matmul(c, adjoint(c)) for c in creators]
    for i in range(2):
        for j in range(2):
            require(matrix_add(matmul(adjoint(creators[i]), creators[j]),
                               matmul(creators[j], adjoint(creators[i]))) ==
                    (ident(4) if i == j else [[ZERO]*4 for _ in range(4)]), "finite CAR witness")
    energy = [[mul((kappa, Q(0)), z) for z in row] for row in matrix_add(ns[1], ns[0], -1)]
    states = []
    for mask in range(4):
        nminus, nplus = int(ns[0][mask][mask][0]), int(ns[1][mask][mask][0])
        e = energy[mask][mask][0]
        require(e+kappa >= 0, "positive sea shifted energy")
        states.append({"occupation_minus_plus": [nminus, nplus], "energy": str(e),
                       "sea_shifted_energy": str(e+kappa), "particle_hole_count": nplus+1-nminus})
    equal(data["CAR_two_mode_sea"], {"single_particle_energies": [str(-kappa), str(kappa)],
          "states": states, "ground_energy": str(-kappa), "energy_shift": str(kappa)}, "CAR sea")
    equal(data["CAR_full_edge_sea"], {"negative_modes_per_internal_channel": 2,
          "occupied_negative_modes": 2*channels, "ground_energy": str(-2*channels*kappa),
          "total_charge": 2*charge}, "full supplied edge sea")
    c2, a2 = integer_matrix([[0, 0], [1, 0]]), integer_matrix([[0, 1], [0, 0]])
    require(matrix_add(matmul(a2, c2), matmul(c2, a2)) == ident(2), "two-dimensional CAR")
    require(sum(matrix_add(matmul(a2, c2), matmul(c2, a2), -1)[i][i][0]
                for i in range(2)) == 0, "finite commutator trace")
    equal(data["finite_dimension_alternative"], {
        "argument": "q=1 requires trace([a,a_dagger])=trace(I); finite matrix traces give 0=dimension",
        "dimension_two_CAR_witness": {"annihilation": [[0, 1], [0, 0]], "creation": [[0, 0], [1, 0]]},
        "finite_trace_controls": [{"dimension": d, "commutator_trace": 0, "identity_trace": d}
                                  for d in [1, 2, 3, 4]],
        "alternative_to_energy_stability_not_an_additional_requirement": True,
        "full_matter_Hilbert_finite_dimension_supplied": True,
        "A1_finite_response_implies_this_matter_realization": False,
    }, "finite matter dimension alternative")
    truncated = data["truncated_boson_control"]
    c, a = integer_matrix(truncated["creation"]), integer_matrix(truncated["annihilation"])
    gram = integer_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    n = matmul(c, a)
    require(matmul(gram, a) == matmul(adjoint(c), gram), "truncated positive metric adjoint")
    require(matrix_add(matmul(n, c), matmul(c, n), -1) == c, "truncated exact phase generator")
    require(matmul(c, c)[2][0] == ONE, "nonzero truncated double creation")
    residual_ccr = matrix_add(matrix_add(matmul(a, c), matmul(c, a), -1), ident(3), -1)
    require([residual_ccr[i][i][0] for i in range(3)] == [0, 0, -3], "truncated q1 top boundary")
    equal(truncated, {
        "occupation_basis": [0, 1, 2], "gram_diagonal": ["1", "1", "2"],
        "creation": [[0, 0, 0], [1, 0, 0], [0, 1, 0]], "annihilation": [[0, 1, 0], [0, 0, 2], [0, 0, 0]],
        "number_diagonal": [0, 1, 2], "negative_energy_diagonal": ["0", str(-kappa), str(-2*kappa)],
        "double_creation_vacuum_norm_squared": "2", "exact_number_generator": True,
        "no_single_real_q_satisfies_scalar_relation": True,
    }, "truncated boson escape")
    equal(data["escapes"], {
        "fixed_total_number_N_energy_lower_bound": "-kappa*N",
        "uniform_number_shift": {"coefficient": str(kappa), "shifted_one_particle_energies": ["0", str(2*kappa)],
            "preserves_fixed_N_continuous_number_conserving_readouts": True, "same_step_Cayley_map_inherited": False},
        "finite_occupation_caps_allow_nonzero_double_creation": True, "boundedness_alone_selects_CAR": False,
    }, "stability escape boundaries")
    equal(data["premises"], {
        "one_common_real_scalar_q_for_all_modes": True,
        "positive_Hilbert_metric_and_normalized_annihilation_vacuum": True,
        "common_invariant_algebraic_operator_domain": True,
        "number_is_a_dagger_a_with_exact_unit_phase_generator": True,
        "linear_creation_map_and_normalized_mode_relation": True,
        "unrestricted_ladder_domain_and_fixed_signed_generator": True,
        "reference_vacuum_is_generator_eigenvector": True,
        "negative_mode_creation_has_fixed_energy_commutator": True,
        "global_energy_lower_bound_on_that_domain": True, "these_premises_are_source_derived": False,
    }, "selector premises")
    equal(data["scope"], {
        "conditional_scalar_q_class_selection": True,
        "same_mode_Pauli_and_polarized_CAR_are_analytic_consequences": True,
        "finite_word_and_edge_witnesses_exact": True, "spin_blind_energy_selector": True,
        "physical_spin_statistics_theorem": False, "A1_A3_derive_CAR": False,
        "all_statistics_or_all_observer_models_excluded": False,
        "native_preparation_or_physical_energy_identified": False,
        "finite_samples_prove_unboundedness": False, "full_gauge_constrained_source_Hamiltonian": False,
    }, "scientific scope")
    return {"verified": True, "q_roots": [-1, 1], "stable_signed_full_ladder_q": -1,
            "finite_dimensional_q": -1, "source_premises_derived": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", nargs="?", default=HERE/"receipt.json")
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True))


if __name__ == "__main__":
    main()
