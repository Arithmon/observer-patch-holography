"""Independent labeled-tensor verification of the capped symmetric controls.

The embedding uses metric 1/2 on degree-two labeled tensors. Thus x_i x_j
maps to |ij>+|ji> with no irrational coefficients, including i=j.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
Z, O, IM = (F(0), F(0)), (F(1), F(0)), (F(0), F(1))


def check(ok, why):
    if not ok:
        raise ValueError(why)


def equal(a, b, why):
    check(json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True), why)


def unique(items):
    out = {}
    for key, value in items:
        check(key not in out, "duplicate key")
        out[key] = value
    return out


def load(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=unique)


def plus(x, y):
    return x[0]+y[0], x[1]+y[1]


def times(x, y):
    return x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]


def conjugate(x):
    return x[0], -x[1]


def sumz(xs):
    out = Z
    for x in xs:
        out = plus(out, x)
    return out


def zeros(n, m=None):
    return [[Z for _ in range(n if m is None else m)] for _ in range(n)]


def eye(n):
    return [[O if i == j else Z for j in range(n)] for i in range(n)]


def decode(a, dimension):
    check(type(a) is list and len(a) == dimension, "matrix dimension")
    result = []
    for row in a:
        check(type(row) is list and len(row) == dimension, "matrix dimension")
        parsed = []
        for z in row:
            check(type(z) is list and len(z) == 2, "canonical complex pair")
            check(all(type(x) is str for x in z), "canonical fraction strings")
            try:
                value = tuple(map(F, z))
            except (ValueError, ZeroDivisionError):
                raise ValueError("canonical fraction strings") from None
            check([str(x) for x in value] == z, "canonical fraction strings")
            parsed.append(value)
        result.append(parsed)
    return result


def encode(a):
    return [[[str(z[0]), str(z[1])] for z in row] for row in a]


def adj(a):
    return [[conjugate(a[j][i]) for j in range(len(a))] for i in range(len(a[0]))]


def mm(a, b):
    out = zeros(len(a), len(b[0]))
    nz = [[(j, z) for j, z in enumerate(row) if z != Z] for row in b]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x != Z:
                for j, y in nz[k]:
                    out[i][j] = plus(out[i][j], times(x, y))
    return out


def linear(terms):
    a = terms[0][1]
    return [[sumz(times(c, m[i][j]) for c, m in terms) for j in range(len(a[0]))]
            for i in range(len(a))]


def comm(a, b):
    return linear([(O, mm(a, b)), ((F(-1), F(0)), mm(b, a))])


def field(v, operators):
    return linear(list(zip(v, operators)))


def labeled_creator(v):
    c = zeros(21)
    for i in range(4):
        c[1+i][0] = v[i]
        for j in range(4):
            c[5+4*i+j][1+j] = plus(c[5+4*i+j][1+j], v[i])
            c[5+4*j+i][1+j] = plus(c[5+4*j+i][1+j], v[i])
    return c


def verify(data):
    equal(sorted(data), sorted(["schema", "source_pins", "basis", "gram_diagonal", "creation", "annihilation",
        "source_edge", "quadratic_hamiltonian", "lifted_cayley", "normalized_mode_controls",
        "occupation_sectors", "scalar_q_failure", "three_outcome_instrument", "binary_instrument_route",
        "contractive_phase_route", "basis_only_control", "source_history_number_obstruction", "scope"]), "schema keys")
    equal(data["schema"], "oph.pauli_source_selection.v1", "schema")
    parent_path = "code/sm_fermion_current/current_receipt.json"
    history = "Lean/InformationProjection/SourceHistoryPacket.lean"
    paths = [parent_path, "code/sm_fermion_current/current.py", history,
             "paper/tex_fragments/PAULI_SOURCE_SELECTION_BOUNDARY.tex",
             "code/pauli_source_selection/build.py", "code/pauli_source_selection/verify.py",
             "code/pauli_source_selection/test_pauli_source_selection.py"]
    equal(data["source_pins"], {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, "source pins")
    parent = load(ROOT/parent_path)
    # Independent carrier constructed from unordered labeled pairs.
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    occupations = [[0]*4]+[[int(i == j) for i in range(4)] for j in range(4)]
    occupations += [[int(k == i)+int(k == j) for k in range(4)] for i, j in pairs]
    equal(data["basis"], occupations, "occupation basis")
    embedding = zeros(21, 15)
    for i in range(5):
        embedding[i][i] = O
    for col, (i, j) in enumerate(pairs, 5):
        embedding[5+4*i+j][col] = plus(embedding[5+4*i+j][col], O)
        embedding[5+4*j+i][col] = plus(embedding[5+4*j+i][col], O)
    weights = [F(1)]*5+[F(1, 2)]*16
    W = [[(weights[i], F(0)) if i == j else Z for j in range(21)] for i in range(21)]
    gram = mm(mm(adj(embedding), W), embedding)
    equal(data["gram_diagonal"], [str(gram[i][i][0]) for i in range(15)], "positive tensor Gram")
    check(all(gram[i][i][0] > 0 for i in range(15)), "positive metric")
    check(type(data["creation"]) is list and type(data["annihilation"]) is list and
          len(data["creation"]) == len(data["annihilation"]) == 4, "four modes")
    cs = [decode(c, 15) for c in data["creation"]]
    aa = [decode(a, 15) for a in data["annihilation"]]
    for i in range(4):
        v = [O if j == i else Z for j in range(4)]
        tc = labeled_creator(v)
        ta = [[times((weights[j]/weights[k], F(0)), conjugate(tc[j][k]))
               for j in range(21)] for k in range(21)]
        check(mm(embedding, cs[i]) == mm(tc, embedding), "tensor creation intertwiner")
        check(mm(embedding, aa[i]) == mm(ta, embedding), "tensor annihilation intertwiner")
        check(mm(gram, aa[i]) == mm(adj(cs[i]), gram), "positive-metric adjoint")
    edge = parent["geometry"]["edges"][0]
    k, step = F(parent["law"]["hopping"]), F(parent["law"]["step"])
    check(k > 0, "positive hopping")
    paulis = [[[Z, O], [O, Z]], [[Z, (F(0), F(-1))], [IM, Z]], [[O, Z], [Z, (F(-1), F(0))]]]
    K = [[times(IM, z) for z in row] for row in paulis[edge["axis"]]]
    h = zeros(4)
    for i in range(2):
        for j in range(2):
            h[i][j+2] = times((k, F(0)), K[i][j])
            h[i+2][j] = times((k, F(0)), conjugate(K[j][i]))
    u = decode(data["source_edge"]["cayley"], 4)
    # Defining Cayley equation independently of producer closed-form factors.
    lhs = linear([(O, eye(4)), ((F(0), step/2), h)])
    rhs = linear([(O, eye(4)), ((F(0), -step/2), h)])
    check(mm(lhs, u) == rhs and mm(adj(u), u) == eye(4), "source Cayley equation")
    expected_edge = {"index": 0, "ends": edge["ends"], "axis": edge["axis"], "hopping": str(k),
                     "step": str(step), "unit_link_phase": True, "hamiltonian": encode(h), "cayley": encode(u)}
    equal(data["source_edge"], expected_edge, "source edge")
    TH, TU = zeros(21), zeros(21)
    TU[0][0] = O
    for i in range(4):
        for j in range(4):
            TH[1+i][1+j], TU[1+i][1+j] = h[i][j], u[i][j]
    for i in range(4):
        for j in range(4):
            for a in range(4):
                for b in range(4):
                    TH[5+4*i+j][5+4*a+b] = sumz([h[i][a] if j == b else Z, h[j][b] if i == a else Z])
                    TU[5+4*i+j][5+4*a+b] = times(u[i][a], u[j][b])
    H, U = decode(data["quadratic_hamiltonian"], 15), decode(data["lifted_cayley"], 15)
    check(mm(embedding, H) == mm(TH, embedding), "tensor quadratic generator")
    check(mm(embedding, U) == mm(TU, embedding), "tensor unitary lift")
    check(mm(gram, H) == mm(adj(H), gram), "self-adjoint capped energy")
    check(mm(mm(adj(U), gram), U) == gram, "unitary covariance lift")
    vectors = [[O, Z, Z, Z], [(F(3, 5), F(0)), (F(4, 5), F(0)), Z, Z],
               [(F(1, 2), F(0)), (F(0), F(1, 2)), (F(1, 2), F(0)), (F(0), F(1, 2))],
               [(F(1, 2), F(0)), (F(1, 2), F(0)), (F(0), F(1, 2)), (F(0), F(1, 2))]]
    controls = []
    for v in vectors:
        c, a = field(v, cs), field(list(map(conjugate, v)), aa)
        n = mm(c, a)
        check(comm(n, c) == c, "exact mode phase law")
        hv = [sumz(times(z, x) for z, x in zip(row, v)) for row in h]
        uv = [sumz(times(z, x) for z, x in zip(row, v)) for row in u]
        check(comm(H, c) == field(hv, cs), "exact source generator lift")
        check(mm(U, c) == mm(field(uv, cs), U), "capped field covariance")
        double = [[row[0]] for row in mm(c, c)]
        norm = mm(mm(adj(double), gram), double)[0][0]
        check(norm == (F(2), F(0)), "double occupation survives")
        controls.append({"mode": [[str(z[0]), str(z[1])] for z in v],
                         "double_creation_norm_squared": "2", "number_phase_identity": True})
    equal(data["normalized_mode_controls"], controls, "normalized mode controls")
    spectrum = [{"energy": str(k*j), "multiplicity": m} for j, m in [(-2, 3), (-1, 2), (0, 5), (1, 2), (2, 3)]]
    powers = [eye(15)]
    for degree in range(1, 6):
        powers.append(mm(powers[-1], H))
    check(linear([(O, powers[5]), ((-5*k*k, F(0)), powers[3]), ((4*k**4, F(0)), powers[1])]) == zeros(15),
          "finite energy spectral polynomial")
    for degree in range(5):
        trace = sumz(powers[degree][i][i] for i in range(15))
        check(trace == (sum(F(row["energy"])**degree*row["multiplicity"] for row in spectrum), F(0)),
              "energy spectral multiplicities")
    equal(data["occupation_sectors"], {"dimensions": [1, 4, 10], "maximum_total_number": 2,
          "full_dimension": 15, "energy_spectrum_multiplicities": spectrum}, "bounded capped sectors")
    equal(data["scalar_q_failure"], {"vacuum_requires": "no restriction", "single_v_occupation_requires_q": "1",
          "double_v_occupation_requires_q": "-1/2", "one_common_q_exists": False}, "scalar-q failure")
    n0, ac = mm(cs[0], aa[0]), mm(aa[0], cs[0])
    check(n0[1][1] == O and ac[1][1] == (F(2), F(0)), "single occupancy forces q1")
    check(n0[5][5] == (F(2), F(0)) and ac[5][5] == Z, "double occupancy forces q minus half")
    effects = linear([((F(1, 4), F(0)), n0), ((F(1, 4), F(0)), ac)])
    idle = [1-effects[i][i][0] for i in range(15)]
    check(all(x > 0 for x in idle), "positive idle instrument effect")
    equal(data["three_outcome_instrument"], {"selected_coordinate_mode": 0,
        "addition_and_removal_amplitude_scale": "1/2", "idle_effect_diagonal": list(map(str, idle)),
        "idle_Kraus": "positive diagonal square root of idle_effect_diagonal",
        "minimum_idle_effect": str(min(idle)), "effect_sum_is_identity": True,
        "exhaustive_two_adjoint_outcomes": False, "unscaled_field_number_law_inherited_by_scaled_branch": False},
        "three-outcome instrument")
    equal(data["binary_instrument_route"], {
        "if_exactly_a_and_a_dagger_are_complete_Kraus_then_diagonal_CAR": True,
        "exact_field_branch_identification_and_two_outcome_completeness_supplied": True, "source_derived": False},
        "conditional binary instrument route")
    scaled_c = linear([((F(1, 2), F(0)), cs[0])])
    scaled_n = linear([((F(1, 4), F(0)), n0)])
    check(comm(scaled_n, scaled_c) == linear([((F(1, 4), F(0)), scaled_c)]) and
          comm(scaled_n, scaled_c) != scaled_c, "scaled branch fails unit phase")
    c2 = [[Z, Z], [O, Z]]
    n2 = mm(c2, adj(c2))
    check(comm(n2, c2) == c2 and mm(c2, c2) == zeros(2), "contractive phase witness")
    check(linear([(O, n2), (O, mm(adj(c2), c2))]) == eye(2), "paired binary completeness")
    equal(data["contractive_phase_route"], {"two_dimensional_creator": [[0, 0], [1, 0]],
        "number_diagonal": [0, 1], "creator_squared_zero": True, "paired_adjoint_effects_sum_to_identity": True,
        "capped_unscaled_creator_norm_squared": "2", "capped_half_scaled_creator_norm_squared": "1/2",
        "capped_half_scaled_number_phase_coefficient": "1/4",
        "same_physical_creation_is_a_contraction_and_unit_phase_raiser_supplied": True,
        "source_derivation_claimed": False}, "contractive phase route")
    bit0 = [[Z, Z, Z, Z], [O, Z, Z, Z], [Z, Z, Z, Z], [Z, Z, O, Z]]
    models = {}
    for name, final_sign in [("commuting_bits", 1), ("signed_CAR", -1)]:
        bit1 = [[Z, Z, Z, Z], [Z, Z, Z, Z], [O, Z, Z, Z], [Z, (F(final_sign), F(0)), Z, Z]]
        for c in [bit0, bit1]:
            n = mm(c, adj(c))
            check(mm(n, n) == n and comm(n, c) == c, "basis contraction and phase")
        c = linear([((F(3, 5), F(0)), bit0), ((F(4, 5), F(0)), bit1)])
        n = mm(c, adj(c))
        check(mm(n, n) == n, "mixed field remains a contraction")
        square = mm(c, c)
        defect = linear([(O, comm(n, c)), ((F(-1), F(0)), c)])
        models[name] = {"creation": [encode(bit0), encode(bit1)],
                        "mixed_double_creation_vacuum_amplitude": [str(x) for x in square[3][0]],
                        "mixed_phase_defect_on_mask1_amplitude": [str(x) for x in defect[3][1]],
                        "sequential_basis_creation_occupation_probability": "1"}
        if name == "signed_CAR":
            for mode in [bit0, bit1, c, linear([((F(3, 5), F(0)), bit0), ((F(0), F(4, 5)), bit1)])]:
                number = mm(mode, adj(mode))
                check(comm(number, mode) == mode and
                      linear([(O, number), (O, mm(adj(mode), mode))]) == eye(4),
                      "CAR real and imaginary pair acceptance")
        else:
            check(square[3][0] == (F(24, 25), F(0)) and defect[3][1] == (F(-72, 125), F(0)),
                  "basis-only false green")
    equal(data["basis_only_control"], {"basis_order": ["00", "10", "01", "11"],
        "real_mixing_coefficients": ["3/5", "4/5"], "models": models,
        "individual_basis_phase_and_contraction_pass_both": True,
        "signed_CAR_basis_and_real_imaginary_pair_acceptance": True,
        "occupation_record_probabilities_select_relative_sign": False}, "basis-only control")
    # Source history is the eight binary triples; action counts adjacent changes.
    actions = [int((i >> 2) != ((i >> 1) & 1))+int(((i >> 1) & 1) != (i & 1)) for i in range(8)]
    source = (ROOT/history).read_text()
    declared = list(map(int, re.search(r"def sourceAction : Fin 8 → ℕ := !\[([^]]+)\]", source).group(1).split(",")))
    check(declared == actions, "source repair count")
    dims = [actions.count(i) for i in range(3)]
    check(dims == [2, 4, 2] and dims[1] > dims[0], "isometry rank obstruction")
    equal(data["source_history_number_obstruction"], {"sourceAction": actions, "eigenspace_dimensions": dims,
        "hypothesis": "repair-count matrix equals a_dagger*a and [N,a_dagger]=a_dagger",
        "lowering_from_energy_one_must_be_isometric": True,
        "required_domain_dimension": dims[1], "available_target_dimension": dims[0],
        "hypothesis_inconsistent_on_this_eight_dimensional_carrier": True,
        "all_source_Hamiltonians_or_enlarged_carriers_excluded": False}, "source-history number obstruction")
    equal(data["scope"], {"finite_positive_capped_symmetric_model": True,
        "linear_creation_with_nonzero_double_occupation": True, "same_supplied_one_particle_source_hopping": True,
        "finite_controls_for_analytic_all_mode_covariance_and_phase": True,
        "all_modes_and_all_unitaries_exhaustively_numerically_tested": False,
        "scalar_q_law_derived_from_source": False, "full_A1_A3_countermodel": False,
        "physical_spin_statistics_theorem": False,
        "native_physical_matter_or_instrument": False, "local_quantum_Gauss_or_relativistic_field_locality": False,
        "same_parent_fifteen_particle_state": False, "positive_binary_instrument_route_is_conditional": True}, "scope")
    return {"verified": True, "finite_source_controls": True, "source_premises_derived": False,
            "physical_spin_statistics": False, "dimension": 15, "normalized_mode_controls": 4}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", nargs="?", default=HERE/"receipt.json")
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True))


if __name__ == "__main__":
    main()
