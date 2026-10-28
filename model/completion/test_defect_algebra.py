"""Independent exact witnesses for the source-charge obstruction.

Endpoint balancing below does not call the native transition implementation.
The tests check the proof objects and scope, rather than a floating nullity.
"""
import json
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys

import pytest

from model.completion import defect_algebra as algebra
from model.observer import native


@pytest.fixture(scope="module")
def receipt():
    return algebra.report()


def _oracle(state, branch):
    a, b = branch["seam"]
    low, parity = divmod(state[a]+state[b], 2)
    output = list(state)
    output[a] = low + branch["coin"]*parity
    output[b] = low + (1-branch["coin"])*parity
    current = state[a]-output[a]
    raw = {0: 0, 1: 1, -1: 2}[current]
    if branch["law"] == "E":
        raw += 4*parity
    return tuple(output), raw


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("domain", ["full78", "binary66"])
def test_projector_products_separate_all_basis_states_from_independent_effects(receipt, law, domain):
    certificate = receipt["C_E"][law][domain]
    states = native.states()
    expected = []
    for i, _ in certificate["effect_signatures"]:
        signature = sum((_oracle(states[i], branch)[1] == branch["raw"]) << bit
                        for bit, branch in enumerate(certificate["separating_effects"]))
        expected.append([i, str(signature)])
    assert expected == certificate["effect_signatures"]
    assert len({signature for _, signature in expected}) == certificate["dimension"]
    assert certificate["separated_basis_pairs"] == len(expected)*(len(expected)-1)//2
    # For each pattern the product of supports/complements is exactly P_i.
    for i, pattern in expected:
        assert [j for j, other in expected if other == pattern] == [i]


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("domain", ["full78", "binary66"])
def test_forest_constructs_all_matrix_units_with_native_or_adjoint_edges(receipt, law, domain):
    certificate = receipt["C_E"][law][domain]
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    links = {i: set() for i, _ in certificate["effect_signatures"]}
    for edge in certificate["spanning_forest"]:
        output, raw = _oracle(states[edge["native_source"]], edge["branch"])
        assert index[output] == edge["native_target"]
        assert raw == edge["branch"]["raw"]
        pair = (edge["native_source"], edge["native_target"])
        if edge["adjoint_for_parent_to_child"]:
            pair = pair[::-1]
        assert pair == (edge["parent"], edge["child"])
        links[edge["parent"]].add(edge["child"])
        links[edge["child"]].add(edge["parent"])
    assert len(certificate["spanning_forest"]) == len(links)-1
    reached, pending = {min(links)}, [min(links)]
    for source in pending:
        for target in links[source]-reached:
            reached.add(target)
            pending.append(target)
    assert reached == set(links)
    assert certificate["generated_star_algebra_complex_dimension"] == len(links)**2
    assert certificate["commutant_complex_dimension"] == 1


@pytest.mark.parametrize("law", ["C", "E"])
def test_invariant_binary_face_is_not_misreported_as_a_conserved_sector(receipt, law):
    control = receipt["C_E"][law]["invariant_not_reducing_control"]
    assert control["binary_to_collision_addresses"] == 0
    assert control["collision_to_binary_addresses"] == 120
    assert control["binary_subspace_invariant"]
    assert not control["binary_projector_commutes_with_every_branch"]
    witness = control["commutator_witness"]
    states = native.states()
    assert max(states[witness["source"]]) == 2
    assert max(states[witness["target"]]) == 1
    output, raw = _oracle(states[witness["source"]], witness["branch"])
    assert output == states[witness["target"]]
    assert raw == witness["branch"]["raw"]


@pytest.mark.parametrize("law", ["C", "E"])
def test_removing_incident_repairs_really_restores_a_nontrivial_charge(receipt, law):
    control = receipt["C_E"][law]["held_port_control"]
    certificate = control["certificate"]
    states = native.states()
    assert certificate["commutant_complex_dimension"] == 3
    assert certificate["generated_star_algebra_complex_dimension"] == 1+11**2+66**2
    assert sorted(map(len, certificate["components"])) == [1, 11, 66]
    for component in certificate["components"]:
        assert len({states[i][0] for i in component}) == 1
    assert control["outside_campaign_positive_exploration"]


def test_diagonal_distinguishability_without_transport_does_not_imply_full_algebra():
    rows = [{"law": "test", "seam": [], "coin": 0, "raw": i, "pairs": [[i, i]]}
            for i in range(3)]
    certificate = algebra.certify_maps(rows, (0, 1, 2))
    assert certificate["components"] == [[0], [1], [2]]
    assert certificate["commutant_complex_dimension"] == 3
    assert certificate["generated_star_algebra_complex_dimension"] == 3


def test_identity_only_control_cannot_supply_the_diagonal_separation_premise():
    row = {"law": "hold", "seam": [], "coin": 0, "raw": 0,
           "pairs": [[0, 0], [1, 1], [2, 2]]}
    with pytest.raises(ValueError, match="do not separate"):
        algebra.certify_maps([row], (0, 1, 2))


def test_full_branch_algebra_does_not_exclude_general_expectation_fixed_observables():
    # Absorbing0,2; middle1 jumps to either with probability1/2. The full
    # Kraus *-algebra is M3, but Q=diag(0,1/2,1) is harmonic in expectation.
    maps = [([[0, 0]], F(1)), ([[2, 2]], F(1)),
            ([[1, 0]], F(1, 2)), ([[1, 2]], F(1, 2))]
    rows = [{"law": "control", "seam": [], "coin": 0, "raw": i, "pairs": pairs}
            for i, (pairs, _) in enumerate(maps)]
    certificate = algebra.certify_maps(rows, (0, 1, 2))
    assert certificate["generated_star_algebra_complex_dimension"] == 9
    assert certificate["commutant_complex_dimension"] == 1
    charge = (F(0), F(1, 2), F(1))
    effects, evolved = [F(0)]*3, [F(0)]*3
    for pairs, weight in maps:
        for source, target in pairs:
            effects[source] += weight
            evolved[source] += weight*charge[target]
    assert effects == [1, 1, 1]
    assert tuple(evolved) == charge
    assert any(charge[source] != charge[target] for pairs, _ in maps for source, target in pairs)


def test_restricted_domain_must_be_invariant():
    row = {"law": "test", "seam": [], "coin": 0, "raw": 0, "pairs": [[0, 1]]}
    with pytest.raises(ValueError, match="not invariant"):
        algebra.certify_maps([row], (0,))


def test_retained_receipt_replays_exactly_and_preserves_scope(receipt):
    path = Path(algebra.__file__).with_name("defect_algebra_results.json")
    assert json.loads(path.read_text()) == receipt
    assert receipt["E_refines_C"]["identities_checked"] == 240
    assert not receipt["family_extension"]["requires_IID_scheduler"]
    assert not receipt["family_extension"]["new_adaptive_mixing_rate_claimed"]
    assert receipt["scope"]["arbitrary_quantum_basis_for_source_charge"]
    assert receipt["scope"]["exact_source_only_reducing_sectors_excluded"]
    assert not receipt["scope"]["arbitrary_Heisenberg_fixed_observable_classification_from_this_certificate"]
    for key in ("all_source_invariant_subspaces_excluded", "controller_or_history_encoded_charges_excluded",
                "metastability_or_controlled_lifetimes_excluded", "refined_or_federated_defects_excluded",
                "full_CP03_gate_accepted", "whole_campaign_or_program_falsified"):
        assert receipt["scope"][key] is False


def test_cli_rejects_scope_promotion_without_rewriting_input(receipt, tmp_path):
    changed = json.loads(json.dumps(receipt))
    changed["scope"]["whole_campaign_or_program_falsified"] = True
    path = tmp_path/"incorrect.json"
    original = json.dumps(changed)
    path.write_text(original)
    process = subprocess.run([sys.executable, "-m", "model.completion.defect_algebra", "--verify", str(path)],
                             cwd=algebra.ROOT, capture_output=True, text=True)
    assert process.returncode != 0
    assert "scope or source custody changed" in process.stderr
    assert path.read_text() == original
