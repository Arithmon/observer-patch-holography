"""Independent native-effect and action checks for the proposed energy map."""
from fractions import Fraction as F

import pytest
import sympy as sp

from . import energy as e
from ..completion import observer


@pytest.mark.parametrize("law", ("C", "E"))
def test_complete_quantum_effects_give_the_reported_activity(law):
    # Different producer: sum the actual complete Kraus effects in H78,
    # then restrict to binary66. No simulated first-exit measurement is used.
    values = [F() for _ in observer.STATES]
    normalization = [F() for _ in observer.STATES]
    for seam in observer.SEAMS:
        for _, _, current, pairs in observer._raw_maps(law, seam):
            for i, _ in pairs:
                values[i] += F(abs(current), 60)
                normalization[i] += F(1, 60)
    assert normalization == [1]*78
    actual = e.activity_table(law)
    for row in actual:
        state = row["state"]
        assert row["absolute_current"] == values[observer.STATES.index(state)]
        occupied = tuple(i for i, n in enumerate(state) if n)
        # Boundary-edge count is an independent graph calculation.
        cut = sum((a in occupied) != (b in occupied) for a, b in observer.SEAMS)
        assert row["absolute_current"] == F(cut, 60)


def test_candidate_reference_subtraction_is_signed_and_unconditioned():
    result = e.maintenance_certificate()
    assert result["carrier_dimension"] == 66 and result["cluster_states"] == 30
    assert result["background_cost_per_attempt"]["absolute_current"] == F(5, 33)
    assert result["cluster_cost_per_attempt"]["absolute_current"] == F(2, 15)
    assert result["signed_excess_cost_per_attempt"]["absolute_current"] == F(-1, 55)
    assert result["signed_excess_cost_per_attempt"]["attempts"] == 0
    # Absolute values after subtracting the reference would reverse the result.
    assert not result["positive_action_scale_times_this_excess_is_positive"]


def test_actual_transported_references_also_have_negative_excess():
    result = e.transported_reference_activity()
    assert not result["adjacency_conditioning_is_executed_measurement"]
    assert [tuple(row["core"]) for row in result["rows"]] == [(0, 1), (0, 8), (1, 5), (5, 8)]
    # Independent geometric charge-weight formula; no history compiler.
    weights = (F(2, 225), F(1, 36), F(2, 45))
    for row in result["rows"]:
        core = row["core"]
        expected = F()
        for state in e.native.states(binary=True):
            weight = weights[sum(state[p] for p in core)]
            cut = sum(state[a] != state[b] for a, b in observer.SEAMS)
            expected += weight*F(cut, 60)
        assert row["reference_activity_per_attempt"] == expected
        assert row["adjacency_conditioned_activity_per_attempt"] == F(2, 15)
        assert row["signed_excess_activity_per_attempt"] == F(2, 15)-expected < 0


@pytest.mark.parametrize("law", ("C", "E"))
def test_closed_source_cycle_keeps_the_record_and_positive_path_activity(law):
    result = e.source_cycle(law)
    assert result["unsigned_activity"] == 2
    assert result["conditional_coin_word_probability"] == F(1, 4)
    assert result["source_state_function_energy_increment"] == 0
    assert not result["entire_source_plus_archive_returns_to_initial_state"]
    state = result["initial_and_final_load"]
    # Direct endpoint oracle, independent of native.native and record recovery.
    for a, b, coin, _ in result["complete_records"]:
        out = list(state)
        total = out[a]+out[b]
        out[a], out[b] = (total+coin)//2, (total+1-coin)//2
        state = tuple(out)
    assert state == result["initial_and_final_load"]


def test_equal_cost_does_not_identify_the_actual_modular_energy():
    witness = e.modular_cost_witness()
    i, j = witness["source_indices"]
    assert witness["reference_weight_ratio_j_over_i"] != 1
    for law in ("C", "E"):
        rows = e.activity_table(law)
        assert rows[i]["absolute_current"] == rows[j]["absolute_current"]
        assert rows[i]["attempts"] == rows[j]["attempts"]
    assert not witness["modular_energy_is_physical_rest_energy_claimed"]


def test_action_coefficients_follow_independent_symbolic_differentiation():
    v = sp.Symbol("v")
    for rest, speed in ((F(7, 3), F(5, 2)), (F(11, 13), F(3, 7))):
        er, cr = sp.Rational(rest.numerator, rest.denominator), sp.Rational(speed.numerator, speed.denominator)
        lagrangian = -er*sp.sqrt(1-v*v/(cr*cr))
        a, b = e.proper_time_coefficients(rest, speed)
        assert sp.diff(lagrangian, v, 2).subs(v, 0)/2 == a
        assert sp.diff(lagrangian, v, 4).subs(v, 0)/24 == b
        result = e.kinetic_consistency(rest, a, speed, b)
        assert result["rest_inertia_match"] and result["quartic_match"]


def test_independently_changed_kinetics_is_not_hidden_by_mass_definition():
    controls = e.action_controls()["controls"]
    assert controls["same_rest_energy_twice_inertia"]["rest_to_inertia_ratio"] == F(1, 2)
    assert not controls["same_rest_energy_twice_inertia"]["rest_inertia_match"]
    weak = controls["rest_inertia_equality_but_wrong_quartic"]
    assert weak["rest_inertia_match"] and not weak["quartic_match"]


@pytest.mark.parametrize("units", ((F(11, 7), F(3, 5), F(17, 13)), (F(1, 1000), F(1000), F(1, 1000))))
def test_unit_choices_cannot_remove_a_kinetic_mismatch(units):
    rest, speed = F(7, 3), F(5, 2)
    a, b = e.proper_time_coefficients(rest, speed)
    original = e.kinetic_consistency(rest, 2*a, speed, b)
    changed = e.scale_coefficients(rest, 2*a, speed, b, energy_unit=units[0], length_unit=units[1], time_unit=units[2])
    scaled = e.kinetic_consistency(*changed[:3], quartic=changed[3])
    for key in ("rest_to_inertia_ratio", "quartic_to_relativistic_ratio", "rest_inertia_match", "quartic_match"):
        assert scaled[key] == original[key]


@pytest.mark.parametrize("bad", (True, 1.0, float("nan"), "1", F(0), F(-1)))
def test_kinetic_inputs_require_positive_exact_normalizations(bad):
    with pytest.raises(ValueError):
        e.kinetic_consistency(bad, F(1), F(1))


def test_corrupt_fraction_and_caller_alias_are_rejected_or_detached():
    corrupt = F(1, 2)
    corrupt._denominator = 0
    with pytest.raises(ValueError):
        e.rational(corrupt)
    original = F(3, 2)
    copied = e.rational(original)
    original._numerator = 9
    assert copied == F(3, 2)


def test_source_bound_receipt_replays():
    e.check_receipt(e.HERE/"energy_results.json", e.report())
