"""Boundary tests prevent silent widening and false program-level promotion."""
from fractions import Fraction as F
import json

import pytest

from .family import CONTRACT, FamilyParameters, load_contract


@pytest.mark.parametrize("weight", [F(0), F(1, 2), F(1)])
def test_both_native_laws_and_mixture_round_trip(weight):
    parameters = FamilyParameters(c_weight=weight)
    assert FamilyParameters.from_dict(parameters.as_dict()) == parameters


@pytest.mark.parametrize("parameters", [
    {"c_weight": 0.5}, {"c_weight": True}, {"c_weight": F(-1, 2)},
    {"c_weight": F(1, 257)}, {"repair_probability": F(0)},
    {"repair_probability": F(1, 17)}, {"exploration_probability": F(0)},
    {"memory_capacity": True}, {"memory_capacity": 9}, {"overflow": "erase_silently"},
])
def test_outside_class_is_rejected(parameters):
    with pytest.raises(ValueError):
        FamilyParameters(**parameters)


def test_retained_sources_and_parameter_contract_are_bound():
    contract = load_contract()
    assert contract["boundary_status"] == "FROZEN_FINITE_CAMPAIGN"
    assert contract["whole_program_frozen"] is False


@pytest.mark.parametrize("damage", ["promotion", "comparison", "floor", "custody"])
def test_contract_damage_fails(tmp_path, damage):
    data = json.loads(CONTRACT.read_text())
    if damage == "promotion":
        data["whole_program_frozen"] = True
    elif damage == "comparison":
        data["physical_comparison_count_added"] = 1
    elif damage == "floor":
        data["parameter_boundary"]["repair_probability_min"] = "0"
    else:
        data["source_bindings"][0]["sha256"] = "0" * 64
    path = tmp_path / "campaign.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        load_contract(path)


@pytest.mark.parametrize("damage", ["duplicate", "lost_C", "lost_E", "boundary"])
def test_named_endpoints_and_finite_boundary_cannot_silently_change(tmp_path, damage):
    data = json.loads(CONTRACT.read_text())
    if damage == "duplicate":
        data["initial_members"].append(data["initial_members"][0])
    elif damage == "boundary":
        data["boundary_status"] = "FROZEN"
    else:
        name = damage.removeprefix("lost_")
        next(row for row in data["initial_members"] if row["name"] == name)["parameters"]["c_weight"] = "1/2"
    path = tmp_path / "campaign.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        load_contract(path)
