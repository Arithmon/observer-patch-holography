"""Exact checks of the admissible record-law classification."""
from fractions import Fraction as F
import json
from pathlib import Path

import pytest

from tools import strict_json

from . import family as f

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def result():
    return f.run()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def test_retained_receipt_replays_exactly(result):
    assert canonical(strict_json.load(HERE/"results.json")) == canonical(result)


def test_surd_arithmetic_and_exact_sign():
    root = f.S(0, 1)
    assert root*root == 2 and (f.ROOT_HALF*f.ROOT_HALF) == F(1, 2)
    assert f.S(F(3, 8), F(-1, 4)).sign() == 1      # 3/8 - sqrt2/4 > 0
    assert f.S(1, -1).sign() == -1 and f.S(-3, 2).sign() == -1 and f.S(-2, 2).sign() == 1
    assert f.S(0, 0, 1)*f.S(0, 0, 1) == -1
    with pytest.raises(ValueError):
        f.S(0, 0, 1).sign()
    with pytest.raises(ValueError):
        f.power(f.table(1, 0, 0, 0, 0), 0)


def test_named_members_are_admissible_and_match_executed_instruments(result):
    members = result["members"]
    assert set(members) == {"C0", "C", "E", "T", "C_signed"}
    for row in members.values():
        assert row["complete"] and row["definite_current_per_branch"]
        assert row["number_phase_covariant"] and row["positive_semidefinite"]
        assert row["channel_equals_table_on_all_matrix_units"]
    for name in ("C", "E"):
        assert members[name]["equals_executed_instrument"]
        assert members[name]["equals_executed_prediction"]
    assert [members[n]["record_values_on_seam"] for n in ("C0", "C", "E", "T")] == [3, 6, 10, 12]


def test_interference_reads_the_waiting_coherence(result):
    plus = {n: result["members"][n]["wait_coherence_plus"] for n in result["members"]}
    assert plus["C"] == {"phase_0": "5/8", "phase_pi": "1/8", "phase_pi_over_2": "3/8",
                         "diagonal_control": "3/8"}
    assert set(plus["E"].values()) == {"3/8"} == set(plus["T"].values())
    assert plus["C0"]["phase_0"] == "3/8+1/4*sqrt2" and plus["C0"]["phase_pi"] == "3/8-1/4*sqrt2"
    assert plus["C_signed"]["phase_0"] == "1/8"
    full = {n: result["members"][n]["full_coherence_plus"]["phase_0"] for n in ("C0", "C", "E", "T")}
    assert full == {"C0": "1", "C": "1", "E": "1", "T": "1/2"}


def test_ladder_completion_and_idempotence(result):
    assert all(result["ladder"].values())
    assert result["idempotent_members"] == ["E", "T"]
    assert {row["limit"] for row in result["completion"].values()} == {"E"}
    assert result["C_power_identity_C^r=E+2^(1-r)(C-E)_for_r_1_to_8"]


def test_composition_law_on_an_arbitrary_point():
    t = f.generic_tables()["C0_quarter_turn"]
    u = f.generic_tables()["C_T_mixture"]
    for x, y in [(f.representative(a), f.representative(b)) for a in f.KEYS for b in f.KEYS]:
        unit = {(x, y): f.S(1)}
        assert f.apply_table(u, f.apply_table(t, unit)) == f.apply_table(f.compose(t, u), unit)


def test_controls_and_bounds(result):
    control = result["prevalue_control"]
    assert control["complete"] and control["same_classical_law_as_C"]
    assert not control["number_phase_covariant"] and not control["definite_current_per_branch"]
    assert result["same_classical_law_for_every_member"]
    partition = result["coin_partition_laws"]
    assert (partition["count"], partition["distinct_channels"], partition["idempotent"],
            partition["with_C_interference"]) == (10, 10, 2, 4)
    assert all(v for k, v in result["bounds"].items() if k != "wait_coherence_bound_squared")
    assert not f.positive_semidefinite(f.table(1, F(1, 2), 0, 0, 0))   # k=1 forces g=h
