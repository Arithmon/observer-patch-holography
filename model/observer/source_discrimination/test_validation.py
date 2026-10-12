"""The replay gate must check every retained record, not only its marginal."""
from copy import deepcopy
import math

import numpy as np
import pytest

from . import circuit, experiment


@pytest.fixture(scope="module")
def diagonal_events():
    return {candidate: experiment.events(candidate, "diagonal_control")
            for candidate in circuit.RECORD_BITS}


@pytest.mark.parametrize("candidate", ("C", "E"))
@pytest.mark.parametrize("mutation", (
    "missing_field", "extra_field", "raw_code", "bool_raw_code", "bool_coin",
    "bool_current", "bool_outcome", "integer_abort", "erased_preparation",
    "wrong_contrast", "drop_empty_branch", "duplicate_empty_branch",
    "split_duplicates", "column_amplitude", "masked_amplitude", "nan_amplitude",
))
def test_terminal_gate_rejects_invalid_records_and_incomplete_branch_grid(
        diagonal_events, candidate, mutation):
    rows = deepcopy(diagonal_events[candidate])
    admitted = next(row for row in rows if not row["abort"] and row["current"] == 0)
    if mutation == "missing_field":
        rows[0].pop("record_code")
    elif mutation == "extra_field":
        rows[0]["undeclared_record"] = 0
    elif mutation == "raw_code":
        for row in rows:
            row["record_code"] = 999
    elif mutation == "bool_raw_code":
        admitted["record_code"] = bool(admitted["record_code"])
    elif mutation == "bool_coin":
        admitted["coin"] = bool(admitted["coin"])
    elif mutation == "bool_current":
        admitted["current"] = False
    elif mutation == "bool_outcome":
        next(row for row in rows if row["outcome"] == 1)["outcome"] = True
    elif mutation == "integer_abort":
        admitted["abort"] = 0
    elif mutation == "erased_preparation":
        for row in rows:
            row["preparation_record"] = None
    elif mutation == "wrong_contrast":
        admitted["contrast"] = 99
    elif mutation in ("drop_empty_branch", "duplicate_empty_branch"):
        index = next(i for i, row in enumerate(rows) if row["abort"])
        assert experiment.probability(rows[index]) < 1e-20
        if mutation == "drop_empty_branch":
            rows.pop(index)
        else:
            rows.append(deepcopy(rows[index]))
    elif mutation == "split_duplicates":
        # Both total mass and every density marginal are unchanged.
        rows = [{**row, "amplitude": row["amplitude"] / math.sqrt(2)}
                for row in rows for _ in range(2)]
    elif mutation == "column_amplitude":
        admitted["amplitude"] = admitted["amplitude"].reshape(512, 1)
    elif mutation == "masked_amplitude":
        admitted["amplitude"] = np.ma.array(admitted["amplitude"], mask=True)
    else:
        admitted["amplitude"][0] = np.nan
    with pytest.raises((ValueError, ArithmeticError)):
        experiment._check_events(candidate, "diagonal_control", rows)


@pytest.mark.parametrize("candidate", ("C", "E"))
def test_diagonal_gate_checks_conditional_densities_not_only_fair_marginals(
        diagonal_events, candidate):
    rows = diagonal_events[candidate]
    lookup = {(row["preparation_record"], row["coin"], row["record_code"],
               row["outcome"]): row for row in rows}
    changed = []
    maximum_conditional_change = 0.
    for row in rows:
        partner = lookup[1 - row["preparation_record"], row["coin"],
                         row["record_code"], row["outcome"]]
        a, b = row["amplitude"], partner["amplitude"]
        amplitude = ((a + b) if row["preparation_record"] == 0 else (b - a)) / math.sqrt(2)
        changed.append({**row, "amplitude": amplitude})
        maximum_conditional_change = max(maximum_conditional_change, float(np.linalg.norm(
            np.outer(amplitude, amplitude.conj()) - np.outer(a, a.conj()))))
        if row["preparation_record"] == 0:
            # A unitary change of the two amplitudes preserves the density
            # after forgetting the preparation record, for EACH old block.
            u, v = (a + b) / math.sqrt(2), (a - b) / math.sqrt(2)
            original = np.outer(a, a.conj()) + np.outer(b, b.conj())
            mixed = np.outer(u, u.conj()) + np.outer(v, v.conj())
            assert np.linalg.norm(original - mixed) < 1e-14
    assert maximum_conditional_change > .12
    for preparation_record in (0, 1):
        marginal = sum(experiment.probability(row) for row in changed
                       if row["preparation_record"] == preparation_record)
        assert abs(marginal - .5) < experiment.TOLERANCE
    with pytest.raises(ArithmeticError):
        experiment._check_events(candidate, "diagonal_control", changed)


@pytest.mark.parametrize("candidate", ("C", "E"))
def test_independent_global_phases_of_committed_branches_are_unobservable(
        diagonal_events, candidate):
    rows = [{**row, "amplitude": row["amplitude"] * np.exp(1j * (index + 1) / 7)}
            for index, row in enumerate(diagonal_events[candidate])]
    checked = experiment._check_events(candidate, "diagonal_control", rows)
    assert checked["maximum_branch_density_frobenius_error"] < experiment.TOLERANCE
    assert checked["maximum_retained_preparation_density_frobenius_error"] < experiment.TOLERANCE


def test_run_rejects_corrupted_raw_record_metadata(monkeypatch):
    original = experiment.events

    def corrupted(candidate, setting):
        return [{**row, "record_code": 999} for row in original(candidate, setting)]

    monkeypatch.setattr(experiment, "events", corrupted)
    with pytest.raises((ValueError, ArithmeticError)):
        experiment.run()


def test_returned_resource_values_cannot_poison_later_receipts():
    first = experiment.run()
    baseline = deepcopy(first["candidates"]["C"]["executions"][0]["resources"])
    external = experiment.resources("C", "phase_0")
    try:
        # A caller may freely annotate or edit its own returned receipt.
        # Neither direct resource results nor receipt dictionaries may be
        # aliases of the private accounting cache used by a later replay.
        external["native_pulses"] = 0
        first["candidates"]["C"]["executions"][0]["resources"]["proper_code_load_events"] = 0
        second = experiment.run()
        assert second["all_pass"]
        assert second["candidates"]["C"]["executions"][0]["resources"] == baseline
        assert experiment.resources("C", "phase_0") == baseline
    finally:
        # Keep even a failing pre-fix regression run from poisoning other tests.
        external.clear()
        external.update(baseline)
