"""The field instrument continues the existing work and every old record."""
import math

import numpy as np
import pytest

from . import response, run
from ..mass_response import experiment as e


@pytest.fixture(scope="module")
def continued():
    original = e.execute(backend="continuous")
    actual = run.append_transfer(original)
    target = run.append_transfer(original, backend="continuous")
    return original, actual, target


def test_complete_continuing_field_outputs_match_independent_instrument(continued):
    before, actual, target = continued
    assert len(actual.branches) == 64
    assert e.quantum_distance(actual, target) < 2e-11
    assert abs(math.fsum(actual.probabilities().values())-1) < 2e-11
    assert actual.resources["work_preparations"] == 1
    assert actual.specification["work_reset_count"] == 0
    assert actual.specification["additional_work_preparations"] == 0
    for record, old in before.probabilities().items():
        new = math.fsum(e.norm(state) for key, state in actual.branches.items() if key[:-1] == record)
        assert abs(new-old) < 2e-12
    assert abs(e.recover_gap_ratio(actual)["finite_single_charged_gap_over_clock"]-
               e.recover_gap_ratio(before)["finite_single_charged_gap_over_clock"]) < 2e-12
    forbidden = np.array([(mask&255).bit_count() != 4 for mask in range(1024)])
    assert math.fsum(e.norm(state[forbidden]) for state in actual.branches.values()) < 1e-24
    # Complete local instruments keep the inaccessible reference marginal.
    old_reference = sum(state.conj().T@state for state in before.branches.values())
    new_reference = sum(state.conj().T@state for state in actual.branches.values())
    np.testing.assert_allclose(old_reference, new_reference, atol=2e-12, rtol=0)


def test_source_field_instrument_does_not_consult_target_matrices(monkeypatch):
    original = e.execute(backend="continuous")
    def reject(*args, **kwargs):
        raise AssertionError("source field read consulted a target matrix")
    monkeypatch.setattr(response, "transfer_operator", reject)
    monkeypatch.setattr(response, "transfer_kraus", reject)
    actual = run.append_transfer(original)
    assert len(actual.branches) == 64


def test_frozen_public_witnesses_are_checked_without_running_writing_producers():
    packets, bindings = run._read_public()
    assert set(packets) == set(run.PUBLIC_NAMES)
    assert len(bindings) > 6
    assert not packets["representation_identification"]["scope"]["laboratory_species_identified"]
