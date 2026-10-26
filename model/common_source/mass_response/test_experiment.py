from dataclasses import replace
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import dispersion, experiment as e, history, model


@pytest.fixture(scope="module")
def source_history():
    return e.execute(steps=32)


@pytest.mark.parametrize("interaction", [False, True])
@pytest.mark.parametrize("interval", [F(1, 4), F(1, 2)])
def test_continuous_record_inverse_matches_independent_charged_spectrum(interaction, interval):
    result = e.execute(backend="continuous", interaction=interaction, interval=interval)
    p = model.parameters()
    exact = dispersion.ring_gap_predictions(coupling=p.interaction if interaction else 0.)
    inferred = e.recover_gap_ratio(result)
    assert abs(inferred["finite_single_charged_gap_over_clock"]-exact["charged_addition_gap"]/p.clock) < 2e-13
    assert len(result.branches) == 32
    assert not result.aborts
    assert abs(sum(result.probabilities().values())-1) < 2e-14


def test_one_history_contains_interaction_dressing_and_later_propagation():
    interacting = e.execute(backend="continuous")
    free = e.execute(backend="continuous", interaction=False)
    no_kick = e.execute(backend="continuous", kick=False)
    pinned = e.execute(backend="continuous", propagation=False)
    late_off = e.execute(backend="continuous", second_interval_interaction=False)
    assert e.recover_gap_ratio(interacting)["finite_single_charged_gap_over_clock"] > e.recover_gap_ratio(free)["finite_single_charged_gap_over_clock"] + .2
    assert abs(e.marginal(interacting, 2)-e.marginal(no_kick, 2)) > .2
    assert abs(e.marginal(interacting, 2)-e.marginal(pinned, 2)) > .2
    for control in (no_kick, pinned):
        assert abs(e.marginal(interacting, 1)-e.marginal(control, 1)) < 1e-14
    assert abs(e.marginal(interacting, 2)-e.marginal(late_off, 2)) > .0005
    assert interacting.prefix_probabilities == late_off.prefix_probabilities


def test_complete_source_outputs_match_independent_continuous_history(source_history):
    target = e.execute(backend="continuous")
    bound = 2*model.split_operator_bound(F(1, 8), 32)+model.split_operator_bound(F(1, 4), 32)
    assert e.quantum_distance(source_history, target) <= float(bound)+2e-8
    assert abs(sum(source_history.probabilities().values())+sum(map(e.norm, source_history.aborts.values()))-1) < 1e-10
    assert len(source_history.branches) == 32 and len(source_history.aborts) == 12
    assert sum(map(e.norm, source_history.aborts.values())) < 1e-20


def test_source_clock_and_charged_response_are_read_from_records(source_history, monkeypatch):
    expected = dispersion.ring_gap_predictions()["charged_addition_gap"]/model.parameters().clock
    def reject(*args, **kwargs):
        raise AssertionError("observer consulted the hidden target spectrum")
    monkeypatch.setattr(model, "ring_hamiltonian", reject)
    hidden_duration = replace(source_history, specification={**source_history.specification,
                             "interval_exact": "not available to the decoder"})
    actual = e.recover_gap_ratio(hidden_duration)
    assert abs(actual["finite_single_charged_gap_over_clock"]-expected) < 6e-6
    assert abs(actual["measured_clock_phase"]-math.log(5)/2) < 2e-10
    assert not actual["infinite_volume_mass_identification"]


def test_later_source_operations_preserve_every_committed_prefix(source_history):
    probabilities = source_history.probabilities()
    for prefix, old in source_history.prefix_probabilities.items():
        new = math.fsum(p for record, p in probabilities.items() if record[:2] == prefix)
        assert abs(new-old) < 2e-11
    assert source_history.resources["work_preparations"] == 1
    assert source_history.specification["work_reset_count"] == 0
    assert source_history.resources["native_pulses"] > 0
    assert source_history.resources["proper_code_service_events"] == (
        source_history.resources["proper_code_load_events"] + source_history.resources["proper_code_export_events"])
    assert source_history.resources["maximum_logical_qubits_including_untouched_reference"] == 15


def test_total_charge_and_untouched_reference_are_retained(source_history):
    forbidden = np.array([(mask&255).bit_count() != 4 for mask in range(1024)])
    assert all(state.shape == (1024, 4) for state in source_history.branches.values())
    leaked = math.fsum(float(np.vdot(state[forbidden], state[forbidden]).real) for state in source_history.branches.values())
    assert leaked < 1e-20


def test_inconclusive_clock_or_failure_is_not_silently_renormalized():
    result = e.execute(backend="continuous")
    altered = replace(result, aborts={("failure",): np.ones((1, 1), complex)/10})
    with pytest.raises(ValueError, match="inconclusive"):
        e.recover_gap_ratio(altered)
    altered = replace(result, branches={record: value*math.sqrt(.8) for record, value in result.branches.items()})
    with pytest.raises(ValueError, match="calibrated circle"):
        e.recover_gap_ratio(altered)


@pytest.mark.parametrize("kwargs", [{"interval": F(0)}, {"interval": F(3, 4)}, {"interval": .25},
                                    {"steps": True}, {"backend": "spectrum-only"}, {"interaction": 1}])
def test_invalid_or_ambiguous_protocols_are_rejected(kwargs):
    with pytest.raises(ValueError):
        e.execute(**kwargs)


def test_cq_distance_includes_phase_invariant_reference_outputs_and_failures():
    original = e.execute(backend="continuous")
    phased = replace(original, branches={k: np.exp(.37j)*v for k, v in original.branches.items()})
    assert e.quantum_distance(original, phased) < 1e-14
    lost = replace(original, branches={k: .9*v for k, v in original.branches.items()})
    assert abs(e.quantum_distance(original, lost)-.095) < 2e-14
