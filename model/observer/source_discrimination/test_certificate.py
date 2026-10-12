"""Independent full-carrier, complement-readout and error-composition checks."""
from fractions import Fraction as F

import pytest

from model.predictions import protocol
from . import certificate, circuit


@pytest.mark.parametrize("candidate,records", (("C", 6), ("E", 10)))
def test_exact_intertwining_checks_every_native_matrix_unit(candidate, records):
    proof = certificate.channel_certificate(candidate)
    assert proof["native_dimension"] == 78
    assert proof["matrix_units_checked"] == 78 * 78
    assert proof["kraus_entries"] == 156
    assert proof["record_count"] == records
    assert proof["floor_ceiling_intertwining_exact"]
    assert proof["arbitrary_reference_extension"]


@pytest.mark.parametrize("candidate", ("C", "E"))
def test_certificate_rejects_an_endpoint_word_mutation(monkeypatch, candidate):
    changed = dict(circuit.REPAIR_SWAPS)
    changed[candidate, 0] = changed[candidate, 0][1:]
    monkeypatch.setattr(circuit, "REPAIR_SWAPS", changed)
    with pytest.raises(AssertionError, match="floor/ceiling"):
        certificate.channel_certificate(candidate)


def test_full_channel_keeps_current_coherence_but_records_contrast():
    # The fixed X,Z witness occupies addresses 0,1. For coin zero, C
    # leaves their shared current record blank, while E records contrast
    # -1 for Z. Consequently the off-diagonal survives only under C.
    images = {}
    for candidate in ("C", "E"):
        columns = []
        for address in (0, 1):
            column = []
            for coin in (0, 1):
                output = certificate._permutation_image(candidate, coin, address)
                column.append(((coin, output // 8), output % 8))
            columns.append(column)
        images[candidate] = certificate._matrix_unit_image(*columns)
    assert images["C"] == {((0, 0), 0, 1): F(1, 2)}
    assert images["E"] == {}


def test_readout_keeps_the_complete_minus_space_and_reference_coherences():
    proof = certificate.readout_certificate()
    assert proof["native_matrix_units_checked"] == 78 * 78
    assert proof["encoded_basis_columns_checked"] == 512
    assert proof["native_minus_rank"] == 77
    assert proof["encoded_minus_rank"] == 511
    # These two distinct basis states are both in the minus subspace. Their
    # off-diagonal matrix unit must survive; a source-basis readout loses it.
    codes = [code for code in circuit.SOURCE_CODES if code not in (0, 1)]
    left, right = codes[:2]
    a = certificate._encoded_readout_column(left, -1)
    b = certificate._encoded_readout_column(right, -1)
    assert {(x, y): u * v for x, u in a.items() for y, v in b.items()} == {(left, right): F(1)}


def test_trace_budget_counts_complete_services_and_retains_frozen_allowance():
    resources = {"native_pulses": 7, "proper_code_load_events": 7,
                 "proper_code_export_events": 7}
    terms = {name: F(1, 64 * 30) for name in certificate.ERROR_TERMS}
    # Three per-service terms contribute 21 units, and nine aggregate
    # channel/numerical terms contribute another nine: exactly eta.
    proof = certificate.apparatus_error_contract(resources, terms)
    assert F(proof["additive_error_upper"]) == F(1, 64)
    assert proof["sufficient_for_frozen_protocol_if_premises_hold"]
    assert protocol.SERVICE_ERROR == certificate.APPARATUS_ALLOWANCE == F(1, 64)
    assert not proof["physical_error_certificate"]
    assert not proof["certified_numerical_execution"]
    terms["diagonal_randomizer"] += F(1, 10**9)
    assert not certificate.apparatus_error_contract(resources, terms)["sufficient_for_frozen_protocol_if_premises_hold"]


def test_complete_instrument_bound_controls_aggregate_effects_not_sum_of_bins():
    resources = {"native_pulses": 0, "proper_code_load_events": 0,
                 "proper_code_export_events": 0}
    terms = {name: F() for name in certificate.ERROR_TERMS}
    terms["record"] = F(1, 128)
    terms["readout"] = F(1, 128)
    proof = certificate.apparatus_error_contract(resources, terms)
    assert proof["every_joint_plus_and_abort_probability_error_upper_if_premises_hold"] == "1/64"
    terms["compilation"] = F(2)
    proof = certificate.apparatus_error_contract(resources, terms)
    assert proof["complete_instrument_trace_distance_upper_if_premises_hold"] == "1"
    assert not proof["sufficient_for_frozen_protocol_if_premises_hold"]


@pytest.mark.parametrize("change", ("missing", "float", "negative", "bool_count"))
def test_error_contract_requires_explicit_exact_nonnegative_premises(change):
    resources = {name: 1 for name in certificate.RESOURCE_TERMS}
    terms = {name: F() for name in certificate.ERROR_TERMS}
    if change == "missing":
        terms.pop("numerical_execution")
    elif change == "float":
        terms["pulse_operator"] = 0.
    elif change == "negative":
        terms["coin"] = F(-1)
    else:
        resources["native_pulses"] = True
    with pytest.raises(ValueError):
        certificate.apparatus_error_contract(resources, terms)
