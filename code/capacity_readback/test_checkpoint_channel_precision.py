"""Exact channel controls for capacity decisions, independent of decoder search."""
from fractions import Fraction

import pytest

from correctable_public_record_capacity import (
    approximate_public_capacity, compound_confusability_graph,
    evaluate_terminal, support_relation_semigroup, tv_robustness_bound,
)
from test_correctable_public_record_capacity import binary_packet


def binary_channel(error):
    return {"rows": {"a": {"a": 1-error, "b": error},
                     "b": {"a": error, "b": 1-error}}}


@pytest.mark.parametrize("error", [Fraction(1, 2**40), Fraction(1, 2**60),
                                  Fraction(1, 10**400)])
def test_every_positive_binary_error_collapses_zero_error_capacity(error):
    channels = [binary_channel(error)]
    assert compound_confusability_graph(["a", "b"], channels) == {"a": {"b"}, "b": {"a"}}
    # Every output is possible for both inputs, so no two-word zero-error code exists.
    assert approximate_public_capacity(["a", "b"], channels, 0)["capacity"] == 1


def test_error_threshold_is_not_enlarged_by_an_acceptance_tolerance():
    error = Fraction(1, 8)
    channels = [binary_channel(error)]
    # The identity decoder's worst error is exactly 1/8. The swapped decoder
    # has error 7/8; either constant decoder has worst error one.
    assert approximate_public_capacity(["a", "b"], channels,
                                       error-Fraction(1, 10**14))["capacity"] == 1
    assert approximate_public_capacity(["a", "b"], channels, error)["capacity"] == 2


@pytest.mark.parametrize("value", [float("nan"), True, "1.0"])
def test_nonprobability_values_do_not_manufacture_a_support_graph(value):
    with pytest.raises(ValueError):
        compound_confusability_graph(["a", "b"], [
            {"rows": {"a": {"shared": value}, "b": {"shared": value}}}])


def test_nan_channel_cannot_certify_saturated_public_capacity():
    packet = binary_packet()
    rows = packet["global_checkpoint_kernels"][0]["rows"]
    for key in rows:
        rows[key] = {"same": float("nan")}
    assert evaluate_terminal(packet)["status"] == "NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING"


def test_duplicate_records_are_not_two_codewords():
    with pytest.raises(ValueError):
        approximate_public_capacity(["a", "a"], [{"rows": {"a": {"x": 1}}}], 0)


def test_output_labels_are_not_silently_stringified():
    with pytest.raises(ValueError):
        compound_confusability_graph(["a", "b"], [
            {"rows": {"a": {0: 1}, "b": {"0": 1}}}])


def test_indefinite_continuation_cannot_discard_positive_escape_mass():
    channel = {"rows": {"a": {"a": Fraction(1, 2), "escape": Fraction(1, 2)},
                        "b": {"b": Fraction(1, 2), "escape": Fraction(1, 2)}}}
    assert compound_confusability_graph(["a", "b"], [channel]) == {"a": {"b"}, "b": {"a"}}
    with pytest.raises(ValueError):
        support_relation_semigroup(["a", "b"], [channel])


@pytest.mark.parametrize(("epsilon", "delta"), [(.3, .6), (.1, 1e-20)])
def test_tv_error_bound_never_rounds_down(epsilon, delta):
    assert Fraction(tv_robustness_bound(epsilon, delta)) >= Fraction(epsilon)+Fraction(delta)
