"""Maintainer audit: compound codes, rational thresholds and parser boundaries."""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
from itertools import product

import pytest

from correctable_public_record_capacity import (
    _channel_rows, approximate_public_capacity, compound_confusability_graph,
    support_relation_semigroup, tv_robustness_bound,
)
from test_checkpoint_decoder_certificates import channel, integer_oracle, read_ratio
from verify_checkpoint_decoder import verify_decoder_witness


def flat_minimax_error(counts, code):
    """Independent integer-count minimization, with no production reductions."""
    return min(max(Fraction(sum(n for y, n in enumerate(counts[x]) if decoder[y] != x),
                            sum(counts[x])) for x in code)
               for decoder in product((-1, *code), repeat=len(counts[0])))


WEIGHTED_ROWS = ((2, 1, 0, 0), (0, 2, 1, 0), (1, 0, 0, 4),
                 (2, 3, 5, 7), (0, 0, 0, 1))


@pytest.mark.parametrize("rows", tuple(product(WEIGHTED_ROWS, repeat=3)))
def test_unequal_denominators_at_both_sides_of_compound_optimum(rows):
    family = [rows, rows[1:]+rows[:1]]
    kernels = [channel(counts) for counts in family]
    records = ["0", "1", "2"]
    threshold = max(flat_minimax_error(counts, (0, 1, 2)) for counts in family)
    step = Fraction(1, 10**30)
    tolerances = {threshold, max(Fraction(0), threshold-step), min(Fraction(1), threshold+step)}
    for epsilon in sorted(tolerances):
        result = approximate_public_capacity(records, kernels, epsilon)
        assert result["capacity"] == integer_oracle(family, epsilon)
        code = tuple(map(int, result["code_witness"]))
        for counts, certificate in zip(family, result["decoder_certificates"]):
            # Compare the actual optimum, not just which side of epsilon it lies on.
            assert read_ratio(certificate["worst_input_error"]) == flat_minimax_error(counts, code)
        assert verify_decoder_witness(records, kernels, epsilon, result)


def test_compound_capacity_requires_one_common_code_not_individual_maxima():
    records = ["0", "1", "2"]
    # Each channel merges one pair. All three pairs occur across the family.
    kernels = [channel(rows) for rows in [((1, 0), (1, 0), (0, 1)),
                                          ((0, 1), (1, 0), (1, 0)),
                                          ((1, 0), (0, 1), (1, 0))]]
    assert [approximate_public_capacity(records, [k], 0)["capacity"] for k in kernels] == [2, 2, 2]
    for epsilon in [0, Fraction(1, 2), 1-Fraction(1, 10**400)]:
        result = approximate_public_capacity(records, kernels, epsilon)
        assert result["capacity"] == 1
        assert verify_decoder_witness(records, kernels, epsilon, result)
    assert approximate_public_capacity(records, kernels, 1)["capacity"] == 3


@pytest.mark.parametrize("sign", [-1, 1])
def test_exact_row_mass_acceptance_boundary_and_disclosure(sign):
    tolerance = Fraction(1, 10**12)
    mass = 1+sign*tolerance
    records = ["0", "1"]
    kernels = [{"rows": {"0": {"0": 3*mass/4, "1": mass/4},
                          "1": {"0": Fraction(1, 4), "1": Fraction(3, 4)}}}]
    result = approximate_public_capacity(records, kernels, Fraction(1, 4))
    assert result["capacity"] == 2
    certificate = result["decoder_certificates"][0]
    assert read_ratio(certificate["input_row_masses"]["0"]) == mass
    assert read_ratio(certificate["worst_input_error"]) == Fraction(1, 4)
    assert verify_decoder_witness(records, kernels, Fraction(1, 4), result)
    outside = deepcopy(kernels)
    escaped_mass = mass+sign*Fraction(1, 10**30)
    outside[0]["rows"]["0"] = {"0": 3*escaped_mass/4, "1": escaped_mass/4}
    with pytest.raises(ValueError, match="unit mass"):
        approximate_public_capacity(records, outside, Fraction(1, 4))
    assert not verify_decoder_witness(records, outside, Fraction(1, 4), result)


def test_exact_decimal_channels_do_not_inherit_context_rounding():
    from decimal import localcontext
    kernels = [{"rows": {"0": {"0": Decimal("0.9"), "1": Decimal("0.1")},
                          "1": {"0": Decimal("0.1"), "1": Decimal("0.9")}}}]
    with localcontext() as context:
        context.prec = 1
        result = approximate_public_capacity(["0", "1"], kernels, Decimal("0.1"))
        assert result["capacity"] == 2
        assert read_ratio(result["decoder_certificates"][0]["worst_input_error"]) == Fraction(1, 10)
        assert verify_decoder_witness(["0", "1"], kernels, Decimal("0.1"), result)


@pytest.mark.parametrize("bad", [Decimal("NaN"), Decimal("sNaN"), Decimal("Infinity"),
                                 float("inf"), -Fraction(1, 10**400), 1+Fraction(1, 10**400)])
def test_invalid_rows_outside_the_certified_code_are_still_rejected(bad):
    kernels = [channel(((1, 0), (0, 1)))]
    result = approximate_public_capacity(["0"], [{"rows": {"0": kernels[0]["rows"]["0"]}}], 0)
    kernels[0]["rows"]["1"] = {"1": bad}
    with pytest.raises(ValueError):
        compound_confusability_graph(["0", "1"], kernels)
    assert not verify_decoder_witness(["0", "1"], kernels, 0, result)


def test_zero_external_entries_do_not_invalidate_same_alphabet_closure():
    kernel = {"rows": {"0": {"0": 1, "outside": 0}, "1": {"1": 1, "outside": 0}}}
    assert support_relation_semigroup(["0", "1"], [kernel]) == {frozenset({("0", "0"), ("1", "1")})}


def test_legacy_packet_adapter_refuses_lossy_conversion_but_general_evaluator_accepts():
    kernels = [channel(((2, 1), (1, 2)))]
    with pytest.raises(ValueError, match="exactly representable"):
        _channel_rows(kernels[0], ["0", "1"])
    result = approximate_public_capacity(["0", "1"], kernels, Fraction(1, 3))
    assert result["capacity"] == 2
    assert verify_decoder_witness(["0", "1"], kernels, Fraction(1, 3), result)


def test_tv_transfer_preserves_a_decoder_at_the_exact_new_boundary():
    original = [channel(((7, 1), (1, 7)))]
    perturbed = [channel(((6, 2), (2, 6)))]
    epsilon = Fraction(1, 8)
    delta = Fraction(1, 8)
    result = approximate_public_capacity(["0", "1"], original, epsilon)
    assert result["capacity"] == 2
    transferred = approximate_public_capacity(["0", "1"], perturbed, tv_robustness_bound(epsilon, delta))
    assert transferred["capacity"] == 2
    assert transferred["decoder_certificates"][0]["decoder"] == result["decoder_certificates"][0]["decoder"]
    assert verify_decoder_witness(["0", "1"], perturbed, epsilon+delta, transferred)
    assert approximate_public_capacity(["0", "1"], perturbed, epsilon+delta-Fraction(1, 10**30))["capacity"] == 1
