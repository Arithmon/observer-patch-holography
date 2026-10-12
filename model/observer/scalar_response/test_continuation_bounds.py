"""Independent arithmetic and custody checks of the two-read certificate."""
from copy import deepcopy
from fractions import Fraction as Q
from itertools import product
import json
from math import isqrt
from pathlib import Path

import pytest

from . import continuation_bounds as subject


@pytest.fixture(scope="module")
def receipt():
    return subject.run()


def rational_sqrt_upper(x):
    """Integer arithmetic independent of the parent's interval package."""
    scale = 10 ** 30
    root = isqrt(x.numerator * scale * scale // x.denominator)
    if Q(root * root, scale * scale) < x:
        root += 1
    return Q(root, scale)


def test_walsh_total_variation_identity_exhaustively():
    # Directly sum four signed atom differences. These are exact rationals,
    # and this identity does not assume that arbitrary coefficients are laws.
    for a, b, c in product((Q(-2), Q(-1, 3), Q(0), Q(1, 4), Q(3)), repeat=3):
        direct = sum(abs(y * a + z * b + y * z * c)
                     for y, z in product((-1, 1), repeat=2)) / 8
        assert subject.walsh_tv_upper(abs(a), abs(b), abs(c)) == direct
    for bad in (-1, .1, True):
        with pytest.raises(ValueError):
            subject.walsh_tv_upper(bad, Q(0), Q(0))


def test_independent_rational_kernel_error_reconstruction(receipt):
    parent = json.loads((subject.ROOT / subject.PARENT).read_text())
    b = Q(parent["continuum"]["coherent_mode_norm"][1])
    tail = rational_sqrt_upper(Q(parent["continuum"]["smooth_comparison_tail_mass"][1]) / 2)
    root2 = rational_sqrt_upper(Q(2))
    # 1/root2 is a lower approximation; adding one final lattice unit gives
    # an upper reciprocal sqrt and avoids relying on binary floating sqrt.
    invroot2 = 1 / root2 + Q(1, 10 ** 30)
    for source, result in zip(parent["levels"], receipt["levels"]):
        db = (Q(source["full_preparation_frequency_residual_upper"][1]) / 2
              + Q(source["sampling_isometry_error_upper"][1])) * invroot2
        eta = Q(source["full_coherent_frequency_residual_upper"][1]) * invroot2
        eg = db + rational_sqrt_upper(Q(source["sampled_detector_tail_mass"][1]) / 2)
        graph_v = eg * (2 * b + eg)
        graph_g = b * (2 * eg + eta / 20) + eg * (eg + eta / 20)
        continuum = tail * (2 * b + tail)
        # The independent 30-place square roots and producer's 50-place
        # roots agree far inside its outward 12-place exported intervals.
        for key, expected in (("variance_difference_upper", graph_v + continuum),
                              ("cross_kernel_difference_upper", graph_g + continuum)):
            reported = Q(result[key])
            assert expected - Q(1, 10 ** 20) <= reported <= expected + Q(1, 10 ** 12)


def test_complete_history_bound_includes_zero_current_branch(receipt):
    for row in receipt["levels"]:
        zero = Q(row["vacuum_conditional_history_tv_upper"])
        nonzero = Q(row["nonzero_current_conditional_history_tv_upper"])
        complete = Q(row["complete_history_tv_upper"])
        # All three exported upper bounds were separately rounded outward;
        # compare with an allowance for those rational export roundings.
        independent = Q(28, 33) * zero + Q(5, 33) * nonzero
        assert abs(complete - independent) <= Q(1, 10 ** 12)
        assert 0 < zero <= 1 and 0 < complete <= 1
    assert Q(receipt["levels"][2]["complete_history_tv_upper"]) < Q(522, 10000)
    assert Q(receipt["levels"][3]["complete_history_tv_upper"]) < Q(2571, 100000)
    assert receipt["settings"] == {
        "times": ["19/20", "1"], "gamma": "4", "detector_couplings": ["1", "1"],
        "field_reset_between_reads": False,
    }


def test_certified_noncommutation_and_two_distinct_witnesses(receipt):
    for witnesses in [receipt["continuum_witnesses"],
                      *(row["finite_field_witnesses"] for row in receipt["levels"][2:])]:
        lo, hi = map(Q, witnesses["commutator_over_i_interval"])
        assert 0 < lo < hi < 1
        assert Q(witnesses["vacuum_bit_correlation_lower"]) > Q(35, 1000)
        assert Q(witnesses["second_plus_disturbance_given_positive_current_lower"]) > 0
    finest = receipt["levels"][-1]["finite_field_witnesses"]
    assert Q(finest["vacuum_bit_correlation_lower"]) > Q(456, 10000)
    assert Q(finest["second_plus_disturbance_given_positive_current_lower"]) > Q(2228, 10 ** 8)


def test_coarse_witnesses_use_signed_fallback_when_sign_is_unresolved(receipt):
    parent = json.loads((subject.ROOT / subject.PARENT).read_text())
    # These coarse enclosures admit a positive definite two-direction
    # covariance with C<0, and a positive unmeasured response r2. Hence
    # exp(-V)*sinh(C)<0 and -r2*(1-cos(c/2))<0: zero is not a lower bound
    # that follows from the available premises. This does not substitute
    # the example for the fixed field; it independently checks the missing
    # sign information that requires the universal signed fallback.
    variance, cross, commutator, response = Q(1, 10), Q(-1, 20), Q(1, 20), Q(1, 10)
    assert variance**2 > cross**2 + commutator**2 / 4
    assert cross < 0 and 0 < commutator < 1 and response > 0
    for source, row in zip(parent["levels"][:2], receipt["levels"][:2]):
        witnesses = row["finite_field_witnesses"]
        for key, value in (("variance_interval", variance),
                           ("real_cross_covariance_interval", cross),
                           ("commutator_over_i_interval", commutator)):
            low, high = map(Q, witnesses[key])
            assert low <= value <= high
        low, high = map(Q, source["graph_compact_induced_response"])
        assert low <= response <= high
        assert Q(witnesses["vacuum_bit_correlation_lower"]) == -1
        assert Q(witnesses["second_plus_disturbance_given_positive_current_lower"]) == -1


@pytest.mark.parametrize("mutation", ("time", "reset", "zero_current", "kernel", "disturbance", "parent", "type"))
def test_replay_rejects_resealed_scientific_mutations(receipt, mutation):
    forged = deepcopy(receipt)
    if mutation == "time":
        forged["settings"]["times"][0] = "39/40"
    elif mutation == "reset":
        forged["settings"]["field_reset_between_reads"] = True
    elif mutation == "zero_current":
        forged["native_first_current_probabilities"]["0"] = "0"
    elif mutation == "kernel":
        forged["levels"][-1]["cross_kernel_difference_upper"] = "0"
    elif mutation == "disturbance":
        forged["levels"][-1]["finite_field_witnesses"]["second_plus_disturbance_given_positive_current_lower"] = "1"
    elif mutation == "parent":
        forged["parent"]["sha256"] = "0" * 64
    else:
        forged["levels"][-1]["q"] = 233.0
    with pytest.raises(ValueError, match="replay mismatch"):
        subject.verify(forged)


@pytest.mark.parametrize("relative", (subject.PARENT, subject.NATIVE,
                                     "code/source_scalar_packet/source_scalar_intervals.py"))
def test_parent_native_and_arithmetic_sources_are_bound(monkeypatch, relative):
    original = Path.read_bytes
    target = subject.ROOT / relative

    def changed(path):
        data = original(path)
        return data + b"\n# injected drift\n" if path == target else data

    monkeypatch.setattr(Path, "read_bytes", changed)
    with pytest.raises(ValueError, match="pin mismatch"):
        subject.run()


def test_no_arbitrary_parent_or_settings_input_and_exact_replay(receipt):
    with pytest.raises(TypeError):
        subject.run({"times": [0, 1]})
    assert subject.verify(receipt)
