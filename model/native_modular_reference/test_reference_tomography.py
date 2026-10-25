"""Native effects, target-blind inversion, continued sampling and audit guards."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import reference_tomography as r
from ..native_interface import construction as interface


def _family_counts(mean=F(2, 5), pair=F(1, 20), success=F(7, 13)):
    weights = r.variation.family(mean, pair)
    basis = r.native.states(binary=True)
    probabilities = {label: F() for label in r.LABELS}
    # Independent native source branches; no production current law or inverse.
    for x, w in zip(basis, weights, strict=True):
        for coin in (0, 1):
            _, (_, current) = r.native.native(x, r.SEAM, coin)
            probabilities[f"success_{current}"] += success*w/2
    probabilities["herald_failure"] = 1-success
    unit = math.lcm(*(p.denominator for p in probabilities.values()))
    scale = max(1, math.ceil((1 << 27)/unit))
    return {label: int(p*unit)*scale for label, p in probabilities.items()}, weights, success


def test_actual_native_effects_and_nonidentifiability_without_completion():
    evidence = r.exact_identification()
    assert evidence["charge_sector_current_probabilities"] == [
        (F(1, 10), F(9, 10), F()), (F(1, 40), F(3, 4), F(9, 40)),
        (F(), F(1, 2), F(1, 2))]
    assert evidence["uncompleted_countercontrol"]["same_current_law"] == (F(), F(1, 2), F(1, 2))
    core, reads = r.population_counts()
    for law in ("C", "E"):
        branches = interface.fixed_maps(law, r.SEAM)
        rebuilt = {j: np.zeros((66, 66), dtype=np.int64) for j in (-1, 0, 1)}
        for branch in branches:
            for before, after in branch.mapping:
                rebuilt[branch.current][after, before] += 1
        for current in rebuilt:
            np.testing.assert_array_equal(reads[current], rebuilt[current])
    assert np.all(core.sum(axis=0) == 2)


def test_cached_native_counts_cannot_be_poisoned():
    core, reads = r.population_counts()
    original = core.copy(), {j: x.copy() for j, x in reads.items()}
    core[:] = 100
    reads[0][:] = 100
    reads.clear()
    fresh_core, fresh_reads = r.population_counts()
    np.testing.assert_array_equal(fresh_core, original[0])
    assert set(fresh_reads) == {-1, 0, 1}
    for j in fresh_reads:
        np.testing.assert_array_equal(fresh_reads[j], original[1][j])


def test_record_contract_matches_existing_preparation_without_reference_lookup(monkeypatch):
    original = r.local.budget()
    r._record_contract.cache_clear()
    monkeypatch.setattr(r.local, "reference", lambda: (_ for _ in ()).throw(AssertionError("target lookup")))
    error, q, cost = r._record_contract()
    assert error == original["complete_flagged_reference_halftrace_upper"]
    assert q == original["uniform_ideal_herald_probability"]
    assert cost == original["total_native_attempts_upper"]


def test_inverse_uses_unknown_mass_and_unfitted_family(monkeypatch):
    counts, weights, success = _family_counts()
    monkeypatch.setattr(r.local, "reference", lambda: (_ for _ in ()).throw(AssertionError("target lookup")))
    inferred = r.infer_counts(counts)
    assert inferred["status"] == "faithful"
    assert inferred["observed_success_mass"] == success
    assert inferred["moments"] == (F(2, 5), F(1, 20))
    by_charge = tuple(next(w for w, q in zip(weights, r.variation.charges(), strict=True) if q == i)
                      for i in range(3))
    assert inferred["central_subnormalized_weights"] == tuple(success*w for w in by_charge)
    assert inferred["central_weights"] == by_charge
    assert inferred["confidence_failure_upper"] == F(3, 512)
    # Point inversion works for other mass/family data; the stated statistical
    # guarantee still concerns the source preparer's complete record contract.


def test_bernstein_radius_and_log_box_are_conservative():
    counts, _, _ = _family_counts()
    row = r.infer_counts(counts)
    n, k = row["trials"], row["confidence_bits"]
    error = row["systematic_halftrace_upper"]
    for i, ((lo, hi), center) in enumerate(zip(row["subnormalized_weight_intervals"], row["central_subnormalized_weights"], strict=True)):
        radius = (hi-lo)/2-r.WIDTHS[i]*error
        assert n*radius**2 >= 2*k*(row["score_conditional_second_moment_upper"][i]+r.WIDTHS[i]*radius/3)
        assert lo <= center <= hi
        assert float(row["logarithm_operator_error_upper"]) >= max(abs(math.log(float(lo/center))), abs(math.log(float(hi/center))))


@pytest.mark.parametrize("mutate", [
    lambda c: c.pop("herald_failure"),
    lambda c: c.update(extra=0),
    lambda c: c.update(success_0=True),
    lambda c: c.update(success_0=0.0),
    lambda c: c.update(success_0=-1),
    lambda c: c.update({key: 0 for key in c}),
])
def test_invalid_or_incomplete_records_fail_closed(mutate):
    counts, _, _ = _family_counts()
    mutate(counts)
    with pytest.raises(ValueError):
        r.infer_counts(counts)


@pytest.mark.parametrize("bias", [True, -1, 2, 0.0, float("nan"), "0"])
def test_unproved_invalid_numeric_bias_forms_are_rejected(bias):
    counts, _, _ = _family_counts()
    with pytest.raises(ValueError):
        r.infer_counts(counts, sampling_bias=bias)


def test_tiny_sample_and_all_failures_are_unresolved():
    counts = {key: 0 for key in r.LABELS}
    counts["success_1"] = 1
    assert r.infer_counts(counts)["status"] == "unresolved"
    counts["success_1"] = 0
    counts["herald_failure"] = 1 << 25
    row = r.infer_counts(counts)
    assert row["status"] == "unresolved"
    assert "central_weights" not in row
    assert "logarithm_operator_error_upper" not in row


def test_residual_sampler_reuses_terminal_state_instead_of_resetting():
    kernel = np.zeros((462, 66), dtype=np.uint64)
    for state in range(66):
        kernel[2*66+(state+1) % 66, state] = r.GRID
    sample = r.sample_kernel(kernel, 137, initial=17, chunk_size=13)
    assert sample["residual_draws"] == 137
    assert sample["final_state"] == (17+137) % 66
    assert sample["prefix"] == [("success_1", (17+i+1) % 66) for i in range(32)]
    assert sample["counts"]["success_1"] == 137
    assert sum(map(sum, sample["terminal_state_counts"])) == 137


def test_common_and_residual_samples_carry_joint_record_state_across_chunks():
    kernel = np.zeros((462, 66), dtype=np.uint64)
    for state in range(66):
        kernel[66+3, state] = r.GRID//2
        kernel[2*66+(state+1) % 66, state] += r.GRID//2
    a = r.sample_kernel(kernel, 1000, seed=77, initial=9, chunk_size=1)
    b = r.sample_kernel(kernel, 1000, seed=77, initial=9, chunk_size=137)
    assert a == b
    assert 400 < a["residual_draws"] < 600
    assert a["common_subkernel_probability"] == F(1, 2)


@pytest.fixture(scope="module")
def receipt():
    return r._json(r.run())


def test_executed_source_kernel_counts_certify_a_nontrivial_reference(receipt):
    row = receipt["continuing_execution"]
    infer, sample, audit = row["inference"], row["execution"], row["sampler_certificate"]
    assert sum(sample["counts"].values()) == 1 << 24
    assert sample["same_carrier_reused"]
    assert not sample["microscopic_attempt_tapes_executed"]
    assert not sample["seeded_PRNG_demonstration_proves_independent_uniform_randomness"]
    assert infer["status"] == "faithful"
    assert F(infer["logarithm_operator_error_upper"]) < F(17, 50)
    assert audit == r._json(r._portable_sampler_certificate())
    assert sample["common_subkernel_probability_lower"] == str(1-r.SAMPLER_GATE)
    assert "common_subkernel_probability" not in sample
    assert row["validation_only_known_target_in_confidence_box"]
    assert infer["source_native_attempts_all_trials_upper"] == (1 << 24)*(r._record_contract()[2]+2)
    assert r.compare_receipt(receipt, receipt)


def test_all_original_population_columns_and_both_laws_are_used():
    c, diagnostic = r.quantized_kernel("C")
    e, _ = r.quantized_kernel("E")
    np.testing.assert_array_equal(c, e)
    assert c.shape == (7*66, 66)
    assert diagnostic["binary_population_mass_residual_diagnostic"] < 5e-10
    audit = r.audit_kernel(c)
    assert audit["exact_dyadic_columns_checked"] == 66
    assert audit["complete_total_variation_upper"] < r.SAMPLER_GATE
    caps = r._portable_sampler_certificate()
    for field in ("score_mean_bias_upper", "score_second_moment_upper"):
        assert all(actual <= bound for actual, bound in zip(audit[field], caps[field], strict=True))


@pytest.mark.parametrize("field", ["complete_total_variation_upper", "score_mean_bias_upper",
                                  "score_second_moment_upper", "common_subkernel_probability"])
def test_export_requires_exact_sampler_audit_to_pass_every_portable_cap(monkeypatch, field):
    audit = r._portable_sampler_certificate()
    sample = {"common_subkernel_probability": F(1)}
    if field == "common_subkernel_probability":
        sample[field] = 1-2*r.SAMPLER_GATE
    elif field == "complete_total_variation_upper":
        audit[field] *= 2
    else:
        audit[field] = (2*audit[field][0], *audit[field][1:])
    monkeypatch.setattr(r, "quantized_kernel", lambda law: (None, {}))
    monkeypatch.setattr(r, "audit_kernel", lambda kernel: audit)
    monkeypatch.setattr(r, "sample_kernel", lambda *args, **kwargs: sample)
    with pytest.raises(ArithmeticError, match="portable source-comparison gates"):
        r.execute_trials(1)


@pytest.mark.parametrize("field, value", [
    ("complete_total_variation_upper", "0"),
    ("complete_total_variation_upper", "-1/100000000000000000000"),
    ("complete_total_variation_upper", "1/100"),
    ("score_mean_bias_upper", ["0", "0", "0"]),
    ("score_mean_bias_upper", ["-1/100000000000000000000", "0", "0"]),
    ("score_second_moment_upper", ["0", "0", "0"]),
    ("score_second_moment_upper", ["0", "0", "-1/100000000000000000000"]),
])
def test_receipt_never_tolerates_false_or_changed_sampler_caps(receipt, field, value):
    changed = deepcopy(receipt)
    changed["continuing_execution"]["sampler_certificate"][field] = value
    assert not r.compare_receipt(changed, receipt)
    assert not r.compare_receipt(receipt, changed)
    assert not r.compare_receipt(changed, changed)


def test_receipt_rejects_false_perfect_sampler_and_only_tolerates_float_mass_diagnostic(receipt):
    changed = deepcopy(receipt)
    row = changed["continuing_execution"]
    row["sampler_certificate"]["complete_total_variation_upper"] = "0"
    row["sampler_certificate"]["score_mean_bias_upper"] = ["0"]*3
    row["sampler_certificate"]["score_second_moment_upper"] = ["0"]*3
    row["execution"]["common_subkernel_probability_lower"] = "1"
    row["numerical_kernel"]["binary_population_mass_residual_diagnostic"] = 0.0
    assert not r.compare_receipt(changed, receipt)
    assert not r.compare_receipt(receipt, changed)
    changed = deepcopy(receipt)
    changed["continuing_execution"]["execution"]["common_subkernel_probability_lower"] = "1"
    assert not r.compare_receipt(changed, receipt)
    assert not r.compare_receipt(changed, changed)
    changed = deepcopy(receipt)
    changed["continuing_execution"]["numerical_kernel"]["binary_population_mass_residual_diagnostic"] = 0.0
    assert r.compare_receipt(changed, receipt)
    assert r.compare_receipt(receipt, changed)


def test_receipt_keeps_counts_and_source_hashes_exact(receipt):
    changed = deepcopy(receipt)
    changed["continuing_execution"]["execution"]["counts"]["herald_failure"] += 1
    assert not r.compare_receipt(changed, receipt)
    changed = deepcopy(receipt)
    first = next(iter(changed["source_sha256"]))
    changed["source_sha256"][first] = "0"*64
    assert not r.compare_receipt(changed, receipt)


@pytest.mark.parametrize("mutation", [
    lambda row: row["execution"]["counts"].update(warmup_timeout=False),
    lambda row: row["execution"]["counts"].update(
        herald_failure=float(row["execution"]["counts"]["herald_failure"])),
    lambda row: row["execution"].update(trials=float(row["execution"]["trials"])),
    lambda row: row["execution"].update(same_carrier_reused=1),
    lambda row: row["numerical_kernel"].update(binary_population_mass_residual_diagnostic=False),
    lambda row: row["sampler_certificate"].update(
        complete_total_variation_upper="+"+row["sampler_certificate"]["complete_total_variation_upper"]),
    lambda row: row["sampler_certificate"].update(
        complete_total_variation_upper=" "+row["sampler_certificate"]["complete_total_variation_upper"]),
    lambda row: row["sampler_certificate"].update(complete_total_variation_upper="{}/{}".format(
        2*F(row["sampler_certificate"]["complete_total_variation_upper"]).numerator,
        2*F(row["sampler_certificate"]["complete_total_variation_upper"]).denominator)),
])
def test_receipt_exact_types_and_canonical_rationals_cannot_alias(receipt, mutation):
    changed = deepcopy(receipt)
    mutation(changed["continuing_execution"])
    assert not r.compare_receipt(changed, receipt)
    assert not r.compare_receipt(receipt, changed)
