"""Finite native-record phase inversion and exact confidence countercontrols."""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
import pytest

from . import phase_inference as p

WEIGHTS = (F(2, 225), F(1, 36), F(2, 45))
N = 1 << 22


def _training(weights=WEIGHTS):
    """Independent exact native current law, used only as a validation fixture."""
    r = p.tomography
    q = r._record_contract()[1]
    probabilities = dict.fromkeys(r.LABELS, F())
    for state, charge in zip(r.native.states(binary=True), r.variation.charges(), strict=True):
        for coin in (0, 1):
            _, (_, current) = r.native.native(state, r.SEAM, coin)
            probabilities[f"success_{current}"] += q*weights[charge]/2
    probabilities["herald_failure"] = 1-q
    denominator = math.lcm(*(x.denominator for x in probabilities.values()))
    multiplier = math.ceil((1 << 24)/denominator)
    return {key: int(value*denominator)*multiplier for key, value in probabilities.items()}


def _source_allowance():
    return p.readout.complete_trial_budget()["complete_trial_halftrace_upper"]+F(1, 10**10)


@lru_cache(maxsize=16)
def _record_fixture(s=F(3, 4), *, law="C", index=0, weights=WEIGHTS, trials=N):
    """Rounded complete ideal-instrument frequencies, not an execution receipt."""
    counts = {}
    q = p.tomography._record_contract()[1]
    density = np.array([float(q*weights[charge]) for charge in p.tomography.variation.charges()])
    for quadrature in p.QUADRATURES:
        keys, kernel = p.readout.population_instrument(weights, s, law=law, index=index, quadrature=quadrature)
        probabilities = (kernel@density).reshape(len(keys), 66).sum(axis=1)
        column = dict.fromkeys(p.complete_alphabet(law, index), 0)
        for key, probability in zip(keys, probabilities, strict=True):
            column[key] = int(probability*trials)
        column[("preparation", "herald_failure")] = trials-sum(column.values())
        counts[quadrature] = column
    return counts


def _records(*args, **kwargs):
    return deepcopy(_record_fixture(*args, **kwargs))


def _infer(counts, training=None, **kwargs):
    return p.infer_counts(counts, _training() if training is None else training,
                         systematic_halftrace_upper=_source_allowance(), **kwargs)


@pytest.mark.parametrize("law,index,s", [
    ("C", 0, F(3, 4)), ("C", 1, F(3, 4)),
    ("E", 0, -F(3, 4)), ("E", 1, F()),
])
def test_complete_native_frequencies_identify_blind_signed_parameter(law, index, s):
    counts = _records(s, law=law, index=index)
    result = _infer(counts, law=law, index=index)
    assert result["status"] == "identified"
    low, high = result["parameter_interval"]
    assert low <= s <= high
    assert high-low < F(1, 5)
    assert result["confidence_failure_upper"] == F(7, 1024)
    assert result["orientation"] == ("positive" if s > 0 else "negative" if s < 0 else "unresolved")
    assert result["unknown_amplitude_fitted"] is False
    assert result["successful_trials_conditioned_on"] is False
    for quadrature, row in result["observed_scores"].items():
        assert row["trials"] == N
        assert len(row["complete_counts"]) == len(counts[quadrature])
        assert row["score_sum"] == sum(p.readout.score(key)*value for key, value in counts[quadrature].items()
                                       if key != ("numerical_abort",))
        assert row["score_square_sum"] == sum(p.readout.score(key)**2*value for key, value in counts[quadrature].items()
                                              if key != ("numerical_abort",))
        assert counts[quadrature][("preparation", "herald_failure")] > N*4//5
    reference_units = result["reference_unit_parameter_interval"]
    assert reference_units[0] <= s <= reference_units[1]


def test_decoder_never_reads_target_parameter_or_reference(monkeypatch):
    counts, training = _records(), _training()
    def forbidden(*args, **kwargs):
        raise AssertionError("the blind decoder looked up a validation target")
    monkeypatch.setattr(p.readout.local, "reference", forbidden)
    monkeypatch.setattr(p.readout, "response_enclosure", forbidden)
    monkeypatch.setattr(p.readout, "stage", forbidden)
    monkeypatch.setattr(p.readout.learned, "clock", forbidden)
    assert _infer(counts, training)["status"] == "identified"
    with pytest.raises(TypeError):
        _infer(counts, training, s=F(3, 4))


@pytest.mark.parametrize("mutate", [
    lambda c: c.pop("sin"),
    lambda c: c.update(other={}),
    lambda c: c["cos"].pop(("numerical_abort",)),
    lambda c: c["sin"].pop(("abort", 0)),
    lambda c: c["cos"].pop(("preparation", "completion_timeout")),
    lambda c: c["cos"].update({("extra",): 0}),
    lambda c: c["cos"].update({("numerical_abort",): True}),
    lambda c: c["cos"].update({("numerical_abort",): 0.0}),
    lambda c: c["cos"].update({("numerical_abort",): -1}),
    lambda c: c["cos"].update({("numerical_abort",): 1 << 41}),
    lambda c: c.update(cos=dict.fromkeys(c["cos"], 0)),
])
def test_all_complete_branches_original_counts_and_totals_required(mutate):
    counts = _records()
    mutate(counts)
    with pytest.raises(ValueError):
        _infer(counts)


def test_equal_but_nonoriginal_boolean_record_is_rejected():
    counts = _records()
    key = next(key for key in counts["cos"] if key[0] == "read" and key[1] == 1)
    value = counts["cos"].pop(key)
    counts["cos"][(key[0], True)+key[2:]] = value
    with pytest.raises(ValueError):
        _infer(counts)


@pytest.mark.parametrize("law,index", [("other", 0), ("C", True), ("C", -1), ("E", 2)])
def test_source_law_and_cap_are_original_admitted_values(law, index):
    with pytest.raises(ValueError):
        _infer(_records(), law=law, index=index)


@pytest.mark.parametrize("value", [True, 0.0, float("nan"), "1/1000", -1, 2, F()])
def test_systematic_error_cannot_be_missing_floating_or_erased(value):
    with pytest.raises(ValueError):
        p.infer_counts(_records(), _training(), systematic_halftrace_upper=value)
    with pytest.raises(TypeError):
        p.infer_counts(_records(), _training())


@pytest.mark.parametrize("value", [True, 3, 65, F(12), 12.0])
def test_confidence_parameter_is_original_bounded_integer(value):
    with pytest.raises(ValueError):
        _infer(_records(), confidence_bits=value)


@pytest.mark.parametrize("window", [(-1.0, 1), [F(-1), F(1)], (0, 0), (1, -1), (False, 1)])
def test_window_requires_explicit_exact_ordered_endpoints(window):
    with pytest.raises(ValueError):
        _infer(_records(), parameter_window=window)


def test_training_with_no_positive_family_does_not_create_clock():
    training = dict.fromkeys(p.tomography.LABELS, 0)
    training["herald_failure"] = 1 << 24
    result = _infer(_records(), training)
    assert result["status"] == "unresolved"
    assert "reference" in result["reason"]
    assert "parameter_interval" not in result


def test_equal_reference_weights_do_not_create_orientation():
    result = _infer(_records(), _training((F(1, 66),)*3))
    assert result["reference_inference"]["status"] == "faithful"
    assert result["status"] == "unresolved"
    assert "positive quadrature orientation" in result["reason"]


def test_wide_window_is_refused_even_for_well_resolved_data():
    result = _infer(_records(), parameter_window=(F(-8), F(8)))
    assert result["status"] == "unresolved"
    assert "right half-plane" in result["reason"]
    lo, hi = result["score_mean_alias_period_interval"]
    exact_diagnostic = 2*math.pi/math.log(float(WEIGHTS[1]/WEIGHTS[0]))
    assert float(lo)-1e-14 <= exact_diagnostic <= float(hi)+1e-14
    assert "equality of complete instruments is not claimed" in result["alias_scope"]


def test_no_successful_phase_records_do_not_invent_a_parameter():
    counts = {q: dict.fromkeys(p.complete_alphabet(), 0) for q in p.QUADRATURES}
    for row in counts.values():
        row[("preparation", "herald_failure")] = N
    result = _infer(counts)
    assert result["status"] == "unresolved"
    assert "parameter_interval" not in result


def test_negative_cosine_or_inconsistent_amplitude_is_refused():
    counts = _records()
    negative = next(key for key in counts["cos"] if key[0] == "read" and p.readout.score(key) == -2)
    counts["cos"] = dict.fromkeys(counts["cos"], 0)
    counts["cos"][negative] = N
    result = _infer(counts)
    assert result["status"] == "unresolved"
    assert "parameter_interval" not in result


def test_observed_sign_cannot_be_approved_outside_declared_window():
    result = _infer(_records(-F(3, 4)), parameter_window=(F(), F(1)))
    assert result["status"] == "unresolved"
    assert "misses" in result["reason"]


def test_statistics_keep_all_trials_and_exact_martingale_exponent():
    result = _infer(_records())
    budget = result["budget"]
    for row in result["observed_scores"].values():
        n, r, k = row["trials"], row["statistical_radius_upper"], budget["confidence_bits"]
        variance = budget["score_conditional_second_moment_upper"]
        assert n*r*r >= 2*k*(variance+4*r/3)
        assert row["total_mean_radius_upper"] == r+4*budget["systematic_halftrace_upper"]
        assert row["empirical_second_moment"] < variance
        assert n == sum(value for _, value in row["complete_counts"])
    assert result["independent_input_copies_required"] is False


def test_systematic_bias_survives_large_trial_counts():
    budget = p.score_budget(_source_allowance())
    assert p._radius(1 << 40, budget) > 0
    assert budget["score_mean_bias_upper"] > 0
    request = F(1, 1000)
    planned = p.trial_budget(request, systematic_halftrace_upper=_source_allowance())
    assert planned["attained_score_radius_upper"] <= request
    assert planned["trials_per_quadrature"] < 1 << 24
    assert planned["complete_trials"] == 2*planned["trials_per_quadrature"]
    with pytest.raises(ValueError):
        p.trial_budget(budget["score_mean_bias_upper"], systematic_halftrace_upper=_source_allowance())


@pytest.mark.parametrize("x", [F(), F(1, 100), F(1, 2), F(1), F(2), F(101), -F(7, 5)])
def test_atan_rational_enclosures_and_odd_symmetry(x):
    lo, hi = p.atan_interval(x)
    nlo, nhi = p.atan_interval(-x)
    assert lo == -nhi and hi == -nlo
    assert hi-lo <= F(1, 1 << 62)
    assert float(lo)-2e-16 <= math.atan(float(x)) <= float(hi)+2e-16
    assert all(type(value) is F for value in (lo, hi))


def test_atan_reduction_and_signed_interval_division_are_exact():
    pi_lo, pi_hi = p.pi_interval()
    lo, hi = p.atan_interval(F(1))
    assert lo <= pi_lo/4 <= pi_hi/4 <= hi
    a, b = p._small_atan(F(1, 3))
    assert a < b
    assert b-a == F(1, 3)**129/129
    assert p._divide((-F(2), F(1)), (F(2), F(4))) == (-F(1), F(1, 2))
    with pytest.raises(ValueError):
        p._divide((F(), F(1)), (F(), F(1)))


def test_extra_failures_are_retained_instead_of_success_renormalization():
    counts = _records()
    first = _infer(counts)
    altered = deepcopy(counts)
    for column in altered.values():
        column[("preparation", "herald_failure")] += N
    second = _infer(altered)
    for quadrature in p.QUADRATURES:
        assert second["observed_scores"][quadrature]["mean"] == first["observed_scores"][quadrature]["mean"]/2
        assert second["observed_scores"][quadrature]["trials"] == 2*N
    assert second["status"] == "unresolved"
    assert "amplitude" in second["reason"]




def test_unresolved_cosine_is_refused_before_any_angle_branch_is_selected():
    result = _infer(_records(F(7, 5)))
    assert result["status"] == "unresolved"
    assert "positive cosine" in result["reason"]
    assert "phase_interval" not in result


def test_native_zero_and_failed_records_remain_zero_score_without_being_discarded():
    counts = _records()
    row = _infer(counts)["observed_scores"]["cos"]
    zero_count = sum(value for key, value in counts["cos"].items()
                     if key == ("numerical_abort",) or p.readout.score(key) == 0)
    assert zero_count > N*4//5
    assert row["trials"] == N
    assert row["trials"] != N-zero_count
