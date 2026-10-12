"""Independent source histories, controller retention and causal controls."""
from copy import deepcopy
from fractions import Fraction as Q
from itertools import combinations, product
import inspect
import json
import math

import pytest

from . import continuation as subject
from .field import ScalarField
from .kernel import WeylKernel


def endpoint(loads, seam, coin):
    a, b = seam
    total = loads[a]+loads[b]
    out = list(loads)
    out[a] = (total+1)//2 if coin else total//2
    out[b] = total-out[a]
    return tuple(out), (coin, loads[a]-out[a])


def independent_paths():
    rows = {}
    for occupied in combinations(range(12), 2):
        initial = tuple(int(port in occupied) for port in range(12))
        for c1, a, c2, b, c3 in product((0, 1), (-1, 1), (0, 1), (-1, 1), (0, 1)):
            first, r1 = endpoint(initial, (0, 1), c1)
            second, r2 = endpoint(first, (0, 5) if a == 1 else (0, 1), c2)
            final, r3 = endpoint(second, (0, 5) if b == 1 else (0, 1), c3)
            rows[initial, (r1, a, r2, b, r3)] = first, second, final
    return rows


def test_all_2112_paths_continue_same_source_and_preserve_the_full_record():
    expected = independent_paths()
    actual = {(row.initial, row.record.prefix): (row.after_first, row.after_second, row.final)
              for row in subject.native_paths()}
    assert len(actual) == 2112
    assert actual == expected
    for row in subject.native_paths():
        assert row.record.halted
        assert subject.native_step(row.final, row.record, 0) == (row.final, row.record)
        assert subject.native_step(row.final, row.record, 1) == (row.final, row.record)
        # Source conservation and unchanged spectator ports are independent
        # of the field probabilities or observation sampling.
        assert sum(row.final) == 2
        assert all(row.final[k] == row.initial[k] for k in range(12) if k not in (0, 1, 5))


@pytest.mark.parametrize("moments", [
    (Q(), Q(), Q(), Q()),
    (Q(-1, 3), Q(-2, 5), Q(1, 20), Q(3, 20)),
    (Q(1, 5), Q(-1, 8), Q(-1, 10), Q(1, 9)),
])
def test_native_consequences_are_exact_for_independent_parity_laws(moments):
    mu1, mu2, k0, k1 = moments
    total = second = third = occupation = Q()
    for (_, word), (_, _, final) in independent_paths().items():
        r1, a, r2, b, r3 = word
        j = r1[1]
        probability = (1+j*mu1*a+j*mu2*b+(k0 if j == 0 else k1)*a*b)/2112
        assert probability > 0
        total += probability
        second += probability*r2[1]
        third += probability*r3[1]
        occupation += probability*final[0]
    assert total == 1
    assert second == 5*mu1/264
    assert third == 5*(mu2-mu1)/352+(k0-k1)/264
    assert occupation == Q(1, 6)-second-third
    expected = subject.exact_statistics()["moment_coefficients"]
    for name, value in (("second_current", second), ("third_current", third), ("final_port0", occupation)):
        coefficients = expected[name]
        assert Q(coefficients["constant"])+sum(Q(coefficients[key])*m for key, m in
                                               zip(("mu1", "mu2", "kappa0", "kappa1"), moments)) == value


def test_reachable_controller_encoding_and_deletion_counts():
    words = {word for _, word in independent_paths()}
    prefixes = {word[:n] for word in words for n in range(6)}
    counts = [sum(len(p) == n for p in prefixes) for n in range(6)]
    assert counts == [1, 4, 8, 24, 48, 124]
    assert len(prefixes) == 209
    assert 2**7 < len(prefixes) <= 2**8
    codes = set()
    for prefix in prefixes:
        record = subject.Record(prefix)
        code = subject.encode_record(record)
        codes.add(code)
        assert subject.decode_record(code) == record
        checkpoint = json.loads(json.dumps(record.checkpoint()))
        assert set(checkpoint) == {"phase", "prefix"}
        assert subject.Record.restore(checkpoint) == record
    assert codes == set(range(209))
    exact = subject.exact_statistics()
    assert exact["source_times_record_dimension_excluding_field"] == 85722
    assert exact["coarse_terminal_record_counts"] == {
        "native_letters_only": 44, "latest_detector_and_native_letter": 8,
        "detector_pair": 4, "delete_first_detector_only": 78, "complete_erasure": 1}


def test_actuator_access_and_invalid_controller_states():
    assert tuple(inspect.signature(subject.selected_seam).parameters) == ("record",)
    for prefix in subject.record_codes():
        record = subject.Record(prefix)
        if len(prefix) in (2, 4):
            assert subject.selected_seam(record) == ((0, 5) if prefix[-1] == 1 else (0, 1))
        elif len(prefix) == 0:
            assert subject.selected_seam(record) == (0, 1)
        else:
            assert subject.selected_seam(record) is None
    for prefix in ([], ((True, 0),), ((0, 0), True), ((0, 0), 0), ((0, 2),), ((0, 0),)*6):
        with pytest.raises(ValueError):
            subject.Record(prefix)
    for code in (True, -1, 209, 1.0):
        with pytest.raises(ValueError):
            subject.decode_record(code)
    with pytest.raises(ValueError):
        subject.Record.restore({"phase": 1, "prefix": [], "loads": []})
    with pytest.raises(ValueError):
        subject.encode_record(subject.Record(((0, -1),)))
    with pytest.raises(ValueError):
        subject.native_step((2,)+(0,)*11, subject.Record(((0, 0),)), 0)
    with pytest.raises(ValueError):
        subject.Record(next(p for p in subject.record_codes() if len(p) == 5)).append(1)
    halted = subject.Record(next(p for p in subject.record_codes() if len(p) == 5))
    for loads in ((1, "bad")+(0,)*10, (3,)+(0,)*11, (True, 1)+(0,)*10):
        with pytest.raises(ValueError):
            subject.native_step(loads, halted, 0)


@pytest.fixture(scope="module")
def q13_laws():
    return subject.field_laws(ScalarField(13))


def test_actual_conditional_instrument_laws_match_independent_weyl_formula(q13_laws):
    laws, checks = q13_laws
    base = WeylKernel.from_field(ScalarField(13), times=(19/20, 1), current=1, gamma=4, couplings=(1, 1))
    for j in (-1, 0, 1):
        expected = base.for_current(j).two_read_law()
        for (a, b), probability in expected.items():
            assert laws["retained"][j, a, b] == pytest.approx(probability, abs=4e-14)
    assert checks["maximum_conditional_norm_error"] < 2e-13
    assert checks["maximum_checkpoint_continuation_probability_error"] < 2e-13
    assert len(checks["terminal_field_branches"]) == 12
    for item in checks["terminal_field_branches"]:
        state = item["conditional_field_descriptor"]
        assert "kernel" not in state
        assert len(state["half_shifts"]) == len(state["coefficients"]) == 4


def test_execution_does_not_use_the_closed_law_as_a_probability_oracle(monkeypatch):
    def forbidden(self):
        raise AssertionError("closed law is a test oracle only")
    monkeypatch.setattr(WeylKernel, "two_read_law", forbidden)
    laws, _ = subject.field_laws(ScalarField(3))
    assert set(laws) == {"retained", "conditional_field_donor", "field_reset",
                         "independent_second_bit", "commuting_gaussian"}


def test_field_donor_keeps_marginals_and_disturbance_but_deletes_retained_correlation(q13_laws):
    laws, _ = q13_laws
    retained, donor = laws["retained"], laws["conditional_field_donor"]
    for j in (-1, 0, 1):
        for a in (-1, 1):
            assert sum(retained[j, a, b] for b in (-1, 1)) == pytest.approx(sum(donor[j, a, b] for b in (-1, 1)), abs=1e-14)
        for b in (-1, 1):
            assert sum(retained[j, a, b] for a in (-1, 1)) == pytest.approx(sum(donor[j, a, b] for a in (-1, 1)), abs=1e-14)
    k0 = sum(a*b*retained[0, a, b] for a, b in product((-1, 1), repeat=2))
    assert k0 > .01
    assert abs(sum(a*b*donor[0, a, b] for a, b in product((-1, 1), repeat=2))) < 1e-14
    actual_prediction = subject.evaluate(retained)[0]
    donor_prediction = subject.evaluate(donor)[0]
    mu1 = sum(a*retained[1, a, b] for a, b in product((-1, 1), repeat=2))
    mu2 = sum(b*retained[1, a, b] for a, b in product((-1, 1), repeat=2))
    k1 = sum(a*b*retained[1, a, b] for a, b in product((-1, 1), repeat=2))
    assert actual_prediction["third_current"]-donor_prediction["third_current"] == pytest.approx(
        (k0-k1+mu1*mu2)/264, abs=2e-14)


def test_correlated_commuting_control_separates_memory_from_commutator_backaction(q13_laws):
    laws, checks = q13_laws
    actual, classical, reset = laws["retained"], laws["commuting_gaussian"], laws["field_reset"]
    cross_imaginary = checks["kernel_statistics"]["covariance_imag"][0][1]
    for j in (-1, 0, 1):
        a_mean = sum(b*actual[j, a, b] for a, b in product((-1, 1), repeat=2))
        c_mean = sum(b*classical[j, a, b] for a, b in product((-1, 1), repeat=2))
        r_mean = sum(b*reset[j, a, b] for a, b in product((-1, 1), repeat=2))
        assert a_mean == pytest.approx(math.cos(cross_imaginary)*c_mean, abs=4e-14)
        assert c_mean == pytest.approx(r_mean, abs=4e-14)
        # The commuting Gaussian alternative retains the same temporal parity
        # correlation. This is why a positive y1*y2 alone is not a quantum test.
        assert sum(a*b*actual[j, a, b] for a, b in product((-1, 1), repeat=2)) == pytest.approx(
            sum(a*b*classical[j, a, b] for a, b in product((-1, 1), repeat=2)), abs=4e-14)
    quantum_second = sum(b*actual[1, a, b] for a, b in product((-1, 1), repeat=2))
    classical_second = sum(b*classical[1, a, b] for a, b in product((-1, 1), repeat=2))
    assert abs(quantum_second-classical_second) > 1e-6


def test_forecast_deletion_is_matched_to_the_same_native_action_and_field_state(q13_laws):
    laws, _ = q13_laws
    stats, full, coarse = subject.evaluate(laws["retained"])
    assert len(full) == 48
    assert stats["native_record_deletion_brier_gain"] > 0
    assert stats["field_record_deletion_brier_gain"] > 0
    assert stats["native_within_latest_pair_shuffle_brier"]-stats["native_full_history_brier"] == pytest.approx(
        2*stats["native_record_deletion_brier_gain"], abs=2e-14)
    for prefix in full:
        assert (prefix[2], prefix[3]) in coarse
        assert subject.selected_seam(subject.Record(prefix)) == subject.decoder_seam(prefix[3])


def test_small_execution_replays_and_accounts_for_every_retained_history():
    result = subject.run(trials=2000, seed=731, grids=(13,))
    assert result == subject.run(trials=2000, seed=731, grids=(13,))
    assert result["all_pass"]
    resources = result["resources_per_history"]
    assert resources["field_preparations"] == resources["source_preparations"] == 1
    assert resources["field_resets"] == resources["source_resets"] == 0
    assert resources["native_events"] == 3 and resources["field_observations"] == 2
    sampled = result["executions"][0]["sample"]
    assert sum(row[1] for row in sampled["record_counts"]) == 2000
    assert all(subject.decode_record(code).halted and count > 0 for code, count in sampled["record_counts"])
    assert len({code for code, _ in sampled["record_counts"]}) == len(sampled["record_counts"])
    assert len(sampled["first_eight_histories"]) == 8
    assert all(trace["halted"] for trace in sampled["first_eight_histories"])
    for trial_count, seed in ((True, 1), (1, 1), (2, True)):
        with pytest.raises(ValueError):
            subject.sample({}, {}, trial_count, seed)
    for grids in ((), (3,), (True,), (13, 13), [233, 3], "13"):
        with pytest.raises(ValueError):
            subject.run(trials=2, grids=grids)


def test_finest_executed_witnesses_reject_false_zero_memory_and_wrong_commutators():
    from . import continuation_bounds
    certificate = continuation_bounds.run()
    laws, checked = subject.field_laws(ScalarField(233))
    actual, commuting = laws["retained"], laws["commuting_gaussian"]
    moments = {"kappa0": sum(a*b*actual[0, a, b] for a, b in product((-1, 1), repeat=2)),
               "commutator_backaction_signed_mean_difference": sum(
                   b*(actual[1, a, b]-commuting[1, a, b]) for a, b in product((-1, 1), repeat=2))}
    stats = checked["kernel_statistics"]
    assert all(subject.certificate_checks(233, stats, moments, certificate).values())
    for wrong_cross in (0., -stats["covariance_imag"][0][1], 2.):
        wrong = deepcopy(stats)
        wrong["covariance_imag"][0][1] = wrong_cross
        wrong["covariance_imag"][1][0] = -wrong_cross
        checks = subject.certificate_checks(233, wrong, moments, certificate)
        assert not checks["commutator_inside_rational_enclosure"]
    for key in ("kappa0", "commutator_backaction_signed_mean_difference"):
        false = dict(moments, **{key: 0.})
        assert not all(subject.certificate_checks(233, stats, false, certificate).values())
    wrong_grid = dict(stats, mode_count=1728)
    assert not subject.certificate_checks(233, wrong_grid, moments, certificate)["kernel_matches_certified_grid_and_times"]
    with pytest.raises(ValueError):
        subject.certificate_checks(3, stats, moments, certificate)
