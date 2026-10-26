"""Independent source, stopping-cost, charge-ambiguity and consensus controls."""
from collections import defaultdict
from fractions import Fraction as F
import math

import pytest

from . import history_entropy as h, common_variation as v
from ..native_modular_reference import local_preparation as local
from ..native_modular_reference import common_caps as common
from ..observer import native


@pytest.mark.parametrize("law", ["C", "E"])
def test_complete_native_history_reversal_current_balance_and_record_decoder(law):
    source = native.states(binary=True)
    paths = h.histories(law)
    totals = [F() for _ in source]
    flag = [F() for _ in source]
    crossings = 0
    feedback_seams = set()
    for path in paths:
        totals[path.initial] += path.chance
        reconstructed = list(path.endpoint)
        for record in reversed(path.tape):
            a, b = record[:2]
            transfer = h.current(record, law)
            reconstructed[a] += transfer
            reconstructed[b] -= transfer
        assert tuple(reconstructed) == source[path.initial]
        # Independent cap current balance on the same trajectory.
        for cap in ((0, 1), (0, 2, 8), (2, 3, 4, 6)):
            delta = sum(path.endpoint[i]-source[path.initial][i] for i in cap)
            flux = sum(h.current(r, law)*((r[1] in cap)-(r[0] in cap)) for r in path.tape)
            assert delta == flux
        crossings += h.current(path.tape[0], law) != 0
        feedback_seams.add(path.tape[1][:2])
        if path.terminal == "success":
            flag[path.initial] += path.chance
            assert h.decode_core(path.tape, 2, law) == source[path.initial][:2]
        else:
            with pytest.raises(ValueError):
                h.decode_core(path.tape, 2, law)
    assert totals == [F(1)]*66
    assert flag == [F(9, 16) if law == "C" else F(1)]*66
    assert crossings > 0 and feedback_seams == {(0, 1), (1, 8)}


def test_actual_record_first_law_exact_and_negative_fisher_curvature():
    sigma = local.reference()
    paths = tuple(p for p in h.histories() if p.terminal == "success")
    reference = h.ledger(paths, sigma)
    assert reference["functional"] == h.source_entropy(sigma)
    # A general source perturbation checks both moment multipliers, rather
    # than merely a tangent engineered to leave the observable unchanged.
    p = tuple(F(i+1, sum(range(1, 67))) for i in range(66))
    moved = h.ledger(paths, p)
    m, c = v.moments(sigma)
    dm, dc = (a-b for a, b in zip(v.moments(p), (m, c), strict=True))
    cert = v.family_certificate(m, c)
    linear = defaultdict(F)
    v._log_add(linear, dm, cert["linear_multiplier_log_argument"])
    v._log_add(linear, dc, cert["pair_multiplier_log_argument"])
    relative = h.relative_entropy(moved["probabilities"], reference["probabilities"])
    assert h._sum_forms(moved["functional"], relative) == h._sum_forms(reference["functional"], linear)
    # A within-q perturbation that is seen by the actual current readout.
    basis = native.states(binary=True)
    seen = basis.index(tuple(int(i in (2, 5)) for i in range(12)))
    unseen = basis.index(tuple(int(i in (2, 3)) for i in range(12)))
    direction = [F()]*66
    direction[seen], direction[unseen] = F(1), F(-1)
    step = F(1, 100000)
    plus = h.ledger(paths, tuple(w+step*x for w, x in zip(sigma, direction, strict=True)))
    minus = h.ledger(paths, tuple(w-step*x for w, x in zip(sigma, direction, strict=True)))
    delta = {key: (plus["probabilities"][key]-minus["probabilities"][key])/(2*step)
             for key in reference["probabilities"]}
    fisher = sum(float(d*d/reference["probabilities"][key]) for key, d in delta.items())
    assert fisher > 0
    numerical = (h.value(plus["functional"])+h.value(minus["functional"])
                 -2*h.value(reference["functional"]))/float(step*step)
    assert abs(numerical+fisher) < 2e-3
    for row in (plus, minus):
        deficit = h.relative_entropy(row["probabilities"], reference["probabilities"])
        assert h._sum_forms(row["functional"], deficit) == reference["functional"]


def test_native_timeouts_are_retained_and_have_nonzero_unresolved_charge_excess():
    weights = local.reference()
    paths = h.histories()
    rows = h.fibres(paths)
    unresolved = [row for row in rows.values() if len({v.charges()[p.initial] for p in row}) > 1]
    assert unresolved and all(row[0].terminal == "timeout" for row in unresolved)
    full = h.ledger(paths, weights)
    excess = h.value(full["functional"])-h.value(h.source_entropy(weights))
    assert 0 < excess <= float(F(7, 16))*math.log(5)
    # No hidden full-state read is allowed to masquerade as a native record.
    with pytest.raises(ValueError, match="charge-resolved"):
        h.conditional_cap(unresolved[0], (0, 1))
    with pytest.raises(ValueError, match="state-independent"):
        h.flag_probability(tuple(p for p in paths if p.initial != 0))


def test_same_history_cap_predictions_agree_on_nested_and_overlapping_regions():
    paths = tuple(p for p in h.histories() if p.terminal == "success")
    caps = ((0, 1, 2, 5), (0, 1, 5, 8))
    overlap = (0, 1, 5)
    nontrivial = 0
    for row in h.fibres(paths).values():
        common = h.conditional_cap(row, overlap)
        nontrivial += len(row) > 1
        for cap in caps:
            bigger = h.conditional_cap(row, cap)
            reduced = defaultdict(F)
            for bits, probability in bigger.items():
                reduced[tuple(bits[cap.index(i)] for i in overlap)] += probability
            assert dict(reduced) == common
    assert nontrivial
    target = h.source_entropy(local.reference())
    for cap in ((), *caps, overlap, tuple(range(12))):
        assert h.ledger(paths, local.reference(), cap)["functional"] == target


def test_stopping_cost_is_retained_and_one_forgotten_record_can_mix_charge():
    paths = tuple(p for p in h.histories() if p.terminal == "success")
    ledger = h.ledger(paths, local.reference())
    # Census stopping lengths genuinely vary. Ignoring the retained coin
    # entropy or the success normalization changes the equality.
    assert {len(p.tape) for p in paths} == {6, 7, 8}
    assert ledger["mean_native_attempts"] == F(20, 3)
    assert h._sum_forms(ledger["bulk"], ledger["edge"]) != h.source_entropy(local.reference())
    forgotten = defaultdict(set)
    for path in paths:
        # Keeping only the fair-coin word preserves no native load response.
        forgotten[tuple(r[2] for r in path.tape)].add(v.charges()[path.initial])
    assert all(charges == {0, 1, 2} for charges in forgotten.values())


def test_record_dependent_stopping_uses_the_perturbed_inputs_own_cost():
    paths = h.histories("E", prefix=3, stop_on_current=True)
    assert h.flag_probability(paths) == 1
    for path in paths:
        count = h.prefix_length(path.tape, "E", 3, True)
        assert h.decode_core(path.tape, count, "E") == native.states(binary=True)[path.initial][:2]
    sigma = local.reference()
    reference = h.ledger(paths, sigma)
    assert reference["functional"] == h.source_entropy(sigma)
    per_input = [sum((p.chance*len(p.tape) for p in paths if p.initial == i), F())
                 for i in range(66)]
    # Keep q and pair moments fixed, but change the stopping distribution.
    q0 = [i for i, q in enumerate(v.charges()) if q == 0]
    low, high = min(q0, key=per_input.__getitem__), max(q0, key=per_input.__getitem__)
    assert per_input[low] < per_input[high]
    moved_weights = list(sigma)
    moved_weights[low] += F(1, 1000)
    moved_weights[high] -= F(1, 1000)
    moved = h.ledger(paths, tuple(moved_weights))
    assert v.moments(tuple(moved_weights)) == v.moments(sigma)
    assert moved["mean_native_attempts"] != reference["mean_native_attempts"]
    divergence = h.relative_entropy(moved["probabilities"], reference["probabilities"])
    assert h._sum_forms(moved["functional"], divergence) == reference["functional"]
    wrong = h._sum_forms(moved["bulk"], moved["edge"], reference["coin_cost"], divergence)
    assert wrong != reference["functional"]


def test_finite_sandwich_bounds_record_tangents_and_unread_directions_stay_flat():
    paths = h.histories("E", prefix=3, stop_on_current=True)
    sigma = local.reference()
    reference = h.ledger(paths, sigma)
    radius = F(1, 10)
    moved = list(sigma)
    basis = native.states(binary=True)
    a = basis.index(tuple(int(i in (2, 5)) for i in range(12)))
    b = basis.index(tuple(int(i in (2, 3)) for i in range(12)))
    moved[a] += radius*sigma[a]
    moved[b] -= radius*sigma[b]
    assert v.moments(tuple(moved)) == v.moments(sigma)
    observed = h.ledger(paths, tuple(moved))
    variation = {key: (observed["probabilities"][key]-reference["probabilities"][key])/radius
                 for key in reference["probabilities"]}
    derivative = -sum(float(delta)*math.log(float(observed["probabilities"][key]
                                                 /reference["probabilities"][key]))
                      for key, delta in variation.items())
    norm = 2*float(sigma[a])
    assert abs(derivative) <= 2*float(radius/(1-radius))*norm
    # Two untouched spectator states are distinguishable by another readout,
    # but are a real null direction of this protocol's Fisher response.
    hidden = list(sigma)
    a = basis.index(tuple(int(i in (9, 10)) for i in range(12)))
    b = basis.index(tuple(int(i in (10, 11)) for i in range(12)))
    hidden[a] += F(1, 1000)
    hidden[b] -= F(1, 1000)
    assert h.ledger(paths, tuple(hidden))["probabilities"] == reference["probabilities"]


@pytest.mark.parametrize("deadline", [2, 3, 4])
def test_deadline_bound_is_a_complete_state_independent_native_instrument(deadline):
    paths = h.histories(prefix=0, deadline=deadline)
    successful = tuple(p for p in paths if p.terminal == "success")
    alpha = 1-F(1, 2**(deadline-1))
    assert h.flag_probability(paths) == 1
    assert h.flag_probability(successful) == alpha**2
    assert max(map(lambda p: len(p.tape), paths)) == 2*deadline
    assert h.ledger(successful, local.reference())["functional"] == h.source_entropy(local.reference())


def test_c_and_e_readout_laws_both_realize_the_same_entropy_result():
    report = h.run()
    for row in report["laws"].values():
        assert row["exact_source_entropy"] == row["exact_success_functional"]
        assert row["successful_completion_dimensions"][-1] > 1
        assert row["same_history_cap_refinements_checked"] == 6
    assert report["bounded_C_readout_service"]["native_attempts_upper"] == 50
    assert report["laws"]["E"]["complete_counting_excess"] == {}
    assert report["unconditional_E_readout_service"]["native_attempts"] == 4


@pytest.mark.parametrize("weights", [common.reference(), v.family(F(3, 5), F(1, 10))])
def test_envelope_uses_the_full_faithful_moment_family(weights):
    paths = h.histories("E", prefix=3, stop_on_current=True)
    assert h.ledger(paths, weights)["functional"] == h.source_entropy(weights)
