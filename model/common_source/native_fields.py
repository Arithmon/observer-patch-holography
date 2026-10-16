"""Exact continuity and closed mean transport of the actual native repair.

This iterates the existing all-seam calibration from slow_geometry.py: its
fixed twelve-port carrier, total-two integer loads, independent uniform seam
choices and fair placement coins. It does not replace the existing E/F
record-controlled production protocol. With B_e=e_a-e_b, the native record j
obeys x'=x-B_e*j pathwise and E[j|x,e]=B_e.T*x/2. Thus the event-count mean
obeys m[n+1]=(I-L/60)m[n], L=sum_e B_e B_e.T. This reuses the existing carrier
mean identity; the full rational Markov iteration checks its attachment to
actual recorded currents, including continuing fluctuations.

The independently specified agreement subspace is constant vertex readings,
ker(B.T); on this connected carrier it is ker(L). This mean law dissipates
nonconstant modes. For its slow scalar mode a[n]=q**n*a[0], 0<q<1,
a[n+1]-2*a[n]+a[n-1]=(q+1/q-2)*a[n] has positive coefficient. An undamped
sampled oscillator has coefficient -4*sin(omega*delta/2)**2 <= 0. No positive
scalar conversion from event count to time equates these nonzero readouts.
This is not an obstruction for other native histories, quantum sectors,
state-dependent interpretations, or OPH as a whole. Loads/means are external
working-state audits; only the current record is already a native readout.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from ..observer import native
from . import slow_geometry as geometry


def laplacian():
    matrix = [[0] * 12 for _ in range(12)]
    for a, b in geometry.source().seams():
        matrix[a][a] += 1
        matrix[b][b] += 1
        matrix[a][b] -= 1
        matrix[b][a] -= 1
    return tuple(map(tuple, matrix))


def mean_operator():
    return tuple(tuple(F(i == j) - F(value, 60) for j, value in enumerate(row))
                 for i, row in enumerate(laplacian()))


def _matvec(matrix, vector):
    return tuple(sum((a * b for a, b in zip(row, vector, strict=True)), F(0)) for row in matrix)


@lru_cache(maxsize=1)
def transition():
    """Sparse row-stochastic 78-state kernel, enumerated from native outputs."""
    states = native.states()
    indices = {state: i for i, state in enumerate(states)}
    rows = []
    for state in states:
        counts = Counter(indices[native.native(state, seam, coin)[0]]
                         for seam in geometry.source().seams() for coin in (0, 1))
        rows.append(tuple((i, F(count, 60)) for i, count in sorted(counts.items())))
    return tuple(rows)


def _distribution(values):
    values = tuple(values)
    if (len(values) != 78 or any(type(value) not in (int, F) or value < 0 for value in values)
            or sum(values) != 1):
        raise ValueError("an exact normalized nonnegative distribution on all 78 states is required")
    return tuple(map(F, values))


def advance(distribution):
    distribution = _distribution(distribution)
    result = [F(0)] * 78
    for mass, row in zip(distribution, transition(), strict=True):
        for index, probability in row:
            result[index] += mass * probability
    return tuple(result)


def mean(distribution):
    distribution = _distribution(distribution)
    return tuple(sum((p * state[i] for p, state in zip(distribution, native.states(), strict=True)), F(0))
                 for i in range(12))


def conditional_path_surprisal(path):
    """Exact state-path probability and its additive, dimensionless surprisal.

    The initial state is conditioned on, so no start-prior term is included.
    A caller who supplies a prior may add -log(pi(x_0)) separately.
    """
    path = tuple(tuple(state) for state in path)
    states = native.states()
    indices = {state: i for i, state in enumerate(states)}
    if len(path) < 2 or any(state not in indices for state in path):
        raise ValueError("a path of at least two valid total-two load states is required")
    terms, probability = [], F(1)
    for source, target in zip(path, path[1:], strict=False):
        row = dict(transition()[indices[source]])
        edge_probability = row.get(indices[target], F(0))
        if edge_probability <= 0:
            raise ValueError("the path contains an impossible native transition")
        multiplicity = edge_probability * 60
        if multiplicity.denominator != 1:
            raise AssertionError("transition probability is not an integer primitive-branch count over 60")
        ratio = F(1, 1) / edge_probability
        potential_drop = F(sum(value*value for value in source)
                           - sum(value*value for value in target), 2)
        terms.append({"conditional_probability": str(edge_probability),
                      "primitive_branch_multiplicity": multiplicity.numerator,
                      "surprisal_nats": f"log({ratio.numerator}/{ratio.denominator})"
                          if ratio.denominator != 1 else f"log({ratio.numerator})",
                      "quadratic_load_potential_drop": str(potential_drop)})
        probability *= edge_probability
    inverse = F(1, 1) / probability
    return {"states": [list(state) for state in path],
            "conditional_path_probability": str(probability),
            "conditional_path_surprisal_nats":
                f"log({inverse.numerator}/{inverse.denominator})"
                if inverse.denominator != 1 else f"log({inverse.numerator})",
            "initial_state_prior_included": False,
            "local_terms": terms}


def stochastic_action_audit():
    """Audit Markov path surprisal against the native quadratic repair cost."""
    multiplicities = Counter()
    for row in transition():
        for _, probability in row:
            count = probability * 60
            if count.denominator != 1:
                raise AssertionError("kernel entry is not a primitive branch multiplicity over 60")
            multiplicities[count.numerator] += 1
    expected = Counter({1: 600, 2: 60, 50: 48, 52: 30})
    if multiplicities != expected:
        raise AssertionError("exact state-kernel multiplicity classes changed")
    concentrated = (2,) + (0,) * 11
    concentrated_index = native.states().index(concentrated)
    concentrated_target = next(native.states()[j] for j, p in transition()[concentrated_index]
                               if F(sum(a*a for a in concentrated)-sum(a*a for a in native.states()[j]), 2) == 1)
    split = (1, 1) + (0,) * 10
    split_index = native.states().index(split)
    split_zero_drop_target = next(native.states()[j] for j, p in transition()[split_index]
                                  if native.states()[j] != split
                                  and F(sum(a*a for a in split)-sum(a*a for a in native.states()[j]), 2) == 0)
    examples = {
        "concentration_repair": conditional_path_surprisal((concentrated, concentrated_target)),
        "split_state_self_transition": conditional_path_surprisal((split, split)),
        "split_state_zero_drop_move": conditional_path_surprisal((split, split_zero_drop_target)),
    }
    split_row = dict(transition()[split_index])
    most_likely_index = max(split_row, key=split_row.get)
    most_likely_state = native.states()[most_likely_index]
    split_mean_after_one = _matvec(mean_operator(), split)
    if most_likely_state != split or split_mean_after_one == split:
        raise AssertionError("state surprisal minimization unexpectedly equals the mean-field update")
    terms = [row["local_terms"][0] for row in examples.values()]
    if (terms[1]["quadratic_load_potential_drop"] != terms[2]["quadratic_load_potential_drop"]
            or terms[1]["surprisal_nats"] == terms[2]["surprisal_nats"]):
        raise AssertionError("the zero-drop transition controls no longer separate surprisal from repair cost")
    return {"conditional_path_probability": "product_k K[x_k,x_(k+1)]",
            "local_surprisal": "-log K[x,y]=log(60/m[x,y]) nats",
            "start_prior": "excluded; a supplied pi(x0) adds the separate boundary term -log pi(x0)",
            "exact_transition_multiplicity_histogram": {str(k): v for k, v in sorted(multiplicities.items())},
            "local_surprisal_classes_nats": {str(k): f"log({F(60, k).numerator}/{F(60, k).denominator})"
                                             if F(60, k).denominator != 1 else f"log({F(60, k).numerator})"
                                             for k in sorted(multiplicities)},
            "fully_labeled_seam_coin_branch_probability": "1/60",
            "fully_labeled_branch_surprisal_nats": "log(60), constant per event",
            "state_only_path_examples": examples,
            "one_step_most_likely_vs_mean": {
                "initial_state": list(split),
                "most_likely_next_state": list(most_likely_state),
                "most_likely_transition_probability": str(split_row[most_likely_index]),
                "conditional_mean_after_one_event": list(map(str, split_mean_after_one)),
                "greedy_surprisal_path_is_not_the_conditional_mean_update": True},
            "not_proportional_to_quadratic_load_drop": True,
            "not_the_logarithm_of_the_kernel_as_a_linear_operator": True,
            "physical_action_units_derived": False,
            "all_pass": True}


def primitive_audit():
    states, seams, operator = native.states(), geometry.source().seams(), mean_operator()
    branches = 0
    for state, row in zip(states, transition(), strict=True):
        if sum(p for _, p in row) != 1 or any(p <= 0 for _, p in row):
            raise AssertionError("native transition is not stochastic")
        for a, b in seams:
            currents = []
            for coin in (0, 1):
                final, record = native.native(state, (a, b), coin)
                current = record[1]
                expected = tuple(state[i] - (int(i == a) - int(i == b)) * current for i in range(12))
                delta = state[a] - state[b]
                norm_drop = F(delta * delta - (abs(delta) % 2), 2)
                if (final != expected or sum(final) != 2 or min(final) < 0
                        or sum(x*x for x in state) - sum(x*x for x in final) != norm_drop or norm_drop < 0):
                    raise AssertionError("native continuity, mass or rounding dissipation failed")
                currents.append(current)
                branches += 1
            if F(sum(currents), 2) != F(state[a] - state[b], 2):
                raise AssertionError("fair native current is not the half difference")
            difference = state[a] - state[b]
            if F(sum(j*j for j in currents), 2) != F(difference*difference + abs(difference) % 2, 4):
                raise AssertionError("the exact native current noise moment changed")
        actual = tuple(sum(p * states[j][i] for j, p in row) for i in range(12))
        if actual != _matvec(operator, state):
            raise AssertionError("full native transition mean differs from I-L/60")
    if (any(sum(row) != 1 for row in operator)
            or any(sum(operator[i][j] for i in range(12)) != 1 for j in range(12))):
        raise AssertionError("the mean operator does not preserve constant readings and total mass")
    return {"primitive_branches_checked": branches, "conditional_current_means_checked": len(states) * len(seams),
            "markov_rows_checked": len(states), "nonzero_markov_entries": sum(map(len, transition())),
            "pathwise_continuity": "x' = x - (e_a-e_b)*j",
            "conditional_current": "E[j | x,(a,b)] = (x_a-x_b)/2",
            "conditional_current_variance": "Var[j | x,(a,b)] = ((x_a-x_b) mod 2)/4",
            "pathwise_squared_load_norm_drop": "((x_a-x_b)^2 - ((x_a-x_b) mod 2))/2",
            "all_pass": True}


def spectral_audit():
    """Verify the entire declared spectrum against incidence, exactly in Q(sqrt5)."""
    c, lap = geometry.source(), laplacian()
    exact_lap = [[c.q5(value) for value in row] for row in lap]
    total = [[c.q5(0) for _ in range(12)] for _ in range(12)]
    rows = []
    for index, (value, multiplicity) in enumerate(c.LAPLACIAN_BANDS):
        projector = c.band_projector_exact(index)
        if (c.q5_matmul(projector, projector) != list(map(list, projector))
                or c.q5_matmul(exact_lap, projector) != [[c.q5_mul(value, x) for x in row] for row in projector]):
            raise AssertionError("declared spectral band fails its exact projector equations")
        trace = c.q5(0)
        for i in range(12):
            trace = c.q5_add(trace, projector[i][i])
            for j in range(12):
                total[i][j] = c.q5_add(total[i][j], projector[i][j])
        q = c.q5_sub(c.q5(1), c.q5_mul(c.q5(F(1, 60)), value))
        if trace != c.q5(multiplicity) or c.q5_sign(q) <= 0 or (index > 0 and c.q5_sign(c.q5_sub(c.q5(1), q)) <= 0):
            raise AssertionError("mean bands are not positive nonconstant contractions")
        rows.append({"laplacian_eigenvalue_Qsqrt5": c.q5_str(value), "multiplicity": multiplicity,
                     "mean_multiplier_Qsqrt5": c.q5_str(q), "mean_multiplier_numeric": c.q5_float(q)})
    if total != [[c.q5(int(i == j)) for j in range(12)] for i in range(12)]:
        raise AssertionError("exact spectral projectors do not resolve the identity")
    slow = c.q5(F(55, 60), F(1, 60))
    coefficient = c.q5_sub(c.q5_add(slow, c.q5_inv(slow)), c.q5(2))
    if c.q5_sign(coefficient) <= 0:
        raise AssertionError("slow mean is not strictly dissipative")
    return {"bands": rows, "slow_multiplier_Qsqrt5": c.q5_str(slow),
            "slow_multiplier_numeric": c.q5_float(slow),
            "slow_symmetric_second_difference_coefficient_Qsqrt5": c.q5_str(coefficient),
            "slow_symmetric_second_difference_coefficient_numeric": c.q5_float(coefficient),
            "undamped_sampled_oscillator_coefficient": "-4*sin(omega*delta/2)^2 <= 0",
            "squared_mean_norm_drop": "m^T*(L/30-L^2/3600)*m",
            "centered_mean_potential_change": "V(Tm)-V(m)=-z^T*L*z/60+z^T*L^2*z/7200; z=m-1/6, V=||z||^2/2",
            "constant_agreement_subspace_dimension": 1, "all_pass": True}


def chain_audit(steps=16):
    if type(steps) is not int or steps < 1:
        raise ValueError("a positive integer event count is required")
    states, operator = native.states(), mean_operator()
    concentrated = (2,) + (0,) * 11
    initial = tuple(F(state == concentrated) for state in states)
    uniform = (F(1, 78),) * 78
    binary = tuple(F(1, 66) if max(state) == 1 else F(0) for state in states)
    if advance(binary) != binary or advance(uniform) == uniform:
        raise AssertionError("the expected finite stationary/nonstationary controls changed")
    result = []
    for name, distribution in (("two_loads_at_port_zero", initial), ("uniform_all_78_states", uniform)):
        expected, previous_norm = mean(distribution), None
        rows = []
        for event in range(steps + 1):
            actual = mean(distribution)
            norm = sum(x*x for x in actual)
            if actual != expected or sum(actual) != 2 or (previous_norm is not None and norm > previous_norm):
                raise AssertionError("independently iterated Markov law violates mean closure/dissipation")
            if event in (0, 1, 2, 4, 8, steps):
                micro_norm = sum(p * sum(x*x for x in state) for p, state in zip(distribution, states, strict=True))
                rows.append({"native_event_count": event, "mean": list(map(str, actual)),
                             "squared_mean_norm": str(norm), "expected_squared_load_norm": str(micro_norm),
                             "centered_mean_norm_squared": str(norm-F(1, 3))})
            previous_norm = norm
            if event < steps:
                distribution, expected = advance(distribution), _matvec(operator, expected)
        result.append({"initial_distribution": name, "mean_equalities_checked": steps + 1,
                       "samples": rows, "all_pass": True})
    return {"chains": result, "uniform_78_distribution_is_stationary": False,
            "uniform_78_mean_is_preserved": True, "uniform_66_binary_distribution_is_stationary": True,
            "stationary_binary_expected_squared_load_norm": "2",
            "stationary_binary_squared_mean_norm": "1/3",
            "stationary_binary_fluctuation_norm_squared": "5/3", "all_pass": True}


def controls():
    state = (1, 0, 1) + (0,) * 9
    outcomes = [native.native(state, native.E, coin) for coin in (0, 1)]
    biased = F(1, 4) * outcomes[0][1][1] + F(3, 4) * outcomes[1][1][1]
    half_difference = F(state[0]-state[1], 2)
    fixed_mean = tuple(F(outcomes[0][0][i]+outcomes[1][0][i], 2) for i in range(12))
    uniform_mean = _matvec(mean_operator(), state)
    if biased == half_difference or fixed_mean == uniform_mean:
        raise AssertionError("schedule/coin premise negative control failed")
    paths = [(p*q, first, last) for p, first in native.tick(native.State(state))
             for q, last in native.tick(first)]
    policy_mean = tuple(sum(weight * last.loads[i] for weight, _, last in paths) for i in range(12))
    calibration_mean = _matvec(mean_operator(), uniform_mean)
    if policy_mean == calibration_mean or sum(weight for weight, _, _ in paths) != 1:
        raise AssertionError("the actual E/F feedback was incorrectly identified with uniform calibration")
    return {"initial_state": list(state), "fixed_seam": list(native.E), "coin_one_probability": "3/4",
            "biased_conditional_current": str(biased), "fair_conditional_current": str(half_difference),
            "fixed_seam_fair_mean": list(map(str, fixed_mean)), "uniform_seams_fair_mean": list(map(str, uniform_mean)),
            "actual_E_F_feedback_paths": [{"probability": str(weight), "first_record": list(first.prefix[0]),
                                           "second_selected_seam": list(native.selected_seam(first)),
                                           "second_record": list(last.prefix[1]), "final_loads": list(last.loads)}
                                          for weight, first, last in paths],
            "actual_E_F_feedback_two_event_mean": list(map(str, policy_mean)),
            "uniform_calibration_two_event_mean": list(map(str, calibration_mean)),
            "all_pass": True}


def run():
    kernel = [[[i, str(p)] for i, p in row] for row in transition()]
    return {"schema": "oph.common_source.native_current_mean_field.v1",
            "source": {"native_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
                       "transition_exact_sha256": hashlib.sha256(json.dumps(kernel, separators=(",", ":")).encode()).hexdigest()},
            "primitive": primitive_audit(), "spectrum": spectral_audit(), "iteration": chain_audit(),
            "premise_controls": controls(), "stochastic_action": stochastic_action_audit(), "all_pass": True,
            "assumptions": {"graph": "existing fixed 12-port 30-seam carrier",
                            "protocol": "existing slow_geometry all-seam calibration, iterated",
                            "load_sector": "nonnegative integer total two",
                            "schedule": "uniform independent seam choice at each event, independent of loads",
                            "coins": "fair placement coins, independent of loads and seam choices",
                            "parameter": "number of native seam events; no observer-clock identification"},
            "scope": {"mean_law": "m[n+1]=(I-L/60)m[n]", "field_values_are_working_state_audits": True,
                      "native_current_records_used": True, "mean_convergence_implies_static_paths": False,
                      "uniform_calibration_is_actual_E_F_feedback_protocol": False,
                      "undamped_wave_excluded_for_this_nonconstant_mean_readout": True,
                      "other_native_history_readouts_excluded": False,
                      "physical_clock_or_spatial_continuum_derived": False,
                      "canonical_uniform_path_law_derived_from_OPH_axioms": False}}


def main():
    from ..observer.scalar_response.experiment import canonical, compare_receipt
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--output", "--write", dest="output", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding="utf-8")), result)
    print(json.dumps({"all_pass": result["all_pass"], "primitive": result["primitive"],
                      "slow_multiplier": result["spectrum"]["slow_multiplier_numeric"]}, indent=2))


if __name__ == "__main__":
    main()
