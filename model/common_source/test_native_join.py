"""Complete native/quantum history composition, including adversarial controls."""
from copy import deepcopy
import math
import unittest
from unittest.mock import patch

from . import native_join as join


def finite_instrument_fixture():
    """A small exact projective instrument for composition tests, not evidence.

    The j=0 detector is deliberately unfair: blindly reusing an old odd/fair
    response formula would give the wrong joint law. Branch vectors are
    orthogonal projections of a normalized declared vector.
    """
    vectors = {-1: [1., 0., 0., 0.], 0: [3/5, 4/5, 0., 0.], 1: [0., 0., 0., 1.]}
    result = {}
    for current, amplitudes in vectors.items():
        branches = []
        for index, record in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
            state = [[0., 0.] for _ in range(16)]
            state[index] = [amplitudes[index], 0.]
            branches.append({"record": list(record),
                             "probability": amplitudes[index]**2,
                             "unnormalized_terminal_state": state})
        result[current] = {"native_current": current, "dimension": 2, "steps_per_segment": 2,
                           "coupling": .25, "feedback_enabled": True,
                           "branches": branches, "all_pass": True}
    return result


def endpoint(loads, seam, coin):
    """Independent floor/ceiling oracle, never calling native.native."""
    a, b = seam
    total = loads[a]+loads[b]
    output = list(loads)
    low, high = total//2, (total+1)//2
    output[a], output[b] = (low, high) if coin == 0 else (high, low)
    return tuple(output), (coin, loads[a]-output[a])


class NativeJoinTests(unittest.TestCase):
    def test_complete_paths_match_independent_endpoint_and_quantum_vectors(self):
        histories = finite_instrument_fixture()
        result = join.compose_histories(histories)
        self.assertEqual(len(result["paths"]), 1056)
        self.assertEqual(len({(row["initial_source"], row["first_native_record"][0],
                               tuple(row["quantum_record"]), row["second_native_record"][0])
                              for row in result["paths"]}), 1056)
        zeros = 0
        for row in result["paths"]:
            initial = join.ALL_STATES[row["initial_source"]]
            after_first, first = endpoint(initial, (0, 1), row["first_native_record"][0])
            self.assertEqual(list(first), row["first_native_record"])
            self.assertEqual(after_first, join.ALL_STATES[row["after_first_source"]])
            a, b = row["quantum_record"]
            expected_seam = (0, 5) if b else (0, 1)
            final, second = endpoint(after_first, expected_seam, row["second_native_record"][0])
            self.assertEqual(row["second_seam"], list(expected_seam))
            self.assertEqual(row["second_native_record"], list(second))
            self.assertEqual(row["final_native_loads"], list(final))
            self.assertEqual(join.ALL_STATES[row["final_source"]], final)
            reference = row["quantum_branch_reference"]
            self.assertEqual(reference["current"], first[1])
            branch = histories[first[1]]["branches"][reference["branch_index"]]
            self.assertEqual(branch["record"], [a, b])
            probability = sum(real**2+imaginary**2
                              for real, imaginary in branch["unnormalized_terminal_state"])
            self.assertEqual(row["terminal_density_prefactor"], "1/264")
            self.assertAlmostEqual(row["history_probability"], probability/264, places=16)
            zeros += row["history_probability"] == 0
        self.assertGreater(zeros, 0)
        self.assertAlmostEqual(result["summary"]["history_probability_sum"], 1, places=14)

    def test_mean_and_matched_independent_control_use_actual_unfair_detector(self):
        result = join.compose_histories(finite_instrument_fixture())
        summary = result["summary"]
        control = result["controls"]["independent_bit_same_marginal"]
        self.assertAlmostEqual(summary["final_scalar_bit_one_probability"],
                               (28/33)*(16/25)+5/66, places=14)
        self.assertNotAlmostEqual(summary["final_scalar_bit_one_probability"], .5)
        self.assertEqual(control["actuator_bit_one_probability"],
                         summary["final_scalar_bit_one_probability"])
        self.assertAlmostEqual(control["mean_second_native_current"], 0, places=14)
        self.assertAlmostEqual(control["mean_final_native_port0_load"], 1/6, places=14)
        self.assertNotAlmostEqual(summary["feedback_current_contrast_against_independent_same_marginal"], 0)
        self.assertAlmostEqual(summary["feedback_current_contrast_against_independent_same_marginal"],
                               -summary["feedback_port0_contrast_against_independent_same_marginal"], places=14)

    def test_disconnected_current_has_zero_feedback_effect_despite_biased_detector(self):
        histories = finite_instrument_fixture()
        for current in (-1, 1):
            histories[current]["branches"] = deepcopy(histories[0]["branches"])
        result = join.compose_histories(histories)
        self.assertAlmostEqual(result["summary"]["final_scalar_bit_one_probability"], 16/25, places=14)
        self.assertAlmostEqual(result["summary"]["feedback_current_contrast_against_independent_same_marginal"],
                               0, places=14)

    def test_first_native_record_survives_and_terminal_states_are_retained_once(self):
        histories = finite_instrument_fixture()
        result = join.compose_histories(histories)
        for history in result["quantum_histories"]:
            self.assertIs(history, histories[history["native_current"]])
        expected = {(0, 0): 14/33, (1, 0): 14/33, (0, 1): 5/66, (1, -1): 5/66}
        for row in result["first_native_record_marginal"]:
            self.assertAlmostEqual(row["probability"], expected[tuple(row["first_native_record"])], places=14)
        self.assertTrue(result["checks"]["controller_uses_only_retained_final_bit"])

    def test_controller_input_is_strictly_the_committed_bit(self):
        self.assertEqual(join.actuator_seam(0), (0, 1))
        self.assertEqual(join.actuator_seam(1), (0, 5))
        for invalid in (True, False, -1, 2, 0., None, (0, 1)):
            with self.assertRaises(ValueError):
                join.actuator_seam(invalid)
        for state in join.ALL_STATES:
            for bit in (0, 1):
                for coin in (0, 1):
                    self.assertEqual(join.native_continuation(state, bit, coin),
                                     endpoint(state, (0, 5) if bit else (0, 1), coin))

    def test_terminal_state_or_complete_record_mutations_are_rejected(self):
        for mutation in ("probability", "amplitude", "duplicate", "missing", "current", "dimension", "failed"):
            histories = finite_instrument_fixture()
            if mutation == "probability":
                histories[0]["branches"][0]["probability"] += .1
            elif mutation == "amplitude":
                histories[0]["branches"][0]["unnormalized_terminal_state"][0][0] = math.nan
            elif mutation == "duplicate":
                histories[0]["branches"][0]["record"] = [0, 1]
            elif mutation == "missing":
                histories[0]["branches"].pop()
            elif mutation == "current":
                histories[0]["native_current"] = 1
            elif mutation == "dimension":
                histories[0]["branches"][0]["unnormalized_terminal_state"].pop()
            elif mutation == "failed":
                histories[0]["all_pass"] = False
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                join.compose_histories(histories)

    def test_run_executes_and_caches_each_current_history(self):
        fixtures = finite_instrument_fixture()
        join._quantum_history.cache_clear()
        with patch.object(join.experiment, "history", side_effect=lambda *args, **kwargs:
                          fixtures[kwargs["current"]]) as call:
            first = join.run(d=2, steps=2)
            second = join.run(d=2, steps=2)
        join._quantum_history.cache_clear()
        self.assertEqual(call.call_count, 3)
        self.assertEqual([invocation.kwargs["current"] for invocation in call.call_args_list], [-1, 0, 1])
        self.assertEqual(first, second)
        for params in ({"d": True}, {"steps": 0}, {"steps": True},
                       {"coupling": math.inf}, {"feedback_enabled": 1}):
            with self.assertRaises(ValueError):
                join.run(**params)

    def test_ideal_temporal_bound_contracts_through_complete_native_history(self):
        histories = finite_instrument_fixture()
        for current in (-1, 0, 1):
            histories[current]["ideal_finite_continuous_comparison"] = {
                "complete_instrument_half_diamond_error_upper": "1/10" if current else "1/20"}
        certificate = join.compose_histories(histories)["finite_temporal_transfer"]
        self.assertTrue(certificate["joint_error_bound_claimed"])
        self.assertEqual(certificate["complete_classical_history_tv_error_upper"], "1/10")
        self.assertEqual(certificate["mean_second_native_current_error_upper"], "1/5")
        histories[0]["ideal_finite_continuous_comparison"][
            "complete_instrument_half_diamond_error_upper"] = "-1/10"
        with self.assertRaises(ValueError):
            join.compose_histories(histories)

    def test_real_common_source_words_drive_native_continuation(self):
        # Small finite action: integration coverage, not a continuum assertion.
        result = join.run(d=2, steps=2)
        self.assertTrue(result["all_pass"])
        self.assertEqual(len(result["paths"]), 1056)
        for history in result["quantum_histories"]:
            self.assertEqual(history["native_current"], history["execution"]["native_current_preparation"]["native_current"])
            self.assertTrue(history["all_pass"])
        self.assertLess(result["maximum_probability_state_norm_discrepancy"], 2e-10)
        self.assertTrue(result["finite_temporal_transfer"]["joint_error_bound_claimed"])


if __name__ == "__main__":
    unittest.main()
