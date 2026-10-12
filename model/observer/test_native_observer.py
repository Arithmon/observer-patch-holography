"""Independent endpoint arithmetic and finite observer controls."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import unittest

from model.observer import native as subject


def direct_event(loads, seam, coin):
    """Independent floor/ceiling oracle; current is reconstructed afterward."""
    a, b = seam
    total = loads[a] + loads[b]
    low, high = total // 2, (total + 1) // 2
    out = list(loads)
    out[a], out[b] = (low, high) if coin == 0 else (high, low)
    return tuple(out), (coin, loads[a] - out[a])


def configurations(binary=False):
    pairs = combinations(range(12), 2) if binary else combinations_with_replacement(range(12), 2)
    for a, b in pairs:
        x = [0] * 12
        x[a] += 1
        x[b] += 1
        yield tuple(x)


def independent_paths():
    for x in configurations(binary=True):
        for c1, c2 in product((0, 1), repeat=2):
            y, r = direct_event(x, (0, 1), c1)
            seam = (0, 5) if r[1] > 0 else (0, 1)
            z, s = direct_event(y, seam, c2)
            yield x, r, seam, s, z


class ObserverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = subject.run()

    def test_generated_state_carrier(self):
        self.assertEqual(subject.states(), tuple(sorted(configurations())))
        self.assertEqual(subject.states(binary=True), tuple(sorted(configurations(True))))
        self.assertEqual(len(subject.states()), 78)
        self.assertEqual(len(subject.states(binary=True)), 66)

    def test_independent_all_native_columns(self):
        for x in configurations():
            for seam in ((0, 1), (0, 5)):
                for coin in (0, 1):
                    self.assertEqual(subject.native(x, seam, coin), direct_event(x, seam, coin))

    def test_all_history_rows_and_preissued_forecasts(self):
        states = subject.states()
        records = {(tuple(states[row["initial_state"]]), tuple(row["first_record"]),
                    tuple(row["second_seam"]), tuple(row["second_record"]),
                    tuple(states[row["final_state"]])): row
                   for row in self.result["main_histories"]}
        oracle = list(independent_paths())
        self.assertEqual(len(records), 264)
        self.assertEqual(set(records), set(oracle))
        counts = defaultdict(int)
        successes = defaultdict(int)
        for x, r, seam, s, z in oracle:
            counts[r] += 1
            successes[r] += s[1] > 0
            row = records[x, r, seam, s, z]
            self.assertEqual(row["probability"], "1/264")
        for row in self.result["forecast"]["table"]:
            r = tuple(row["record"])
            self.assertEqual(Q(row["forecast"]), Q(successes[r], counts[r]))
        for key, row in records.items():
            self.assertEqual(Q(row["issued_before_event_2"]), Q(successes[key[1]], counts[key[1]]))

    def test_forecast_cannot_read_working_state_or_future(self):
        # Deliberately give a duck-typed state no loads or future at all.
        class PrefixOnly:
            phase = 1
            def __init__(self, record):
                self.prefix = (record,)
        for row in self.result["forecast"]["table"]:
            state = PrefixOnly(tuple(row["record"]))
            self.assertEqual(subject.issue_forecast(state), Q(row["forecast"]))
            self.assertEqual(subject.selected_seam(state),
                             (0, 5) if state.prefix[0][1] > 0 else (0, 1))

    def test_independent_predictive_risks_and_action_matched_shuffle(self):
        paths = list(independent_paths())
        nr, ny, ns, nys = (defaultdict(int) for _ in range(4))
        for _, r, seam, s, _ in paths:
            nr[r] += 1
            ny[r] += s[1] > 0
            ns[seam] += 1
            nys[seam] += s[1] > 0
        forecast = {r: Q(ny[r], nr[r]) for r in nr}
        risk = sum((forecast[r] - (s[1] > 0)) ** 2 for _, r, _, s, _ in paths) / 264
        prior = Q(sum(ny.values()), 264)
        marginal_risk = prior * (1 - prior)
        action_prior = sum((Q(nys[seam], ns[seam]) - (s[1] > 0)) ** 2
                           for _, _, seam, s, _ in paths) / 264
        global_shuffle = sum(Q(nr[r], 264 * 264) * (forecast[r] - (s[1] > 0)) ** 2
                             for _, _, _, s, _ in paths for r in nr)
        action_shuffle = sum(Q(nr[r], 264 * ns[seam]) * (forecast[r] - (s[1] > 0)) ** 2
                             for _, _, seam, s, _ in paths for r in nr
                             if ((0, 5) if r[1] > 0 else (0, 1)) == seam)
        expected = {"brier_risk": risk, "prior_only_brier_risk": marginal_risk,
                    "independent_shuffle_brier_risk": global_shuffle,
                    "action_known_prior_brier_risk": action_prior,
                    "within_action_shuffle_brier_risk": action_shuffle}
        for key, value in expected.items():
            self.assertEqual(Q(self.result["forecast"][key]), value)
        self.assertLess(risk, action_prior)
        self.assertLess(action_prior, action_shuffle)
        self.assertEqual(action_shuffle - risk, 2 * (action_prior - risk))

    def test_feedback_is_not_just_changed_seam_frequencies(self):
        paths = list(independent_paths())
        real = Q(sum(z[0] for _, _, _, _, z in paths), 264)
        # Draw a replacement first record independently from an actual donor.
        # This independent 132 by 132 construction preserves its marginal.
        firsts = [direct_event(x, (0, 1), c) for x in configurations(True) for c in (0, 1)]
        shuffled = Q(0)
        disconnected = Q(0)
        for y, _ in firsts:
            for coin in (0, 1):
                disconnected += Q(direct_event(y, (0, 1), coin)[0][0], 264)
            for _, donor_r in firsts:
                seam = (0, 5) if donor_r[1] > 0 else (0, 1)
                for coin in (0, 1):
                    shuffled += Q(direct_event(y, seam, coin)[0][0], 132 * 132 * 2)
        self.assertEqual(real, Q(35, 264))
        self.assertEqual(shuffled, disconnected)
        self.assertEqual(disconnected, Q(1, 6))
        self.assertEqual(real - shuffled, -Q(3, 88))

    def test_checkpoint_roundtrip_all_first_states_and_halt(self):
        for x in configurations():
            for r in subject.LETTERS:
                state = subject.State(x, (r,))
                restored = subject.State.restore(json.loads(json.dumps(state.checkpoint())))
                self.assertEqual(subject.tick(state), subject.tick(restored))
                self.assertEqual(subject.issue_forecast(state), subject.issue_forecast(restored))
                self.assertEqual(state.controller_checkpoint(), {"phase": 1, "prefix": [list(r)]})
                local = subject.State.restore_controller(state.controller_checkpoint(), x)
                self.assertEqual(subject.tick(local), subject.tick(state))
                self.assertEqual(subject.issue_forecast(local), subject.issue_forecast(state))
                for _, last in subject.tick(state):
                    self.assertEqual(subject.tick(last), ((Q(1), last),))
        bad = subject.State(tuple(next(configurations())), ((0, 0),)).checkpoint()
        bad["phase"] = 0
        with self.assertRaises(ValueError):
            subject.State.restore(bad)

    def test_generic_complex_coherence_against_dense_small_stinespring(self):
        # Independent complex-state check, including arbitrary record inputs.
        # The exact producer proves all columns; this catches prefix-trace mistakes.
        inputs = [(r, x) for r in subject.LETTERS for x in configurations()]
        amplitudes = {key: complex((k % 7) - 3, (k % 5) - 2)
                      for k, key in enumerate(inputs)}
        norm = sum(abs(a) ** 2 for a in amplitudes.values())
        out = {}
        for (r, x), amplitude in amplitudes.items():
            seam = (0, 5) if r[1] > 0 else (0, 1)
            for coin in (0, 1):
                y, s = direct_event(x, seam, coin)
                key = (r, s, y)
                self.assertNotIn(key, out)
                out[key] = amplitude / (2 ** 0.5)
        self.assertAlmostEqual(sum(abs(a) ** 2 for a in out.values()), norm, places=9)
        before = defaultdict(float)
        after = defaultdict(float)
        for (r, _), amplitude in amplitudes.items():
            before[r] += abs(amplitude) ** 2
        for (r, _, _), amplitude in out.items():
            after[r] += abs(amplitude) ** 2
        for r in subject.LETTERS:
            self.assertAlmostEqual(before[r], after[r], places=9)
        witness = self.result["quantum_prefix_counterexample"]
        self.assertEqual((witness["before"], witness["after"]), ("1/4", "1/8"))

    def test_canonical_replay(self):
        path = Path(__file__).with_name("results") / "native.json"
        self.assertEqual(path.read_text(), subject.canonical(self.result))


if __name__ == "__main__":
    unittest.main()
