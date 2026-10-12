"""Independent endpoint, channel, record and interference checks."""
from collections import defaultdict
from fractions import Fraction as Q
import math
import unittest

from model.predictions import observer_instrument as m


def dense_vector_oracle(candidate, setting):
    """Direct floor/ceiling construction; never calls the native endpoint map."""
    phase = {"phase_0": 1, "phase_pi": -1, "phase_pi_over_2": 1j}[setting]
    vector = {m.X: 1 / math.sqrt(2), m.Z: phase / math.sqrt(2)}
    branches = defaultdict(lambda: [0j] * 78)
    for source, amplitude in vector.items():
        state = m.STATES[source]
        for coin in (0, 1):
            out = list(state)
            total = state[0] + state[1]
            low, high = total // 2, (total + 1) // 2
            out[0], out[1] = (low, high) if coin == 0 else (high, low)
            current = state[0] - out[0]
            label = (coin, current if candidate == "C" else state[0] - state[1])
            branches[label][m.INDEX[tuple(out)]] += amplitude / math.sqrt(2)
    return {label: (abs((vector[m.X] + vector[m.Z]) / math.sqrt(2))**2,
                   sum(abs(amplitude)**2 for amplitude in vector))
            for label, vector in branches.items()}


class ObserverInstrumentTests(unittest.TestCase):
    def test_exact_predictions_and_reversal(self):
        self.assertEqual(m.prediction_table(), {
            "C": {"phase_0": "5/8", "phase_pi": "1/8",
                  "phase_pi_over_2": "3/8", "diagonal_control": "3/8"},
            "E": {setting: "3/8" for setting in m.SETTINGS}})

    def test_kraus_gram_identity_on_full_carrier(self):
        for candidate in m.CANDIDATES:
            gram = [[Q() for _ in m.STATES] for _ in m.STATES]
            for entries in m.kraus_support(candidate).values():
                for a, out_a in entries:
                    for b, out_b in entries:
                        if out_a == out_b:
                            gram[a][b] += Q(1, 2)
            self.assertTrue(all(gram[a][b] == int(a == b)
                                for a in range(78) for b in range(78)))
            receipt = m.cptp_certificate(candidate)
            self.assertEqual(receipt["nonzero_kraus_entries"], 156)
            self.assertEqual(receipt["record_count"], 6 if candidate == "C" else 10)

    def test_complete_events_and_independent_coherent_oracle(self):
        for candidate in m.CANDIDATES:
            for setting in m.SETTINGS[:-1]:
                oracle = dense_vector_oracle(candidate, setting)
                events = m.event_outputs(candidate, setting)
                self.assertEqual(len(events), 12 if candidate == "C" else 20)
                self.assertEqual(sum(Q(event["probability"]) for event in events), 1)
                for event in events:
                    record = event["record"]
                    label = (record["coin"], record[
                        "current" if candidate == "C" else "contrast"])
                    plus, total = oracle.get(label, (0, 0))
                    expected = plus if event["outcome"] == 1 else total - plus
                    self.assertAlmostEqual(float(Q(event["probability"])), expected, places=14)
                    trace = sum(Q(entry["real"]) for entry in event["source_density"]
                                if entry["row"] == entry["column"])
                    self.assertEqual(trace, Q(event["probability"]))

    def test_readout_preserves_every_committed_record(self):
        for candidate in m.CANDIDATES:
            for setting in m.SETTINGS:
                for block in m.instrument_density(candidate, m.prepared_density(setting)).values():
                    total = sum((m.trace(m.readout_density(block, outcome))
                                 for outcome in (-1, 1)), Q())
                    self.assertEqual(total, m.trace(block))
                    for outcome in (-1, 1):
                        after = m.readout_density(block, outcome)
                        self.assertEqual(m.readout_density(after, outcome), after)
                        self.assertEqual(m.readout_density(after, -outcome), {})

    def test_all_classical_basis_inputs_are_observationally_identical(self):
        for source in range(78):
            density = {(source, source): (Q(1), Q())}
            channels = []
            for candidate in m.CANDIDATES:
                reduced = {}
                for record, block in m.instrument_density(candidate, density).items():
                    coin, value = record
                    current = value if candidate == "C" else (
                        value - (2 * coin - 1) * (abs(value) % 2)) // 2
                    for pair, entry in block.items():
                        key = (coin, current, pair)
                        reduced[key] = tuple(a + b for a, b in zip(
                            reduced.get(key, (Q(), Q())), entry))
                channels.append(reduced)
            self.assertEqual(*channels)

    def test_phase_controls_have_identical_diagonal_preparations(self):
        expected = {(m.X, m.X): (Q(1, 2), Q()),
                    (m.Z, m.Z): (Q(1, 2), Q())}
        for setting in m.SETTINGS:
            density = m.prepared_density(setting)
            self.assertEqual({key: value for key, value in density.items()
                              if key[0] == key[1]}, expected)
            for (a, b), (real, imaginary) in density.items():
                self.assertEqual(density[b, a], (real, -imaginary))
            self.assertEqual(m.trace(density), 1)

    def test_contrast_record_coarse_graining_does_not_restore_interference(self):
        def recorded_plus(candidate):
            result = defaultdict(Q)
            for event in m.event_outputs(candidate, "phase_0"):
                if event["outcome"] == 1:
                    record = event["record"]
                    result[record["coin"], record["current"]] += Q(event["probability"])
            return result
        self.assertEqual(recorded_plus("C")[0, 0], Q(1, 2))
        self.assertEqual(recorded_plus("E")[0, 0], Q(1, 4))
        for candidate in m.CANDIDATES:
            distribution = m.public_distribution(candidate, "phase_0")
            self.assertEqual(len(distribution), 12)
            self.assertEqual(sum(distribution.values()), 1)
            self.assertEqual(sum(probability for (_, _, outcome), probability
                                 in distribution.items() if outcome == 1),
                             m.prediction(candidate, "phase_0"))

    def test_invalid_symbols_rejected(self):
        for candidate in ("c", "unknown", None):
            with self.assertRaises(ValueError):
                m.prediction(candidate, "phase_0")
        with self.assertRaises(ValueError):
            m.prediction("C", "phase_fitted")
        with self.assertRaises(ValueError):
            m.readout_density({}, True)
        with self.assertRaises(ValueError):
            m.instrument_density("C", {(True, 0): (Q(1), Q())})


if __name__ == "__main__":
    unittest.main()
