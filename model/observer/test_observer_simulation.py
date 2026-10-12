"""Sequential execution must preserve the complete exact finite experiment."""
import json
from pathlib import Path
import unittest

from model.observer import simulate


class ObserverSimulationTests(unittest.TestCase):
    def test_published_seeded_run_replays(self):
        receipt = Path(__file__).with_name("results") / "run.json"
        expected = json.loads(receipt.read_text())
        result = simulate.run(expected["trials"], expected["seed"])
        self.assertEqual(receipt.read_text(), simulate.canonical(result))
        self.assertEqual(result["execution"]["exact_path_reference_matches"], result["trials"])
        self.assertEqual(result["execution"]["native_events"], 2 * result["trials"])
        self.assertEqual(result["execution"]["explicit_halt_checks"], result["trials"])

    def test_invalid_trial_count_is_rejected(self):
        for value in (0, -1, True, 0.5):
            with self.assertRaises(ValueError):
                simulate.run(value)


if __name__ == "__main__":
    unittest.main()
