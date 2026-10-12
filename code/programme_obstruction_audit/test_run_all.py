"""Portable replay must leave its input package and frozen evidence untouched."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent


def snapshot(root):
    return {
        str(path.relative_to(root)): (
            hashlib.sha256(path.read_bytes()).hexdigest(),
            path.stat().st_mtime_ns,
        )
        for path in root.rglob('*') if path.is_file()
    }


class PackagePreservationTest(unittest.TestCase):
    def test_two_fresh_replays_preserve_package(self):
        with tempfile.TemporaryDirectory(prefix='obstruction-replay-test-') as tmp:
            temporary = Path(tmp)
            package = temporary / 'package'
            shutil.copytree(ROOT, package, ignore=shutil.ignore_patterns(
                '__pycache__', 'run_receipts'))
            before = snapshot(package)
            for run in range(2):
                output = temporary / f'replay-{run}'
                destination_args = ['--output-dir', str(output)] if run == 0 else []
                result = subprocess.run(
                    [sys.executable, '-B', '-W', 'error',
                     str(package / 'run_all.py'), *destination_args],
                    cwd=temporary, capture_output=True, text=True, timeout=900,
                    env={**os.environ, 'PYTHONWARNINGS': 'error', 'TMPDIR': str(temporary)},
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                verdict = json.loads(result.stdout)
                self.assertTrue(verdict['all_seven_controls_passed'])
                summary_path = Path(verdict['full_summary'])
                self.assertFalse(package == summary_path or package in summary_path.parents)
                summary = json.loads(summary_path.read_text())
                self.assertEqual(len(summary['runs']), 7)
                self.assertTrue(all(row['exit_code'] == 0 for row in summary['runs']))
                self.assertEqual(snapshot(package), before,
                                 'Replay changed an input or frozen evidence file')

    def test_unsafe_destinations_preserve_existing_files(self):
        with tempfile.TemporaryDirectory(prefix='obstruction-output-test-') as tmp:
            temporary = Path(tmp)
            package = temporary / 'package'
            shutil.copytree(ROOT, package, ignore=shutil.ignore_patterns(
                '__pycache__', 'run_receipts'))
            retained = temporary / 'retained'
            retained.mkdir()
            (retained / 'summary.json').write_text('{"previous_success": true}\n')
            before = snapshot(temporary)
            for destination, message in [
                (package, 'outside the input package'),
                (package / 'receipts', 'outside the input package'),
                (package / 'new-results', 'outside the input package'),
                (retained, 'must be empty'),
            ]:
                with self.subTest(destination=destination):
                    result = subprocess.run(
                        [sys.executable, '-B', str(package / 'run_all.py'),
                         '--output-dir', str(destination)],
                        cwd=temporary, capture_output=True, text=True, timeout=30,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(message, result.stderr)
                    self.assertEqual(snapshot(temporary), before)
            # A caller-selected temporary root must not bypass the package guard.
            result = subprocess.run(
                [sys.executable, '-B', str(package / 'run_all.py')],
                cwd=temporary, capture_output=True, text=True, timeout=30,
                env={**os.environ, 'TMPDIR': str(package / 'receipts'),
                     'TEMP': str(package / 'receipts'), 'TMP': str(package / 'receipts')},
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('outside the input package', result.stderr)
            self.assertEqual(snapshot(temporary), before)


if __name__ == '__main__':
    unittest.main()
