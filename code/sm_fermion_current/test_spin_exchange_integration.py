"""Keep the current's statistics boundary in its mandatory test collection."""
from pathlib import Path
import subprocess
import sys


def test_independent_spin_exchange_certificate_and_mutation_controls():
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q",
         "code/spin_exchange", "code/pauli_stability", "code/pauli_source_selection"],
        cwd=root, text=True, capture_output=True, timeout=300,
    )
    assert result.returncode == 0, result.stdout + result.stderr
