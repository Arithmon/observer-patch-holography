"""A scientifically broken implementation must fail an ordinary assertion.

Each mutant runs in an isolated copy. Import errors and collection failures
do not count as successful detection. These controls run in both CI systems.
"""
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


@pytest.mark.parametrize("module,before,after,test", [
    ("correctable_public_record_capacity.py", "if error > epsilon:",
     "if error > epsilon + Fraction(1, 10**12):",
     "test_checkpoint_channel_precision.py::test_error_threshold_is_not_enlarged_by_an_acceptance_tolerance"),
    ("checkpoint_channels.py", "probability(p) for out, p in row.items()",
     "probability(float(p)) for out, p in row.items()",
     "test_checkpoint_channel_precision.py::test_every_positive_binary_error_collapses_zero_error_capacity"),
    ("correctable_public_record_capacity.py", "if p > 0} for source, row in rows.items()",
     "if p > Fraction(1, 10**12)} for source, row in rows.items()",
     "test_checkpoint_channel_precision.py::test_every_positive_binary_error_collapses_zero_error_capacity"),
    ("correctable_public_record_capacity.py", "return directed_float(bound, upward=True)",
     "return float(bound)", "test_checkpoint_channel_precision.py::test_tv_error_bound_never_rounds_down"),
    ("correctable_public_record_capacity.py",
     "if any(p > 0 and target not in universe for row in rows.values() for target, p in row.items()):",
     "if False:",
     "test_checkpoint_channel_precision.py::test_indefinite_continuation_cannot_discard_positive_escape_mass"),
    ("checkpoint_channels.py", "error+(p if i != owner else 0)", "error+0",
     "test_checkpoint_decoder_certificates.py::test_every_three_input_three_output_half_count_channel"),
    ("verify_checkpoint_decoder.py", "return False", "return True",
     "test_checkpoint_decoder_certificates.py::test_independent_replay_rejects_forged_certificates"),
    ("checkpoint_channels.py", "out: p/mass for out, p in exact.items()",
     "out: p for out, p in exact.items()",
     "test_checkpoint_audit_boundaries.py::test_exact_row_mass_acceptance_boundary_and_disclosure"),
    ("checkpoint_channels.py", "for owner in reversed(owners):", "for owner in reversed(owners[:1]):",
     "test_checkpoint_audit_boundaries.py::test_unequal_denominators_at_both_sides_of_compound_optimum"),
    ("verify_checkpoint_decoder.py", "for source, row in raw.items():",
     "for source, row in ((x, r) for x, r in raw.items() if x in code):",
     "test_checkpoint_audit_boundaries.py::test_invalid_rows_outside_the_certified_code_are_still_rejected"),
])
def test_semantic_implementation_mutations_are_detected(tmp_path, module, before, after, test):
    source = Path(__file__).parent
    for name in ["correctable_public_record_capacity.py", "checkpoint_channels.py",
                 "public_record_csp.py", "verify_checkpoint_decoder.py",
                 "test_correctable_public_record_capacity.py",
                 "test_checkpoint_channel_precision.py", "test_checkpoint_decoder_certificates.py",
                 "test_checkpoint_audit_boundaries.py"]:
        shutil.copy2(source/name, tmp_path/name)
    target = tmp_path/module
    original = target.read_text(encoding="utf-8")
    assert before in original
    target.write_text(original.replace(before, after), encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", test],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
    )
    output = completed.stdout+completed.stderr
    assert completed.returncode == 1, output
    assert "failed" in output and "ERROR" not in output, output
    assert "AssertionError" in output or "DID NOT RAISE" in output, output
