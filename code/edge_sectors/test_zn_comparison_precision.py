"""Keep unresolved Z_n cancellation out of a successful held-out report."""

from pathlib import Path
import re
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge_sectors import heat_kernel_holdout_validation as reporting
from edge_sectors.test_electric_sector_controls import _original_incidence_oracle


def test_nearly_uniform_comparison_refuses_unresolved_mismatch_sign():
    h = 5.708163265306122e-9
    ctx, reference = _original_incidence_oracle(5, h, dps=65, return_high_precision=True)
    eigenvalues = [4*ctx.sin(ctx.pi*q/5)**2 for q in range(5)]
    fitted = ctx.log(reference[0]/reference[1])/eigenvalues[1]
    predicted_gap = fitted*eigenvalues[2]
    mismatch = (predicted_gap-ctx.log(reference[0]/reference[2]))/abs(predicted_gap)
    assert float(mismatch) == pytest.approx(-1.09016217468756855e-9, rel=1e-13, abs=0)
    # On reviewed head 8bf11dbb these accurate binary64 probabilities led to
    # the opposite-sign +8.867257980040205e-9 mismatch and printed +0.00%.
    # The fitted log gap exceeded the existing 1e-8 gate. The independent
    # probability observable must remain available when its comparison fails.
    measured = reporting.zn_edge_distribution(5, h)
    assert np.all(measured > 0)
    assert measured == pytest.approx([float(p) for p in reference], rel=2e-14, abs=0)
    with pytest.raises(ValueError, match="unresolved"):
        reporting.report_zn(5, [h])


def test_ordinary_overconstrained_comparison_remains_available(capsys):
    reporting.report_zn(5, [.05])
    output = capsys.readouterr().out
    assert "q=2:" in output
    assert "res -0.96%" in output


def test_small_resolved_comparison_does_not_print_zero_percent(capsys):
    reporting.report_zn(5, [1e-5])
    output = capsys.readouterr().out
    assert re.search(r"res -[1-9]\.\d{3}e-\d{2}%", output)
    assert "0.00%" not in output


def test_roundoff_sized_mismatch_is_unresolved_for_supplied_heat_weights(monkeypatch):
    # These positive conjugate Z4 weights have exact binary ratios 2 and 4.
    # With the exact Laplacian ratio 2, the held-out mismatch vanishes.
    # Float sine/log arithmetic cannot certify that zero from its residual.
    probabilities = np.array([4/9, 2/9, 1/9, 2/9])
    assert probabilities[0] == 2*probabilities[1] == 4*probabilities[2]
    monkeypatch.setattr(reporting, "zn_edge_distribution", lambda n, h: probabilities)
    with pytest.raises(ValueError, match="unresolved"):
        reporting.report_zn(4, [1.])
