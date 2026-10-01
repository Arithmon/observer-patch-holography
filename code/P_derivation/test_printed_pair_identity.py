#!/usr/bin/env python3
"""Check exact output pairs and expose the archived approximate counterexample.

A converged runtime artifact that claims an exact closure pair must satisfy

    alpha_root = (P - phi) / sqrt(pi),    alpha_root = 1 / alpha_inv,

because ``P`` is defined by the outer equation ``P = phi + alpha*sqrt(pi)``.
The full report has a converged pair. The frozen compressed trunk instead
retains an 18-digit, 12-iteration approximate candidate with a nonzero
fixed-point residual. Its bytes are source-pinned by later evidence; test
its residual accounting and non-promotion rather than falsely certifying it.

Stated tolerance: relative defect <= 1e-30 (at least 30 significant digits).
"""

from __future__ import annotations

from decimal import Decimal, localcontext
import json
from pathlib import Path

from paper_math import decimal_pi


RUNTIME = Path(__file__).resolve().parent / "runtime"
TRUNK = RUNTIME / "p_closure_trunk_current.json"
FULL_REPORT = RUNTIME / "full_p_alpha_report_current.json"

RELATIVE_TOLERANCE = Decimal("1e-30")
WORK_PRECISION = 160


def _constants() -> tuple[Decimal, Decimal]:
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        pi = +decimal_pi(WORK_PRECISION)
        sqrt_pi = pi.sqrt()
        phi = (Decimal(1) + Decimal(5).sqrt()) / Decimal(2)
        return +phi, +sqrt_pi


def _relative_identity_defect(p: str, alpha_inv: str) -> Decimal:
    phi, sqrt_pi = _constants()
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        alpha_root = Decimal(1) / Decimal(alpha_inv)
        rhs = (Decimal(p) - phi) / sqrt_pi
        return +(abs(alpha_root - rhs) / alpha_root)


def test_decimal_pi_matches_independent_mpmath_pi() -> None:
    """The chain's Chudnovsky pi agrees with an independent mpmath pi."""
    from mpmath.ctx_mp import MPContext

    mp = MPContext()
    mp.dps = 140
    independent = mp.nstr(mp.pi, 120, strip_zeros=False)
    ours = format(decimal_pi(150), "f")[:121]
    assert ours.startswith(independent[:118])


def test_archived_trunk_discloses_approximate_pair_and_matching_residual() -> None:
    payload = json.loads(TRUNK.read_text(encoding="utf-8"))
    fixed_point = payload["fixed_point_candidate"]
    defect = _relative_identity_defect(fixed_point["P"], fixed_point["alpha_inv"])
    assert payload["claim_status"] == "compressed_candidate_trunk_not_final_particle_root"
    assert payload["consumer_policy"]["may_feed_live_particle_predictions"] is False
    assert payload["source_report_precision"] == 18
    assert Decimal("5e-6") < defect < Decimal("6e-6")

    # The nonzero map-minus-probe residual must explain the pair mismatch.
    # 1e-28 covers the decimal display rounding of this archived artifact;
    # it is not an accepted fixed-point error or a certification tolerance.
    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        alpha = Decimal(fixed_point["alpha"])
        alpha_root = Decimal(1) / Decimal(fixed_point["alpha_inv"])
        residual = Decimal(fixed_point["alpha_fixed_point_residual"])
        assert Decimal("3e-8") < residual < Decimal("4e-8")
        assert abs((alpha_root - alpha) - residual) < Decimal("1e-28")
        assert abs(defect - residual / alpha_root) < Decimal("1e-25")
        # Check the archive's own constants at one unit of their last printed
        # decimal place, without promoting those displays to exact constants.
        phi, sqrt_pi = _constants()
        for key, independent in (("phi", phi), ("sqrt_pi", sqrt_pi)):
            displayed = Decimal(payload["closed_form_candidate"][key])
            last_place = Decimal(1).scaleb(displayed.as_tuple().exponent)
            assert last_place <= Decimal("1e-28")
            assert abs(displayed - independent) <= last_place


def test_full_report_printed_pair_identity() -> None:
    payload = json.loads(FULL_REPORT.read_text(encoding="utf-8"))
    defect = _relative_identity_defect(payload["p"], payload["alpha_inv"])
    assert defect <= RELATIVE_TOLERANCE

    with localcontext() as ctx:
        ctx.prec = WORK_PRECISION
        residual = abs(Decimal(payload["alpha_fixed_point_residual"]))
        alpha = Decimal(payload["alpha"])
        assert residual / alpha <= RELATIVE_TOLERANCE
