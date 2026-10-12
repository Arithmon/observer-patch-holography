"""Independent source binding, fixed-angle enclosures and fail-closed controls."""
from fractions import Fraction
import json
from pathlib import Path

import mpmath
from mpmath.ctx_mp import MPContext
import numpy as np
import pytest

from native_response import certificate
from native_response.certificate import IntervalCalculation, certify
from native_response.memory import closure_data


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def fixed_certificates():
    """Actually recompute both cases while unrelated contexts are poisoned."""
    old_mp, old_iv = mpmath.mp.dps, mpmath.iv.dps
    mpmath.mp.dps, mpmath.iv.dps = 7, 6
    try:
        rows = {d: certify(d) for d in (16, 8)}
        assert (mpmath.mp.dps, mpmath.iv.dps) == (7, 6)
    finally:
        mpmath.mp.dps, mpmath.iv.dps = old_mp, old_iv
    return rows


def bounds(row):
    return (Fraction(row["lower_numerator"], row["denominator"]),
            Fraction(row["upper_numerator"], row["denominator"]))


def test_both_fixed_angles_have_recomputed_positive_enclosures(fixed_certificates):
    lower = {16: 186661676, 8: 3778430107}
    for denominator, row in fixed_certificates.items():
        assert row["theta_numerator"] == 1
        assert row["theta_denominator"] == denominator
        assert row["representative_edges"] == [[0, 2], [0, 4], [0, 8]]
        assert row["off_diagonal_multiplicities"] == [4, 4]
        assert row["taylor_degree"] == 24
        assert row["entry_tail_denominator"] == 2**50
        assert row["interval_decimal_precision"] == 35
        assert row["enclosures"]["gershgorin_margin"]["lower_numerator"] == lower[denominator]
        assert bounds(row["enclosures"]["gershgorin_margin"])[0] > 0
        # This direct thirty-seam numerical calculation imports no interval
        # producer and evaluates all overlaps rather than three representatives.
        data = closure_data(1 / denominator)
        gram = data["edge_gram"]
        values = {"beta": 1 - data["alpha"], "diagonal": gram[0, 0],
                  "cross_first": gram[0, 1], "cross_second": gram[0, 3],
                  "gershgorin_margin": np.min(np.diag(gram)
                      - np.sum(abs(gram - np.diag(np.diag(gram))), axis=1))}
        for name, value in values.items():
            lo, hi = bounds(row["enclosures"][name])
            assert float(lo) <= value <= float(hi)
        assert json.loads(json.dumps(row)) == row


@pytest.mark.parametrize("value", [True, False, None, 8.0, "8", 0, -8, 1, 9, 32,
                                   np.array([8]), np.ma.masked])
def test_certificate_rejects_undeclared_or_malformed_angles(value):
    with pytest.raises(ValueError, match="only theta"):
        certify(value)


def test_interval_source_contains_independent_authoritative_entries():
    """Decode exact recorded coefficients using a separate 100-digit oracle."""
    packet = json.loads((ROOT / "code/source_selection_model/response.json").read_text())
    calculation = IntervalCalculation()
    oracle = MPContext()
    oracle.dps = 100
    root = oracle.sqrt(5)

    def rational(value):
        number = Fraction(value)
        return oracle.mpf(number.numerator) / number.denominator

    for actual, recorded in zip(calculation.generators, packet["generators"]):
        for block, raw in zip(actual, recorded):
            for i in range(3):
                for j in range(3):
                    a, b, c, d = map(rational, raw[i][j])
                    for interval, expected in ((block[i, j].real, a + b * root),
                                               (block[i, j].imag, c + d * root)):
                        assert oracle.mpf(interval.a) <= expected <= oracle.mpf(interval.b)


@pytest.mark.parametrize("interval, expected", [
    (("0.1234567890122", "0.1234567890123"), (123456789012, 123456789013)),
    (("-0.1234567890123", "-0.1234567890122"), (-123456789013, -123456789012)),
    (("-0.0000000000001", "0.0000000000001"), (-1, 1)),
    (("0", "0"), (0, 0)),
])
def test_serialization_rounds_outward(interval, expected):
    calculation = IntervalCalculation()
    result = calculation.rational_enclosure(calculation.iv.mpf(interval))
    assert (result["lower_numerator"], result["upper_numerator"]) == expected
    assert result["denominator"] == 10**12


def test_declared_precision_contexts_are_isolated():
    old_mp, old_iv = mpmath.mp.dps, mpmath.iv.dps
    try:
        mpmath.mp.dps, mpmath.iv.dps = 4, 4
        calculation = IntervalCalculation()
        assert calculation.iv.dps == 35
        assert calculation.mp.dps == 80
        assert (mpmath.mp.dps, mpmath.iv.dps) == (4, 4)
    finally:
        mpmath.mp.dps, mpmath.iv.dps = old_mp, old_iv


def test_no_leakage_control_cannot_pass_positivity(monkeypatch):
    def zero_local(self, edge, theta):
        return self.iv.mpf(0), [[self.iv.zeros(6), self.iv.zeros(6)] for _ in (0, 1)]

    monkeypatch.setattr(IntervalCalculation, "local", zero_local)
    with pytest.raises(ArithmeticError, match="positive full-network"):
        certify(8)


@pytest.mark.parametrize("name, value", [("DEGREE", 0), ("TAIL_DENOMINATOR", 2**100)])
def test_unproved_taylor_tail_is_rejected(monkeypatch, name, value):
    calculation = IntervalCalculation()
    monkeypatch.setattr(certificate, name, value)
    with pytest.raises(ArithmeticError, match="Taylor tail"):
        calculation.exponential(calculation.iv.zeros(18), calculation.iv.mpf(1) / 8)


def test_unproved_reference_norm_is_rejected():
    calculation = IntervalCalculation()
    with pytest.raises(ArithmeticError, match="Taylor norm"):
        calculation.exponential(17 * calculation.iv.j * calculation.iv.eye(18),
                                calculation.iv.mpf(1) / 8)


def test_small_row_norm_does_not_hide_a_large_column_norm():
    calculation = IntervalCalculation()
    generator = calculation.iv.zeros(18)
    for i in range(18):
        generator[i, 0] = 8
    with pytest.raises(ArithmeticError, match="Taylor norm"):
        calculation.exponential(generator, calculation.iv.mpf(1) / 8)


def test_polynomial_exponential_encloses_an_exact_zero_generator():
    calculation = IntervalCalculation()
    result = calculation.exponential(calculation.iv.zeros(18), calculation.iv.mpf(1) / 8)
    tail = calculation.iv.mpf(1) / 2**50
    for i in range(18):
        for j in range(18):
            assert result[i, j].real.a <= int(i == j) - tail
            assert result[i, j].real.b >= int(i == j) + tail
            assert result[i, j].imag.a <= -tail
            assert result[i, j].imag.b >= tail


def test_polynomial_exponential_contains_independent_nonzero_oracle():
    calculation = IntervalCalculation()
    iv = calculation.iv
    generator = iv.j * iv.sqrt(5) * iv.eye(18)
    result = calculation.exponential(generator, iv.mpf(1) / 8)
    oracle = MPContext()
    oracle.dps = 100
    diagonal = oracle.exp(oracle.j * oracle.sqrt(5) / 8)
    for i in range(18):
        for j in range(18):
            expected = diagonal if i == j else oracle.mpc(0)
            assert oracle.mpf(result[i, j].real.a) <= expected.real <= oracle.mpf(result[i, j].real.b)
            assert oracle.mpf(result[i, j].imag.a) <= expected.imag <= oracle.mpf(result[i, j].imag.b)
