from __future__ import annotations

import json
import math
import cmath
from fractions import Fraction
from pathlib import Path

import pytest

from verify_collar_gap_certificate import verify


HERE = Path(__file__).resolve().parent


def test_exact_rational_contract_witness() -> None:
    payload = json.loads(
        (HERE / "certificates" / "issue_306_theorem_contract_witness.json").read_text()
    )
    result = verify(payload)
    assert result["valid"] is True
    assert result["physical_clay_receipt"] is False
    assert result["c_floor"] == "3/4"
    assert result["eta_upper"] == "1/2"
    assert result["gap_lower"] == "3/8"


def test_no_mixing_local_countermodel_has_zero_gap() -> None:
    # On support {00, 11}, conditioning either spin on the other fixes it.
    states = [(0, 0), (1, 1)]
    law = [Fraction(1, 2), Fraction(1, 2)]
    generator = [[Fraction(0) for _ in states] for _ in states]
    for site in range(2):
        for i, state in enumerate(states):
            fiber = [j for j, other in enumerate(states) if other[1-site] == state[1-site]]
            mass = sum(law[j] for j in fiber)
            generator[i][i] += 1
            for j in fiber:
                generator[i][j] -= law[j]/mass
    centered_sign = [-1, 1]
    assert sum(p*f for p, f in zip(law, centered_sign)) == 0
    assert sum(p*f*f for p, f in zip(law, centered_sign)) == 1
    assert [sum(a*f for a, f in zip(row, centered_sign)) for row in generator] == [0, 0]

    # Exact conditional rows have TV=1; the positive-gap gate must refuse.
    payload = json.loads((HERE / "certificates" / "issue_306_theorem_contract_witness.json").read_text())
    for entry in payload["types"]:
        entry["influences"] = [{"target_type": entry["id"], "upper": "1",
                                "conditional_rows": [["1", "0"], ["0", "1"]]}]
    payload.pop("expected")
    with pytest.raises(ValueError, match="< 1"):
        verify(payload)


def test_product_mixing_nonlocal_gray_cycle_gap_vanishes() -> None:
    gaps = []
    for m in range(3, 9):
        size = 2**m
        mode = [cmath.exp(2j * math.pi * index / size) for index in range(size)]
        # (I-E_0)+(I-E_1) is one half of the cycle graph Laplacian.
        applied = [
            value - (mode[(index - 1) % size] + mode[(index + 1) % size]) / 2
            for index, value in enumerate(mode)
        ]
        gap = 1 - math.cos(2 * math.pi / size)
        assert max(abs(lhs - gap * rhs) for lhs, rhs in zip(applied, mode)) < 1e-12
        gaps.append(gap)
    assert all(later < earlier for earlier, later in zip(gaps, gaps[1:]))
    assert gaps[-1] < 0.001
