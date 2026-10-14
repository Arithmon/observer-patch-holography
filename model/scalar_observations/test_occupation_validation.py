"""Reject inputs which could otherwise turn a skipped evolution into a certificate."""
from fractions import Fraction as F

import pytest

from .occupation_exact import Generator, Vector, certificate, flow, polynomial_step, quantize


@pytest.fixture
def blank():
    return Vector({(0,)*7: (1, 0)})


@pytest.mark.parametrize("substeps", [0, -1, True, 2.0])
def test_flow_cannot_silently_skip_requested_evolution(blank, substeps):
    with pytest.raises(ValueError, match="substeps"):
        flow(Generator(4), blank, F(1, 8), substeps=substeps)


@pytest.mark.parametrize("duration", [0.125, True, "1/8", float("nan")])
def test_residual_certificate_requires_exact_duration(blank, duration):
    for operation in (flow, polynomial_step):
        with pytest.raises(ValueError, match="exact rational"):
            operation(Generator(4), blank, duration)


@pytest.mark.parametrize("degree", [0, -1, True, 1.5])
def test_polynomial_degree_must_be_positive_integer(blank, degree):
    with pytest.raises(ValueError, match="degree"):
        polynomial_step(Generator(4), blank, F(1, 8), degree=degree)


@pytest.mark.parametrize("bits", [0, -1, True, 80.0])
def test_rounding_precision_must_be_positive_integer(blank, bits):
    with pytest.raises(ValueError, match="rounding bits"):
        quantize(blank, bits)


@pytest.mark.parametrize("dimensions", [None, True, 1, (4,)*6, (4,)*6+(False,)])
def test_invalid_occupation_box_is_rejected(dimensions):
    with pytest.raises(ValueError, match="occupation dimensions"):
        Generator(dimensions)


@pytest.mark.parametrize("kwargs", [
    {"coupling": True}, {"coupling": 2}, {"sign": 0}, {"sign": True},
    {"interval": -F(1, 8)}, {"interval": 0.125},
    {"prefinal_only": 1}, {"return_terminal_state": "yes"},
])
def test_history_certificate_rejects_undeclared_parameters(kwargs):
    with pytest.raises(ValueError):
        certificate(4, **kwargs)


def test_zero_duration_is_exact_identity_with_zero_budget(blank):
    output, budget = flow(Generator(4), blank, F(0), degree=2, substeps=1)
    assert output.coefficients == {(0,)*7: (output.denominator, 0)}
    assert all(value == 0 for value in budget.values())
