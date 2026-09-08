"""Serialize finite mpmath interval endpoints with exact directed rounding.

Binary endpoints are converted as integers, without passing through the
independently configured scalar mpmath context. Decimal division then rounds
the lower endpoint down and the upper endpoint up.
"""

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from mpmath.libmp import to_rational


def endpoint_decimal(endpoint: tuple, digits: int, *, upper: bool) -> str:
    if type(digits) is not int or digits < 1:
        raise ValueError("decimal precision must be a positive integer")
    # mpmath encodes NaN and both infinities with negative bit counts.
    # to_rational itself rejects NaN only; infinities would become zero.
    if endpoint[3] < 0:
        raise ValueError("interval endpoints must be finite")
    numerator, denominator = to_rational(endpoint)
    with localcontext() as context:
        context.prec = digits
        context.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        value = Decimal(numerator) / Decimal(denominator)
    return str(value)


def interval_json(value, digits: int = 40) -> dict[str, str]:
    lower, upper = value._mpi_
    return {
        "lo": endpoint_decimal(lower, digits, upper=False),
        "hi": endpoint_decimal(upper, digits, upper=True),
    }


def scalar_bound(value, digits: int, *, upper: bool) -> str:
    """Round an already computed finite binary bound in its declared direction."""
    return endpoint_decimal(value._mpf_, digits, upper=upper)
