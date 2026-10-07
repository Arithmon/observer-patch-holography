"""Independent replay of a code's attainability, without importing its producer.

This checks explicit decoders against the supplied channels, not the claim
that no larger code exists. Maximality is established by the exact search;
independent exhaustive integer-channel controls test that separate claim.
"""
from collections.abc import Mapping, Sequence
from decimal import Decimal
from fractions import Fraction
import math


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float, Fraction, Decimal)):
        raise ValueError("not a numerical probability")
    value = Fraction(value)
    if not 0 <= value <= 1:
        raise ValueError("not a probability")
    return value


def _wire(value):
    if not isinstance(value, dict) or set(value) != {"numerator", "denominator"}:
        raise ValueError("not an exact ratio")
    n, d = value["numerator"], value["denominator"]
    if type(n) is not str or type(d) is not str:
        raise ValueError("ratio components must be decimal integer strings")
    result = Fraction(int(n), int(d))
    if str(result.numerator) != n or str(result.denominator) != d:
        raise ValueError("ratio is not canonical")
    return result


def _labels(values):
    if (not isinstance(values, Sequence) or isinstance(values, (str, bytes))
            or not values or any(type(x) is not str or not x for x in values)
            or len(set(values)) != len(values)):
        raise ValueError("invalid alphabet")
    return set(values)


def verify_decoder_witness(reachable, channels, epsilon, result):
    """Return False for malformed, unbound, or scientifically false witnesses."""
    try:
        universe = _labels(reachable)
        tolerance = _number(epsilon)
        if not isinstance(result, Mapping):
            return False
        code = _labels(result["code_witness"])
        if (not code <= universe or type(result["capacity"]) is not int
                or result["capacity"] != len(code)
                or _wire(result["epsilon_exact"]) != tolerance
                or result["probability_contract"] != "exact_relative_weights_with_disclosed_near_unit_row_masses"
                or result["certificate_scope"] != "deterministic_code_attainability"):
            return False
        if (not isinstance(channels, Sequence) or isinstance(channels, (str, bytes))
                or not channels):
            return False
        certificates = result["decoder_certificates"]
        successes = result["worst_input_success_by_channel"]
        errors = result["worst_input_error_by_channel"]
        if any(not isinstance(x, list) or len(x) != len(channels)
               for x in (certificates, successes, errors)):
            return False
        for channel, certificate, success_bound, error_bound in zip(channels, certificates, successes, errors):
            if not isinstance(channel, Mapping) or not isinstance(certificate, Mapping):
                return False
            raw = channel["rows"]
            if not isinstance(raw, Mapping) or set(raw) != universe:
                return False
            rows, masses = {}, {}
            for source, row in raw.items():
                if (not isinstance(row, Mapping) or not row
                        or any(type(y) is not str or not y for y in row)):
                    return False
                cells = {y: _number(p) for y, p in row.items()}
                mass = sum(cells.values(), Fraction(0))
                if mass == 0 or abs(mass-1) > Fraction(1, 10**12):
                    return False
                rows[source] = {y: p/mass for y, p in cells.items()}
                masses[source] = mass
            decoder = certificate["decoder"]
            outputs = {y for x in code for y in rows[x]}
            if (not isinstance(decoder, Mapping) or set(decoder) != outputs
                    or any(type(x) is not str or x not in code for x in decoder.values())):
                return False
            exact_errors = {x: sum((p for y, p in rows[x].items() if decoder[y] != x), Fraction(0))
                            for x in code}
            reported = certificate["error_by_input"]
            reported_masses = certificate["input_row_masses"]
            if not isinstance(reported, Mapping) or not isinstance(reported_masses, Mapping):
                return False
            if set(reported) != code or set(reported_masses) != code:
                return False
            if any(_wire(reported[x]) != exact_errors[x] or _wire(reported_masses[x]) != masses[x] for x in code):
                return False
            worst = max(exact_errors.values())
            if worst > tolerance or _wire(certificate["worst_input_error"]) != worst:
                return False
            if not math.isfinite(success_bound) or not math.isfinite(error_bound):
                return False
            if _number(success_bound) > 1-worst or _number(error_bound) < worst:
                return False
        return True
    except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError):
        return False
