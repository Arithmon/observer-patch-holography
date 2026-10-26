"""Rational enclosure of an interaction-generated Hubbard charge threshold.

This is an alternative coherent architecture, not the finite dimerized walk.
The native preparation fixes two dimensionless coefficients; assigning them
to a homogeneous half-filled Hubbard chain is an explicit architecture choice.
Lieb--Wu supplies the thermodynamic theorem. No empirical mass, quadrature,
asymptotic equality, fitted chemical potential or relaxation rate enters the
certificate. The elementary enclosure uses only standard-library arithmetic.

Primary sources:
  Lieb and Wu, https://arxiv.org/html/cond-mat/0207529, section 7, eqs. 66--67.
  Jeckelmann, Gebhard and Essler, https://arxiv.org/html/cond-mat/9911281,
  eq. 2 (exact gap), eqs. 6--8 (the separate low-energy field-theory limit).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import json
import math


PRECISION_BITS = 160


def _bits(value: int) -> int:
    if type(value) is not int or not 64 <= value <= 512:
        raise ValueError("precision_bits must be an integer from 64 to 512")
    return value


def _rational(value: int | F) -> F:
    if type(value) not in (int, F):
        raise TypeError("exact integer or Fraction input required")
    return F(value)


@dataclass(frozen=True)
class RationalInterval:
    """Closed rational interval, rounded outward onto a fixed dyadic grid."""

    lower: F
    upper: F
    precision_bits: int = PRECISION_BITS

    def __post_init__(self) -> None:
        lo, hi, bits = _rational(self.lower), _rational(self.upper), _bits(self.precision_bits)
        if lo > hi:
            raise ValueError("reversed interval")
        scale = 1 << bits
        lower = F((lo.numerator * scale) // lo.denominator, scale)
        upper = F(-((-hi.numerator * scale) // hi.denominator), scale)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @classmethod
    def enclose(cls, value: int | F, precision_bits: int = PRECISION_BITS) -> RationalInterval:
        value = _rational(value)
        return cls(value, value, precision_bits)

    def _coerce(self, value: RationalInterval | int | F) -> RationalInterval:
        if isinstance(value, RationalInterval):
            if value.precision_bits != self.precision_bits:
                raise ValueError("interval precisions must agree")
            return value
        return RationalInterval.enclose(value, self.precision_bits)

    def __add__(self, other: RationalInterval | int | F) -> RationalInterval:
        other = self._coerce(other)
        return RationalInterval(self.lower + other.lower, self.upper + other.upper, self.precision_bits)

    __radd__ = __add__

    def __neg__(self) -> RationalInterval:
        return RationalInterval(-self.upper, -self.lower, self.precision_bits)

    def __sub__(self, other: RationalInterval | int | F) -> RationalInterval:
        return self + (-self._coerce(other))

    def __rsub__(self, other: RationalInterval | int | F) -> RationalInterval:
        return self._coerce(other) - self

    def __mul__(self, other: RationalInterval | int | F) -> RationalInterval:
        other = self._coerce(other)
        products = tuple(a * b for a in (self.lower, self.upper)
                         for b in (other.lower, other.upper))
        return RationalInterval(min(products), max(products), self.precision_bits)

    __rmul__ = __mul__

    def __truediv__(self, other: RationalInterval | int | F) -> RationalInterval:
        other = self._coerce(other)
        if other.lower <= 0 <= other.upper:
            raise ZeroDivisionError("denominator interval contains zero")
        reciprocal = RationalInterval(1 / other.upper, 1 / other.lower, self.precision_bits)
        return self * reciprocal

    def __rtruediv__(self, other: RationalInterval | int | F) -> RationalInterval:
        return self._coerce(other) / self

    def __pow__(self, power: int) -> RationalInterval:
        if type(power) is not int or power < 0:
            raise ValueError("nonnegative integer power required")
        result = RationalInterval.enclose(1, self.precision_bits)
        base = self
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def as_json(self) -> dict:
        """Exact strings are authoritative; float endpoints also round outward."""
        lo, hi = float(self.lower), float(self.upper)
        if not math.isfinite(lo) or not math.isfinite(hi):
            raise ValueError("finite float display endpoints required")
        if F.from_float(lo) > self.lower:
            lo = math.nextafter(lo, -math.inf)
        if F.from_float(hi) < self.upper:
            hi = math.nextafter(hi, math.inf)
        return {"lower_exact": str(self.lower), "upper_exact": str(self.upper),
                "lower": lo, "upper": hi,
                "midpoint": float((self.lower + self.upper) / 2),
                "width": float(self.upper - self.lower)}


def log_bounds(value: int | F, *, precision_bits: int = PRECISION_BITS) -> RationalInterval:
    """Use log(x)=2 atanh((x-1)/(x+1)) and a positive geometric tail."""
    x, bits = _rational(value), _bits(precision_bits)
    if x <= 0:
        raise ValueError("positive logarithm argument required")
    if x < 1:
        return -log_bounds(1 / x, precision_bits=bits)
    q = (x - 1) / (x + 1)
    q2, power, total = q * q, q, F(0)
    for k in range(bits):
        total += power / (2 * k + 1)
        power *= q2
    tail = 2 * power / ((2 * bits + 1) * (1 - q2))
    return RationalInterval(2 * total, 2 * total + tail, bits)


def pi_bounds(*, precision_bits: int = PRECISION_BITS) -> RationalInterval:
    """Machin's identity with exact alternating-arctangent remainders."""
    bits = _bits(precision_bits)
    terms = bits // 3 + 16

    def atan(q: F) -> RationalInterval:
        power, total = q, F(0)
        for k in range(terms):
            total += (-1 if k & 1 else 1) * power / (2 * k + 1)
            power *= q * q
        next_term = (-1 if terms & 1 else 1) * power / (2 * terms + 1)
        return RationalInterval(min(total, total + next_term), max(total, total + next_term), bits)

    return 16 * atan(F(1, 5)) - 4 * atan(F(1, 239))


def sqrt_bounds(value: RationalInterval) -> RationalInterval:
    """Enclose endpoint square roots using integer isqrt only."""
    if value.lower < 0:
        raise ValueError("nonnegative square-root interval required")
    scale = 1 << value.precision_bits

    def floor_root(x: F) -> F:
        return F(math.isqrt((x.numerator * scale * scale) // x.denominator), scale)

    lo, hi = floor_root(value.lower), floor_root(value.upper)
    if hi * hi != value.upper:
        hi += F(1, scale)
    return RationalInterval(lo, hi, value.precision_bits)


def exp_negative_bounds(value: RationalInterval) -> RationalInterval:
    """Reciprocal positive Taylor sums, range reduction, and exact tail bounds."""
    if value.lower < 0:
        raise ValueError("nonnegative exponent interval required")
    bits, subdivisions = value.precision_bits, 0
    while value.upper > (1 << subdivisions):
        subdivisions += 1
    factor = 1 << subdivisions
    terms = bits // 3 + 16

    def scalar(x: F) -> tuple[F, F]:
        total, term = F(1), F(1)
        for k in range(1, terms + 1):
            term *= x / k
            total += term
        first_omitted = term * x / (terms + 1)
        ratio_upper = x / (terms + 2)
        tail = first_omitted / (1 - ratio_upper)
        return 1 / (total + tail), 1 / total

    lo = scalar(value.upper / factor)[0]
    hi = scalar(value.lower / factor)[1]
    result = RationalInterval(lo, hi, bits)
    for _ in range(subdivisions):
        result = result * result
    return result


def _k1_bounds(a: RationalInterval, exp_minus_a: RationalInterval, multiple: int,
               even_order: int, pi: RationalInterval) -> tuple[RationalInterval, RationalInterval]:
    z = multiple * a
    reciprocal = 1 / z
    power = RationalInterval.enclose(1, a.precision_bits)
    total, coefficient = power, F(1)
    lower = total
    for k in range(1, even_order + 2):
        coefficient *= F(4 - (2 * k - 1) ** 2, 8 * k)
        power = power * reciprocal
        total = total + coefficient * power
        if k == even_order:
            lower = total
    base = sqrt_bounds(pi / (2 * z)) * exp_minus_a ** multiple
    return base * lower, base * total


def hubbard_charge_gap_bounds(hopping: RationalInterval, repulsion: RationalInterval,
                              *, bessel_even: int = 8) -> RationalInterval:
    """Bound the exact thermodynamic half-filled 1D charge gap.

    This quantity is mu_plus-mu_minus, twice the symmetric single-holon
    rest-energy threshold. It is not the global many-body gap or an electron
    pole mass. The zero-repulsion return is the exact gapless limit.
    """
    if hopping.precision_bits != repulsion.precision_bits:
        raise ValueError("interval precisions must agree")
    if hopping.lower <= 0:
        raise ValueError("strictly positive hopping required")
    if type(bessel_even) is not int or bessel_even % 2 or not 0 <= bessel_even <= 32:
        raise ValueError("even Bessel bound order from 0 to 32 required")
    if repulsion.lower == repulsion.upper == 0:
        return RationalInterval.enclose(0, hopping.precision_bits)
    if repulsion.lower <= 0:
        raise ValueError("strictly positive repulsion or exact zero required")
    pi = pi_bounds(precision_bits=hopping.precision_bits)
    a = 2 * pi * hopping / repulsion
    exponential = exp_negative_bounds(a)
    if (1 - exponential ** 2).lower <= 0:
        raise ValueError("precision cannot resolve the positive exponential-tail denominator")
    lo1, hi1 = _k1_bounds(a, exponential, 1, bessel_even, pi)
    lo3, hi3 = _k1_bounds(a, exponential, 3, bessel_even, pi)
    lower = 16 * hopping / pi * (lo1 + lo3 / 3)
    tail = hi1 * exponential ** 4 / (1 - exponential ** 2)
    upper = 16 * hopping / pi * (hi1 + hi3 / 3 + tail)
    return RationalInterval(max(F(0), lower.lower), upper.upper, hopping.precision_bits)


def native_reference_data() -> dict:
    """Read all native source weights; never import fitted particle values."""
    from ...native_modular_reference import local_preparation
    from ...observer import native

    basis = native.states(binary=True)
    weights = local_preparation.reference()
    levels = [set() for _ in range(3)]
    multiplicities = [0, 0, 0]
    for state, weight in zip(basis, weights, strict=True):
        q = sum(state[p] for p in local_preparation.CORE)
        levels[q].add(F(weight))
        multiplicities[q] += 1
    if any(len(level) != 1 for level in levels):
        raise ArithmeticError("reference is not constant on each core-charge sector")
    w = tuple(next(iter(level)) for level in levels)
    if min(w) <= 0 or sum(n * p for n, p in zip(multiplicities, w, strict=True)) != 1:
        raise ArithmeticError("reference is not faithful and normalized")
    return {"weights": w, "multiplicities": tuple(multiplicities),
            "clock_ratio": w[2] / w[0], "connected_ratio": w[1] ** 2 / (w[0] * w[2])}


def certificate(*, precision_bits: int = PRECISION_BITS, bessel_even: int = 8) -> dict:
    """Return JSON-safe exact bounds and the theorem/architecture boundary."""
    source = native_reference_data()
    omega = log_bounds(source["clock_ratio"], precision_bits=precision_bits)
    repulsion = log_bounds(source["connected_ratio"], precision_bits=precision_bits)
    hopping = omega / 2
    gap = hubbard_charge_gap_bounds(hopping, repulsion, bessel_even=bessel_even)
    zero = hubbard_charge_gap_bounds(hopping, RationalInterval.enclose(0, precision_bits))
    return {
        "schema": "oph.mass_response.hubbard_charge_gap.v1",
        "source": {
            "reference": "model/native_modular_reference/local_preparation.reference",
            "weights_per_configuration": list(map(str, source["weights"])),
            "charge_sector_multiplicities": list(source["multiplicities"]),
            "clock_ratio": str(source["clock_ratio"]),
            "connected_ratio": str(source["connected_ratio"]),
            "empirical_mass_inputs": [],
            "ideal_native_reference_contract": True,
            "finite_preparation_error_propagated": False,
        },
        "omega": omega.as_json(), "hopping": hopping.as_json(),
        "repulsion": repulsion.as_json(), "repulsion_over_hopping": (repulsion / hopping).as_json(),
        "charge_gap": gap.as_json(), "gap_over_clock": (gap / omega).as_json(),
        "holon_energy": (gap / 2).as_json(), "holon_over_clock": (gap / (2 * omega)).as_json(),
        "zero_interaction_control": zero.as_json(),
        "architecture": {
            "law": "homogeneous repulsive spinful one-dimensional Hubbard chain",
            "hopping_assignment": "t=Omega/2; Omega=log(w2/w0)",
            "repulsion_assignment": "U=log(w1^2/(w0*w2))",
            "half_filling_and_ground_state_are_contract_inputs": True,
            "native_reference_is_not_the_Hubbard_vacuum": True,
            "no_bare_staggered_or_Dirac_mass": True,
            "architecture_assignment_forced_by_native_repair": False,
            "same_as_primary_finite_dimerized_experiment": False,
        },
        "proof": {
            "thermodynamic_input": "Lieb-Wu half-filled repulsive chain charge-gap theorem",
            "gap_formula": "16*t^2/U * integral_1^infinity sqrt(y^2-1)/sinh(2*pi*t*y/U) dy",
            "positive_series": "16*t/pi * sum_r>=0 K1((2*r+1)*a)/(2*r+1); a=2*pi*t/U",
            "precision_bits": precision_bits,
            "log": "finite atanh series with geometric tail",
            "pi": "Machin identity and alternating arctangent remainders",
            "sqrt": "integer dyadic endpoint bounds",
            "exp": "reciprocal positive Taylor sums, geometric remainder and repeated squaring",
            "K1": "global sqrt(1+x) Taylor lower/upper polynomials integrated exactly",
            "K1_lower_order": bessel_even, "K1_upper_order": bessel_even + 1,
            "explicit_positive_series_terms": [1, 3],
            "remaining_series_upper": "K1(a)*exp(-4*a)/(1-exp(-2*a))",
            "quadrature_used": False, "all_arithmetic_outward_rational": True,
            "strictly_positive_certified_lower_endpoint": gap.lower > 0,
        },
        "interpretation": {
            "quantity": "thermodynamic charge chemical-potential jump and optical threshold",
            "single_holon_rest_energy_is_half_charge_gap": True,
            "spin_sector_gapless": True, "positive_global_many_body_gap_claimed": False,
            "isolated_electron_pole_claimed": False,
            "relativistic_dispersion": "separate 1+1 low-energy/scaling-limit theorem only",
            "fixed_lattice_Lorentz_error_bound_supplied": False,
            "three_plus_one_identification": False, "observed_Standard_Model_mass": False,
            "native_relaxation_rate_identified_as_mass": False,
            "physical_clock_calibration_supplied": False,
            "finite_experiment_observed_thermodynamic_threshold": False,
        },
        "references": [
            "https://arxiv.org/html/cond-mat/0207529",
            "https://arxiv.org/html/cond-mat/9911281",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2, sort_keys=True, allow_nan=False))
