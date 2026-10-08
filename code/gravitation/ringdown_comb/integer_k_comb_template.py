"""Integer-k Kerr ringdown comb template. Build-stage numeric instrument.

Historical input (not modified here):
`falsification/frozen_targets/fz01_2026-07-17/frozen_target_integer_k_comb_2026-07-17.md`
with companion statement
`proof/epic_wins/ringdown_comb/INTEGER_K_COMB_STATEMENT.md`, under the
scientific-status erratum
`falsification/frozen_targets/fz01_2026-07-17/SCIENTIFIC_STATUS_ERRATUM_2026-07-29.md`.

Imported continuation law. Candidate spectral features above the rotation line
satisfy, under the comb hypothesis,

    (f_a - m*Omega_H/(2*pi)) / (f_b - m*Omega_H/(2*pi)) = ln(k_a)/ln(k_b)

for integers k >= 2, equivalently universal-coordinate positions
x_k = ln(k)/(8*pi) with x = (G*M/(c^3*g(chi))) * (omega - m*Omega_H).
Secondary structure: the within-line (k-1)/k KMS net-response factor and the
mass-independent but spin-dependent linewidth-to-spacing ratio
64*pi^2*p_0/(a*g(chi)^2*ln(k)) with declared a in [1, 10]. This factor follows
from P = p_0*hbar*c^6/(G^2*M^2), mean emitted energy
a*hbar*kappa/(2*pi), spacing hbar*kappa*ln(k)/(2*pi), and
g(chi)=4*G*M*kappa/c^3. The superseded expression omitted g(chi)^(-2).

Entropy sign and divisibility. For an emission with record dimensions
d_before and d_after, the integer-division premise is
d_before = k*d_after. It therefore requires k to divide d_before. With
S_BH = ln(d), the signed black-hole entropy change is
ln(d_after)-ln(d_before) = -ln(k); the positive entropy loss entering the
emitted-frequency formula is +ln(k). The rule that such division steps dominate
is an imported continuation premise, not a source-derived OPH theorem.

Kerr horizon functions (geometric input, SI output). With
s(chi) = sqrt(1 - chi^2):

    r_plus      = (G*M/c^2) * (1 + s(chi))
    Omega_H     = c^3 * chi     / (2*G*M*(1 + s(chi)))      [rad/s]
    kappa       = c^3 * s(chi)  / (2*G*M*(1 + s(chi)))      [1/s]

Derivation of g(chi) from the imported law's universal coordinate (the
KMS/temperature reading pinned by the companion statement). The leading
semiclassical first-law template of the statement is

    hbar*(omega - m*Omega_H) = k_B*T_H * (S_before-S_after)
                              = k_B*T_H * ln(k),

and with the Hawking temperature k_B*T_H = hbar*kappa/(2*pi) (kappa in
1/s as above) this reads

    omega - m*Omega_H = (kappa/(2*pi)) * ln(k).

The imported universal coordinate x = (G*M/(c^3*g(chi)))*(omega - m*Omega_H)
takes the value ln(k)/(8*pi) on these lines exactly when

    kappa/(2*pi) = (c^3*g(chi)/(G*M)) / (8*pi),

that is g(chi) = 4*G*M*kappa/c^3 = 2*sqrt(1-chi^2)/(1+sqrt(1-chi^2)),
which is the explicit g(chi) of the companion statement. This relation treats
the Kerr background as fixed during one transition; a finite-step completion
must account for its change. g(chi) is
therefore statement-pinned within this imported continuation, not a free
numeric selection; the receipt records the
identity g(chi) = 4*G*M*kappa/c^3 as its provenance. Note that the
alternative normalization g = G*M*kappa/c^3 is inconsistent with
x_k = ln(k)/(8*pi) by a factor of 4 and is not used.

Factor-of-2*pi bookkeeping for tooth frequencies in Hz. With
omega = 2*pi*f,

    omega_k - m*Omega_H = (c^3*g(chi)/(G*M)) * x_k
                        = (c^3*g(chi)/(G*M)) * ln(k)/(8*pi),
    Delta_f_k = (omega_k - m*Omega_H)/(2*pi)
              = c^3*g(chi)*ln(k) / (16*pi^2*G*M),
    f_{k,m}   = m*Omega_H/(2*pi) + Delta_f_k.

Frame rule. The formulas return source-frame hertz when given source-frame
mass. Observed detector-frame hertz require the redshifted mass
M_det=(1+z)M_source in both Omega_H and the tooth offset. Equivalently, every
source-frame frequency is divided by 1+z. The offset-subtracted ratio is
redshift invariant only under this consistent transformation.

Constants. Only exact definitional constants enter: the SI defined
c = 299792458 m/s and the IAU 2015 Resolution B3 nominal solar mass
parameter (GM)_sun = 1.3271244e20 m^3/s^2. Masses are parameterized in
nominal solar masses so that G*M = mass_solar * (GM)_sun with no separate
measured G; this parameterization is a declared selection recorded in the
receipt.

Declared selections (all recorded in the receipt): the imported
integer-division rule; the linewidth nuisance
range a in [1, 10] and the display endpoints a in {1, 10}; the Page
emission coefficient p_0 = 2e-4 as pinned in the companion statement; the
tooth range k in {2, ..., 12} matching the historical draft's proposed
finite ladder; the mass parameterization above; and the synthetic reference
point M = 62 nominal solar masses, chi = 0.67, m = 2, which is a
synthetic reference, not an event fit, and matches no published remnant
posterior.

Numeric discipline. Receipt values request 50 significant Decimal digits;
public helpers request the caller's precision. Outward interval arithmetic
refines complete expressions, with at most 4096 additional working digits.
Logarithms and square roots use neighboring values around the correctly
rounded Decimal primitives; Machin pi uses exact alternating-series bounds.
Results use HALF_EVEN rounding and are checked against their enclosures for
error below one output ulp and relative error at most 10^(1-p), where p is
the requested precision. This is not a correct-rounding promise. Unresolved
cancellation or output range is refused. Inputs named pi are exact supplied
finite values, not replacements by mathematical pi. Caller Decimal context
settings and flags are preserved; receipt construction owns a fixed context.
Rendered strings carry 40 significant digits; float renderings live under
`derived_for_display`.
The canonical serialization is sorted-key, separator-minimal JSON with a
trailing newline, no timestamps, and no machine paths, so the receipt is
checksum-stable. The independent verifier
`verify_integer_k_comb_independent.py` re-derives every number through
different numeric primitives and byte-compares the canonical JSON.

What is not proved here. This module is a build-stage template
instrument, target-blind by construction: no gravitational-wave event
data, no remnant posterior, no detector likelihood, and no comparison
dataset is read, fetched, or evaluated anywhere in this directory. The
integer-division selection and physical reading of the template scale (that
G*M/(c^3*g(chi)) belongs to a Kerr remnant) are imported continuation inputs,
not derived here. The KMS factor is not a cross-k transition prior, and the
small-transition first-law formula is not an exact finite-step result.
Posterior samples alone do not supply a detector likelihood or evidence.
Nothing in this module or its receipt is a registered,
frozen, or scored prediction; the registration contract in this
directory is a draft pending the owner's freeze. The strain likelihood/evidence
interface, prior normalization, source-transition derivation, event selection,
and trials accounting demanded by the scientific-status erratum are open.
"""

from __future__ import annotations

import hashlib
import json
import os
from decimal import (
    Context, Decimal, DecimalException, DivisionByZero, InvalidOperation,
    MAX_EMAX, MAX_PREC, MIN_EMIN, Overflow, ROUND_CEILING, ROUND_FLOOR,
    ROUND_HALF_EVEN, getcontext, localcontext,
)
from fractions import Fraction

# Requested receipt precision; interval evaluation adds guard digits.
# Rendered strings carry SIG_DIGITS significant digits.
WORKING_PRECISION = 50
SIG_DIGITS = 40

# Exact definitional constants (strings; parsed into Decimal).
C_LIGHT_M_PER_S = "299792458"  # SI defined value, exact.
GM_SUN_NOMINAL_M3_PER_S2 = "1.3271244E20"  # IAU 2015 B3 nominal, exact.

# Declared selections.
DECLARED_A_RANGE = ("1", "10")  # linewidth nuisance a in [1, 10], declared.
DECLARED_A_DISPLAY = ("1", "10")  # display endpoints, declared.
DECLARED_P0 = "2E-4"  # Page emission coefficient, statement-pinned, declared.
DECLARED_K_MIN = 2
DECLARED_K_MAX = 12  # KILL-condition ladder set {2, ..., 12}, declared.
DECLARED_REFERENCE_MASS_SOLAR = "62"  # synthetic reference, declared.
DECLARED_REFERENCE_CHI = "0.67"  # synthetic reference, declared.
DECLARED_REFERENCE_M_AZIMUTHAL = 2  # synthetic reference, declared.

RECEIPT_BASENAME = "integer_k_comb_template_receipt.json"


MAX_GUARD_DIGITS = 4096


class NumericalResolutionError(ArithmeticError):
    """The requested precision/range could not enclose a reportable value."""


def _context(precision: int, rounding=ROUND_HALF_EVEN, *, emin=MIN_EMIN, emax=MAX_EMAX) -> Context:
    # Do not inherit ambient rounding, flags, exponent bounds or trap settings.
    return Context(
        prec=precision, rounding=rounding, Emin=emin, Emax=emax,
        capitals=1, clamp=0, flags=[],
        traps=[InvalidOperation, DivisionByZero, Overflow],
    )


def _power10(exponent: int) -> Decimal:
    return Decimal((0, (1,), exponent))


def _scalar(value: Decimal | int, name: str, *, positive=False, nonnegative=False) -> Decimal:
    if type(value) is int:
        value = Decimal(value)
    if not isinstance(value, Decimal):
        raise TypeError(f"{name} must be an exact Decimal or integer")
    if not value.is_finite():
        raise ValueError(f"{name} must be finite")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _integer(value: int, name: str, minimum: int | None = None) -> int:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer, not a coerced numeric value")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _chi(value: Decimal) -> Decimal:
    value = _scalar(value, "chi")
    # copy_abs is exact; ambient Decimal abs/unary-plus can round the source.
    if value.copy_abs() > 1:
        raise ValueError("chi must lie in [-1, 1]")
    return value


class _Interval:
    """Closed Decimal enclosure; every arithmetic operation rounds outward."""

    def __init__(self, arithmetic, lo, hi=None):
        self.arithmetic = arithmetic
        self.lo = lo
        self.hi = lo if hi is None else hi

    def __add__(self, other):
        other = self.arithmetic(other)
        return _Interval(self.arithmetic,
                         self.arithmetic.down.add(self.lo, other.lo),
                         self.arithmetic.up.add(self.hi, other.hi))

    __radd__ = __add__

    def __neg__(self):
        return _Interval(self.arithmetic, self.hi.copy_negate(), self.lo.copy_negate())

    def __sub__(self, other):
        return self + -self.arithmetic(other)

    def __rsub__(self, other):
        return self.arithmetic(other) + -self

    def __mul__(self, other):
        other = self.arithmetic(other)
        pairs = ((a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi))
        pairs = tuple(pairs)
        return _Interval(self.arithmetic,
                         min(self.arithmetic.down.multiply(a, b) for a, b in pairs),
                         max(self.arithmetic.up.multiply(a, b) for a, b in pairs))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.arithmetic(other)
        if other.lo <= 0 <= other.hi:
            raise _Refine()
        pairs = tuple((a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi))
        return _Interval(self.arithmetic,
                         min(self.arithmetic.down.divide(a, b) for a, b in pairs),
                         max(self.arithmetic.up.divide(a, b) for a, b in pairs))

    def __rtruediv__(self, other):
        return self.arithmetic(other) / self


class _Refine(Exception):
    pass


class _Arithmetic:
    def __init__(self, precision):
        self.precision = precision
        self.down = _context(precision, ROUND_FLOOR)
        self.up = _context(precision, ROUND_CEILING)
        self.nearest = _context(precision)

    def __call__(self, value):
        return value if isinstance(value, _Interval) else _Interval(self, Decimal(value))

    def sqrt(self, value):
        value = self(value)
        if value.lo < 0:
            raise _Refine()
        if value.lo == value.hi == 0:
            return self(0)
        # Decimal sqrt/ln are correctly rounded HALF_EVEN. Adjacent values
        # enclose the exact result regardless of the endpoint's rounding side.
        low = self.nearest.sqrt(value.lo)
        high = self.nearest.sqrt(value.hi)
        return _Interval(self, self.nearest.next_minus(low) if low else low,
                         self.nearest.next_plus(high))

    def ln(self, value):
        value = self(value)
        if value.lo <= 0:
            raise _Refine()
        if value.lo == value.hi == 1:
            return self(0)
        return _Interval(self, self.nearest.next_minus(self.nearest.ln(value.lo)),
                         self.nearest.next_plus(self.nearest.ln(value.hi)))

    def arctan_inverse(self, denominator):
        # Exact rational partial sums with the alternating-series remainder.
        # This also encloses pi; no agreement-between-precisions assumption.
        power = Fraction(1, denominator)
        total = Fraction(0)
        limit = Fraction(1, 10 ** (self.precision + 4))
        index, sign = 0, 1
        while True:
            total += sign * power / (2 * index + 1)
            power /= denominator * denominator
            index += 1
            tail = power / (2 * index + 1)
            sign = -sign
            if tail < limit:
                low, high = sorted((total, total + sign * tail))
                return _Interval(
                    self,
                    self.down.divide(Decimal(low.numerator), Decimal(low.denominator)),
                    self.up.divide(Decimal(high.numerator), Decimal(high.denominator)),
                )


def _evaluate(expression) -> Decimal:
    """Return at caller precision, HALF_EVEN, with a checked error < 1 ulp.

    A second relative bound, 10**(1-p), prevents coarse subnormal rounding
    from silently discarding the requested significant precision. The caller's
    exponent range is respected. At most 4096 extra working digits are tried;
    unresolved cancellation/range is refused, not replaced by a numerical zero.
    """
    caller = getcontext()
    precision = caller.prec
    output = _context(precision, emin=caller.Emin, emax=caller.Emax)
    ceiling = min(MAX_PREC, precision + MAX_GUARD_DIGITS)
    work = min(ceiling, precision + 16)
    relative_width = _power10(-precision - 3)
    relative_error = _power10(1 - precision)
    while True:
        arithmetic = _Arithmetic(work)
        try:
            interval = expression(arithmetic)
            lo, hi = interval.lo, interval.hi
            if lo == hi == 0:
                return Decimal(0)
            if lo <= 0 <= hi:
                raise _Refine()
            scale = min(lo.copy_abs(), hi.copy_abs())
            width = arithmetic.up.subtract(hi, lo)
            if width > arithmetic.down.multiply(scale, relative_width):
                raise _Refine()
            midpoint = arithmetic.nearest.add(lo, arithmetic.nearest.divide(width, Decimal(2)))
            result = output.plus(midpoint)
            if not result.is_finite() or result == 0:
                raise NumericalResolutionError("nonzero result is outside the requested Decimal output range")
            error = max(
                arithmetic.up.subtract(max(result, endpoint), min(result, endpoint))
                for endpoint in (lo, hi)
            )
            ulp = _power10(max(result.adjusted(), output.Emin) - precision + 1)
            if error < ulp and error <= arithmetic.down.multiply(scale, relative_error):
                return result
        except _Refine:
            pass
        except DecimalException as exc:
            raise NumericalResolutionError("Decimal exponent range exhausted during evaluation") from exc
        if work == ceiling:
            raise NumericalResolutionError("requested accuracy unresolved within the guard-digit budget")
        work = min(ceiling, 2 * work)


def _spin(arithmetic, chi):
    spin = arithmetic(chi)
    root = arithmetic.sqrt((1 - spin) * (1 + spin))
    return root, 1 + root


def _kerr(arithmetic, chi):
    root, horizon = _spin(arithmetic, chi)
    c = int(C_LIGHT_M_PER_S)
    # Keep the mass outside this scale: G*M or the separate dimensional
    # frequency terms can overflow while the requested result is finite.
    scale = arithmetic(c ** 3) / (2 * arithmetic(dec(GM_SUN_NOMINAL_M3_PER_S2)) * horizon)
    return root, scale


def _frequency(arithmetic, mass, chi, m, k, pi):
    root, scale = _kerr(arithmetic, chi)
    pi = arithmetic(pi)
    return scale * (m * arithmetic(chi) / (2 * pi)
                    + root * arithmetic.ln(k) / (4 * pi * pi)) / mass


def _arctan_inv(x: int) -> Decimal:
    _integer(x, "arctan denominator", 2)
    return _evaluate(lambda arithmetic: arithmetic.arctan_inverse(x))


def compute_pi() -> Decimal:
    """Machin pi, enclosed by exact alternating-series remainder bounds."""
    return _evaluate(lambda arithmetic: 16 * arithmetic.arctan_inverse(5) - 4 * arithmetic.arctan_inverse(239))


def dec(value: str | int | Decimal) -> Decimal:
    """Parse an exact finite decimal; never import a binary float silently."""
    if isinstance(value, str):
        # Parsing invalid text must not change the caller's Decimal flags.
        try:
            with localcontext(_context(WORKING_PRECISION)):
                value = Decimal(value)
        except DecimalException as exc:
            raise ValueError("value must be an exact finite decimal") from exc
    return _scalar(value, "value")


def integer_division_after(d_before: int, k: int) -> int:
    """Return d_after for the imported rule d_before = k*d_after.

    The function rejects nonpositive dimensions, k < 2, and nondivisible
    pairs. It validates the continuation premise; it does not derive or select
    a physical transition.
    """
    _integer(d_before, "d_before")
    _integer(k, "k", 2)
    if d_before <= 0:
        raise ValueError("d_before must be positive")
    d_after, remainder = divmod(d_before, k)
    if remainder != 0:
        raise ValueError("k must divide d_before")
    if d_after <= 0:
        raise ValueError("d_after must be positive")
    return d_after


def transition_entropy_nats(d_before: int, k: int) -> tuple[Decimal, Decimal]:
    """Signed black-hole entropy change and positive entropy loss in nats."""
    integer_division_after(d_before, k)
    entropy_loss = _evaluate(lambda arithmetic: arithmetic.ln(k))
    return entropy_loss.copy_negate(), entropy_loss


def detector_frame_mass_solar(
    source_frame_mass_solar: Decimal, redshift: Decimal
) -> Decimal:
    """M_det=(1+z)M_source for an observed-frequency template."""
    mass = _scalar(source_frame_mass_solar, "source-frame mass", positive=True)
    redshift = _scalar(redshift, "redshift", nonnegative=True)
    return _evaluate(lambda arithmetic: arithmetic(mass) * (1 + arithmetic(redshift)))


def sqrt_one_minus_chi_squared(chi: Decimal) -> Decimal:
    """s(chi) = sqrt(1 - chi^2), the Kerr root factor."""
    chi = _chi(chi)
    return _evaluate(lambda arithmetic: _spin(arithmetic, chi)[0])


def r_plus_hat(chi: Decimal) -> Decimal:
    """Outer horizon radius in units of G*M/c^2: 1 + sqrt(1 - chi^2)."""
    chi = _chi(chi)
    return _evaluate(lambda arithmetic: _spin(arithmetic, chi)[1])


def gm_si(mass_solar: Decimal) -> Decimal:
    """G*M in m^3/s^2 from the nominal solar mass parameter."""
    mass = _scalar(mass_solar, "mass", positive=True)
    return _evaluate(lambda arithmetic: arithmetic(mass) * dec(GM_SUN_NOMINAL_M3_PER_S2))


def r_plus_si(mass_solar: Decimal, chi: Decimal) -> Decimal:
    """Outer horizon radius in meters: (G*M/c^2)*(1 + sqrt(1 - chi^2))."""
    mass, chi = _scalar(mass_solar, "mass", positive=True), _chi(chi)
    return _evaluate(lambda arithmetic: (
        arithmetic(mass) * dec(GM_SUN_NOMINAL_M3_PER_S2)
        * _spin(arithmetic, chi)[1] / int(C_LIGHT_M_PER_S) ** 2
    ))


def omega_h_si(mass_solar: Decimal, chi: Decimal) -> Decimal:
    """Horizon angular frequency in rad/s:
    c^3*chi / (2*G*M*(1 + sqrt(1 - chi^2)))."""
    mass, chi = _scalar(mass_solar, "mass", positive=True), _chi(chi)
    if chi == 0:
        return Decimal(0)
    return _evaluate(lambda arithmetic: _kerr(arithmetic, chi)[1] * chi / mass)


def kappa_si(mass_solar: Decimal, chi: Decimal) -> Decimal:
    """Surface gravity in 1/s:
    c^3*sqrt(1 - chi^2) / (2*G*M*(1 + sqrt(1 - chi^2)))."""
    mass, chi = _scalar(mass_solar, "mass", positive=True), _chi(chi)
    if chi.copy_abs() == 1:
        return Decimal(0)
    def expression(arithmetic):
        root, scale = _kerr(arithmetic, chi)
        return scale * root / mass
    return _evaluate(expression)


def g_of_chi(chi: Decimal) -> Decimal:
    """Statement-pinned spin factor
    g(chi) = 2*sqrt(1 - chi^2)/(1 + sqrt(1 - chi^2)) = 4*G*M*kappa/c^3."""
    chi = _chi(chi)
    def expression(arithmetic):
        root, horizon = _spin(arithmetic, chi)
        return 2 * root / horizon
    return _evaluate(expression)


def base_spacing_hz_per_nat(mass_solar: Decimal, chi: Decimal, pi: Decimal) -> Decimal:
    """Tooth spacing per nat of ln(k): c^3*g(chi) / (16*pi^2*G*M) in Hz."""
    mass, chi = _scalar(mass_solar, "mass", positive=True), _chi(chi)
    pi = _scalar(pi, "pi", positive=True)
    if chi.copy_abs() == 1:
        return Decimal(0)
    def expression(arithmetic):
        root, scale = _kerr(arithmetic, chi)
        return scale * root / (4 * arithmetic(pi) * pi) / mass
    return _evaluate(expression)


def rotation_line_hz(mass_solar: Decimal, chi: Decimal, m: int, pi: Decimal) -> Decimal:
    """Rotation line m*Omega_H/(2*pi) in Hz."""
    mass, chi = _scalar(mass_solar, "mass", positive=True), _chi(chi)
    m, pi = _integer(m, "m"), _scalar(pi, "pi", positive=True)
    if m == 0 or chi == 0:
        return Decimal(0)
    return _evaluate(lambda arithmetic: (
        m * _kerr(arithmetic, chi)[1] * chi / (2 * arithmetic(pi)) / mass
    ))


def universal_position(k: int, pi: Decimal) -> Decimal:
    """Frozen universal-coordinate tooth position x_k = ln(k)/(8*pi)."""
    k, pi = _integer(k, "k", 2), _scalar(pi, "pi", positive=True)
    return _evaluate(lambda arithmetic: arithmetic.ln(k) / (8 * arithmetic(pi)))


def ladder_ratio(k: int) -> Decimal:
    """Offset-subtracted ratio against the k = 2 tooth: ln(k)/ln(2)."""
    k = _integer(k, "k", 2)
    if k == 2:
        return Decimal(1)
    return _evaluate(lambda arithmetic: arithmetic.ln(k) / arithmetic.ln(2))


def kms_weight(k: int) -> Decimal:
    """Within-line KMS net-response factor (k-1)/k.

    The legacy function name does not make this a normalized transition
    probability or prior across different k.
    """
    k = _integer(k, "k", 2)
    return _evaluate(lambda arithmetic: arithmetic(k - 1) / k)


def tooth_offset_hz(mass_solar: Decimal, chi: Decimal, k: int, pi: Decimal) -> Decimal:
    """Delta_f_k = c^3*g(chi)*ln(k) / (16*pi^2*G*M) in Hz."""
    mass, chi = _scalar(mass_solar, "mass", positive=True), _chi(chi)
    k, pi = _integer(k, "k", 2), _scalar(pi, "pi", positive=True)
    if chi.copy_abs() == 1:
        return Decimal(0)
    def expression(arithmetic):
        root, scale = _kerr(arithmetic, chi)
        return scale * root * arithmetic.ln(k) / (4 * arithmetic(pi) * pi) / mass
    return _evaluate(expression)


def tooth_frequency_hz(
    mass_solar: Decimal, chi: Decimal, m: int, k: int, pi: Decimal
) -> Decimal:
    """Frequency in the frame of ``mass_solar``.

    Pass source-frame mass for source-frame hertz and detector-frame mass for
    observed detector-frame hertz.
    """
    mass, chi = _scalar(mass_solar, "mass", positive=True), _chi(chi)
    m, k = _integer(m, "m"), _integer(k, "k", 2)
    pi = _scalar(pi, "pi", positive=True)
    if m == 0 and chi.copy_abs() == 1:
        return Decimal(0)
    return _evaluate(lambda arithmetic: _frequency(arithmetic, mass, chi, m, k, pi))


def detector_frame_tooth_frequency_hz(
    source_frame_mass_solar: Decimal,
    redshift: Decimal,
    chi: Decimal,
    m: int,
    k: int,
    pi: Decimal,
) -> Decimal:
    """Observed tooth frequency using M_det=(1+z)M_source."""
    mass = _scalar(source_frame_mass_solar, "source-frame mass", positive=True)
    redshift, chi = _scalar(redshift, "redshift", nonnegative=True), _chi(chi)
    m, k = _integer(m, "m"), _integer(k, "k", 2)
    pi = _scalar(pi, "pi", positive=True)
    if m == 0 and chi.copy_abs() == 1:
        return Decimal(0)
    return _evaluate(lambda arithmetic: _frequency(
        arithmetic, arithmetic(mass) * (1 + arithmetic(redshift)), chi, m, k, pi
    ))


def linewidth_fraction(a: Decimal, chi: Decimal, k: int, pi: Decimal) -> Decimal:
    """Linewidth-to-spacing ratio 64*pi^2*p_0/(a*g(chi)^2*ln(k)).

    The mass cancels, but the Kerr spin factor does not. The constant-p_0
    approximation is a declared template nuisance model, not a controlled
    near-extremal Page calculation."""
    a, chi = _scalar(a, "a", positive=True), _chi(chi)
    k, pi = _integer(k, "k", 2), _scalar(pi, "pi", positive=True)
    if chi.copy_abs() == 1:
        raise ValueError("linewidth is singular at extremal chi")
    def expression(arithmetic):
        root, horizon = _spin(arithmetic, chi)
        g_chi = 2 * root / horizon
        return (64 * arithmetic(pi) * pi * dec(DECLARED_P0)
                / (a * g_chi * g_chi * arithmetic.ln(k)))
    return _evaluate(expression)


def sig40(x: Decimal) -> str:
    """Render exactly SIG_DIGITS significant digits in scientific form."""
    x = _scalar(x, "value")
    with localcontext(_context(WORKING_PRECISION)):
        return format(x, ".%dE" % (SIG_DIGITS - 1))


def build_receipt() -> dict:
    """Assemble the full receipt dictionary (pure; no I/O)."""
    with localcontext(_context(WORKING_PRECISION, emin=-999999, emax=999999)):
        return _build_receipt()


def _build_receipt() -> dict:
    pi = compute_pi()

    mass = dec(DECLARED_REFERENCE_MASS_SOLAR)
    chi = dec(DECLARED_REFERENCE_CHI)
    m_az = DECLARED_REFERENCE_M_AZIMUTHAL
    ks = list(range(DECLARED_K_MIN, DECLARED_K_MAX + 1))

    display: dict[str, float] = {}

    ladder = []
    for k in ks:
        xk = universal_position(k, pi)
        rk = ladder_ratio(k)
        wk = kms_weight(k)
        key = "k%02d" % k
        ladder.append(
            {
                "k": k,
                "x_exact": "ln(%d)/(8*pi)" % k,
                "x_sig40": sig40(xk),
                "ratio_to_k2_exact": "ln(%d)/ln(2)" % k,
                "ratio_to_k2_sig40": sig40(rk),
                "kms_weight_exact": "%d/%d" % (k - 1, k),
                "kms_weight_sig40": sig40(wk),
            }
        )
        display["universal_ladder.%s.x" % key] = float(xk)
        display["universal_ladder.%s.ratio_to_k2" % key] = float(rk)
        display["universal_ladder.%s.kms_weight" % key] = float(wk)

    omega_h = omega_h_si(mass, chi)
    kappa = kappa_si(mass, chi)
    g_chi = g_of_chi(chi)
    rhat = r_plus_hat(chi)
    r_plus_m = r_plus_si(mass, chi)
    rot = rotation_line_hz(mass, chi, m_az, pi)
    base = base_spacing_hz_per_nat(mass, chi, pi)

    display["reference.omega_h_rad_per_s"] = float(omega_h)
    display["reference.kappa_per_s"] = float(kappa)
    display["reference.g_chi"] = float(g_chi)
    display["reference.r_plus_hat"] = float(rhat)
    display["reference.r_plus_m"] = float(r_plus_m)
    display["reference.rotation_line_hz"] = float(rot)
    display["reference.base_spacing_hz_per_nat"] = float(base)

    teeth = []
    a_lo = dec(DECLARED_A_DISPLAY[0])
    a_hi = dec(DECLARED_A_DISPLAY[1])
    for k in ks:
        dfk = tooth_offset_hz(mass, chi, k, pi)
        fk = tooth_frequency_hz(mass, chi, m_az, k, pi)
        lw_lo = linewidth_fraction(a_hi, chi, k, pi)  # a = 10: narrow end
        lw_hi = linewidth_fraction(a_lo, chi, k, pi)  # a = 1: wide end
        key = "k%02d" % k
        teeth.append(
            {
                "k": k,
                "delta_f_hz_sig40": sig40(dfk),
                "f_hz_sig40": sig40(fk),
                "linewidth_fraction_a10_sig40": sig40(lw_lo),
                "linewidth_fraction_a1_sig40": sig40(lw_hi),
            }
        )
        display["reference.teeth.%s.delta_f_hz" % key] = float(dfk)
        display["reference.teeth.%s.f_hz" % key] = float(fk)
        display["reference.teeth.%s.linewidth_fraction_a10" % key] = float(lw_lo)
        display["reference.teeth.%s.linewidth_fraction_a1" % key] = float(lw_hi)

    receipt = {
        "schema": "oph.ringdown.integer_k_comb_template.v3",
        "status": (
            "build-stage instrument, target-blind draft; not a registered, "
            "frozen, or scored prediction"
        ),
        "imported_continuation_law": {
            "source_of_record": (
                "falsification/frozen_targets/fz01_2026-07-17/"
                "frozen_target_integer_k_comb_2026-07-17.md"
            ),
            "companion_statement": (
                "proof/epic_wins/ringdown_comb/INTEGER_K_COMB_STATEMENT.md"
            ),
            "erratum": (
                "falsification/frozen_targets/fz01_2026-07-17/"
                "SCIENTIFIC_STATUS_ERRATUM_2026-07-29.md"
            ),
            "ratio_law": (
                "(f_a - m*Omega_H/(2*pi)) / (f_b - m*Omega_H/(2*pi)) "
                "= ln(k_a)/ln(k_b), integers k >= 2"
            ),
            "transition_rule": (
                "d_before = k*d_after with positive integer dimensions and "
                "k >= 2; therefore k must divide d_before"
            ),
            "transition_status": (
                "imported continuation premise; not source-derived by OPH"
            ),
            "signed_black_hole_entropy_change": (
                "ln(d_after)-ln(d_before) = -ln(k)"
            ),
            "positive_entropy_loss": (
                "ln(d_before)-ln(d_after) = ln(k)"
            ),
            "universal_coordinate": (
                "x = (G*M/(c^3*g(chi))) * (omega - m*Omega_H); "
                "x_k = ln(k)/(8*pi)"
            ),
            "g_of_chi": (
                "g(chi) = 2*sqrt(1-chi^2)/(1+sqrt(1-chi^2)) "
                "= 4*G*M*kappa(M,chi)/c^3; statement-pinned via the KMS "
                "reading omega - m*Omega_H = (kappa/(2*pi))*ln(k)"
            ),
            "tooth_offset": "Delta_f_k = c^3*g(chi)*ln(k)/(16*pi^2*G*M)",
            "kms_weight": "(k-1)/k",
            "kms_scope": (
                "within-line net absorption-minus-stimulated-emission factor; "
                "not a normalized cross-k transition probability or prior"
            ),
            "linewidth_fraction": "64*pi^2*p_0/(a*g(chi)^2*ln(k))",
            "linewidth_scope": "Gamma/Delta_E_k linewidth-to-spacing ratio",
            "transition_approximation": (
                "leading small-transition Kerr first-law template; finite-step "
                "background corrections are not derived"
            ),
            "computed_k_scope": (
                "finite display/submodel candidate k in {2,...,12}; not an "
                "unrestricted integer-family likelihood"
            ),
            "kerr_functions": {
                "r_plus": "(G*M/c^2)*(1 + sqrt(1-chi^2))",
                "omega_h": "c^3*chi/(2*G*M*(1 + sqrt(1-chi^2)))",
                "kappa": "c^3*sqrt(1-chi^2)/(2*G*M*(1 + sqrt(1-chi^2)))",
            },
        },
        "constants_exact": {
            "c_m_per_s": C_LIGHT_M_PER_S,
            "c_provenance": "SI defined value, exact",
            "gm_sun_nominal_m3_per_s2": GM_SUN_NOMINAL_M3_PER_S2,
            "gm_sun_provenance": (
                "IAU 2015 Resolution B3 nominal solar mass parameter, "
                "exact nominal value"
            ),
        },
        "declared_selections": {
            "a_range": {
                "value": list(DECLARED_A_RANGE),
                "status": "declared; historical-draft nuisance interval",
            },
            "a_display_endpoints": {
                "value": list(DECLARED_A_DISPLAY),
                "status": "declared; display endpoints of the a range",
            },
            "p_0": {
                "value": DECLARED_P0,
                "status": (
                    "declared; Page emission coefficient as pinned in the "
                    "historical companion statement"
                ),
            },
            "k_range": {
                "value": [DECLARED_K_MIN, DECLARED_K_MAX],
                "status": (
                    "declared finite display/submodel candidate matching the "
                    "historical draft's ladder set {2, ..., 12}; not an "
                    "unrestricted integer-family support"
                ),
            },
            "mass_parameterization": {
                "value": "G*M = mass_solar * (GM)_sun nominal",
                "status": (
                    "declared; nominal solar mass parameter in place of a "
                    "separate measured G"
                ),
            },
            "g_of_chi_provenance": {
                "value": "algebraically fixed within the imported KMS reading",
                "status": (
                    "not a free selection; the identity "
                    "g = 4*G*M*kappa/c^3 is recorded in "
                    "imported_continuation_law.g_of_chi"
                ),
            },
            "reference_point": {
                "value": {
                    "mass_nominal_solar": DECLARED_REFERENCE_MASS_SOLAR,
                    "chi": DECLARED_REFERENCE_CHI,
                    "m_azimuthal": DECLARED_REFERENCE_M_AZIMUTHAL,
                },
                "status": (
                    "declared; synthetic-reference, not an event fit; "
                    "matches no published remnant posterior"
                ),
            },
        },
        "universal_ladder": ladder,
        "reference_point_synthetic": {
            "label": "synthetic-reference, not an event fit",
            "mass_nominal_solar": DECLARED_REFERENCE_MASS_SOLAR,
            "chi": DECLARED_REFERENCE_CHI,
            "m_azimuthal": DECLARED_REFERENCE_M_AZIMUTHAL,
            "frame": "source frame; observed hertz require M_det=(1+z)M_source",
            "omega_h_rad_per_s_sig40": sig40(omega_h),
            "kappa_per_s_sig40": sig40(kappa),
            "g_chi_sig40": sig40(g_chi),
            "r_plus_hat_sig40": sig40(rhat),
            "r_plus_m_sig40": sig40(r_plus_m),
            "rotation_line_hz_sig40": sig40(rot),
            "base_spacing_hz_per_nat_sig40": sig40(base),
            "teeth": teeth,
        },
        "derived_for_display": display,
        "numerics": {
            "working_precision_decimal_digits": WORKING_PRECISION,
            "rendered_significant_digits": SIG_DIGITS,
            "pi_method": "Machin: 16*arctan(1/5) - 4*arctan(1/239)",
        },
        "frame_contract": {
            "detector_mass": "M_det = (1+z)*M_source",
            "observed_frequency_rule": (
                "evaluate Omega_H and every tooth with M_det; equivalently "
                "divide all source-frame frequencies by 1+z"
            ),
            "ratio_boundary": (
                "redshift cancellation holds only when tooth and rotation "
                "offset use the same frame"
            ),
        },
        "comparison_contract_boundary": (
            "Published posterior samples alone are not a detector likelihood "
            "or model evidence. A future comparison needs a common "
            "strain/readout likelihood or sufficient likelihood product, "
            "sampling-prior and normalization metadata, prospectively "
            "normalized model priors, and convergence diagnostics."
        ),
        "boundary": (
            "Target-blind build-stage instrument. No gravitational-wave "
            "event data, remnant posterior, detector likelihood, or "
            "comparison dataset was read, fetched, or evaluated. The "
            "integer-division transition and physical reading of the template "
            "scale as a Kerr remnant quantity are imported continuation "
            "premises. Source derivation and registration remain open."
        ),
    }
    return receipt


def canonical_bytes(receipt: dict) -> bytes:
    """Canonical serialization: sorted keys, minimal separators, ASCII,
    trailing newline. No timestamps and no machine paths appear."""
    return (
        json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode("ascii")
        + b"\n"
    )


def receipt_path() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "runtime", RECEIPT_BASENAME)


def main() -> int:
    receipt = build_receipt()
    payload = canonical_bytes(receipt)
    path = receipt_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(payload)
    digest = hashlib.sha256(payload).hexdigest()
    print("wrote %s" % os.path.relpath(path, os.path.dirname(os.path.abspath(__file__))))
    print("sha256 %s" % digest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
