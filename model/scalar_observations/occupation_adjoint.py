"""Observation-specific residual bounds using the complete reverse instrument.

The complete-output occupation certificate is often much larger than the error
of a classical event. Exterior forward residuals are orthogonal to an in-code
reverse trajectory, so their contribution to an event probability has a second
residual factor. Both reads are reversed as full pointer unitaries; their
outcomes are never reset, normalized, or treated as coherent Weyl alternatives.
"""
from fractions import Fraction as F
from math import lcm

from .bounds import sqrt_upper
from .occupation_exact import (Generator, Vector, _add, combine, flow,
                               norm_squared)


def _reverse_read_combine(positive, negative, bit):
    """Apply the raw inverse pointer rotation, keeping both existing bits.

    R = [[K0,-K1],[K1,K0]], K0=(C-S)/sqrt(2), K1=(C+S)/sqrt(2).
    This returns sqrt(2)*R* on the input; its caller halves the squared
    state weight. The two arguments are exp(+/-i Q/2) applied to that input.
    """
    denominator = 2*lcm(positive.denominator, negative.denominator)
    out = {}
    for vector, phase in ((positive, 1), (negative, -1)):
        scale = denominator//(2*vector.denominator)
        for key, (real, imaginary) in vector.coefficients.items():
            old = key[7+bit]
            factors = ((1, phase), (-1, phase)) if old == 0 else ((1, -phase), (1, phase))
            for new, (a, b) in enumerate(factors):
                destination = key[:7+bit]+(new,)+key[8+bit:]
                _add(out, destination,
                     scale*(a*real-b*imaginary), scale*(b*real+a*imaginary))
    return Vector(out, denominator)


class _PointerCoordinate:
    """Exact full generator Q tensor Y on one existing pointer register."""

    def __init__(self, dimensions, bit):
        self.coordinate = Generator(tuple(dimensions), "Q")
        self.bit = bit

    def apply(self, vector, exterior=False):
        applied = self.coordinate.apply(vector, exterior=exterior)
        out = {}
        for key, (real, imaginary) in applied.coefficients.items():
            old = key[7+self.bit]
            destination = key[:7+self.bit]+(1-old,)+key[8+self.bit:]
            # Y|0>=i|1>, Y|1>=-i|0>.
            phase = 1 if old == 0 else -1
            _add(out, destination, -phase*imaginary, phase*real)
        return Vector(out, applied.denominator)


def _inverse_pointer_phase(vector, bit):
    """Return sqrt(2)*exp(+i*pi*Y/4) vector using exact integer sums."""
    out = {}
    for key, (real, imaginary) in vector.coefficients.items():
        old = key[7+bit]
        for new in (0, 1):
            phase = -1 if old == 0 and new == 1 else 1
            destination = key[:7+bit]+(new,)+key[8+bit:]
            _add(out, destination, phase*real, phase*imaginary)
    return Vector(out, vector.denominator)


def reverse_history(vector, weight, dimensions, coupling, sign, interval, *,
                    degree=8, substeps=8, bits=100):
    """Bound every reverse prefix from the supplied unnormalized terminal input.

    Both pointer bits exist throughout this reversed unitary dilation. A
    retained event projector is applied before this call, not between stages.
    Sum-of-local-error bounds are monotone and control partial polynomial
    intervals as well as endpoints.
    """
    hamiltonian = Generator(tuple(dimensions), "H", coupling)
    coordinate = Generator(tuple(dimensions), "Q", coupling)
    total = F(0)
    stages = []

    def advance(name, generator, time, selected=None):
        nonlocal vector, total
        current = vector if selected is None else Vector(
            {key: value for key, value in vector.coefficients.items() if selected(key)},
            vector.denominator)
        updated, budget = flow(generator, current, time, degree, substeps, bits)
        if selected is not None:
            kept = Vector({key: value for key, value in vector.coefficients.items()
                           if not selected(key)}, vector.denominator)
            updated = combine([updated, kept])
        vector = updated
        error = sqrt_upper(weight)*budget["unrestricted_vector_error_upper"]
        total += error
        stages.append({"operation": name, "raw_state_squared_weight": weight,
                       "local_vector_error_upper": error,
                       "reverse_prefix_vector_error_upper": total})

    def read(name, bit):
        nonlocal vector, weight, total
        # R* = exp(+i Q Y/2) exp(+i pi Y/4). The constant pointer
        # rotation is exact; its raw norm sqrt(2) is charged to the weight.
        vector = _inverse_pointer_phase(vector, bit)
        weight /= 2
        vector, budget = flow(_PointerCoordinate(dimensions, bit), vector,
                              -F(1, 2), degree, substeps, bits)
        error = sqrt_upper(weight)*budget["unrestricted_vector_error_upper"]
        total += error
        stages.append({"operation": name, "raw_state_squared_weight": weight,
                       "local_vector_error_upper": error,
                       "reverse_prefix_vector_error_upper": total,
                       "existing_pointer_bit": bit,
                       "generator": "Q tensor Y after exact inverse pi/4 pointer rotation",
                       "full_inverse_pointer_unitary": True})

    read("inverse_second_read", 1)
    advance("inverse_second_evolution", hamiltonian, -F(interval))
    advance("inverse_conditional_feedback", coordinate, F(1, 4),
            lambda key: key[7] == 1)
    read("inverse_first_read", 0)
    advance("inverse_first_evolution", hamiltonian, -F(interval))
    advance("inverse_initial_kick", coordinate, F(sign, 4))
    return {"state": vector, "raw_state_squared_weight": weight,
            "all_reverse_prefix_error_upper": total, "stages": stages}


def integrated_forward_budgets(forward):
    """Enclose integrated full/inside defects, not only endpoint errors.

    In a read, the two Weyl evolutions run coherently and simultaneously.
    The integral of their direct-sum residual norm is at most the sum of
    their two integral bounds, which is at most sqrt(2) times the recorded
    Euclidean combination. The recombination rounding is charged too. Other
    stages already retain integrated-residual and actual-jump bounds.
    """
    full, inside = F(0), F(0)
    for stage in forward["stages"]:
        factor = sqrt_upper(F(2)) if stage.get("coherent_weyl_branches_retained") else F(1)
        full += factor*stage["unrestricted_vector_error_upper"]
        inside += factor*stage["finite_vector_error_upper"]
    return full, inside


def event_certificate(forward, event):
    """Enclose an unrestricted event from an exact forward polynomial receipt.

    ``forward`` must come from occupation_exact.certificate with
    return_terminal_state=True. ``event`` is a nonempty set of two-bit
    classical histories, selected without using an experimental dataset.
    The shared forward approximation is never normalized.
    """
    if not isinstance(forward, dict) or not isinstance(forward.get("terminal_state"), Vector):
        raise ValueError("a forward certificate with its exact terminal state is required")
    if not isinstance(event, (list, tuple, set, frozenset)) or not event:
        raise ValueError("a nonempty set of two-bit outcomes is required")
    outcomes = []
    for outcome in event:
        if (not isinstance(outcome, (tuple, list)) or len(outcome) != 2
                or any(type(bit) is not int or bit not in (0, 1) for bit in outcome)):
            raise ValueError("event outcomes must each contain two binary integer records")
        outcomes.append(tuple(outcome))
    if len(set(outcomes)) != len(outcomes):
        raise ValueError("event outcomes must be distinct")
    outcomes = tuple(sorted(outcomes))
    terminal = forward["terminal_state"]
    selected = Vector({key: value for key, value in terminal.coefficients.items()
                       if key[7:] in outcomes}, terminal.denominator)
    weight = forward["raw_state_squared_weight"]
    probability = weight*norm_squared(selected)
    reverse = reverse_history(selected, weight, forward["local_dimensions"],
                              forward["coupling"], forward["preparation_sign"],
                              forward["interval"], degree=forward["degree"],
                              substeps=forward["substeps_per_flow"],
                              bits=forward["rounding_bits"])
    integrated_full, integrated_inside = integrated_forward_budgets(forward)
    backward = reverse["all_reverse_prefix_error_upper"]
    endpoint = forward["unrestricted_vector_error_upper"]
    inside_term = 2*sqrt_upper(probability)*integrated_inside
    exterior_term = 2*integrated_full*backward
    quadratic_term = endpoint*endpoint
    radius = inside_term+exterior_term+quadratic_term
    ordinary_radius = 2*sqrt_upper(probability)*endpoint+endpoint*endpoint
    radius = min(radius, ordinary_radius)
    return {"schema": "oph.scalar_observations.adjoint_event.v1",
            "event": [list(outcome) for outcome in outcomes],
            "local_dimensions": forward["local_dimensions"],
            "coupling": forward["coupling"],
            "preparation_sign": forward["preparation_sign"],
            "interval": forward["interval"],
            "degree": forward["degree"],
            "substeps_per_flow": forward["substeps_per_flow"],
            "rounding_bits": forward["rounding_bits"],
            "raw_polynomial_event_probability": probability,
            "forward_endpoint_vector_error_upper": endpoint,
            "forward_integrated_full_defect_and_jumps_upper": integrated_full,
            "forward_integrated_inside_defect_and_jumps_upper": integrated_inside,
            "all_reverse_prefix_error_upper": backward,
            "inside_defect_probability_error_upper": inside_term,
            "exterior_forward_reverse_probability_error_upper": exterior_term,
            "quadratic_endpoint_probability_error_upper": quadratic_term,
            "ordinary_probability_error_upper": ordinary_radius,
            "event_probability_error_upper": radius,
            "unrestricted_event_probability_interval": [max(F(0), probability-radius),
                                                         min(F(1), probability+radius)],
            "reverse_stages": reverse["stages"],
            "scope": {"exact_gradient": True,
                      "same_supplied_seven_mode_action": True,
                      "every_record_retained": True,
                      "event_selection_is_postselection": False,
                      "rational_residual_and_rounding_arithmetic": True,
                      "source_execution_or_empirical_measurement": False}}


def adjoint_certificate(forward, event="bits_differ"):
    """Named principal event, or an explicit nonempty set of two-bit outcomes."""
    if isinstance(event, str):
        if event != "bits_differ":
            raise ValueError("the only named event is bits_differ")
        event = ((0, 1), (1, 0))
    return event_certificate(forward, event)
