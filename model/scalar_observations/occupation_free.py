"""Certified free comparison after tracing out six unobserved oscillator modes.

At g=0 the seven-mode Hamiltonian is a tensor sum. All preparations, reads and
feedback act on cos_x alone, and the blank is a product state. The spectator
unitaries therefore factor out of every joint record probability, including
conditional continuation. The following certificate concerns that exact
reduced observation law; it does not approximate the seven-mode output state.
"""
from fractions import Fraction as F

from .bounds import sqrt_upper
from .occupation_exact import (Generator, Vector, _add, _positive_integer,
                               _read_combine, _signed_time, combine, flow, quantize)
from .occupation_readout import terminal_probabilities


class ObservedFreeHamiltonian:
    """Centered observed-mode H=3N/2-(a**2+a†**2)/4, with Omega=2.

    The physical frequency is sqrt(2): H=P**2/2+Q**2. Its blank energy 3/4
    is a common scalar. All spectator coordinates here are placeholders fixed
    to zero after the tensor-factorization proof, not simulated field states.
    """

    def __init__(self, dimension):
        if type(dimension) is not int or dimension < 2:
            raise ValueError("observed occupation dimension must be an integer of at least two")
        self.dimension = dimension

    def apply(self, vector, exterior=False):
        out = {}
        for key, (real, imaginary) in vector.coefficients.items():
            if any(n != 0 for i, n in enumerate(key[:7]) if i != 1):
                raise ValueError("the reduced free state must have zero spectator placeholders")
            n = key[1]
            if not 0 <= n < self.dimension:
                raise ValueError("the input must lie inside the observed occupation code")
            for destination, coefficient in ((n, 6*n), (n+2, -1), (n-2, -n*(n-1))):
                if coefficient and destination >= 0 and (destination < self.dimension) != exterior:
                    _add(out, key[:1]+(destination,)+key[2:], coefficient*real, coefficient*imaginary)
        return Vector(out, 4*vector.denominator)


def free_prefinal_certificate(dimension=12, *, degree=6, substeps=4, bits=100,
                              sign=1, interval=F(1, 8), return_state=False):
    """Principal free record law, with exact unrestricted terminal effects.

    Every residual and rounding budget is rational. The occupation cutoff
    affects only the observed oscillator. The other six free oscillator
    factors have been traced out exactly, without postselection.
    """
    hamiltonian = ObservedFreeHamiltonian(dimension)
    degree = _positive_integer(degree, "Taylor degree")
    substeps = _positive_integer(substeps, "flow substeps")
    bits = _positive_integer(bits, "rounding bits")
    interval = _signed_time(interval)
    if interval < 0:
        raise ValueError("the observation interval must be nonnegative")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("the native preparation sign must be minus or plus one")
    if type(return_state) is not bool:
        raise ValueError("return_state must be Boolean")
    dimensions = (2, dimension, 2, 2, 2, 2, 2)
    coordinate = Generator(dimensions, "Q", 0)
    vector, weight = Vector({(0,)*7: (1, 0)}), F(1)
    errors = {"unrestricted_vector_error_upper": F(0), "finite_vector_error_upper": F(0)}
    stages = []

    def advance(name, generator, time, selected=None):
        nonlocal vector
        current = vector if selected is None else Vector(
            {key: value for key, value in vector.coefficients.items() if selected(key)},
            vector.denominator)
        updated, budget = flow(generator, current, time, degree, substeps, bits)
        if selected is not None:
            kept = Vector({key: value for key, value in vector.coefficients.items() if not selected(key)},
                          vector.denominator)
            updated = combine([updated, kept])
        vector = updated
        physical = {key: sqrt_upper(weight)*budget[key] for key in errors}
        for key in errors:
            errors[key] += physical[key]
        stages.append({"operation": name, **physical})

    advance("initial_kick", coordinate, -F(sign, 4))
    advance("first_evolution", hamiltonian, interval)
    positive, ep = flow(coordinate, vector, -F(1, 2), degree, substeps, bits)
    negative, em = flow(coordinate, vector, F(1, 2), degree, substeps, bits)
    vector, rounding = quantize(_read_combine(positive, negative), bits)
    weight /= 2
    read_budget = {key: sqrt_upper(weight*(ep[key]**2+em[key]**2))+sqrt_upper(weight)*rounding
                   for key in errors}
    for key in errors:
        errors[key] += read_budget[key]
    stages.append({"operation": "first_read", **read_budget})
    advance("conditional_feedback_kick", coordinate, -F(1, 4), lambda key: key[7] == 1)
    advance("second_evolution", hamiltonian, interval)
    terminal = terminal_probabilities(vector, weight, errors["unrestricted_vector_error_upper"])
    result = {"schema": "oph.scalar_observations.factorized_free.v1",
              "observed_mode": "cos_x", "observed_occupation_dimension": dimension,
              "coupling": 0, "preparation_sign": sign, "interval": interval,
              "degree": degree, "substeps_per_flow": substeps, "rounding_bits": bits,
              "raw_state_squared_weight": weight, **errors, "stages": stages,
              "terminal_effect_certificate": terminal,
              "scope": {"exact_gradient": True, "same_supplied_seven_mode_record_law": True,
                        "unobserved_modes_traced_out_exactly": 6,
                        "product_blank_and_observed_mode_only_operations_required": True,
                        "seven_mode_quantum_output_state_certified": False,
                        "every_record_retained": True, "branch_normalization_or_postselection": False,
                        "rational_residual_and_rounding_arithmetic": True,
                        "source_execution_or_empirical_measurement": False}}
    if return_state:
        result["state"] = vector
    return result
