"""Independent high-precision checks of every complex generator entry.

These tests use exact binary input values and high-precision square roots,
without BLAS or an additive comparison tolerance. They validate the interval
implementation on representative configurations; the finite tests themselves
are not a universal arithmetic proof.
"""
from fractions import Fraction

import mpmath
import numpy as np
import pytest

from .solver import generator


@pytest.mark.parametrize('sectors,flip,imaginary,pointer', [
    (2, False, False, False),
    (2, True, False, False),
    (2, False, True, False),
    (4, False, False, True),
    (4, True, True, True),
])
@pytest.mark.parametrize('duration', [.037, Fraction(1, 37)])
def test_generator_intervals_enclose_independent_high_precision_entries(
        sectors, flip, imaginary, pointer, duration):
    ctx = mpmath.mp.clone()
    ctx.dps = 100

    def exact(value):
        # Fraction(float) preserves its exact binary rational, including the
        # signed state entries. Rational cell durations are not rounded first.
        rational = Fraction(value)
        return ctx.mpf(rational.numerator)/rational.denominator

    dimension = 3
    size = dimension*sectors
    rng = np.random.default_rng(129+sectors+int(flip))
    raw = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
    rho = (raw+raw.conj().T)/2
    delta, force, spin_x, spin_z, gamma = 1.3, -.7, .2, .3, .1
    operation = generator(
        dimension, sectors=sectors, duration=duration, delta=delta,
        force=force, force_spin_flip=flip, force_imaginary=imaginary,
        spin_x=spin_x, spin_z=spin_z, gamma=gamma, pointer_sign=pointer)

    # Rebuild H from its occupation/spin definition, independently of the
    # sparse operation.terms and of operation.apply_float.
    hamiltonian = ctx.zeros(size)
    state = ctx.matrix([[ctx.mpc(exact(value.real), exact(value.imag))
                         for value in row] for row in rho])
    for i in range(size):
        sector, occupation = divmod(i, dimension)
        hamiltonian[i, i] = exact(delta)*occupation+exact(spin_z)*(1-2*(sector % 2))
        hamiltonian[i, (sector ^ 1)*dimension+occupation] += exact(spin_x)
        for j in range(size):
            source_sector, source_occupation = divmod(j, dimension)
            if (source_sector == (sector ^ int(flip))
                    and abs(source_occupation-occupation) == 1):
                value = exact(force)*ctx.sqrt(max(occupation, source_occupation))
                if pointer:
                    value *= 1-2*(sector//2)
                if imaginary:
                    value *= ctx.j*(occupation-source_occupation)
                hamiltonian[i, j] += value

    direct = -ctx.j*(hamiltonian*state-state*hamiltonian)*exact(duration)
    for i in range(size):
        for j in range(size):
            difference = (i % dimension)-(j % dimension)
            direct[i, j] -= exact(gamma)*exact(duration)/2*difference**2*state[i, j]

    real, imag = operation.apply_interval(rho)
    for i in range(size):
        for j in range(size):
            assert exact(real[0][i, j]) <= direct[i, j].real <= exact(real[1][i, j])
            assert exact(imag[0][i, j]) <= direct[i, j].imag <= exact(imag[1][i, j])
