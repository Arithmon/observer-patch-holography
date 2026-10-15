"""Sparse polynomial paths with integrated, outward-enclosed full defects.

The numerical path need not be positive. Trace-norm contractivity of the
reference dynamics, not stability or positivity of a numerical integrator,
controls its error. Coefficients are treated as exact binary matrices when
their residual is enclosed. Every cell has its own endpoint-jump budget.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import numpy as np

from . import intervals as iv


def _fraction(value):
    if isinstance(value, bool) or not math.isfinite(value):
        raise ValueError('finite real parameter required')
    return Fraction(value)


@dataclass
class Generator:
    dimension: int
    sectors: int
    terms: list
    noise: tuple
    boundary: tuple

    def apply_float(self, state):
        left = np.zeros_like(state)
        for rows, cols, real, imag in self.terms:
            value = (np.asarray(real[0])+np.asarray(real[1]))/2
            value = value+1j*(np.asarray(imag[0])+np.asarray(imag[1]))/2
            left[rows] += value[:, None]*state[cols]
        noise = (self.noise[0]+self.noise[1])/2
        return -1j*(left-left.conj().T)-noise*state

    def apply_interval(self, state):
        re, im = iv.complex_exact(state)
        shape = state.shape
        lr = (np.zeros(shape), np.zeros(shape))
        li = (np.zeros(shape), np.zeros(shape))
        for rows, cols, hr, hi in self.terms:
            xr, xi = (re[0][cols], re[1][cols]), (im[0][cols], im[1][cols])
            hr, hi = tuple(np.asarray(x)[:, None] for x in hr), tuple(np.asarray(x)[:, None] for x in hi)
            ar = iv.sub(iv.mul(hr, xr), iv.mul(hi, xi))
            ai = iv.add(iv.mul(hr, xi), iv.mul(hi, xr))
            for accumulator, value in ((lr, ar), (li, ai)):
                lo, upper = iv.add((accumulator[0][rows], accumulator[1][rows]), value)
                accumulator[0][rows], accumulator[1][rows] = lo, upper
        real = iv.sub(iv.add(li, (li[0].T, li[1].T)), iv.mul(self.noise, re))
        imag = iv.sub(iv.sub((lr[0].T, lr[1].T), lr), iv.mul(self.noise, im))
        return real, imag


def generator(dimension, *, sectors=2, duration=1., delta=0., force=0.,
              force_spin_flip=False, force_imaginary=False, spin_x=0.,
              spin_z=0., gamma=0., pointer_sign=False):
    """H=delta*N+force*Q (+ spin controls), scaled by exact duration.

For force_imaginary, the oscillator term is i*force*(a†-a).
For pointer_sign, sectors are [pointer0/spin0, pointer0/spin1,
pointer1/spin0, pointer1/spin1] and the force changes sign by pointer.
All other inputs denote exact binary rationals. The finite ladder square
roots and the time scale are enclosed, not silently rounded into the target.
"""
    if type(dimension) is not int or dimension < 2 or sectors not in (2, 4):
        raise ValueError('dimension >=2 and two or four finite sectors required')
    if any(type(flag) is not bool for flag in (force_spin_flip, force_imaginary, pointer_sign)):
        raise ValueError('boolean force selector flags required')
    duration = _fraction(duration)
    values = [delta, force, spin_x, spin_z, gamma]
    if duration < 0 or any(not math.isfinite(x) for x in values) or gamma < 0:
        raise ValueError('finite parameters and nonnegative duration/noise required')
    def scaled(value):
        return iv.rational(_fraction(value)*duration)
    size = dimension*sectors
    n = np.tile(np.arange(dimension), sectors)
    sector = np.repeat(np.arange(sectors), dimension)
    rows = np.arange(size)
    terms = []
    diagonal = iv.add(iv.scale(scaled(delta), n),
                      iv.scale(scaled(spin_z), 1-2*(sector % 2)))
    terms.append((rows, rows, diagonal, iv.exact(np.zeros(size))))
    if spin_x:
        columns = (sector ^ 1)*dimension+n
        terms.append((rows, columns, iv.scale(scaled(spin_x), np.ones(size)), iv.exact(np.zeros(size))))
    for offset in (-1, 1):
        valid = (n+offset >= 0) & (n+offset < dimension)
        dest = rows[valid]
        source_sector = sector[valid] ^ int(force_spin_flip)
        source = source_sector*dimension+n[valid]+offset
        root = iv.square_root(iv.exact(np.maximum(n[valid], n[valid]+offset)))
        value = iv.mul(root, scaled(force))
        if pointer_sign:
            value = iv.scale(value, 1-2*(sector[valid]//2))
        if force_imaginary:
            value = iv.scale(value, -offset)
            terms.append((dest, source, iv.exact(np.zeros(len(dest))), value))
        else:
            terms.append((dest, source, value, iv.exact(np.zeros(len(dest)))))
    differences = n[:, None]-n[None, :]
    noise = iv.scale(iv.rational(_fraction(gamma)*duration/2), differences**2)
    boundary = iv.mul(iv.rational(abs(_fraction(force))*duration),
                      iv.square_root(iv.exact(dimension)))
    return Generator(dimension, sectors, terms, noise, boundary)


def cell(state, operation, *, order=12):
    """Advance one normalized cell, bounding its complete infinite defect."""
    state = np.asarray(state, dtype=complex)
    if (type(order) is not int or order < 2 or state.shape !=
            (operation.dimension*operation.sectors,)*2 or
            not np.isfinite(state).all() or not np.array_equal(state, state.conj().T)):
        raise ValueError('finite exactly Hermitian state, valid shape and order >=2 required')
    coefficients = [state]
    for k in range(1, order+1):
        value = operation.apply_float(coefficients[-1])/k
        # Exact Hermiticity is a requirement for the contraction inequality.
        coefficients.append((value+value.conj().T)/2)
    finite, boundary = [], []
    top = np.arange(operation.dimension-1, len(state), operation.dimension)
    sum_interval = iv.complex_exact(np.zeros_like(state))
    for k, value in enumerate(coefficients):
        image = operation.apply_interval(value)
        derivative = (iv.complex_scale(iv.complex_exact(coefficients[k+1]), iv.exact(k+1))
                      if k < order else iv.complex_exact(np.zeros_like(state)))
        residual = iv.complex_sub(derivative, image)
        finite.append(float(iv.up(iv.trace_upper(*residual)/(2*(k+1)))))
        row = iv.complex_exact(value[top])
        product_upper = iv.up(operation.boundary[1]*iv.trace_upper(*row))
        boundary.append(float(iv.up(product_upper/(k+1))))
        sum_interval = (iv.add(sum_interval[0], iv.exact(value.real)),
                        iv.add(sum_interval[1], iv.exact(value.imag)))
    endpoint = sum(coefficients)
    endpoint = (endpoint+endpoint.conj().T)/2
    jump = iv.complex_sub(iv.complex_exact(endpoint), sum_interval)
    budget = {'finite_defect': iv.sum_up(finite),
              'occupation_defect': iv.sum_up(boundary),
              'endpoint_jump': float(iv.up(iv.trace_upper(*jump)/2))}
    return endpoint, budget


def evolve(state, *, cells, order=12, **parameters):
    if type(cells) is not int or cells < 1:
        raise ValueError('positive integer cell count required')
    duration = _fraction(parameters.pop('duration', 1.))
    operation = generator(len(state)//parameters.get('sectors', 2),
                          duration=duration/cells, **parameters)
    budgets = []
    for _ in range(cells):
        state, budget = cell(state, operation, order=order)
        budgets.append(budget)
    return state, {name: iv.sum_up([budget[name] for budget in budgets]) for name in budgets[0]}
