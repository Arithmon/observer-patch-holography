"""Seven-dimensional contact control, with an exact rational pole enclosure.

This reduction only covers the contact-generated cyclic subspace at L=4,
K=0. It does not replace the full field, its other inputs or local detectors.
"""
from fractions import Fraction as F
from math import comb
import numpy as np


def _parameters(speed, interaction):
    if type(speed) is not F or type(interaction) is not F:
        raise ValueError('exact Fraction couplings required')
    if speed <= 0 or interaction <= 0 or 12*speed**2 >= interaction**2:
        raise ValueError('positive couplings with U^2 > 12 c^2 required')
    return speed, interaction


def contact_resolvent(energy, speed=F(1, 8)):
    """Exact <d,(E-T)^-1 d> above the complete free L=4, K=0 spectrum."""
    if type(energy) is not F or type(speed) is not F:
        raise ValueError('exact Fraction energy and speed required')
    if speed <= 0 or energy <= 0 or energy**2 <= 12*speed**2:
        raise ValueError('positive energy above 2 sqrt(3) c required')
    return F(1, 8)/energy + sum(
        (F(comb(3, r), 8)*energy/(energy**2-4*speed**2*r)
         for r in (1, 2, 3)), F())


def pole_interval(speed=F(1, 8), interaction=F(1), iterations=64):
    """Exact bisection certificate; no floating eigensolver or square root.

The resolvent is strictly decreasing above the free spectrum. At U its
value exceeds 1/U. At U+4c every denominator E-lambda exceeds U, so its
value is less than 1/U. Each returned endpoint retains the corresponding
exact secular sign.
"""
    c, u = _parameters(speed, interaction)
    if type(iterations) is not int or not 0 <= iterations <= 256:
        raise ValueError('integer bisection count 0..256 required')
    lo, hi = u, u+4*c
    for _ in range(iterations):
        middle = (lo+hi)/2
        if u*contact_resolvent(middle, c) >= 1:
            lo = middle
        else:
            hi = middle
    return lo, hi


def matrices(speed=F(1, 8), interaction=F(1)):
    """Numerical representation of the analytically exact seven-state model.

Free atoms ordered by energy have weights (1,3,3,2,3,3,1)/16.
The on-site contact vector is their positive square root. All seven states,
including the six interacting scattering eigenstates, remain present.
"""
    c, u = _parameters(speed, interaction)
    free = 2*float(c)*np.array([-np.sqrt(3), -np.sqrt(2), -1, 0,
                               1, np.sqrt(2), np.sqrt(3)])
    weights = np.array([1, 3, 3, 2, 3, 3, 1], float)/16
    contact = np.sqrt(weights)
    hamiltonian = np.diag(free)+float(u)*np.outer(contact, contact)
    return free, contact, hamiltonian
