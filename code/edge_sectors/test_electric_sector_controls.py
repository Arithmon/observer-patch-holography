"""Independent original-incidence controls for small and rare electric sectors."""

import itertools
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge_sectors.abelian_ground_state import zn_edge_distribution


# Oriented boundaries in (00x,00y,01x,01y,10x,10y,11x,11y) order.
_FACES = np.array([[1, -1, -1, 0, 0, 1, 0, 0],
                   [-1, 0, 1, -1, 0, 0, 0, 1],
                   [0, 1, 0, 0, 1, -1, -1, 0],
                   [0, 0, 0, 1, -1, 0, 1, -1]])


def _original_incidence_oracle(n, h, *, dps=100, return_high_precision=False):
    """Diagonalize the independently enumerated electric-loop Hamiltonian.

    All eight original links contribute their electric energy. Magnetic
    hopping adds or subtracts each oriented face boundary, and the edge
    charge uses its original three-link restricted star. No producer's
    potential chart, trigonometric simplification or sparse solver is used.
    """
    ctx = mpmath.mp.clone()
    ctx.dps = dps
    fields = sorted({tuple(np.array([0, *potentials]) @ _FACES % n)
                     for potentials in itertools.product(range(n), repeat=3)})
    indices = {field: row for row, field in enumerate(fields)}
    matrix = ctx.zeros(len(fields))
    original_h = ctx.mpf(h)
    for row, field in enumerate(fields):
        matrix[row, row] = sum(2*original_h*ctx.sin(ctx.pi*int(q)/n)**2
                               for q in field)
        for face in _FACES:
            for sign in (-1, 1):
                column = indices[tuple((np.array(field)+sign*face) % n)]
                matrix[row, column] -= ctx.mpf(".5")
    _, vectors = ctx.eigsy(matrix)
    probabilities = [ctx.mpf(0)]*n
    for row, field in enumerate(fields):
        charge = int((field[0]+field[1]-field[3]) % n)
        probabilities[charge] += vectors[row, 0]**2
    if return_high_precision:
        return ctx, tuple(probabilities)
    return np.array([float(value) for value in probabilities])


@pytest.mark.parametrize("n,h", [(3, 1e-40), (3, 1e6), (4, 10.)])
def test_each_sector_matches_independent_original_link_incidence(n, h):
    expected = _original_incidence_oracle(n, h)
    measured = zn_edge_distribution(n, h)
    assert np.all(expected > 0)
    assert measured == pytest.approx(expected, rel=2e-12, abs=0)
