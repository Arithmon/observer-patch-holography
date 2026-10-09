"""Original-input controls for rare finite-group ground-state sectors."""

import mpmath
import numpy as np
import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from edge_sectors.heat_kernel_holdout_validation import (
    report_s3, s3_edge_distribution, zn_edge_distribution,
)


def _s3_secular_probabilities(h):
    """Independent positive root, with no matrix eigensolver or fit reused.

    The supplied Hamiltonian has ground energy -x. Its first two rows
    give positive amplitudes (1, x/(24h+x), x); the last row gives the
    scalar equation below. Perron-Frobenius selects its unique positive
    root, which lies in (0, 2] for h >= 0.
    """
    ctx = mpmath.mp.clone()
    ctx.dps = 120
    original_h = ctx.mpf(h)
    lower, upper = ctx.mpf(0), ctx.mpf(2)
    for _ in range(500):
        x = (lower+upper)/2
        residual = (12*original_h-1+x)*x-1-x/(24*original_h+x)
        if residual > 0:
            upper = x
        else:
            lower = x
    x = (lower+upper)/2
    amplitudes = (ctx.mpf(1), x/(24*original_h+x), x)
    norm = sum(value*value for value in amplitudes)
    return {name: float(value*value/norm) for name, value in
            zip(("triv", "sign", "std"), amplitudes)}


def _z2_electric_loop_probabilities(h):
    """Eight-state physical loop Hamiltonian, separately assembled.

    Bits label (00x,00y,01x,01y,10x,10y,11x,11y). Plaquette boundaries
    generate eight distinct loops, with their product the identity.
    The nonnegative ground state occupies the trivial Gauss and winding
    sectors. Removing the common -8h-20 energy leaves 2h times the number
    of occupied links on the diagonal and -1 for each plaquette flip.
    Edge charge is the parity on links (00x,00y,01y), measured as a sum
    of squared electric amplitudes, not a difference of unit overlaps.
    """
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    faces = (0b00100111, 0b10001101, 0b01110010, 0b11011000)
    loops = sorted({(faces[0] if index & 1 else 0)
                    ^ (faces[1] if index & 2 else 0)
                    ^ (faces[2] if index & 4 else 0)
                    for index in range(8)})
    hamiltonian = ctx.zeros(8)
    for row, loop in enumerate(loops):
        hamiltonian[row, row] = 2*ctx.mpf(h)*loop.bit_count()
        for face in faces:
            hamiltonian[row, loops.index(loop ^ face)] -= 1
    _, vectors = ctx.eigsy(hamiltonian)
    probabilities = [ctx.mpf(0), ctx.mpf(0)]
    for row, loop in enumerate(loops):
        charge = (loop & 0b00001011).bit_count() % 2
        probabilities[charge] += vectors[row, 0]**2
    return np.array([float(value) for value in probabilities])


@pytest.mark.parametrize("h", [0., .5, 12., 100., 1e12, 1e15])
def test_s3_retains_every_representable_sector(h):
    expected = _s3_secular_probabilities(h)
    measured = s3_edge_distribution(h)
    for sector in ("triv", "sign", "std"):
        assert measured[sector] == pytest.approx(expected[sector], rel=2e-11, abs=0)
    assert sum(measured.values()) == pytest.approx(1., abs=2e-15)


def test_s3_report_retains_a_normal_range_held_out_probability(capsys):
    h = 1e15
    expected = _s3_secular_probabilities(h)
    assert expected["sign"] > np.finfo(float).tiny
    report_s3([h])
    rows = [line.split() for line in capsys.readouterr().out.splitlines()
            if line.startswith(f"{h:.2f}")]
    assert len(rows) == 1
    assert float(rows[0][2]) == pytest.approx(expected["sign"], rel=5e-5, abs=0)


@pytest.mark.parametrize("h", [.2, 2., 1e6, 1e7, 1e8])
def test_z2_readout_matches_positive_electric_sector_norms(h):
    expected = _z2_electric_loop_probabilities(h)
    measured = zn_edge_distribution(2, h)
    assert np.all(expected > 0)
    assert measured == pytest.approx(expected, rel=2e-8, abs=0)
