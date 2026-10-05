"""Analytic controls for Gibbs energy origins and sector normalization (#1033)."""

import numpy as np
import pytest

from collar_alignment.msa_characterizations import (
    entropic_alignment_defect, gibbs_blocks, is_ec_aligned,
    modular_splitting_defect, takesaki_defect,
)


X = np.array([[0., 1.], [1., 0.]])
XX = np.kron(X, X)
DIMS = (1, 2, 2, 1)


@pytest.mark.parametrize("origin", ["central", "matrix"])
def test_energy_origin_cannot_hide_a_cross_cut_interaction(origin):
    # XX has eigenvalues +/-1, each twice. The state and its mutual
    # information follow exactly from XX**2=I, without a Gibbs oracle.
    beta = .4
    ham = XX + (1e20*np.eye(4) if origin == "matrix" else 0)
    energy = 1e20 if origin == "central" else 0
    blocks = gibbs_blocks([(ham, DIMS)], [energy], beta)
    expected = (np.eye(4)-np.tanh(beta)*XX)/4
    np.testing.assert_allclose(blocks[0][1], expected, atol=2e-15, rtol=0)
    assert entropic_alignment_defect(blocks) == pytest.approx(
        beta*np.tanh(beta)-np.log(np.cosh(beta)), abs=2e-14)
    assert modular_splitting_defect(blocks) == pytest.approx(beta, abs=2e-14)
    assert takesaki_defect(blocks) > .1
    assert not is_ec_aligned(blocks)


def test_common_central_origin_retains_distinct_sector_partition_functions():
    sectors = [(np.diag([0., 1.]), (1, 1, 1, 2)),
               (np.zeros((1, 1)), (1, 1, 1, 1))]
    blocks = gibbs_blocks(sectors, [1e20, 1e20])
    z = 2+np.exp(-1)
    np.testing.assert_allclose([p for p, _, _ in blocks],
                               [(1+np.exp(-1))/z, 1/z], atol=2e-15, rtol=0)
    np.testing.assert_allclose(blocks[0][1],
                               np.diag([1., np.exp(-1)])/(1+np.exp(-1)),
                               atol=2e-15, rtol=0)


def test_cancelling_sector_offsets_retain_small_trace_energy():
    # The small relative trace energy must survive cancellation of the
    # second sector's large matrix and central offsets.
    sectors = [(np.diag([0., 2.]), (1, 1, 1, 2)),
               (1e20*np.eye(2)+X, (1, 1, 1, 2))]
    blocks = gibbs_blocks(sectors, [0., -1e20])
    # Spectra are (0,2) and (-1,1); Z_0/Z_1=exp(-1).
    assert blocks[0][0]/blocks[1][0] == pytest.approx(np.exp(-1), rel=2e-14)
    np.testing.assert_allclose(blocks[1][1], (np.eye(2)-np.tanh(1)*X)/2,
                               atol=2e-15, rtol=0)


def test_unobservable_common_origin_may_exceed_float_energy_range():
    # beta*e overflows, but beta*H and every output probability are moderate.
    blocks = gibbs_blocks([(1e-200*X, (1, 1, 1, 2))], [1e200], beta=1e200)
    np.testing.assert_allclose(blocks[0][1], (np.eye(2)-np.tanh(1)*X)/2,
                               atol=2e-15, rtol=0)


def test_small_units_cannot_hide_a_nonhermitian_hamiltonian():
    with pytest.raises(ValueError, match="Hermitian"):
        gibbs_blocks([(np.array([[0., 1e-200], [0., 0.]]),
                       (1, 1, 1, 2))], [0], beta=1e200)


@pytest.mark.parametrize("bad", [[[False, 1.], [1., 0.]],
                                 [["0", "1"], ["1", "0"]]])
def test_sector_operator_cannot_be_silently_coerced(bad):
    with pytest.raises(ValueError):
        gibbs_blocks([(bad, (1, 1, 1, 2))], [0])
