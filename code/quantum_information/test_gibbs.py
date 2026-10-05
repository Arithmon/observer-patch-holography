"""Independent full-matrix and variational controls for sector Gibbs states."""

import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import block_diag, expm, logm

from quantum_information import gibbs_sectors
from maxent.information_projection import gibbs_state


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


def embedded(blocks):
    return block_diag(*(p*rho for p, rho in blocks))


def matrix_gibbs(ham, beta):
    # Independent scaling-and-squaring exponential, no production eigensolve.
    unnormalized = expm(-beta*ham)
    return unnormalized/np.trace(unnormalized)


@pytest.mark.parametrize("seed", range(10))
def test_unequal_sectors_equal_one_full_matrix_exponential(seed):
    rng = np.random.default_rng(seed)
    operators = []
    for d in [1, 2, 3, 4]:
        a = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
        operators.append((a+a.conj().T)/4)
    energies = rng.normal(size=4)
    beta = rng.uniform(-.6, .6)
    full = block_diag(*(h+e*np.eye(len(h)) for h, e in zip(operators, energies)))
    actual = embedded(gibbs_sectors(operators, energies, beta))
    np.testing.assert_allclose(actual, matrix_gibbs(full, beta), atol=3e-15, rtol=0)


@pytest.mark.parametrize("scale", [1e-300, 1e-200, 1., 1e200, 1e300, -1e200])
def test_temperature_and_energy_units_transform_together(scale):
    operators = [scale*(X+.5*Y), scale*Z]
    actual = embedded(gibbs_sectors(operators, [0., .25*scale], beta=.4/scale))
    expected = matrix_gibbs(block_diag(X+.5*Y, Z+.25*I), .4)
    np.testing.assert_allclose(actual, expected, atol=3e-15, rtol=0)


@pytest.mark.parametrize("scale", [1e-308, 1e-310])
def test_subnormal_operator_units_retain_a_resolved_thermal_interaction(scale):
    # A complex reciprocal 1/scale would overflow before beta compensates.
    beta = 1e308
    (p, rho), = gibbs_sectors([scale*X], [0.], beta)
    expected = -.5*np.tanh(beta*scale)
    assert p == 1
    assert rho[0, 1].real == pytest.approx(expected, rel=2e-14, abs=0)
    assert rho[0, 1].real != 0


def test_zero_temperature_multiplier_weights_dimension_not_sector_count():
    operators = [X, np.diag([1., 2., 3.]), np.zeros((1, 1))]
    actual = embedded(gibbs_sectors(operators, [1e308, -1e308, 0.], beta=0.))
    np.testing.assert_allclose(actual, np.eye(6)/6, atol=1e-16, rtol=0)


def test_opposite_offsets_cancel_before_subtracting_relative_log_weights():
    # The first sector retains its off-diagonal spectrum after +1e20/-1e20.
    actual = embedded(gibbs_sectors([1e20*I+X, Y+Z], [-1e20, .5], .7))
    expected = matrix_gibbs(block_diag(X, Y+Z+.5*I), .7)
    np.testing.assert_allclose(actual, expected, atol=3e-15, rtol=0)


def test_small_mean_is_retained_when_large_central_offsets_cancel():
    # diag(0,2) has trace mean 1. Adding that mean to 1e20 too early would
    # erase it and change this partition ratio by e.
    blocks = gibbs_sectors([np.diag([0., 2.]), X], [1e20, 1e20])
    assert blocks[0][0]/blocks[1][0] == pytest.approx(np.exp(-1), rel=2e-14)


def test_representable_independent_sector_rechart_preserves_full_state():
    operators = [X+.5*Z, .7*Y-.25*Z]
    reference = matrix_gibbs(block_diag(operators[0]+.75*I, operators[1]-.5*I), .8)
    for offsets in [[2., -4.], [-32., 16.]]:
        actual = embedded(gibbs_sectors([h+s*I for h, s in zip(operators, offsets)],
                                        [.75-offsets[0], -.5-offsets[1]], .8))
        np.testing.assert_allclose(actual, reference, atol=3e-15, rtol=0)


def test_relative_central_energy_changes_weight_but_not_conditional_state():
    blocks = gibbs_sectors([X, X], [0., 2.], beta=.3)
    assert blocks[0][0]/blocks[1][0] == pytest.approx(np.exp(.6), rel=2e-14)
    for _, rho in blocks:
        np.testing.assert_allclose(rho, (I-np.tanh(.3)*X)/2, atol=2e-15, rtol=0)


def test_splitting_and_permuting_a_direct_sum_sector_preserves_full_state():
    operators = [X, .5*Y, .7*Z]
    separate = gibbs_sectors(operators, [0., 0., .3], .8)
    joined = gibbs_sectors([block_diag(X, .5*Y), .7*Z], [0., .3], .8)
    np.testing.assert_allclose(embedded(separate), embedded(joined), atol=3e-15, rtol=0)
    permutation = [2, 0, 1]
    permuted = gibbs_sectors([operators[i] for i in permutation], [.3, 0., 0.], .8)
    for actual, i in zip(permuted, permutation):
        assert actual[0] == pytest.approx(separate[i][0], abs=3e-15)
        np.testing.assert_allclose(actual[1], separate[i][1], atol=3e-15, rtol=0)


def test_maxent_and_sector_constructor_share_the_same_finite_gibbs_family():
    h = .3*X+.2*Y-.4*Z
    (p, sector), = gibbs_sectors([1e20*I+h], [1e20], .5)
    # A huge offset added to a diagonal term already destroys that term in
    # the caller's literal array. This test uses the actual supplied operator.
    stored_h = (1e20*I+h)-1e20*I
    expected = matrix_gibbs(stored_h, .5)
    single, _ = gibbs_state([1e20*I+h], [.5])
    assert p == 1
    np.testing.assert_allclose(sector, expected, atol=2e-15, rtol=0)
    np.testing.assert_allclose(single, expected, atol=2e-15, rtol=0)


@pytest.mark.parametrize("seed", range(4))
def test_relative_entropy_splits_into_center_and_conditional_excess(seed):
    rng = np.random.default_rng(seed)
    operators, energies, beta = [X+.2*Z, .5*Y], [.3, -.2], .7
    full_h = block_diag(*(h+e*I for h, e in zip(operators, energies)))
    optimum = gibbs_sectors(operators, energies, beta)
    target = embedded(optimum)
    assert np.linalg.norm(target-matrix_gibbs(full_h, beta)) < 3e-15
    weights = rng.dirichlet([2., 3.])
    trial = [matrix_gibbs(rng.uniform(-1, 1)*Z, 1.), matrix_gibbs(.3*X, 1.)]
    candidate = block_diag(*(p*rho for p, rho in zip(weights, trial)))
    entropy = lambda a, b: float(np.trace(a@(logm(a)-logm(b))).real)
    center_excess = sum(p*np.log(p/q) for p, (q, _) in zip(weights, optimum))
    conditional_excess = sum(p*entropy(rho, tau)
                             for p, rho, (_, tau) in zip(weights, trial, optimum))
    total_excess = entropy(candidate, target)
    assert center_excess >= 0
    assert conditional_excess > .01
    assert total_excess == pytest.approx(center_excess+conditional_excess, abs=1e-14)
    free_energy = lambda a: float(np.trace(a@(beta*full_h+logm(a))).real)
    assert free_energy(candidate)-free_energy(target) == pytest.approx(total_excess, abs=1e-14)


@pytest.mark.parametrize("bad", [
    np.diag(np.array([2**53, 2**53+1], dtype=np.int64)),
    [[2**60, 0.], [0., 2**60+1]],
])
def test_lossy_integer_conversion_cannot_erase_a_level_gap(bad):
    for call in [lambda: gibbs_sectors([bad], [0.]),
                 lambda: gibbs_state([bad], [1.])]:
        with pytest.raises(ValueError, match="conversion"):
            call()


@pytest.mark.parametrize("which", ["central", "beta"])
def test_lossy_integer_parameter_cannot_change_relative_energies(which):
    with pytest.raises(ValueError, match="conversion"):
        gibbs_sectors([X, X], [2**53, 2**53+1] if which == "central" else [0., 0.],
                       2**53+1 if which == "beta" else 1.)


@pytest.mark.parametrize("case", ["spectrum", "sector", "joint", "dense"])
def test_positive_mass_cannot_be_clipped_away(case):
    hams, energies = {
        "spectrum": ([np.diag([0., 1000.])], [0.]),
        "sector": ([I, I], [0., 1000.]),
        # Both sector probabilities and conditional eigenvalues survive;
        # one product does not. Local support checks alone are insufficient.
        "joint": ([np.zeros((1, 1)), np.diag([0., 400.])], [0., 400.]),
        "dense": ([30*X], [0.]),
    }[case]
    with pytest.raises(ValueError, match="underflow|faithful"):
        gibbs_sectors(hams, energies)


def test_huge_relative_energy_does_not_return_nan_or_drop_a_sector():
    with pytest.raises(ValueError, match="underflow"):
        gibbs_sectors([np.zeros((1, 1))]*2, [1e308, -1e308], beta=1e308)


@pytest.mark.parametrize("bad", [[], [np.eye(2), np.eye(3)], [[np.nan]],
                                 [[np.inf]], [[1j]], np.ones((2, 3))])
def test_invalid_operator_fails_even_at_beta_zero(bad):
    with pytest.raises(ValueError):
        gibbs_sectors([bad], [0.], beta=0.)


def test_validation_does_not_depend_on_python_assertions(tmp_path):
    source = '''
import numpy as np
from quantum_information import gibbs_sectors
bad = [lambda: gibbs_sectors([], []),
       lambda: gibbs_sectors([[[0., 1e-200], [0., 0.]]], [0.], 1e200),
       lambda: gibbs_sectors([np.eye(2)], [True]),
       lambda: gibbs_sectors([np.diag([0.,1000.])], [0.])]
for call in bad:
    try: call()
    except ValueError: pass
    else: raise RuntimeError('invalid input accepted with assertions disabled')
'''
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    subprocess.run([sys.executable, "-O", "-c", source], cwd=tmp_path, env=env,
                   check=True, capture_output=True, text=True)
