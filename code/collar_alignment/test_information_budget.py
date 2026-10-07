"""Independent witnesses and variational checks for the full alignment budget."""

import itertools

import numpy as np
import pytest
from scipy.linalg import logm

from collar_alignment.msa_characterizations import (
    alignment_information_budget, bell_counterexample, collar_cmi,
    entropic_alignment_defect, is_ec_aligned, random_generic_blocks,
)


@pytest.mark.parametrize("pair,term", (((0, 3), "collar"), ((1, 2), "middle"),
                                     ((0, 2), "left_leakage"), ((1, 3), "right_leakage")))
def test_each_of_the_four_obstructions_can_be_present_alone(pair, term):
    p = [1/8 if bits[pair[0]] == bits[pair[1]] else 0
         for bits in itertools.product((0, 1), repeat=4)]
    result = alignment_information_budget([(1., np.diag(p), (2, 2, 2, 2))])
    weighted = result["weighted"]
    for key in ("collar", "middle", "left_leakage", "right_leakage"):
        assert weighted[key] == pytest.approx(np.log(2) if key == term else 0, abs=1e-14)
    assert weighted["alignment"] == pytest.approx(np.log(2), abs=1e-14)


def test_bell_counterexample_resolves_into_two_endpoint_leakages():
    result = alignment_information_budget(bell_counterexample())
    values = result["weighted"]
    assert values["collar"] == values["middle"] == 0
    assert values["left_leakage"] == pytest.approx(2*np.log(2), abs=1e-14)
    assert values["right_leakage"] == pytest.approx(2*np.log(2), abs=1e-14)
    assert values["alignment"] == pytest.approx(4*np.log(2), abs=1e-14)


@pytest.mark.parametrize("seed", range(4))
def test_budget_is_the_exact_relative_entropy_projection_cost(seed):
    blocks = random_generic_blocks(np.random.default_rng(700+seed), [(2, 2, 2, 2)])
    rho = blocks[0][1]
    left = np.einsum('arbr->ab', rho.reshape(4, 4, 4, 4))
    right = np.einsum('aras->rs', rho.reshape(4, 4, 4, 4))
    projected = np.kron(left, right)/np.trace(rho).real
    cost = np.trace(rho@(logm(rho)-logm(projected))).real
    result = alignment_information_budget(blocks)
    assert result["weighted"]["alignment"] == pytest.approx(cost, abs=2e-14)
    assert abs(result["sectors"][0]["chain_rule_residual"]) < 2e-14
    distance = np.abs(np.linalg.eigvalsh(rho-projected)).sum()/2
    assert distance <= result["trace_distance_bound"]+1e-14


def test_variational_identity_includes_wrong_marginals_and_wrong_center():
    # Direct independent logarithms; the minimizer is not supplied to the
    # evaluator. Perturb both sector weights and product marginals.
    states = [np.diag([.375, .125, .125, .375]), np.diag([.25, .25, .375, .125])]
    p, q = [.25, .75], [.625, .375]
    trial_l = [np.diag([.375, .625]), np.diag([.75, .25])]
    trial_r = [np.diag([.625, .375]), np.diag([.25, .75])]
    kl = lambda a, b: float(np.trace(a@(logm(a)-logm(b))).real)
    result = alignment_information_budget([(w, r, (2, 1, 1, 2)) for w, r in zip(p, states)])
    correction = sum(w*np.log(w/v) for w, v in zip(p, q))
    direct = correction
    for w, rho, a, b in zip(p, states, trial_l, trial_r):
        left = np.diag(np.diag(rho).reshape(2, 2).sum(axis=1))
        right = np.diag(np.diag(rho).reshape(2, 2).sum(axis=0))
        correction += w*(kl(left, a)+kl(right, b))
        direct += w*kl(rho, np.kron(a, b))
    assert direct == pytest.approx(result["weighted"]["alignment"]+correction, abs=2e-14)
    assert correction > 0


def test_small_sector_average_does_not_erase_worst_sector_misalignment():
    correlated = np.diag([.5, 0., 0., .5])
    blocks = [(1., np.eye(4)/4, (2, 1, 1, 2)), (1e-100, correlated, (2, 1, 1, 2))]
    result = alignment_information_budget(blocks)
    assert result["weighted"]["alignment"] == pytest.approx(1e-100*np.log(2), rel=1e-14, abs=0)
    assert result["maximum_alignment"] == entropic_alignment_defect(blocks)
    assert not is_ec_aligned(blocks)
    assert collar_cmi(blocks) == result["weighted"]["collar"]


def test_collar_wrapper_does_not_erase_subnormal_hermitian_roundoff():
    rho = np.eye(4)/4
    rho[0, 3] = np.nextafter(0., 1.)
    with pytest.raises(ValueError, match="precision|underflow|unresolved"):
        is_ec_aligned([(1., rho, (2, 1, 1, 2))])


@pytest.mark.parametrize("target", ("state", "weight"))
def test_collar_validation_does_not_strip_missing_data(target):
    rho = np.ma.array(np.eye(4)/4, mask=False)
    weight = 1.
    if target == "state":
        rho.mask[0, 1] = True
    else:
        weight = np.ma.masked
    with pytest.raises(ValueError, match="masked|missing"):
        alignment_information_budget([(weight, rho, (2, 1, 1, 2))])
