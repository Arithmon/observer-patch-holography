"""Exact common-reference characterization and native mixed-protocol checks."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from . import common_caps as c


@pytest.mark.parametrize("n", (3, 4, 5, 12))
def test_constructive_factorization_and_independent_spectator_ratios(n):
    relative = tuple(F(i+2, i+1) for i in range(n))
    weights = {(i, j): F(3, 7)*relative[i]*relative[j] for i, j in combinations(range(n), 2)}
    row = c.factor_characterization(n, weights)
    assert row["all_bare_two_port_algebras_invariant"]
    scale = F(row["factor_scale"]); factors = tuple(map(F, row["relative_factors"]))
    assert all(scale*factors[i]*factors[j] == w for (i, j), w in weights.items())
    for i, j in combinations(range(n), 2):
        ratios = {weights[tuple(sorted((i, k)))]/weights[tuple(sorted((j, k)))]
                  for k in range(n) if k not in (i, j)}
        assert ratios == {relative[i]/relative[j]}
    if n >= 4:
        changed = dict(weights); changed[0, 1] *= F(101, 100)
        assert not c.factor_characterization(n, changed)["all_bare_two_port_algebras_invariant"]


def test_arbitrary_three_port_weights_and_exact_native_countercontrols():
    row = c.factor_characterization(3, {(0, 1): F(2), (0, 2): F(7), (1, 2): F(11)})
    assert row["all_bare_two_port_algebras_invariant"]
    for sigma in (c.regional.source()[2], c.cut.completed_palm(c.CAPS[0]), c.reference()):
        row = c.factor_characterization(12, c.native_weights(sigma))
        assert not row["all_bare_two_port_algebras_invariant"]
        quad = row["noninvariance_four_port_witness"]
        assert len(set(quad["matching_products"])) > 1
    assert c.factor_characterization(12, c.native_weights((F(1, 66),)*66))["all_bare_two_port_algebras_invariant"]


def test_nonuniform_compatible_product_is_not_selected_by_connected_native_repair():
    # Mathematical compatible-clock control, explicitly not a sourced field.
    relative = (F(2), *(F(1) for _ in range(11)))
    states = c.regional.source()[0].states
    weights = tuple(relative[i]*relative[j] for x in states
                    for i, j in [tuple(p for p, bit in enumerate(x) if bit)])
    prior = np.array([float(w/sum(weights)) for w in weights])
    kernel = c.regional.source()[0]
    after = np.array(kernel.counts).T@prior/60
    assert np.linalg.norm(after-prior, 1) > .005
    for i, j in kernel.seams:
        if relative[i] != relative[j]:
            spectator = next(k for k in range(12) if k not in (i, j))
            assert relative[i]*relative[spectator] != relative[j]*relative[spectator]
            break
    else:
        raise AssertionError("nonconstant field has no native seam witness")


def test_actual_mixed_cap_source_and_all_1024_common_algebras():
    row = c.exact_audit()
    assert row["common_reference_cap_algebras_checked"] == 1024
    assert row["protected_intersection_modular_action"] == "identity"
    core, weights, schedule = c.cap_schedule()
    assert core == (0, 1) and sum(weights.values()) == 1
    assert schedule["native_seam_union"] == 20
    assert schedule["uniform_component_weight"] == "2/3"
    assert schedule["gap_lower_per_attempt"] == "1/600"
    assert all(e not in c.CORE for edge in weights for e in edge)
    assert set(c.reference()) == {F(17, 1125), F(23, 1500), F(1, 75)}


def test_mixed_native_complete_record_branch_covariance_and_cap_marginal_control():
    sigma = c.reference()
    kernel = c.regional.source()[0]
    index = {x: i for i, x in enumerate(kernel.states)}
    # Every branch commutes with one same reference, not different cap Gibbs
    # matrices. Labels C/E can split branches but cannot change this identity.
    for seam in c.cap_schedule()[1]:
        for coin in (0, 1):
            for i, x in enumerate(kernel.states):
                y, _ = c.cut.native.native(x, seam, coin)
                assert sigma[i] == sigma[index[y]]
    # The other cap's repair can change the first cap marginal. Only the
    # intersection is preserved throughout the mixed preparation.
    counts, scale = c.mixed_counts("C")
    cap_unit = None
    for unit in c.regional.basis_supports(c.CAPS[0], completed=False):
        vector = np.zeros(4356, np.int64)
        for i, j in unit:
            vector[66*i+j] = 1
        if np.any(counts.T@vector != scale*vector):
            cap_unit = unit; break
    assert cap_unit is not None
    for unit in c.regional.basis_supports(c.CORE, completed=False):
        vector = np.zeros(4356, np.int64)
        for i, j in unit:
            vector[66*i+j] = 1
        np.testing.assert_array_equal(counts.T@vector, scale*vector)


def test_unprotected_native_control_has_connected_tracial_selection():
    control = c.unprotected_control()
    assert control["schedule"]["protected_intersection"] == []
    assert control["schedule"]["native_seam_union"] == 28
    assert control["schedule"]["gap_lower_per_attempt"] == "1/1080"
    assert control["unique_fixed_binary_sector_state"] == "I66/66"


def test_finite_common_preparation_floor_and_reference_error():
    row = c.preparation_budget()
    assert row["mixed_cap_native_attempts_after_success"] == 18600
    assert row["actual_common_reference_modular_half_diamond_upper"] < F(99, 100000)
    assert row["minimum_work_eigenvalue_lower"] > F(1333, 100000)
    assert row["expected_native_preparation_attempts_upper"] < F(21042015, 1000)
    compact = c.preparation._compact(row)
    for key, value in row.items():
        if isinstance(value, F) and key.endswith("_lower"):
            assert F(compact[key]) <= value
        elif isinstance(value, F) and key.endswith("_upper"):
            assert F(compact[key]) >= value
    with pytest.raises(ValueError, match="faithful"):
        c.preparation_budget(attempts=0)


def test_one_complete_source_word_and_two_cap_Ramsey_readouts():
    s = F(16)
    bound = c.checked_bound(s)
    assert len(c.word(s)) == 2042
    assert bound["operator_error_upper"] == F(1, 10**9)
    outputs = [c.witness(i, s) for i in (0, 1)]
    assert [row["source_pulses"] for row in outputs] == [2056, 2068]
    for row in outputs:
        assert .98 < row["complete_bit_probabilities"][0] < .99
        assert abs(sum(row["complete_bit_probabilities"])-1) < 2e-11
    assert abs(outputs[0]["complete_bit_probabilities"][0]-outputs[1]["complete_bit_probabilities"][0]) < 2e-11
    rng = np.random.default_rng(447)
    probe = rng.normal(size=(512, 3))+1j*rng.normal(size=(512, 3))
    probe /= np.linalg.norm(probe)
    phase = np.array([float(w)**(-1j*16) for w in c.target_ratios()])[:, None]
    assert np.linalg.norm(c.matrix(s)@probe-phase*probe) < float(c.WORD_ERROR_GATE)
    # Full matrix certificate catches a terminal padding error, with no
    # native-shell projection or renormalization allowed.
    wrong = np.diag(phase[:, 0]); wrong[511, 511] += .125
    assert c.certify(s, wrong) >= F(1, 8)


def test_shared_clock_nontrivial_cap_action_and_trivial_intersection():
    sigma = c.reference()
    for cap in c.CAPS:
        ratios = {sigma[i]/sigma[j] for unit in c.regional.basis_supports(cap, completed=False) for i, j in unit}
        assert F(69, 68) in ratios and F(23, 20) in ratios
    for unit in c.regional.basis_supports(c.CORE, completed=False):
        assert {sigma[i]/sigma[j] for i, j in unit} == {F(1)}


@pytest.mark.parametrize("bad", ({(0, 1): F(1)}, {(0, 1): 1., (0, 2): F(1), (1, 2): F(1)},
                                   {(0, 1): F(0), (0, 2): F(1), (1, 2): F(1)}))
def test_factor_contract_rejects_missing_float_or_nonfaithful_weights(bad):
    with pytest.raises(ValueError, match="exact positive"):
        c.factor_characterization(3, bad)
