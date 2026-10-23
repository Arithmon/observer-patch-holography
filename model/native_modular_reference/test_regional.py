"""Independent native-region, CP, minimality and source-word controls."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import block_diag

from . import regional as r


def _matrix(unit):
    out = np.zeros((66, 66), complex)
    for i, j in unit:
        out[i, j] = 1
    return out


def _global_modular(array, s):
    phases = np.array([float(w)**(1j*s) for w in r.source()[2]])
    return phases[:, None]*array*phases.conj()[None, :]


def test_actual_degree_generator_and_all_region_factorizations():
    kernel, degree, sigma, _ = r.source()
    assert len(kernel.seams) == 30 and {d: degree.count(d) for d in set(degree)} == {8: 30, 10: 36}
    # Every port subset has a complete, faithful native decomposition.
    for mask in range(4096):
        ports = tuple(i for i in range(12) if mask & (1 << i))
        bs = r.blocks(ports)
        assert sum(b.probability for b in bs) == 1
        for b in bs:
            assert sum(b.local_weights) == sum(b.outside_weights) == 1
            for a, row in enumerate(b.indices):
                for j, index in enumerate(row):
                    assert sigma[index] == b.probability*b.local_weights[a]*b.outside_weights[j]


@pytest.mark.parametrize("ports", ((0, 1), (0, 2, 3), r.CAP, (0, 1, 2, 3, 4, 5)))
def test_minimal_context_partition_from_independent_all_ratio_classes(ports):
    # Independently classify external configurations by ALL matrix-unit
    # modular eigenvalues, not by the normalized-column implementation.
    kernel, _, sigma, index = r.source()
    assert len(kernel.states) == 66
    outside = tuple(p for p in range(12) if p not in ports)
    for q in range(3):
        local = tuple(combinations(ports, q))
        groups = {}
        for y in combinations(outside, 2-q):
            values = [sigma[index[tuple(sorted(x+y))]] for x in local]
            signature = tuple(a/b for a in values for b in values)
            groups.setdefault(signature, set()).add(y)
        actual = {frozenset(b.outside) for b in r.blocks(ports) if b.charge == q}
        assert actual == {frozenset(g) for g in groups.values()}
    # Every local matrix unit has one common modular ratio within its context.
    for unit in r.basis_supports(ports):
        assert len({sigma[i]/sigma[j] for i, j in unit}) == 1


def test_bare_record_clock_is_trivial_and_bare_coherent_clock_is_not_restriction():
    for port in range(12):
        assert len(r.basis_supports((port,))) == 2
        for unit in r.basis_supports((port,)):
            assert all(i == j for i, j in unit)
    # The adjacent two-port marginal is tracial in each charge block: its
    # intrinsic regional modular action is identity, but ambient action moves
    # its off-diagonal unit outside the bare algebra.
    assert r.bare_reference((0, 1)) == ((0, (F(17, 25),)),
                                     (1, (F(23, 150), F(23, 150))),
                                     (2, (F(1, 75),)))
    witness = r.noninvariance_witness((0, 1))
    assert witness["first_ratio"] == F(4, 5) and witness["second_ratio"] == 1
    unit = _matrix(witness["matrix_unit_support"])
    evolved = _global_modular(unit, 1)
    first, second = witness["first_pair"], witness["second_pair"]
    assert abs(evolved[first]-evolved[second]) > .2
    np.testing.assert_allclose(r.conditional_expectation(evolved, (0, 1)), evolved, atol=1e-15)


def test_state_preserving_expectation_unital_CP_bimodule_and_covariance():
    rng = np.random.default_rng(691)
    sigma = np.diag(list(map(float, r.source()[2])))
    a = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    e = lambda x: r.conditional_expectation(x, r.CAP)
    np.testing.assert_allclose(e(np.eye(66)), np.eye(66), atol=2e-15)
    np.testing.assert_allclose(e(e(a)), e(a), atol=2e-15)
    np.testing.assert_allclose(np.trace(sigma@e(a)), np.trace(sigma@a), atol=2e-15)
    # The Palm-preserving expectation is generally NOT tracial. It therefore
    # cannot instantiate CPRegionalNet.expect_trace just from this formula.
    first = r.blocks(r.CAP)[0]
    chosen = next(j for j, w in enumerate(first.outside_weights)
                  if w != F(1, len(first.outside)))
    projector = _matrix(((first.indices[0][chosen], first.indices[0][chosen]),))
    assert abs(np.trace(e(projector))-np.trace(projector)) > .01
    # An independent Kraus calculation on every context: sqrt(nu_z) I|z><y|.
    expected = np.zeros_like(a)
    for block in r.blocks(r.CAP):
        n, m = block.shape
        for z in range(m):
            for y in range(m):
                rows = [row[z] for row in block.indices]
                cols = [row[y] for row in block.indices]
                expected[np.ix_(cols, cols)] += float(block.outside_weights[z])*a[np.ix_(rows, rows)]
    np.testing.assert_allclose(expected, e(a), atol=2e-15)
    unit = next(u for u in r.basis_supports(r.CAP) if any(i != j for i, j in u))
    b = _matrix(unit)
    np.testing.assert_allclose(e(b@a@b.conj().T), b@e(a)@b.conj().T, atol=2e-15)
    np.testing.assert_allclose(e(_global_modular(a, 1.3)), _global_modular(e(a), 1.3), atol=3e-15)
    positive = e(a@a.conj().T)
    assert np.linalg.eigvalsh(positive).min() > 0


def test_restriction_complete_contexts_and_arbitrary_reference():
    rng = np.random.default_rng(834)
    amplitude = rng.normal(size=(66, 3))+1j*rng.normal(size=(66, 3))
    amplitude /= np.linalg.norm(amplitude)
    kraus = r.restriction_kraus(r.CAP)
    np.testing.assert_array_equal(sum(k.conj().T@k for k in kraus), np.eye(66))
    vector = amplitude.ravel()
    reference_state = amplitude.T@amplitude.conj()
    out = sum(np.outer((k@amplitude).ravel(), (k@amplitude).ravel().conj()) for k in kraus)
    dimension = kraus[0].shape[0]
    out_reference = np.trace(out.reshape(dimension, 3, dimension, 3), axis1=0, axis2=2)
    np.testing.assert_allclose(out_reference, reference_state, atol=3e-16)
    assert abs(np.trace(out)-np.vdot(vector, vector)) < 3e-16
    assert np.linalg.eigvalsh(out).min() > -2e-16
    source_state = amplitude@amplitude.conj().T
    reduced = sum(k@source_state@k.conj().T for k in kraus)
    np.testing.assert_allclose(reduced, block_diag(*r.restrict_state(source_state, r.CAP)), atol=3e-16)


def test_isotony_intersection_agreement_and_disjoint_locality_control():
    # Minimal modular completion is isotone. Independent sparse membership
    # checks cover nontrivial nested coherent algebras.
    for small, large in (((0,), (0, 1)), ((0, 1), r.CAP),
                         ((0, 5), r.CAP), (r.CAP, (0, 1, 2, 5, 7)),
                         ((2, 8), (0, 1, 2, 8))):
        assert all(r.membership(unit, large) for unit in r.basis_supports(small))
    common = (0, 1)
    for unit in r.basis_supports(common):
        assert r.membership(unit, r.CAP) and r.membership(unit, (0, 1, 2, 8))
        evolved = _global_modular(_matrix(unit), 2.1)
        np.testing.assert_allclose(r.conditional_expectation(evolved, r.CAP), evolved, atol=2e-15)
        np.testing.assert_allclose(r.conditional_expectation(evolved, (0, 1, 2, 8)), evolved, atol=2e-15)
    counter = r.locality_counterexample()
    a, b = _matrix(counter["left_matrix_unit"]), _matrix(counter["right_matrix_unit"])
    np.testing.assert_array_equal(a@b-b@a, _matrix(((64, 58),)))
    # The original disjoint local algebras do commute; refinement loses that.
    for u in r.basis_supports((0, 1), completed=False):
        for v in r.basis_supports((2, 8), completed=False):
            left, right = _matrix(u), _matrix(v)
            np.testing.assert_array_equal(left@right, right@left)


def test_nested_state_preserving_expectations_compose():
    rng = np.random.default_rng(738)
    a = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    for small, large in (((0,), (0, 1)), ((0, 1), r.CAP)):
        e_small = lambda x: r.conditional_expectation(x, small)
        e_large = lambda x: r.conditional_expectation(x, large)
        np.testing.assert_allclose(e_small(e_large(a)), e_small(a), atol=2e-15)
        np.testing.assert_allclose(e_large(e_small(a)), e_small(a), atol=2e-15)


def test_entropy_completion_and_modular_completion_are_distinct_actual_maps():
    assert len(r.basis_supports(r.CAP, completed=False)) == 53
    assert len(r.basis_supports(r.CAP)) == 149
    assert r.tracial_completion_distance(r.CAP) == F(8, 175)
    assert [sum(b.probability for b in r.blocks(r.CAP) if b.charge == q) for q in range(3)] == [F(5, 12), F(1, 2), F(1, 12)]


def test_same_source_word_realizes_all_regional_flows_retaining_padding():
    s = F(1)
    budget = r.regional_budget(r.CAP, s)
    rng = np.random.default_rng(448)
    amplitude = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    amplitude /= np.linalg.norm(amplitude)
    codes = r.full_flow.support()[0]
    encoded = np.zeros((512, 2), complex)
    encoded[list(codes)] = amplitude
    actual = r.full_flow.matrix(s)@encoded
    complete = actual@actual.conj().T
    for ports in ((), (0,), (0, 1), r.CAP, tuple(range(12))):
        actual_blocks = r.restrict_encoded_state(complete, ports)
        target_blocks = []
        for block, state in zip(r.blocks(ports), r.restrict_state(amplitude@amplitude.conj().T, ports), strict=True):
            phase = np.array([float(w)**(-1j*float(s)) for w in block.local_weights])
            target_blocks.append(phase[:, None]*state*phase.conj()[None, :])
        target_blocks.append(np.zeros((1, 1)))
        distance = sum(np.linalg.svd(a-b, compute_uv=False).sum()/2
                       for a, b in zip(actual_blocks, target_blocks, strict=True))
        assert distance < float(r.full_flow.checked_bound(s)["half_diamond_upper"])
        assert abs(sum(np.trace(a) for a in actual_blocks)-np.trace(complete)) < 2e-15
    assert budget["source_then_restrict_vs_actual_regional_modular_half_diamond_upper"] < F(61, 10**6)
    # Loss outside the source shell is an explicit output, never normalized away.
    control = np.zeros((512, 512), complex); control[511, 511] = .25
    control[codes[0], codes[0]] = .75
    reduced = r.restrict_encoded_state(control, r.CAP)
    assert reduced[-1][0, 0] == .25 and sum(np.trace(a) for a in reduced) == 1


def test_finite_preparation_sandwich_survives_regional_partial_trace():
    # Independent coherent perturbation with a known relative operator bound.
    rng = np.random.default_rng(41)
    sigma = np.diag(list(map(float, r.source()[2])))
    root = np.sqrt(sigma)
    h = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    h = h+h.conj().T; h /= np.linalg.norm(h, 2)
    t = .01
    unnormalized = root@(np.eye(66)+t*h)@root
    for actual, ideal in zip(r.restrict_state(unnormalized, r.CAP), r.restrict_state(sigma, r.CAP), strict=True):
        assert np.linalg.eigvalsh(actual-(1-t)*ideal).min() >= -2e-16
        assert np.linalg.eigvalsh((1+t)*ideal-actual).min() >= -2e-16
        assert np.linalg.eigvalsh(actual).min() > 0
    with pytest.raises(ValueError, match="faithful"):
        r.regional_budget(r.CAP, F(1), attempts=0)


@pytest.mark.parametrize("bad", ((0, 0), (12,), (True,), (1.0,), "01"))
def test_region_contract_rejects_ambiguous_port_inputs(bad):
    with pytest.raises(ValueError, match="ports"):
        r.blocks(bad)
