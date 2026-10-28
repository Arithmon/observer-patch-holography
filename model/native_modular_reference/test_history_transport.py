"""Independent endpoints, complex spectators, support and control failures."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from . import history_transport as h
from ..observer import native, record_recovery as rr


def independent_pairs(law, tape):
    """Floor arithmetic, not the producer's recorded_step or native.native."""
    basis = tuple(sorted(tuple(int(i in pair) for i in range(12))
                         for pair in combinations(range(12), 2)))
    positions = {state: i for i, state in enumerate(basis)}
    pairs = []
    for i, initial in enumerate(basis):
        state = list(initial)
        for a, b, coin, value in tape:
            before_a, before_b = state[a], state[b]
            total = before_a+before_b
            state[a], state[b] = (total+coin)//2, (total+1-coin)//2
            label = before_a-state[a] if law == 'C' else before_a-before_b
            if label != value:
                break
        else:
            pairs.append((i, positions[tuple(state)]))
    return tuple(pairs)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_every_seam_both_orientations_complete_source_columns(law):
    for seam in rr.seams():
        for oriented in (seam, seam[::-1]):
            branches = h.word_instrument(law, (oriented,))
            counts = [F()]*66
            for branch in branches:
                assert tuple(zip(branch.inputs, branch.outputs)) == independent_pairs(law, branch.history)
                assert branch.weight == F(1, 2)
                for i in branch.inputs:
                    counts[i] += branch.weight
            assert counts == [1]*66


@pytest.mark.parametrize('kind', ['local', 'palm', 'common', 'tracial'])
def test_source_references_are_exact_faithful_and_detached(kind):
    sigma = h.reference(kind)
    assert len(sigma) == 66 and sum(sigma) == 1 and min(sigma) > 0
    original = F(sigma[0].numerator, sigma[0].denominator)
    sigma[0]._numerator = 0
    assert h.reference(kind)[0] == original
    assert min(h.reference(kind)) > 0


@pytest.mark.parametrize('law,number,units', [('C', 14, 9524), ('E', 32, 4704)])
def test_adaptive_history_complete_support_and_phase_cocycle(law, number, units):
    branches = h.instrument(law, h.adaptive_policy(law))
    sigma, omega = h.reference(), h.reference('palm')
    report = h.phase_ratio_certificate(law, branches, sigma, omega)
    assert report['branches'] == number and report['matrix_units_checked'] == units
    blocks = h.transport(law, branches, sigma)
    assert sum(block.mass for block in blocks) == 1
    for block in blocks:
        b = block.branch
        assert tuple(zip(b.inputs, b.outputs)) == independent_pairs(law, b.history)
        assert block.support == b.outputs
        assert block.mass == sum(block.unnormalized)
        assert sum(block.normalized) == 1
        assert min(block.normalized) > 0
        for a, i in enumerate(b.inputs):
            assert block.unnormalized[a] == F(1, 4)*sigma[i]
            assert block.normalized[a] == sigma[i]/sum(sigma[j] for j in b.inputs)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_modular_conjugation_on_complex_entangled_untouched_reference(law):
    sigma = h.reference()
    blocks = h.transport(law, h.instrument(law, h.adaptive_policy(law)), sigma)
    rng = np.random.default_rng(713)
    vectors = rng.normal(size=(198, 4))+1j*rng.normal(size=(198, 4))
    rho = vectors@vectors.conj().T
    rho /= np.trace(rho)
    rho = rho.reshape(66, 3, 66, 3)
    t = .83
    phase = np.exp(1j*t*np.log(np.array(sigma, float)))
    rotated = phase[:, None, None, None]*rho*phase.conj()[None, None, :, None]
    for block in blocks:
        ix = np.array(block.branch.inputs)
        original = float(block.branch.weight)*rho[ix][:, :, ix, :]
        expected = float(block.branch.weight)*rotated[ix][:, :, ix, :]
        for weights in (block.unnormalized, block.normalized):
            p = np.exp(1j*t*np.log(np.array(weights, float)))
            actual = p[:, None, None, None]*original*p.conj()[None, None, :, None]
            assert np.max(abs(actual-expected)) < 1e-16


@pytest.mark.parametrize('law', ['C', 'E'])
def test_sequential_reference_updates_are_the_composed_actual_history(law):
    sigma = h.reference()
    first = h.word_instrument(law, ((0, 5),))
    combined = h.word_instrument(law, ((0, 5), (0, 1)))
    for block in h.transport(law, combined, sigma):
        prefix = next(b for b in first if b.history == block.branch.history[:1])
        mid = h.branch_reference(law, prefix, sigma)
        mid_weights = dict(zip(mid.support, mid.unnormalized))
        for initial, final, expected in zip(block.branch.inputs, block.support, block.unnormalized):
            middle = prefix.outputs[prefix.inputs.index(initial)]
            assert expected == mid_weights[middle]/2
            assert final == dict(independent_pairs(law, block.branch.history))[initial]


@pytest.mark.parametrize('law', ['C', 'E'])
def test_fixed_reference_and_erased_record_counterexamples(law):
    sigma = h.reference()
    branches = h.word_instrument(law, ((0, 5),))
    fixed = h.covariance_witness(law, branches, sigma)
    assert fixed['input_ratio'] == '5/8'
    assert fixed['comparison_output_ratio'] == '8/25'
    assert fixed['ratio_mismatch'] == '125/64'
    lost = h.covariance_witness(law, branches, sigma, forget_records=True)
    assert lost and F(lost['ratio_mismatch']) != 1
    assert h.covariance_witness(law, branches, h.reference('tracial'), forget_records=True) is None


def test_erasing_only_coin_already_loses_modular_transport():
    sigma = h.reference()
    basis = native.states(binary=True)
    state = lambda pair: tuple(int(i in pair) for i in range(12))
    u, v = basis.index(state((0, 1))), basis.index(state((1, 2)))
    zero = tuple(b for b in h.word_instrument('C', ((0, 5),)) if b.history[-1][-1] == 0)
    output = [F()]*66
    coefficient = F()
    for b in zero:
        for i, j in zip(b.inputs, b.outputs):
            output[j] += b.weight*sigma[i]
        if u in b.inputs and v in b.inputs:
            assert b.outputs[b.inputs.index(u)] == u and b.outputs[b.inputs.index(v)] == v
            coefficient += b.weight
    assert coefficient == F(1, 2)
    assert sigma[u]/sigma[v] == F(8, 5)
    assert output[u]/output[v] == F(4, 5)


def test_noncommuting_reference_is_not_covered_by_exact_identity():
    b = next(b for b in h.word_instrument('C', ((0, 5),)) if b.history[-1][-1] == 1)
    inside, outside = b.inputs[0], next(i for i in range(66) if i not in b.inputs)
    sigma = np.diag(np.array(h.reference(), float)).astype(complex)
    off = .4*np.sqrt(sigma[inside, inside]*sigma[outside, outside])
    sigma[inside, outside] = 1j*off
    sigma[outside, inside] = -1j*off
    eig, v = np.linalg.eigh(sigma)
    assert min(eig) > 0
    modular = (v*np.exp(.7j*np.log(eig)))@v.conj().T
    # K annihilates this input; K sigma^(it) does not. No scalar fixes it.
    assert np.linalg.norm(modular[np.array(b.inputs), outside]) > .01
    with pytest.raises(ValueError):
        h.branch_reference('C', b, sigma)


def test_sharp_hadamard_domain_counterexample_and_compact_controls():
    report=h.noncommuting_reference_control()
    u, v=report['rotated_coordinates']
    b=next(b for b in h.word_instrument('C', ((0, 5),)) if b.history == tuple(map(tuple, report['history'])))
    sigma=np.array(h.reference(), float)
    rotation=np.eye(66, dtype=complex)
    rotation[np.ix_((u, v), (u, v))]=np.array([[1, 1], [1, -1]])/np.sqrt(2)
    s=np.pi/np.log(sigma[u]/sigma[v])
    rho=rotation@np.diag(sigma)@rotation.conj().T
    flow=(rotation*np.exp(1j*s*np.log(sigma)))@rotation.conj().T
    k=np.zeros((len(b.inputs), 66),complex)
    k[np.arange(len(b.inputs)),b.inputs]=np.sqrt(float(b.weight))
    tau=k@rho@k.conj().T
    vals, vecs=np.linalg.eigh(tau)
    actual=(vecs*np.exp(1j*s*np.log(vals)))@vecs.conj().T
    defect=actual@k-np.exp(1j*s*np.log(float(b.weight)))*k@flow
    assert np.linalg.norm(defect, 2) == pytest.approx(1., abs=2e-14)
    assert np.linalg.norm(defect[:, v])**2 == pytest.approx(.5, abs=2e-14)
    assert h.lost_coin_control()['matrix_unit_entry_defect_at_s_pi_over_log2'] == '1'


def test_poisoned_upstream_reference_cache_is_not_certificate_authority(monkeypatch):
    from . import local_preparation
    sigma=list(h.reference());sigma[0],sigma[-1]=sigma[-1],sigma[0]
    monkeypatch.setattr(local_preparation, 'reference', lambda: tuple(sigma))
    with pytest.raises(ArithmeticError, match='source-selected'):
        h.reference()


def test_record_only_stopping_retains_variable_coin_cost_and_all_flags():
    first = h.word_instrument('C', ((0, 5),))
    policy = {(): (0, 5)}
    for branch in first:
        if branch.history[-1][-1] != 0:
            policy[branch.history] = None
        else:
            policy[branch.history] = (0, 1)
            for i in branch.outputs:
                for coin in (0, 1):
                    _, record = rr.recorded_step('C', native.states(True)[i], (0, 1), coin)
                    policy[(*branch.history, record)] = None
    branches = h.instrument('C', policy)
    assert set(b.weight for b in branches) == {F(1, 2), F(1, 4)}
    assert sum(b.mass for b in h.transport('C', branches, h.reference())) == 1


@pytest.mark.parametrize('mutation', ['endpoint', 'input', 'coin_weight', 'drop', 'duplicate', 'hidden_control'])
def test_forged_branch_and_incomplete_or_hidden_control_rejected(mutation):
    branches = h.word_instrument('C', ((0, 5),))
    if mutation == 'endpoint':
        b = branches[0]; changed = (replace(b, outputs=b.outputs[::-1]), *branches[1:])
    elif mutation == 'input':
        b = branches[0]; changed = (replace(b, inputs=b.inputs[:-1], outputs=b.outputs[:-1]), *branches[1:])
    elif mutation == 'coin_weight':
        changed = (replace(branches[0], weight=F(1)), *branches[1:])
    elif mutation == 'drop':
        changed = branches[:-1]
    elif mutation == 'duplicate':
        changed = (*branches, branches[0])
    else:
        # Two valid but incompatible root decisions cannot form one causal policy.
        changed = (*branches, *h.word_instrument('C', ((0, 1),)))
    with pytest.raises(ValueError):
        h.validate_instrument('C', changed)


@pytest.mark.parametrize('policy', [lambda history: (0, 5), {(): 'Hadamard'}, {(): np.eye(2)},
                                    {(): (0, 5)}, {(): None, ((0, 5, 0, 1),): None}])
def test_unimplemented_controls_and_missing_or_unreachable_decisions_rejected(policy):
    with pytest.raises(ValueError):
        h.instrument('C', policy)


@pytest.mark.parametrize('bad', [True, 1., F(0), F(-1, 66)])
def test_reference_rejects_nonexact_nonfaithful_inputs(bad):
    sigma = list(h.reference());sigma[0] = bad
    with pytest.raises(ValueError):
        h.transport('C', h.word_instrument('C', ()), tuple(sigma))


@pytest.mark.parametrize('slot,value', [('_denominator', -66), ('_denominator', 0),
                                      ('_numerator', True), ('_numerator', 2)])
def test_malformed_fraction_rejected_before_arithmetic(slot, value):
    bad = F(1, 66);setattr(bad, slot, value)
    with pytest.raises(ValueError):
        h.sandwich_bound(bad, F(1))


def test_zero_branch_and_unsupported_carrier_refused():
    with pytest.raises(ValueError):
        h.branch_reference('C', h.Branch(((0, 5, 0, -1),), (), (), F(1, 2)), h.reference())
    with pytest.raises(ValueError):
        h.transport('C', h.word_instrument('C', ()), (F(1, 78),)*78)


def test_reference_block_mutation_does_not_change_later_results():
    branches = h.word_instrument('C', ((0, 5),))
    block = h.transport('C', branches, h.reference())[0]
    block.unnormalized[0]._numerator = -1
    block.branch.weight._denominator = -2
    assert min(h.transport('C', branches, h.reference())[0].unnormalized) > 0
    with pytest.raises(ValueError):
        h.branch_reference('C', block.branch, h.reference())


def test_exact_sandwich_bounds_zero_and_near_one_without_float_rounding():
    assert h.sandwich_bound(F(0), F(9))['complete_intertwining_half_diamond_upper'] == 0
    result = h.sandwich_bound(F(1, 10), F(-3, 2))
    assert result['generator_error_up_to_scalar'] == F(1, 9)
    assert result['complete_intertwining_half_diamond_upper'] == F(1, 3)
    assert h.sandwich_bound(1-F(1, 10**80), 1)['complete_intertwining_half_diamond_upper'] == 1
    for invalid in (F(1), F(2), F(-1), .1, True):
        with pytest.raises(ValueError):h.sandwich_bound(invalid, 1)


def test_noncommuting_sandwich_survives_rare_branches_and_normalization():
    sigma = np.diag(np.array(h.reference(), float))
    rng = np.random.default_rng(1729)
    x = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66));x=(x+x.conj().T)/2
    x /= np.linalg.norm(x, 2)
    epsilon = .02;root=np.diag(np.sqrt(np.diag(sigma)))
    rho = root@(np.eye(66)+epsilon*x)@root
    def log(a):
        vals, vecs = np.linalg.eigh(a)
        return (vecs*np.log(vals))@vecs.conj().T
    for b in h.word_instrument('E', ((0, 5), (0, 1), (1, 8))):
        ix=np.array(b.inputs);a=float(b.weight)*rho[np.ix_(ix, ix)];s=float(b.weight)*sigma[np.ix_(ix, ix)]
        assert np.linalg.eigvalsh(a-(1-epsilon)*s).min() >= -1e-17
        assert np.linalg.norm(log(a)-log(s), 2) <= epsilon/(1-epsilon)+1e-13
        # Normalization changes only the generator scalar, including rare ranks1.
        difference=log(a/np.trace(a))-log(s/np.trace(s))+np.log(np.trace(a)/np.trace(s))*np.eye(len(ix))
        assert np.linalg.norm(difference, 2) <= epsilon/(1-epsilon)+1e-13


def test_compact_certificate_is_source_bound_and_not_a_waiting_or_locality_claim():
    report=h.exact_certificate()
    assert report['finite_preparation']['preparation_fine_tapes_accessible_to_later_controller'] is False
    assert F(report['finite_preparation']['complete_intertwining_half_diamond_upper']) < F(1, 1000)
    assert report['domain']['coherent_interventions_accepted'] is False
    assert report['domain']['fixed_bare_cap_locality_claimed'] is False
    assert report['domain']['autonomous_modular_waiting_claimed'] is False
    assert 'model/observer/record_recovery.py' in h.source_bindings()
