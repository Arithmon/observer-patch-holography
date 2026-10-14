"""Independent finite-code, active-memory and passive-erasure controls.

The endpoint oracle below does not use the implementation's branch images,
probabilities or recovery targets. Exact counting checks the established code
on all native edges; source executions check its concrete new recoveries.
"""
from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from ..observer import record_recovery
from ..observer.source_discrimination import circuit
from . import encoded_memory as memory
from . import protocol, quantum_native as joint


PAIRS = (((0, 2), (1, 3)), ((0, 3), (1, 2)))
SUPPORTS = tuple(tuple(tuple(int(port in pair) for port in range(12))
                       for pair in logical) for logical in PAIRS)


def _event(state, seam, coin, law):
    a, b = seam
    low, parity = divmod(state[a]+state[b], 2)
    final = list(state)
    final[a], final[b] = low+coin*parity, low+(1-coin)*parity
    current = state[a]-final[a]
    raw = {0: 0, 1: 1, -1: 2}[current]
    if law == 'E':
        raw += 4*parity
    return tuple(final), raw


def _images(law, seam, coin, raw):
    return tuple(tuple(circuit.encode(final) for initial in logical
                       for final, actual in [_event(initial, seam, coin, law)]
                       if actual == raw) for logical in SUPPORTS)


def _weight(law, seam, coin, raw):
    left, right = _images(law, seam, coin, raw)
    assert len(left) == len(right)
    assert len(set(left)) == len(left) and len(set(right)) == len(right)
    assert set(left).isdisjoint(right)
    return F(len(left), 4)


def _code():
    result = np.zeros((512, 2), complex)
    for logical, support in enumerate(SUPPORTS):
        for state in support:
            result[circuit.encode(state), logical] = 1/math.sqrt(2)
    return result


@pytest.mark.parametrize('law', ['C', 'E'])
def test_original_endpoint_effect_is_exactly_scalar_on_the_code_for_all_native_edges(law):
    edges = record_recovery.seams()
    assert len(edges) == 30
    for seam in edges:
        weights = []
        for coin, raw in product((0, 1), range(4 if law == 'C' else 8)):
            # This checks the exact Gram matrix: the common diagonal is k/4,
            # and disjoint support makes the cross-logical entry exactly zero.
            weight = _weight(law, seam, coin, raw)
            assert weight in (F(0), F(1, 4), F(1, 2))
            weights.append(weight)
        assert sum(weights, F(0)) == 1


@pytest.mark.parametrize('law', ['C', 'E'])
def test_randomized_record_adaptive_seam_policy_retains_exact_branch_weights(law):
    # A two-event corrected policy. The second random choice depends on the
    # first retained native coin, so it cannot be replaced by one global
    # seam probability. Selection is classical, independent of the qubit.
    # Each corrected branch multiplies every logical matrix unit by its
    # native probability times the explicitly retained policy probability.
    histories = [((), F(1), F(1))]
    for stage in range(2):
        extended = []
        for prefix, complete_weight, native_weight in histories:
            first_probability = (F(1, 3) if stage == 0 or prefix[-1][1] == 0
                                 else F(2, 3))
            policy = (((0, 1), first_probability), ((2, 3), 1-first_probability))
            assert sum((p for _, p in policy), F(0)) == 1
            for seam, probability in policy:
                for coin, raw in product((0, 1), range(4 if law == 'C' else 8)):
                    native = _weight(law, seam, coin, raw)
                    if native:
                        extended.append((prefix+((seam, coin, raw),),
                                         complete_weight*probability*native,
                                         native_weight*native))
        histories = extended
    assert len(histories) == 64
    assert sum((weight for _, weight, _ in histories), F(0)) == 1
    # Dropping the sqrt(policy probability) from either branch operator
    # incorrectly counts the two seam choices as two complete instruments.
    assert sum((weight for _, _, weight in histories), F(0)) == 4
    assert {weight/native for _, weight, native in histories} == {
        F(1, 9), F(2, 9), F(4, 9)}
    # This exact coefficient applies also to off-diagonal logical/reference
    # matrix units, rather than only the diagonal input probabilities.
    for entry in (F(1), F(-1), F(2, 7)):
        assert sum((weight*entry for _, weight, _ in histories), F(0)) == entry


def test_source_encoder_and_inverse_preserve_unknown_complex_reference_columns():
    rng = np.random.default_rng(1101001)
    logical = rng.normal(size=(2, 3))+1j*rng.normal(size=(2, 3))
    logical /= np.linalg.norm(logical)
    actual = memory.encoding(logical)
    np.testing.assert_allclose(actual, _code()@logical, atol=4e-11, rtol=0)
    decoded = memory.decoding(actual)
    target = np.zeros((512, 3), complex)
    target[:2] = logical
    assert decoded.shape == target.shape  # No syndrome/padding compression.
    np.testing.assert_allclose(decoded, target, atol=8e-11, rtol=0)
    np.testing.assert_allclose(memory.code_matrix(), _code(), atol=4e-11, rtol=0)


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('seam', [(0, 1), (0, 5), (2, 3)])
def test_every_executed_corrected_branch_has_the_logical_identity_intertwiner(law, seam):
    branches = memory.repair_branches(law, seam, memory.code_matrix())
    assert [(b['coin'], b['record_code']) for b in branches] == list(
        product((0, 1), range(4 if law == 'C' else 8)))
    gram = np.zeros((2, 2), complex)
    for branch in branches:
        coin, raw = branch['coin'], branch['record_code']
        weight = _weight(law, seam, coin, raw)
        assert branch['code_branch_weight'] == weight
        assert memory.branch_weight(law, seam, coin, raw) == weight
        actual = branch['amplitude']
        assert actual.shape == (512, 2)
        np.testing.assert_allclose(actual, math.sqrt(float(weight))*_code(),
                                   atol=8e-11, rtol=0)
        gram += actual.conj().T@actual
    np.testing.assert_allclose(gram, np.eye(2), atol=2e-10, rtol=0)


def _passive_histories(law, word):
    # Exact amplitudes have one common factor for this fixed word. Integer
    # supports suffice to test whether any history preserves logical phase.
    histories = defaultdict(lambda: [[], []])
    for logical, support in enumerate(SUPPORTS):
        for original in support:
            for coins in product((0, 1), repeat=len(word)):
                state, record = original, ()
                for seam, coin in zip(word, coins, strict=True):
                    state, raw = _event(state, seam, coin, law)
                    record += ((coin, raw),)
                histories[record][logical].append(circuit.encode(state))
    return histories


@pytest.mark.parametrize('law', ['C', 'E'])
def test_passive_crossed_seams_erase_phase_but_the_old_seams_do_not(law):
    crossed = _passive_histories(law, ((0, 1), (2, 3)))
    assert len(crossed) == 16
    # Every history has only one original code support point, belonging to
    # one logical basis state. Thus its logical effect is diag(1/8,0) or
    # diag(0,1/8), and |+> / |-> differ only by an unobservable branch phase.
    counts = Counter()
    for left, right in crossed.values():
        assert (len(left), len(right)) in ((1, 0), (0, 1))
        counts[bool(right)] += 1
        plus = {address: F(1, 4) for address in (*left, *right)}
        minus = {address: F(1, 4) if left else F(-1, 4)
                 for address in (*left, *right)}
        assert {(x, y): a*b for x, a in plus.items() for y, b in plus.items()} == {
            (x, y): a*b for x, a in minus.items() for y, b in minus.items()}
    assert counts == {False: 8, True: 8}
    old = _passive_histories(law, ((0, 1), (0, 5)))
    assert all(len(left) == len(right) for left, right in old.values())
    assert all(len(set(left)) == len(left) and len(set(right)) == len(right)
               and set(left).isdisjoint(right) for left, right in old.values())


@pytest.mark.parametrize('law', ['C', 'E'])
def test_actual_two_event_active_memory_preserves_bell_reference_and_every_raw_branch(law):
    # Columns are the untouched reference; this is a complex Bell pair.
    logical = np.diag([1., 1j])/math.sqrt(2)
    encoded = memory.encoding(logical)
    histories = []
    first = memory.repair_branches(law, (0, 1), encoded)
    for a in first:
        if a['abort']:
            histories.append((a['amplitude'], F(0)))
            continue
        for b in memory.repair_branches(law, (2, 3), a['amplitude']):
            weight = _weight(law, (0, 1), a['coin'], a['record_code'])*_weight(
                law, (2, 3), b['coin'], b['record_code'])
            histories.append((b['amplitude'], weight))
    assert len(histories) == (50 if law == 'C' else 166)
    assert sum(weight > 0 for _, weight in histories) == 16
    assert sum((weight for _, weight in histories), F(0)) == 1
    reference = np.zeros((2, 2), complex)
    fidelity = 0.
    ideal = _code()@logical
    for actual, weight in histories:
        np.testing.assert_allclose(actual, math.sqrt(float(weight))*ideal,
                                   atol=2e-10, rtol=0)
        reference += actual.conj().T@actual
        fidelity += abs(np.vdot(ideal, actual))**2
    np.testing.assert_allclose(reference, np.eye(2)/2, atol=3e-10, rtol=0)
    assert abs(fidelity-1) < 3e-10


@pytest.fixture(scope='module')
def work_histories():
    return {current: protocol.execute(protocol.continuation_program(
                2, 8, native_current=current), backend='product')
            for current in (-1, 0, 1)}


@pytest.fixture(scope='module')
def protected_history():
    return joint.execute('C', 2, 8, native_mode='protected_qubit')


def test_integrated_source_branch_intertwiner_on_entangled_logical_work_and_reference(
        protected_history, work_histories):
    rng = np.random.default_rng(1101002)
    initial = rng.normal(size=(2, 16, 2))+1j*rng.normal(size=(2, 16, 2))
    initial /= np.linalg.norm(initial)
    targets = {j: {b.record: b.kraus for b in history.branches}
               for j, history in work_histories.items()}
    errors, mass = [], []
    expected_marginal = np.einsum('lwr,mws->lrms', initial, initial.conj())
    marginal = np.zeros((2, 2, 2, 2), complex)
    coded_marginal = np.zeros((2, 2, 2, 2), complex)
    outside_errors = []
    for record, actual in joint.output_on_joint_input(protected_history, initial):
        c1, raw1, bits, c2, raw2 = record
        if bits:
            current = {0: 0, 1: 1, 2: -1}[raw1 % 4]
            seam = (2, 3) if bits[2] == 0 else (0, 1)
            weight = _weight('C', (0, 1), c1, raw1)*_weight('C', seam, c2, raw2)
            work = np.einsum('ab,lbr->lar', targets[current][bits], initial)
            expected = math.sqrt(float(weight))*np.einsum('il,lwr->iwr', _code(), work)
            marginal += float(weight)*np.einsum('lwr,mws->lrms', work, work.conj())
        else:
            expected = np.zeros_like(actual)
        errors.append(float(np.linalg.norm(actual-expected)**2))
        mass.append(float(np.vdot(actual, actual).real))
        # These are bounded encoded observable tests, not a replacement of
        # the complete output by its compression. Separately retain the full
        # off-code amplitude error, with no renormalization or postselection.
        coded = np.einsum('il,iwr->lwr', _code().conj(), actual)
        coded_marginal += np.einsum('lwr,mws->lrms', coded, coded.conj())
        outside_errors.append(float(np.linalg.norm(actual-np.einsum('il,lwr->iwr', _code(), coded))**2))
    np.testing.assert_allclose(marginal, expected_marginal, atol=3e-11, rtol=0)
    np.testing.assert_allclose(coded_marginal, expected_marginal, atol=3e-10, rtol=0)
    assert math.sqrt(math.fsum(outside_errors)) < 4e-10
    assert math.sqrt(math.fsum(errors)) < 4e-10
    assert abs(math.fsum(mass)-1) < 3e-10


def _phase_distance(history):
    plus, minus = np.array([1., 1.])/math.sqrt(2), np.array([1., -1.])/math.sqrt(2)
    pieces = []
    for branch in history.branches:
        a, b = branch.native_map@plus, branch.native_map@minus
        # Reduce uu* - vv* by QR before diagonalizing. The closed-form square
        # root subtracts almost equal squares on erased branches and can turn
        # ordinary rounding into a fabricated O(sqrt(epsilon)) distance.
        _, triangle = np.linalg.qr(np.column_stack((a, b)), mode='reduced')
        reduced = triangle@np.diag([1., -1.])@triangle.conj().T
        distance = float(np.sum(abs(np.linalg.eigvalsh(reduced))))/2
        pieces.append(distance*float(np.vdot(branch.work_state, branch.work_state).real))
    return math.fsum(pieces)


def test_logical_memory_survives_a_deliberately_different_complete_work_channel(protected_history):
    # A measure-and-erase instrument K_r=|0><r| replaces the Hamiltonian work
    # program. It is deliberately inaccurate as that program, but is CPTP.
    # This is an independent matrix-instrument control, not a claim that a
    # new reset operation was added to the source execution grammar.
    branches = []
    for branch in protected_history.branches:
        bits = branch.record[2]
        matrix = np.eye(16, dtype=complex) if not bits else np.zeros((16, 16), complex)
        if bits:
            matrix[0, sum(bit << index for index, bit in enumerate(bits))] = 1
        branches.append(replace(branch, work_map=matrix,
                                work_state=matrix[:, :1]))
    altered = replace(protected_history, branches=tuple(branches))
    rng = np.random.default_rng(1101003)
    initial = rng.normal(size=(2, 16, 2))+1j*rng.normal(size=(2, 16, 2))
    initial /= np.linalg.norm(initial)
    target = np.einsum('lwr,mws->lrms', initial, initial.conj())
    marginal = np.zeros_like(target)
    outside, work_outside_zero = [], []
    for _, actual in joint.output_on_joint_input(altered, initial):
        coded = np.einsum('il,iwr->lwr', _code().conj(), actual)
        marginal += np.einsum('lwr,mws->lrms', coded, coded.conj())
        outside.append(float(np.linalg.norm(actual-np.einsum('il,lwr->iwr', _code(), coded))**2))
        work_outside_zero.append(float(np.linalg.norm(actual[:, 1:])**2))
    np.testing.assert_allclose(marginal, target, atol=3e-10, rtol=0)
    assert math.sqrt(math.fsum(outside)) < 4e-10
    assert math.fsum(work_outside_zero) < 1e-18
    assert np.linalg.norm(initial[:, 1:]) > .5  # The work input really changed.


@pytest.mark.parametrize('law', ['C', 'E'])
def test_integrated_passive_phase_distance_and_active_unit_contrast(law, work_histories):
    q = {j: math.fsum(b.probability for b in h.branches if b.record[2] == 0)
         for j, h in work_histories.items()}
    qbar = (q[-1]+2*q[0]+q[1])/4
    active = joint.execute(law, 2, 8, native_mode='protected_qubit')
    passive = joint.execute(law, 2, 8, native_mode='passive_qubit')
    assert abs(_phase_distance(active)-1) < 3e-10
    assert abs(_phase_distance(passive)-(1-qbar)) < 3e-10
    assert 0 < qbar/2 < .5  # Strict generic-decoder obstruction for this run.
    # The executed decoder and source basis read, rather than a substituted
    # code projector, read out the protected phase.
    if law == 'C':
        p = joint.readout(active, np.array([1., 1.])/math.sqrt(2))
        m = joint.readout(active, np.array([1., -1.])/math.sqrt(2))
        assert abs(p['plus_code_zero_probability']-1) < 3e-10
        assert m['plus_code_zero_probability'] < 3e-10
        assert abs(p['probability_mass']-1) < 3e-10
        assert abs(m['probability_mass']-1) < 3e-10


def test_uniform_passive_obstruction_constants_have_independent_rational_certificates():
    from . import encoded_observer
    # d2: |sin(lambda)| <= |lambda| < 3/4 for lambda^2=1/2.
    assert F(1, 2) < F(3, 4)**2
    assert (1-F(3, 4))/2 == F(1, 8)
    # d4: positive roots are sqrt((3 +/- sqrt(6))/2). Bounds on sqrt(6)
    # establish their location; 3<pi<22/7 is the classical rational pi bound.
    assert F(12, 5)**2 < 6 < 3**2 < 5**2
    assert F(27, 10) == (3+F(12, 5))/2 > F(41, 25)**2
    assert F(41, 25)-F(11, 7) == F(12, 175) > F(1, 16)
    assert 1**2 < 6  # lambda_minus<1, hence pi/2-lambda_minus>1/2.
    x = F(1, 16)
    assert 1-x*x/2+x**4/24 <= 1-x*x/3
    assert x*x/6 == F(1, 1536)
    # Derive the characteristic polynomial in exact Q[sqrt(2),sqrt(3)]
    # through the Jacobi determinant recurrence, whose squared links are
    # 1/2,1,3/2. Entries are polynomial coefficients in ascending order.
    polys = [[F(1)], [F(0), F(1)]]
    for square in (F(1, 2), F(1), F(3, 2)):
        lifted = [F(0), *polys[-1]]
        for k, coefficient in enumerate(polys[-2]):
            lifted[k] -= square*coefficient
        polys.append(lifted)
    assert polys[-1] == [F(3, 4), F(0), F(-3), F(0), F(1)]
    for dimension, effect, obstruction in ((2, F(1, 8), F(1, 16)),
                                            (4, F(1, 1536), F(1, 3072))):
        report = encoded_observer.passive_obstruction(dimension)
        assert F(report['uniform_second_smooth_plus_lower_exact']) == effect
        assert F(report['uniform_passive_worst_case_decoder_error_lower_exact']) == obstruction == effect/2
