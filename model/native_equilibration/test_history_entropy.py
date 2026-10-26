"""Independent source, stopping-cost, charge-ambiguity and consensus controls."""
from collections import defaultdict
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from . import history_entropy as h, common_variation as v
from ..native_modular_reference import local_preparation as local
from ..native_modular_reference import common_caps as common
from ..observer import native


@pytest.mark.parametrize("law", ["C", "E"])
def test_complete_native_history_reversal_current_balance_and_record_decoder(law):
    source = native.states(binary=True)
    paths = h.histories(law)
    totals = [F() for _ in source]
    flag = [F() for _ in source]
    crossings = 0
    feedback_seams = set()
    for path in paths:
        totals[path.initial] += path.chance
        reconstructed = list(path.endpoint)
        for record in reversed(path.tape):
            a, b = record[:2]
            transfer = h.current(record, law)
            reconstructed[a] += transfer
            reconstructed[b] -= transfer
        assert tuple(reconstructed) == source[path.initial]
        # Independent cap current balance on the same trajectory.
        for cap in ((0, 1), (0, 2, 8), (2, 3, 4, 6)):
            delta = sum(path.endpoint[i]-source[path.initial][i] for i in cap)
            flux = sum(h.current(r, law)*((r[1] in cap)-(r[0] in cap)) for r in path.tape)
            assert delta == flux
        crossings += h.current(path.tape[0], law) != 0
        feedback_seams.add(path.tape[1][:2])
        if path.terminal == "success":
            flag[path.initial] += path.chance
            assert h.decode_core(path.tape, 2, law) == source[path.initial][:2]
        else:
            with pytest.raises(ValueError):
                h.decode_core(path.tape, 2, law)
    assert totals == [F(1)]*66
    assert flag == [F(9, 16) if law == "C" else F(1)]*66
    assert crossings > 0 and feedback_seams == {(0, 1), (1, 8)}


def test_actual_record_first_law_exact_and_negative_fisher_curvature():
    sigma = local.reference()
    paths = tuple(p for p in h.histories() if p.terminal == "success")
    reference = h.ledger(paths, sigma)
    assert reference["functional"] == h.source_entropy(sigma)
    # A general source perturbation checks both moment multipliers, rather
    # than merely a tangent engineered to leave the observable unchanged.
    p = tuple(F(i+1, sum(range(1, 67))) for i in range(66))
    moved = h.ledger(paths, p)
    m, c = v.moments(sigma)
    dm, dc = (a-b for a, b in zip(v.moments(p), (m, c), strict=True))
    cert = v.family_certificate(m, c)
    linear = defaultdict(F)
    v._log_add(linear, dm, cert["linear_multiplier_log_argument"])
    v._log_add(linear, dc, cert["pair_multiplier_log_argument"])
    relative = h.relative_entropy(moved["probabilities"], reference["probabilities"])
    assert h._sum_forms(moved["functional"], relative) == h._sum_forms(reference["functional"], linear)
    # A within-q perturbation that is seen by the actual current readout.
    basis = native.states(binary=True)
    seen = basis.index(tuple(int(i in (2, 5)) for i in range(12)))
    unseen = basis.index(tuple(int(i in (2, 3)) for i in range(12)))
    direction = [F()]*66
    direction[seen], direction[unseen] = F(1), F(-1)
    step = F(1, 100000)
    plus = h.ledger(paths, tuple(w+step*x for w, x in zip(sigma, direction, strict=True)))
    minus = h.ledger(paths, tuple(w-step*x for w, x in zip(sigma, direction, strict=True)))
    delta = {key: (plus["probabilities"][key]-minus["probabilities"][key])/(2*step)
             for key in reference["probabilities"]}
    fisher = sum(float(d*d/reference["probabilities"][key]) for key, d in delta.items())
    assert fisher > 0
    numerical = (h.value(plus["functional"])+h.value(minus["functional"])
                 -2*h.value(reference["functional"]))/float(step*step)
    assert abs(numerical+fisher) < 2e-3
    for row in (plus, minus):
        deficit = h.relative_entropy(row["probabilities"], reference["probabilities"])
        assert h._sum_forms(row["functional"], deficit) == reference["functional"]


def test_native_timeouts_are_retained_and_have_nonzero_unresolved_charge_excess():
    weights = local.reference()
    paths = h.histories()
    rows = h.fibres(paths)
    unresolved = [row for row in rows.values() if len({v.charges()[p.initial] for p in row}) > 1]
    assert unresolved and all(row[0].terminal == "timeout" for row in unresolved)
    full = h.ledger(paths, weights)
    excess = h.value(full["functional"])-h.value(h.source_entropy(weights))
    assert 0 < excess <= float(F(7, 16))*math.log(5)
    # No hidden full-state read is allowed to masquerade as a native record.
    with pytest.raises(ValueError, match="charge-resolved"):
        h.conditional_cap(unresolved[0], (0, 1))
    with pytest.raises(ValueError, match="state-independent"):
        h.flag_probability(tuple(p for p in paths if p.initial != 0))


def test_same_history_cap_predictions_agree_on_nested_and_overlapping_regions():
    paths = tuple(p for p in h.histories() if p.terminal == "success")
    caps = ((0, 1, 2, 5), (0, 1, 5, 8))
    overlap = (0, 1, 5)
    nontrivial = 0
    for row in h.fibres(paths).values():
        common = h.conditional_cap(row, overlap)
        nontrivial += len(row) > 1
        for cap in caps:
            bigger = h.conditional_cap(row, cap)
            reduced = defaultdict(F)
            for bits, probability in bigger.items():
                reduced[tuple(bits[cap.index(i)] for i in overlap)] += probability
            assert dict(reduced) == common
    assert nontrivial
    target = h.source_entropy(local.reference())
    for cap in ((), *caps, overlap, tuple(range(12))):
        assert h.ledger(paths, local.reference(), cap)["functional"] == target


def test_stopping_cost_is_retained_and_one_forgotten_record_can_mix_charge():
    paths = tuple(p for p in h.histories() if p.terminal == "success")
    ledger = h.ledger(paths, local.reference())
    # Census stopping lengths genuinely vary. Ignoring the retained coin
    # entropy or the success normalization changes the equality.
    assert {len(p.tape) for p in paths} == {6, 7, 8}
    assert ledger["mean_native_attempts"] == F(20, 3)
    assert h._sum_forms(ledger["bulk"], ledger["edge"]) != h.source_entropy(local.reference())
    forgotten = defaultdict(set)
    for path in paths:
        # Keeping only the fair-coin word preserves no native load response.
        forgotten[tuple(r[2] for r in path.tape)].add(v.charges()[path.initial])
    assert all(charges == {0, 1, 2} for charges in forgotten.values())


def test_record_dependent_stopping_uses_the_perturbed_inputs_own_cost():
    paths = h.histories("E", prefix=3, stop_on_current=True)
    assert h.flag_probability(paths) == 1
    for path in paths:
        count = h.prefix_length(path.tape, "E", 3, True)
        assert h.decode_core(path.tape, count, "E") == native.states(binary=True)[path.initial][:2]
    sigma = local.reference()
    reference = h.ledger(paths, sigma)
    assert reference["functional"] == h.source_entropy(sigma)
    per_input = [sum((p.chance*len(p.tape) for p in paths if p.initial == i), F())
                 for i in range(66)]
    # Keep q and pair moments fixed, but change the stopping distribution.
    q0 = [i for i, q in enumerate(v.charges()) if q == 0]
    low, high = min(q0, key=per_input.__getitem__), max(q0, key=per_input.__getitem__)
    assert per_input[low] < per_input[high]
    moved_weights = list(sigma)
    moved_weights[low] += F(1, 1000)
    moved_weights[high] -= F(1, 1000)
    moved = h.ledger(paths, tuple(moved_weights))
    assert v.moments(tuple(moved_weights)) == v.moments(sigma)
    assert moved["mean_native_attempts"] != reference["mean_native_attempts"]
    divergence = h.relative_entropy(moved["probabilities"], reference["probabilities"])
    assert h._sum_forms(moved["functional"], divergence) == reference["functional"]
    wrong = h._sum_forms(moved["bulk"], moved["edge"], reference["coin_cost"], divergence)
    assert wrong != reference["functional"]


def test_finite_sandwich_bounds_record_tangents_and_unread_directions_stay_flat():
    paths = h.histories("E", prefix=3, stop_on_current=True)
    sigma = local.reference()
    reference = h.ledger(paths, sigma)
    radius = F(1, 10)
    moved = list(sigma)
    basis = native.states(binary=True)
    a = basis.index(tuple(int(i in (2, 5)) for i in range(12)))
    b = basis.index(tuple(int(i in (2, 3)) for i in range(12)))
    moved[a] += radius*sigma[a]
    moved[b] -= radius*sigma[b]
    assert v.moments(tuple(moved)) == v.moments(sigma)
    observed = h.ledger(paths, tuple(moved))
    variation = {key: (observed["probabilities"][key]-reference["probabilities"][key])/radius
                 for key in reference["probabilities"]}
    derivative = -sum(float(delta)*math.log(float(observed["probabilities"][key]
                                                 /reference["probabilities"][key]))
                      for key, delta in variation.items())
    norm = 2*float(sigma[a])
    assert abs(derivative) <= 2*float(radius/(1-radius))*norm
    # Two untouched spectator states are distinguishable by another readout,
    # but are a real null direction of this protocol's Fisher response.
    hidden = list(sigma)
    a = basis.index(tuple(int(i in (9, 10)) for i in range(12)))
    b = basis.index(tuple(int(i in (10, 11)) for i in range(12)))
    hidden[a] += F(1, 1000)
    hidden[b] -= F(1, 1000)
    assert h.ledger(paths, tuple(hidden))["probabilities"] == reference["probabilities"]


@pytest.mark.parametrize("deadline", [2, 3, 4])
def test_deadline_bound_is_a_complete_state_independent_native_instrument(deadline):
    paths = h.histories(prefix=0, deadline=deadline)
    successful = tuple(p for p in paths if p.terminal == "success")
    alpha = 1-F(1, 2**(deadline-1))
    assert h.flag_probability(paths) == 1
    assert h.flag_probability(successful) == alpha**2
    assert max(map(lambda p: len(p.tape), paths)) == 2*deadline
    assert h.ledger(successful, local.reference())["functional"] == h.source_entropy(local.reference())


def test_c_and_e_readout_laws_both_realize_the_same_entropy_result():
    report = h.run()
    for row in report["laws"].values():
        assert row["exact_source_entropy"] == row["exact_success_functional"]
        assert row["successful_completion_dimensions"][-1] > 1
        assert row["same_history_cap_refinements_checked"] == 6
    assert report["bounded_C_readout_service"]["native_attempts_upper"] == 50
    assert report["laws"]["E"]["complete_counting_excess"] == {}
    assert report["unconditional_E_readout_service"]["native_attempts"] == 4


@pytest.mark.parametrize("weights", [common.reference(), v.family(F(3, 5), F(1, 10))])
def test_envelope_uses_the_full_faithful_moment_family(weights):
    paths = h.histories("E", prefix=3, stop_on_current=True)
    assert h.ledger(paths, weights)["functional"] == h.source_entropy(weights)


def _coin_word_oracle(law):
    """Independent complete coin-word coupling, using integer endpoints only."""
    horizon = 8 if law == "C" else 4
    result = []
    for word in product((0, 1), repeat=horizon):
        outputs = []
        for original in native.states(binary=True):
            state, tape, offset, transfer = list(original), [], 0, 0

            def step(edge):
                nonlocal offset, transfer
                a, b = edge
                coin = word[offset]
                total = state[a]+state[b]
                left = (total+coin)//2
                transfer = state[a]-left
                label = transfer if law == "C" else state[a]-state[b]
                state[a], state[b] = left, total-left
                tape.append((*edge, coin, label))
                offset += 1

            step((0, 5))
            step((0, 1) if transfer else (1, 8))
            flag = "success"
            for edge in ((0, 1), (0, 5)):
                first = word[offset]
                step(edge)
                if law == "C":
                    for _ in range(2):
                        flip = word[offset] != first
                        step(edge)
                        if flip:
                            break
                    else:
                        flag = "timeout"
                        break
            outputs.append((tuple(tape), flag))
        result.append(outputs)
    return result


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("success_only", [False, True])
def test_every_coherent_matrix_unit_matches_independent_shared_coin_source(law, success_only):
    paths = h.histories(law)
    if success_only:
        paths = tuple(path for path in paths if path.terminal == "success")
    response = h.quantum_response(paths, law=law)
    oracle = _coin_word_oracle(law)
    for i in range(66):
        for j in range(66):
            matches = sum(row[i] == row[j] and (not success_only or row[i][1] == "success")
                          for row in oracle)
            assert response.overlap[i][j] == F(matches, len(oracle))/response.flag_probability
    assert sorted(map(len, response.recovery_classes)) == [1]*6+[8]*4+[28]
    occupation = defaultdict(set)
    for i, state in enumerate(native.states(binary=True)):
        occupation[tuple(state[port] for port in (0, 1, 5, 8))].add(i)
    assert {frozenset(row) for row in response.recovery_classes} == {frozenset(row) for row in occupation.values()}
    assert sum(len(row)**2 for row in response.recovery_classes) == 1046
    expected_zero = (2814 if success_only else 558) if law == "C" else 3278
    assert sum(x == 0 for row in response.overlap for x in row) == expected_zero
    # Equality classes are the actual primitive instrument's membership
    # classes, not eigensolver ranks or separately stipulated memory codes.
    for group in response.recovery_classes:
        assert all(response.overlap[i][j] == 1 for i in group for j in group)


def _positive_entropy(matrix):
    values = np.linalg.eigvalsh(matrix)
    assert min(values) > 0  # All test inputs and supported blocks are faithful.
    return -float(values@np.log(values))


def _quantum_functional(response, rho):
    total, cost = 0., 0.
    for branch in response.branches:
        block = float(branch.weight)*rho[np.ix_(branch.inputs, branch.inputs)]
        total += _positive_entropy(block)
        cost += float(np.trace(block).real)*len(branch.tape)*math.log(2)
    return total-cost-math.log(float(response.flag_probability))


@pytest.mark.parametrize("law", ["C", "E"])
def test_actual_Petz_decoder_preserves_recoverable_entanglement_and_every_diagonal(law):
    response = h.quantum_response(h.histories(law), law=law)
    rng = np.random.default_rng(910662)
    raw = rng.normal(size=(132, 5))+1j*rng.normal(size=(132, 5))
    rho = raw@raw.conj().T
    rho /= np.trace(rho)
    sigma = np.diag(list(map(float, local.reference())))
    tau = np.array([[.6, .1+.07j], [.1-.07j, .4]])
    reference = np.kron(sigma, tau)
    sqrt_reference = np.kron(np.diag(np.sqrt(np.diag(sigma))),
                             np.linalg.eigh(tau)[1]@np.diag(np.sqrt(np.linalg.eigvalsh(tau)))
                             @np.linalg.eigh(tau)[1].conj().T)
    pulled = np.zeros_like(rho)
    for branch in response.branches:
        # Build each actual endpoint partial injection in endpoint order.
        # Evaluate the defining Petz expression; G is used only afterward.
        ordered = sorted(zip(branch.outputs, branch.inputs, strict=True))
        kraus = np.zeros((len(ordered), 66), complex)
        for row, (_, initial) in enumerate(ordered):
            kraus[row, initial] = math.sqrt(float(branch.weight))
        kraus = np.kron(kraus, np.eye(2))
        target = kraus@reference@kraus.conj().T
        eigenvalues, eigenvectors = np.linalg.eigh(target)
        inverse = (eigenvectors/np.sqrt(eigenvalues))@eigenvectors.conj().T
        decoder = sqrt_reference@kraus.conj().T@inverse
        pulled += decoder@(kraus@rho@kraus.conj().T)@decoder.conj().T
    gram = np.asarray(response.overlap, float)
    expected = rho*np.repeat(np.repeat(gram, 2, axis=0), 2, axis=1)
    np.testing.assert_allclose(pulled, expected, atol=2e-15, rtol=0)
    np.testing.assert_allclose(np.diag(pulled), np.diag(rho), atol=2e-15, rtol=0)
    # A genuinely entangled state inside the 28-dimensional spectator block
    # is recovered with its reference, without selecting a successful record.
    a, b = max(response.recovery_classes, key=len)[:2]
    vector = np.zeros(132, complex)
    vector[2*a], vector[2*b+1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    bell = np.outer(vector, vector.conj())
    np.testing.assert_array_equal(bell*np.repeat(np.repeat(gram, 2, axis=0), 2, axis=1), bell)


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("weights", [local.reference(), v.family(F(3, 5), F(1, 10))])
def test_coherent_stationarity_hessian_and_information_loss_are_different(law, weights):
    response = h.quantum_response(h.histories(law), law=law)
    sigma = np.diag(list(map(float, weights))).astype(complex)
    # Exercise a fully invisible, partially visible, and lossless coherence.
    pairs = [next((i, j) for i in range(66) for j in range(i+1, 66)
                  if response.overlap[i][j] == target) for target in (0, F(1, 2), 1)]
    epsilon = min(map(float, weights))/8
    for index, (a, b) in enumerate(pairs):
        x = np.zeros((66, 66), complex)
        x[a, b], x[b, a] = .6+.8j, .6-.8j
        report = h.quantum_stationarity(response, weights, x)
        assert report["fixed_moment_first_derivative"] == 0
        assert report["exact_structural_flat_direction"] == (index == 0)
        assert report["exact_structural_lossless_direction"] == (index == 2)
        assert report["negative_hessian_lower"] <= -report["output_entropy_hessian"] <= report["negative_hessian_upper"]
        assert report["information_loss_hessian"] >= 0
        f0 = _quantum_functional(response, sigma)
        fp = _quantum_functional(response, sigma+epsilon*x)
        fm = _quantum_functional(response, sigma-epsilon*x)
        assert abs(f0-_positive_entropy(sigma)) < 3e-13
        observed = (fp+fm-2*f0)/epsilon**2
        assert abs(observed-report["output_entropy_hessian"]) < .01*max(1., -observed)
        if index == 0:
            assert abs(fp-f0) < 3e-13
            assert _positive_entropy(sigma+epsilon*x) < _positive_entropy(sigma)
        else:
            assert fp < f0
        # Equality in source->record relative entropy concerns G=1, whereas
        # equality in the maximum-output-entropy envelope concerns G=0.
        source_relative = _positive_entropy(sigma)-_positive_entropy(sigma+epsilon*x)
        output_relative = f0-fp
        assert source_relative+3e-13 >= output_relative >= -3e-13
        if index == 2:
            assert abs(source_relative-output_relative) < 3e-13


@pytest.mark.parametrize("law", ["C", "E"])
def test_common_reference_cap_intersections_match_independent_completion_matrices(law):
    from . import cut_completion as cut
    response = h.quantum_response(h.histories(law), law=law)
    gram = np.asarray(response.overlap, float).ravel()
    for cap in ((0, 1), (0, 1, 2, 5), (0, 1, 5, 8)):
        completion, scale = cut.completion_counts(cap)
        # Every E_R range basis operator has a disjoint set of matrix-unit
        # rows. Group equal sparse-column signatures to count its range
        # and intersections independently of quantum_region_response.
        columns = completion.tocsc()
        supports = {tuple(columns.indices[columns.indptr[i]:columns.indptr[i+1]])
                    for i in range(4356) if columns.indptr[i] != columns.indptr[i+1]}
        expected = h.quantum_region_response(response, cap)
        assert expected["cap_algebra_dimension"] == len(supports)
        assert expected["record_output_image_dimension"] == sum(any(gram[list(row)] > 0) for row in supports)
        assert expected["relative_entropy_equality_algebra_dimension"] == sum(all(gram[list(row)] == 1) for row in supports)
        sigma = np.diag(list(map(float, local.reference()))).ravel()
        np.testing.assert_allclose(completion@sigma/scale, sigma, atol=2e-17, rtol=0)
    core = h.quantum_region_response(response, (0, 1))
    assert (core["cap_algebra_dimension"], core["record_output_image_dimension"],
            core["relative_entropy_equality_algebra_dimension"]) == (6, 4, 4)


def test_quantum_response_rejects_hidden_control_wrong_endpoints_and_state_dependent_flags():
    paths = h.histories("E")
    with pytest.raises(ValueError, match="duplicate"):
        h.quantum_response((*paths, paths[0]), law="E")
    with pytest.raises(ValueError, match="fair-coin"):
        h.quantum_response((replace(paths[0], chance=F(1)), *paths[1:]), law="E")
    endpoint = native.states(binary=True)[(paths[0].initial+1) % 66]
    with pytest.raises(ValueError, match="endpoint"):
        h.quantum_response((replace(paths[0], endpoint=endpoint), *paths[1:]), law="E")
    with pytest.raises(ValueError, match="state-independent"):
        h.quantum_response(tuple(path for path in paths if path.initial != 0), law="E")
    # Every input gets a normalized unitary history, but its first seam is
    # chosen using a hidden initial occupation rather than a committed tape.
    hidden = []
    for initial, state in enumerate(native.states(binary=True)):
        edge = (0, 1) if state[0] else (0, 5)
        for coin in (0, 1):
            out, record = h._step(state, (), edge, coin, "E")
            hidden.append(h.Path(initial, out, record, F(1, 2), "success"))
    with pytest.raises(ValueError, match="committed records"):
        h.quantum_response(tuple(hidden), law="E")
    # A constant beta is not enough: hiding the coin branch selected using
    # an unobserved occupation yields beta=1/2 for every source input, yet
    # drops compatible inputs from individual observable tape domains.
    concealed = []
    for initial, state in enumerate(native.states(binary=True)):
        coin = state[4]
        out, tape = h._step(state, (), (0, 1), coin, "C")
        concealed.append(h.Path(initial, out, tape, F(1, 2), "success"))
    assert h.flag_probability(tuple(concealed)) == F(1, 2)
    with pytest.raises(ValueError, match="every compatible initial state"):
        h.quantum_response(tuple(concealed), law="C")


@pytest.mark.parametrize("bad", [
    True, np.zeros((66, 66), bool), [[False]*66 for _ in range(66)],
    np.ma.array(np.zeros((66, 66)), mask=False), np.zeros((65, 65)),
    np.full((66, 66), math.nan), np.eye(66),
])
def test_coherent_tangent_domain_rejects_coercion_nonfinite_and_changed_constraints(bad):
    response = h.quantum_response(h.histories("E"), law="E")
    with pytest.raises(ValueError):
        h.quantum_stationarity(response, local.reference(), bad)


def test_quantum_paths_refuse_nonoriginal_or_malformed_primitive_records():
    paths = h.histories("E")
    for tape in ((list(paths[0].tape[0]),), ((0, 5, True, 0),), ((0, 5, 0, 17),)):
        with pytest.raises(ValueError):
            h.quantum_response((replace(paths[0], tape=tape, chance=F(1, 2)),), law="E")
    response = h.quantum_response(paths, law="E")
    for cap in ((0, 1, []), (False, 1), (1, 0), (2, 3)):
        with pytest.raises(ValueError):
            h.quantum_region_response(response, cap)


def test_coherent_stationarity_refuses_each_changed_moment_and_unrepresentable_reference():
    response = h.quantum_response(h.histories("E"), law="E")
    indices = [v.charges().index(q) for q in range(3)]
    for direction in ((1, -1, 0), (1, -2, 1)):
        x = np.zeros((66, 66))
        x[indices, indices] = direction
        assert np.trace(x) == 0
        with pytest.raises(ValueError, match="both initial moments"):
            h.quantum_stationarity(response, local.reference(), x)
    x = np.zeros((66, 66), complex)
    x[0, 1] = 1j
    with pytest.raises(ValueError, match="Hermitian"):
        h.quantum_stationarity(response, local.reference(), x)
    tiny = v.family(F(1, 3), F(1, 10**1000))
    with pytest.raises(ValueError, match="execution precision"):
        h.quantum_stationarity(response, tiny, np.zeros((66, 66)))


def test_coherent_entropy_divided_difference_survives_an_underflowed_exact_gap():
    response = h.quantum_response(h.histories("E"), law="E")
    midpoint = (F(1e-308)+F(math.nextafter(1e-308, math.inf)))/2
    epsilon = F(1, 10**340)
    a, b = midpoint-epsilon, midpoint+epsilon
    assert float(a) != float(b) and float(b-a) == 0
    weights = v.family(20*a+2*b, b)
    result = h.quantum_stationarity(response, weights, np.zeros((66, 66)))
    assert result["source_entropy_hessian"] == result["output_entropy_hessian"] == 0
    assert result["information_loss_hessian"] == 0
    assert all(math.isfinite(value) for value in result.values() if type(value) is float)


def test_faithful_rare_populations_refuse_nonfinite_curvature_diagnostics():
    response = h.quantum_response(h.histories("E"), law="E")
    epsilon = F(1, 10**308)
    weights = v.family(22*epsilon, epsilon)
    pair = next((i, j) for group in response.recovery_classes for i in group for j in group
                if i < j and native.states(binary=True)[i][0] == 1)
    assert response.overlap[pair[0]][pair[1]] == 1
    tangent = np.zeros((66, 66))
    tangent[pair] = tangent[pair[::-1]] = 1
    # Each reciprocal population is representable, yet their two coherent
    # contributions overflow. A diagnostic must refuse rather than return
    # (-inf)-(-inf)=NaN as an allegedly lossless curvature.
    assert math.isfinite(1/float(min(weights)))
    with pytest.raises(ValueError, match="diagnostics must remain finite"):
        h.quantum_stationarity(response, weights, tangent)


def test_coherent_output_gradient_inherits_the_same_finite_preparation_sandwich():
    response = h.quantum_response(h.histories("E", prefix=3, stop_on_current=True), law="E")
    weights = local.reference()
    sigma = np.diag(list(map(float, weights))).astype(complex)
    a, b = next((i, j) for i in range(66) for j in range(i+1, 66)
                if response.overlap[i][j] == F(1, 2))
    x = np.zeros((66, 66), complex)
    x[a, b], x[b, a] = 1j, -1j
    radius = .1
    rho = sigma+radius*float(min(weights))*x
    assert np.linalg.eigvalsh(rho-(1-radius)*sigma).min() > -1e-16
    assert np.linalg.eigvalsh((1+radius)*sigma-rho).min() > -1e-16
    derivative = 0.
    for branch in response.branches:
        block = float(branch.weight)*rho[np.ix_(branch.inputs, branch.inputs)]
        direction = float(branch.weight)*x[np.ix_(branch.inputs, branch.inputs)]
        values, vectors = np.linalg.eigh(block)
        logarithm = (vectors*np.log(values))@vectors.conj().T
        derivative -= float(np.trace(direction@logarithm).real)
        derivative -= float(np.trace(direction).real)*(1+len(branch.tape)*math.log(2))
    assert 0 < abs(derivative) <= 2*radius/(1-radius)*2  # ||X||_1=2.
