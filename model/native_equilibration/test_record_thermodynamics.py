"""Independent full-state, complete-tape and failed-extension controls."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest

from . import record_thermodynamics as r
from . import cut_completion as c
from ..observer import native, quantum_mixing


def entropy(a):
    values = np.linalg.eigvalsh(a)
    assert values.min() > -2e-13
    positive = values[values > 1e-15]
    return float(-np.dot(positive, np.log(positive)))


def relative(a, b):
    values, vectors = np.linalg.eigh(b)
    assert values.min() > 0
    logb = (vectors*np.log(values))@vectors.conj().T
    return -entropy(a)-float(np.trace(a@logb).real)


def expand(indices, ref):
    return np.array([ref*i+j for i in indices for j in range(ref)])


def reference_marginal(rho, ref):
    return np.trace(rho.reshape(66, ref, 66, ref), axis1=0, axis2=2)


def source_marginal(rho, ref):
    return np.trace(rho.reshape(66, ref, 66, ref), axis1=1, axis2=3)


def apply_native(rho, law, ref=1, cap=c.CAP):
    p, ps = c.native_cut_counts(cap, law)
    matrix_units = rho.reshape(66, ref, 66, ref).transpose(0, 2, 1, 3).reshape(4356, ref*ref)
    return ((p@matrix_units)/ps).reshape(66, 66, ref, ref).transpose(0, 2, 1, 3).reshape(66*ref, 66*ref)


def coherent_state(ref=2):
    # Coherent/entangled full-rank input, including all charge sectors.
    rng = np.random.default_rng(413)
    a = rng.normal(size=(66*ref, 3))+1j*rng.normal(size=(66*ref, 3))
    a = a@a.conj().T
    return .72*a/np.trace(a)+.28*np.eye(66*ref)/(66*ref)


def tape_ledger(rho, tau, laws, word_alphabet):
    """Independent block-cq computation, never constructing an enormous tape."""
    ref = len(rho)//66
    choices_count = len(word_alphabet)
    output = np.zeros_like(rho)
    tape_entropy = tape_relative = pinched_entropy = record_entropy = 0.
    record_conditioned_entropy = 0.
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    for word in word_alphabet:
        groups = defaultdict(list)
        for i, initial in enumerate(basis):
            state, record = initial, []
            for (seam, coin), law in zip(word, laws, strict=True):
                new, (_, current) = native.native(state, seam, coin)
                record.append(current if law == "C" else state[seam[0]]-state[seam[1]])
                state = new
            groups[tuple(record)].append((i, positions[state]))
        for pairs in groups.values():
            inputs, outputs = zip(*pairs, strict=True)
            ix, iy = expand(inputs, ref), expand(outputs, ref)
            a, b = rho[np.ix_(ix, ix)], tau[np.ix_(ix, ix)]
            block, blockref = a/choices_count, b/choices_count
            output[np.ix_(iy, iy)] += block
            tape_entropy += entropy(block)
            tape_relative += relative(block, blockref)
            pinched_entropy += entropy(a)/choices_count
            probability = np.trace(block).real
            record_entropy -= probability*np.log(probability)
            record_conditioned_entropy += entropy(block)+probability*np.log(probability)
    return output, tape_entropy, tape_relative, pinched_entropy, record_entropy, record_conditioned_entropy


def test_exact_source_and_history_receipt():
    report = r.run()
    assert json.loads(r.OUTPUT.read_text()) == report
    for law in ("C", "E"):
        audit = report["laws"][law]
        assert audit["detailed_balance"]["choice_alphabet"] == 30
        assert audit["detailed_balance"]["record_branches"] == (60 if law == "C" else 90)
        assert audit["complete_history"]["all_input_choice_pairs"] == 59400
        assert audit["complete_history"]["all_choice_words"] == 900


@pytest.mark.parametrize("law", ["C", "E"])
def test_coherent_complete_tape_clausius_and_relative_entropy_chain(law):
    rho = coherent_state()
    rho_ref = reference_marginal(rho, 2)
    tau = np.kron(np.diag([float(x) for x in c.completed_palm()]), rho_ref)
    words = tuple((choice,) for choice in r.choices())
    out, s_tape, d_tape, s_pinched, h_j, mean_s = tape_ledger(rho, tau, (law,), words)
    np.testing.assert_allclose(out, apply_native(rho, law, 2), atol=2e-16)
    np.testing.assert_allclose(apply_native(tau, law, 2), tau, atol=2e-16)
    np.testing.assert_allclose(reference_marginal(out, 2), rho_ref, atol=2e-16)
    d_in, d_out = relative(rho, tau), relative(out, tau)
    q, delta = s_pinched-entropy(rho), entropy(out)-entropy(rho)
    assert 1e-3 < q < np.log(66)
    assert 1e-3 < delta-q < np.log(30)
    assert abs(delta-(d_in-d_out)) < 2e-13
    assert abs(s_tape-(np.log(30)+entropy(rho)+q)) < 2e-13
    assert abs(d_tape-(d_in-q)) < 2e-13
    assert abs(s_tape-(h_j+mean_s)) < 2e-13
    conditional_tape = s_tape-entropy(out)
    chi = entropy(out)-mean_s
    assert 0 < conditional_tape < h_j
    assert 0 < chi < h_j
    assert abs(conditional_tape-(np.log(30)+q-delta)) < 2e-13
    # Ignoring initial source-reference correlations changes the entropy ledger.
    source_delta = entropy(source_marginal(out, 2))-entropy(source_marginal(rho, 2))
    assert delta-source_delta > .01
    assert abs((delta-source_delta)-(relative(rho, np.kron(source_marginal(rho, 2), rho_ref))
                -relative(out, np.kron(source_marginal(out, 2), rho_ref)))) < 2e-13


@pytest.mark.parametrize("law", ["C", "E"])
def test_two_step_coherent_tape_and_noncommuting_seam_schedule(law):
    rho = coherent_state(ref=1)
    tau = np.diag([float(x) for x in c.completed_palm()])
    seams = [s for s in quantum_mixing.seams() if not set(s).intersection(c.CAP)]
    first, second = next((a, b) for a, b in product(seams, repeat=2)
                         if len(set(a).intersection(b)) == 1)
    # These actual two-seam compositions need not be self-adjoint; the
    # entropy identities require reference preservation, not that upgrade.
    words = tuple(((first, a), (second, b)) for a, b in product((0, 1), repeat=2))
    out, s_tape, d_tape, s_pinched, _, _ = tape_ledger(rho, tau, (law, law), words)
    q = s_pinched-entropy(rho)
    delta = entropy(out)-entropy(rho)
    assert 0 < q <= delta
    assert delta-q < np.log(4)
    assert abs(s_tape-(np.log(4)+entropy(rho)+q)) < 2e-13
    assert abs(d_tape-(relative(rho, tau)-q)) < 2e-13
    assert abs(delta-(relative(rho, tau)-relative(out, tau))) < 2e-13


@pytest.mark.parametrize("law", ["C", "E"])
def test_retained_population_tape_has_exact_likelihood_ratio_and_joint_entropy(law):
    sigma = c.completed_palm()
    initial = tuple(F(i+1, 2211) for i in range(66))
    output = [F(0)]*66
    s_joint = d_joint = 0.
    for word in product(r.choices(), repeat=2):
        for _, pairs in r.history_domains(word, law):
            for i, j in pairs:
                p, q = initial[i]/900, sigma[i]/900
                assert p/q == initial[i]/sigma[i]
                output[j] += p
                s_joint -= float(p)*np.log(float(p))
                d_joint += float(p)*np.log(float(p/q))
    h = lambda p: -sum(float(x)*np.log(float(x)) for x in p)
    d_initial = sum(float(p)*np.log(float(p/q)) for p, q in zip(initial, sigma, strict=True))
    assert sum(output) == 1
    assert abs(s_joint-h(initial)-2*np.log(30)) < 3e-11
    assert abs(d_joint-d_initial) < 3e-12
    assert h(output) > h(initial)
    conditional_tape = s_joint-h(output)
    assert 0 < conditional_tape < 2*np.log(30)


@pytest.mark.parametrize("law", ["C", "E"])
def test_gns_kms_weighted_symmetry_with_complex_observables(law):
    rng = np.random.default_rng(755)
    x, y = [rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66)) for _ in range(2)]
    sigma = np.diag([float(x) for x in c.completed_palm()])
    square = np.sqrt(sigma)
    px, py = apply_native(x, law), apply_native(y, law)
    for inner in (lambda a, b: np.trace(sigma@a.conj().T@b),
                  lambda a, b: np.trace(square@a.conj().T@square@b)):
        assert abs(inner(x, py)-inner(px, y)) < 5e-14


@pytest.mark.parametrize("law", ["C", "E"])
def test_selected_record_and_wrong_reference_do_not_inherit_clausius(law):
    sigma = np.diag([float(x) for x in c.completed_palm()])
    branch = next(b for b in r.branches(law=law) if any(i != j for i, j in b.pairs))
    inputs, outputs = zip(*branch.pairs, strict=True)
    selected = np.zeros((66, 66))
    selected[np.ix_(outputs, outputs)] = sigma[np.ix_(inputs, inputs)]
    selected /= np.trace(selected)
    assert entropy(selected) < entropy(sigma)-1
    assert relative(selected, sigma) > 1
    degree = np.array([sum(x[a] != x[b] for a, b in quantum_mixing.seams())
                       for x in native.states(binary=True)], float)
    palm = np.diag(degree/degree.sum())
    out = apply_native(palm, law)
    assert entropy(out) > entropy(palm)
    assert relative(out, palm) > 1e-5  # D descent against unpreserved Palm fails.


def test_cross_cut_repair_does_not_preserve_the_same_reference():
    sigma = c.completed_palm()
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    witnesses = []
    for seam in quantum_mixing.seams():
        if len(set(seam).intersection(c.CAP)) != 1:
            continue
        for i, x in enumerate(basis):
            for coin in (0, 1):
                y, _ = native.native(x, seam, coin)
                if sigma[i] != sigma[positions[y]]:
                    witnesses.append((seam, coin, i))
    assert witnesses


@pytest.mark.parametrize("law", ["C", "E"])
def test_all_initial_domains_including_empty_history_and_api_guards(law):
    assert r.history_domains((), law) == (((), tuple((i, i) for i in range(66))),)
    assert r.exact_history_audit(law=law, attempts=0)["all_labelled_words"] == 1
    with pytest.raises(ValueError):
        r.exact_history_audit(law=law, attempts=3)
    with pytest.raises(ValueError):
        r.branches(law="unknown")
    with pytest.raises(ValueError):
        r.history_domains((((0, 0), 0),), law)


@pytest.mark.parametrize("law", ["C", "E"])
def test_singular_untouched_reference_uses_its_support_without_extra_entropy(law):
    rho = coherent_state(ref=2)
    # Embed a correlated two-dimensional reference into a three-dimensional
    # one. The reference has a genuine zero eigenvalue; the source proof
    # must restrict to its support, never add a faithful-state regularizer.
    u = np.array([[1, 0], [0, 1/np.sqrt(2)], [0, 1j/np.sqrt(2)]], complex)
    embedding = np.kron(np.eye(66), u)
    larger = embedding@rho@embedding.conj().T
    out = apply_native(larger, law, 3)
    expected = embedding@apply_native(rho, law, 2)@embedding.conj().T
    np.testing.assert_allclose(out, expected, atol=2e-16)
    np.testing.assert_allclose(reference_marginal(out, 3), reference_marginal(larger, 3), atol=2e-16)
    assert abs(entropy(out)-entropy(larger)
               -(entropy(apply_native(rho, law, 2))-entropy(rho))) < 2e-13
