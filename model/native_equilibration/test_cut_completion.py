"""Native Kraus, coherent/reference, variational and failed-extension controls."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest
from scipy.sparse import coo_matrix

from . import cut_completion as c
from ..observer import native, quantum_mixing


def entropy(rho):
    values = np.linalg.eigvalsh(rho)
    assert values.min() > -2e-12
    positive = values[values > 1e-14]
    return -np.dot(positive, np.log(positive))


def log_matrix(rho):
    values, vectors = np.linalg.eigh(rho)
    assert values.min() > 0
    return (vectors*np.log(values))@vectors.conj().T


def apply_with_reference(op, rho, ref=1):
    data = rho.reshape(66, ref, 66, ref).transpose(0, 2, 1, 3).reshape(4356, ref*ref)
    return (op@data).reshape(66, 66, ref, ref).transpose(0, 2, 1, 3).reshape(66*ref, 66*ref)


def independent_partial_trace_completion(rho):
    result = np.zeros_like(rho)
    for sector in c.cut_sectors():
        ix = np.array(sector.indices)
        a, b = ix.shape
        reduced = np.zeros((a, a), complex)
        for j in range(b):
            reduced += rho[np.ix_(ix[:, j], ix[:, j])]
        for j in range(b):
            result[np.ix_(ix[:, j], ix[:, j])] = reduced/b
    return result


def test_native_factorization_graph_certificate_and_exact_projection():
    report = c.run()
    cert = report["native_cut"]
    assert [(s["bulk_dimension"], s["complement_dimension"]) for s in cert["sector_dimensions"]] == [(1, 28), (4, 8), (6, 1)]
    assert cert["fixed_algebra_dimension"] == 53
    assert cert["internal_native_seams"] == 15
    assert F(cert["gap_lower_per_attempt"]) == F(1, 240)
    assert cert["preparation"]["attempts"] == 4080
    assert F(cert["preparation"]["trace_distance_upper"]) == F(33, 131072)
    assert F(cert["source_entropy_deficit_upper_nats"]) == F(4719, 1310720)
    assert report["all_matrix_unit_audit"]["completion_rank"] == 53
    assert report["all_matrix_unit_audit"]["completion_denominator"] == 56
    assert json.loads(c.OUTPUT.read_text()) == report


@pytest.mark.parametrize("precision", [1, 10, 30, 1000])
def test_entropy_budget_is_uniform_without_a_reference_population_floor(precision):
    import mpmath as mp
    from math import factorial
    assert sum(F(7, 10)**k/factorial(k) for k in range(5)) > 2
    cert = c.cut_certificate(precision_bits=precision)
    epsilon = F(cert["preparation"]["trace_distance_upper"])
    assert 0 < epsilon < F(1, 2)
    with mp.workdps(80):
        x = mp.mpf(epsilon.numerator)/epsilon.denominator
        independent = x*mp.log(65)-x*mp.log(x)-(1-x)*mp.log1p(-x)
        bound = F(cert["source_entropy_deficit_upper_nats"])
        assert independent < mp.mpf(bound.numerator)/bound.denominator


@pytest.mark.parametrize("law", ["C", "E"])
def test_every_actual_kraus_matrix_unit_matches_outside_signed_difference_formula(law):
    # Independent analytic expression: this does not invoke native repair.
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    edges = tuple((a, b) for a, b in quantum_mixing.seams() if a not in c.CAP and b not in c.CAP)
    rows, cols, values = [], [], []
    for i, x in enumerate(basis):
        for j, y in enumerate(basis):
            z = tuple(a-b for a, b in zip(x, y, strict=True))
            for a, b in edges:
                col = i*66+j
                if z[a] == z[b] == 0 and x[a] != x[b]:
                    xx, yy = list(x), list(y)
                    xx[a], xx[b] = xx[b], xx[a]
                    yy[a], yy[b] = yy[b], yy[a]
                    rows.extend((col, 66*positions[tuple(xx)]+positions[tuple(yy)]))
                    cols.extend((col, col)); values.extend((1, 1))
                else:
                    loss = abs(z[a]-z[b]) if law == "C" else 2*int(z[a] != z[b])
                    rows.append(col); cols.append(col); values.append(2-loss)
    predicted = coo_matrix((values, (rows, cols)), shape=(4356, 4356), dtype=np.int64).tocsr()
    actual, denominator = c.native_cut_counts(law=law)
    assert denominator == 30
    assert c._zero(actual-predicted)
    # Numerical spectrum diagnostics complement the analytic exact lower bound.
    groups = defaultdict(list)
    complement = tuple(i for i in range(12) if i not in c.CAP)
    for i, x in enumerate(basis):
        for j, y in enumerate(basis):
            key = tuple(x[v] for v in c.CAP), tuple(y[v] for v in c.CAP), tuple(x[v]-y[v] for v in complement)
            groups[key].append(66*i+j)
    eigenvalues = np.concatenate([np.linalg.eigvalsh(actual[ix][:, ix].toarray()/denominator) for ix in groups.values()])
    assert np.count_nonzero(abs(eigenvalues-1) < 1e-11) == 53
    assert eigenvalues.min() > -1e-12
    assert eigenvalues[eigenvalues < 1-1e-11].max() <= 1-1/240+1e-12


@pytest.mark.parametrize("law", ["C", "E"])
def test_actual_finite_repair_retains_cap_reference_and_completes_environment(law):
    p, ps = c.native_cut_counts(law=law)
    e, es = c.completion_counts()
    # The source is entangled with a two-dimensional untouched reference;
    # three cap charges occur and one cap charge has a nontrivial coherence.
    ket = np.zeros((66, 2), complex)
    sectors = c.cut_sectors()
    ket[sectors[0].indices[0][0], 0] = 1/2
    ket[sectors[1].indices[0][0], 0] = 1j/2
    ket[sectors[1].indices[1][1], 1] = -1/2
    ket[sectors[2].indices[0][0], 1] = 1/2
    rho = np.outer(ket.ravel(), ket.ravel().conj())
    expected = apply_with_reference(e/es, rho, 2)
    initial_reference = np.trace(rho.reshape(66, 2, 66, 2), axis1=0, axis2=2)
    budget = c.cut_certificate()["preparation"]
    channel = p/ps
    for _ in range(budget["attempts"]):
        rho = apply_with_reference(channel, rho, 2)
    error = np.abs(np.linalg.eigvalsh(rho-expected)).sum()/2
    assert error <= float(F(budget["trace_distance_upper"]))+2e-12
    assert abs(np.trace(rho)-1) < 3e-12
    np.testing.assert_allclose(np.trace(rho.reshape(66, 2, 66, 2), axis1=0, axis2=2), initial_reference, atol=2e-12)
    np.testing.assert_allclose(apply_with_reference(e/es, rho, 2), expected, atol=2e-12)
    # A complete computational-basis readout, including every outcome.
    assert np.abs(np.diag(rho-expected)).sum()/2 <= error+1e-13


def test_source_entropy_maximum_modular_split_and_all_fibre_tangents():
    rng = np.random.default_rng(82)
    ket = rng.normal(size=66)+1j*rng.normal(size=66)
    ket /= np.linalg.norm(ket)
    rho = 0.6*np.eye(66)/66+0.4*np.outer(ket, ket.conj())
    e, es = c.completion_counts()
    sigma = apply_with_reference(e/es, rho)
    np.testing.assert_allclose(sigma, independent_partial_trace_completion(rho), atol=1e-15)
    ls = log_matrix(sigma)
    divergence = np.trace(rho@(log_matrix(rho)-ls)).real
    assert divergence > 0
    assert abs((entropy(sigma)-entropy(rho))-divergence) < 1e-12
    bulk_entropy, edge_entropy = 0., 0.
    bulk_generator, center_generator = np.zeros((66, 66), complex), np.zeros((66, 66), complex)
    for s in c.cut_sectors():
        ix = np.array(s.indices)
        a, b = ix.shape
        cap_block = sum(sigma[np.ix_(ix[:, j], ix[:, j])] for j in range(b))
        probability = np.trace(cap_block).real
        bulk_entropy += -probability*np.log(probability)+probability*entropy(cap_block/probability)
        edge_entropy += probability*np.log(b)
        for j in range(b):
            bulk_generator[np.ix_(ix[:, j], ix[:, j])] = -log_matrix(cap_block)
            center_generator[np.ix_(ix[:, j], ix[:, j])] = np.eye(a)*np.log(b)
    np.testing.assert_allclose(-ls, bulk_generator+center_generator, atol=2e-14)
    assert abs(entropy(sigma)-bulk_entropy-edge_entropy) < 2e-12
    # Because E is HS-self-adjoint and fixes log(sigma), this is a complete
    # dual witness that every Hermitian X with E(X)=0 has zero derivative.
    np.testing.assert_allclose(apply_with_reference(e/es, ls), ls, atol=2e-14)
    y = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    y = (y+y.conj().T)/2
    x = y-apply_with_reference(e/es, y)
    x /= np.linalg.norm(x, 2)
    assert abs(np.trace(x@ls)) < 3e-13
    for step in (-1e-4, 1e-4):
        varied = sigma+step*x
        assert np.linalg.eigvalsh(varied).min() > 0
        assert entropy(varied) < entropy(sigma)
        np.testing.assert_allclose(independent_partial_trace_completion(varied), sigma, atol=1e-15)


@pytest.mark.parametrize("law", ["C", "E"])
def test_full_empty_complement_coherence_blocks_unqualified_extension(law):
    # B=(0,1) is connected, but both empty and full are possible at k=2.
    cap = tuple(range(2, 12))
    p, ps = c.native_cut_counts(cap, law)
    e, es = c.completion_counts(cap)
    basis = native.states(binary=True)
    x = basis.index((1, 1)+(0,)*10)
    y = basis.index((0, 0, 1, 1)+(0,)*8)
    unit = np.zeros(4356); unit[66*x+y] = 1
    np.testing.assert_array_equal(p@unit, ps*unit)
    np.testing.assert_array_equal(e@unit, np.zeros(4356))
    with pytest.raises(ValueError, match="empty/full coherence"):
        c.cut_certificate(cap)


def test_palm_reference_and_conditioned_records_are_different_contracts():
    basis = native.states(binary=True)
    degrees = np.array([sum(x[a] != x[b] for a, b in quantum_mixing.seams()) for x in basis])
    assert set(degrees) == {8, 10} and degrees.sum() == 600
    e, es = c.completion_counts()
    palm = np.diag(degrees/600)
    completed = apply_with_reference(e/es, palm)
    np.testing.assert_allclose(np.diag(completed), np.asarray(c.completed_palm(), float), atol=1e-16)
    assert min(c.completed_palm()) == F(1, 75)
    assert [sum(c.completed_palm()[i] for row in s.indices for i in row) for s in c.cut_sectors()] == [F(5, 12), F(1, 2), F(1, 12)]
    for s in c.cut_sectors():
        for row in s.indices:
            assert len({c.completed_palm()[i] for i in row}) == 1
    assert c.completed_palm()[c.cut_sectors()[1].indices[2][0]]/c.completed_palm()[c.cut_sectors()[1].indices[0][0]] == F(38, 37)
    assert abs(np.abs(np.diag(completed-palm)).sum()/2-8/175) < 1e-15
    # An accepted movement flag from the uniform state breaks conditional
    # expectation completion: records cannot silently be retained/conditioned.
    seam = next((a, b) for a, b in quantum_mixing.seams() if a not in c.CAP and b not in c.CAP)
    plus = np.zeros(66)
    for state, coin in product(basis, (0, 1)):
        output, (_, current) = native.native(state, seam, coin)
        if current == 1:
            plus[basis.index(output)] += 1/132
    assert abs(plus.sum()-5/66) < 1e-15
    posterior = np.diag(plus/plus.sum())
    assert entropy(apply_with_reference(e/es, posterior))-entropy(posterior) > 0.1


@pytest.mark.parametrize("cap", [(), (0, 0), (2, 1), (True,), (12,), list(c.CAP), tuple(range(12))])
def test_malformed_cap_rejected(cap):
    with pytest.raises(ValueError):
        c.cut_certificate(cap)
    with pytest.raises(ValueError):
        c.cut_sectors(cap)
