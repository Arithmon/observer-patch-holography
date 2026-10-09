"""Original-link controls for the electric plaquette-potential reduction."""
import itertools
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np
import mpmath
import pytest
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge_sectors import abelian_ground_state as reduced


# Independent incidence tables in documented 00x,00y,01x,01y,10x,10y,11x,11y
# order. Neither the old producer nor the new reduced generator is used here.
PLAQUETTES = np.array([
    [1, -1, -1, 0, 0, 1, 0, 0], [-1, 0, 1, -1, 0, 0, 0, 1],
    [0, 1, 0, 0, 1, -1, -1, 0], [0, 0, 0, 1, -1, 0, 1, -1],
])
STARS = np.array([
    [1, 1, 0, -1, -1, 0, 0, 0], [0, -1, 1, 1, 0, 0, -1, 0],
    [-1, 0, 0, 0, 1, 1, 0, -1], [0, 0, -1, 0, 0, -1, 1, 1],
])


def original_link_distribution(n, h):
    """Build all n**8 link states directly; project the full electric vector."""
    states = np.array(list(itertools.product(range(n), repeat=8)))
    diagonal = -np.cos(2 * np.pi * (states @ PLAQUETTES.T) / n).sum(axis=1)
    sources = np.arange(n**8)
    row, column, values = [], [], []
    for shift, strength in [(s, h) for s in np.eye(8, dtype=int)] + [(s, 5.0) for s in STARS]:
        for sign in (-1, 1):
            row.append(np.ravel_multi_index(((states + sign * shift) % n).T, (n,) * 8))
            column.append(sources)
            values.append(np.full(n**8, -strength / 2))
    operator = diags(diagonal) + coo_matrix(
        (np.concatenate(values), (np.concatenate(row), np.concatenate(column))),
        shape=(n**8, n**8),
    ).tocsr()
    if n == 2:
        _, vectors = eigh(operator.toarray(), subset_by_index=[0, 0])
    else:
        _, vectors = eigsh(operator, k=1, which="SA", tol=1e-13,
                           v0=np.ones(n**8), maxiter=20000)
    electric = np.fft.ifftn(vectors[:, 0].reshape((n,) * 8), norm="ortho").ravel()
    charges = (states[:, 0] + states[:, 1] - states[:, 3]) % n
    return np.bincount(charges, weights=abs(electric)**2, minlength=n)


def original_electric_sector_distribution(n, h):
    """Restrict the original electric-link action using explicit plaquette moves.

    Unlike the implementation's potential-coordinate rolls, this builds and
    looks up eight-link fields and changes them by the original incidence rows.
    """
    fields = np.array([[0, *abc] for abc in itertools.product(range(n), repeat=3)]) @ PLAQUETTES % n
    lookup = {tuple(field): i for i, field in enumerate(fields)}
    operator = np.diag(-h * np.cos(2 * np.pi * fields / n).sum(axis=1) - 20)
    for source, field in enumerate(fields):
        for plaquette in PLAQUETTES:
            for sign in (-1, 1):
                operator[lookup[tuple((field + sign * plaquette) % n)], source] -= 0.5
    _, vectors = eigh(operator, subset_by_index=[0, 0])
    charges = (fields[:, 0] + fields[:, 1] - fields[:, 3]) % n
    return np.bincount(charges, weights=vectors[:, 0]**2, minlength=n)


def high_precision_z2_probabilities(h):
    """Independently diagonalize the original electric-link operator at 800 dps."""
    ctx = mpmath.mp.clone()
    ctx.dps = 800
    fields = np.array([[0, *abc] for abc in itertools.product(range(2), repeat=3)]) @ PLAQUETTES % 2
    lookup = {tuple(field): i for i, field in enumerate(fields)}
    operator = ctx.zeros(8)
    for source, field in enumerate(fields):
        operator[source, source] = 2 * ctx.mpf(h) * int(field.sum())
        for plaquette in PLAQUETTES:
            operator[lookup[tuple((field + plaquette) % 2)], source] -= 1
    _, vectors = ctx.eigsy(operator)
    probabilities = [ctx.mpf(0), ctx.mpf(0)]
    for row, field in enumerate(fields):
        charge = int((field[0] + field[1] - field[3]) % 2)
        probabilities[charge] += vectors[row, 0]**2
    return np.array([float(probability) for probability in probabilities])


@pytest.mark.parametrize("h", [0.01, 0.2, 0.5, 1.0, 2.0, 10.0])
def test_z2_reduction_matches_independent_full_link_dense_ground(h):
    np.testing.assert_allclose(reduced.zn_edge_distribution(2, h),
                               original_link_distribution(2, h), rtol=2e-11, atol=1e-15)


@pytest.mark.parametrize("h", [0.2, 0.5, 1.0, 2.0])
def test_z3_reduction_matches_independent_full_link_sparse_ground(h):
    np.testing.assert_allclose(reduced.zn_edge_distribution(3, h),
                               original_link_distribution(3, h), rtol=2e-10, atol=1e-14)


def test_composite_modulus_zero_flux_space_is_exactly_plaquette_boundaries():
    n = 4
    potentials = np.array([[0, *abc] for abc in itertools.product(range(n), repeat=3)])
    boundaries = potentials @ PLAQUETTES % n
    all_fields = np.array(list(itertools.product(range(n), repeat=8)))
    gauss = np.all((all_fields @ STARS.T) % n == 0, axis=1)
    winding_x = (all_fields[:, 0] + all_fields[:, 2]) % n
    winding_y = (all_fields[:, 1] + all_fields[:, 5]) % n
    allowed = all_fields[gauss & (winding_x == 0) & (winding_y == 0)]
    assert len(allowed) == n**3
    assert set(map(tuple, boundaries)) == set(map(tuple, allowed))
    restricted_star = (boundaries[:, 0] + boundaries[:, 1] - boundaries[:, 3]) % n
    np.testing.assert_array_equal(restricted_star, (potentials[:, 2] - potentials[:, 3]) % n)


@pytest.mark.parametrize("h", [0.05, 0.1, 0.2, 0.5])
def test_z5_default_cases_resolve_positive_conjugate_sectors(h):
    probabilities = reduced.zn_edge_distribution(5, h)
    assert np.all(probabilities > 0)
    assert probabilities.sum() == pytest.approx(1, abs=3e-16)
    np.testing.assert_allclose(probabilities[1:], probabilities[:0:-1], rtol=1e-11, atol=0)
    np.testing.assert_allclose(probabilities, original_electric_sector_distribution(5, h),
                               rtol=2e-11, atol=0)


@pytest.mark.parametrize("n", [1, 0, -2, True, 2.0])
def test_reject_invalid_modulus(n):
    with pytest.raises(ValueError, match="integer >= 2"):
        reduced.zn_edge_distribution(n, 1)


@pytest.mark.parametrize("h", [float("nan"), float("inf"), -1, True, "1",
                             np.ma.array(1.0, mask=True), np.ma.masked])
def test_reject_invalid_coupling(h):
    with pytest.raises(ValueError, match="finite"):
        reduced.zn_edge_distribution(2, h)


@pytest.mark.parametrize("h", [2**53 + 1, Fraction(1, 3), Fraction(1, 2**1075)])
def test_original_coupling_cannot_round_during_binary64_conversion(h):
    with pytest.raises(ValueError, match="exactly representable in binary64"):
        reduced.zn_edge_distribution(2, h)


def test_exact_dyadic_coupling_is_accepted():
    np.testing.assert_array_equal(reduced.zn_edge_distribution(2, Fraction(1, 2)),
                                  reduced.zn_edge_distribution(2, 0.5))


def test_zero_coupling_does_not_silently_select_a_degenerate_topological_sector():
    with pytest.raises(ValueError, match="degenerate topological sectors"):
        reduced.zn_edge_distribution(2, 0)


@pytest.mark.parametrize("corruption", ["zero", "negative", "nonfinite", "constant"])
def test_unresolved_eigenvector_cannot_become_a_successful_distribution(monkeypatch, corruption):
    vector = np.ones(8)
    if corruption == "zero":
        vector[3] = 0
    elif corruption == "negative":
        vector[3] = -1
    elif corruption == "nonfinite":
        vector[3] = np.nan
    # A constant positive vector is the wrong ground state at h=.2 and must
    # fail its local eigen-equations, despite looking normalized and positive.
    monkeypatch.setattr(reduced, "eigsh", lambda *a, **kw: (np.array([-4.0]), vector[:, None]))
    with pytest.raises(ValueError, match="precision unresolved"):
        reduced.zn_edge_distribution(2, 0.2)


def test_corrupted_positive_schur_amplitudes_fail_componentwise_check(monkeypatch):
    monkeypatch.setattr(reduced, "_schur_ground", lambda *a: np.ones(8))
    with pytest.raises(ValueError, match="componentwise"):
        reduced.zn_edge_distribution(2, 1e6)


@pytest.mark.parametrize("h", [1e80, 1e100, 1e154])
def test_representable_sectors_survive_individual_weight_underflow(h):
    expected = high_precision_z2_probabilities(h)
    measured = reduced.zn_edge_distribution(2, h)
    assert np.all(expected > 0)
    np.testing.assert_allclose(measured, expected, rtol=1e-12, atol=0)


def test_underresolved_subnormal_sector_conversion_is_explicit():
    with pytest.raises(ValueError, match="sector probability conversion"):
        reduced.zn_edge_distribution(2, 1e156)


def test_unresolved_schur_underflow_is_explicit_and_never_zero_probability():
    with pytest.raises(ValueError, match="precision unresolved"):
        reduced.zn_edge_distribution(2, 1e200)


@pytest.mark.parametrize("amplitude", [1e-155, 1e-160, 1e-200])
def test_final_probability_conversion_respects_original_amplitude_value(amplitude):
    exact_weight = Fraction(amplitude)**2
    exact_probability = exact_weight / (1 + exact_weight)
    expected = float(exact_probability)
    relative_loss = abs(Fraction(expected) - exact_probability) / exact_probability
    if relative_loss > Fraction(1, 10**12):
        with pytest.raises(ValueError, match="sector probability (conversion|underflow)"):
            reduced._sector_probabilities([1.0, amplitude], [0, 1], 2)
    else:
        measured = reduced._sector_probabilities([1.0, amplitude], [0, 1], 2)
        assert 0 < measured[1] < np.finfo(float).tiny
        assert measured[1] == expected
