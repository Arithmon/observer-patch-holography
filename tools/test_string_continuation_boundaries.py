"""Exact controls for the string paper's continuum inference boundaries."""

from collections import defaultdict
from pathlib import Path


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def test_odd_square_relation_is_already_inhabited_by_two_state_mechanics():
    q = ((0, 1), (1, 0))
    parity = ((1, 0), (0, -1))
    identity = ((1, 0), (0, 1))
    assert matmul(q, q) == identity
    qf, fq = matmul(q, parity), matmul(parity, q)
    assert all(qf[i][j] + fq[i][j] == 0 for i in range(2) for j in range(2))
    # No local current or conformal modes occur in the premise. In a finite
    # Hilbert space tr([L1,L-1])=0 also prevents [L1,L-1]=2 H for this H.
    assert 2 * sum(identity[i][i] for i in range(2)) == 4


def convolve(a, b, cutoff):
    out = defaultdict(int)
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= cutoff:
                out[i + j] += x * y
    return {i: x for i, x in out.items() if x}


def power(series, degree, cutoff):
    result = {0: 1}
    for _ in range(degree):
        result = convolve(result, series, cutoff)
    return result


def test_supersymmetric_spin_sum_cancels_even_with_nonzero_characters():
    # Formal variable x=q^(1/8): theta3=sum x^(4n²), theta4 has (-1)^n,
    # theta2=sum x^((2n+1)²). Retain all coefficients through q^10 exactly.
    cutoff = 80
    theta2, theta3, theta4 = (defaultdict(int) for _ in range(3))
    for n in range(-cutoff, cutoff + 1):
        if 4 * n * n <= cutoff:
            theta3[4 * n * n] += 1
            theta4[4 * n * n] += 1 if n % 2 == 0 else -1
        if (2 * n + 1) ** 2 <= cutoff:
            theta2[(2 * n + 1) ** 2] += 1
    fourth2, fourth3, fourth4 = (power(s, 4, cutoff)
                                for s in (theta2, theta3, theta4))
    assert fourth3[0] == fourth4[0] == 1
    assert fourth2[4] == 16
    assert all(fourth3.get(i, 0) - fourth4.get(i, 0) - fourth2.get(i, 0) == 0
               for i in range(cutoff + 1))
    # This is a finite coefficient regression; the all-order identity is
    # Jacobi's identity cited in the paper, not inferred from the truncation.


def test_string_paper_requires_local_ope_and_full_sector_identification():
    paper = (Path(__file__).resolve().parents[1]
             / "extra/observer_patch_holography_as_string_vacuum_selector.tex").read_text()
    assert "The integrated square relation alone does not establish this local OPE" in paper
    assert "Ordinary, unprojected vacuum characters" in paper
    assert "They are not equivalent to the full edge-VOA identification" in paper
    assert "field-content, charge-lattice, and sector-extension receipts" in paper
