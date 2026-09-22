"""Additional semantic and high-precision checks of the analytic matrix proof."""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F

import pytest

from . import check_quantum, codec, quantum


@pytest.mark.parametrize("spec", codec.QUANTUM_CASES)
def test_exact_noncommuting_certificate(spec):
    codec.equal(quantum.case(spec), check_quantum.segment(spec), "matrix certificate")


def decimal_matrix(a):
    return [[D(F(x).numerator)/D(F(x).denominator) for x in row] for row in a]


def roots(a):
    trace = a[0][0]+a[1][1]
    delta = ((a[0][0]-a[1][1])**2+4*a[0][1]*a[1][0]).sqrt()
    return (trace+delta)/2, (trace-delta)/2


def entropy(a):
    return sum(x*x.ln() if x else D(0) for x in roots(a))


def relative(a, reference):
    high, low = roots(reference)
    assert high >= low > 0
    if high == low:
        logarithm = [[high.ln(), D(0)], [D(0), high.ln()]]
    else:
        beta = (high.ln()-low.ln())/(high-low)
        alpha = (high*low.ln()-low*high.ln())/(high-low)
        logarithm = [[beta*reference[i][j]+(alpha if i==j else 0) for j in range(2)] for i in range(2)]
    return entropy(a)-sum(a[i][j]*logarithm[j][i] for i in range(2) for j in range(2))


@pytest.mark.parametrize("t", ("1/16", "1/4", "1/2"))
def test_noncommuting_entropy_identity_and_boundary_gap(t):
    # This is an independent numerical regression of the analytic identities,
    # not the proof or a numerical optimizer used to certify the receipt.
    with localcontext() as context:
        context.prec = 85
        item = check_quantum.segment((3, 4, 5))
        left, right = decimal_matrix(item["left"]), decimal_matrix(item["right"])
        reference = decimal_matrix([[F(2, 3), F(1, 6)], [F(1, 6), F(1, 3)]])
        t = D(F(t).numerator)/D(F(t).denominator)
        mix = [[(1-t)*left[i][j]+t*right[i][j] for j in range(2)] for i in range(2)]
        gap = (1-t)*relative(left, reference)+t*relative(right, reference)-relative(mix, reference)
        donald = (1-t)*relative(left, mix)+t*relative(right, mix)
        assert abs(gap-donald) < D("1e-70")
        bad = right[1][1]
        binary_entropy = -bad*bad.ln()-(1-bad)*(1-bad).ln()
        assert gap >= -t*bad*t.ln()-t*binary_entropy-D("1e-70")


def test_explicit_improvement_has_correct_sign_and_size():
    with localcontext() as context:
        context.prec = 85
        item = check_quantum.segment((3, 4, 5))
        left, right = decimal_matrix(item["left"]), decimal_matrix(item["right"])
        reference = decimal_matrix([[F(1, 2), 0], [0, F(1, 2)]])
        bad = right[1][1]
        h = -bad*bad.ln()-(1-bad)*(1-bad).ln()
        cost = relative(right, reference)+h
        t = (-1-cost/bad).exp()
        mix = [[(1-t)*left[i][j]+t*right[i][j] for j in range(2)] for i in range(2)]
        actual = relative(left, reference)-relative(mix, reference)
        assert 0 < t < 1 and actual >= bad*t > 0


def test_classical_record_face_does_not_select_identity_quantum_channel():
    identity, depolarizing, dephasing, flip = check_quantum.expected()["channels"]
    assert identity["uniform_bit_joint"] == dephasing["uniform_bit_joint"]
    assert identity["choi"] != dephasing["choi"]
    assert identity["matrix_unit_images"][1] != dephasing["matrix_unit_images"][1]
    assert [row["mismatch"] for row in (identity, depolarizing, dephasing, flip)] == ["0", "1/2", "0", "1"]
