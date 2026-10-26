"""Cross-sector CAR signs, native records and invisible modular extensions."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import charged_sector as charged
from ...native_modular_reference import local_preparation
from ...observer import native, quantum_mixing


def _compose(mask, first, second):
    middle = first(mask)
    if middle is None:
        return None
    final = second(middle[0])
    return None if final is None else (final[0], middle[1]*final[1])


def test_adjacent_sectors_recover_car_on_the_actual_total_two_domain():
    # The source shell alone annihilates every c_i. Adjacent sectors restore
    # both CAR contractions, including signs for non-neighbouring ports.
    for mask in charged.masks(2):
        for i, j in product(range(12), repeat=2):
            terms = defaultdict(int)
            for term in (
                _compose(mask, lambda m: charged.creation(m, j),
                         lambda m: charged.annihilation(m, i)),
                _compose(mask, lambda m: charged.annihilation(m, i),
                         lambda m: charged.creation(m, j)),
            ):
                if term is not None:
                    terms[term[0]] += term[1]
            assert {m: sign for m, sign in terms.items() if sign} == ({mask: 1} if i == j else {})


def test_extension_agrees_with_every_existing_binary_source_branch():
    assert charged.primitive_certificate() == 2*30*2*66
    assert [len(charged.masks(load)) for load in range(4)] == [1, 12, 66, 220]


@pytest.mark.parametrize("law", ["C", "E"])
def test_full_integer_primitive_preserves_difference_with_identical_records(law):
    # This includes nonbinary loads and different total loads; it tests the
    # primitive formula independently of the truncated binary completion.
    seam = quantum_mixing.seams()[0]
    for a, b, c, d, coin in product(range(4), range(4), range(4), range(4), (0, 1)):
        x, y = [0]*12, [0]*12
        x[seam[0]], x[seam[1]] = a, b
        y[seam[0]], y[seam[1]] = c, d
        xp, rx = native.native(x, seam, coin)
        yp, ry = native.native(y, seam, coin)
        same_record = rx == ry if law == "C" else a-b == c-d
        if same_record:
            assert tuple(u-v for u, v in zip(xp, yp, strict=True)) == tuple(u-v for u, v in zip(x, y, strict=True))


@pytest.mark.parametrize("law", ["C", "E"])
def test_complete_proposal_has_no_off_port_charged_response(law):
    accepted, terminal = charged.response(law)
    timeout = local_preparation.selectors(1)["completion_timeout_probability"]
    assert 0 < timeout < 1
    for i, j in product(range(12), repeat=2):
        assert charged.read_field(charged.field_insertion(j), i) == (i == j)
        assert terminal[i][j] == (timeout if i == j else 0)
        if i != j:
            assert accepted[i][j] == 0
    # Genuine attenuation exists. Omitting all accepted outcomes must fail.
    assert any(0 < accepted[i][i] < 1-timeout for i in range(12))


def test_distinct_faithful_reference_extensions_preserve_the_conditioned_source():
    a, b = charged.reference_extension(1), charged.reference_extension(F(1, 2))
    assert a["partition_function_exact"] != b["partition_function_exact"]
    for witness in (a, b):
        probabilities = list(map(F, witness["number_sector_probabilities_exact"]))
        assert all(value > 0 for value in probabilities)
        assert sum(probabilities) == 1
        assert witness["conditional_total_two_reference_identical_exact"]
    # A neutral two-copy (3,1)<->(2,2) interferometer sees the discrete
    # second difference: a common chemical potential cancels, curvature does
    # not. This is exact arithmetic on the specified witness family.
    for kappa in (F(0), F(1, 7), F(2)):
        for mu in (F(-11, 3), F(0), F(5, 9)):
            f = lambda n: kappa*(n-2)**2+mu*n
            assert f(3)+f(1)-2*f(2) == 2*kappa


def test_coherent_hopping_control_has_cross_port_charged_response():
    # Exact two-mode CAR control: the response evaluator is not restricted
    # to diagonal propagation by definition. This continuation is supplied;
    # it is outside the current-record native channel tested above.
    operators = []
    for port in (0, 1):
        operator = s.zeros(4)
        for mask in range(4):
            removed = charged.annihilation(mask, port)
            if removed is not None:
                operator[removed[0], mask] = removed[1]
        operators.append(operator)
    c0, c1 = operators
    hopping = c0.H*c1+c1.H*c0
    unitary = s.eye(4)-hopping**2-s.I*hopping  # exp(-i*pi*hopping/2)
    assert unitary.H*unitary == s.eye(4)
    sigma = s.eye(4)/4
    insertion = c1.H*sigma+sigma*c1.H
    assert s.trace(c0*unitary*insertion*unitary.H) == -s.I
    assert s.trace(c1*unitary*insertion*unitary.H) == 0


@pytest.mark.parametrize("eta", [0, -1, F(-1, 2)])
def test_reference_witness_requires_faithfulness(eta):
    with pytest.raises(ValueError):
        charged.reference_extension(eta)
