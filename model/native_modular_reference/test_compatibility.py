"""Exact all-reference compatibility, independent algebra and hostile controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations

import pytest

from . import compatibility as c
from . import common_caps, regional


def diagonal(values):
    return {(i, i): value for i, value in enumerate(values)}


def test_native_exact_constraint_certificate():
    result = c.exact_certificate()
    assert result['shell_dimension'] == 66
    assert result['off_diagonal_ordered_entries_forced_zero'] == 4290
    assert result['spectator_equations_checked'] == 594
    assert result['diagonal_exact_rank'] == 54
    assert result['hermitian_real_constraints'] == 4344
    assert result['hermitian_generator_real_dimension'] == 12
    assert result['number_operator_intersection_dimensions'] == [2]*12
    assert result['all_faithful_references_covered']
    assert not result['nonfaithful_references_classified']
    assert not result['numerical_rank_used']


@pytest.mark.parametrize('vertices', range(3, 10))
def test_exact_rank_and_three_port_exception(vertices):
    result = c.exact_certificate(vertices)
    assert result['diagonal_exact_rank'] == vertices*(vertices-1)//2-vertices
    assert result['number_operator_intersection_dimensions'] == ([3]*3 if vertices == 3 else [2]*vertices)


def test_every_native_offdiagonal_entry_has_actual_forbidden_support():
    states = c.basis()
    witnesses = c.off_diagonal_witnesses()
    assert {(w.row, w.column) for w in witnesses} == {
        (i, j) for i in range(66) for j in range(66) if i != j}
    for w in witnesses:
        supports = c.bare_basis(12, (w.number_port, w.algebra_partner))
        assert all((w.row, w.column) not in support for support in supports)
        expected = int(w.number_port in states[w.column])-int(w.number_port in states[w.row])
        assert expected == w.commutator_coefficient != 0


def test_bare_supports_match_existing_native_regional_algebras():
    states = regional.source()[0].states
    canonical = {pair: i for i, pair in enumerate(c.basis())}
    reindex = {i: canonical[tuple(p for p, occupied in enumerate(state) if occupied)]
               for i, state in enumerate(states)}
    for ports in combinations(range(12), 2):
        independent = {frozenset((reindex[i], reindex[j]) for i, j in support)
                       for support in regional.basis_supports(ports, completed=False)}
        assert independent == {frozenset(support) for support in c.bare_basis(12, ports)}


@pytest.mark.parametrize('vertices', [3, 4, 5, 12])
def test_additive_generators_normalize_actual_algebras(vertices):
    onsite = [F((-1)**i*(i+1), i+2) for i in range(vertices)]
    operator = diagonal([onsite[i]+onsite[j] for i, j in c.basis(vertices)])
    before = deepcopy(operator)
    result = c.analyze_generator(vertices, operator)
    assert result['compatible']
    assert list(map(F, result['onsite_energies'])) == onsite
    assert c.directly_normalizes(vertices, operator)
    assert operator == before


@pytest.mark.parametrize('vertices', [3, 4, 12])
@pytest.mark.parametrize('imaginary', [False, True])
def test_nondiagonal_faithful_gibbs_references_are_excluded(vertices, imaginary):
    # exp(-H)/Tr(exp(-H)) is faithful for either exact Hermitian H.
    d = len(c.basis(vertices))
    operator = diagonal([F(i+1, 7) for i in range(d)])
    operator[0, d-1] = (0, F(2, 3)) if imaginary else F(2, 3)
    operator[d-1, 0] = (0, -F(2, 3)) if imaginary else F(2, 3)
    assert not c.analyze_generator(vertices, operator)['compatible']
    assert not c.directly_normalizes(vertices, operator)


def test_all_three_port_diagonal_references_are_allowed():
    assert c.directly_normalizes(3, diagonal([F(2, 7), -3, 11]))
    assert c.analyze_generator(3, diagonal([F(2, 7), -3, 11]))['compatible']
    weights = dict(zip(c.basis(3), (F(2, 7), F(3, 11), F(5, 13)), strict=True))
    assert common_caps.factor_characterization(3, weights)['all_bare_two_port_algebras_invariant']


def test_diagonal_obstruction_and_existing_product_weight_test_agree():
    vertices = 6
    energies = [F(i-2) for i in range(vertices)]
    values = [energies[i]+energies[j] for i, j in c.basis(vertices)]
    weights = {pair: F(2)**(-int(value)) for pair, value in zip(c.basis(vertices), values, strict=True)}
    assert c.analyze_generator(vertices, diagonal(values))['compatible']
    assert common_caps.factor_characterization(vertices, weights)['all_bare_two_port_algebras_invariant']
    values[0] += 1
    weights[c.basis(vertices)[0]] /= 2
    result = c.analyze_generator(vertices, diagonal(values))
    assert result['obstruction'] == 'spectator-dependent modular difference'
    assert not c.directly_normalizes(vertices, diagonal(values))
    assert not common_caps.factor_characterization(vertices, weights)['all_bare_two_port_algebras_invariant']


def test_compatible_generator_need_not_be_scalar():
    operator = diagonal([F(i+j) for i, j in c.basis(5)])
    assert len(set(operator.values())) > 1
    assert c.directly_normalizes(5, operator)
    # Compatibility does not assert stationarity under connected symmetric repair.


@pytest.mark.parametrize('vertices', [True, False, 2, 17, 3., '3'])
def test_vertex_domain_is_original_and_bounded(vertices):
    with pytest.raises(ValueError):
        c.exact_certificate(vertices)


@pytest.mark.parametrize('operator', [
    [], {(True, 0): 1}, {(0, 3): 1}, {(0, 0): True}, {(0, 0): .5},
    {(0, 0): float('nan')}, {(0, 0): complex(1)}, {(0, 0): (1, False)},
    {(0, 0): (1, 1)}, {(0, 1): 1}, {(0, 1): (0, 1), (1, 0): (0, 1)},
])
def test_hostile_exact_operator_inputs_fail_closed(operator):
    with pytest.raises(ValueError):
        c.analyze_generator(3, operator)
    with pytest.raises(ValueError):
        c.directly_normalizes(3, operator)


def test_returned_certificate_cannot_mutate_cached_constraints():
    first = c.exact_certificate(4)
    expected = deepcopy(first)
    first['diagonal_pivot_columns'].append(999)
    first['number_operator_intersection_dimensions'][0] = 999
    assert c.exact_certificate(4) == expected


def test_exact_sparse_algebra_membership_rejects_hidden_outside_dependence():
    support = c.bare_basis(5, (0, 1))[1]
    complete = {address: F(2, 3) for address in support}
    assert c.in_bare_algebra(5, (0, 1), complete)
    complete[support[0]] += F(1, 11)
    assert not c.in_bare_algebra(5, (0, 1), complete)
