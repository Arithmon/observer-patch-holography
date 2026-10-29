"""Independent finite-matrix and source controls for CP05 classifications."""
from fractions import Fraction
from itertools import product
import json

import pytest
import sympy as s

from . import classification as c
from model.local_repair_family import local_bias, modular_clock


@pytest.fixture(scope='module')
def receipt():
    return c.report()


def test_complete_current_trace_tensor_and_unbounded_integer_kernel(receipt):
    result = receipt['anomaly_classification']
    equations = result['all_trace_equations']
    assert len(equations) == 376
    assert sum(row['kind'] == 'symmetric_gauge' for row in equations) == 364
    matrix = s.Matrix([row['coefficients'] for row in equations])
    assert matrix.rank() == 4
    assert matrix.nullspace() == [s.ones(5, 1)]
    assert result['binary_submodule_controls']['admitted'] == [[0]*5, [1]*5]
    assert len(result['binary_submodule_controls']['excluded_with_witness']) == 30
    for row in result['binary_submodule_controls']['excluded_with_witness']:
        actual = (matrix*s.Matrix(row['retained_blocks']))[row['equation_index']]
        assert actual == s.Rational(row['nonzero_residual']) != 0


def test_occupation_charge_oracle_without_particle_names(receipt):
    rows = receipt['anomaly_classification']['multiplet_order']
    assert [row['dimension'] for row in rows] == [1, 6, 3, 3, 2]
    for row in rows:
        assert int(row['primitive_q6Y']) == -2*row['color_occupations']+3*row['weak_occupations']
        assert all(mask.bit_count() in (2, 4) for mask in row['source_masks'])
    # Independent standard invariant contractions from the occupation counts.
    charges = [s.Rational(row['primitive_q6Y']) for row in rows]
    dimensions = [row['dimension'] for row in rows]
    grav = [q*d for q, d in zip(charges, dimensions)]
    cubic = [q**3*d for q, d in zip(charges, dimensions)]
    color_mixed = [q*d*s.Rational(4, 3) if row['color_occupations'] in (1, 2) else 0
                   for q, d, row in zip(charges, dimensions, rows)]
    weak_mixed = [q*d*s.Rational(3, 4) if row['weak_occupations'] == 1 else 0
                  for q, d, row in zip(charges, dimensions, rows)]
    color_cubic = [0, 2, -1, -1, 0]
    independent = s.Matrix([grav, cubic, color_mixed, weak_mixed, color_cubic])
    assert independent.rank() == 4
    assert independent.nullspace() == [s.ones(5, 1)]


def test_arbitrary_copies_and_opposite_chiralities():
    assert c.validate_multiplicities((10**40,)*5)['complete_signed_generation_count'] == 10**40
    assert c.validate_multiplicities((2, 3, 4, 5, 6), (1, 2, 3, 4, 5))['all_local_traces_zero']
    assert c.validate_multiplicities((3, 1, 4, 1, 5), (3, 1, 4, 1, 5))['all_local_traces_zero']
    for i in range(5):
        values = [1]*5; values[i] = 0
        assert not c.validate_multiplicities(tuple(values))['all_local_traces_zero']
    assert c.validate_multiplicities((3,)*5)['original_SU2_doublet_parity'] == 0
    assert c.validate_multiplicities((0, 0, 0, 0, 1))['original_SU2_doublet_parity'] == 1


@pytest.mark.parametrize('bad', [(True,)*5, (1.0,)*5, (Fraction(1),)*5, (-1,)*5, (1,)*4, [1]*5])
def test_multiplicity_domain_rejects_lossy_or_negative_values(bad):
    with pytest.raises(ValueError):
        c.validate_multiplicities(bad)


def test_fixed_reference_group_and_scheduler_controls(receipt):
    result = receipt['clock_family_symmetry']
    p = modular_clock.parameters()
    assert p.c_weight == Fraction(1, 2)
    assert len(c.fixed_reference_stabilizer(local_bias.Parameters())) == 60
    assert len(c.fixed_reference_stabilizer(p)) == 2
    g = tuple(result['fixed_source_proper_stabilizer'][1])
    edge = next(e for e in local_bias.SEAMS if tuple(sorted(g[i] for i in e)) != e)
    weights = tuple(2 if e == edge else 1 for e in local_bias.SEAMS)
    changed = local_bias.Parameters(p.port_weights, weights, p.c_weight)
    assert len(c.fixed_reference_stabilizer(changed)) == 1
    with pytest.raises(ValueError):
        c.branch_covariance(changed, g)
    assert result['complete_C_E_branch_columns_covariant_with_record_recharting'] == 7920


def test_source_reference_breaks_nonresidual_rotation():
    p = modular_clock.parameters()
    states = local_bias.basis()
    index = {state: j for j, state in enumerate(states)}
    reference = local_bias.reference(p)
    admitted = {g for g, _ in c.fixed_reference_stabilizer(p)}
    for g, _ in c.rotation_representation():
        preserved = all(reference[index[tuple(state[g.index(j)] for j in range(12))]] == reference[i]
                        for i, state in enumerate(states))
        assert preserved == (g in admitted)


def test_commutant_and_positive_mass_controls(receipt):
    result = receipt['clock_family_symmetry']
    r = s.Matrix(result['residual_family_involution'])
    matrix = s.Matrix(3, 3, s.symbols('a:9'))
    solved = s.linsolve(list(r*matrix-matrix*r), list(matrix))
    pattern = next(iter(solved))
    assert [pattern[i] for i in (1, 3, 5, 7)] == [0]*4
    control = result['distinct_mass_and_two_family_mixing_control']
    yu, yd = s.Matrix(control['Y_u']), s.Matrix(control['Y_d'])
    assert yu.is_positive_definite and yd.is_positive_definite
    assert sorted((yu*yu.H).eigenvals()) == [1, 4, 9]
    assert sorted((yd*yd.H).eigenvals()) == [1, 9, 25]
    hu, hd = yu*yu.H, yd*yd.H
    bracket = hu*hd-hd*hu
    assert bracket != s.zeros(3)
    assert bracket.det() == s.trace(bracket**3) == 0


def test_complex_yukawas_and_invariant_noncanonical_kinetic_metrics():
    r = s.diag(-1, 1, -1)
    yu = s.Matrix([[1, 0, s.I], [0, 3, 0], [2, 0, 4]])
    yd = s.Matrix([[2, 0, 1+s.I], [0, 2, 0], [1, 0, 3-s.I]])
    left = s.Matrix([[2, 0, s.I], [0, 4, 0], [-s.I, 0, 3]])
    right_u, right_d = s.diag(2, 5, 3), s.diag(7, 2, 4)
    assert all(x*r == r*x for x in (yu, yd, left, right_u, right_d))
    assert left.is_positive_definite
    # Generalized mass maps in a common left kinetic metric. The common
    # isolated family is unaffected by a different right kinetic metric.
    hu = left.inv()*yu*right_u.inv()*yu.H
    hd = left.inv()*yd*right_d.inv()*yd.H
    bracket = hu*hd-hd*hu
    assert bracket[1, :] == s.zeros(1, 3)
    assert bracket[:, 1] == s.zeros(3, 1)
    assert bracket.det() == s.trace(bracket**3) == 0


def test_supplied_symmetry_breaking_can_escape(receipt):
    row = receipt['clock_family_symmetry']['symmetry_breaking_control']
    yd = s.Matrix(row['Y_d'])
    assert yd.is_positive_definite
    r = s.diag(-1, 1, -1)
    assert yd*r != r*yd
    hu, hd = s.diag(1, 4, 9), yd*yd.H
    bracket = hu*hd-hd*hu
    assert s.expand(bracket.det()) == 46080*s.I
    assert s.expand(s.trace(bracket**3)) == 138240*s.I


def test_retained_receipt_matches_sources_and_no_gate_promotion(receipt):
    assert json.loads(c.OUTPUT.read_text()) == receipt
    assert not receipt['scope']['integrated_gate_accepted']
    assert receipt['scope']['physical_comparison_count_added'] == 0
    assert not receipt['clock_family_symmetry']['family_identification_or_Yukawa_selection_derived']


def test_mutated_frozen_parameter_slots_are_revalidated():
    bad = modular_clock.parameters()
    object.__setattr__(bad, 'port_weights', (0,)*12)
    with pytest.raises(ValueError):
        c.fixed_reference_stabilizer(bad)
    other = modular_clock.parameters()
    object.__setattr__(other, 'repair_probability', Fraction(-1, 10**400))
    with pytest.raises(ValueError):
        c.branch_covariance(other, tuple(range(12)))


def test_actual_stationary_reference_supplies_three_distinct_response_axes(receipt):
    row = receipt['clock_family_symmetry']['derived_source_anisotropy']
    occupation = [s.Rational(x) for x in row['stationary_one_port_occupations']]
    assert occupation == [s.Rational(525, 1501)]*2 + [s.Rational(976, 7505)]*10
    assert sum(occupation) == 2
    # Independent endpoint-pair partition function and two-class mean.
    z = 25**2 + 2*10*25*8 + 45*8**2
    assert s.Rational(25*(25+10*8), z) == occupation[0]
    assert s.Rational(8*(2*25+9*8), z) == occupation[2]
    a = s.Matrix(row['matrix'])
    assert len(a.eigenvals()) == 3
    assert all(x.is_positive for x in a.diagonal())
    assert s.trace(a) == s.Rational(1, 2)
    assert a*a**2 == a**2*a
    assert a*s.diag(-1, 1, -1) == s.diag(-1, 1, -1)*a
