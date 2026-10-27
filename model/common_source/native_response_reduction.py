"""Exact limits on reducing the native response by channels or quotients.

The canonical source algebra is su(3)+su(3)+iR I6. Its scalar center is
invisible under conjugation. The specified current interface instead has
u(3)+so(3), with a visible relative phase. These are different algebras.

The channel statement is deliberately narrow. A unital completely positive
Heisenberg map fixing the right block projector and all three spin-one
generators fixes every operator on that block. Their Schwarz defects are
positive and sum to zero because sum J_i**2=2I; hence each generator belongs
to the multiplicative domain. Since E_ji=delta_ij I-J_i J_j, every matrix
unit is fixed. This covers a discarded-record recovery channel satisfying
these hypotheses, not arbitrary instruments, approximate protection, or a
new dynamics which disturbs the retained rotations.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
from pathlib import Path

import sympy as s


def _rational(value):
    if type(value) is int:
        return s.Rational(value)
    if type(value) is Fraction:
        return s.Rational(value.numerator, value.denominator)
    if isinstance(value, s.Rational) and not isinstance(value, s.Float):
        return value
    raise ValueError('an original exact rational coefficient is required')


def _unit(i, j, dimension=3):
    value = s.zeros(dimension)
    value[i, j] = 1
    return s.ImmutableMatrix(value)


def spin_generators():
    """Hermitian spin-one generators in an orthonormal real geometric frame."""
    return tuple(s.ImmutableMatrix(s.I*value) for value in (
        s.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]]),
        s.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]]),
        s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])))


def _multipole_action(a, b, matrix):
    # Complex linear transpose, not the anti-linear operation of conjugation.
    return s.ImmutableMatrix((a+b)*matrix/2+(b-a)*matrix.T/2
                             +(1-b)*s.trace(matrix)*s.eye(3)/3)


def multipole_choi(a, b):
    """Unnormalized Choi matrix, including non-CP diagnostic controls.

    Under the actual real three-dimensional A5 representation, End(C3)
    decomposes as 1+3+5 with multiplicity one. A covariant unital Hermiticity
    preserving map therefore has real multipliers (1,a,b). Complete
    positivity is a separate condition; this function does not assume it.
    """
    a, b = _rational(a), _rational(b)
    return s.ImmutableMatrix(sum((s.kronecker_product(_unit(i, j),
        _multipole_action(a, b, _unit(i, j))) for i in range(3) for j in range(3)), s.zeros(9)))


def multipole_eigenvalues(a, b):
    """Exact Choi eigenvalues on scalar, symmetric-traceless and skew sectors."""
    a, b = _rational(a), _rational(b)
    return ((s.factor((1+3*a+5*b)/3), 1),
            (s.factor((2-3*a+b)/6), 5),
            (s.factor((2+3*a-5*b)/6), 3))


def _choi_projectors():
    vector = s.eye(3).reshape(9, 1)
    scalar = vector*vector.T/3
    swap = sum((s.kronecker_product(_unit(i, j), _unit(j, i))
                for i in range(3) for j in range(3)), s.zeros(9))
    return tuple(s.ImmutableMatrix(x) for x in
                 (scalar, (s.eye(9)+swap)/2-scalar, (s.eye(9)-swap)/2))


def covariant_channel(a, b):
    """Return the exact row-vectorized qutrit channel, only when it is CPTP."""
    a, b = _rational(a), _rational(b)
    if any(value < 0 for value, _ in multipole_eigenvalues(a, b)):
        raise ValueError('the proposed rotation/quadrupole multipliers are not completely positive')
    return s.ImmutableMatrix.hstack(*(_multipole_action(a, b, _unit(i, j)).reshape(9, 1)
                                      for i in range(3) for j in range(3)))


def _su3():
    skew = tuple(_unit(i, j)-_unit(j, i) for i in range(3) for j in range(i+1, 3))
    symmetric = tuple(s.I*(_unit(i, j)+_unit(j, i)) for i in range(3) for j in range(i+1, 3))
    diagonal = (s.I*s.diag(1, -1, 0), s.I*s.diag(0, 1, -1))
    return tuple(s.ImmutableMatrix(x) for x in (*skew, *symmetric, *diagonal))


def canonical_bases():
    """Real anti-Hermitian bases; no complex-rank substitution is made."""
    left = tuple(s.ImmutableMatrix(s.diag(x, s.zeros(3))) for x in _su3())
    right = tuple(s.ImmutableMatrix(s.diag(s.zeros(3), x)) for x in _su3())
    rotations = tuple(s.ImmutableMatrix(s.diag(s.zeros(3), -s.I*j)) for j in spin_generators())
    scalar = s.ImmutableMatrix(s.I*s.eye(6))
    relative = s.ImmutableMatrix(s.diag(s.I, s.I, s.I, 0, 0, 0))
    return {'source': (*left, *right, scalar), 'target': (*left, *rotations, relative),
            'scalar_center': scalar, 'target_center': relative}


def _real_rank(matrices):
    return s.Matrix.hstack(*(s.Matrix([s.re(x) for x in value]+[s.im(x) for x in value])
                             for value in matrices)).to_DM(extension=True).rank()


def _native_covariance_certificate():
    """Prove multiplicity-free 1+3prime+5 using all actual port rotations."""
    from . import native_gauge_object as objects
    rows = objects.geometric_projectors()['rotation_rows']
    if len(rows) != 60:
        raise ArithmeticError('the complete native proper-rotation group is required')
    characters = []
    for row in rows:
        record = row['rotation3']
        rotation = s.zeros(3)
        if record['shape'] != [3, 3]:
            raise ArithmeticError('actual native rotation has the wrong carrier')
        for i, j, value in record['entries']:
            rotation[i, j] = s.sympify(value).xreplace({s.sqrt(5): -s.sqrt(5)})
        if (not objects._zero(rotation.T*rotation-s.eye(3))
                or s.simplify(rotation.det()) != 1):
            raise ArithmeticError('the right native triplet has lost its real orthogonal action')
        trace, square_trace = s.trace(rotation), s.trace(rotation*rotation)
        characters.append(s.Matrix([1, (trace**2-square_trace)/2,
                                     (trace**2+square_trace)/2-1]))
    gram = objects._clean(sum((row*row.T for row in characters), s.zeros(3))/60)
    if gram != s.eye(3):
        raise ArithmeticError('actual native scalar, rotation and quadrupole representations are not distinct irreducibles')
    return {'actual_proper_native_rotations_checked': 60,
            'right_EndC3_irreducible_dimensions': [1, 3, 5],
            'exact_character_inner_product_matrix': [[int(x) for x in row] for row in gram.tolist()],
            'every_A5_covariant_linear_map_is_scalar_on_these_three_summands': True}


def report():
    """Bind the canonical theorem to the owned actual-source comparison."""
    from . import native_gauge_object as objects, native_response_comparison as comparison
    actual = comparison.exact_data()
    canonical = canonical_bases()
    frame = actual['geometric_frame']
    if frame.shape != (36, 6) or not objects._zero(frame.conjugate().T*frame-s.eye(6)):
        raise ArithmeticError('the actual geometric source chart is not an isometry')
    right_frame = frame[:, 3:]
    if not objects._zero(right_frame*right_frame.conjugate().T-objects.geometric_projectors()['three_prime']):
        raise ArithmeticError('the right source block is not the native geometric three-prime sector')
    for key, expected, dimension in (('native_basis', canonical['source'], 17),
                                      ('declared_basis', canonical['target'], 12)):
        values = actual[key]
        if (len(values) != dimension or _real_rank(values) != dimension
                or _real_rank((*values, *expected)) != dimension):
            raise ArithmeticError('the actual source/current span differs from the canonical theorem domain')
    result = canonical_certificate()
    result['actual_source_binding'] = {
        'owning_comparison': 'native_response_comparison.exact_data',
        'complete_native_source_basis_dimension': 17,
        'complete_declared_current_basis_dimension': 12,
        'same_actual_orthonormal_geometric_frame': True,
        'right_block_is_actual_native_three_prime': True,
        'both_exact_real_spans_independently_checked': True,
        **_native_covariance_certificate()}
    return result


def source_bindings():
    from . import native_response_comparison as comparison
    root = Path(__file__).resolve().parents[2]
    bindings = comparison.source_bindings()
    for name in ('model/common_source/native_response_reduction.py',
                 'model/common_source/test_native_response_reduction.py'):
        bindings[name] = hashlib.sha256((root/name).read_bytes()).hexdigest()
    return bindings


def canonical_certificate():
    """Finite exact witnesses for the algebraic and physical-channel boundaries."""
    basis = canonical_bases()
    source, target = basis['source'], basis['target']
    source_rank, target_rank = _real_rank(source), _real_rank(target)
    intersection = source_rank+target_rank-_real_rank((*source, *target))
    observable = lambda values: _real_rank(tuple(x-s.trace(x)*s.eye(6)/6 for x in values))
    if (source_rank, target_rank, intersection, observable(source), observable(target)) != (17, 12, 11, 16, 12):
        raise ArithmeticError('source and target algebra/center distinction failed')
    probe = _unit(0, 3, 6)
    if (basis['scalar_center']*probe-probe*basis['scalar_center'] != s.zeros(6)
            or basis['target_center']*probe-probe*basis['target_center'] != s.I*probe):
        raise ArithmeticError('relative-center observation control failed')
    js = spin_generators()
    if sum((j*j for j in js), s.zeros(3)) != 2*s.eye(3):
        raise ArithmeticError('the right spin Casimir differs')
    for i in range(3):
        for j in range(3):
            if int(i == j)*s.eye(3)-js[i]*js[j] != _unit(j, i):
                raise ArithmeticError('spin products do not reconstruct every matrix unit')
    # The projector hypothesis matters for a map on the whole M6. This valid
    # Schrödinger channel replaces the left block by |B0><B0| and preserves
    # the right block. Its Heisenberg dual fixes right rotations, but moves
    # both the right projector and a quadrupole. It does not preserve u3_A.
    left, right = s.diag(1, 1, 1, 0, 0, 0), s.diag(0, 0, 0, 1, 1, 1)
    leak_kraus = (right, *(_unit(3, i, 6) for i in range(3)))
    if sum((k.conjugate().T*k for k in leak_kraus), s.zeros(6)) != s.eye(6):
        raise ArithmeticError('the scope countercontrol is not trace preserving')
    heisenberg = lambda x: sum((k.conjugate().T*x*k for k in leak_kraus), s.zeros(6))
    quadrupole = s.diag(0, 0, 0, 1, -1, 0)
    if (heisenberg(right) != s.eye(6) or heisenberg(quadrupole) != quadrupole+left
            or any(heisenberg(s.diag(s.zeros(3), j)) != s.diag(s.zeros(3), j) for j in js)):
        raise ArithmeticError('the full-carrier projector hypothesis countercontrol failed')
    # An explicit attempt to keep rotations but erase all five quadrupoles
    # is not even positive. The witness has unit trace and rank one.
    rho = s.Matrix([[1, -s.I, 0], [s.I, 1, 0], [0, 0, 0]])/2
    negative = s.Matrix([1, -s.I, 0])/s.sqrt(2)
    expectation = (negative.conjugate().T*_multipole_action(1, 0, rho)*negative)[0]
    if expectation != -s.Rational(1, 6):
        raise ArithmeticError('false multipole projection lost its negative-state witness')
    # A physical quintet-erasing alternative necessarily damages the rotations.
    # This boundary point is attained by Kraus (I,Jx,Jy,Jz)/sqrt(3).
    kraus = tuple(s.ImmutableMatrix(k/s.sqrt(3)) for k in (s.eye(3), *js))
    if sum((k.conjugate().T*k for k in kraus), s.zeros(3)) != s.eye(3):
        raise ArithmeticError('the physical erasure control is not trace preserving')
    for i in range(3):
        for j in range(3):
            x = _unit(i, j)
            actual = sum((k*x*k.conjugate().T for k in kraus), s.zeros(3))
            if actual != _multipole_action(s.Rational(2, 3), 0, x):
                raise ArithmeticError('the attained CP tradeoff differs on a matrix unit')
    covariant_channel(s.Rational(2, 3), 0)
    projectors = _choi_projectors()
    if (sum(projectors, s.zeros(9)) != s.eye(9)
            or tuple(s.trace(p) for p in projectors) != (1, 5, 3)
            or any(p*q != (p if i == j else s.zeros(9))
                   for i, p in enumerate(projectors) for j, q in enumerate(projectors))):
        raise ArithmeticError('the complete Choi eigenspaces differ')
    # Both sides depend affinely on (a,b); three affine basis points prove
    # the eigenvalue formula for every real coefficient, not just a sample.
    for a, b in ((0, 0), (1, 0), (0, 1)):
        predicted = sum((value*p for (value, _), p in
                         zip(multipole_eigenvalues(a, b), projectors, strict=True)), s.zeros(9))
        if multipole_choi(a, b) != predicted:
            raise ArithmeticError('the symbolic affine Choi decomposition failed')
    return {'schema': 'oph.native_response_reduction.v1',
        'canonical_algebras': {'source': 'su3_A + su3_B + iR I6',
            'source_real_dimension': source_rank, 'source_observable_dimension': observable(source),
            'target': 'u3_A + so3_B', 'target_real_dimension': target_rank,
            'target_observable_dimension': observable(target), 'intersection_real_dimension': intersection,
            'target_relative_center_in_source': False,
            'source_scalar_center_visible_on_EndM6': False,
            'target_center_cross_block_commutator': 'i E_03'},
        'exact_unitary_or_quotient': {
            'unitary_reencoding_preserves_dimension_and_scalar_center': True,
            'surjective_real_Lie_quotient_dimensions': [0, 1, 8, 9, 16, 17],
            'reason': 'each su3 summand is a real simple ideal; the scalar center has dimension one',
            'target_dimension12_is_exact_homomorphic_quotient': False},
        'unital_CP_fixed_observable_obstruction': {
            'hypotheses': 'fix the right block projector and all three Hermitian spin-one generators',
            'Casimir': 'sum_i J_i^2=2 I3',
            'Schwarz_defects': 'positive and sum to zero, hence each vanishes',
            'matrix_units': 'E_ji=delta_ij I3-J_i J_j',
            'fixed_right_operator_algebra_complex_dimension': 9,
            'right_quadrupoles_cannot_be_deleted_under_these_hypotheses': True,
            'arbitrary_untouched_reference_domain': True,
            'projector_hypothesis_countercontrol':
                'Kraus P_B and |B0><Ai| fix right J_i but map P_B to I6 and Q_B to Q_B+P_A in the Heisenberg picture',
            'countercontrol_preserves_left_u3_observables': False,
            'applies_to_general_irreversible_or_postselected_dynamics': False},
        'invalid_positive_or_CP_controls': {
            'rotation_preserving_quintet_erasure_negative_state_expectation': str(expectation),
            'rotation_preserving_quintet_erasure_Choi_minimum': '-1/6',
            'transpose_Choi_minimum': '-1',
            'transpose_average_Choi_minimum': '-1/2',
            'entrywise_real_part_is_complex_linear': False},
        'covariant_CP_tradeoff': {
            'domain': 'one qutrit, unital A5-covariant complex-linear Hermiticity-preserving channel',
            'multipliers': 'scalar1, rotation a, symmetric-traceless quadrupole b',
            'Choi_eigenvalues_and_multiplicities': [
                ['(1+3a+5b)/3', 1], ['(2-3a+b)/6', 5], ['(2+3a-5b)/6', 3]],
            'affine_Choi_identity_checked_on_a_complete_parameter_basis': True,
            'exact_CP_conditions': ['1+3a+5b>=0', '2-3a+b>=0', '2+3a-5b>=0'],
            'a_equals_one_forces_b_equals_one': True,
            'b_equals_zero_maximum_a': '2/3',
            'maximum_attained_by_Kraus': '(I3,Jx,Jy,Jz)/sqrt(3)',
            'SO3_twirl_also_erases_rotations': True,
            'new_physical_selection_mechanism_derived': False}}
