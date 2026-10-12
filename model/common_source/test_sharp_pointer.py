"""Independent source-matrix controls and hostile eigensolver candidates."""
from fractions import Fraction as F
import math
from pathlib import Path
import sys
from types import ModuleType

import mpmath as mp
import numpy as np
import pytest

from . import sharp_pointer as pointer


@pytest.mark.parametrize('d', (2, 4, 16, 32))
def test_projector_action_matches_independent_normalized_fock_eigensolve(d):
    coordinate = np.diag(np.sqrt(np.arange(1, d)), 1)
    coordinate += coordinate.T
    eigenvalues, eigenvectors = np.linalg.eigh(coordinate)
    exact = (eigenvectors[:, eigenvalues > 0] @ eigenvectors[:, eigenvalues > 0].T)
    # Different spectator/gauge labels are independent reference components.
    source = {(0, 0, 0): (F(2), F(-3)), (1, 0, 0): (F(1, 7), F(5, 4)),
              (0, 1, 3): (F(-4, 3), F(2, 5))}
    if d >= 4:
        source[3, 0, 0] = F(1, 9), F(-1, 13)
        source[2, 1, 3] = F(-1, 7), F(3, 8)
    for outcome in (0, 1):
        actual = pointer.apply_pointer(source, d, outcome)
        operator = exact if outcome else np.eye(d)-exact
        for spectator, gauge in ((0, 0), (1, 3)):
            incoming = np.array([complex(*map(float, source.get((n, spectator, gauge), (0, 0))))
                                 * math.sqrt(math.factorial(n)) for n in range(d)])
            outgoing = np.array([complex(*map(float, actual.get((n, spectator, gauge), (0, 0))))
                                 * math.sqrt(math.factorial(n)) for n in range(d)])
            np.testing.assert_allclose(outgoing, operator @ incoming, atol=2e-14, rtol=2e-13)
    eta = F(pointer.pointer_certificate(d)['positive_projector_operator_error_upper'])
    assert 0 < eta < F(1, 10**25)


def test_two_level_exact_projector_is_inside_reported_operator_enclosure():
    # Independent exact source sign: (I+X)/2. Do not compare only probabilities.
    cert = pointer.pointer_certificate(2)
    eta = F(cert['positive_projector_operator_error_upper'])
    entries = []
    for column in range(2):
        result = pointer.apply_pointer({(column, 0, 0): (F(1), F(0))}, 2, 1)
        entries.extend((result[row, 0, 0][0]-F(1, 2))**2 for row in range(2))
    assert sum(entries, F(0)) <= eta**2


def test_each_residual_matches_independent_high_precision_original_source():
    basis, scale, diagonal, value_scale = pointer._seed(4)
    certificate = pointer.certify_eigenbasis(basis, scale, diagonal, value_scale)
    with mp.workdps(100):
        source = mp.matrix(4)
        for n in range(1, 4):
            source[n-1, n] = source[n, n-1] = mp.sqrt(n)
        columns = mp.matrix([[mp.mpf(basis[n][j])/scale*mp.sqrt(math.factorial(n))
                              for j in range(4)] for n in range(4)])
        values = mp.diag([mp.mpf(x)/value_scale for x in diagonal])
        for key, residual in (
            ('basis_orthogonality_frobenius_error_upper', columns.T*columns-mp.eye(4)),
            ('source_eigenvector_residual_frobenius_upper', source*columns-columns*values),
        ):
            bound = F(certificate[key])
            reported = mp.mpf(bound.numerator)/bound.denominator
            independent = mp.sqrt(sum(abs(x)**2 for x in residual))
            assert independent <= reported+mp.mpf('1e-95')
            assert reported-independent < mp.mpf('1e-59')


@pytest.mark.parametrize('mutation', ('zero', 'wrong_spectrum', 'wrong_vectors', 'wrong_metric'))
def test_independent_certificate_rejects_plausible_bad_eigensolver_outputs(mutation):
    basis, scale, diagonal, value_scale = pointer._seed(4)
    if mutation == 'zero':
        basis = ((0,)*4,)*4
    elif mutation == 'wrong_spectrum':
        diagonal = tuple(-x for x in diagonal)
    elif mutation == 'wrong_vectors':
        basis = tuple(tuple(reversed(row)) for row in basis)
    else:
        # Treating normalized Fock components as Bargmann coefficients loses
        # the factorial metric. Ordinary Euclidean orthogonality must not pass.
        basis = tuple(tuple(x*(2 if n == 3 else 1) for x in row)
                      for n, row in enumerate(basis))
    with pytest.raises(ArithmeticError):
        pointer.certify_eigenbasis(basis, scale, diagonal, value_scale)


def test_valid_column_permutation_and_sign_changes_preserve_certification():
    basis, scale, diagonal, value_scale = pointer._seed(4)
    changed = tuple(tuple(-x for x in reversed(row)) for row in basis)
    first = pointer.certify_eigenbasis(basis, scale, diagonal, value_scale)
    second = pointer.certify_eigenbasis(changed, scale, tuple(reversed(diagonal)), value_scale)
    assert first == second


def test_seed_uses_and_restores_its_own_precision():
    before = mp.mp.dps
    try:
        mp.mp.dps = 9
        low = pointer._seed.__wrapped__(2)
        assert mp.mp.dps == 9
        mp.mp.dps = 83
        high = pointer._seed.__wrapped__(2)
        assert mp.mp.dps == 83
        assert low == high
    finally:
        mp.mp.dps = before


@pytest.mark.parametrize('d', (True, 2., 0, -2, 3, 130))
def test_rejects_bad_dimension_without_eigensolve(d):
    with pytest.raises(ValueError):
        pointer.pointer_certificate(d)


def test_rejects_bad_basis_and_out_of_code_input_and_preserves_cached_evidence():
    basis, scale, diagonal, value_scale = pointer._seed(2)
    with pytest.raises(ValueError):
        pointer.certify_eigenbasis(basis, True, diagonal, value_scale)
    with pytest.raises(ValueError):
        pointer.certify_eigenbasis(basis, scale, (False, diagonal[1]), value_scale)
    with pytest.raises(ValueError):
        pointer.apply_pointer({(0, 2, 0): (F(1), F(0))}, 2, 1)
    with pytest.raises(ValueError):
        pointer.apply_pointer({}, 2, True)
    changed = pointer.pointer_certificate(2)
    changed['positive_projector_operator_error_upper'] = '0'
    assert F(pointer.pointer_certificate(2)['positive_projector_operator_error_upper']) > 0


def test_complete_comparison_preserves_initial_norm_all_flags_and_both_error_terms(monkeypatch):
    from . import occupation
    # Two raw complex vacuum branches have mass 1/2 after division by the
    # one initial norm 4. Use a closed distance formula, not the producer.
    history = {'initial_polynomial_norm_squared': '4', 'configurations': [{
        'feedback_enabled': True,
        'complete_output_half_trace_norm_error_upper': '1/100',
        'retained_first_records': [
            {'first_record': a, 'pre_readout_polynomial': [[0, 0, a, '1', '1']]}
            for a in (0, 1)],
    }]}
    # Isolate budget composition; the live receipt separately evaluates the
    # unmodified interacting occupation theorem at dimensions 32,64,128.
    monkeypatch.setattr(occupation, 'occupation_cutoff_bound',
                        lambda d: {'pre_final_readout_trace_distance_upper': '1/50'})
    row = pointer.complete_cutoff_comparison(history, 2)['configurations'][0]
    assert [(r['first_record'], r['sign_record']) for r in row['retained_record_comparisons']] == [
        (0, 0), (0, 1), (1, 0), (1, 1)]
    with mp.workdps(70):
        overlap = mp.mpf(1)/4+1/(2*mp.sqrt(2*mp.pi))
        distance = mp.sqrt(1-4*overlap**2)
        def real(value):
            value = F(value)
            return mp.mpf(value.numerator)/value.denominator
        lo, hi = map(real, row['polynomial_instrument_distance_interval'])
        assert lo <= distance <= hi
        assert hi-lo < mp.mpf('2e-12')
        lo, hi = map(real, row['complete_finite_to_infinite_trace_distance_interval'])
        assert abs(lo-(distance-mp.mpf('0.04'))) < mp.mpf('2e-12')
        assert abs(hi-(distance+mp.mpf('0.04'))) < mp.mpf('2e-12')
    eta = F(pointer.pointer_certificate(2)['positive_projector_operator_error_upper'])
    assert F(row['pointer_seed_density_error_upper']) == 2*eta+eta**2


@pytest.mark.parametrize('d', (2, 32))
def test_valid_rotated_candidate_bound_contains_exact_parity_witness(d):
    # Rotate only the smallest positive/negative candidate eigenvectors by
    # the exact rational orthogonal matrix with cos=9999/10001, sin=200/10001.
    # This is a valid, deliberately imperfect seed, rather than a candidate
    # that should be rejected by a malformed-input or usefulness guard.
    basis, scale, diagonal, value_scale = pointer._seed(d)
    negative = max((j for j in range(d) if diagonal[j] < 0), key=lambda j: diagonal[j])
    positive = min((j for j in range(d) if diagonal[j] > 0), key=lambda j: diagonal[j])
    changed = []
    for row in basis:
        output = [10001*x for x in row]
        output[negative] = 9999*row[negative]+200*row[positive]
        output[positive] = -200*row[negative]+9999*row[positive]
        changed.append(tuple(output))
    changed, changed_scale = tuple(changed), 10001*scale
    certificate = pointer.certify_eigenbasis(changed, changed_scale, diagonal, value_scale)

    # Independent exact source oracle: R=(-1)^N anticommutes with a+a*,
    # whose even-dimensional Jacobi matrix is invertible. Thus R B_d R=I-B_d
    # and every parity-pure vector has B_d expectation exactly 1/2. Neither
    # the eigensolver nor its residual certificate supplies this expectation.
    # Use the even part of a seed column merely to select a rational vector.
    vector = tuple(basis[n][positive] if n % 2 == 0 else 0 for n in range(d))
    weights = tuple(math.factorial(n) for n in range(d))
    norm = sum(weight*x*x for weight, x in zip(weights, vector, strict=True))
    assert norm > 0
    expectation_numerator = sum(
        sum(weights[n]*vector[n]*changed[n][j] for n in range(d))**2
        for j in range(d) if diagonal[j] > 0)
    expectation = F(expectation_numerator, changed_scale**2*norm)
    exact_error_lower = abs(expectation-F(1, 2))
    reported = F(certificate['positive_projector_operator_error_upper'])
    assert F(1, 100) < exact_error_lower < F(1, 40)
    assert 0 < reported < F(1, 10)
    assert exact_error_lower <= reported, 'certificate understates the exact parity-witness error'


@pytest.mark.parametrize('mutation,d', (('omit_gap', 32), ('quarter_bound', 2)))
def test_exact_parity_witness_kills_understated_projector_bounds(monkeypatch, mutation, d):
    code = Path(pointer.__file__).read_text(encoding='utf-8')
    original = 'projection_error = operator_residual/gap+epsilon*(2+epsilon)'
    assert code.count(original) == 1
    replacement = ('projection_error = operator_residual+epsilon*(2+epsilon)'
                   if mutation == 'omit_gap' else
                   'projection_error = (operator_residual/gap+epsilon*(2+epsilon))/4')
    mutant = ModuleType('model.common_source._understated_pointer')
    mutant.__package__ = 'model.common_source'
    mutant.__file__ = pointer.__file__
    exec(compile(code.replace(original, replacement), pointer.__file__, 'exec'), mutant.__dict__)
    monkeypatch.setattr(sys.modules[__name__], 'pointer', mutant)
    with pytest.raises(AssertionError, match='understates the exact parity-witness error'):
        test_valid_rotated_candidate_bound_contains_exact_parity_witness(d)
