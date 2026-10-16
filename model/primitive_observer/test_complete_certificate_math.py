"""Independent composition and exact phase-loss checks of the new certificate.

The composition check is numerical falsification on synthetic complete output
arrays, including a mixed entangled reference and non-Hermitian work blocks.
Its dense ideal square roots are not inputs to the exact bound. The phase
check is exact rational algebra on the independently enumerated native maps.
Neither test substitutes synthetic arrays for a source execution receipt.
"""
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from ..common_source import quantum_native as native_source
from . import complete_certificate as certificate
from . import exact_native


BITS = tuple(product((0, 1), repeat=4))
FIRST = {
    'C': {(0, 1), (1, 0), (0, 0), (1, 2)},
    'E': {(0, 5), (1, 4), (0, 4), (1, 6)},
}


def _phase_distance(matrix):
    """Exact half norm of M X M* when its two columns are orthogonal.

    X is the difference between normalized logical plus/minus densities.
    Disjoint column support gives eigenvalues +/-sqrt(n0*n1). This does
    not call the certificate's effect or completeness implementation.
    """
    supports = [{row for row, column in matrix if column == c} for c in (0, 1)]
    assert supports[0].isdisjoint(supports[1])
    norms = [sum((weight for (_, column), weight in matrix.items() if column == c), F(0))
             for c in (0, 1)]
    squared = norms[0]*norms[1]
    numerator, denominator = math.isqrt(squared.numerator), math.isqrt(squared.denominator)
    assert numerator*numerator == squared.numerator
    assert denominator*denominator == squared.denominator
    return F(numerator, denominator)


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('mode', ('protected_qubit', 'passive_qubit'))
def test_exact_native_maps_preserve_or_erase_phase_on_the_selected_branch(law, mode):
    tree = exact_native.ideal_tree(law, mode)
    assert {key[:2] for key in tree} == FIRST[law]
    conditional = {}
    for prefix, selector in product(sorted(FIRST[law]), (0, 1)):
        conditional[prefix, selector] = sum((
            _phase_distance(matrix) for key, matrix in tree.items()
            if key[:3] == (*prefix, selector)), F(0))
        expected = F(1, 4) if mode == 'protected_qubit' or selector == 1 else F(0)
        assert conditional[prefix, selector] == expected

    # Use different exact work probabilities at each current. All later
    # work outcomes remain grouped by the actual third-bit seam selector.
    q = {-1: F(1, 4), 0: F(1, 2), 1: F(3, 4)}
    current = {0: 0, 1: 1, 2: -1}
    complete = sum((q[current[raw % 4]]*conditional[((coin, raw), 0)]
                    +(1-q[current[raw % 4]])*conditional[((coin, raw), 1)]
                    for coin, raw in FIRST[law]), F(0))
    qbar = (q[-1]+2*q[0]+q[1])/4
    assert complete == (1 if mode == 'protected_qubit' else 1-qbar)


def _dense(matrix):
    value = np.zeros((512, 2), complex)
    for entry, square in matrix.items():
        value[entry] = math.sqrt(float(square))
    return value


def _nuclear(value):
    return float(np.linalg.svd(value, compute_uv=False).sum())


def _work_blocks():
    state = np.array([[.75, .25], [.25, .25]], complex)
    values = {}
    errors = []
    for index, current in enumerate((-1, 0, 1), 1):
        q, xi = F(index, 4), F(index, 64)
        errors.append(xi/2)
        for bits in BITS:
            weight = q if bits == (0, 0, 0, 0) else 1-q if bits == (1, 0, 1, 1) else F(0)
            target = float(weight)*state
            actual = target.copy()
            actual[0, 1] += 1j*float(weight*xi)
            # Every error is exactly i*weight*xi*|0><1|, of nuclear norm
            # weight*xi. The entries are dyadic and introduce no rounding.
            values[current, bits] = actual, target
    return state, values, tuple(errors)


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('mode', ('protected_qubit', 'passive_qubit'))
def test_composition_bounds_complete_mixed_reference_and_nonhermitian_outputs(law, mode):
    tree = exact_native.ideal_tree(law, mode)
    menu = native_source._menu(law)
    assert len(menu) == (770 if law == 'C' else 2566)
    padding = next(row for row in range(512) if row not in native_source.CODES)
    nonzero_key = sorted(tree)[0]
    first_abort_key = next(record for record in menu if not record[2])
    final_abort_key = (*nonzero_key[:2], 0, 0, 3)

    # Mixed, exactly normalized, entangled logical/reference state:
    # (3/4)|Phi+><Phi+| + (1/4)|01><01|. All entries are exact dyadics.
    tau = np.zeros((4, 4), complex)
    tau[np.ix_([0, 3], [0, 3])] = 3/8
    tau[1, 1] = 1/4
    assert np.trace(tau) == 1
    assert np.linalg.eigvalsh(tau.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4))[0] < 0

    prepared, work, work_errors = _work_blocks()
    native_matrices, norm_cache = {}, {}
    reduced_rows = sorted({padding, *(row for matrix in tree.values() for row, _ in matrix)})
    native_bound = F(0)
    observed = 0.
    visited = set()
    contributions = {'first_abort': F(0), 'final_abort': F(0)}
    for record, abort in menu.items():
        c1, raw1, bits, c2, raw2 = record
        key = (c1, raw1, bits[2], c2, raw2) if bits else record
        if key not in native_matrices:
            ideal = tree.get(key, {}) if bits else {}
            target = _dense(ideal)
            actual = target.copy()
            if key == nonzero_key:
                actual[padding, 0] += 1/16
            if key == first_abort_key:
                actual[padding, 1] += 1j/32
            if key == final_abort_key:
                actual[padding, 0] += 1/64
            exact = exact_native.compare_map(actual, ideal)['distance_upper']
            native_matrices[key] = actual, target, exact
        actual_native, target_native, native_error = native_matrices[key]
        if bits:
            current = {0: 0, 1: 1, 2: -1}[raw1 % 4]
            actual_work, target_work = work[current, bits]
            norm_key = current, bits
        else:
            actual_work, target_work = prepared, np.zeros_like(prepared)
            norm_key = None
        if norm_key not in norm_cache:
            norm_cache[norm_key] = certificate.nuclear_upper(actual_work)
        contribution = native_error*norm_cache[norm_key]
        native_bound += contribution
        if abort:
            contributions['final_abort' if bits else 'first_abort'] += contribution
        # Removing rows which are identically zero on both sides is an
        # exact common isometry; the padding perturbation stays included.
        a = np.kron(actual_native[reduced_rows], np.eye(2))
        b = np.kron(target_native[reduced_rows], np.eye(2))
        difference = np.kron(a@tau@a.conj().T, actual_work)-np.kron(b@tau@b.conj().T, target_work)
        observed += _nuclear(difference)/2
        visited.add(record)
    bound = certificate.compose(native_bound, work_errors)
    assert visited == set(menu)
    assert contributions['first_abort'] > 0 and contributions['final_abort'] > 0
    assert all(not np.array_equal(actual, actual.conj().T)
               for (current, bits), (actual, _) in work.items() if bits == (0, 0, 0, 0))
    assert observed > 0
    assert bound < F(1, 4)
    assert observed <= float(bound)+2e-12


def test_work_weighted_lift_has_no_native_record_multiplier():
    """Exact scalar errors saturate the weighted bound, not a branch count."""
    errors = (F(1, 16), F(1, 8), F(3, 16))
    assert certificate.compose(F(0), errors) == F(1, 8)
    # The four admitted first records carry 1/4 each; the two current-zero
    # records must both be present. This checks the sum independently.
    for law in ('C', 'E'):
        by_current = dict(zip((-1, 0, 1), errors, strict=True))
        exact_sum = sum((by_current[{0: 0, 1: 1, 2: -1}[raw % 4]]/4
                         for _, raw in FIRST[law]), F(0))
        assert exact_sum == certificate.compose(F(0), errors)
