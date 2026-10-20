"""Source words for the actual finite, successfully prepared edge reference.

The declared sampler and herald select M=sum_n w_n (I-L/120)^(2n).
Its positive real logarithm fixes H=-log(M/m(0)), without spectral fitting or
integer phase branches.  The source clock is exp(-i H), so its Heisenberg
action obeys sigma_s=alpha_{-s} for rho=M/Tr(M).  The separate sqrt(L) clock
is a comparison, not the modular generator of this finite prepared state.

All routines act on the same sixteen-address protected register.  Four
padding addresses are retained.  Matrices returned as ``raw`` are executed
source words or products of those executed routines, never inserted targets
or polar-corrected matrices.  Fixed numerical guards are explicitly not
directed-rounding certificates or physical apparatus error measurements.
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from . import inverse_clock
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*inverse_clock.SOURCE_PATHS,
    'model/predictions/prepared_clock_source.py',
    'model/observer/source_scalar/compiler.py'))))
DIMENSION = 16
QUBITS = 4
SYNTHESIS_BITS = 44
BASIS_GUARD = F(1, 10**8)
PHASE_GUARD = F(1, 10**10)
PAIR_GUARD = F(1, 10**9)
BAND_ADDRESSES = ((0,), (1, 3, 5), (2, 6, 7, 8, 9), (4, 10, 11))
PAIR_ADDRESSES = (1, 2, 4)


def _dyadic(expression):
    """Evaluate symbolic inputs before one reproducible source-grid rounding."""
    return F(int(sp.floor(sp.N(expression, 90)*(1 << SYNTHESIS_BITS)
                          + sp.Rational(1, 2))), 1 << SYNTHESIS_BITS)


def _quantize(word):
    return [compiler.pulse(event['kind'],
            round(event['angle']*(1 << SYNTHESIS_BITS))/(1 << SYNTHESIS_BITS),
            event['targets']) for event in word]


def _float_matrix(matrix):
    return np.array([[float(sp.N(x, 80)) for x in row]
                     for row in matrix.tolist()])


def _basis_word(exact_basis):
    """Gray elimination at fixed precision before rounding emitted angles.

    The public source compiler's primitive Gray rotations and Walsh phases
    are reused.  Deriving the elimination angles here from the exact basis
    avoids making source-word identities depend on platform libm residuals
    at the many zero pivots of a degenerate eigenspace.
    """
    gray = [j ^ (j >> 1) for j in range(16)]
    rotations = []
    with mp.workdps(90):
        residual = [[mp.mpf(str(sp.N(exact_basis[a, b], 100)))
                     for b in gray] for a in gray]
        for column in range(15):
            for row in range(15, column, -1):
                x, y = residual[row-1][column], residual[row][column]
                theta = mp.mpf(0) if x == 0 and y == 0 else mp.atan2(y, x)
                c, s = mp.cos(theta), mp.sin(theta)
                first, second = residual[row-1][:], residual[row][:]
                residual[row-1] = [c*a+s*b for a, b in zip(first, second, strict=True)]
                residual[row] = [-s*a+c*b for a, b in zip(first, second, strict=True)]
                angle = F(int(mp.floor(theta*(1 << SYNTHESIS_BITS)+mp.mpf('0.5'))),
                          1 << SYNTHESIS_BITS)
                rotations.append((gray[row-1], gray[row], angle))
        signs = [1 if residual[j][j] >= 0 else -1 for j in range(16)]
        if max(abs(residual[a][b]-(signs[a] if a == b else 0))
               for a in range(16) for b in range(16)) > mp.mpf('1e-75'):
            raise ArithmeticError('fixed-precision Gray elimination lost its complete basis')
    phases = [0.]*16
    for j, sign in enumerate(signs):
        phases[gray[j]] = 0. if sign == 1 else float(_dyadic(sp.pi))
    word = compiler.compile_diagonal(phases)
    for a, b, angle in reversed(rotations):
        word.extend(compiler._gray_rotation(a, b, -float(angle), 4))
    return _quantize(word)


@lru_cache(maxsize=1)
def _spectral_data():
    """Exact, deterministic Gram-Schmidt inside each native eigenspace."""
    laplacian = sp.zeros(12)
    seams = tuple(inverse_clock.continuation.quantum_mixing.seams())
    if len(seams) != 30 or len(set(seams)) != 30:
        raise ArithmeticError('the declared native seam census changed')
    for a, b in seams:
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    bands = (sp.Integer(0), 5-sp.sqrt(5), sp.Integer(6), 5+sp.sqrt(5))
    basis = sp.zeros(16)
    values = [sp.Integer(0)]*16
    for eigenvalue, addresses in zip(bands, BAND_ADDRESSES, strict=True):
        vectors = (laplacian-eigenvalue*sp.eye(12)).nullspace()
        if len(vectors) != len(addresses):
            raise ArithmeticError('the native eigenspace multiplicity changed')
        orthogonal = []
        for address, vector in zip(addresses, vectors, strict=True):
            for previous in orthogonal:
                vector = (vector-previous.dot(vector)/previous.dot(previous)*previous
                          ).applyfunc(sp.simplify)
            first = next(x for x in vector if x != 0)
            if float(sp.N(first, 80)) < 0:
                vector = -vector
            if any(sp.simplify(x) != 0 for x in laplacian*vector-eigenvalue*vector):
                raise ArithmeticError('an exact basis vector left its native eigenspace')
            squared = sp.simplify(vector.dot(vector))
            if squared.is_positive is not True:
                raise ArithmeticError('an exact native mode has no positive norm')
            for previous in orthogonal:
                if sp.simplify(previous.dot(vector)) != 0:
                    raise ArithmeticError('exact native mode orthogonality failed')
            orthogonal.append(vector)
            basis[:12, address] = vector/sp.sqrt(squared)
            values[address] = eigenvalue
    for address in range(12, 16):
        basis[address, address] = 1
    array = _float_matrix(basis)
    if (np.linalg.norm(array.T@array-np.eye(16), 'fro') > 2e-13
            or any(sp.simplify(x-1/sp.sqrt(12)) != 0 for x in basis[:12, 0])):
        raise ArithmeticError('the normalized native basis or constant mode changed')
    return {'laplacian': sp.ImmutableMatrix(laplacian),
        'basis_exact': sp.ImmutableMatrix(basis), 'basis': array,
        'eigenvalues': tuple(values), 'band_eigenvalues': bands,
        'band_addresses': BAND_ADDRESSES, 'pair_addresses': PAIR_ADDRESSES,
        'seams': seams,
        'exact_basis_certificate': {
            'symmetric_source_laplacian': laplacian == laplacian.T,
            'all_twelve_eigenvector_equations_checked_exactly': True,
            'within_each_band_orthogonality_and_positive_norms_checked_exactly': True,
            'cross_band_orthogonality': 'distinct eigenvalues of the exact real symmetric L',
            'full_space_completeness': '1+3+5+3=12 mutually orthogonal nonzero vectors',
            'padding_extension': 'four orthogonal standard coordinate vectors'}}


def native_spectral_data():
    """Detached exact source L and basis; no caller supplies an eigenvalue."""
    return deepcopy(_spectral_data())


@lru_cache(maxsize=1)
def _reference_data():
    spectrum = _spectral_data()
    weights = inverse_clock.subordinated_filter_weights()
    if (any(type(w) is not F or w < 0 for w in weights)
            or not 0 < weights[0] <= sum(weights) <= 1
            or not any(weights[1:])):
        raise ArithmeticError('the declared finite sampler is not positive and nontrivial')
    lazy = sp.eye(12)-spectrum['laplacian']/120
    polynomial = sp.zeros(12)
    for n, weight in enumerate(weights):
        polynomial += sp.Rational(weight)*lazy**(2*n)
    mass0 = sum(weights)
    masses = tuple(sp.simplify(sum((sp.Rational(w)*(1-lam/120)**(2*n)
                                  for n, w in enumerate(weights)), sp.Integer(0)))
                   for lam in spectrum['eigenvalues'][:12])
    gaps = tuple(sp.log(sp.Rational(mass0)/mass) for mass in masses)
    normalization = sp.trace(polynomial)
    density = polynomial/normalization
    if (normalization <= 0 or any(not 0 < float(sp.N(mass, 80)) <= float(mass0)
                                  for mass in masses)
            or sp.simplify(sum(masses)-normalization) != 0):
        raise ArithmeticError('finite prepared reference lost its positive spectral law')
    extended_density = np.zeros((16, 16), complex)
    extended_density[:12, :12] = _float_matrix(density)
    return {'weights': weights, 'lazy_filter': sp.ImmutableMatrix(lazy),
        'polynomial': sp.ImmutableMatrix(polynomial),
        'normalization': normalization, 'm0': mass0,
        'density_exact': sp.ImmutableMatrix(density), 'density': extended_density,
        'masses': masses, 'modular_gaps_exact': gaps,
        'modular_gaps': tuple(float(sp.N(gap, 80)) for gap in gaps),
        'reference_algebra': 'B(C^12); the four retained padding addresses are outside its support',
        'modular_convention': 'U(t)=exp(-i*t*H); sigma_s=alpha_(-s), H=-log(M/m0)',
        'logarithm_choice': 'positive real logarithm of the positive prepared polynomial',
        'temperature_rescaled_to_match_sqrt_L': False,
        'padding_clock_extension': 'identity',
        'selection': 'fixed supplied subordination sampler and all-success herald'}


def prepared_reference():
    """Actual finite sampler reference, not its asymptotic sqrt(L) target."""
    return deepcopy(_reference_data())


def _band(band):
    if type(band) is not int or band not in (1, 2, 3):
        raise ValueError('an original nonconstant native band 1, 2 or 3 is required')
    return PAIR_ADDRESSES[band-1]


def _pair_spectral_word(band, axis, angle):
    """Rank-two X/Y pulse; the other fourteen addresses are left untouched."""
    address = _band(band)
    if type(axis) is not str or axis not in ('X', 'Y'):
        raise ValueError('a declared X or Y pair axis is required')
    target = address.bit_length()-1
    result = []
    for mask in range(16):
        if mask & address:
            continue
        label = [axis if j == target else 'Z' if mask >> j & 1 else 'I'
                 for j in range(4)]
        result.extend(compiler.pauli_rotation(label, float(angle)/8))
    return _quantize(result)


def _pair_target(band, axis, angle):
    address = _band(band)
    result = np.eye(16, dtype=complex)
    c, s = math.cos(float(angle)), math.sin(float(angle))
    result[0, 0] = result[address, address] = c
    result[0, address] = 1j*s if axis == 'X' else s
    result[address, 0] = 1j*s if axis == 'X' else -s
    return result


def _check_matrix(raw, target, guard, name):
    for label, matrix in (('raw source', raw), ('target', target)):
        if (type(matrix) is not np.ndarray or matrix.shape != (16, 16)
                or matrix.dtype.kind not in 'iufc' or not np.all(np.isfinite(matrix))):
            raise ArithmeticError(f'{label} {name} is not a complete finite register matrix')
    if (type(guard) is not F or type(guard.numerator) is not int
            or type(guard.denominator) is not int or guard.numerator < 0
            or guard.denominator <= 0 or math.gcd(guard.numerator, guard.denominator) != 1):
        raise ArithmeticError(f'source {name} requires a finite nonnegative rational guard')
    try:
        error_limit = float(guard)
        defect_limit = float(2*guard+guard*guard)
    except OverflowError as exc:
        raise ArithmeticError(
            f'source {name} requires a finite nonnegative rational guard') from exc
    if not math.isfinite(error_limit) or not math.isfinite(defect_limit):
        raise ArithmeticError(f'source {name} requires a finite nonnegative rational guard')
    with np.errstate(over='ignore', invalid='ignore'):
        # Integer matrix products can wrap to an apparently unitary matrix;
        # lower-precision floating norms can also underflow before comparison.
        numeric_raw = raw.astype(np.result_type(raw.dtype, np.complex128), copy=False)
        numeric_target = target.astype(np.result_type(target.dtype, np.complex128), copy=False)
        difference = float(np.linalg.norm(numeric_raw-numeric_target, 'fro'))
        defect = float(np.linalg.norm(numeric_raw.conj().T@numeric_raw-np.eye(16), 'fro'))
    if (not math.isfinite(difference) or not math.isfinite(defect)
            or difference > error_limit or defect > defect_limit):
        raise ArithmeticError(f'raw source {name} exceeds its complete-register numerical guard')
    return {'complete_target_Frobenius_error': difference,
        'raw_unitarity_Frobenius_defect': defect,
        'complete_target_error_guard': guard,
        'interval_or_physical_error_certificate': False}


def _execute(word, target, guard, name):
    raw = compiler.execute(word, 4)
    return raw, _check_matrix(raw, target, guard, name)


def _combined_guard(*guards):
    result = F(1)
    for guard in guards:
        result *= 1+guard
    return result-1


@lru_cache(maxsize=1)
def _compiled_data():
    spectrum, reference = _spectral_data(), _reference_data()
    basis = spectrum['basis']
    words, raw, targets, guards, diagnostics = {}, {}, {}, {}, {}
    words['basis'] = _basis_word(spectrum['basis_exact'])
    words['unbasis'] = compiler.inverse(words['basis'])
    targets.update(basis=basis, unbasis=basis.T)
    for name in ('basis', 'unbasis'):
        guards[name] = BASIS_GUARD
        raw[name], diagnostics[name] = _execute(words[name], targets[name], guards[name], name)
    for name, values in (
            ('modular_phase', reference['modular_gaps_exact']),
            ('sqrt_phase', tuple(sp.sqrt(x) for x in spectrum['eigenvalues'][:12]))):
        # Padding has exactly zero requested phase.  Targets independently use
        # high-precision unrounded symbolic phases and all sixteen columns.
        phases = tuple(-_dyadic(value) for value in values)+(F(0),)*4
        words[name] = _quantize(compiler.compile_diagonal([float(x) for x in phases]))
        wanted = np.array([-float(sp.N(value, 80)) for value in values]+[0.]*4)
        targets[name] = np.diag(np.exp(1j*wanted))
        guards[name] = PHASE_GUARD
        raw[name], diagnostics[name] = _execute(words[name], targets[name], guards[name], name)
    for name, phase in (('clock1', 'modular_phase'), ('sqrt_clock1', 'sqrt_phase')):
        words[name] = words['unbasis']+words[phase]+words['basis']
        raw[name] = raw['basis']@raw[phase]@raw['unbasis']
        targets[name] = basis@targets[phase]@basis.T
        guards[name] = _combined_guard(BASIS_GUARD, PHASE_GUARD, BASIS_GUARD)
        diagnostics[name] = _check_matrix(raw[name], targets[name], guards[name], name)
    for band in (1, 2, 3):
        settings = (
            (f'kick_plus_{band}', 'X', sp.pi/4, True),
            (f'kick_minus_{band}', 'X', -sp.pi/4, True),
            (f'detector_x_{band}', 'Y', sp.pi/4, False),
            (f'detector_y_{band}', 'X', -sp.pi/4, False))
        for name, axis, angle, return_to_site in settings:
            pair_word = _pair_spectral_word(band, axis, _dyadic(angle))
            target = _pair_target(band, axis, float(sp.N(angle, 80)))
            pair_raw, pair_diagnostic = _execute(pair_word, target, PAIR_GUARD, name+'_pair')
            words[name] = words['unbasis']+pair_word
            raw[name] = pair_raw@raw['unbasis']
            targets[name] = target@basis.T
            guards[name] = _combined_guard(BASIS_GUARD, PAIR_GUARD)
            if return_to_site:
                words[name] += words['basis']
                raw[name] = raw['basis']@raw[name]
                targets[name] = basis@targets[name]
                guards[name] = _combined_guard(BASIS_GUARD, PAIR_GUARD, BASIS_GUARD)
            diagnostics[name] = _check_matrix(raw[name], targets[name], guards[name], name)
            diagnostics[name]['pair_routine'] = pair_diagnostic
    digests = {name: hashlib.sha256(json.dumps(word, sort_keys=True,
                separators=(',', ':'), allow_nan=False).encode()).hexdigest()
               for name, word in words.items()}
    # Exact-algebra identity, numerically checked on the actual finite density;
    # the source implementation allowance is retained separately above.
    state = reference['density']
    stationarity = float(np.linalg.norm(targets['clock1']@state
                                       @targets['clock1'].conj().T-state, 'fro'))
    if stationarity > 2e-13:
        raise ArithmeticError('finite reference and its positive-log clock disagree')
    return {'words': words, 'raw': raw, 'targets': targets, 'guards': guards,
        'diagnostics': diagnostics, 'word_sha256': digests,
        'pulse_counts': {name: len(word) for name, word in words.items()},
        'pair_addresses': {band: _band(band) for band in (1, 2, 3)},
        'reference': reference, 'spectrum': spectrum,
        'ideal_reference_stationarity_error': stationarity,
        'source_synthesis_grid': F(1, 1 << SYNTHESIS_BITS),
        'raw_source_routines_used_without_polar_correction': True,
        'complete_sixteen_address_maps_checked_with_numerical_guards': True,
        'modular_clock_is_sqrt_L_clock': False,
        'physical_seconds_or_autonomous_clock_claimed': False}


def compiled_source():
    """Fresh words and raw source maps for the fixed declared experiment.

    ``clock1`` and ``sqrt_clock1`` act in the site basis.  ``kick_plus_b``
    and ``kick_minus_b`` are exp(+/-i*pi*X_0j/4) in the same register.
    ``detector_x_b`` and ``detector_y_b`` end in the spectral address basis;
    address-zero minus address-j has expectation X_0j or Y_0j.  All other
    addresses, including padding, remain outcomes.  Band 2 is lambda=6.

    Reusing the executed raw routines is linear source execution, not a
    replacement of the prepared density.  Untouched reference columns can
    be carried through each returned map without modification.
    """
    return deepcopy(_compiled_data())


def clock_word(kind='modular'):
    if type(kind) is not str or kind not in ('modular', 'sqrt'):
        raise ValueError('the declared modular or sqrt comparison clock is required')
    return deepcopy(_compiled_data()['words']['clock1' if kind == 'modular' else 'sqrt_clock1'])


def pair_word(band, sign=1):
    _band(band)
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('an original signed pair kick is required')
    return deepcopy(_compiled_data()['words'][f'kick_{"plus" if sign == 1 else "minus"}_{band}'])


def detector_word(band, axis='X'):
    _band(band)
    if type(axis) is not str or axis not in ('X', 'Y'):
        raise ValueError('a declared X or Y detector is required')
    return deepcopy(_compiled_data()['words'][f'detector_{axis.lower()}_{band}'])
