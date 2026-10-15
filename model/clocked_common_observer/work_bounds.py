"""Finite all-input source-work bounds, with exact binary-array arithmetic.

The supplied work state is arbitrary: the sixteen normalized basis columns
replace the vacuum premise of smooth_history's prepared-state certificate.
The certified maps are the stored complete branch operators interpreted
exactly. Later floating applications are charged by exact_density_residual.
Only dimension two, the declared coupling and an empty environment are covered.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import hashlib
import math

import numpy as np

from ..common_source import exact_output as exact
from ..common_source import infinite_history as ih, protocol, smooth_history as sh
from ..primitive_observer import common, runtime


DIMENSION = 16
SOURCE_PATHS = (
    'model/clocked_common_observer/work_bounds.py',
    'model/common_source/exact_output.py', 'model/common_source/smooth_history.py',
    'model/common_source/infinite_history.py', 'model/common_source/protocol.py',
    'model/primitive_observer/common.py', 'model/primitive_observer/runtime.py',
)


def _matrix(array):
    return exact._array(array, (DIMENSION, DIMENSION))


@lru_cache(maxsize=128)
def _integer_snapshot(dtype, payload):
    array = np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(DIMENSION, DIMENSION)
    values = tuple(exact._component(x) for x in array.flat)
    denominator = max((x.denominator for pair in values for x in pair), default=1)
    if denominator & (denominator-1):
        raise ArithmeticError('binary arrays require power-of-two denominators')
    real = np.array([x.numerator*(denominator//x.denominator) for x, _ in values], dtype=object)
    imag = np.array([x.numerator*(denominator//x.denominator) for _, x in values], dtype=object)
    real, imag = real.reshape(array.shape), imag.reshape(array.shape)
    real.flags.writeable = imag.flags.writeable = False
    return real, imag, denominator


def exact_matrix(array):
    """Detach original 16x16 binary entries into exact integer lanes.

    Returned lanes are copies, so caller mutation cannot poison the private
    cache. No float matrix multiplication, norm or normalization occurs.
    """
    array = _matrix(array)
    real, imag, denominator = _integer_snapshot(array.dtype.str, array.tobytes())
    return real.copy(), imag.copy(), denominator


def _multiply(left, right):
    lr, li, ld = left
    rr, ri, rd = right
    return lr@rr-li@ri, lr@ri+li@rr, ld*rd


def _adjoint(value):
    real, imag, denominator = value
    return real.T, -imag.T, denominator


def _difference_squared(left, right):
    lr, li, ld = left
    rr, ri, rd = right
    denominator = math.lcm(ld, rd)
    real = lr*(denominator//ld)-rr*(denominator//rd)
    imag = li*(denominator//ld)-ri*(denominator//rd)
    return F(sum(int(a)*int(a)+int(b)*int(b) for a, b in zip(real.flat, imag.flat, strict=True)),
             denominator*denominator)


def exact_product_residual(computed, left, right):
    """Bound ||computed-left@right||_F using the original arrays exactly."""
    squared = _difference_squared(exact_matrix(computed),
                                  _multiply(exact_matrix(left), exact_matrix(right)))
    return exact._sqrt_interval(squared)[1]


def exact_density_residual(computed, factors, initial):
    """Half-nuclear error of a returned density versus exact binary factors.

    Factors are in execution order: (P,K) means K P initial P* K*. The
    bound applies to general matrices, without positivity or Hermiticity.
    It charges floating execution, not a change to any ideal source action.
    """
    if type(factors) is not tuple or len(factors) > 64:
        raise ValueError('at most 64 immutable execution-ordered factors required')
    target = exact_matrix(initial)
    for factor in factors:
        matrix = exact_matrix(factor)
        target = _multiply(_multiply(matrix, target), _adjoint(matrix))
    squared = _difference_squared(exact_matrix(computed), target)
    return exact._sqrt_interval(DIMENSION*squared)[1]/2


def nuclear_upper(array):
    real, imag, denominator = exact_matrix(array)
    squared = F(sum(int(a)*int(a)+int(b)*int(b) for a, b in zip(real.flat, imag.flat, strict=True)),
                denominator*denominator)
    return exact._sqrt_interval(DIMENSION*squared)[1]


def certify_preparation(initial):
    """Actual d2 preparation against its exact normalized finite target.

    The existing finite preparation is ih._initial(d-2), hence |0,0,0>
    at d=2. This is not a certificate against the infinite free vacuum.
    """
    declared = ih._initial(0)
    if declared.coefficients != {(0, 0, 0): (1, 0)} or declared.denominator != 1:
        raise ArithmeticError('the declared finite preparation changed')
    target = np.zeros((DIMENSION, DIMENSION), dtype=int)
    target[0, 0] = 1
    return exact_density_residual(initial, (), target)


def _validate_program(program):
    program = protocol.validate_program(program)
    if (program.dimension != 2 or program.coupling != F(1, 4)
            or not 1 <= program.steps <= 4096
            or any(op.kind == 'sharp_read' for op in program.operations)
            or sum(op.kind in protocol.READS for op in program.operations) > 4
            or sum(op.kind == 'smooth_read' for op in program.operations) > 2
            or sum((op.time for op in program.operations), F(0)) > 4):
        raise ValueError('bounded d2 finite smooth work program at coupling 1/4 required')
    return program


def _basis(index):
    gauge, occupation = index % 4, index//4
    return ih._Polynomial({(occupation % 2, occupation//2, gauge): (1, 0)}, 1)


@lru_cache(maxsize=16)
def _column_targets(operations, current):
    """Exact finite actions and carried error on every orthonormal column.

    At d2 all factorial Fock norms are one. The existing error induction
    uses complete isometric record dilations, so it applies to every basis
    vector without a vacuum or initial-moment premise.
    """
    resolution = sh.Resolution()
    columns, errors = [], []
    for index in range(DIMENSION):
        branches = {(): sh._kick(_basis(index), current)}
        total_error = F(0)
        for operation in operations:
            following, local_errors = {}, []
            for record, vector in branches.items():
                local = {}
                if operation.kind == 'evolve':
                    value, local, _ = sh._evolve(vector, operation.time, F(1), resolution, 2)
                    following[record] = value
                elif operation.kind == 'smooth_read':
                    values, local, _ = sh._pointer(vector, F(1), resolution, 2)
                    following.update({(*record, bit): value for bit, value in enumerate(values)})
                elif operation.kind == 'gauge_read':
                    following.update({(*record, bit): ih._first_record(vector, bit, False)
                                      for bit in (0, 1)})
                elif operation.kind == 'feedback':
                    following[record] = sh._kick(vector, int(operation.enabled and record[operation.record]))
                else:
                    raise ValueError('unsupported exact finite operation')
                local_errors.append(sum(local.values(), F(0)))
            total_error += exact._sqrt_interval(sum((x*x for x in local_errors), F(0)))[1]
            branches = following
        columns.append(branches)
        errors.append(total_error)
    return tuple(columns), tuple(errors)


def _source_maps(compiled, include_kick):
    events = compiled.events if include_kick else compiled.events[1:]
    maps = {(): np.eye(DIMENSION, dtype=complex)}
    for event in events:
        if event.sequence:
            raise ValueError('the finite pure source-map certificate has no hidden sequences')
        following = {}
        for record, previous in maps.items():
            instrument = event.instruments[0 if event.control is None else record[event.control]]
            if instrument.schur is not None or any(len(group) != 1 for group in instrument.kraus):
                raise ValueError('one source Kraus map per complete outcome required')
            for bit, group in enumerate(instrument.kraus):
                following[(*record, bit) if event.visible else record] = group[0]@previous
        maps = following
    return tuple((record, runtime._immutable(matrix)) for record, matrix in maps.items())


@dataclass(frozen=True)
class WorkCertificate:
    program: protocol.Program
    includes_current_kick: bool
    maps: tuple
    stinespring_error_upper: F
    half_diamond_upper: F
    channel_norm_upper: F
    polynomial_residual_upper: F
    polynomial_remainder_upper: F
    source_maps_sha256: str
    compiled_program: runtime.Program


def certify_program(program, *, include_kick=True):
    """Certify actual composed source maps against the exact finite process.

    The stored rounded matrix products are compared directly with the exact
    target, so their internal product rounding is included. Later uses of
    these matrices need exact_density_residual or another arithmetic bound.
    """
    program = _validate_program(program)
    if type(include_kick) is not bool or (not include_kick and program.native_current != 0):
        raise ValueError('an explicit kick policy consistent with the current required')
    compiled = common.compile_program(program)
    maps = _source_maps(compiled, include_kick)
    columns, errors = _column_targets(program.operations, program.native_current if include_kick else 0)
    expected_records = tuple(product((0, 1), repeat=sum(op.kind in protocol.READS for op in program.operations)))
    if tuple(record for record, _ in maps) != expected_records:
        raise ArithmeticError('the complete source record stack changed')
    squared = F(0)
    digest = hashlib.sha256()
    for record, matrix in maps:
        digest.update(repr(record).encode())
        digest.update(matrix.tobytes())
        for row, column in product(range(DIMENSION), repeat=2):
            value = columns[column][record]
            key = (row//4 % 2, row//8, row % 4)
            real, imag = value.coefficients.get(key, (0, 0))
            ar, ai = exact._component(matrix[row, column])
            squared += (ar-F(real, value.denominator))**2+(ai-F(imag, value.denominator))**2
    residual = exact._sqrt_interval(squared)[1]
    remainder = exact._sqrt_interval(sum((e*e for e in errors), F(0)))[1]
    delta = residual+remainder
    # ||V||=1 for the ideal complete record dilation, ||A||<=1+delta.
    # Tensor extension and final record dephasing do not change the estimate.
    return WorkCertificate(program, include_kick, maps, delta, delta*(2+delta)/2,
                           (1+delta)**2, residual, remainder, digest.hexdigest(), compiled)


def certify_step(time=F(1, 4), steps=128):
    if type(time) is not F or time not in (F(1, 4), F(1, 2)):
        raise ValueError('the declared exact quarter or half model-time step required')
    return certify_program(protocol.Program(2, steps, operations=(protocol.Operation('evolve', time),)),
                           include_kick=False)


def certify_suffix(current, steps=128):
    if type(current) is not int or current not in (-1, 1):
        raise ValueError('a signed first-current outcome required')
    return certify_program(protocol.continuation_program(2, steps, native_current=current))


def summary(certificate):
    if type(certificate) is not WorkCertificate:
        raise ValueError('a computed finite work certificate required')
    return {'dimension': 2, 'work_dimension': DIMENSION,
            'canonical_program': certificate.program.specification(),
            'vacuum_preparation_included': False,
            'includes_current_kick': certificate.includes_current_kick,
            'retained_records': len(certificate.maps),
            'stinespring_error_upper': str(certificate.stinespring_error_upper),
            'half_diamond_upper': str(certificate.half_diamond_upper),
            'channel_norm_upper': str(certificate.channel_norm_upper),
            'polynomial_residual_upper': str(certificate.polynomial_residual_upper),
            'polynomial_remainder_upper': str(certificate.polynomial_remainder_upper),
            'source_maps_sha256': certificate.source_maps_sha256,
            'arbitrary_work_reference_input': True,
            'subsequent_float_applications_included': False,
            'physical_service_errors_included': False,
            'infinite_occupation_comparison': False}
