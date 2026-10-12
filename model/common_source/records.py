"""Exact checkpoint record algebra and controlled-instrument adapter.

This executes the postcommit construction in model/observer/DERIVATION.md.
Channel completeness and Kraus intertwiners imply the dual fixed-observable
criterion of Lean/Dynamics/ProtectedCharge.lean. Writing a new record is a
separate event: raw record projectors are not fixed on the blank sector.

Dense matrices use immutable tuples of Gaussian rationals (real, imaginary),
compatible with the common-source gauge/matter matrix interface. Sparse native
operators retain the actual 3354-dimensional source/tape basis. Arbitrary
operators, including off-diagonal matrix units, are accepted; these routines
never normalize a branch or infer positivity from a numerical tolerance.
"""
from __future__ import annotations

from fractions import Fraction as Q
from itertools import product

from model.observer import native

ZERO = (Q(0), Q(0))
ONE = (Q(1), Q(0))


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _multiply(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])


def _conjugate(a):
    return (a[0], -a[1])


def _scalar(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(part) is not Q for part in value)):
        raise ValueError("Expected a pair of exact Fraction components")
    return value


def _matrix(value):
    if (type(value) is not tuple or not value
            or any(type(row) is not tuple or len(row) != len(value)
                   for row in value)):
        raise ValueError("Expected a nonempty immutable square matrix")
    for row in value:
        for entry in row:
            _scalar(entry)
    return value


def _identity(dimension):
    return tuple(tuple(ONE if row == column else ZERO
                       for column in range(dimension))
                 for row in range(dimension))


def _dagger(matrix):
    return tuple(tuple(_conjugate(matrix[column][row])
                       for column in range(len(matrix)))
                 for row in range(len(matrix)))


def _matmul(left, right):
    size = len(left)
    return tuple(tuple(_sum(_multiply(left[row][k], right[k][column])
                            for k in range(size))
                       for column in range(size))
                 for row in range(size))


def _sum(values):
    result = ZERO
    for value in values:
        result = _add(result, value)
    return result


def _shape(operator, record_count):
    operator = _matrix(operator)
    if (type(record_count) is not int or record_count < 1
            or len(operator) % record_count):
        raise ValueError("Positive record count must divide the matrix size")
    return operator, len(operator) // record_count


def block_readback(operator, record_count):
    """Pinch old record labels, retaining EVERY operator within each block.

    This conditional expectation fixes the generally noncommutative algebra
    direct_sum_r B(work). It is distinct from a classical public readback.
    """
    operator, dimension = _shape(operator, record_count)
    return tuple(tuple(value if row//dimension == column//dimension else ZERO
                       for column, value in enumerate(entries))
                 for row, entries in enumerate(operator))


def public_record_readback(operator, record_count):
    """Conditional expectation onto sum_r C Q_r using normalized block trace."""
    operator, dimension = _shape(operator, record_count)
    means = [_sum(operator[r*dimension+k][r*dimension+k]
                  for k in range(dimension)) for r in range(record_count)]
    means = [(a/dimension, b/dimension) for a, b in means]
    return tuple(tuple(means[row//dimension] if row == column else ZERO
                       for column in range(len(operator)))
                 for row in range(len(operator)))


def _controlled_kraus(record_kraus):
    """Build block diagonal Kraus maps; every record retains every outcome."""
    if (type(record_kraus) is not tuple or not record_kraus
            or any(type(branches) is not tuple or not branches
                   for branches in record_kraus)):
        raise ValueError("Expected record-indexed tuples of Kraus matrices")
    dimension = len(_matrix(record_kraus[0][0]))
    outcomes = len(record_kraus[0])
    for branches in record_kraus:
        if len(branches) != outcomes:
            raise ValueError("Every record must retain the complete outcome menu")
        for matrix in branches:
            if len(_matrix(matrix)) != dimension:
                raise ValueError("Every branch must act on the same work space")
        products = [_matmul(_dagger(matrix), matrix) for matrix in branches]
        gram = tuple(tuple(_sum(matrix[a][b] for matrix in products)
                           for b in range(dimension))
                     for a in range(dimension))
        if gram != _identity(dimension):
            raise ValueError("Each record's Kraus family must be exactly complete")
    size = len(record_kraus)*dimension
    return tuple(tuple(tuple(
        record_kraus[a//dimension][outcome][a % dimension][b % dimension]
        if a//dimension == b//dimension else ZERO
        for b in range(size)) for a in range(size)) for outcome in range(outcomes))


def controlled_instrument(operator, record_kraus):
    """Apply the complete old-record-controlled instrument without readback.

    record_kraus[r][y] is K_(r,y). Each returned unnormalized operator is the
    retained outcome y, including zero-probability outcomes. Old-record
    off-diagonal inputs are acted on as K_(r,y) X_rs K_(s,y)^dagger.
    """
    operator = _matrix(operator)
    kraus = _controlled_kraus(record_kraus)
    if len(operator) != len(kraus[0]):
        raise ValueError("Operator and controlled instrument dimensions differ")
    return tuple(_matmul(_matmul(matrix, operator), _dagger(matrix))
                 for matrix in kraus)


def controlled_certificate(record_kraus):
    """Exact Gram and record-projector intertwiners on the WHOLE carrier."""
    kraus = _controlled_kraus(record_kraus)
    count = len(record_kraus)
    size = len(kraus[0])
    dimension = size//count
    products = [_matmul(_dagger(matrix), matrix) for matrix in kraus]
    gram = tuple(tuple(_sum(matrix[a][b] for matrix in products)
                       for b in range(size)) for a in range(size))
    if gram != _identity(size):
        raise AssertionError("Controlled Gram identity failed")
    for matrix in kraus:
        if any(matrix[a][b] != ZERO for a in range(size) for b in range(size)
               if a//dimension != b//dimension):
            raise AssertionError("An old record is changed")
    return {"record_count": count, "work_dimension": dimension,
            "retained_outcomes": len(kraus), "gram_identity_exact": True,
            "old_record_kraus_intertwiners_exact": True,
            "dual_fixes_all_old_record_projectors": True,
            "arbitrary_operator_and_reference_scope": True,
            "whole_work_block_algebra_preservation_claimed": False}


def _state(state):
    if (type(state) is not native.State or type(state.loads) is not tuple
            or len(state.loads) != 12 or any(type(v) is not int or v < 0
                                           for v in state.loads)
            or sum(state.loads) != 2 or type(state.prefix) is not tuple
            or len(state.prefix) > 2
            or any(type(r) is not tuple or len(r) != 2
                   or any(type(v) is not int for v in r) or r not in native.LETTERS
                   for r in state.prefix)):
        raise ValueError("Expected a state of the actual native source/tape")
    return state


def _sparse(operator):
    if type(operator) is not dict:
        raise ValueError("Expected a sparse native operator dictionary")
    for pair, coefficient in operator.items():
        if type(pair) is not tuple or len(pair) != 2:
            raise ValueError("Expected two source/tape basis states")
        for state in pair:
            _state(state)
        _scalar(coefficient)
    return operator


def _accumulate(output, pair, coefficient):
    updated = _add(output.get(pair, ZERO), coefficient)
    if updated == ZERO:
        output.pop(pair, None)
    else:
        output[pair] = updated


def native_tick_operator(operator):
    """The existing coherent T=A0(.)A0*+A1(.)A1*+P2(.)P2, exactly.

    Equal-amplitude native write columns make every matrix-unit coefficient
    rational: nonfull columns have amplitude 1/sqrt(2), giving product 1/2.
    This is NOT an incoherent coin or current measurement. Different phase
    sectors are separate Kraus branches, as in the original construction.
    """
    output = {}
    for (left, right), coefficient in _sparse(operator).items():
        if left.phase != right.phase:
            continue
        if left.phase == 2:
            _accumulate(output, (left, right), coefficient)
            continue
        for weight_left, after_left in native.tick(left):
            for weight_right, after_right in native.tick(right):
                if weight_left != Q(1, 2) or weight_right != Q(1, 2):
                    raise AssertionError("Native write amplitudes changed")
                _accumulate(output, (after_left, after_right),
                            (coefficient[0]/2, coefficient[1]/2))
    return output


def native_block_readback(operator):
    """Pinch the first committed letter, with blank as a separate block."""
    return {pair: coefficient for pair, coefficient in _sparse(operator).items()
            if pair[0].prefix[:1] == pair[1].prefix[:1]}


def native_first_record_matrix(operator):
    """Trace work/phase/second-slot data on the already committed carrier."""
    result = [[ZERO for _ in native.LETTERS] for _ in native.LETTERS]
    index = {record: i for i, record in enumerate(native.LETTERS)}
    for (left, right), coefficient in _sparse(operator).items():
        if not left.prefix or not right.prefix:
            raise ValueError("First record has not been committed")
        if (left.phase == right.phase and left.loads == right.loads
                and left.prefix[1:] == right.prefix[1:]):
            a, b = index[left.prefix[0]], index[right.prefix[0]]
            result[a][b] = _add(result[a][b], coefficient)
    return tuple(tuple(row) for row in result)


def native_record_weight(operator, record):
    """Linear expectation of a first-letter projector, zero on blank tape."""
    if record not in native.LETTERS:
        raise ValueError("Unknown native record")
    return _sum(coefficient for (left, right), coefficient in _sparse(operator).items()
                if left == right and left.prefix[:1] == (record,))


def native_checkpoint_certificate():
    """Check all write columns and their checkpoint intertwiners exactly.

    Orthogonal column supports and norms prove the full Gram identity, not
    just classical conservation. The phase Kraus split and retained first
    letter then imply the dual projector identity on every coherent input.
    """
    counts = {}
    for phase in (0, 1, 2):
        seen = set()
        columns = 0
        for prefix in product(native.LETTERS, repeat=phase):
            for loads in native.states():
                before = native.State(loads, prefix)
                images = native.tick(before)
                if sum(weight for weight, _ in images) != 1:
                    raise AssertionError("Native column norm differs from one")
                support = {after for _, after in images}
                if len(support) != len(images) or seen.intersection(support):
                    raise AssertionError("Native columns are not orthogonal")
                seen.update(support)
                if phase and any(after.prefix[:1] != prefix[:1]
                                 for _, after in images):
                    raise AssertionError("A committed first record changed")
                columns += 1
        counts[str(phase)] = columns
    blank = native.State((1, 1) + (0,)*10)
    operator = {(blank, blank): ONE}
    record = (0, 0)
    before = native_record_weight(operator, record)
    after = native_record_weight(native_tick_operator(operator), record)
    if before != ZERO or after != (Q(1, 2), Q(0)):
        raise AssertionError("Record-creation boundary changed")
    return {"source_tape_dimension": sum(counts.values()),
            "all_phase_column_counts": counts,
            "postcommit_dimension": counts["1"] + counts["2"],
            "all_phase_gram_identities_exact": True,
            "postcommit_first_record_intertwiners_exact": True,
            "postcommit_dual_projectors_fixed": True,
            "arbitrary_operator_and_reference_scope": True,
            "global_blank_projector_counterexample": {
                "first_record": list(record), "before": "0", "after": "1/2"},
            "full_operational_observer_inhabitant_claimed": False}
