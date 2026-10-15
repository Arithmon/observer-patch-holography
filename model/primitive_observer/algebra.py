"""Observer algebras of the actual retained classical/quantum history.

Tensor factors and access are supplied apparatus choices.  Restriction keeps
every public record and its unnormalized mass, while discarding the untouched
reference.  It does not provide autonomous dynamics of an interacting region.
GNS and modular computations reuse the existing finite matrix implementation;
their numerical rank tolerance is not an eigenvalue certificate or a clock.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from collections.abc import Mapping

import numpy as np

from ..common_source.operational import MatrixGNS, modular_flow as _modular_flow


TOLERANCE = 2e-9
MAX_MATRIX_ENTRIES = 2_000_000


def _matrix(value, size, name='operator'):
    if np.ma.isMaskedArray(value):
        raise ValueError(f'{name} cannot be masked')
    def original(item):
        if isinstance(item, (bool, np.bool_)):
            raise ValueError('Boolean quantum entries are inadmissible')
        if isinstance(item, (list, tuple)):
            for child in item:
                original(child)
    original(value)
    raw = np.asarray(value)
    if (raw.shape != (size, size) or raw.dtype.kind not in 'iufc'
            or raw.size > MAX_MATRIX_ENTRIES or not np.isfinite(raw).all()):
        raise ValueError(f'{name} must be a finite numeric {size} by {size} matrix')
    result = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(result).all():
        raise ValueError(f'{name} must remain finite in binary64')
    return result


def _immutable(value):
    """Detached bytes prevent callers from re-enabling writes on retained states."""
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


def _layout(factors, region):
    if (not isinstance(factors, tuple) or not factors
            or any(type(n) is not int or n < 1 for n in factors)):
        raise ValueError('supply a nonempty tuple of positive tensor dimensions')
    if (not isinstance(region, tuple)
            or any(type(i) is not int or i not in range(len(factors)) for i in region)
            or len(set(region)) != len(region)):
        raise ValueError('region must contain distinct supplied tensor factor indices')
    kept = tuple(i for i in range(len(factors)) if i in region)
    rest = tuple(i for i in range(len(factors)) if i not in kept)
    size = math.prod(factors)
    if size*size > MAX_MATRIX_ENTRIES:
        raise ValueError('declared tensor algebra exceeds the dense access budget')
    return kept, rest, math.prod(factors[i] for i in kept), math.prod(factors[i] for i in rest)


def partial_trace(density, factors, region):
    """Restrict a state; trace and every record weight are preserved."""
    kept, rest, local, complement = _layout(factors, region)
    size = math.prod(factors)
    rho = _matrix(density, size, 'density')
    order = kept+rest
    tensor = rho.reshape(factors+factors).transpose(order+tuple(i+len(factors) for i in order))
    with np.errstate(over='ignore', invalid='ignore'):
        result = np.trace(tensor.reshape(local, complement, local, complement), axis1=1, axis2=3)
    if not np.isfinite(result).all():
        raise ValueError('restricted operator must remain finite')
    return result


def embed_operator(operator, factors, region):
    """The multiplicative inclusion A -> A tensor I in the declared order."""
    kept, rest, local, complement = _layout(factors, region)
    matrix = _matrix(operator, local)
    order = kept+rest
    joint = np.kron(matrix, np.eye(complement))
    tensor = joint.reshape(tuple(factors[i] for i in order)*2)
    inverse = tuple(order.index(i) for i in range(len(factors)))
    return tensor.transpose(inverse+tuple(i+len(factors) for i in inverse)).reshape(joint.shape)


def conditional_expectation(operator, factors, region):
    """Normalized partial trace followed by inclusion: CP, generally not multiplicative."""
    _, _, _, complement = _layout(factors, region)
    return embed_operator(partial_trace(operator, factors, region)/complement, factors, region)


@dataclass(frozen=True)
class CQView:
    """State on direct-sum_r M_dimension, with all committed records retained."""
    records: tuple[tuple[int, ...], ...]
    blocks: tuple[np.ndarray, ...]
    dimension: int

    def __post_init__(self):
        if type(self.dimension) is not int or self.dimension < 1:
            raise ValueError('a positive integer algebra dimension is required')
        if (not isinstance(self.records, tuple) or not self.records
                or any(not isinstance(record, tuple) or any(type(bit) is not int or bit < 0 for bit in record)
                       for record in self.records)
                or len(set(self.records)) != len(self.records)
                or len({len(record) for record in self.records}) != 1
                or not isinstance(self.blocks, tuple) or len(self.records) != len(self.blocks)):
            raise ValueError('distinct equally long records and every corresponding block are required')
        blocks = []
        for value in self.blocks:
            block = _matrix(value, self.dimension, 'density block')
            scale = float(np.max(np.abs(block)))
            mass = float(np.trace(block).real)
            if (np.linalg.norm(block-block.conj().T, 2) > TOLERANCE*scale
                    or np.linalg.eigvalsh((block+block.conj().T)/2)[0] < -TOLERANCE*max(scale, mass)
                    or mass < 0):
                raise ValueError('each unnormalized block must be positive and Hermitian')
            blocks.append(_immutable(block))
        if abs(math.fsum(float(np.trace(block).real) for block in blocks)-1) > TOLERANCE:
            raise ValueError('the complete observer state must have unit total mass')
        object.__setattr__(self, 'blocks', tuple(blocks))

    @property
    def probabilities(self):
        return dict(zip(self.records, (float(np.trace(block).real) for block in self.blocks), strict=True))


def observer_view(history, factors=None, region=None):
    """Restrict a validated execution, without granting reference access to its observer."""
    from .runtime import validate_history
    history = validate_history(history)
    if factors is None:
        factors = (history.work_dimension,)
    if region is None:
        region = tuple(range(len(factors)))
    _, _, local, _ = _layout(factors, region)
    if math.prod(factors) != history.work_dimension:
        raise ValueError('supplied factorization must describe the actual work register')
    work = tuple(partial_trace(branch.density,
                               (history.work_dimension, history.reference_dimension), (0,))
                 for branch in history.branches)
    return CQView(tuple(branch.record for branch in history.branches),
                  tuple(partial_trace(block, factors, region) for block in work), local)


def restrict(view, factors, region):
    """Keep the same classical record algebra and restrict its quantum blocks."""
    _, _, local, _ = _layout(factors, region)
    if math.prod(factors) != view.dimension:
        raise ValueError('factorization does not match the accessible algebra')
    return CQView(view.records, tuple(partial_trace(block, factors, region) for block in view.blocks), local)


def _observables(view, observable):
    if isinstance(observable, Mapping):
        if set(observable) != set(view.records):
            raise ValueError('one observable is required for every retained record')
        return tuple(_matrix(observable[record], view.dimension) for record in view.records)
    matrix = _matrix(observable, view.dimension)
    return (matrix,)*len(view.records)


def expectation(view, observable):
    """Positive functional sum_r Tr(rho_r A_r), with no branch normalization."""
    with np.errstate(over='ignore', invalid='ignore'):
        result = sum((np.trace(block@matrix) for block, matrix in
                      zip(view.blocks, _observables(view, observable), strict=True)), 0j)
    if not np.isfinite(result):
        raise ValueError('expectation must remain finite')
    return result


def public_expectation(view, values):
    """Read a function of committed records, including every zero-weight outcome."""
    if not isinstance(values, Mapping) or set(values) != set(view.records):
        raise ValueError('public readout must specify every retained record')
    scalars = {}
    for record, value in values.items():
        raw = np.asarray(value)
        if raw.shape != () or raw.dtype.kind not in 'iufc' or not np.isfinite(raw):
            raise ValueError('public values must be finite numeric scalars')
        scalars[record] = complex(value)
        if not (math.isfinite(scalars[record].real) and math.isfinite(scalars[record].imag)):
            raise ValueError('public values must remain finite in binary64')
    result = sum((probability*scalars[record] for record, probability in view.probabilities.items()), 0j)
    if not np.isfinite(result):
        raise ValueError('public expectation must remain finite')
    return result


def _normalized(block, mass):
    # Complex division may form 1/mass first and overflow for a valid tiny
    # branch. Componentwise real division preserves its normalized state.
    result = np.empty_like(block)
    result.real = block.real/mass
    result.imag = block.imag/mass
    return result


@dataclass(frozen=True)
class CQGNS:
    """Direct sum of the existing finite GNS spaces, weighted by actual masses."""
    view: CQView
    components: tuple[MatrixGNS | None, ...]

    @property
    def dimension(self):
        return sum(component.dimension for component in self.components if component is not None)

    def lift(self, observable):
        matrices = _observables(self.view, observable)
        values = [math.sqrt(mass)*component.lift(matrix).ravel()
                  for mass, component, matrix in zip(self.view.probabilities.values(), self.components,
                                                      matrices, strict=True) if component is not None]
        return np.concatenate(values)

    def inner(self, first, second):
        return np.vdot(self.lift(first), self.lift(second))

    def expectation(self, observable):
        return self.inner(np.eye(self.view.dimension), observable)


def gns(view):
    """Numerical GNS quotient; zero-mass blocks remain in the view with zero representation."""
    components = []
    for block, mass in zip(view.blocks, view.probabilities.values(), strict=True):
        if mass > 0:
            component = MatrixGNS.from_density(_normalized(block, mass))
            components.append(MatrixGNS(_immutable(component.density), _immutable(component.root_columns)))
        else:
            components.append(None)
    return CQGNS(view, tuple(components))


def gns_embedding(view, factors, region):
    """Return J, regional GNS and common GNS, with J[A]=[embed(A)].

    The coordinate construction also admits nonfaithful local states by using
    the support inverse. Independent numerical rank cuts can disagree between
    the common and regional states; the resulting map is rejected unless its
    isometry and cyclic-vector identities pass. These are numerical checks,
    not eigenvalue enclosures. The dense map has an explicit allocation budget.
    """
    regional = restrict(view, factors, region)
    local, common = gns(regional), gns(view)
    if local.dimension*common.dimension > MAX_MATRIX_ENTRIES:
        raise ValueError('explicit GNS inclusion exceeds the dense access budget')
    inclusion = np.zeros((common.dimension, local.dimension), complex)
    row_start = column_start = 0
    for small, large in zip(local.components, common.components, strict=True):
        if small is None:
            continue
        inverse = np.linalg.pinv(small.root_columns)
        for column in range(small.dimension):
            coordinate = np.zeros_like(small.root_columns)
            coordinate.flat[column] = 1
            observable = coordinate@inverse
            inclusion[row_start:row_start+large.dimension, column_start+column] = large.lift(
                embed_operator(observable, factors, region)).ravel()
        row_start += large.dimension
        column_start += small.dimension
    gram = inclusion.conj().T@inclusion
    cyclic_error = np.linalg.norm(inclusion@local.lift(np.eye(regional.dimension))
                                  -common.lift(np.eye(view.dimension)))
    if (not np.isfinite(inclusion).all() or not np.isfinite(gram).all()
            or not np.isfinite(cyclic_error)
            or np.linalg.norm(gram-np.eye(local.dimension), 2) > TOLERANCE
            or cyclic_error > TOLERANCE):
        raise ValueError('numerical GNS support cuts do not give a compatible regional isometry')
    return inclusion, local, common


def modular_flow(view, observable, parameter):
    """Block modular automorphism for faithful blocks; parameter is dimensionless.

    Zero-weight records make the state nonfaithful on this full algebra and
    are rejected.  No reference-state padding or physical clock is inserted.
    """
    raw = np.asarray(parameter)
    if raw.shape != () or raw.dtype.kind not in 'iuf' or not np.isfinite(raw):
        raise ValueError('a finite real modular parameter is required')
    matrices = _observables(view, observable)
    result = {}
    for record, block, matrix in zip(view.records, view.blocks, matrices, strict=True):
        mass = view.probabilities[record]
        if mass <= 0:
            raise ValueError('modular flow on the full record algebra requires a faithful state')
        result[record] = _modular_flow(_normalized(block, mass), matrix, float(parameter))
    return result


def forecast(view, kraus):
    """Joint old-record/new-outcome law for a complete instrument on this algebra.

    A restricted view supports only instruments on its accessible factors.
    Interacting futures require the complete work view and their full maps;
    this function never supplies a missing environment or assignment map.
    """
    if not isinstance(kraus, tuple) or not kraus or any(not isinstance(group, tuple) or not group for group in kraus):
        raise ValueError('supply every outcome and its nonempty Kraus tuple')
    groups = tuple(tuple(_matrix(matrix, view.dimension, 'Kraus map') for matrix in group) for group in kraus)
    with np.errstate(over='ignore', invalid='ignore'):
        effects = tuple(sum((matrix.conj().T@matrix for matrix in group),
                            np.zeros((view.dimension, view.dimension), complex)) for group in groups)
        total = sum(effects)
    if not np.isfinite(total).all() or any(not np.isfinite(effect).all() for effect in effects):
        raise ValueError('Kraus effects must remain finite')
    if np.linalg.norm(total-np.eye(view.dimension), 2) > TOLERANCE:
        raise ValueError('forecast requires a complete trace-preserving instrument')
    return {(*record, outcome): float(np.trace(block@effect).real)
            for record, block in zip(view.records, view.blocks, strict=True)
            for outcome, effect in enumerate(effects)}
