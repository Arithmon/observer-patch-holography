"""Finite operational attachment to twelve primitive central address ports.

States are twelve subnormalized d-by-d blocks, not a dephased native state.
The logical algebra M_d embeds diagonally in direct_sum_12 M_d. Address
refresh and star repair reproduce the supplied uniform seam channel exactly;
no native or private quantum factor is reset. This finite adapter is not a
full A1 source, a response-completeness proof, or a physical clock selection.
Floating matrix-product evaluators conservatively refuse subnormal operands;
their ordinary roundoff is not enclosed by an error certificate.
"""
from __future__ import annotations

from collections.abc import Mapping
from numbers import Integral

import numpy as np

from quantum_information.states import ATOL, density_matrix

from .source import checked_edge, edges, vertices


def _numeric(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} must not be masked")
    try:
        raw = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric array") from exc
    if raw.dtype.kind not in "iufc" or not np.all(np.isfinite(raw)):
        raise ValueError(f"{name} must be finite and numeric")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        numeric = raw.astype(complex)
    if (np.any((raw.real != 0) & (numeric.real == 0))
            or np.any((raw.imag != 0) & (numeric.imag == 0))):
        raise ValueError(f"{name} underflows the numerical representation")
    if (not np.all(np.isfinite(numeric))
            or np.max(np.abs(numeric.real), initial=0) > 1 + ATOL
            or np.max(np.abs(numeric.imag), initial=0) > 1 + ATOL):
        raise ValueError(f"{name} has entries outside its normalized domain")
    # Entries of states, projectors and CPTP Kraus operators have modulus <=1.
    if np.max(np.abs(numeric), initial=0) > 1 + ATOL:
        raise ValueError(f"{name} has entries outside its normalized domain")
    return numeric


def _port(port):
    if (isinstance(port, (bool, np.bool_)) or not isinstance(port, Integral)
            or not 0 <= port < 12):
        raise ValueError("port must be an integer from 0 to 11")
    return int(port)


def _matrix_operand(value, name):
    """Refuse subnormal operands before BLAS can silently erase them.

    This is a conservative floating boundary, not an error certificate for
    arbitrary matrix products. Direct attachment and address moves keep
    their separate component-preservation checks.
    """
    tiny = np.finfo(float).tiny
    for component in (value.real, value.imag):
        magnitude = np.abs(component)
        if np.any((magnitude != 0) & (magnitude < tiny)):
            raise ValueError(f"{name} has subnormal matrix-product components")


def antipode(port):
    """Source-chart antipode, used for contextual binary outcome storage."""
    p = _port(port)
    frame = vertices()
    return int(np.argmin(np.linalg.norm(frame + frame[p], axis=1)))


def cq_state(value):
    """Validate twelve PSD blocks with total trace one; never renormalize.

    Positivity uses the shared numerical state tolerance, applied separately
    to every populated conditional block. These are floating diagnostics.
    """
    blocks = _numeric(value, "cq state")
    if (blocks.ndim != 3 or blocks.shape[0] != 12 or not blocks.shape[1]
            or blocks.shape[1] != blocks.shape[2]):
        raise ValueError("cq state needs shape (12, d, d), d positive")
    weights = np.trace(blocks, axis1=1, axis2=2)
    if abs(weights.sum() - 1) > ATOL:
        raise ValueError("cq state must have total trace one")
    for index, weight in enumerate(weights):
        if weight.real < 0 or abs(weight.imag) > ATOL:
            raise ValueError("cq block must have nonnegative real trace")
        if weight.real == 0:
            if np.any(blocks[index] != 0):
                raise ValueError("zero-trace cq block must be zero")
        else:
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                conditional = blocks[index] / weight.real
            conditional = _numeric(conditional, "conditional cq state")
            blocks[index] = weight.real * density_matrix(conditional)
    return blocks


def decode(blocks):
    """Forget only the central address, retaining every native coherence."""
    return cq_state(blocks).sum(axis=0)


def encode(rho):
    """Uniform attachment rho -> (rho/12,...,rho/12); reject erased entries.

    The exact map is recoverable. Floating evaluation is not a certificate
    and rejects component underflow instead of silently deleting input data.
    """
    state = density_matrix(_numeric(rho, "native state"))
    with np.errstate(under="ignore"):
        block = state / 12
    if (np.any((state.real != 0) & (block.real == 0))
            or np.any((state.imag != 0) & (block.imag == 0))):
        raise ValueError("uniform attachment underflows a native state component")
    return np.repeat(block[None, :, :], 12, axis=0)


def refresh_address(blocks):
    """Replace the address by independent uniform randomness only."""
    return encode(decode(blocks))


def select_address(blocks, port):
    """Set a requested central address without changing the native marginal."""
    p = _port(port)
    state = decode(blocks)
    result = np.zeros((12, *state.shape), dtype=complex)
    result[p] = state
    return result


class SeamChannels:
    """A fixed CPTP Kraus family for every source seam, copied on entry.

    Supplying these thirty maps is an explicit input. Kraus validation does
    not derive their physical realization or the response grammar.
    """

    def __init__(self, families):
        if not isinstance(families, Mapping) or set(families) != set(edges()):
            raise ValueError("one Kraus family per canonical source seam required")
        for edge in families:
            if not isinstance(edge, tuple) or checked_edge(edge) != edge:
                raise ValueError("canonical integer seam keys required")
        stored, dimension = {}, None
        for edge in edges():
            kraus = _numeric(families[edge], "Kraus family")
            if (kraus.ndim != 3 or not kraus.shape[0] or not kraus.shape[1]
                    or kraus.shape[1] != kraus.shape[2]):
                raise ValueError("Kraus family needs shape (k, d, d)")
            if dimension is None:
                dimension = kraus.shape[1]
            if kraus.shape[1] != dimension:
                raise ValueError("all seam channels must have the same dimension")
            _matrix_operand(kraus, "Kraus family")
            complete = sum(k.conj().T @ k for k in kraus)
            if np.linalg.norm(complete - np.eye(dimension), ord="fro") > ATOL:
                raise ValueError("seam Kraus family must be trace preserving")
            kraus.setflags(write=False)
            stored[edge] = kraus
        self._families = stored
        self._dimension = dimension

    @property
    def dimension(self):
        return self._dimension

    def _apply(self, edge, value):
        _matrix_operand(value, "channel input")
        return sum(k @ value @ k.conj().T for k in self._families[edge])

    def apply(self, edge, rho):
        """Apply one supplied seam to a normalized native state."""
        seam = tuple(sorted(checked_edge(edge)))
        state = density_matrix(_numeric(rho, "native state"))
        if state.shape != (self.dimension, self.dimension):
            raise ValueError("native state and seam dimensions differ")
        return self._apply(seam, state)


def star_repair(blocks, channels):
    """The actual address p selects a uniform one-of-five incident seam.

    Without a fresh independent address this generally differs from the
    global uniform seam mixture, including on correlated cq inputs.
    """
    state = cq_state(blocks)
    if not isinstance(channels, SeamChannels):
        raise ValueError("validated SeamChannels required")
    if state.shape[1] != channels.dimension:
        raise ValueError("cq state and seam dimensions differ")
    result = np.zeros_like(state)
    for a, b in edges():
        result[a] += channels._apply((a, b), state[a])
        result[b] += channels._apply((a, b), state[b])
    with np.errstate(under="ignore"):
        average = result / 5
    if (np.any((result.real != 0) & (average.real == 0))
            or np.any((result.imag != 0) & (average.imag == 0))):
        raise ValueError("star averaging underflows a state component")
    return average


def tick(blocks, channels):
    """One native global attempt: fresh address, then addressed star repair.

    decode(tick(sigma)) = (sum_e Phi_e/30)(decode(sigma)) for every cq state.
    This intertwining iterates without intermediate native measurements.
    """
    return star_repair(refresh_address(blocks), channels)


def binary_query(blocks, occupied_projectors):
    """Read the native binary projector selected by the central address.

    Occupied writes p, empty writes antipode(p), retaining the full native
    Lueders branch. Decoding the outcome requires the previously retained
    query address. The center alone does not store both all twelve possible
    query addresses and their two outcomes simultaneously.
    """
    state = cq_state(blocks)
    projectors = _numeric(occupied_projectors, "binary projectors")
    if projectors.shape != state.shape:
        raise ValueError("binary projectors need shape (12, d, d)")
    _matrix_operand(state, "query input")
    _matrix_operand(projectors, "binary projectors")
    identity = np.eye(state.shape[1])
    result = np.zeros_like(state)
    for p, projector in enumerate(projectors):
        if (np.linalg.norm(projector - projector.conj().T, ord="fro") > ATOL
                or np.linalg.norm(projector @ projector - projector, ord="fro") > ATOL):
            raise ValueError("binary query requires orthogonal projectors")
        complement = identity - projector
        result[p] += projector @ state[p] @ projector
        result[antipode(p)] += complement @ state[p] @ complement
    return result
