"""Finite block expectations and an exact boundary-hopping certificate.

The equivalence theorem is proved in README.md. The block routines are
numerical controls of that theorem; they do not decide exact factorization
of an arbitrary floating density. The field witness uses Gaussian-integer
entries and Python integer accumulation, without diagonalizing a Gibbs state.
"""
from fractions import Fraction as F
from math import prod
import numpy as np
from ..interacting_field.lattice import WeylLattice
from ..field_stress.energy import GaussianInteger


def _dimensions(dimensions):
    if not isinstance(dimensions, (tuple, list)) or not dimensions:
        raise ValueError('a nonempty list of local/exterior dimension pairs is required')
    result = []
    for pair in dimensions:
        if (not isinstance(pair, (tuple, list)) or len(pair) != 2
                or any(type(n) is not int or n < 1 for n in pair)):
            raise ValueError('positive integer local/exterior dimensions are required')
        result.append(tuple(pair))
    return tuple(result)


def _matrix(value, size):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked matrices are not admitted')
    value = np.asarray(value)
    if (value.shape != (size, size) or value.dtype.kind not in 'ficu'
            or not np.isfinite(value).all()):
        raise ValueError('a finite matrix on the full declared domain is required')
    return np.asarray(value, complex)


def block_expectation(matrix, dimensions, exterior_states):
    """Heisenberg UCP retraction onto direct_sum_n (M_a[n] tensor I_b[n]).

    ``exterior_states`` are declared faithful normalized density matrices.
    Their Hermiticity and unit trace are checked literally, with no clipping
    or normalization. Positive eigenvalues are floating diagnostics only.
    The exact theorem and complete Kraus identity concern exact states.
    """
    dimensions = _dimensions(dimensions)
    if not isinstance(exterior_states, (tuple, list)) or len(exterior_states) != len(dimensions):
        raise ValueError('one exterior state per nonzero sector is required')
    size = sum(prod(pair) for pair in dimensions)
    matrix = _matrix(matrix, size)
    result = np.zeros_like(matrix)
    start = 0
    for (a, b), supplied in zip(dimensions, exterior_states, strict=True):
        tau = _matrix(supplied, b)
        if (not np.array_equal(tau, tau.conj().T) or np.trace(tau) != 1
                or np.linalg.eigvalsh(tau)[0] <= 0):
            raise ValueError('a faithful Hermitian unit-trace exterior state is required')
        stop = start+a*b
        block = matrix[start:stop, start:stop].reshape(a, b, a, b)
        local = np.einsum('aibj,ji->ab', block, tau)
        result[start:stop, start:stop] = np.kron(local, np.eye(b))
        start = stop
    return result


def crossing_certificate(lattice=None, kept_modes=None):
    """Exact ||[H,N_U]||_HS^2 on the actual canonical full N=2 field.

    ``ordered_crossing_entries`` counts both matrix orientations p,q. With
    m modes, the exterior-square identity multiplies the one-body squared
    norm by m-2. The contact action commutes with N_U. The witness excludes
    a Gibbs-preserving UCP retraction only for beta != 0 on the stated local
    number-preserving algebra; it does not invalidate ordinary state traces.
    """
    lattice = WeylLattice(1) if lattice is None else lattice
    if lattice.speed != 1/8 or lattice.interaction not in (0., 1.):
        raise ValueError('this integer certificate covers c=1/8 and U=0 or 1')
    if kept_modes is None:
        kept_modes = (lattice.origin, lattice.origin+1,
                      lattice.neighbor, lattice.neighbor+1)
    if (not isinstance(kept_modes, (tuple, list))
            or any(type(p) is not int or not 0 <= p < lattice.modes for p in kept_modes)
            or tuple(sorted(set(kept_modes))) != tuple(kept_modes)):
        raise ValueError('sorted distinct retained modes are required')
    inside = np.zeros(lattice.modes, dtype=np.int64)
    inside[list(kept_modes)] = 1
    occupation = inside[lattice.first]+inside[lattice.second]
    one = GaussianInteger.from_exact(16*lattice.single_h)
    full = GaussianInteger.from_exact(16*lattice.hamiltonian)

    def squared_commutator(operator, diagonal):
        total = 0
        for part in (operator.real, operator.imag):
            part = part.tocoo()
            total += sum((int(z)*(int(diagonal[row])-int(diagonal[col])))**2
                         for row, col, z in zip(part.row, part.col, part.data, strict=True))
        return total

    one_squared = squared_commutator(one, inside)
    full_squared = squared_commutator(full, occupation)
    if full_squared != (lattice.modes-2)*one_squared:
        raise ArithmeticError('actual N=2 action fails the boundary-current identity')
    h = lattice.single_h.tocoo()
    count = sum(bool(inside[p] != inside[q]) and z != 0
                for p, q, z in zip(h.row, h.col, h.data, strict=True))
    return {
        'modes': lattice.modes,
        'N2_dimension': lattice.dimension,
        'kept_modes': list(kept_modes),
        'ordered_crossing_entries': int(count),
        'one_body_commutator_HS_squared': str(F(one_squared, 16**2)),
        'full_N2_commutator_HS_squared': str(F(full_squared, 16**2)),
        'exterior_square_multiplicity': lattice.modes-2,
        'exact_integer_witness': True,
        'Gibbs_preserving_regional_retraction_excluded_for_beta_nonzero': full_squared != 0,
        'ordinary_regional_CP_restriction_excluded': False,
        'physical_clock_mass_or_continuum_claimed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(crossing_certificate(), indent=2, sort_keys=True))
