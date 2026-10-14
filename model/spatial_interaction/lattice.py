"""A supplied cubic Bose field on its exact two-particle invariant sector.

The spatial domain is the open cube {-R,...,R}^3 at fixed lattice spacing.
The finite Hamiltonian is the graph-Laplacian kinetic energy plus repulsive
onsite contact energy.  The common 2*m phase is omitted throughout this
fixed-number execution.  No spatial continuum or relativistic claim follows.

Local two-site operators use the same sixteen-dimensional occupation code
as the source compiler: index = n_x + 4*n_y.  Only the sectors containing
zero, one and two particles are accessed; their exact invariant compression
is not a truncation of the interacting dynamics.  The rest of the global
two-particle state is retained, including when a local Kraus map has a
nontrivial vacuum amplitude.
"""
from __future__ import annotations

from itertools import product
from math import isfinite, sqrt
from numbers import Integral, Real

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


LOCAL_SECTORS = ((0,), (1, 4), (2, 5, 8))
SECTOR_TOLERANCE = 2e-8


def _real(value, name, *, nonnegative=False):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not isfinite(value) or (nonnegative and value < 0)):
        raise ValueError(f"{name} must be a finite {'nonnegative ' if nonnegative else ''}real")
    return float(value)


def local_sector_leakage(gate):
    """Maximum entry leaving an accessed local number sector.

    Numerical source replay can leave roundoff-sized forbidden entries.  A
    separate numerical error account must cover their compression; this
    diagnostic is not a mathematical replacement for exact conservation.
    """
    if np.ma.isMaskedArray(gate):
        raise ValueError("a local operator cannot be masked")
    raw = np.asarray(gate)
    if raw.shape != (16, 16) or raw.dtype.kind not in "fiuc" or not np.isfinite(raw).all():
        raise ValueError("a finite numeric sixteen-dimensional local operator is required")
    return max(float(np.max(abs(raw[np.ix_([i for i in range(16) if i not in sector],
                                          sector)]))) for sector in LOCAL_SECTORS)


class CubicBosonLattice:
    """Exact N=2 Fock basis with six nearest-neighbor edge matchings.

    ``basis[k]=(i,j)`` means one particle at each site when i<j, and a
    normalized double-occupation state when i=j.  Site coordinates are
    supplied; their index order is lexicographic.  Every undirected edge
    appears once in ``edge_groups[2*axis+parity]``.  Parity is taken in the
    absolute integer coordinate, so the schedule is consistent across boxes.
    """

    def __init__(self, radius, hopping=1.0, interaction=1.0):
        if (isinstance(radius, (bool, np.bool_)) or not isinstance(radius, Integral)
                or radius < 1):
            raise ValueError("the cubic radius must be a positive integer")
        self.radius = int(radius)
        self.hopping = _real(hopping, "hopping", nonnegative=True)
        self.interaction = _real(interaction, "interaction", nonnegative=True)
        self.sites = tuple(product(range(-self.radius, self.radius+1), repeat=3))
        self.site_index = {site: i for i, site in enumerate(self.sites)}
        self.volume = len(self.sites)
        first, second = np.triu_indices(self.volume)
        self._first, self._second = first, second
        self.basis = tuple(zip(map(int, first), map(int, second)))
        self.index = {pair: i for i, pair in enumerate(self.basis)}
        self.dimension = len(self.basis)
        self.origin = self.site_index[(0, 0, 0)]
        self.neighbor = self.site_index[(1, 0, 0)]
        self.ramsey_edge = (self.origin, self.neighbor)
        self.collision_mask = first == second
        groups = [[] for _ in range(6)]
        degree = np.zeros(self.volume, dtype=int)
        for i, site in enumerate(self.sites):
            for axis in range(3):
                neighbor = list(site)
                neighbor[axis] += 1
                j = self.site_index.get(tuple(neighbor))
                if j is not None:
                    groups[2*axis+(site[axis] % 2)].append((i, j))
                    degree[i] += 1
                    degree[j] += 1
        self.edge_groups = tuple(tuple(group) for group in groups)
        self.edges = tuple(edge for group in self.edge_groups for edge in group)
        self.degree = degree
        self._local_indices = {}
        self.hamiltonian = self._hamiltonian()

    def _site(self, value):
        if isinstance(value, (bool, np.bool_)):
            raise ValueError("a site index or cubic coordinate is required")
        if isinstance(value, Integral):
            if 0 <= value < self.volume:
                return int(value)
            raise ValueError("site index lies outside the cube")
        if (isinstance(value, (tuple, list)) and len(value) == 3
                and all(isinstance(x, Integral) and not isinstance(x, (bool, np.bool_))
                        for x in value)):
            try:
                return self.site_index[tuple(map(int, value))]
            except KeyError as exc:
                raise ValueError("site coordinate lies outside the cube") from exc
        raise ValueError("a site index or cubic coordinate is required")

    def state_index(self, first, second):
        """Return the normalized Fock-basis index for two sites."""
        return self.index[tuple(sorted((self._site(first), self._site(second))))]

    def _edge(self, edge):
        if not isinstance(edge, (tuple, list)) or len(edge) != 2:
            raise ValueError("two distinct local sites are required")
        x, y = map(self._site, edge)
        if x == y:
            raise ValueError("two distinct local sites are required")
        return x, y

    def _indices(self, edge):
        x, y = self._edge(edge)
        if (x, y) not in self._local_indices:
            spectators = [z for z in range(self.volume) if z not in (x, y)]
            one = np.array([[self.index[tuple(sorted((x, z)))] for z in spectators],
                            [self.index[tuple(sorted((y, z)))] for z in spectators]], dtype=int)
            two = np.array([self.index[(x, x)], self.index[tuple(sorted((x, y)))],
                            self.index[(y, y)]], dtype=int)
            self._local_indices[(x, y)] = (one, two)
        return self._local_indices[(x, y)]

    def _hamiltonian(self):
        diagonal = self.hopping*(self.degree[self._first]+self.degree[self._second])
        diagonal = diagonal + self.interaction*self.collision_mask
        rows, cols, values = [np.arange(self.dimension)], [np.arange(self.dimension)], [diagonal]
        for edge in self.edges:
            one, two = self._indices(edge)
            rows.extend((one[0], one[1], two[[0, 1, 1, 2]]))
            cols.extend((one[1], one[0], two[[1, 0, 2, 1]]))
            values.extend((np.full(self.volume-2, -self.hopping),
                           np.full(self.volume-2, -self.hopping),
                           np.full(4, -sqrt(2)*self.hopping)))
        return sparse.coo_matrix((np.concatenate(values),
                                  (np.concatenate(rows), np.concatenate(cols))),
                                 shape=(self.dimension, self.dimension)).tocsr()

    def _state(self, state):
        if np.ma.isMaskedArray(state):
            raise ValueError("a state cannot be masked")
        raw = np.asarray(state)
        if (raw.ndim not in (1, 2) or raw.shape[0] != self.dimension
                or raw.dtype.kind not in "fiuc" or not np.isfinite(raw).all()):
            raise ValueError("finite amplitude columns on the declared two-particle basis are required")
        return np.asarray(raw, dtype=complex)

    @staticmethod
    def _gate(gate):
        leakage = local_sector_leakage(gate)
        if leakage > SECTOR_TOLERANCE:
            raise ValueError("the local operator must preserve each accessed particle-number sector")
        return np.asarray(gate, dtype=complex)

    def _apply_local_inplace(self, state, edge, gate):
        one, two = self._indices(edge)
        # Advanced indexing copies each affected block before any overwrite.
        # The same rule works for vector states and spectator/reference columns.
        block = state[one]
        transformed_one = np.einsum("ij,j...->i...", gate[np.ix_((1, 4), (1, 4))], block)
        transformed_two = gate[np.ix_((2, 5, 8), (2, 5, 8))] @ state[two]
        if gate[0, 0] != 1:
            state *= gate[0, 0]
        state[one] = transformed_one
        state[two] = transformed_two

    def apply_local(self, state, edge, gate):
        """Apply a local number-preserving operator, including a Kraus map.

        No normalization or outcome selection is performed.  Tiny forbidden
        numerical entries are discarded only within SECTOR_TOLERANCE; callers
        retain their source-replay residual in the numerical error account.
        """
        result = self._state(state).copy()
        gate = self._gate(gate)
        self._apply_local_inplace(result, edge, gate)
        return result

    def apply_matching(self, state, group, gate):
        """Apply one disjoint edge matching with a single full-state copy."""
        if isinstance(group, Integral) and not isinstance(group, (bool, np.bool_)):
            if not 0 <= group < 6:
                raise ValueError("matching index must lie from zero through five")
            group = self.edge_groups[int(group)]
        if not isinstance(group, (tuple, list)):
            raise ValueError("a matching index or a sequence of disjoint edges is required")
        edges = tuple(self._edge(edge) for edge in group)
        used = [site for edge in edges for site in edge]
        if len(set(used)) != len(used):
            raise ValueError("matching edges must be disjoint")
        result = self._state(state).copy()
        gate = self._gate(gate)
        # Source roundoff may give the vacuum a tiny phase.  Pulling that
        # scalar outside the complete matching is an algebraic factorization,
        # not a state normalization or replacement by one.  The conservative
        # range avoids magnifying tiny Kraus vacuum factors in floating point.
        vacuum = gate[0, 0]
        factor_vacuum = vacuum != 1 and .5 <= abs(vacuum) <= 2
        if factor_vacuum:
            gate = gate/vacuum
            gate[0, 0] = 1
        for edge in edges:
            self._apply_local_inplace(result, edge, gate)
        if factor_vacuum:
            result *= vacuum**len(edges)
        return result

    def doublon_state(self, site=None):
        """Two particles at one supplied site, normalized in the Fock basis."""
        site = self.origin if site is None else self._site(site)
        state = np.zeros(self.dimension, complex)
        state[self.index[(site, site)]] = 1
        return state

    def initial_state(self):
        """Local positive Ramsey superposition of doublon and split pair."""
        state = self.doublon_state()/sqrt(2)
        state[self.state_index(self.origin, self.neighbor)] = 1/sqrt(2)
        return state

    def occupation(self, site=None):
        site = self.origin if site is None else self._site(site)
        return (self._first == site).astype(int)+(self._second == site)

    def origin_occupation_mask(self, threshold=1):
        if (isinstance(threshold, (bool, np.bool_)) or not isinstance(threshold, Integral)
                or threshold not in (1, 2)):
            raise ValueError("occupation threshold must be one or two")
        return self.occupation() >= threshold

    def phase_feedback(self, state, angle, site=None):
        angle = _real(angle, "local feedback angle")
        state = self._state(state)
        phase = np.exp(1j*angle*self.occupation(site))
        return phase*state if state.ndim == 1 else phase[:, None]*state

    def ramsey_branches(self, state, sign=1):
        """Return (no, yes) for a local +/-i Ramsey projector, unnormalized."""
        if isinstance(sign, (bool, np.bool_)) or not isinstance(sign, Integral) or sign not in (-1, 1):
            raise ValueError("the Ramsey quadrature sign must be minus or plus one")
        state = self._state(state)
        detector = self.doublon_state()/sqrt(2)
        detector[self.state_index(self.origin, self.neighbor)] = sign*1j/sqrt(2)
        amplitude = detector.conj() @ state
        yes = detector*amplitude if state.ndim == 1 else detector[:, None]*amplitude[None, :]
        return state-yes, yes

    def reference_evolve(self, state, time):
        """Independent sparse exponential of the full interacting generator.

        This floating numerical oracle is separate from the executed local
        product word; it is not an interval certificate or an OPH derivation.
        The parameter is the supplied Hamiltonian's model-coordinate time.
        """
        time = _real(time, "model time")
        state = self._state(state)
        return expm_multiply((-1j*time)*self.hamiltonian, state)
