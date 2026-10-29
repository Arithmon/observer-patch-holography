"""Peierls coupling of the existing Weyl/contact action, with all N=2 states.

The link field is an explicitly supplied external U(1) background. This
constructs neither dynamical gauge degrees of freedom nor physical charge.
Durations are Hamiltonian model coordinates with the convention hbar=1.
"""
from itertools import product
from numbers import Real

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from ..interacting_field.lattice import SIGMA, TwoFermions, real

MAX_RESOLVENT_DIMENSION = 8192


def _angle(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError('finite real link or site angle required')
    result = float(value)
    if not np.isfinite(result):
        raise ValueError('finite real link or site angle required')
    return result


class GaugedField(TwoFermions):
    """Finite induced cubic-site action; absent links are not boundary wraps.

    Site order is canonical lexicographic order. A key (x,y), x<y, denotes
    h_xy=e^(i A_xy) i c sigma_axis/2, the coefficient of c_x* c_y.
    Omitted nearest-neighbor keys have zero angle. All supplied keys must
    actually be bonds of the declared induced site set.
    """
    def __init__(self, sites=None, *, angles=None, speed=.125, interaction=1.):
        if sites is None:
            sites = tuple(product(range(-1, 2), repeat=3))
        if not isinstance(sites, (tuple, list)) or not 2 <= len(sites) <= 350:
            raise ValueError('2..350 distinct cubic sites required')
        if any(not isinstance(x, (tuple, list)) or len(x) != 3 or
               any(type(c) is not int for c in x) for x in sites):
            raise ValueError('three integer coordinates per site required')
        self.sites = tuple(sorted(tuple(x) for x in sites))
        if len(set(self.sites)) != len(self.sites):
            raise ValueError('duplicate site')
        self.site_index = {x:i for i,x in enumerate(self.sites)}
        self.volume = len(self.sites)
        self.speed = real(speed, 'coordinate speed', nonnegative=True)
        self.interaction = real(interaction, 'interaction', nonnegative=True)
        super().__init__(2*self.volume)
        edges = {}
        for x, coord in enumerate(self.sites):
            for axis in range(3):
                other = list(coord); other[axis] += 1
                y = self.site_index.get(tuple(other))
                if y is not None:
                    edges[x,y] = axis
        if angles is None:
            angles = {}
        if type(angles) is not dict:
            raise ValueError('bond-angle dictionary required')
        for edge, angle in angles.items():
            if (type(edge) is not tuple or len(edge) != 2 or
                    any(type(i) is not int for i in edge) or edge not in edges):
                raise ValueError('angle key is not an oriented field bond')
            _angle(angle)
        self.axes = edges
        self.angles = {edge:_angle(angles.get(edge, 0.)) for edge in edges}
        groups = [[] for _ in range(6)]
        rows, cols, values = [], [], []
        for (x,y), axis in edges.items():
            matrix = np.exp(1j*self.angles[x,y])*1j*self.speed*SIGMA[axis]/2
            for a,b in product(range(2), repeat=2):
                z = complex(matrix[a,b])
                if z:
                    p,q = 2*x+a,2*y+b
                    groups[2*axis+self.sites[x][axis] % 2].append((p,q,z))
                    rows.extend((p,q)); cols.extend((q,p)); values.extend((z,z.conjugate()))
        self.groups = tuple(tuple(g) for g in groups)
        self.single_h = sparse.coo_matrix((values,(rows,cols)),
            shape=(self.modes,self.modes),dtype=complex).tocsr()
        self.kinetic = self.lift_one_body(self.single_h)
        self.contact_indices = np.array([self.index(2*x,2*x+1) for x in range(self.volume)])
        contact = np.zeros(self.dimension); contact[self.contact_indices] = self.interaction
        self.hamiltonian = self.kinetic+sparse.diags(contact)

    def local_numbers(self):
        return tuple(sparse.diags((self.first//2 == x).astype(int)+
                                  (self.second//2 == x).astype(int),dtype=complex)
                     for x in range(self.volume))

    def current(self, edge):
        if type(edge) is not tuple or len(edge) != 2 or any(type(x) is not int for x in edge) or edge not in self.axes:
            raise ValueError('oriented field bond required')
        x,y = edge
        h = self.single_h.tocoo()
        take = (h.row//2 == x)&(h.col//2 == y)
        one = sparse.coo_matrix((h.data[take],(h.row[take],h.col[take])),shape=h.shape)
        forward = self.lift_one_body(one)
        return (1j*(forward-forward.conj().T)).tocsr()

    def gauge_transform(self, site_angles):
        if (not isinstance(site_angles, (tuple,list,np.ndarray)) or
                np.ndim(site_angles) != 1 or len(site_angles) != self.volume):
            raise ValueError('one angle per site required')
        alpha = np.array([_angle(a) for a in site_angles])
        changed = type(self)(self.sites, speed=self.speed, interaction=self.interaction,
            angles={(x,y):a+alpha[x]-alpha[y] for (x,y),a in self.angles.items()})
        phases = np.exp(1j*(alpha[self.first//2]+alpha[self.second//2]))
        return changed, phases

    def kinetic_norm_upper(self):
        """Schur bound 2 max_p sum_q |h_pq|, valid even at nonzero flux.

        On an induced cubic site set it is at most 6c. The sharper flat-field
        anticommuting-symbol bound is not assumed for an arbitrary link field.
        """
        return 2*float(np.max(np.asarray(abs(self.single_h).sum(axis=1))))


def contact_blocks(field):
    """P projects onto every onsite pair, Q retains every separated pair."""
    p = field.contact_indices
    q = np.setdiff1d(np.arange(field.dimension),p)
    if field.kinetic[p][:,p].nnz:
        raise ArithmeticError('one-body hopping connects two contact states')
    return p,q,field.kinetic[p][:,q].tocsr(),field.kinetic[q][:,q].tocsr()


def second_order_pair(field):
    """U I + PTQ QTP/U, a leading approximation, not projected exact dynamics."""
    if field.interaction <= 0:
        raise ValueError('positive contact interaction required')
    _,_,b,_ = contact_blocks(field)
    return (field.interaction*sparse.eye(field.volume,dtype=complex)+(b@b.conj().T)/field.interaction).tocsr()


def determinant_pair_action(field):
    """Independent local formula for the same second-order matrix.

    Off-diagonal pair transfer is 2 det(h_xy)/U. The diagonal contribution
    of each incident bond is Tr(h_xy* h_xy)/U.
    """
    if field.interaction <= 0:
        raise ValueError('positive contact interaction required')
    matrix = np.eye(field.volume,dtype=complex)*field.interaction
    for x,y in field.axes:
        h = field.single_h[2*x:2*x+2,2*y:2*y+2].toarray()
        value = 2*(h[0,0]*h[1,1]-h[0,1]*h[1,0])/field.interaction
        matrix[x,y] += value; matrix[y,x] += value.conjugate()
        shift = float(np.vdot(h,h).real)/field.interaction
        matrix[x,x] += shift; matrix[y,y] += shift
    return matrix


def feshbach(field, energy):
    """Exact contact Schur complement above a proved Q kinetic bound.

    The bounded numerical evaluator retains every Q state and solves its
    resolvent. The analytic identity does not require materializing a matrix.
    """
    if field.dimension > MAX_RESOLVENT_DIMENSION:
        raise ValueError('field exceeds the host resolvent dimension limit')
    energy = _angle(energy)
    bound = field.kinetic_norm_upper()
    if energy <= bound:
        raise ValueError('energy must exceed the declared kinetic norm bound')
    _,_,b,d = contact_blocks(field)
    if not d.shape[0]:
        return np.eye(field.volume,dtype=complex)*field.interaction
    solved = splu((energy*sparse.eye(d.shape[0],dtype=complex)-d).tocsc()).solve(b.conj().T.toarray())
    return np.eye(field.volume,dtype=complex)*field.interaction+b@solved


def pair_current_response(field, energy, edge):
    """Return fixed-E link derivative and the full-state norm metric.

    For a contact eigenvector u, its separated component is X u with
    X=(E-QTQ)^(-1)QTP. Thus its physical full-sector norm is u* M u,
    M=I+X*X. The exact Hellmann--Feynman current is (u* dF u)/(u* M u),
    not u* dF u divided only by the contact norm.
    """
    if field.dimension > MAX_RESOLVENT_DIMENSION:
        raise ValueError('field exceeds the host resolvent dimension limit')
    energy = _angle(energy)
    if energy <= field.kinetic_norm_upper():
        raise ValueError('energy must exceed the declared kinetic norm bound')
    p,q,b,d = contact_blocks(field)
    x = splu((energy*sparse.eye(len(q),dtype=complex)-d).tocsc()).solve(b.conj().T.toarray())
    current = field.current(edge)
    cross = current[p][:,q]@x
    derivative = current[p][:,p].toarray()+cross+cross.conj().T+x.conj().T@(current[q][:,q]@x)
    metric = np.eye(len(p),dtype=complex)+x.conj().T@x
    return derivative,metric,x


def feshbach_remainder_bound(field, energy):
    """Uniform operator bound against U I+BB*/U for E>||QTQ|| upper.

    BB* has its own computable Schur bound. This compares operators at the
    declared energy; it does not silently substitute U for an unknown pole.
    """
    energy = _angle(energy)
    bound = field.kinetic_norm_upper()
    if field.interaction <= 0 or energy <= bound:
        raise ValueError('positive U and energy above kinetic bound required')
    _,_,b,_ = contact_blocks(field)
    bb = b@b.conj().T
    square = float(np.max(np.asarray(abs(bb).sum(axis=1))))
    return square*(abs(1/energy-1/field.interaction)+bound/(energy*(energy-bound)))
