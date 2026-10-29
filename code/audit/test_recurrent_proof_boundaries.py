"""Exact witnesses for proof domains and noncommutative transport conventions.

The metric-family tests check finite instances and exact witness formulas;
the infinite completeness and nonattainment arguments are in the normal-form
manuscript. Finite truncations are not reported as proofs of an infinite limit.
"""

from fractions import Fraction as Q
import itertools


def test_approximate_descent_allows_nonstationary_periodic_limit_points():
    # On r=1, Phi(r,theta)=sin(theta)/2+(r-1)*sqrt(2-x/2-x*x/4),
    # x=cos(theta), and the smooth flow is theta'=1, r'=0.
    # Its tangential and radial gradient components are orthogonal.
    for x in (Q(-1), Q(-3, 5), Q(0), Q(3, 5), Q(1)):
        tangent_squared = x*x/4
        radial_squared = 2-x/2-x*x/4
        gradient_squared = tangent_squared + radial_squared
        derivative = x/2
        assert radial_squared > 0
        assert derivative == -gradient_squared + 2
    # At theta=pi each revolution, the gradient exceeds epsilon/c=2.
    assert Q(5, 2) > 2


def test_energy_nonincrease_does_not_supply_scheduler_fairness():
    # Independent spins with positive fields: (-1) at coordinate 1 is
    # unstable forever if the scheduler only revisits coordinate 0.
    state = (1, -1)
    def update(s, i):
        return tuple(1 if j == i else value for j, value in enumerate(s))
    assert update(state, 0) == state
    assert update(state, 1) != state


def test_tie_preserving_hopfield_flips_strictly_lower_energy():
    weights = ((0, 2, -1), (2, 0, 1), (-1, 1, 0))
    thresholds = (1, 0, -1)
    def energy(s):
        return -sum(weights[i][j]*s[i]*s[j] for i in range(3)
                    for j in range(i)) + sum(t*x for t, x in zip(thresholds, s))
    saw_tie = False
    for s in itertools.product((-1, 1), repeat=3):
        for i in range(3):
            field = sum(w*x for w, x in zip(weights[i], s)) - thresholds[i]
            spin = s[i] if field == 0 else (1 if field > 0 else -1)
            after = s[:i] + (spin,) + s[i+1:]
            if field == 0:
                saw_tie = True
                assert after == s
            if after != s:
                assert energy(after) - energy(s) == -2*abs(field) < 0
    assert saw_tie


def test_equilibrium_propagation_sign_on_exact_quadratic_branch():
    # E=(s-theta)^2/2, C=(s-y)^2/2, s_beta=(theta+beta*y)/(1+beta).
    theta, y = Q(3), Q(1)
    for beta in (Q(1, 2), Q(1, 10), Q(1, 100)):
        s_beta = (theta+beta*y)/(1+beta)
        contrast = (theta-s_beta)/beta
        assert contrast == (theta-y)/(1+beta)
        assert 0 < contrast < theta-y
        # Finite-nudge contrast is biased; the limit has the positive sign.
        assert (theta-y)-contrast == (theta-y)*beta/(1+beta)


_IDENTITY = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
_SWAP01 = ((0, 1, 0), (1, 0, 0), (0, 0, 1))
_SWAP12 = ((1, 0, 0), (0, 0, 1), (0, 1, 0))


def _matrix_product(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def _inverse_permutation_matrix(a):
    return tuple(zip(*a))


def _chronological_transport(edges, path):
    result = _IDENTITY
    for source, target in zip(path, path[1:]):
        result = _matrix_product(result, edges[source, target])
    return result


def test_strict_noncommuting_triangle_requires_chronological_comparisons():
    frames = (_IDENTITY, _SWAP01, _SWAP12)
    edges = {(i, j): _matrix_product(frames[i], _inverse_permutation_matrix(frames[j]))
             for i, j in itertools.product(range(3), repeat=2)}
    for i, j, k in itertools.product(range(3), repeat=3):
        assert _matrix_product(edges[i, j], edges[j, k]) == edges[i, k]
        assert _chronological_transport(edges, (i, j, k, i)) == _IDENTITY
    # Every triangle defect is exactly trivial, but reversing only the
    # matrix product gives spurious holonomy on the filled triangle.
    reverse_order = _matrix_product(
        _matrix_product(edges[2, 0], edges[1, 2]), edges[0, 1])
    assert reverse_order == ((0, 0, 1), (1, 0, 0), (0, 1, 0))
    assert reverse_order != _IDENTITY
    comparison = _chronological_transport(edges, (0, 1, 2))
    assert comparison == edges[0, 2]
    # A comparison maps terminal-chart column components to initial ones;
    # its adjoint is forward column transport.
    vector = (1, 2, 3)
    forward = _inverse_permutation_matrix(comparison)
    assert tuple(sum(forward[i][j]*vector[j] for j in range(3))
                 for i in range(3)) == (1, 3, 2)


def _nonflat_edges():
    edges = {(i, i): _IDENTITY for i in range(4)}
    for i, j in itertools.combinations(range(4), 2):
        value = _SWAP01 if (i, j) == (1, 2) else (
            _SWAP12 if (i, j) == (2, 3) else _IDENTITY)
        edges[i, j] = value
        edges[j, i] = _inverse_permutation_matrix(value)
    return edges


def test_path_ratios_reverse_paths_and_frame_changes_use_one_order():
    edges = _nonflat_edges()
    p, q = (0, 1, 2, 3), (0, 2, 3)
    up = _chronological_transport(edges, p)
    uq = _chronological_transport(edges, q)
    ratio = _matrix_product(up, _inverse_permutation_matrix(uq))
    loop = p+tuple(reversed(q))[1:]
    assert _chronological_transport(edges, loop) == ratio != _IDENTITY
    assert _chronological_transport(edges, tuple(reversed(p))) == (
        _inverse_permutation_matrix(up))
    # The loop at the terminal chart uses the opposite ratio order. The
    # two based holonomies are conjugate, not generally equal.
    terminal_ratio = _matrix_product(_inverse_permutation_matrix(uq), up)
    terminal_loop = tuple(reversed(q))+p[1:]
    assert _chronological_transport(edges, terminal_loop) == terminal_ratio
    assert terminal_ratio != ratio
    assert terminal_ratio == _matrix_product(
        _matrix_product(_inverse_permutation_matrix(uq), ratio), uq)

    frames = (_SWAP12, _IDENTITY, _SWAP01, _matrix_product(_SWAP01, _SWAP12))
    changed = {(i, j): _matrix_product(
        _matrix_product(frames[i], value), _inverse_permutation_matrix(frames[j]))
        for (i, j), value in edges.items()}
    changed_up = _chronological_transport(changed, p)
    changed_uq = _chronological_transport(changed, q)
    assert changed_up == _matrix_product(
        _matrix_product(frames[0], up), _inverse_permutation_matrix(frames[3]))
    assert _matrix_product(changed_up, _inverse_permutation_matrix(changed_uq)) == (
        _matrix_product(_matrix_product(frames[0], ratio),
                        _inverse_permutation_matrix(frames[0])))
    assert _matrix_product(_inverse_permutation_matrix(changed_uq), changed_up) == (
        _matrix_product(_matrix_product(frames[3], terminal_ratio),
                        _inverse_permutation_matrix(frames[3])))


def test_crossed_module_quadruple_coherence_follows_associativity():
    # H=G=S3, boundary=id and action=conjugation form a crossed module.
    # Arbitrary edge comparisons determine h_ijk=g_ij*g_jk*g_ik^-1.
    edges = _nonflat_edges()

    def h(i, j, k):
        return _matrix_product(_matrix_product(edges[i, j], edges[j, k]), edges[k, i])

    def action(g, value):
        return _matrix_product(_matrix_product(g, value), _inverse_permutation_matrix(g))

    for i, j, k, ell in itertools.product(range(4), repeat=4):
        assert _matrix_product(h(i, j, k), h(i, k, ell)) == (
            _matrix_product(action(edges[i, j], h(j, k, ell)), h(i, j, ell)))
    # Permuting the triple indices changes AB into BA on four distinct
    # vertices and fails even though the triangle law holds exactly.
    wrong_left = _matrix_product(h(1, 2, 3), h(0, 1, 3))
    wrong_right = _matrix_product(action(edges[0, 1], h(0, 2, 3)), h(0, 1, 2))
    assert wrong_left != wrong_right


def _star_state_distance(p, q):
    if p == q:
        return Q(0)

    def length(point):
        return Q(0) if point == 'x' else Q(4) if point == 'y' else 1+Q(1, point)

    return length(p)+length(q)


def _star_observation_distance(p, q):
    if p == q or (p in ('x', 'y') and q in ('x', 'y')):
        return Q(0)
    return sum(Q(0) if point in ('x', 'y') else 1+Q(1, point) for point in (p, q))


def test_complete_star_family_needs_attained_consistency_witness():
    # The formulas describe infinite stars; the manuscript proves every
    # distinct distance is >=1, hence completeness, and inf_n(1+1/n)=1.
    points = ('x', 'y', *range(1, 9))
    for p, q, r in itertools.product(points, repeat=3):
        assert _star_state_distance(p, r) <= (
            _star_state_distance(p, q)+_star_state_distance(q, r))
        assert _star_observation_distance(p, r) <= (
            _star_observation_distance(p, q)+_star_observation_distance(q, r))
    for p, q in itertools.product(points, repeat=2):
        assert _star_observation_distance(p, q) <= _star_state_distance(p, q)

    # The limit radius one is strictly below every candidate distance,
    # including exact rational candidates far below binary64 resolution.
    for n in (1, 2, 8, 10**20, 10**400):
        assert _star_state_distance('x', n) > 1
        assert _star_observation_distance('y', n) > 1
    for epsilon in (Q(1), Q(1, 7), Q(1, 10**20), Q(1, 10**400)):
        n = epsilon.denominator//epsilon.numerator+1
        assert 1 < _star_state_distance('x', n) < 1+epsilon

    # Phi(x)=1, Bx=By=0, and N(x)=y give a (delta,beta)=(1,0)
    # receipt. Every distinct consistent observation is >1 apart, so
    # omega(1)=0. A bare distance-to-C bound would falsely give 4<=1.
    claimed_radius, inverse_at_radius = Q(1), Q(0)
    assert _star_state_distance('x', 'y') > claimed_radius+inverse_at_radius
    # A singleton fine level restricting to x has the same false one-step
    # naturality bound unless the coarse-level witness premise is imposed.
    assert _star_state_distance('x', 'y') == 4


def test_finite_star_cutoffs_supply_a_witness_and_a_valid_receipt_bound():
    for size in (1, 2, 7, 16):
        consistent = ('y', *range(1, size+1))
        eta = min(_star_state_distance('x', point) for point in consistent)
        witness = min(consistent, key=lambda point: _star_state_distance('x', point))
        assert eta == 1+Q(1, size) > 1
        assert _star_state_distance('x', witness) == eta
        omega = max(_star_state_distance(p, q)
                    for p, q in itertools.product(consistent, repeat=2)
                    if _star_observation_distance(p, q) <= eta)
        assert omega == 4+eta
        assert _star_state_distance('x', 'y') <= eta+omega
