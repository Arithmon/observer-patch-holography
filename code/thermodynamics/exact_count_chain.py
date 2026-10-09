"""Finite Markov algebra on supplied rational or represented real weights.

Binary floating-point weights mean their exact represented values. This does
not reconstruct precision lost while acquiring or accumulating those weights.
There is no convergence tolerance, minimum edge size or equilibrium iteration.
"""

from fractions import Fraction as F
from math import gcd
from numbers import Integral, Rational, Real


def rational(value):
    if isinstance(value, bool):
        raise ValueError('Boolean weights are not measured masses')
    if isinstance(value, Integral):
        return F(int(value))
    if isinstance(value, Rational):
        return F(int(value.numerator), int(value.denominator))
    if isinstance(value, Real) and hasattr(value, 'as_integer_ratio'):
        try:
            numerator, denominator = value.as_integer_ratio()
            return F(int(numerator), int(denominator))
        except (ValueError, OverflowError) as exc:
            raise ValueError('weights must be finite') from exc
    raise ValueError('weights must be unmasked real rational or binary scalars')


def weight_matrix(weights):
    try:
        rows = [[rational(value) for value in row] for row in weights]
    except TypeError as exc:
        raise ValueError('weights must be a nonempty square matrix') from exc
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError('weights must be a nonempty square matrix')
    if any(value < 0 for row in rows for value in row):
        raise ValueError('weights must be nonnegative')
    return rows


def kernel_from_weights(weights, empty_rows='reject'):
    """Normalize supplied masses; absorbing completion must be explicit."""
    if empty_rows not in ('reject', 'absorbing'):
        raise ValueError('unknown empty-row policy')
    rows = weight_matrix(weights)
    result = []
    for i, row in enumerate(rows):
        mass = sum(row)
        if mass:
            result.append([value/mass for value in row])
        elif empty_rows == 'absorbing':
            result.append([F(i == j) for j in range(len(rows))])
        else:
            raise ValueError('uncounted row has no transition law')
    return result


def stochastic_matrix(kernel):
    rows = weight_matrix(kernel)
    if any(sum(row) != 1 for row in rows):
        raise ValueError('kernel rows must sum to one exactly; supply weights for normalization')
    return rows


def support_reachability(kernel):
    rows = stochastic_matrix(kernel)
    n = len(rows)
    reach = [[i == j or rows[i][j] > 0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] |= reach[k][j]
    return reach


def closed_classes(kernel):
    rows = stochastic_matrix(kernel)
    reach = support_reachability(rows)
    unseen, result = set(range(len(rows))), []
    while unseen:
        seed = min(unseen)
        component = sorted(j for j in unseen if reach[seed][j] and reach[j][seed])
        unseen.difference_update(component)
        if all(rows[i][j] == 0 for i in component
               for j in range(len(rows)) if j not in component):
            result.append(component)
    return result


def irreducible_period(kernel):
    """GCD of directed edge level differences in a strongly connected graph."""
    rows = stochastic_matrix(kernel)
    if not all(all(row) for row in support_reachability(rows)):
        raise ValueError('period requires an irreducible kernel')
    distance, queue = {0: 0}, [0]
    for source in queue:
        for target, weight in enumerate(rows[source]):
            if weight and target not in distance:
                distance[target] = distance[source]+1
                queue.append(target)
    period = 0
    for source, row in enumerate(rows):
        for target, weight in enumerate(row):
            if weight:
                period = gcd(period, distance[source]+1-distance[target])
    return period


def stationary_distribution(kernel):
    """Unique stationary probability, including transient states when present.

    A finite chain has a unique stationary law iff it has one closed class.
    Solve the independent balance equations on that class, replace one with
    normalization, and verify *all* original balance equations exactly.
    Multiple closed classes require a separate choice of their masses.
    """
    rows = stochastic_matrix(kernel)
    classes = closed_classes(rows)
    if len(classes) != 1:
        raise ValueError('stationary law is not unique; closed-class masses are unspecified')
    indices = classes[0]
    n = len(indices)
    equations = [[rows[indices[j]][indices[i]]-F(i == j) for j in range(n)]+[F(0)]
                 for i in range(n-1)]
    equations.append([F(1)]*(n+1))
    for col in range(n):
        pivot = next((i for i in range(col, n) if equations[i][col]), None)
        if pivot is None:
            raise ValueError('stationary balance system is singular')
        equations[col], equations[pivot] = equations[pivot], equations[col]
        divisor = equations[col][col]
        equations[col] = [value/divisor for value in equations[col]]
        for i in range(n):
            if i != col and equations[i][col]:
                multiplier = equations[i][col]
                equations[i] = [a-multiplier*b for a, b in zip(equations[i], equations[col])]
    pi = [F(0)]*len(rows)
    for i, index in enumerate(indices):
        pi[index] = equations[i][-1]
    if (sum(pi) != 1 or any(value < 0 for value in pi)
            or any(sum(pi[i]*rows[i][j] for i in range(len(rows))) != pi[j]
                   for j in range(len(rows)))):
        raise ValueError('stationary solution does not satisfy original equations')
    return pi


def detailed_balance_defect(kernel, pi):
    rows = stochastic_matrix(kernel)
    law = [rational(value) for value in pi]
    if len(law) != len(rows) or min(law) < 0 or sum(law) != 1:
        raise ValueError('reference must be a probability on the complete state set')
    return max(abs(law[i]*rows[i][j]-law[j]*rows[j][i])
               for i in range(len(rows)) for j in range(len(rows)))


def lumpability_defect(kernel, blocks):
    rows = stochastic_matrix(kernel)
    try:
        blocks = [list(block) for block in blocks]
        indices = [index for block in blocks for index in block]
    except TypeError as exc:
        raise ValueError('blocks must partition the state set') from exc
    if (not blocks or any(not block for block in blocks)
            or any(isinstance(i, bool) or not isinstance(i, Integral) for i in indices)
            or sorted(indices) != list(range(len(rows)))):
        raise ValueError('blocks must partition the state set exactly once')
    defect = F(0)
    for source in blocks:
        for target in blocks:
            masses = [sum(rows[i][j] for j in target) for i in source]
            defect = max(defect, max(masses)-min(masses))
    return defect
