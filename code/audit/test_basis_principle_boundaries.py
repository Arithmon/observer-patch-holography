"""Exact finite boundaries of the closure interpretation.

These countermodels test only the disputed implications. They are not models
of all OPH premises or realizations of an inhabited universe. In particular,
internal availability of a recharting is weaker than the response-endogeneity
premise used by the compact Lie-type theorem.
"""

from fractions import Fraction
from itertools import permutations
import unittest


def compose(left, right):
    """Composition as functions: apply right first."""
    return tuple(left[right[i]] for i in range(len(left)))


def even(permutation):
    return sum(permutation[i] > permutation[j]
               for i in range(len(permutation))
               for j in range(i + 1, len(permutation))) % 2 == 0


def twelve_point_action():
    """A5 acts on the twelve left cosets of one order-five subgroup."""
    identity = tuple(range(5))
    group = frozenset(p for p in permutations(range(5)) if even(p))
    cycle = (1, 2, 3, 4, 0)
    subgroup, power = set(), identity
    for _ in range(5):
        subgroup.add(power)
        power = compose(cycle, power)
    cosets = {frozenset(compose(g, h) for h in subgroup) for g in group}
    cosets = tuple(sorted(cosets, key=lambda coset: min(coset)))
    index = {coset: i for i, coset in enumerate(cosets)}
    action = {g: tuple(index[frozenset(compose(g, h) for h in coset)]
                       for coset in cosets) for g in group}
    return group, subgroup, cosets, action


class TestBasisPrincipleBoundaries(unittest.TestCase):
    def test_fixed_point_equation_does_not_select_the_information_projection(self):
        uniform = (Fraction(1, 2), Fraction(1, 2))
        biased = (Fraction(3, 4), Fraction(1, 4))
        closure = lambda state: state
        self.assertNotEqual(uniform, biased)
        self.assertEqual(closure(uniform), uniform)
        self.assertEqual(closure(biased), biased)
        # Relative to the uniform reference, D(biased || uniform) is
        # (1/4) log(27/16)>0, while D(uniform || uniform)=0. Thus the same
        # fixed-point equation permits a state excluded by the A3 minimum.
        ratio = Fraction(3, 2) ** 3 * Fraction(1, 2)
        self.assertEqual(ratio, Fraction(27, 16))
        self.assertGreater(ratio, 1)

    def test_internal_twelve_port_action_need_not_be_inner_on_diagonal_current(self):
        group, subgroup, cosets, action = twelve_point_action()
        self.assertEqual(len(group), 60)
        self.assertEqual(len(subgroup), 5)
        self.assertEqual(len(cosets), 12)
        self.assertTrue(all(len(coset) == 5 for coset in cosets))
        self.assertEqual(set().union(*cosets), set(group))
        self.assertEqual(len(set(action.values())), 60)  # faithful
        self.assertEqual({p[0] for p in action.values()}, set(range(12)))
        for g in group:
            for h in group:
                self.assertEqual(action[compose(g, h)], compose(action[g], action[h]))

        # Current generators are i*E_ii. They commute and exponentiate to the
        # connected diagonal torus. Any matrix commuting with every E_ii is
        # diagonal: each off-diagonal matrix unit E_ij has a nonzero commutator
        # with E_ii. Hence both the torus and its unitary centralizer fix every
        # current generator under conjugation.
        # Nevertheless a permitted internal permutation maps E_ii to
        # E_(g.i,g.i). The finite permutation group normalizes the diagonal
        # torus, so the semidirect product is a closed internal operation
        # system with a nontrivial action that is not of the g*c inner form.
        nontrivial = next(p for p in action.values() if p[0] != 0)
        self.assertNotEqual(nontrivial[0], 0)
        self.assertTrue(all(tuple(sorted(p)) == tuple(range(12))
                            for p in action.values()))


if __name__ == "__main__":
    unittest.main()
