"""Exact counterexamples to omitted hypotheses in the particle papers.

These are mathematical witnesses, not simulated observations or source
selection mechanisms. They keep the reason for each added premise executable.
"""

from fractions import Fraction

import sympy as sp


def test_contraction_on_open_interval_need_not_have_a_fixed_point():
    x = sp.Symbol("x", real=True)
    domain = sp.Interval.open(0, 1)
    fixed = sp.solveset(x / 2 - x, x, domain=sp.S.Reals)
    assert fixed == sp.FiniteSet(0)
    assert fixed.intersect(domain) == sp.EmptySet
    assert sp.diff(x / 2, x) == sp.Rational(1, 2)


def test_zero_cost_distortion_cannot_replace_a_lift_of_the_coarse_optimum():
    coarse = {"A": Fraction(0), "B": Fraction(1)}
    fine = {"b": Fraction(1)}
    contraction = {"b": "B"}
    coarse_optimum = min(coarse, key=coarse.get)
    fine_optimum = min(fine, key=fine.get)
    eta = max(abs(cost - coarse[contraction[name]]) for name, cost in fine.items())
    gap = coarse["B"] - coarse["A"]
    assert gap > 2 * eta
    assert contraction[fine_optimum] != coarse_optimum
    # Supply a lift and the same exact distortion bound: the conclusion holds.
    fine["a"] = Fraction(0)
    contraction["a"] = "A"
    assert contraction[min(fine, key=fine.get)] == coarse_optimum


def test_repeated_stieltjes_nodes_collapse_the_moment_rank():
    def gram(nodes):
        n = len(nodes)
        return sp.Matrix(n, n, lambda i, j: sum(node ** (i + j) for node in nodes) / n)

    repeated = gram([sp.Integer(1)] * 3)
    distinct = gram([sp.Integer(1), sp.Integer(2), sp.Integer(3)])
    assert repeated.rank() == 1
    assert repeated.det() == 0
    assert distinct.det() > 0
    # m=0 makes the resolvent difference vanish, so normalization by D would
    # divide by zero; the construction requires m>0 as well as positive nodes.
    s, m = sp.symbols("s m", positive=True)
    difference = 1 / s - 1 / (s + m**2)
    assert sp.simplify(difference - m**2 / (s * (s + m**2))) == 0
    assert difference.subs(m, 0) == 0


def test_register_count_and_total_weight_do_not_bound_a_thomson_moment():
    s = sp.Symbol("s", positive=True)
    kernel = 1 / (s * (s + 1))
    assert sp.limit(kernel, s, 0, dir="+") == sp.oo
    assert sp.limit(kernel, s, sp.oo) == 0
    assert sp.cancel(sp.diff(kernel, s) + (2 * s + 1) / (s**2 * (s + 1)**2)) == 0
