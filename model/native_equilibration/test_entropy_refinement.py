"""Independent dimension, tangent, native-channel and conditioning controls."""
from collections import Counter
from fractions import Fraction as F
import json
from math import comb

import numpy as np
import pytest
import sympy as sp

from . import entropy_refinement as e
from .construction import exact_native_channel, graph_certificate
from ..observer import native


@pytest.mark.parametrize("dimension", [1, 2, 6, 66, comb(3840, 1280)])
@pytest.mark.parametrize("level,tangent_bits", [(0, 0), (2, 0), (7, 9)])
def test_integer_schedule_pays_dimension_and_all_declared_tangent_norms(dimension, level, tangent_bits):
    gamma = F(1, 123456789)
    row = e.entropy_budget(dimension, gamma, level, tangent_bits=tangent_bits)
    if dimension == 1:
        assert row["discarded_native_attempts"] == 0
        assert F(row["joint_counting_gradient_upper"]) == 0
        return
    h, blocks = row["log2_dimension_upper"], row["blocks"]
    assert (1 << (h-1)) < dimension <= 1 << h
    assert row["half_contraction_block"]*gamma >= 1
    assert row["discarded_native_attempts"] == row["half_contraction_block"]*blocks
    # The direct d/2 * 2^-blocks bound, no floating logarithms or n-th powers.
    epsilon = F(dimension, 1 << (blocks+1))
    theta = dimension*epsilon
    assert theta <= F(row["relative_density_error_upper"])
    speed, ell = 1 << tangent_bits, F(1, 1 << level)
    assert speed*(theta/(1-theta)+theta) <= F(row["joint_counting_gradient_upper"])
    assert F(row["joint_counting_gradient_upper"]) <= F(3, 32)*ell**5
    assert F(row["quantum_state_gradient_upper"]) <= F(1, 16)*ell**5
    assert 2*dimension*epsilon**2 <= F(row["entropy_deficit_upper"]) <= ell**10/512


def test_all_level_schedule_uses_only_connected_source_census(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("the family bound must not allocate a graph or sector")
    monkeypatch.setattr(e, "comb", forbidden)
    monkeypatch.setattr(e.mixing.source, "federation", forbidden)
    for level in (0, 3, 5, 20):
        row = e.screen_family_budget(level)
        v, m = 240*4**level, 630*4**level
        assert (row["vertices"], row["seams"]) == (v, m)
        assert F(row["gap_lower_per_attempt"]) == F(2, m*v*(v-1))
        assert row["log2_dimension_upper"] == v
        assert row["blocks"] == 2*v+4+5*level
        assert row["discarded_native_attempts"] == (m*v*(v-1)//2)*row["blocks"]
        assert F(row["gradient_over_ell4_upper"]) == F(3, 32*2**level)


def test_actual_retained_screen_levels_supply_the_connectivity_gap_and_density():
    report = e.run()
    for row in report["retained_actual_screen_certificates"]:
        graph = e.mixing.source.federation(row["level"])
        independent = graph_certificate(graph.ports, graph.seams)
        assert row["particles"] == graph.ports//3
        assert row["diameter_upper_bound"] == independent["diameter_upper_bound"]
        assert row["quantum_gap_lower_per_attempt"] == independent["quantum_gap_lower_per_attempt"]
        assert F(row["gap_lower_per_attempt"]) >= F(e.screen_family_budget(row["level"])["gap_lower_per_attempt"])
        d = comb(graph.ports, row["particles"])
        assert row["log2_dimension_upper"] == (d-1).bit_length()
        assert row["discarded_native_attempts"] < e.screen_family_budget(row["level"])["discarded_native_attempts"]


def test_exact_entropy_derivative_witness_and_dimension_free_error_failure():
    t = sp.symbols("t", real=True)
    d = sp.symbols("d", positive=True)
    theta = sp.Rational(1, 2)
    # All unlisted d-2 singleton probabilities are 1/d and stay fixed.
    p = ((1+theta)/d-t/2, (1-theta)/d+t/2)
    entropy = -sum(x*sp.log(x) for x in p)
    derivative = sp.expand_log(sp.diff(entropy, t).subs(t, 0), force=True).simplify()
    assert sp.simplify(derivative-sp.log(3)/2) == 0
    # TV=1/(2d) ->0 and entropy deficit ->0, but this unit-L1 gradient
    # remains log(3)/2. Preparation error must be charged against dimension.
    deficit = sum(x.subs(t, 0)*sp.log(d*x.subs(t, 0)) for x in p)
    assert sp.limit(deficit, d, sp.oo) == 0
    # Joint count variation attains both range terms, with fixed total d.
    capacities = (1+d*t/2, 1-d*t/2)
    joint = -sum(x*sp.log(x/y) for x, y in zip(p, capacities, strict=True))
    assert sp.simplify(sp.diff(joint, t).subs(t, 0)-sp.log(3)/2-theta) == 0


@pytest.mark.parametrize("law", ["C", "E"])
def test_actual_small_graph_quantum_preparation_and_noncommuting_tangents(law):
    # Independent source Kraus enumeration checks the inherited theorem on
    # a different graph, not just the named elementary 66-state carrier.
    edges = ((0, 1), (1, 2))
    counts, scale, basis = exact_native_channel(3, edges, 1, law)
    d = len(basis)
    gamma = F(graph_certificate(3, edges)["quantum_gap_lower_per_attempt"])
    row = e.entropy_budget(d, gamma, 0)
    ket = np.array([1, 1j, 0])/np.sqrt(2)
    rho = np.outer(ket, ket.conj())
    channel = counts/scale
    for _ in range(row["discarded_native_attempts"]):
        rho = (channel@rho.ravel()).reshape(d, d)
    values, vectors = np.linalg.eigh(rho)
    assert values.min() > 0
    epsilon = np.abs(np.linalg.eigvalsh(rho-np.eye(d)/d)).sum()/2
    theta = float(F(row["relative_density_error_upper"]))
    assert d*epsilon <= theta
    # These off-diagonal and diagonal directions span the traceless space.
    log_rho = (vectors*np.log(values))@vectors.conj().T
    directions = []
    for i in range(d):
        for j in range(i):
            for phase in (1, 1j):
                x = np.zeros((d, d), complex)
                x[i, j], x[j, i] = phase/2, np.conj(phase)/2
                directions.append(x)
    for j in range(1, d):
        x = np.zeros((d, d))
        x[0, 0], x[j, j] = F(1, 2), -F(1, 2)
        directions.append(x)
    for x in directions:
        assert abs(np.trace(x@log_rho)) <= float(F(row["quantum_state_gradient_upper"]))
    deficit = np.log(d)+np.sum(values*np.log(values))
    assert -1e-12 <= deficit <= float(F(row["entropy_deficit_upper"]))


def test_retained_current_and_conserved_sector_are_not_erased_by_preparation():
    states = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    unread, plus = Counter(), Counter()
    for state in states:
        for coin in (0, 1):
            output, (_, current) = native.native(state, (0, 1), coin)
            unread[output] += 1
            if current == 1:
                plus[output] += 1
    assert set(unread.values()) == {2} and len(unread) == 3
    assert sum(plus.values()) == 1 and len(plus) == 1  # success 1/6, pure posterior
    assert sp.log(len(unread))-sp.log(len(plus)) == sp.log(3)
    for state in ((0, 0, 0), (1, 1, 1)):
        for seam in ((0, 1), (1, 2)):
            for coin in (0, 1):
                assert native.native(state, seam, coin)[0] == state


@pytest.mark.parametrize("args,kwargs", [
    ((True, F(1, 2), 0), {}), ((0, F(1, 2), 0), {}),
    ((2, 0.5, 0), {}), ((2, F(0), 0), {}), ((2, F(2), 0), {}),
    ((2, F(1, 2), True), {}), ((2, F(1, 2), -1), {}),
    ((2, F(1, 2), 0), {"tangent_bits": True}),
    ((2, F(1, 2), 0), {"tangent_bits": -1}),
])
def test_undeclared_dimension_gap_refinement_or_tangent_scale_rejected(args, kwargs):
    with pytest.raises(ValueError):
        e.entropy_budget(*args, **kwargs)


def test_retained_refinement_receipt_replays():
    assert json.loads(e.OUTPUT.read_text(encoding="utf-8")) == e.run()
