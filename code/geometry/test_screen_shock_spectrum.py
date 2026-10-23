#!/usr/bin/env python3
"""Tests for the finite-screen shock-sign and graph-spectrum receipt."""

from __future__ import annotations

import copy
from collections import Counter
from fractions import Fraction
import json
import pathlib
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT = HERE / "derive_screen_shock_spectrum.py"
ARTIFACT = HERE / "runs" / "screen_shock_spectrum.json"
sys.path.insert(0, str(HERE))

from derive_screen_shock_spectrum import (  # noqa: E402
    ReceiptValidationError,
    _icosahedron_graphs,
    _laplacian,
    build_receipt,
    counting_stationarity_bound,
    finite_native_stationarity_receipt,
    receipt_failures,
    validate_receipt,
)


def _run_to(path: pathlib.Path) -> dict:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(path)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(path.read_text(encoding="utf-8"))


def test_screen_shock_receipt_is_byte_exact_and_lf_only(
    tmp_path: pathlib.Path,
) -> None:
    rebuilt = tmp_path / ARTIFACT.name
    _run_to(rebuilt)
    assert rebuilt.read_bytes() == ARTIFACT.read_bytes()
    assert b"\r\n" not in rebuilt.read_bytes()
    subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )


def test_pure_de_sitter_mu_squared_equals_lambda_l1_for_d3_through_d6() -> None:
    rows = build_receipt()["pure_de_sitter_identity"]["checks"]
    assert [row["spacetime_dimension"] for row in rows] == [3, 4, 5, 6]
    assert [row["mu_squared_exact"] for row in rows] == ["1", "2", "3", "4"]
    assert all(row["mu_squared_exact"] == row["lambda_l1_exact"] for row in rows)
    assert all(row["kappa_times_r_c_exact"] == "1" for row in rows)


def test_entropy_maximum_and_area_hessian_have_distinct_statuses() -> None:
    capacity = build_receipt()["entropy_and_capacity"]
    entropy_maximum = capacity["finite_generalized_entropy_maximum"]
    area = capacity["area_analytic_relaxation"]

    assert entropy_maximum["identity"] == "S_gen(p,d)=log(M)-D_KL(p || d/M)"
    assert entropy_maximum["unique_simplex_maximum"] is True
    assert entropy_maximum["maximum_value"] == "log(M)"

    assert area["domain"] == "positive-real relaxation of integer sector dimensions d_i"
    assert area["logarithmic_area_coordinate"] == (
        "F(d)=sum_i d_i*log(d_i)/M"
    )
    assert area["gradient_component_exact"] == "1/(n*d)"
    assert area["gradient_component"] > 0.0
    assert area["unconstrained_stationary"] is False
    assert area["hessian_exact"] == "(1/(n*d^2))*(I-(2/n)*J)"
    assert area["fixed_M_tangent_eigenvalue"] > 0.0
    assert area["fixed_M_tangent_classification"] == "strict_local_minimum"
    assert area["homogeneous_eigenvalue"] < 0.0
    assert area["homogeneous_directional_second_derivative"] < 0.0


@pytest.mark.parametrize("dimensions", [(2, 3, 5), (1, 4, 5), (4, 4, 2)])
def test_fixed_integer_capacity_entropy_has_the_exact_relative_entropy_deficit(dimensions):
    # Distinct actual integer partitions of one total capacity, no derivative
    # or dimension-changing repair step is invoked.
    total = sum(dimensions)
    probability = tuple(map(sp.Rational, ("1/2", "1/3", "1/6")))
    reference = tuple(sp.Rational(d, total) for d in dimensions)
    entropy = -sum(p*sp.log(p/d) for p, d in zip(probability, dimensions, strict=True))
    deficit = sum(p*sp.log(p/q) for p, q in zip(probability, reference, strict=True))
    assert sp.expand_log(entropy-sp.log(total)+deficit, force=True).simplify() == 0
    maximum = -sum(q*sp.log(q/d) for q, d in zip(reference, dimensions, strict=True))
    assert sp.simplify(maximum-sp.log(total)) == 0


@pytest.mark.parametrize("probability_velocity", [
    ("1/20", "-1/10", "1/20"), ("1/10", "-1/5", "1/10"), ("0", "0", "0"),
])
def test_joint_relaxed_state_capacity_derivatives_equal_the_constrained_score_square(probability_velocity):
    t = sp.symbols("t", real=True)
    base = tuple(map(sp.Rational, (2, 3, 5)))
    v = tuple(map(sp.Rational, (1, -2, 1)))
    u = tuple(map(sp.Rational, probability_velocity))
    # Nonzero acceleration also checks cancellation for curved constrained paths.
    acceleration = tuple(map(sp.Rational, (2, -1, -1)))
    dimensions = tuple(d+t*x+t*t*a for d, x, a in zip(base, v, acceleration, strict=True))
    probability = tuple(d/10+t*x-t*t*a/100 for d, x, a in zip(base, u, acceleration, strict=True))
    assert sum(dimensions) == 10 and sum(probability) == 1
    entropy = -sum(p*sp.log(p/d) for p, d in zip(probability, dimensions, strict=True))
    assert sp.simplify(sp.diff(entropy, t).subs(t, 0)) == 0
    second = sp.simplify(sp.diff(entropy, t, 2).subs(t, 0))
    score_square = -sum((up-vp/10)**2/(d/10) for up, vp, d in zip(u, v, base, strict=True))
    assert second == score_square <= 0
    assert (second == 0) == all(up == vp/10 for up, vp in zip(u, v, strict=True))


def test_native_uniform_microstates_supply_the_actual_sector_counts():
    from model.physical_space_time.occupation_field import native_counts

    states, remote, incident, _, _ = native_counts()
    kernel = (remote+incident)/60
    uniform = sp.ones(66, 1)/66
    assert kernel*uniform == uniform
    # The observed sector is the load in four named ports. Its multiplicities
    # come from the native states, not an assigned thermodynamic weight.
    counts = [sum(sum(state[:4]) == load for state in states) for load in range(3)]
    assert counts == [28, 32, 6]
    masses = [sum(uniform[i] for i, state in enumerate(states) if sum(state[:4]) == load)
              for load in range(3)]
    assert masses == [sp.Rational(count, 66) for count in counts]


def test_bulk_dimensions_must_enter_the_sector_count_and_the_entropy_deficit():
    edge = (2, 3)
    bulk = ((sp.Rational(1, 3), sp.Rational(2, 3)), (sp.Rational(1, 1),))
    probability = (sp.Rational(2, 5), sp.Rational(3, 5))
    weights = tuple(d*len(rho) for d, rho in zip(edge, bulk, strict=True))
    total = sum(weights)
    eigenvalues = [p*r/d for p, d, rho in zip(probability, edge, bulk, strict=True)
                   for r in rho for _ in range(d)]
    entropy = -sum(x*sp.log(x) for x in eigenvalues)
    sector_deficit = sum(p*sp.log(p*total/w) for p, w in zip(probability, weights, strict=True))
    bulk_deficit = sum(p*sum(r*sp.log(r*len(rho)) for r in rho)
                       for p, rho in zip(probability, bulk, strict=True))
    assert sp.expand_log(entropy-sp.log(total)+sector_deficit+bulk_deficit, force=True).simplify() == 0
    assert total == 7 and sum(edge) == 5
    assert sp.simplify(bulk_deficit) != 0


def test_relative_entropy_stationarity_does_not_select_the_counting_entropy_maximum():
    t = sp.symbols("t", real=True)
    reference = sp.Matrix([[sp.Rational(1, 4), sp.Rational(3, 4)]])
    reset = sp.ones(2, 1)*reference
    assert reference*reset == reference and reset*reset == reset
    probability = (sp.Rational(1, 4)+t, sp.Rational(3, 4)-t)
    relative = sum(p*sp.log(p/q) for p, q in zip(probability, reference, strict=True))
    entropy = -sum(p*sp.log(p) for p in probability)  # Two sectors of count one.
    assert sp.simplify(sp.diff(relative, t).subs(t, 0)) == 0
    assert sp.simplify(sp.diff(entropy, t).subs(t, 0)) == sp.log(3)


def test_finite_joint_bound_has_a_sharp_complete_tangent_witness():
    # Differentiate the entropy itself, independently of the bound's formula.
    t = sp.symbols("t", real=True)
    epsilon = sp.Rational(1, 6600)
    p = (sp.Rational(1, 66)+epsilon-t/2, sp.Rational(1, 66)-epsilon+t/2)
    d = (1+33*t, 1-33*t)
    entropy = -sum(x*sp.log(x/y) for x, y in zip(p, d, strict=True))
    derivative = sp.expand_log(sp.diff(entropy, t).subs(t, 0), force=True).simplify()
    expected = sp.log(sp.Rational(101, 99))/2+sp.Rational(1, 100)
    assert sp.simplify(derivative-expected) == 0
    bound = counting_stationarity_bound(Fraction(1, 6600), Fraction(1, 66))
    assert float(derivative) <= float(bound["joint_derivative_upper"])
    assert bound["relative_population_error"] == Fraction(1, 100)
    # A coherent two-dimensional block attains the quantum-state bound too.
    rho = sp.Matrix([[sp.Rational(1, 66), epsilon], [epsilon, sp.Rational(1, 66)]])
    tangent = sp.Matrix([[0, -sp.Rational(1, 2)], [-sp.Rational(1, 2), 0]])
    entropy_gradient = -sp.trace(tangent*rho.log()).simplify()
    assert sp.simplify(entropy_gradient-sp.log(sp.Rational(101, 99))/2) == 0


@pytest.mark.parametrize("arguments", [
    (True, Fraction(1, 66)), (0.001, Fraction(1, 66)),
    (Fraction(1, 66), Fraction(1, 66)), (0, 0), (-1, 1), (0, 2),
    (0, Fraction(1, 66), -1), (0, Fraction(1, 66), 1, -1),
])
def test_stationarity_bound_rejects_missing_floor_and_ambiguous_units(arguments):
    with pytest.raises(ValueError):
        counting_stationarity_bound(*arguments)


def test_native_finite_preparation_bounds_all_partition_tangents():
    from model.observer import native, quantum_mixing
    from model.record_laws import burn_in
    from mpmath import MPContext

    states = tuple(state for state in native.states() if max(state) == 1)
    index = {state: i for i, state in enumerate(states)}
    # Independently compile integer transition columns from native coin calls.
    columns = [Counter(index[native.native(state, seam, coin)[0]]
                       for seam in quantum_mixing.seams() for coin in (0, 1))
               for state in states]
    counts, denominator = [1]+[0]*65, 1
    report = finite_native_stationarity_receipt()
    row = report["finite_budgets"][-1]
    for _ in range(row["repair_attempts"]):
        following = [0]*66
        for amount, column in zip(counts, columns, strict=True):
            for output, multiplicity in column.items():
                following[output] += amount*multiplicity
        counts, denominator = following, denominator*60
    p = [Fraction(value, denominator) for value in counts]
    actual_tv = sum(abs(x-Fraction(1, 66)) for x in p)/2
    assert actual_tv <= burn_in.binary_bound(row["repair_attempts"])
    assert Fraction(row["joint_derivative_upper"]) <= Fraction(row["joint_derivative_target"])
    ctx = MPContext()
    ctx.dps = 80
    for partition in ([i for i in range(66)], [sum(s[:4]) for s in states],
                      [int(s[0]) for s in states]):
        labels = sorted(set(partition))
        q = [Fraction(partition.count(label), 66) for label in labels]
        coarse = [sum((p[i] for i, label in enumerate(partition) if label == wanted), Fraction())
                  for wanted in labels]
        ratios = [ctx.mpf(x.numerator)/x.denominator/(ctx.mpf(y.numerator)/y.denominator)
                  for x, y in zip(coarse, q, strict=True)]
        logs = list(map(ctx.log, ratios))
        # Exact dual norms over all zero-sum, unit-L1 state/count tangents.
        derivative_sup = (max(logs)-min(logs)+max(ratios)-min(ratios))/2
        bound = counting_stationarity_bound(actual_tv, min(q))
        rational = bound["joint_derivative_upper"]
        assert derivative_sup <= ctx.mpf(rational.numerator)/rational.denominator
        deficit = sum((ctx.mpf(x.numerator)/x.denominator)*log
                      for x, log in zip(coarse, logs, strict=True))
        rational = bound["entropy_deficit_upper"]
        assert 0 <= deficit <= ctx.mpf(rational.numerator)/rational.denominator


def test_retained_record_and_full_carrier_expose_stationarity_domain():
    from model.observer import native, quantum_mixing

    states = native.states()
    binary = tuple(state for state in states if max(state) == 1)
    plus = Counter()
    unread = Counter()
    for state in binary:
        for coin in (0, 1):
            output, (_, current) = native.native(state, (0, 1), coin)
            unread[output] += 1
            if current == 1:
                plus[output] += 1
    assert set(unread.values()) == {2} and len(unread) == 66
    assert len(plus) == 10 and set(plus.values()) == {1}
    assert Fraction(sum(plus.values()), 132) == Fraction(5, 66)
    assert sp.simplify(sp.log(66)-sp.log(len(plus))-sp.log(sp.Rational(33, 5))) == 0
    # The full 78-state channel is not unital: from its entropy maximum,
    # one repair strictly decreases entropy. The binary theorem cannot cover it.
    full = Counter(native.native(state, seam, coin)[0]
                   for state in states for seam in quantum_mixing.seams() for coin in (0, 1))
    assert len(full) == 78 and sum(full.values()) == 78*60
    assert set(full.values()) != {60}


def test_native_stationarity_receipt_rejects_a_promoted_or_mutated_bound():
    for section, key, value in [
        ("scope", "physical_fixed_volume_or_area_coefficient_derived", True),
        ("scope", "full_78_state_entropy_gradient_covered", True),
        (None, "record_access", "retain the tape"),
    ]:
        receipt = build_receipt()
        row = receipt["entropy_and_capacity"]["finite_native_stationarity"]
        (row[section] if section else row)[key] = value
        assert "finite_native_stationarity_mismatch" in receipt_failures(receipt)


def _sector_partial_traces(rho, states, region):
    """Independent partial traces in the actual occupied basis, no 2^N matrix."""
    complement = tuple(i for i in range(len(states[0])) if i not in region)
    labels = [tuple(tuple(state[i] for i in part) for state in states)
              for part in (region, complement)]
    indices = [{label: i for i, label in enumerate(sorted(set(side)))} for side in labels]
    partial = [np.zeros((len(index), len(index)), dtype=complex) for index in indices]
    for i in range(len(states)):
        for j in range(len(states)):
            if labels[1][i] == labels[1][j]:
                partial[0][indices[0][labels[0][i]], indices[0][labels[0][j]]] += rho[i, j]
            if labels[0][i] == labels[0][j]:
                partial[1][indices[1][labels[1][i]], indices[1][labels[1][j]]] += rho[i, j]
    return tuple(partial)


def _state_entropy(rho):
    np.testing.assert_allclose(rho, rho.conj().T, atol=2e-13)
    assert abs(np.trace(rho)-1) < 2e-13
    spectrum = np.linalg.eigvalsh(rho)
    assert min(spectrum) >= -2e-13
    positive = spectrum[spectrum > 0]
    return -float(np.dot(positive, np.log(positive)))


def test_native_region_entropy_is_shape_blind_while_crossing_seams_differ():
    from collections import Counter
    from model.primitive_observer.transfer import native_counts

    source = native_counts()
    regions = ((0, 1, 5, 7), (7, 8, 10, 11))
    tau = np.eye(66)/66
    entropy_pairs, boundaries = [], []
    for region in regions:
        reached = {region[0]}
        for _ in region:
            reached |= {b for a, b in source.seams if a in reached and b in region}
            reached |= {a for a, b in source.seams if b in reached and a in region}
        assert reached == set(region)
        boundaries.append(sum((a in region) != (b in region) for a, b in source.seams))
        counts = Counter(sum(state[i] for i in region) for state in source.states)
        assert [counts[j] for j in range(3)] == [28, 32, 6]
        left, right = _sector_partial_traces(tau, source.states, region)
        entropy_pairs.append((_state_entropy(left), _state_entropy(right)))
        q = np.array([28, 32, 6])/66
        center_entropy = -float(np.dot(q, np.log(q)))
        assert abs(sum(entropy_pairs[-1])-np.log(66)-center_entropy) < 2e-13
    assert boundaries == [10, 14]
    np.testing.assert_allclose(entropy_pairs[0], entropy_pairs[1], atol=2e-13)


@pytest.mark.parametrize('kind', ('cross_charge_coherence', 'nonuniform_coherent'))
def test_quantum_region_mutual_information_bound_uses_full_nondiagonal_states(kind):
    from model.primitive_observer.transfer import native_counts

    states = native_counts().states
    region = (0, 1, 5, 7)
    charges = np.array([sum(state[i] for i in region) for state in states])
    dimension = len(states)
    tau = np.eye(dimension, dtype=complex)/dimension
    if kind == 'cross_charge_coherence':
        i, j = int(np.flatnonzero(charges == 0)[0]), int(np.flatnonzero(charges == 1)[0])
        rho = tau.copy()
        rho[i, j] = 1j/(2*dimension); rho[j, i] = -1j/(2*dimension)
    else:
        ket = np.array([(index+1)*np.exp(1j*index/3) for index in range(dimension)])
        ket /= np.linalg.norm(ket)
        rho = (tau+np.outer(ket, ket.conj()))/2
    assert np.linalg.norm(rho-np.diag(np.diag(rho))) > 0.001
    left, right = _sector_partial_traces(rho, states, region)
    entropy = _state_entropy(rho)
    information = _state_entropy(left)+_state_entropy(right)-entropy
    deficit = np.log(dimension)-entropy
    p = np.array([np.real(np.diag(rho))[charges == charge].sum() for charge in range(3)])
    q = np.array([28, 32, 6])/66
    center_entropy = -float(np.dot(p, np.log(p)))
    center_deficit = float(np.dot(p, np.log(p/q)))
    assert 0 <= information <= center_entropy-center_deficit+deficit+3e-13
    assert information <= np.log(3)+deficit+3e-13
    if kind == 'cross_charge_coherence':
        # Both partial traces are unchanged; nonzero coherence entropy is
        # essential, and the sharper quantum bound is attained.
        assert deficit > 0
        np.testing.assert_allclose(p, q, atol=2e-14)
        assert abs(information-center_entropy-deficit) < 3e-13
        assert information > center_entropy


@pytest.mark.parametrize('particles', (2, 4))
def test_native_return_reference_entropy_deficit_for_ideal_and_all_retained_caps(particles):
    from fractions import Fraction as F
    from model.observer_time_fields import continuum_dynamics as dynamics
    from model.observer_time_fields import wave_reference as reference

    states, counts, edges = dynamics._native_trace_counts(particles)
    kernel = np.array(counts, dtype=float)/(2*len(edges))
    spectrum, vectors = np.linalg.eigh(kernel)
    assert spectrum[0] >= -2e-13 and spectrum[-1] <= 1+2e-13
    spectrum = np.clip(spectrum, 0, 1)
    kappa = np.trace(kernel)/len(states)
    dimension, denominator = len(states), 2*len(edges)
    exact_mean = F(sum(counts[i][i] for i in range(dimension)), dimension*denominator)
    exact_second = F(sum(value*value for row in counts for value in row), dimension*denominator**2)
    variance = exact_second-exact_mean**2
    assert abs(float(variance)-np.sum((kernel-kappa*np.eye(dimension))**2)/dimension) < 2e-14
    chi = F(particles*(12-particles), 12*11)
    upper = float(25*variance/chi**2)
    ideal = np.exp(-np.arccos(spectrum))
    for cap in (None, 0, 1, 4, 32):
        weights = ideal if cap is None else np.exp(-np.pi/2)*sum(
            float(c)*spectrum**n for n, c in enumerate(reference.coefficients(cap)))
        ratio = weights/np.mean(weights)
        deficit = float(np.dot(ratio, np.log(ratio))/len(states))
        squared = float(np.mean((ratio-1)**2))
        assert -2e-13 <= deficit <= squared+2e-13 <= upper+2e-13
        if cap == 0:
            assert abs(deficit) < 2e-13
        elif cap == 1 and particles == 2:
            rho = (vectors*(weights/weights.sum()))@vectors.T
            assert np.linalg.norm(rho-np.diag(np.diag(rho))) > 0.001
            left, right = _sector_partial_traces(rho, states, (0, 1, 5, 7))
            information = _state_entropy(left)+_state_entropy(right)-_state_entropy(rho)
            assert information <= np.log(3)+deficit+3e-13


@pytest.mark.parametrize('level', (0, 1, 2, 8, 100))
def test_all_level_reference_entropy_bound_is_absolute_and_cap_independent(level):
    from fractions import Fraction as F
    from model.observer_time_fields.continuum_dynamics import reference_state_bounds

    row = reference_state_bounds(level)
    ports, particles = row['ports'], row['particles']
    chi = F(particles*(ports-particles), ports*(ports-1))
    assert chi >= F(2, 9)
    assert 25*row['transfer_trace_variance']/chi**2 <= F(2025, 64*ports)
    assert F(2025, 64*ports) == F(135, 1024*4**level)


@pytest.mark.parametrize('level', (0, 1, 2))
def test_actual_screen_patch_counts_distinguish_cells_from_boundary_seams(level):
    from model.federated_observer.source import federation

    graph = federation(level)
    patch_ports = 12*4**level
    crossing = sum((a < patch_ports) != (b < patch_ports) for a, b in graph.seams)
    assert patch_ports == graph.ports//20
    assert crossing == 3*2**level
    assert crossing == sum((a < patch_ports) != (b < patch_ports) for a, b in graph.inter_seams)


def test_counting_entropy_marginal_bounds_follow_exact_binomial_mode_control():
    import math
    from fractions import Fraction as F

    for population in range(2, 13):
        for occupied in range(1, population):
            rho = F(occupied, population)
            mode_mass = math.comb(population, occupied)*rho**occupied*(1-rho)**(population-occupied)
            assert mode_mass >= F(1, population+1)
            h = -float(rho)*np.log(float(rho))-float(1-rho)*np.log(float(1-rho))
            for size in range(1, population):
                low, high = max(0, occupied-population+size), min(size, occupied)
                region_entropy = 0.0
                for number in range(low, high+1):
                    multiplicity = math.comb(size, number)
                    eigenvalue = F(math.comb(population-size, occupied-number), math.comb(population, occupied))
                    region_entropy -= multiplicity*float(eigenvalue)*np.log(float(eigenvalue))
                assert size*h-np.log(population+1)-2e-13 <= region_entropy <= size*h+2e-13


def test_area_gradient_and_hessian_match_finite_differences() -> None:
    n = 12
    d = 7.0
    point = np.full(n, d)

    def coordinate(dimensions: np.ndarray) -> float:
        return float(np.dot(dimensions, np.log(dimensions)) / np.sum(dimensions))

    epsilon = 1.0e-3
    basis0 = np.zeros(n)
    basis0[0] = 1.0
    numerical_gradient = (
        coordinate(point + epsilon * basis0)
        - coordinate(point - epsilon * basis0)
    ) / (2.0 * epsilon)
    assert numerical_gradient == pytest.approx(1.0 / (n * d), rel=2.0e-8)

    tangent = np.zeros(n)
    tangent[0] = 1.0
    tangent[1] = -1.0
    tangent_second = (
        coordinate(point + epsilon * tangent)
        - 2.0 * coordinate(point)
        + coordinate(point - epsilon * tangent)
    ) / epsilon**2
    assert tangent_second == pytest.approx(2.0 / (n * d**2), rel=1.0e-6)

    homogeneous = np.ones(n)
    homogeneous_second = (
        coordinate(point + epsilon * homogeneous)
        - 2.0 * coordinate(point)
        + coordinate(point - epsilon * homogeneous)
    ) / epsilon**2
    assert homogeneous_second == pytest.approx(-1.0 / d**2, rel=1.0e-6)


def test_one_sided_capacity_transfer_is_monotone_and_matches_budget_rows() -> None:
    transfer = build_receipt()["entropy_and_capacity"]["one_sided_capacity_transfer"]
    finite = transfer["finite_integer_statement"]
    analytic = transfer["positive_real_analytic_interpolation"]
    assert "positive integer" in finite["domain"]
    assert finite["f_equals_zero_boundary_maximum"] is True
    assert finite["monotone_on_admissible_depletions"] is True
    assert analytic["f_equals_zero_one_sided_boundary_maximum"] is True
    assert analytic["right_derivative_at_zero"] == -1.0
    assert analytic["monotone_coordinate_loss"] is True
    rows = finite["budget_checks"]
    assert [row["observer_capacity_fraction"] for row in rows] == [0.02, 0.05, 0.1]
    assert [row["coordinate_loss"] for row in rows] == [
        pytest.approx(0.020202707318),
        pytest.approx(0.051293294388),
        pytest.approx(0.105360515658),
    ]
    assert all(
        row["coordinate_change"] < 0.0 < row["coordinate_loss"] for row in rows
    )
    assert [row["initial_integer_sector_dimension"] for row in rows] == [100, 100, 100]
    assert [row["depleted_integer_sector_dimension"] for row in rows] == [98, 95, 90]


def test_port_edge_and_face_graph_spectra_are_frozen() -> None:
    raw = build_receipt()["graph_spectra"]["raw_laplacian_spectra"]
    assert [
        (row["raw_laplacian_eigenvalue_exact"], row["multiplicity"])
        for row in raw["ports"]
    ] == [
        ("0", 1),
        ("5 - sqrt(5)", 3),
        ("6", 5),
        ("5 + sqrt(5)", 3),
    ]
    assert [
        (row["raw_laplacian_eigenvalue_exact"], row["multiplicity"])
        for row in raw["edge_sectors"]
    ] == [
        ("0", 1),
        ("5 - sqrt(5)", 3),
        ("6", 5),
        ("5 + sqrt(5)", 3),
        ("10", 18),
    ]
    assert [
        (row["raw_laplacian_eigenvalue_exact"], row["multiplicity"])
        for row in raw["faces"]
    ] == [
        ("0", 1),
        ("3 - sqrt(5)", 3),
        ("2", 5),
        ("3", 4),
        ("5", 4),
        ("3 + sqrt(5)", 3),
    ]


def test_line_graph_theorem_is_checked_by_exact_integer_incidence_identities() -> None:
    theorem = build_receipt()["graph_spectra"]["line_graph_theorem"]
    assert theorem["port_identity_max_integer_residual"] == 0
    assert theorem["edge_identity_max_integer_residual"] == 0
    assert theorem["incidence_rank"] == 12
    assert theorem["extra_eigenvalue_exact"] == "2k=10"
    assert theorem["extra_multiplicity_exact"] == "m-n=18"
    assert theorem["port_edge_low_spectra_identical"] is True


def test_line_graph_negative_control_breaks_after_one_edge_is_deleted() -> None:
    graphs = _icosahedron_graphs()
    mutated = graphs["edge_adjacency"].copy()
    pair = np.argwhere(np.triu(mutated, k=1) == 1)[0]
    mutated[pair[0], pair[1]] = 0
    mutated[pair[1], pair[0]] = 0
    original_eigenvalues = np.linalg.eigvalsh(_laplacian(graphs["edge_adjacency"]))
    mutated_eigenvalues = np.linalg.eigvalsh(_laplacian(mutated))
    assert not np.allclose(original_eigenvalues, mutated_eigenvalues)
    assert not np.all(np.sum(mutated, axis=1) == 8)


def test_normalized_spectrum_is_conditional_on_both_named_premises() -> None:
    receipt = build_receipt()
    boundary = receipt["claim_boundary"]
    conditional = receipt["graph_spectra"][
        "conditional_normalized_port_shock_spectrum"
    ]
    assert boundary["assumption_tokens"] == ["DS-GAUGE", "DS-LAPLACIAN"]
    assert receipt["assumptions"]["DS-GAUGE"]["discharged"] is False
    assert receipt["assumptions"]["DS-LAPLACIAN"]["discharged"] is False
    assert boundary["physical_shock_spectrum_established"] is False
    assert boundary["promotion_allowed"] is False
    assert conditional["status"] == (
        "DERIVED_CONDITIONAL_ON_DS-GAUGE_AND_DS-LAPLACIAN"
    )
    assert conditional["scope_boundary"]["physical_prediction"] is False
    assert len(conditional["scope_boundary"]["escape_routes"]) == 3


def test_primary_source_and_physical_interpretation_boundary_are_explicit() -> None:
    receipt = build_receipt()
    source = receipt["primary_source"]
    interpretation = receipt["physical_interpretation"]
    assert source["arxiv_id"] == "2607.14042v1"
    assert "Eq. (C.9)" in " ".join(
        source["equation_map"]["higher_dimensional_shock_operator"]
    )
    assert interpretation["status"] == (
        "CONDITIONAL_INTERPRETATION_WITH_OPEN_PHYSICAL_PREMISES"
    )
    assert interpretation["premise_tokens"] == [
        "capacity_to_horizon_area",
        "capacity_ledger_to_observer_mass",
        "einstein_branch_realization",
        "coefficient_and_physical_scale",
    ]
    assert all(
        row["discharged"] is False
        for row in interpretation["premises"].values()
    )
    assert interpretation["physical_time_advance_explained"] is False
    assert interpretation["promotion_allowed"] is False
    assert len(interpretation["scope_boundary"]["escape_routes"]) == 4


def test_invariant_shock_mode_is_exactly_minus_two() -> None:
    conditional = build_receipt()["graph_spectra"][
        "conditional_normalized_port_shock_spectrum"
    ]
    spectrum = conditional["spectrum"]
    assert spectrum[0]["block"] == "1"
    assert spectrum[0]["shock_eigenvalue_exact"] == "-2"
    assert spectrum[0]["shock_eigenvalue"] == -2.0
    assert spectrum[1]["block"] == "3"
    assert spectrum[1]["shock_eigenvalue"] == 0.0
    assert conditional["invariant_minus_two_guard"]["passed"] is True


def test_repair_generator_domain_does_not_contain_sector_dimension_shock() -> None:
    boundary = build_receipt()["repair_generator_domain_boundary"]
    assert boundary["repair_domain"] == "K_r=L^2(X_r,pi_r) at fixed sector structure"
    assert boundary["shock_deformation_in_repair_domain"] is False
    assert boundary["repair_gap_constrains_shock_mode"] is False
    assert boundary["repair_gap_forbids_shock_mode"] is False
    assert boundary["repair_gap_protects_shock_mode"] is False
    assert boundary["repair_kernel_exactness_established_by_this_receipt"] is False
    assert boundary["shock_generator_supplied_by_repair_operator"] is False
    assert boundary["classification"] == "domain_separation_with_shock_generator_open"


@pytest.mark.parametrize(
    ("path", "value", "failure"),
    [
        (
            ("claim_boundary", "promotion_allowed"),
            True,
            "promotion_allowed_must_be_false",
        ),
        (
            ("claim_boundary", "physical_shock_spectrum_established"),
            True,
            "physical_shock_spectrum_established_must_be_false",
        ),
        (
            ("assumptions", "DS-GAUGE", "discharged"),
            True,
            "DS-GAUGE_must_remain_undischarged",
        ),
        (
            ("assumptions", "DS-LAPLACIAN", "discharged"),
            True,
            "DS-LAPLACIAN_must_remain_undischarged",
        ),
        (
            ("entropy_and_capacity", "area_analytic_relaxation", "unconstrained_stationary"),
            True,
            "symmetric_area_point_must_be_nonstationary",
        ),
        (
            (
                "entropy_and_capacity",
                "area_analytic_relaxation",
                "fixed_M_tangent_classification",
            ),
            "strict_local_maximum",
            "fixed_M_tangent_sign_mismatch",
        ),
        (
            (
                "graph_spectra",
                "conditional_normalized_port_shock_spectrum",
                "invariant_minus_two_guard",
                "observed",
            ),
            -2.76,
            "minus_two_observed_mismatch",
        ),
        (
            (
                "repair_generator_domain_boundary",
                "shock_deformation_in_repair_domain",
            ),
            True,
            "shock_deformation_in_repair_domain_must_be_false",
        ),
        (
            (
                "repair_generator_domain_boundary",
                "repair_kernel_exactness_established_by_this_receipt",
            ),
            True,
            "repair_kernel_exactness_established_by_this_receipt_must_be_false",
        ),
        (
            ("physical_interpretation", "physical_time_advance_explained"),
            True,
            "physical_interpretation_physical_time_advance_explained_must_be_false",
        ),
        (
            (
                "physical_interpretation",
                "premises",
                "capacity_to_horizon_area",
                "discharged",
            ),
            True,
            "capacity_to_horizon_area_must_remain_undischarged",
        ),
    ],
)
def test_adversarial_status_and_sign_mutations_fail_closed(
    path: tuple[str, ...],
    value: object,
    failure: str,
) -> None:
    mutated = copy.deepcopy(build_receipt())
    cursor = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    assert failure in receipt_failures(mutated)
    with pytest.raises(ReceiptValidationError, match=failure):
        validate_receipt(mutated)


def test_adversarial_capacity_loss_and_line_graph_mutations_fail_closed() -> None:
    mutated_loss = copy.deepcopy(build_receipt())
    mutated_loss["entropy_and_capacity"]["one_sided_capacity_transfer"][
        "finite_integer_statement"
    ][
        "budget_checks"
    ][0]["coordinate_loss"] *= -1.0
    assert "capacity_budget_losses_mismatch" in receipt_failures(mutated_loss)

    mutated_line = copy.deepcopy(build_receipt())
    mutated_line["graph_spectra"]["line_graph_theorem"][
        "edge_identity_max_integer_residual"
    ] = 1
    assert (
        "edge_identity_max_integer_residual_must_be_zero"
        in receipt_failures(mutated_line)
    )
