"""Exact-envelope replay, independent geometry, and false-green controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_quantum_history as producer
import verify_whitney_quantum_history as audit


@pytest.fixture(scope="module")
def packet():
    return audit.load()


@pytest.fixture(scope="module")
def elements():
    return audit.exact_geometry()[0]


def test_committed_receipt_replays_fresh_and_is_canonical(packet):
    assert producer.OUTPUT.read_bytes() == producer.canonical(producer.build())
    result = audit.verify(deepcopy(packet))
    assert result["accepted"] is True
    assert result["horizon"] == "1/36028797018963968"
    assert Q(result["squared_norm_error_upper"]) < Q(1, 100)
    assert result["real_configuration_dimension"] == 56
    assert result["trial_history_computed"] is True
    assert result["exact_Hamiltonian_history_computed"] is False


def test_verifier_runs_without_producer_or_sibling_import_path():
    verifier = Path(audit.__file__).resolve()
    script = """
import importlib.abc, importlib.util, json, pathlib, sys
class BlockProducer(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.rsplit('.', 1)[-1] in {
            'whitney_quantum_history', 'whitney_interacting_quantum',
            'whitney_quantum_state'}:
            raise AssertionError('verifier imported a producer: '+fullname)
sys.meta_path.insert(0, BlockProducer())
source=pathlib.Path(sys.argv[1])
sys.path=[p for p in sys.path if p and pathlib.Path(p).resolve()!=source.parent]
spec=importlib.util.spec_from_file_location('isolated_history_audit', source)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
print(json.dumps(module.verify(module.load()), sort_keys=True))
"""
    result = subprocess.run([sys.executable, "-c", script, str(verifier)],
        cwd=producer.ROOT, env={**os.environ, "PYTHONPATH": ""},
        capture_output=True, text=True, check=True, timeout=30)
    assert json.loads(result.stdout)["accepted"] is True


@pytest.mark.parametrize("a,b,sign", [
    (0, 0, 0), (0, 1, 1), (0, -1, -1), (3, -1, 1),
    (2, -1, -1), (-3, 1, -1), (-2, 1, 1),
    (Q(161, 72), -1, 1), (Q(682, 305), -1, -1),
    (Q(-161, 72), 1, -1), (Q(-682, 305), 1, 1),
])
def test_exact_irrational_order_has_no_decimal_rounding(a, b, sign):
    assert audit.Golden(a, b).sign() == sign


def test_exact_simplex_geometry_and_monomials(elements):
    volume = sum((v for _, v, _ in elements), audit.Golden())
    assert volume.pair() == ["10", "10/3"]
    assert len(elements) == 20
    for _, cell_volume, gradients in elements:
        assert (cell_volume-Q(5, 6)).sign() >= 0
        assert all((audit.Golden(4)-audit.dot(g, g)).sign() >= 0 for g in gradients)
        assert all(sum((g[c] for g in gradients), audit.Golden()).sign() == 0 for c in range(3))
    assert audit.simplex_average((1, 0, 0, 0)) == Q(1, 4)
    assert audit.simplex_average((2, 0, 0, 0)) == Q(1, 10)
    assert audit.simplex_average((1, 1, 0, 0)) == Q(1, 20)
    assert audit.simplex_average((4, 0, 0, 0)) == Q(1, 35)


@pytest.mark.parametrize("dimension", [26, 30])
@pytest.mark.parametrize("power", range(6))
def test_radial_moments_match_independent_gaussian_pairing(dimension, power):
    # Wick expansion of (sum_i Z_i^2)^power uses independent standard normal
    # even moments. Its generating series is multiplied across coordinates.
    from math import factorial
    series = [Q(1)]
    even = [Q(factorial(2*k), 2**k*factorial(k)*factorial(k)) for k in range(power+1)]
    for _ in range(dimension):
        series = [sum((series[j]*even[k-j] for j in range(min(k, len(series)-1)+1)), Q(0))
                  for k in range(power+1)]
    wick = factorial(power)*series[power]/4**power
    assert producer.gaussian_radial_moment(dimension, power) == wick == audit.moment(dimension, power)


def test_constants_derive_from_envelopes_and_full_dimension():
    polys, bounds, ell0, ell1 = audit.independent_bounds()
    assert (ell0, ell1) == (1896, 11286)
    assert polys == producer.envelopes()
    assert bounds == producer.energy_bounds()
    assert bounds["kinetic_upper"] == 135829605095437
    assert bounds["potential_upper"] == Q(619353, 16)
    # A five-coordinate replacement would change the Gaussian moments and
    # cannot serve as the certified full-configuration calculation.
    wrong = sum((c*producer.gaussian_radial_moment(3, x)*producer.gaussian_radial_moment(2, y)
                 for (x, y), c in polys["initial_kinetic"].items()), Q(0))
    assert wrong != bounds["kinetic_upper"]


def test_exact_phase_replay_uses_full_simplex_action(elements):
    configs = producer.phase_configurations()
    for row in configs:
        value = audit.sampled_potential([Q(x) for x in row["coordinates_exact"]], elements)
        assert value.pair() == row["potential_in_Qsqrt5"]
    assert configs[0]["potential_in_Qsqrt5"] == ["0", "0"]
    assert configs[1]["potential_in_Qsqrt5"] == configs[2]["potential_in_Qsqrt5"]
    # A spatially constant complex field has zero gradient and exact scalar
    # potential Vol*(m²*|psi|²+g*|psi|⁴/2), including the quartic factor1/2.
    q = [Q(0)]*30+[Q(1, 3)]*13+[Q(1, 4)]*13
    norm = Q(1, 9)+Q(1, 16)
    expected = audit.Golden(10, Q(10, 3))*(norm/2+norm**2/8)
    assert audit.sampled_potential(q, elements).pair() == expected.pair()


def test_dyadic_error_bound_is_two_sided_and_certifies_whole_interval():
    bounds = producer.energy_bounds()
    power, time = producer.selected_horizon(bounds)
    assert power == 55 and time < Q(1, 10**16)
    assert producer.squared_error_bound(time, bounds) <= Q(1, 100)
    assert producer.squared_error_bound(2*time, bounds) > Q(1, 100)
    for factor in [Q(0), Q(1, 7), Q(2, 3), Q(1)]:
        assert producer.squared_error_bound(-factor*time, bounds) == producer.squared_error_bound(factor*time, bounds)
        assert producer.squared_error_bound(factor*time, bounds) <= producer.squared_error_bound(time, bounds)


def change(packet, path, value):
    row = packet
    for key in path[:-1]:
        row = row[key]
    row[path[-1]] = value


@pytest.mark.parametrize("path,value", [
    (("state", "dimension"), 5), (("state", "dimension"), 56.0),
    (("state", "trial_history_computed"), 1),
    (("state", "exact_Hamiltonian_history_computed"), True),
    (("state", "configuration_density_moves"), True),
    (("state", "physical_state_preparation"), True),
    (("state", "observer_history"), True),
    (("state", "ordinary_physics_benchmark"), True),
    (("state", "classical_five_coordinate_restriction"), True),
    (("state", "measure"), "Lebesgue dq"),
    (("geometry_envelope", "minimum_node_multiplicity"), 5.0),
    (("geometry_envelope", "Maxwell_mass_lower"), "1/72"),
    (("geometry_envelope", "dressed_mass_lower"), "1/64"),
    (("coefficient_bounds", "log_density_gradient_upper"), "1896+11285*sqrt(Y)"),
    (("bounds", "kinetic_upper"), "135829605095436"),
    (("bounds", "force_moment_upper"), "24441816695556086/4"),
    (("integration", "Gaussian_tail_truncation"), True),
    (("integration", "Monte_Carlo_used"), True),
    (("integration", "quadrature_used_for_bound"), True),
    (("integration", "floating_point_used_for_bound"), True),
    (("error_certificate", "dyadic_horizon_power"), 54),
    (("error_certificate", "horizon"), "1/18014398509481984"),
    (("error_certificate", "target_norm_error"), "1/100"),
    (("error_certificate", "trivial_norm_error_upper"), "1"),
    (("error_certificate", "coverage"), "stored samples only"),
    (("error_certificate", "time_units"), "seconds"),
    (("history", 1, "squared_norm_error_upper"), "0"),
    (("history", 1, "time"), "1/4"),
])
def test_mutated_scientific_claims_fail_closed(packet, path, value):
    mutated = deepcopy(packet)
    change(mutated, path, value)
    with pytest.raises(ValueError):
        audit.verify(mutated)


def test_coherent_wrong_phase_sign_cannot_pass_exact_replay(packet):
    mutated = deepcopy(packet)
    for row in mutated["history"]:
        for sample in row["phase_samples"]:
            sample["angle_in_Qsqrt5"] = [str(-Q(x)) for x in sample["angle_in_Qsqrt5"]]
    with pytest.raises(ValueError, match="history phase"):
        audit.verify(mutated)


def test_coherent_wrong_potential_and_corresponding_phases_fail(packet):
    mutated = deepcopy(packet)
    for row in mutated["phase_configurations"]:
        row["potential_in_Qsqrt5"] = [str(2*Q(x)) for x in row["potential_in_Qsqrt5"]]
    for row in mutated["history"]:
        for sample in row["phase_samples"]:
            sample["angle_in_Qsqrt5"] = [str(2*Q(x)) for x in sample["angle_in_Qsqrt5"]]
    with pytest.raises(ValueError, match="independent exact sample potential"):
        audit.verify(mutated)


def test_removed_sample_and_source_pin_fail(packet):
    mutated = deepcopy(packet)
    mutated["history"].pop()
    with pytest.raises(ValueError, match="schedule"):
        audit.verify(mutated)
    mutated = deepcopy(packet)
    mutated["source_pins"][next(iter(mutated["source_pins"]))] = "0"*64
    with pytest.raises(ValueError, match="source pin"):
        audit.verify(mutated)


@pytest.mark.parametrize("value", [True, 1.0, "2/2", "01", "1/0", "NaN", "1e-9"])
def test_noncanonical_rationals_rejected(value):
    with pytest.raises(ValueError):
        audit.rational(value)


@pytest.mark.parametrize("data", ['{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":1e9999}'])
def test_json_duplicate_and_nonfinite_values_rejected(tmp_path, data):
    path = tmp_path/"bad.json"
    path.write_text(data, encoding="utf-8")
    with pytest.raises(ValueError):
        audit.load(path)


@pytest.mark.parametrize("packet", [None, [], 0, {"schema": "oph.whitney_quantum_trial_history.v1"}])
def test_malformed_packet_containers_rejected(packet):
    with pytest.raises(ValueError):
        audit.verify(packet)


@pytest.mark.parametrize("time", [True, "0", complex(1, 0), float("nan"), float("inf"), Q(10**1000)])
def test_numerical_trial_rejects_invalid_time(time):
    with pytest.raises(ValueError, match="time"):
        producer.trial_log_state(np.zeros(56), time)


def test_numerical_trial_keeps_full_density_and_exact_phase_law():
    mesh = producer.quantum.geometry(5)
    q = np.arange(56, dtype=float)/300
    times = [0, Q(1, 7), Q(2, 7), Q(3, 7)]
    states = [producer.trial_log_state(q, t, mesh) for t in times]
    assert all(s["log_magnitude_numeric"] == states[0]["log_magnitude_numeric"] for s in states)
    assert states[0]["phase_numeric"] == 0
    assert abs(states[3]["phase_numeric"]-states[1]["phase_numeric"]-states[2]["phase_numeric"]) < 1e-13
    assert states[1]["phase_numeric"] < 0
    assert "not used to certify" in states[1]["scope"]
    angle = .713
    scalar = (q[30:43]+1j*q[43:])*np.exp(1j*angle)
    rotated = np.concatenate([q[:30], scalar.real, scalar.imag])
    other = producer.trial_log_state(rotated, times[1], mesh)
    assert abs(other["phase_numeric"]-states[1]["phase_numeric"]) < 1e-12
    assert abs(other["log_magnitude_numeric"]-states[1]["log_magnitude_numeric"]) < 1e-12
