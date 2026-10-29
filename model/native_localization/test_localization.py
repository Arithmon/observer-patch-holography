"""Independent physical-interface and exact-algebra checks for localization."""
from fractions import Fraction as F
from itertools import combinations
import math

import numpy as np
import pytest
import sympy as sp

from . import localization as loc
from ..common_source import record_code
from ..local_repair_family import local_bias as local
from ..observer.source_discrimination import circuit


def test_faithful_same_reference_diagonal_and_full_exterior():
    rho = loc.initial_density()
    pi = loc.initial_density(F())
    assert all(rho[i][i] == pi[i][i] for i in range(66))
    assert loc.exterior_restriction(rho) == loc.exterior_restriction(pi)
    assert np.linalg.eigvalsh(np.array(rho, float)).min() > 0
    assert np.trace(np.array(rho, float)) == pytest.approx(1)
    for i, j in loc.TRANSFERS:
        assert rho[i][i]*rho[j][j]-rho[i][j]**2 > 0
    assert loc.contrast(0) == F(776, 7505)
    assert np.linalg.norm(np.array(rho, float)-np.array(pi, float), ord="nuc")/2 == pytest.approx(float(loc.contrast(0)))


def test_unphysical_coherence_rejected_not_projected():
    with pytest.raises(ValueError):
        loc.initial_density(F(15))
    with pytest.raises(ValueError):
        loc.initial_density(.5)
    assert loc.initial_density(F(-8))[loc.TRANSFERS[0][0]][loc.TRANSFERS[0][1]] < 0


def test_all_raw_branches_preserve_difference_and_exterior_including_record():
    report = loc.trajectory_certificate()
    assert report["raw_slots"] == 720
    assert report["matrix_units_examined"] == 7200
    assert report["raw_abort_slots"] == 240
    # Independently check the exterior trace condition on every surviving
    # endpoint image without using trajectory_certificate's transfer index.
    for law in ("C", "E"):
        for seam in local.SEAMS:
            for coin, _, mapping in local.raw_maps(law, seam):
                mapping = dict(mapping)
                for i, j in loc.TRANSFERS:
                    if i in mapping and j in mapping:
                        x, y = loc.BASIS[mapping[i]], loc.BASIS[mapping[j]]
                        assert tuple(x[k] for k in loc.CORE) != tuple(y[k] for k in loc.CORE)
                        assert tuple(a-b for a, b in zip(x, y, strict=True)) == tuple(
                            a-b for a, b in zip(loc.BASIS[i], loc.BASIS[j], strict=True))


@pytest.mark.parametrize("law", ["C", "E", "mixed"])
def test_exact_killed_kernel_against_full_density_channel(law):
    p = loc.parameters(law)
    rho = np.array(loc.initial_density(), float)
    output = np.zeros((66, 66))
    for name, mixture in (("C", p.c_weight), ("E", 1-p.c_weight)):
        for seam in local.SEAMS:
            for coin, _, mapping in local.raw_maps(name, seam):
                source, target = zip(*mapping) if mapping else ((), ())
                output[np.ix_(target, target)] += float(mixture*local.coin_probability(p, seam, coin)/30)*rho[np.ix_(source, source)]
    v = [F(8*p.port_weights[k], 7505) for k in loc.SPECTATORS]
    expected = loc._step(v, loc.coherence_kernel(law))
    for value, (i, j) in zip(expected, loc.TRANSFERS, strict=True):
        assert output[i, j] == pytest.approx(float(value), abs=1e-15)
    assert np.diag(output) == pytest.approx(np.diag(rho), abs=1e-15)
    assert np.trace(output) == pytest.approx(1)
    assert np.linalg.norm(output-np.diag(np.diag(rho)), ord="nuc")/2 == pytest.approx(float(loc.contrast(1, law)))


@pytest.mark.parametrize("law", ["C", "E", "mixed"])
def test_exact_resolvent_lifetime_and_positive_real_spectrum(law):
    report = loc.spectral_certificate(law)
    q = sp.Matrix(loc.coherence_kernel(law))
    h = sp.Matrix([sp.Rational(x) for x in report["potential"]])
    assert h == sp.ones(10, 1)+q*h
    assert float(F(report["integrated_relative_visibility_attempts"])) > 1
    low, high = map(F, report["all_eigenvalues_rational_enclosure"])
    assert 0 < low < high < 1
    eigen = np.linalg.eigvals(np.array(q, float))
    assert np.max(abs(eigen.imag)) < 1e-13
    assert min(eigen.real) >= float(low)-1e-13
    assert max(eigen.real) <= float(high)+1e-13
    assert loc.contrast(32, law) < loc.contrast(16, law) < loc.contrast(0, law)


def test_E_lifetime_and_deleting_incident_seams_countercontrol():
    assert F(loc.spectral_certificate("E")["integrated_relative_visibility_attempts"]) == F(10, 3)
    for n in (0, 1, 2, 16, 32):
        assert loc.contrast(n, "E") == loc.contrast(0, "E")*F(7, 10)**n
    for seam in local.SEAMS:
        if not set(seam) & set(loc.CORE):
            for i, j in loc.TRANSFERS:
                surviving = sum(local.coin_probability(loc.parameters(), seam, coin)
                                for coin, _, mapping in local.raw_maps("E", seam)
                                if i in dict(mapping) and j in dict(mapping))
                assert surviving == 1
    # If all incident seams are removed, the positive all-seam contract is
    # lost and this coherence has no killing. The report excludes that escape.


def test_exterior_obstruction_beyond_sharp_q_and_its_diagonal_attainer():
    report = loc.exact_code_certificate()
    pi = local.reference(loc.parameters())
    seams = set(local.SEAMS)
    minima = []
    for row in report["rows"]:
        core = tuple(row["core"])
        outside = set(range(12))-set(core)
        groups = {}
        for i, pair in enumerate(loc.PAIRS):
            groups.setdefault(tuple(k for k in pair if k in outside), []).append(i)
        target = {key: sum((pi[i] for i in indices), F()) for key, indices in groups.items()}
        allowed = {key: [i for i in indices if loc.PAIRS[i] not in seams] for key, indices in groups.items()}
        missing = sum((target[key] for key, indices in allowed.items() if not indices), F())
        assert missing == F(row["universal_exterior_trace_distance_lower"])
        # A diagonal nonedge state attains this lower bound. It need not be a code.
        actual = {key: target[key]/(1-missing) if indices else F() for key, indices in allowed.items()}
        assert sum(actual.values()) == 1
        assert sum(abs(actual[key]-target[key]) for key in target)/2 == missing
        minima.append(missing)
    assert len(minima) == 66
    assert min(minima) == F(report["minimum_any_two_port_core"]) > 0


def test_arbitrary_six_level_code_has_exact_scalar_clock_and_predicted_leakage():
    code = record_code.ideal_code_matrix()[list(loc.SOURCE_CODES)]
    heavy = np.array([x[0]+x[1] for x in loc.BASIS])
    compressed = code.conj().T@(heavy[:, None]*code)
    nu = float(np.trace(compressed).real/6)
    assert compressed == pytest.approx(nu*np.eye(6), abs=1e-15)
    assert nu == pytest.approx(F(1, 3))
    parameter = .713
    phase = np.exp(-1j*parameter*math.log(25/8)*heavy)
    overlap = code.conj().T@(phase[:, None]*code)
    expected = (1-nu+nu*np.exp(-1j*parameter*math.log(25/8)))*np.eye(6)
    assert overlap == pytest.approx(expected, abs=1e-15)
    leakage = 1-np.linalg.norm(overlap[:, 0])**2
    assert leakage == pytest.approx(4*nu*(1-nu)*math.sin(parameter*math.log(25/8)/2)**2)


def test_preparation_factor_and_zero_energy_are_not_an_arbitrary_qubit_code():
    factor, _ = loc.preparation_factor()
    on_shell = factor[list(loc.SOURCE_CODES)]
    assert on_shell@on_shell.conj().T == pytest.approx(np.array(loc.initial_density(), float), abs=1e-15)
    assert np.linalg.norm(factor)**2 == pytest.approx(1)
    state = loc.state_certificate()
    assert state["modular_energy_excess"] == "0"
    assert state["entropy_deficit_float_diagnostic"] > 0
    assert state["physical_topological_charge_or_particle_identified"] is False


def test_compiled_source_preparation_clock_and_every_two_attempt_raw_history():
    report = loc.source_execution()
    assert report["complete_raw_history_count"] == 128
    assert report["full_history_error_sum_frobenius"] < 1e-8
    assert report["complete_history_mass"] == pytest.approx(1, abs=1e-10)
    assert report["raw_abort_mass"] < 1e-20
    assert report["formal_outward_or_hardware_certificate"] is False


def test_current_response_uses_actual_native_detector_and_all_spectators():
    rho = np.array(loc.initial_density(), complex)
    phase = np.diag(np.exp(-1j*.73*np.log(np.array(local.reference(loc.parameters()), float))))
    rotated = phase@rho@phase.conj().T
    analyzer = np.eye(66)
    for i, j in loc.TRANSFERS:
        analyzer[np.ix_((i, j), (i, j))] = np.array([[1., 1.], [-1., 1.]])/math.sqrt(2)
    analyzed = analyzer@rotated@analyzer.T
    expected = loc.current_response(0, "mixed", .73)
    for name in ("C", "E"):
        probabilities = {str(j): 0. for j in (-1, 0, 1)}
        for coin, value, mapping in local.raw_maps(name, loc.CORE):
            current = local.records.record_current(name, (*loc.CORE, coin, value))
            probabilities[str(current)] += float(local.coin_probability(loc.parameters(), loc.CORE, coin))*sum(analyzed[i, i].real for i, _ in mapping)
        assert probabilities == pytest.approx(expected["current_probabilities"], abs=1e-15)
    assert loc.current_response(8, "mixed", 1.)["cosine_coefficient"] != expected["cosine_coefficient"]


def test_energy_carrying_unitary_control_does_not_retain_reference_exterior():
    report = loc.unitary_preparation_control()
    assert report["rows"][0]["exterior_trace_distance"] == "0"
    assert F(report["rows"][1]["exterior_trace_distance"]) > 0
    assert report["rows"][0]["modular_energy_excess_coefficient_of_log25over8"] == "1649/15010"
    assert F(report["rows"][-1]["modular_energy_excess_coefficient_of_log25over8"]) < F(1649, 15010)


@pytest.fixture(scope="module")
def retained_packet():
    from . import experiment
    return experiment.report()


@pytest.mark.parametrize("mutation", ["probability", "current_mean", "source_count", "source_budget", "eigenvalue", "entropy", "nonfinite", "physical_boundary"])
def test_retained_receipt_rejects_forged_observables_budgets_and_nonfinite(retained_packet, mutation):
    from copy import deepcopy
    from . import experiment
    changed = deepcopy(retained_packet)
    if mutation == "probability":
        changed["current_response"][0]["current_probabilities"]["1"] += .01
    elif mutation == "current_mean":
        changed["current_response"][0]["mean_current"] += .01
    elif mutation == "source_count":
        changed["source_execution"]["retained_preparation_reference_columns"] = 65
    elif mutation == "source_budget":
        changed["source_execution"]["preparation_pulses"] += 1
    elif mutation == "eigenvalue":
        changed["native_lifetime"]["C"]["eigenvalues_float_diagnostic"][0] += .01
    elif mutation == "entropy":
        changed["state"]["entropy_deficit_float_diagnostic"] += .01
    elif mutation == "nonfinite":
        changed["source_execution"]["complete_history_mass"] = float("nan")
    else:
        changed["canonical_physical_premises_discharged"] = 1
    with pytest.raises(ValueError):
        experiment.verify(changed, retained_packet)


def test_retained_receipt_full_verification_and_small_float_tolerance(retained_packet):
    from copy import deepcopy
    from . import experiment
    experiment.verify(retained_packet, retained_packet)
    changed = deepcopy(retained_packet)
    changed["current_response"][0]["mean_current"] += 1e-12
    experiment.verify(changed, retained_packet)


@pytest.mark.parametrize("mutation", ["probability", "source_count", "nonfinite"])
def test_check_cli_rejects_hostile_receipt_path(retained_packet, tmp_path, monkeypatch, mutation):
    from copy import deepcopy
    import json
    import sys
    from . import experiment
    changed = deepcopy(retained_packet)
    if mutation == "probability":
        changed["current_response"][0]["current_probabilities"]["1"] += .01
    elif mutation == "source_count":
        changed["source_execution"]["preparation_pulses"] += 1
    else:
        changed["current_response"][0]["mean_current"] = float("nan")
    path = tmp_path/"hostile.json"
    path.write_text(json.dumps(changed))
    monkeypatch.setattr(experiment, "report", lambda: retained_packet)
    monkeypatch.setattr(sys, "argv", ["native_localization", "--check", str(path)])
    with pytest.raises(ValueError):
        experiment.main()
