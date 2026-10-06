"""Analytic recovery controls; no Qiskit installation or provider call needed."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "programs" / "stage1_markov_fingerprint.py"
)
SPEC = importlib.util.spec_from_file_location("stage1_markov_fingerprint", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
stage1 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage1)


@pytest.mark.parametrize(
    ("occupied", "cmi_bits", "optimal_squared_fidelity"),
    [
        ([0], 0.0, 1.0),  # Product |000>.
        ([0, 7], 1.0, 0.5),  # GHZ: AB loses the coherence between A=0 and A=1.
        ([0, 5], 2.0, 0.25),  # Bell pair AC, with independent B=0.
    ],
)
def test_recovery_bound_matches_analytic_optimal_squared_fidelity(
    occupied, cmi_bits, optimal_squared_fidelity
):
    psi = np.zeros(8, dtype=complex)
    psi[occupied] = 1 / np.sqrt(len(occupied))
    rho = np.outer(psi, psi.conj())
    recovered = stage1.petz_recovery(rho)
    row = stage1.analyze_state(rho)

    # Independent control: squared fidelity to a pure state is <psi|sigma|psi>.
    # For GHZ, every B-only recovery remains block diagonal in A and therefore
    # has overlap <= 1/2. For Bell_AC, AB = I_A/2 tensor |0><0|_B, so every
    # recovered state is I_A/2 tensor sigma_BC and has overlap <= 1/4.
    # The computed Petz states attain these bounds; the product case attains 1.
    overlap = float(np.vdot(psi, recovered @ psi).real)
    assert overlap == pytest.approx(optimal_squared_fidelity, abs=1e-14)
    assert row["petz_fidelity"] == pytest.approx(overlap, abs=1e-7)
    assert row["cmi_bits"] == pytest.approx(cmi_bits, abs=1e-14)
    assert row["fawzi_renner_fidelity_lower_bound"] == pytest.approx(
        optimal_squared_fidelity, abs=1e-14
    )
    assert row["fidelity_convention"] == "squared_uhlmann"
    assert row["fawzi_renner_bound_scope"] == "optimal_recovery_over_B_to_BC_channels"


def test_diagonal_markov_state_is_recovered_exactly():
    # p(a,b,c) = p(b) p(a|b) p(c|b), with asymmetric full-support marginals.
    p_b = np.array([0.3, 0.7])
    p_a_given_b = np.array([[0.2, 0.8], [0.6, 0.4]])
    p_c_given_b = np.array([[0.9, 0.1], [0.3, 0.7]])
    probabilities = np.array([
        p_b[b] * p_a_given_b[b, a] * p_c_given_b[b, c]
        for a in range(2) for b in range(2) for c in range(2)
    ])
    rho = np.diag(probabilities)
    np.testing.assert_allclose(stage1.petz_recovery(rho), rho, atol=1e-15)
    row = stage1.analyze_state(rho)
    assert row["cmi_bits"] == pytest.approx(0.0, abs=1e-14)
    assert row["petz_fidelity"] == pytest.approx(1.0, abs=1e-14)
    assert row["fawzi_renner_fidelity_lower_bound"] == pytest.approx(1.0, abs=1e-14)


@pytest.mark.parametrize("mass", (2.**-35, 2.**-40, 1e-100, 1e-310))
def test_petz_retains_every_resolved_markov_sector(mass):
    rho = np.diag([1-mass, 0, 0, 0, 0, 0, 0, mass])
    recovered = stage1.petz_recovery(rho)
    assert recovered[7, 7].real == pytest.approx(mass, rel=2e-13, abs=0)
    np.testing.assert_allclose(recovered, rho, atol=0, rtol=2e-13)


@pytest.mark.parametrize("bad", (
    np.zeros((8, 8)), np.eye(8), np.diag([1.1, -.1, 0, 0, 0, 0, 0, 0]),
    np.ma.array(np.eye(8)/8, mask=np.eye(8, dtype=bool)),
))
@pytest.mark.parametrize("operation", ("analyze_state", "petz_recovery", "state_fidelity"))
def test_invalid_states_cannot_become_recovery_evidence(bad, operation):
    with pytest.raises(ValueError):
        if operation == "state_fidelity":
            stage1.state_fidelity(bad, np.eye(8)/8)
        else:
            getattr(stage1, operation)(bad)


def test_entropy_retains_positive_spectral_mass():
    mass = 1e-100
    assert stage1.von_neumann_entropy(np.diag([1., mass])) == pytest.approx(
        -mass*np.log2(mass), rel=1e-14, abs=0)


@pytest.mark.parametrize("counts", ({}, {"XYZ": {"000": 100}},
                                      {b: {} for b in stage1.measurement_bases(3)}))
def test_incomplete_tomography_does_not_invent_a_state(counts):
    with pytest.raises(ValueError):
        stage1.reconstruct_density_matrix(counts, 3)


@pytest.mark.parametrize("counts", ({}, {"0": -1, "1": 2}, {"0": True},
                                      {"0": .5, "1": .5}, {"2": 1}, {"00": 1}))
def test_invalid_counts_are_not_zero_expectations(counts):
    with pytest.raises(ValueError):
        stage1.expectation_from_counts(counts, "Z")


def test_tomography_uses_shot_weights_for_repeated_pauli_estimates():
    counts = {b: {"00": 1} for b in stage1.measurement_bases(2)}
    counts["XY"] = {"01": 100}  # q0 is the rightmost count bit.
    _, expectations = stage1.reconstruct_density_matrix(counts, 2)
    assert expectations["XI"] == pytest.approx(-98/102, abs=1e-15)
