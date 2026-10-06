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


@pytest.mark.parametrize("occupied", ([0], [0,7], [0,5]))
def test_complete_counts_reconstruct_analytic_states_in_q0_order(occupied):
    # Independent Born probabilities of the product measurement projectors.
    psi = np.zeros(8)
    psi[occupied] = 1/np.sqrt(len(occupied))
    rho = np.outer(psi,psi)
    counts = {}
    for basis in stage1.measurement_bases(3):
        row = {}
        for outcome in range(8):
            bits = f"{outcome:03b}"
            projector = np.ones((1,1))
            for bit, pauli in zip(bits,basis):
                projector = np.kron(projector, (np.eye(2)+(-1)**int(bit)*stage1.PAULI_MATRICES[pauli])/2)
            shots = 256*np.trace(rho @ projector).real
            assert abs(shots-round(shots)) < 1e-12
            row[bits[::-1]] = round(shots)
        counts[basis] = row
    recovered, expectations, diagnostics = stage1.reconstruct_density_matrix(
        counts,3,return_diagnostics=True)
    np.testing.assert_allclose(recovered,rho,atol=1e-15)
    for label,value in expectations.items():
        assert value == pytest.approx(stage1.pauli_expectation(rho,label),abs=5e-16)
    assert diagnostics["complete_settings"] == 27
    assert set(diagnostics["shots_by_basis"].values()) == {256}
    assert diagnostics["state_correction_frobenius"] < 2e-15
    assert diagnostics["statistical_error_certified"] is False


def test_tomography_discloses_the_nonphysical_linear_estimate_correction():
    # Simultaneous deterministic +X,+Y,+Z is an impossible qubit state.
    rho, _, report = stage1.reconstruct_density_matrix(
        {b:{"0":100} for b in "XYZ"},1,return_diagnostics=True)
    negative = (np.sqrt(3)-1)/2
    assert report["raw_min_eigenvalue"] == pytest.approx(-negative,abs=2e-16)
    assert report["negative_spectral_mass"] == pytest.approx(negative,abs=2e-16)
    assert report["state_correction_frobenius"] == pytest.approx(np.sqrt(2)*negative,abs=4e-16)
    assert np.linalg.eigvalsh(rho)[0] >= -1e-16


def test_one_missing_setting_can_hide_a_full_bit_of_conditional_information():
    uniform = np.eye(8)/8
    parity = np.diag([.25 if i.bit_count()%2 == 0 else 0 for i in range(8)])
    assert stage1.conditional_mutual_information(uniform) == pytest.approx(0,abs=1e-15)
    assert stage1.conditional_mutual_information(parity) == pytest.approx(1,abs=1e-15)
    incomplete = {}
    for basis in stage1.measurement_bases(3):
        if basis == "ZZZ":
            continue
        for i in range(8):
            projector = np.ones((1,1))
            for bit,pauli in zip(f"{i:03b}",basis):
                projector = np.kron(projector,(np.eye(2)+(-1)**int(bit)*stage1.PAULI_MATRICES[pauli])/2)
            assert np.trace(uniform @ projector) == np.trace(parity @ projector) == .125
        incomplete[basis] = {f"{i:03b}":1 for i in range(8)}
    with pytest.raises(ValueError,match="complete"):
        stage1.reconstruct_density_matrix(incomplete,3)


@pytest.mark.parametrize("bad", (-1.,float("nan"),float("inf"),True,[0.],1j))
def test_recovery_bound_requires_valid_cmi(bad):
    with pytest.raises(ValueError):
        stage1.fawzi_renner_fidelity_lower_bound(bad)


@pytest.mark.parametrize("bad", ({"0":0}, {"":1}, {0:1}, {"0":np.int64(1)}))
def test_counts_have_one_unambiguous_json_integer_contract(bad):
    with pytest.raises(ValueError):
        stage1.expectation_from_counts(bad,"Z")


def test_metrics_never_call_tomographic_state_repair(monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("analysis attempted to manufacture a replacement state")
    monkeypatch.setattr(stage1,"project_to_physical_density_matrix",forbidden)
    result = stage1.analyze_state(np.eye(8)/8)
    assert result["petz_fidelity"] == pytest.approx(1.,abs=2e-15)


def test_invalid_evidence_is_rejected_under_optimized_python(tmp_path):
    import subprocess
    import sys
    script = f"""
import importlib.util
import numpy as np
s=importlib.util.spec_from_file_location('stage1', {str(MODULE_PATH)!r})
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)
for call in (lambda: m.analyze_state(np.zeros((8,8))),
             lambda: m.state_fidelity(np.eye(8),np.eye(8)/8),
             lambda: m.reconstruct_density_matrix({{}},3),
             lambda: m.expectation_from_counts({{'0':True}},'Z')):
    try:
        call()
    except ValueError:
        continue
    raise SystemExit('invalid evidence accepted')
if abs(m.analyze_state(np.eye(8)/8)['petz_fidelity']-1)>1e-14:
    raise SystemExit('valid control failed')
"""
    result = subprocess.run([sys.executable,"-O","-c",script],cwd=tmp_path,
                            capture_output=True,text=True,timeout=30)
    assert result.returncode == 0,result.stderr


def test_full_rank_counterexample_keeps_petz_below_the_optimal_recovery_bound():
    # Zhang, Table I and Eq. (20), epsilon=10^-4 (2026-09-12):
    # https://yuxuanzhang1995.github.io/assets/pdf/agentic/petz-cmi.pdf
    # These published coefficients are data, not an OPH fit. Reorder the
    # author's (A,B,C) to our (B,C,A), so our middle slot is the conditioner.
    import mpmath as mp
    rows = (
        ("-.6104357345", "0", "-.0465530650", "0"),
        (".0063927953", "-.0044244357", "-.0377544806", "-.0978816076"),
        (".0026075606", "-.0150331220", ".8582239838", "-.4826969158"),
        ("-.0001180574", "-.0002823919", "-.0132291583", ".0104067043"),
        ("-.0118114830", ".0054943442", ".0018681142", ".0014701336"),
        ("-.0003861017", "-.0098540922", "-.0012199090", ".0109247234"),
        (".7842448379", "-.1079189498", "-.0488402257", "-.0085549482"),
        ("-.0061499192", ".0091457466", "-.0582671635", "-.1049688866"),
    )
    with mp.workdps(65):
        v = mp.matrix([[mp.mpc(row[0],row[1]),mp.mpc(row[2],row[3])] for row in rows])
        original = v*mp.diag([mp.mpf(".0383477362"),mp.mpf(".9616522638")])*v.H
        original /= sum(original[i,i] for i in range(8))
        original = mp.mpf(".9999")*original+mp.mpf(".0001")*mp.eye(8)/8
        permutation = [4*c+2*a+b for a in range(2) for b in range(2) for c in range(2)]
        rho = mp.matrix([[original[i,j] for j in permutation] for i in permutation])
        # Independent index contractions, with no production reduction/root code.
        ab = mp.matrix([[sum(rho[2*i+c,2*j+c] for c in range(2)) for j in range(4)] for i in range(4)])
        bc = mp.matrix([[sum(rho[4*a+i,4*a+j] for a in range(2)) for j in range(4)] for i in range(4)])
        b = mp.matrix([[sum(bc[2*i+c,2*j+c] for c in range(2)) for j in range(2)] for i in range(2)])
        inverse = mp.sqrtm(b)**-1
        whiten = mp.matrix([[inverse[(i//2)%2,(j//2)%2] if i//4 == j//4 and i%2 == j%2 else 0
                            for j in range(8)] for i in range(8)])
        lift_ab = mp.matrix([[ab[i//2,j//2] if i%2 == j%2 else 0 for j in range(8)] for i in range(8)])
        sqrt_bc = mp.sqrtm(bc)
        outer = mp.matrix([[sqrt_bc[i%4,j%4] if i//4 == j//4 else 0 for j in range(8)] for i in range(8)])
        recovered = outer*whiten*lift_ab*whiten*outer
        entropy = lambda a: -mp.fsum(x*mp.log(x,2) for x in mp.eighe(a,eigvals_only=True))
        cmi = entropy(ab)+entropy(bc)-entropy(b)-entropy(rho)
        root = mp.sqrtm(recovered)
        fidelity = mp.fsum(mp.sqrt(x) for x in mp.eighe(root*rho*root,eigvals_only=True))**2
        gap = cmi+mp.log(fidelity,2)
        assert mp.mpf("-.001803") < gap < mp.mpf("-.001802")
        supplied = np.array(rho.tolist(),complex)
        expected_recovery = np.array(recovered.tolist(),complex)
        expected_cmi, expected_fidelity = float(cmi), float(fidelity)
    np.testing.assert_allclose(stage1.petz_recovery(supplied),expected_recovery,atol=1e-14)
    row = stage1.analyze_state(supplied)
    assert row["cmi_bits"] == pytest.approx(expected_cmi,abs=2e-14)
    assert row["petz_fidelity"] == pytest.approx(expected_fidelity,abs=2e-13)
    assert row["petz_fidelity"] < row["fawzi_renner_fidelity_lower_bound"]-.001
    assert row["fawzi_renner_bound_scope"] == "optimal_recovery_over_B_to_BC_channels"


@pytest.mark.parametrize("seed", range(4))
def test_noisy_tomography_returns_an_actually_positive_matrix(seed):
    """Complete finite-shot counts must remain analyzable with exact PSD checks."""
    import json
    import sympy as sp

    rng = np.random.default_rng(seed)
    psi = rng.normal(size=8)+1j*rng.normal(size=8)
    psi /= np.linalg.norm(psi)
    counts = {}
    for basis in stage1.measurement_bases(3):
        probabilities = []
        for bits in range(8):
            projector = np.ones((1, 1))
            for bit, axis in zip(f"{bits:03b}", basis):
                projector = np.kron(projector, (
                    np.eye(2)+(-1)**int(bit)*stage1.PAULI_MATRICES[axis])/2)
            probabilities.append(float(np.vdot(psi, projector@psi).real))
        samples = rng.multinomial(8192, np.array(probabilities)/sum(probabilities))
        counts[basis] = {f"{i:03b}"[::-1]: int(v) for i, v in enumerate(samples)}
    rho, _, diagnostic = stage1.reconstruct_density_matrix(
        counts, 3, return_diagnostics=True)
    # Exact characteristic coefficients independently certify the matrix
    # actually returned, not the intended eigenspectrum before rounding.
    exact = sp.Matrix([[sp.Rational(float(z.real))+sp.I*sp.Rational(float(z.imag))
                        for z in row] for row in rho])
    assert exact == exact.H
    assert all((-1)**k*x >= 0 for k, x in enumerate(exact.charpoly().all_coeffs()))
    assert 0 <= diagnostic["gram_rounding_trace_bound"] < 2e-15
    assert diagnostic["trace_normalization_tolerance"] == 1e-12
    metrics = stage1.analyze_state(rho)
    assert np.isfinite(metrics["cmi_bits"]) and metrics["cmi_bits"] > 0
    replay, _, replay_diagnostic = stage1.reconstruct_density_matrix(
        json.loads(json.dumps(counts)), 3, return_diagnostics=True)
    np.testing.assert_array_equal(replay, rho)
    assert replay_diagnostic == diagnostic


def test_exact_tomographic_states_do_not_receive_a_positive_floor():
    for state in [np.diag([1., 0.]), np.array([[.5, .5j], [-.5j, .5]])]:
        actual, bound = stage1.project_to_physical_density_matrix(state, return_rounding_bound=True)
        assert bound is None
        assert np.linalg.matrix_rank(actual) == 1
        np.testing.assert_allclose(actual, state, atol=2e-16, rtol=0)


def _stub_circuit_statevector(monkeypatch, amplitudes):
    """Supply the SDK's numerical array contract without a cloud dependency."""
    import sys
    from types import ModuleType, SimpleNamespace
    sdk = ModuleType("qiskit")
    info = ModuleType("qiskit.quantum_info")
    info.Statevector = SimpleNamespace(from_instruction=lambda circuit:
                                      SimpleNamespace(data=circuit.amplitudes))
    info.DensityMatrix = lambda state: SimpleNamespace(
        data=np.outer(state.data, state.data.conj()))
    sdk.quantum_info = info
    monkeypatch.setitem(sys.modules, "qiskit", sdk)
    monkeypatch.setitem(sys.modules, "qiskit.quantum_info", info)
    return SimpleNamespace(num_qubits=3, amplitudes=amplitudes)


@pytest.mark.parametrize("case", [0, 1, 2, 3, "theta_.6", "theta_1"])
def test_circuit_state_constructor_returns_exact_psd_entries(monkeypatch, case):
    import sympy as sp
    if isinstance(case, int):
        rng = np.random.default_rng(case)
        psi = rng.normal(size=8)+1j*rng.normal(size=8)
        psi /= np.linalg.norm(psi)
    else:
        theta = .6 if case == "theta_.6" else 1.
        psi = np.zeros(8, complex)
        psi[[0, 3, 7]] = np.array([1., np.cos(theta/2), np.sin(theta/2)])/np.sqrt(2)
    circuit = _stub_circuit_statevector(monkeypatch, psi)
    rho = stage1.circuit_density_q0_order(circuit)
    exact = sp.Matrix([[sp.Rational(float(z.real))+sp.I*sp.Rational(float(z.imag))
                        for z in row] for row in rho])
    assert exact == exact.H
    assert all((-1)**k*x >= 0 for k, x in enumerate(exact.charpoly().all_coeffs()))
    # Independent bit reversal: Qiskit little-endian becomes the q0-first order.
    q0 = psi[[0, 4, 2, 6, 1, 5, 3, 7]]
    np.testing.assert_allclose(rho, np.outer(q0, q0.conj()), atol=1e-15, rtol=0)
    metrics = stage1.analyze_state(rho)
    assert np.isfinite(metrics["cmi_bits"])
