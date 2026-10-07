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
    # LAPACK/platform rounding may already produce an exactly PSD candidate.
    # In that case the documented bound is None; the independent PSD check
    # above still applies. Do not require a particular fallback decision.
    bound = diagnostic["gram_rounding_trace_bound"]
    if bound is not None:
        assert 0 <= bound < 2e-15
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
    actual, diagnostic = stage1.circuit_density_q0_order(circuit, return_diagnostics=True)
    np.testing.assert_array_equal(actual, rho)
    from sympy.polys.matrices import DomainMatrix
    v = sp.Matrix([sp.Rational(float(z.real))+sp.I*sp.Rational(float(z.imag)) for z in q0])
    exact_v = DomainMatrix.from_Matrix(v).convert_to(sp.QQ_I).to_dense()
    adjoint = DomainMatrix.from_Matrix(v.H).convert_to(sp.QQ_I).to_dense()
    error = exact-exact_v.matmul(adjoint).to_Matrix()
    assert all((-1)**k*x >= 0 for k, x in enumerate(error.charpoly().all_coeffs()))
    assert 0 <= sp.trace(error) <= sp.Rational(diagnostic["gram_rounding_trace_bound"])
    assert diagnostic["gram_rounding_trace_bound"] < 2e-15
    assert diagnostic["rounding_bound_scope"] == "supplied_numerical_statevector"
    assert diagnostic["circuit_simulation_error_certified"] is False


@pytest.mark.parametrize("psi", [np.zeros(8), np.ones(8), np.ones(4),
                               np.full(8, np.nan), np.full(8, np.inf)])
def test_invalid_circuit_amplitudes_are_not_normalized_or_replaced(monkeypatch, psi):
    circuit = _stub_circuit_statevector(monkeypatch, psi)
    with pytest.raises(ValueError):
        stage1.circuit_density_q0_order(circuit)


def test_exact_product_circuit_gets_no_rounding_floor(monkeypatch):
    psi = np.zeros(8, complex)
    psi[1] = 1j
    rho, diagnostic = stage1.circuit_density_q0_order(
        _stub_circuit_statevector(monkeypatch, psi), return_diagnostics=True)
    expected = np.zeros((8, 8))
    expected[4, 4] = 1
    np.testing.assert_array_equal(rho, expected)
    assert diagnostic["gram_rounding_trace_bound"] == 0


def test_full_report_retains_circuit_rounding_and_count_replay(monkeypatch, tmp_path):
    import json
    import sys
    from types import ModuleType, SimpleNamespace
    _stub_circuit_statevector(monkeypatch, np.zeros(8))
    def circuit(name, amplitudes):
        return SimpleNamespace(name=name, num_qubits=3, amplitudes=amplitudes)
    def structured(theta):
        psi = np.zeros(8, complex)
        psi[[0, 3, 7]] = np.array([1., np.cos(theta/2), np.sin(theta/2)])/np.sqrt(2)
        return circuit(f"structured_theta_{theta:.2f}", psi)
    def random(seed, depth):
        rng = np.random.default_rng(seed)
        psi = rng.normal(size=8)+1j*rng.normal(size=8)
        return circuit(f"random_seed_{seed}", psi/np.linalg.norm(psi))
    monkeypatch.setattr(stage1, "build_structured_family", structured)
    monkeypatch.setattr(stage1, "build_random_control", random)
    monkeypatch.setattr(stage1, "build_ghz", lambda: circuit(
        "ghz_control", np.array([1., 0, 0, 0, 0, 0, 0, 1.])/np.sqrt(2)))
    monkeypatch.setattr(stage1, "add_measurement_basis", lambda c, b:
                        SimpleNamespace(name=f"{c.name}__{b}"))
    monkeypatch.setattr(stage1, "parse_args", lambda: SimpleNamespace(
        local_testing=True, mode="local", outdir=tmp_path, random_depth=3,
        random_seeds=[0, 1], shots=256, transpile_seed=7,
        credentials_file=tmp_path/'unused', backend=None))
    counts = {f"{i:03b}": 32 for i in range(8)}
    monkeypatch.setattr(stage1, "run_sampler", lambda **kwargs: (
        {"counts_by_name": {c.name: counts for c in kwargs["circuits"]},
         "run_metadata": {"test_counts": True}}, None))
    common = ModuleType("ibm_runtime_common")
    common.ensure_dir = lambda path: path
    common.write_json = lambda path, data: path.write_text(json.dumps(data))
    monkeypatch.setitem(sys.modules, "ibm_runtime_common", common)
    assert stage1.main() == 0
    report = json.loads((tmp_path/'summary.json').read_text())
    assert len(report["exact_analysis"]) == 5
    for metrics in report["exact_analysis"].values():
        assert metrics["circuit_state"]["rounding_bound_scope"] == "supplied_numerical_statevector"
        assert metrics["circuit_state"]["circuit_simulation_error_certified"] is False
    selection = report["catalog"]["random_control_selection"]
    assert all("circuit_state" in c for c in selection["candidate_summary"])
    assert selection["circuit_state"] == next(c["circuit_state"] for c in
        selection["candidate_summary"] if c["seed"] == selection["seed"])
    for name, settings in report["tomography_counts_by_state"].items():
        rho, _ = stage1.reconstruct_density_matrix(settings, 3)
        assert stage1.analyze_state(rho)["cmi_bits"] == report["reconstructed_analysis"][name]["cmi_bits"]


def _review_counts():
    # Review #1052: independent A,C are uniform; B has Bloch vector (.6,.8,0).
    # All 27 settings have 40 shots. No counts were filtered for recovery.
    return {b: {f"{k:03b}": {"X": (8, 2), "Y": (9, 1), "Z": (5, 5)}[b[1]][(k >> 1) & 1]
                for k in range(8)} for b in stage1.measurement_bases(3)}


def _assert_unresolved_report(row):
    assert row["petz_status"] == "unresolved_support"
    assert row["petz_unavailable_reason"] == "reference support is numerically unresolved"
    for key in ("petz_fidelity", "petz_trace_distance", "petz_observable_mismatch"):
        assert row[key] is None
    assert row["cmi_bits"] == pytest.approx(0, abs=2e-13)
    assert row["fawzi_renner_fidelity_lower_bound"] == pytest.approx(1, abs=2e-13)
    assert row["fawzi_renner_bound_scope"] == "optimal_recovery_over_B_to_BC_channels"


def test_review_complete_counts_preserve_cmi_with_unresolved_recovery():
    import sympy as sp
    rho, _ = stage1.reconstruct_density_matrix(_review_counts(), 3)
    # Certify the actual returned entries independently of the producer's Schur test.
    exact = sp.Matrix([[sp.Rational(float(z.real))+sp.I*sp.Rational(float(z.imag))
                        for z in row] for row in rho])
    assert all((-1)**k*x >= 0 for k, x in enumerate(exact.charpoly().all_coeffs()))
    b = np.array([[.5, .3-.4j], [.3+.4j, .5]])
    np.testing.assert_allclose(rho, np.kron(np.kron(np.eye(2)/2, b), np.eye(2)/2),
                               atol=2e-16, rtol=0)
    assert stage1.conditional_mutual_information(rho) == pytest.approx(0, abs=2e-13)
    with pytest.raises(ValueError, match="support is numerically unresolved"):
        stage1.petz_recovery(rho)
    _assert_unresolved_report(stage1.analyze_state(rho))


def test_review_seed_one_reference_retains_valid_information(monkeypatch):
    # Unmodified output of Statevector.from_instruction(random_circuit(
    # 3, depth=3, max_operands=2, measure=False, seed=1)), Qiskit 2.5.2.
    psi = np.array([
        .012583908990800391-.08544448948589424j,
        -.06695570762093886+.009860958104918296j,
        .08544448948589424+.012583908990800391j,
        -.009860958104918296-.06695570762093886j,
        -.54396834895709-.08011339570683529j,
        .06277817483306489+.4262625471456357j,
        .08011339570683529-.54396834895709j,
        -.4262625471456357+.06277817483306489j,
    ])
    circuit = _stub_circuit_statevector(monkeypatch, psi)
    rho = stage1.circuit_density_q0_order(circuit)
    with pytest.raises(ValueError, match="support is numerically unresolved"):
        stage1.petz_recovery(rho)
    _assert_unresolved_report(stage1.analyze_state(rho))


def _review_main_fixture(monkeypatch, tmp_path, *, sdk=False, unresolved_counts=True):
    import json
    import sys
    from types import ModuleType, SimpleNamespace
    if sdk:
        pytest.importorskip("qiskit")
    else:
        psi = np.zeros(8)
        psi[0] = 1
        _stub_circuit_statevector(monkeypatch, psi)
        def circuit(name):
            return SimpleNamespace(name=name, num_qubits=3, amplitudes=psi)
        monkeypatch.setattr(stage1, "build_structured_family", lambda theta:
                            circuit(f"structured_theta_{theta:.2f}"))
        monkeypatch.setattr(stage1, "build_ghz", lambda: circuit("ghz_control"))
        monkeypatch.setattr(stage1, "build_random_control", lambda seed, depth:
                            circuit(f"random_seed_{seed}"))
        monkeypatch.setattr(stage1, "add_measurement_basis", lambda c, b:
                            SimpleNamespace(name=f"{c.name}__{b}"))
    # Use actual CLI parsing, including the maintainer's --random-seeds 1.
    monkeypatch.setattr(sys, "argv", [str(MODULE_PATH), "--local-testing", "--shots", "40",
                                      "--random-seeds", "1", "--outdir", str(tmp_path)])
    settings = _review_counts() if unresolved_counts else {
        b: {f"{k:03b}": 5 for k in range(8)} for b in stage1.measurement_bases(3)}
    acquired = {}
    def sampler(**kwargs):
        acquired.update(counts_by_name={c.name: settings[c.name.rsplit("__", 1)[1]]
                                       for c in kwargs["circuits"]},
                        run_metadata={"mode": "local", "fixture": True, "shots": 40})
        return acquired, "local_fixture"
    monkeypatch.setattr(stage1, "run_sampler", sampler)
    common = ModuleType("ibm_runtime_common")
    common.ensure_dir = lambda path: path
    common.write_json = lambda path, data: path.write_text(json.dumps(data, allow_nan=False))
    monkeypatch.setitem(sys.modules, "ibm_runtime_common", common)
    return acquired


@pytest.mark.parametrize("sdk,unresolved_counts", [(False, True), (True, True), (True, False)])
def test_review_main_retains_evidence_and_partial_results(monkeypatch, tmp_path, sdk, unresolved_counts):
    import json
    acquired = _review_main_fixture(monkeypatch, tmp_path, sdk=sdk,
                                   unresolved_counts=unresolved_counts)
    assert stage1.main() == 0
    raw = json.loads((tmp_path/"acquired_counts.json").read_text())
    report = json.loads((tmp_path/"summary.json").read_text())
    assert raw["counts_by_name"] == acquired["counts_by_name"]
    assert raw["run_metadata"] == acquired["run_metadata"]
    assert raw["backend"] == "local_fixture"
    assert raw["random_seeds"] == [1]
    assert len(raw["counts_by_name"]) == 135
    assert raw["catalog"] == report["catalog"]
    assert raw["timestamp_utc"] == report["timestamp_utc"]
    for name, mapping in raw["measured_index"].items():
        settings = {b: raw["counts_by_name"][c] for b, c in mapping.items()}
        assert settings == report["tomography_counts_by_state"][name]
        rho, _, diagnostic = stage1.reconstruct_density_matrix(settings, 3, return_diagnostics=True)
        row = report["reconstructed_analysis"][name]
        assert row["tomography"] == diagnostic
        assert all(row[k] == v for k, v in stage1.analyze_state(rho).items())
        if unresolved_counts:
            _assert_unresolved_report(row)
        else:
            assert row["petz_status"] == "available"
            assert row["petz_unavailable_reason"] is None
            assert row["petz_fidelity"] == pytest.approx(1, abs=2e-14)
    if sdk:
        _assert_unresolved_report(report["exact_analysis"]["random_seed_1"])
        assert report["exact_analysis"]["random_seed_1"]["circuit_state"][
            "circuit_simulation_error_certified"] is False
    checks = report["fingerprint_checks"]
    assert checks["recovery_improves_as_cmi_drops"] is (None if unresolved_counts else True)
    assert checks["structured_theta_0.00_lt_random_control"] is False
    assert checks["structured_theta_0.00_lt_ghz"] is False
    assert json.loads((tmp_path/"summary_pretty.txt").read_text()) == report


@pytest.mark.parametrize("failing_operation", ["analyze_state", "reconstruct_density_matrix"])
def test_review_counts_are_saved_before_analysis_errors(monkeypatch, tmp_path, failing_operation):
    import json
    acquired = _review_main_fixture(monkeypatch, tmp_path)
    def fail(*args, **kwargs):
        raw = json.loads((tmp_path/"acquired_counts.json").read_text())
        assert raw["counts_by_name"] == acquired["counts_by_name"]
        assert raw["run_metadata"] == acquired["run_metadata"]
        assert len(raw["measured_index"]) == 5
        raise ValueError("unexpected analysis failure")
    monkeypatch.setattr(stage1, failing_operation, fail)
    with pytest.raises(ValueError, match="unexpected analysis failure"):
        stage1.main()
    assert not (tmp_path/"summary.json").exists()


def test_review_does_not_classify_unrelated_errors_by_message(monkeypatch):
    def fail(*args):
        raise ValueError("reference support is numerically unresolved")
    monkeypatch.setattr(stage1, "petz_recovery", fail)
    with pytest.raises(ValueError, match="support is numerically unresolved"):
        stage1.analyze_state(np.eye(8)/8)


@pytest.mark.parametrize("unavailable", ["structured_theta_0.00", "structured_theta_0.60",
                                       "structured_theta_1.00", "random_seed_1"])
def test_review_only_dependent_comparisons_become_unavailable(monkeypatch, tmp_path, unavailable):
    import json
    _review_main_fixture(monkeypatch, tmp_path, unresolved_counts=False)
    sampler = stage1.run_sampler
    def mixed_sampler(**kwargs):
        result, backend = sampler(**kwargs)
        result["counts_by_name"].update({f"{unavailable}__{b}": counts
                                         for b, counts in _review_counts().items()})
        return result, backend
    monkeypatch.setattr(stage1, "run_sampler", mixed_sampler)
    assert stage1.main() == 0
    report = json.loads((tmp_path/"summary.json").read_text())
    for name, row in report["reconstructed_analysis"].items():
        if name == unavailable:
            _assert_unresolved_report(row)
        else:
            assert row["petz_status"] == "available"
            assert row["petz_unavailable_reason"] is None
            assert row["petz_fidelity"] == pytest.approx(1, abs=2e-14)
    expected = None if unavailable.startswith("structured") else True
    assert report["fingerprint_checks"]["recovery_improves_as_cmi_drops"] is expected


def test_review_incomplete_acquisition_is_saved_but_never_reported_as_success(monkeypatch, tmp_path):
    import json
    acquired = _review_main_fixture(monkeypatch, tmp_path)
    sampler = stage1.run_sampler
    missing = "structured_theta_0.00__XYZ"
    def incomplete_sampler(**kwargs):
        result, backend = sampler(**kwargs)
        del result["counts_by_name"][missing]
        return result, backend
    monkeypatch.setattr(stage1, "run_sampler", incomplete_sampler)
    with pytest.raises(KeyError, match=missing):
        stage1.main()
    raw = json.loads((tmp_path/"acquired_counts.json").read_text())
    assert raw["counts_by_name"] == acquired["counts_by_name"]
    assert len(raw["counts_by_name"]) == 134
    assert raw["measured_index"]["structured_theta_0.00"]["XYZ"] == missing
    assert not (tmp_path/"summary.json").exists()


def test_review_support_exception_is_only_handled_for_recovery(monkeypatch):
    from quantum_information.recovery import UnresolvedPetzSupport
    def fail(*args):
        raise UnresolvedPetzSupport("injected information failure")
    monkeypatch.setattr(stage1, "conditional_mutual_information", fail)
    with pytest.raises(UnresolvedPetzSupport, match="injected information failure"):
        stage1.analyze_state(np.eye(8)/8)


@pytest.mark.parametrize("prior_file", ["acquired_counts.json", "summary.json", "summary_pretty.txt"])
def test_review_existing_evidence_is_preserved_before_any_new_acquisition(monkeypatch, tmp_path, prior_file):
    acquired = _review_main_fixture(monkeypatch, tmp_path)
    previous = b"prior evidence must not be overwritten, even if incomplete"
    (tmp_path/prior_file).write_bytes(previous)
    with pytest.raises(FileExistsError, match="fresh --outdir"):
        stage1.main()
    assert acquired == {}
    assert (tmp_path/prior_file).read_bytes() == previous
    assert {p.name for p in tmp_path.iterdir()} == {prior_file}


def test_review_rerun_cannot_pair_new_counts_with_an_old_success(monkeypatch, tmp_path):
    _review_main_fixture(monkeypatch, tmp_path, unresolved_counts=False)
    assert stage1.main() == 0
    previous = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    acquired = _review_main_fixture(monkeypatch, tmp_path, unresolved_counts=True)
    original = stage1.run_sampler
    def incomplete_sampler(**kwargs):
        result, backend = original(**kwargs)
        del result["counts_by_name"]["structured_theta_0.00__XYZ"]
        return result, backend
    monkeypatch.setattr(stage1, "run_sampler", incomplete_sampler)
    with pytest.raises(FileExistsError, match="fresh --outdir"):
        stage1.main()
    assert acquired == {}
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == previous
