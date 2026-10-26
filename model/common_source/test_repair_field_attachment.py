"""Independent native/CAR/source/record controls for finite field attachment."""
from copy import deepcopy
from fractions import Fraction as F
import json
import math
from unittest.mock import patch

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import repair_field_attachment as f


@pytest.fixture(scope="module")
def execution():
    return f.execute()


def _pauli_annihilator(site):
    identity, z = s.eye(2), s.diag(1, -1)
    ladder = s.Matrix([[0, 1], [0, 0]])
    factors = [ladder if j == site else z if j < site else identity for j in reversed(range(3))]
    return s.kronecker_product(*factors)


def test_independent_pauli_tensor_CAR_and_full_interacting_operator():
    a = [_pauli_annihilator(j) for j in range(3)]
    for j in range(3):
        assert a[j] == f.annihilator(j)
    p = f.logical_parameters()
    for sector in range(3):
        h = 3*p["onsite_baselines"][sector]*s.eye(8)
        h += p["positive_occupation_gaps"][sector]*sum((x.T*x for x in a), s.zeros(8))
        for j in range(2):
            n, m = a[j].T*a[j], a[j+1].T*a[j+1]
            h += p["exchange_strength"]*(n+m-2*n*m-a[j].T*a[j+1]-a[j+1].T*a[j])
        assert h == f.chain_hamiltonian(sector)
    assert f.global_car_certificate()["all_basis_anticommutators_checked"] == 82944


def test_three_chains_are_distinct_carriers_with_positive_inherited_thresholds():
    result = f.exact_field_checks()
    assert result["total_independent_CAR_modes"] == 9
    assert result["inter_chain_mixing"] == "0"
    assert f.logical_parameters()["positive_occupation_gaps"] == (s.Rational(1, 256), s.Rational(1, 1024), s.Rational(1, 3072))
    for row in result["per_chain"]:
        assert row["minimum_energy_above_vacuum_by_charge"][0] == "0"
    assert not result["two_particle_is_free_sum_of_one_particle_poles"]
    # Fully occupied exchange energy is0, whereas a free one-body lift has
    # the nonzero trace4J. The density interaction cannot be dropped.
    h = f.chain_hamiltonian(0)
    baseline = 3*f.logical_parameters()["onsite_baselines"][0]
    assert h[7, 7]-baseline == 3*s.Rational(1, 256)


def test_nonadjacent_loop_has_fermion_parity_obstruction():
    a, b = _pauli_annihilator(0), _pauli_annihilator(2)
    n, m = a.T*a, b.T*b
    naive = n+m-2*n*m-a.T*b-b.T*a
    assert ((s.eye(8)-f.carrier_swap(0, 2))-naive)[6, 3] == -2
    assert not f.exact_field_checks()["relativistic_Weyl_Dirac_or_3D_local_CAR_dispersion"]


@pytest.mark.parametrize("ports", ((0, 0), (0, 7), (11, 3)))
def test_native_two_port_cross_copy_swap_is_scalar_not_hopping(ports):
    # Independent original configuration enumeration; changed per-copy
    # occupation leaves the native N=2 code and contributes zero.
    supports = f.code.SUPPORTS
    membership = {tuple(state): logical for logical, row in enumerate(supports) for state in row}
    compressed = s.zeros(36)
    for left, row in enumerate(supports):
        for right, other in enumerate(supports):
            for x in row:
                for y in other:
                    xx, yy = list(x), list(y)
                    xx[ports[0]], yy[ports[1]] = yy[ports[1]], xx[ports[0]]
                    i, j = membership.get(tuple(xx)), membership.get(tuple(yy))
                    if i is not None and j is not None:
                        compressed[6*i+j, 6*left+right] += s.Rational(1, 36)
    assert compressed == s.Rational(13, 18)*s.eye(36)


def test_full_native_orbit_average_fixes_padding_and_preserves_each_code_column():
    original, T, averaged = f._orbit_data()
    native_addresses = {f.circuit.encode(state) for state in f.native.states()}
    for address in range(512):
        current = address
        total = F()
        for _ in range(6):
            total += original[current]
            current = T[current]
        assert current == address and averaged[address] == total/6
        if address not in native_addresses:
            assert T[address] == address and original[address] == averaged[address] == 0
    for logical, row in enumerate(f.code.SUPPORTS):
        assert {averaged[f.circuit.encode(state)] for state in row} == {F(f._parameters()[0][logical])}


def test_native_source_controls_and_complete_charge_readout(execution):
    report, arrays = execution
    f.verify_retained_arrays(arrays)
    assert arrays["native_orbit_averaged_outputs"].shape == (512, 6)
    assert arrays["native_charge_pointer_outputs"].shape == (2, 512, 6)
    assert report["collective_native_source"]["retained_full_native_output_rows"] == 2**19
    assert report["collective_native_source"]["unknown_Fock_columns"] == 4
    assert report["native_charge_source"]["operator_error_diagnostic"] < f.SOURCE_GATE
    assert not report["operational_charge"]["electric_color_weak_charge_or_SM_generation_identified"]


def test_native_collective_swap_phase_is_ordinary_not_fermionic(execution):
    _, arrays = execution
    actual = arrays["collective_logical_projection"]
    J = float(f.logical_parameters()["exchange_strength"])
    expected = expm(-1j*J*(np.eye(4)-np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])))
    np.testing.assert_allclose(actual, expected, atol=f.SOURCE_GATE, rtol=0)
    assert actual[3, 3] == pytest.approx(1, abs=f.SOURCE_GATE)
    fermion = f.gaussian.local_gate("fswap").fock_unitary
    assert abs(fermion[3, 3]-actual[3, 3]) > 1.9


def test_little_endian_macro_factorization_matches_literal_direct_source():
    gate = f.fredkin_source()
    targets = (1, 4, 0)
    direct = f.compiler.execute(f.compiler.remap(gate["word"], targets), 5)
    lifted = f._little_apply(gate["matrix"], targets, 5, np.eye(32, dtype=complex))
    np.testing.assert_allclose(lifted, direct, atol=2e-13, rtol=0)


def test_source_field_all_sectors_and_continuity(execution):
    report, arrays = execution
    number = np.asarray(f.occupation(), complex)
    for sector, U in enumerate(arrays["field_source_unitaries"]):
        row = report["finite_field_source"][sector]
        bound = F(row["splitting_operator_bound_exact_upper"])
        assert row["dense_target_difference_diagnostic"] < float(bound)
        np.testing.assert_allclose(U@number, number@U, atol=1e-10, rtol=0)
        assert np.linalg.norm(U.conj().T@U-np.eye(8), 2) < 1e-10
        h = f.chain_hamiltonian(sector)
        for site in range(3):
            a = f.annihilator(site); n = a.T*a
            expected = (f.current(site-1) if site else s.zeros(8))-(f.current(site) if site < 2 else s.zeros(8))
            assert s.I*(h*n-n*h) == expected


def test_uniform_mode_and_vacuum_are_fixed_by_each_exchange_factor(execution):
    report, arrays = execution
    q = np.zeros(8, complex); q[[1, 2, 4]] = 1/math.sqrt(3)
    for sector, U in enumerate(arrays["uniform_gap_source_unitaries"]):
        p = f.logical_parameters()
        expected = np.exp(-1j*float(f.GAP_TIME*(3*p["onsite_baselines"][sector]+p["positive_occupation_gaps"][sector])))*q
        np.testing.assert_allclose(U@q, expected, atol=3e-13, rtol=0)
    assert report["uniform_charged_gap_records"]["uniform_mode_splitting_error"] == "0"
    prior = report["uniform_charged_gap_records"]["non_alias_prior_from_native_architecture"]
    assert prior["positive_phase_prior"] == ["0", "3"]
    assert f.pi_interval()[0] > 3


def test_count_only_inverse_recovers_all_three_unfitted_native_gaps(execution):
    report, arrays = execution
    with patch.object(f, "logical_parameters", side_effect=AssertionError("native gap access")), \
         patch.object(f, "chain_hamiltonian", side_effect=AssertionError("Hamiltonian access")), \
         patch.object(f, "_parameters", side_effect=AssertionError("source data access")):
        boxes = f.decode_uniform_records(arrays["uniform_gap_complete_simulated_counts"])
    for box, gap in zip(boxes, (F(1, 256), F(1, 1024), F(1, 3072)), strict=True):
        assert 0 < box[0] <= gap <= box[1]
    records = report["uniform_charged_gap_records"]
    assert records["simulated_all_attempt_records"] == 12000000
    assert F(records["simultaneous_sampling_certificate"]["familywise_error_upper"]) < F(1, 100)
    assert not report["source_contract"]["full_native_field_numerical_composition_error_certified"]


def test_two_time_density_observers_share_one_backaction_history(execution):
    report, arrays = execution
    for sector, vectors in enumerate(arrays["complete_two_time_density_outputs"]):
        probabilities = np.array([np.vdot(v, v).real for v in vectors]).reshape(2, 2)
        row = report["continuing_density_observers"][sector]
        np.testing.assert_allclose(probabilities.sum(axis=1), row["observer_A_first_record_marginal"], atol=1e-13)
        np.testing.assert_allclose(probabilities.sum(axis=0), row["observer_B_second_record_marginal"], atol=1e-13)
        assert probabilities.sum() == pytest.approx(1, abs=1e-10)
        assert row["first_measurement_backaction_retained"]


@pytest.mark.parametrize("value", (s.S.Half, s.S.One, s.S.Zero, F(3, 2), 2))
def test_exact_rational_strength_accepts_sympy_singletons(value):
    assert f.chain_hamiltonian(strength=value).shape == (8, 8)


@pytest.mark.parametrize("value", (True, .5, float("nan"), -1))
def test_bad_original_strength_refuses(value):
    with pytest.raises(ValueError):
        f.chain_hamiltonian(strength=value)


def test_complete_checker_rejects_phase_zero_missing_arrays_and_budget_forgery(execution):
    report, arrays = execution
    for kind in ("phase", "zero", "missing"):
        forged = dict(arrays)
        if kind == "missing":
            forged.pop("native_charge_pointer_outputs")
        else:
            forged["native_charge_pointer_outputs"] = forged["native_charge_pointer_outputs"].copy()
            forged["native_charge_pointer_outputs"][0] *= -1 if kind == "phase" else 0
        with pytest.raises(ValueError):
            f.compare_execution(report, forged, report, arrays)
    forged = deepcopy(report)
    forged["source_contract"]["numerical_source_operator_admission"] = "1"
    with pytest.raises(ValueError):
        f.compare_execution(forged, arrays, report, arrays)


def test_independent_retained_operator_gate_catches_distributed_small_forgery(execution):
    _, arrays = execution
    forged = dict(arrays)
    forged["native_encoded_six_columns"] = arrays["native_encoded_six_columns"].copy()+1e-9
    with pytest.raises(ValueError, match="operator gate"):
        f.verify_retained_arrays(forged)
    forged = dict(arrays)
    forged["collective_fredkin_macro"] = arrays["collective_fredkin_macro"].copy()+1e-9
    with pytest.raises(ValueError, match="local macro"):
        f.verify_retained_arrays(forged)


def test_count_schema_zero_visibility_and_JSON_types(execution):
    report, arrays = execution
    counts = arrays["uniform_gap_complete_simulated_counts"]
    for bad in (counts.astype(float), np.ma.array(counts), counts.tolist(), counts[:2]):
        with pytest.raises(ValueError):
            f.decode_uniform_records(bad)
    with pytest.raises(ValueError):
        f.decode_uniform_records(np.full((3, 2, 2), f.TRIALS//2, np.int64))
    f.selection._compare(report, json.loads(json.dumps(report)))


def test_custody_rejects_duplicate_nonfinite_or_stale_JSON_before_execution(tmp_path):
    path = tmp_path/"receipt.json"
    path.write_text('{"source_bindings":{},"source_bindings":{}}')
    with pytest.raises(ValueError, match="duplicate JSON"):
        f.check(path)
    path.write_text('{"bad":NaN}')
    with pytest.raises(ValueError, match="nonfinite"):
        f.check(path)
    path.write_text('{"source_bindings":{}}')
    with patch.object(f, "source_bindings", return_value={"own": "sha"}), patch.object(f, "execute", side_effect=AssertionError("must not execute")):
        with pytest.raises(ValueError, match="source bindings"):
            f.check(path)
