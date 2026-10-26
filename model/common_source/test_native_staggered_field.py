"""Independent echo, parity, Clifford, source and continuum countercontrols."""
from copy import deepcopy
from fractions import Fraction as F
import json
import math
from unittest.mock import patch

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import native_staggered_field as n


@pytest.fixture(scope="module")
def execution():
    return n.execute()


def test_exact_native_echo_retains_double_occupation_and_absolute_phase():
    result = n.echo_certificate()
    assert result["density_term_cancelled_on_double_occupation"]
    np.testing.assert_allclose(result["unitary"], result["target"], atol=2e-13, rtol=0)
    assert result["unitary"][3, 3] == pytest.approx(1, abs=2e-13)
    # Replacing the inverse coherent leg by a forward leg changes the law.
    first = n.attachment.exchange_word(F(3, 10))
    z = n._phase_word(math.pi)
    wrong = n.compiler.execute(list(first)+z+list(first)+z, 2)
    assert np.linalg.norm(wrong-result["target"], 2) > .05


@pytest.mark.parametrize("eta", (-1, 1))
def test_nonadjacent_source_bond_has_full_CAR_parity_and_oriented_phase(eta):
    factors = n.bond_factors(0, 3, eta, F(2, 5))
    word = [event for local, sites in factors for event in n.compiler.remap(local, sites)]
    actual = n.compiler.execute(word, 4)
    h = np.zeros((4, 4), complex)
    h[0, 3] = -1j*eta*float(n.native_parameters()["exchange"])
    h[3, 0] = h[0, 3].conjugate()
    target = expm(-2j/5*n._fock_hamiltonian(h))
    np.testing.assert_allclose(actual, target, atol=3e-13, rtol=0)
    # Empty spectators hide the missing parity service; occupied ones expose it.
    bare = [(local, sites) for local, sites in factors if local != n._cz_word()]
    wrong = n.compiler.execute([event for local, sites in bare for event in n.compiler.remap(local, sites)], 4)
    one = [1 << j for j in range(4)]
    np.testing.assert_allclose(wrong[np.ix_(one, one)], target[np.ix_(one, one)], atol=3e-13, rtol=0)
    assert np.linalg.norm(wrong-target, 2) > .05


def test_derived_Clifford_lorentz_signature_spin_and_two_tastes():
    result = n.clifford_certificate()
    assert result["Dirac_tastes_per_native_sector"] == 2
    assert result["derived_continuum_gamma_signature"] == [1, -1, -1, -1]
    assert result["spin_Casimir"] == "3/4"
    I, X, Z = s.eye(2), s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
    expected = (s.kronecker_product(Z, I, I), s.kronecker_product(X, Z, I),
                s.kronecker_product(X, X, Z), s.kronecker_product(X, X, X))
    assert n.clifford_matrices() == expected


def test_actual_periodic_operator_matches_all_corner_blocks_and_both_signs():
    result = n.periodic_fourier_certificate()
    assert len(result["corner_blocks"]) == 8
    assert max(row["operator_error_diagnostic"] for row in result["corner_blocks"]) < 1e-13
    h = np.asarray(n.lattice_operator(), complex)
    assert np.count_nonzero(np.linalg.eigvalsh(h) < 0) == 32
    with pytest.raises(ValueError, match="2-site periodic"):
        n.lattice_operator((2, 2, 2))


def test_removing_staggered_flux_or_mass_sign_breaks_the_Dirac_relation():
    correct = np.asarray(n.lattice_operator((2, 2, 2), periodic=False), complex)
    mass = float(n.native_parameters()["mass"])
    with patch.object(n, "_eta", return_value=1):
        no_flux = np.asarray(n.lattice_operator((2, 2, 2), periodic=False), complex)
    with patch.object(n, "_epsilon", return_value=1):
        no_mass_sign = np.asarray(n.lattice_operator((2, 2, 2), periodic=False), complex)
    E2 = mass*mass+3*float(n.native_parameters()["exchange"])**2
    np.testing.assert_allclose(correct@correct, E2*np.eye(8), atol=1e-15, rtol=0)
    assert np.linalg.norm(no_flux@no_flux-E2*np.eye(8), 2) > .01
    assert np.linalg.norm(no_mass_sign@no_mass_sign-E2*np.eye(8), 2) > 1e-4


def test_literal_small_source_tape_equals_factor_execution_on_all256_inputs():
    cube = n.source_cube(steps=1)
    word = [event for block, sites in cube["factor_order"]
            for event in n.compiler.remap(cube["catalogue"][block]["word"], sites)]
    literal = n.compiler.execute(word, 8)
    assert len(word) == cube["resources"]["native_pulses"]
    np.testing.assert_allclose(literal, cube["unitary"], atol=2e-12, rtol=0)


def test_source_cube_retains_full_Fock_operator_charge_and_splitting(execution):
    report, arrays = execution
    U = arrays["source_cube_unitary"]
    h = arrays["source_cube_one_body_H"]
    target = expm(-1j*n._fock_hamiltonian(h))
    bound = F(report["source_cube"]["splitting_operator_upper_exact"])
    assert np.linalg.norm(U-target, 2) < float(bound)
    number = np.diag([mask.bit_count() for mask in range(256)])
    assert np.linalg.norm(U@number-number@U, 2) < 1e-10
    assert np.linalg.norm(U.conj().T@U-np.eye(256), 2) < 1e-10
    assert report["source_cube"]["all_Fock_columns"] == 256


def test_prepared_vacuum_is_half_filled_and_global_phase_is_not_discarded(execution):
    report, arrays = execution
    state = arrays["source_filled_vacuum"]
    assert sum(abs(state[j])**2 for j in range(256) if j.bit_count() != 4) < 1e-24
    H = n._fock_hamiltonian(arrays["source_cube_one_body_H"])
    ground = np.linalg.eigvalsh(H)[0]
    assert np.linalg.norm(H@state-ground*state) < 1e-10
    assert ground < 0 and H[0, 0] == 0
    forged = dict(arrays)
    forged["source_vacuum_frame"] = -arrays["source_vacuum_frame"]
    forged["source_filled_vacuum"] = -state
    with pytest.raises(ValueError, match="phase gate"):
        n.verify_retained_arrays(report, forged)


def test_complete_continuing_instrument_keeps_backaction_and_reference_columns(execution):
    report, arrays = execution
    K = arrays["source_history_kraus"]
    assert K.shape == (4, 256, 256)
    np.testing.assert_allclose(sum(k.conj().T@k for k in K), np.eye(256), atol=1e-10, rtol=0)
    p = [float(np.linalg.norm(k@arrays["source_filled_vacuum"])**2) for k in K]
    np.testing.assert_allclose(p, report["continuing_observer"]["record_probabilities"], atol=1e-13)
    assert report["continuing_observer"]["first_read_backaction_and_record_feedback_retained"]


@pytest.mark.parametrize("spacing", (F(1, 2), F(1, 4), F(1, 8)))
def test_continuum_band_and_propagator_bounds_follow_actual_staggered_symbol(spacing):
    momentum = (F(1), F(1, 2), F(-2, 3))
    h = np.asarray(n.corner_symbol(momentum, spacing=spacing).evalf(), complex)
    matrices = tuple(np.asarray(g, complex) for g in n.clifford_matrices())
    p = n.native_parameters()
    continuum = float(p["mass"])*matrices[3]+sum((float(p["coordinate_speed"]*k)*g
        for k, g in zip(momentum, matrices[:3], strict=True)), np.zeros((8, 8), complex))
    contract = n.continuum_contract(spacing=spacing, momentum_cutoff=F(1), time=F(2))
    assert np.linalg.norm(h-continuum, 2) <= float(F(contract["one_body_operator_error_upper"]))
    assert np.linalg.norm(expm(-2j*h)-expm(-2j*continuum), 2) <= float(F(contract["propagator_operator_error_upper"]))
    assert not contract["infinite_volume_vacuum_vector_norm_convergence_claimed"]


def test_excitation_identification_is_explicit_and_not_a_local_detector_claim():
    result = n.fiber_transport_certificate()
    assert result["transported_generator_norm_diagnostic"] == pytest.approx(result["energy_difference_diagnostic"], abs=1e-14)
    assert not result["unchanged_local_detectors_covered_by_excitation_bound"]
    with pytest.raises(ValueError, match="Brillouin"):
        n.continuum_contract(spacing=F(2), momentum_cutoff=F(1))


def test_free_field_normalization_and_common_geometry_cannot_repair_the_native_ratio():
    result = n.hierarchy_certificate()
    assert result["ascending_mass_ratio"] == [1, 3, 12]
    assert not result["constant_positive_field_rescaling_changes_pole_zeros"]
    for spacing in (F(1), F(1, 2), F(1, 4)):
        for sector in range(3):
            H = n.corner_symbol((F(),)*3, spacing=spacing, sector=sector)
            assert H*H == s.Rational(n.native_parameters(sector)["mass"]**2)*s.eye(8)


@pytest.mark.parametrize("bad", (True, .5, float("nan"), -1))
def test_original_exact_parameter_domains_refuse(bad):
    with pytest.raises(ValueError):
        n.source_cube(spacing=bad)


def test_warm_wrong_syndrome_and_wrong_factor_order_requests_refuse():
    word = n.echo_word(F(1, 3))
    n.decoded_occupation_word(word)
    original = n.compiler.remap
    def syndrome(tape, sites):
        literal = original(tape, sites)
        return [{**event, "targets": [target+1 for target in event["targets"]]} for event in literal]
    with patch.object(n.compiler, "remap", side_effect=syndrome):
        with pytest.raises(ArithmeticError, match="syndrome"):
            n.decoded_occupation_word(word)


def test_retained_operator_gate_and_replay_reject_dense_phase_zero_and_schema_forgery(execution):
    report, arrays = execution
    for name in ("source_cube_unitary", "source_vacuum_frame", "source_history_kraus"):
        forged = dict(arrays); forged[name] = arrays[name]*-1
        with pytest.raises(ValueError):
            n.compare_execution(report, forged, report, arrays)
    dense = dict(arrays); dense["source_cube_unitary"] = arrays["source_cube_unitary"]+1e-10
    with pytest.raises(ValueError, match="operator differs"):
        n.compare_execution(report, dense, report, arrays)
    changed = deepcopy(report)
    changed["source_contract"]["local_complete_macro_operator_admission"] = "1"
    with pytest.raises(ValueError):
        n.compare_execution(changed, arrays, report, arrays)
    changed = deepcopy(report)
    changed["backend_diagnostics"]["maximum_macro_operator_residual"] = 1e-9
    with pytest.raises(ValueError, match="own original"):
        n.compare_execution(changed, arrays, report, arrays)


def test_JSON_roundtrip_is_exact_schema_and_custody_refuses_duplicates_nonfinite(tmp_path, execution):
    report, arrays = execution
    n.compare_execution(json.loads(json.dumps(report)), arrays, report, arrays)
    path = tmp_path/"receipt.json"
    path.write_text('{"source_bindings":{},"source_bindings":{}}')
    with pytest.raises(ValueError, match="duplicate JSON"):
        n.check(path)
    path.write_text('{"bad":NaN}')
    with pytest.raises(ValueError, match="nonfinite"):
        n.check(path)
