"""Independent exterior, charge-pairing and tree-spectrum checks."""
from collections import Counter
import json

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import action


@pytest.fixture(scope="module")
def packet():
    return action.report()


def _numeric(matrix):
    return np.array(matrix.evalf(),dtype=complex)


def test_actual_all_left_charge_conjugate_order_and_no_mirror_addition(packet):
    rows = packet["ordering"]["all_left_components"]
    assert Counter(row["internal_type"] for row in rows) == {"Q":6,"u^c":3,"d^c":3,"L":2,"e^c":1}
    assert {row["electric_Q"] for row in rows if row["internal_type"]=="u^c"} == {"-2/3"}
    assert {row["electric_Q"] for row in rows if row["internal_type"]=="e^c"} == {"1"}
    assert packet["ordering"]["massless_neutrino_internal_family_indices"] == [30,31,32]
    assert not packet["interpretation"]["right_handed_mirror_fields_added_to_make_mass"]


def test_exterior_derivative_matches_independent_finite_determinant_representation():
    data = action.matrices()
    theta = .173
    u5 = expm(theta*_numeric(data["generators5"][0]))
    subsets = [row["exterior_subset"] for row in action.ordering()["all_left_components"]]
    finite = np.zeros((15,15),complex)
    for i,left in enumerate(subsets):
        for j,right in enumerate(subsets):
            if len(left)==len(right):
                finite[i,j] = np.linalg.det(u5[np.ix_(left,right)])
    target = expm(theta*_numeric(data["generators15"][0]))
    np.testing.assert_allclose(finite,target,rtol=0,atol=4e-16)
    np.testing.assert_allclose(finite.conj().T@finite,np.eye(15),rtol=0,atol=6e-16)


def test_neutral_mass_pairs_opposite_electric_charges_and_leaves_neutrinos_empty():
    k,charge = action.mass_kernel(),action.electric_charge()
    assert k.T == k
    assert charge.T*k+k*charge == s.zeros(45)
    for (i,j),value in k.todok().items():
        assert charge[i,i]+charge[j,j] == 0
        assert i//3 != j//3  # No artificial self-mass in a charged Weyl slot.
    for i in (30,31,32):
        assert k.row(i) == s.zeros(1,45)
    # A constant identity mass across this chiral table violates electric charge.
    assert charge.T+charge != s.zeros(45)
    # A symmetric all-left mass is a bilinear, rather than a charge-commuting endomorphism.
    assert charge*k-k*charge != s.zeros(45)


def test_creation_pairing_has_car_sign_and_full_nambu_electric_symmetry():
    delta = action.pairing_kernel()
    assert delta.shape == (90,90) and delta.T == -delta
    q = s.kronecker_product(s.eye(2),action.electric_charge())
    assert q*delta+delta*q.T == s.zeros(90)
    # At p=0 the stated BdG convention uses both pairing and its adjoint.
    d = _numeric(delta)
    h = np.block([[np.zeros_like(d),d],[d.conj().T,np.zeros_like(d)]])
    qn = np.block([[_numeric(q),np.zeros_like(d)],[np.zeros_like(d),-_numeric(q).T]])
    np.testing.assert_allclose(h,h.conj().T,rtol=0,atol=0)
    np.testing.assert_allclose(h@qn,qn@h,rtol=0,atol=4e-17)


def test_independent_dense_singular_values_have_expected_color_pairing_multiplicities(packet):
    k = _numeric(action.mass_kernel())
    singular = np.sort(np.linalg.svd(k,compute_uv=False))
    expected = [0.]*3
    for row in packet["tree_parameters"]:
        multiplicity = 2 if row["channel"]=="charged_lepton" else 6
        for expression in row["tree_mass_parameters_exact"]:
            expected.extend([float(s.sympify(expression))]*multiplicity)
    assert len(expected)==45
    np.testing.assert_allclose(singular,np.sort(expected),rtol=0,atol=3e-16)
    assert np.linalg.matrix_rank(k,tol=1e-12)==42
    assert np.max(singular)/np.min(singular[singular>1e-12]) < 15/4


def test_higgs_off_and_exact_scale_controls_are_explicit_architecture_changes():
    assert action.mass_kernel(0) == s.zeros(45)
    assert action.pairing_kernel(0) == s.zeros(90)
    assert action.mass_kernel(s.Rational(2)) == 2*action.mass_kernel(s.Rational(1))
    for value in (-1,s.I,s.Symbol("v",positive=True),s.oo):
        with pytest.raises(ValueError):
            action.mass_kernel(value)


def test_exact_current_covariance_anomaly_and_hierarchy_bound_are_scoped(packet):
    identities = packet["identities"]
    assert identities["one_generation_current_Higgs_identities"]==72
    assert identities["native_family_current_Higgs_identities"]==72
    assert set(identities["anomaly_traces_one_generation"].values()) == {"0"}
    assert identities["weak_doublets_one_generation"]==4
    assert packet["mass_hierarchy"]["same_family_singular_value_ratios_in_all_three_channels"]
    assert not packet["mass_hierarchy"]["observed_charged_fermion_hierarchy_reproduced"]
    assert not packet["interpretation"]["physical_Standard_Model_pole_masses_derived"]
    assert not packet["interpretation"]["Yukawa_architecture_forced_by_repair"]
    json.dumps(packet,allow_nan=False)


def test_public_matrices_and_coefficient_inputs_cannot_be_corrupted():
    first = action.source_coefficients()
    first["a"] = s.Integer(0)
    assert action.source_coefficients()["a"] == s.log(s.Rational(25,8))
    k = action.mass_kernel()
    with pytest.raises(TypeError):
        k[0,0]=1


def test_full_spin_mode_labels_select_only_three_of_six_lepton_pair_blocks():
    labels = action.mode_labels()
    charge,delta = action.spin_electric_charge(),action.pairing_kernel()
    assert [row["index"] for row in labels] == list(range(90))
    probes = [row["index"] for row in labels if row["spin"]==0 and charge[row["index"],row["index"]]==-1]
    assert probes == [33,34,35]
    assert [[j for j in range(90) if delta[i,j]!=0] for i in probes] == [[72],[73],[74]]
    assert sum(row["internal_type"] in ("L","e^c") and row["electric_Q"]!="0" for row in labels)==12
    assert sum(row["electric_Q"]=="0" for row in labels)==6
