"""Full-complement, derivative and uniform-scheduler controls."""
import json

import numpy as np
import pytest
import sympy as s

from . import coupling as c
from ...native_modular_reference import local_preparation


def numeric(matrix):
    return np.array(matrix.evalf(17), dtype=float)


def test_source_D_is_complete_weighted_seam_energy_with_protected_kernel():
    data = c._data()
    selection = local_preparation.selectors(32)
    independent = s.zeros(12)
    for (a, b), weight in selection["completion_accepted_seam_weights"].items():
        edge = s.zeros(12, 1)
        edge[a], edge[b] = 1, -1
        independent += s.Rational(weight)/2*edge*edge.T
    assert data["D"] == independent
    assert data["D"].rank() == 9
    assert data["D"]*s.eye(12)[:, 0] == s.zeros(12, 1)
    assert data["D"]*s.eye(12)[:, 1] == s.zeros(12, 1)
    assert data["D"]*s.Matrix([0, 0]+[1]*10) == s.zeros(12, 1)


def test_exact_harmonic_extension_keeps_all_port_constraints_and_energy():
    data = c._data()
    b, u, d, g, static = (data[key] for key in ("B", "harmonic", "D", "G", "static"))
    assert c._zero(b.T*u-g*s.eye(3))
    assert c._zero(s.ones(1, 12)*u)
    assert c._zero(d*u-b*static)
    assert c._zero(u.T*d*u/g-static)
    assert static.rank() == 1
    assert static[0, 0] == static[1, 1] == 0
    assert float(static[2, 2]) == pytest.approx(.03717449835323584)
    # Independent admissible perturbation: its energy cross term vanishes,
    # and the retained full seam form adds only its nonnegative own energy.
    complement = (b.T.col_join(s.ones(1, 12))).nullspace()
    delta = complement[0]+2*complement[-1]
    vector = s.Matrix([2, -3, 1])
    baseline, perturbed = u*vector, u*vector+delta
    difference = c._exact(perturbed.T*d*perturbed-baseline.T*d*baseline-delta.T*d*delta)
    assert difference == s.zeros(1)
    assert float((delta.T*d*delta)[0]) > 0


def test_frequency_derivative_is_full_harmonic_gram_and_cannot_lift_zeros():
    data = c._data()
    z, b, u, g = (data[key] for key in ("derivative", "B", "harmonic", "G"))
    assert c._zero(z-s.eye(3)-(u-b).T*(u-b)/g)
    assert np.linalg.eigvalsh(numeric(z))[0] >= 1-1e-12
    # Direct canonical normalization independently preserves two mass zeros.
    eig, basis = np.linalg.eigh(numeric(z))
    invroot = basis@np.diag(1/np.sqrt(eig))@basis.T
    canonical = invroot@numeric(data["static"])@invroot
    eigenvalues = np.linalg.eigvalsh(canonical)
    np.testing.assert_allclose(eigenvalues[:2], 0, atol=1e-14)
    assert eigenvalues[2] == pytest.approx(float(data["canonical_gap"]), rel=1e-12)


def test_finite_frequency_full_inverse_confirms_the_first_derivative():
    data = c._data()
    errors = []
    for z in (s.Rational(1, 1000), s.Rational(1, 2000)):
        full_inverse = c.family_inverse_response(z)
        remainder = numeric((full_inverse-data["static"])/z-data["derivative"])
        errors.append(np.linalg.norm(remainder))
    assert errors[1] < .6*errors[0]
    assert errors[1] > 0
    for invalid in (0, -1, .1, s.sqrt(2)):
        with pytest.raises(ValueError):
            c.family_inverse_response(invalid)


def test_uniform_native_continuation_has_scalar_static_response_and_unit_kinetic():
    uniform = c.uniform_repair_control()
    data = c._data()
    assert c._zero(uniform["D12"]*data["B"]-data["B"]*uniform["coefficient"])
    assert uniform["static"] == (5-s.sqrt(5))/60*s.eye(3)
    assert uniform["derivative"] == s.eye(3)
    assert uniform["D12"].rank() == 11


def test_fixed_Higgs_tensor_attachment_has_correct_current_and_CAR_ranks():
    packet = c.higgs_attachment_report()
    assert packet["retained_internal_current_Higgs_identities_rechecked_exact"] == 72
    assert packet["K45_rank_exact_by_disjoint_nonzero_matching"] == 14
    assert packet["K45_zero_dimension"] == 31
    assert packet["Delta90_rank_exact"] == 28
    assert packet["massive_Dirac_fields_counting_color"] == 7
    assert packet["massless_charged_Dirac_fields_counting_color"] == 14
    kernel = c.higgs_mass_kernel()
    np.testing.assert_array_equal(np.count_nonzero(numeric(kernel), axis=1) != 0,
                                  np.count_nonzero(numeric(kernel), axis=0) != 0)
    assert np.linalg.matrix_rank(numeric(kernel)) == 14


def test_report_separates_euclidean_form_from_declared_Higgs_assignment():
    packet = c.report()
    json.dumps(packet, allow_nan=False)
    harmonic, assignment = packet["harmonic_elimination"], packet["declared_Higgs_coupling_control"]
    assert harmonic["static_rank"] == 1
    assert harmonic["positive_derivative_at_least_identity_exact"]
    assert not harmonic["frequency_derivative_lifts_static_zero_modes"]
    assert all(gaps[:2] == ["0", "0"] for gaps in assignment["source_predicted_family_gaps_over_clock"].values())
    assert not assignment["assignment_forced_by_native_repair"]
    assert not assignment["all_three_charged_family_ratios_defined"]
    assert not packet["interpretation"]["native_coherent_left_right_Higgs_vertex_derived"]
    assert not packet["interpretation"]["Dirichlet_or_Laplace_coefficients_are_automatically_Lorentzian_masses"]
