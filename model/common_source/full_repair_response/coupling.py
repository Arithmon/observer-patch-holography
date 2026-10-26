"""Harmonic full-repair elimination and a declared Higgs-coupling control.

The source fixes a positive one-proposal Dirichlet form and its complete port
complement. Eliminating that complement produces a rank-one static family
kernel and a positive frequency derivative. These are not automatically a
Lorentzian kinetic term or a Yukawa: assigning those roles is explicitly a
tested architecture prescription, with two exact zero charged family gaps.
"""
from __future__ import annotations

from functools import lru_cache

import sympy as s
from sympy.polys.matrices import DomainMatrix

from ..chiral_candidate import action
from ..hierarchy_diagnosis import native_family_dynamics
from .. import slow_geometry


_FIELD = s.QQ.algebraic_field(s.sqrt(5))


def _exact(matrix):
    return s.ImmutableMatrix(DomainMatrix.from_Matrix(s.Matrix(matrix)).convert_to(_FIELD).to_Matrix())


def _inverse(matrix):
    return _exact(DomainMatrix.from_Matrix(s.Matrix(matrix)).convert_to(_FIELD).inv().to_Matrix())


def _zero(matrix):
    return _exact(matrix) == s.zeros(*matrix.shape)


def _entries(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


def _scalar(value):
    value = _exact(s.Matrix([[value]]))[0, 0]
    return {"exact": str(value), "float_diagnostic": float(value)}


@lru_cache(maxsize=2, typed=True)
def _data(selector_rounds=32):
    source = native_family_dynamics._data(selector_rounds)
    transition, basis = source["T"], source["basis"]
    gram = _exact(basis.T*basis)[0, 0]
    d = _exact(s.eye(12)-transition)
    if not _zero(d-d.T) or not _zero(d*s.ones(12, 1)):
        raise ArithmeticError("a symmetric constant-preserving source repair is required")
    if not _zero(d[:2, :]) or not _zero(d[:, :2]):
        raise ArithmeticError("the actual completion must protect both core ports")
    outer = s.Matrix(d[2:, 2:])
    constant_projector = s.ones(10)/10
    outer_pinv = _exact((outer+constant_projector).inv(method="DM")-constant_projector)
    if not _zero(outer*outer_pinv-(s.eye(10)-constant_projector)):
        raise ArithmeticError("the retained outer graph has an additional zero direction")
    # Impose B^T u=G*x and 1^T u=0; minimize u^T D u over the entire
    # complementary port space. Core ports can freely satisfy x/y constraints.
    phi = (1+s.sqrt(5))/2
    bz = basis[2:, 2]
    susceptibility = _exact(bz.T*outer_pinv*bz)[0, 0]
    gap = _exact(s.Matrix([[gram/susceptibility]]))[0, 0]
    outer_z = _exact(gap*outer_pinv*bz)
    mx, my = (_exact(basis[2:, j].T*outer_z)[0, 0] for j in (0, 1))
    harmonic = s.zeros(12, 3)
    harmonic[0, 0], harmonic[1, 0] = -gram/2, gram/2
    harmonic[0, 1] = harmonic[1, 1] = 5*gram/(12*phi)
    for i in range(2, 12):
        harmonic[i, 1] = -gram/(12*phi)
    shift = my/(12*phi)
    harmonic[0, 2], harmonic[1, 2] = -5*my/(12*phi)+mx/2, -5*my/(12*phi)-mx/2
    for i in range(2, 12):
        harmonic[i, 2] = outer_z[i-2]+shift
    harmonic = _exact(harmonic)
    static = _exact(harmonic.T*d*harmonic/gram)
    derivative = _exact(harmonic.T*harmonic/gram)
    outside = _exact(harmonic-basis)
    if not _zero(basis.T*harmonic-gram*s.eye(3)) or not _zero(s.ones(1, 12)*harmonic):
        raise ArithmeticError("the full harmonic extension does not preserve its family coordinates")
    if static != s.diag(0, 0, gap) or not _zero(d*harmonic-basis*static):
        raise ArithmeticError("the full static Schur kernel did not have the predicted rank")
    if not _zero(derivative-s.eye(3)-outside.T*outside/gram):
        raise ArithmeticError("the exact positive frequency-derivative Gram identity failed")
    # Positive canonical normalization preserves rank. The one nonzero
    # generalized coefficient of S0 v = lambda Z v is gap*(Z^-1)_zz.
    inverse_derivative = _inverse(derivative)
    canonical_gap = _exact(s.Matrix([[gap*inverse_derivative[2, 2]]]))[0, 0]
    return {"source": source, "T": transition, "D": d, "B": basis, "G": gram,
            "outer_pinv": outer_pinv, "susceptibility": susceptibility,
            "harmonic": harmonic, "static": static, "derivative": derivative,
            "outside": outside, "canonical_gap": canonical_gap}


def static_family_kernel(selector_rounds=32):
    return _data(selector_rounds)["static"]


def frequency_derivative(selector_rounds=32):
    """Z=S'(0)=U_harm^T U_harm/G, positive and at least identity."""
    return _data(selector_rounds)["derivative"]


def family_inverse_response(laplace_parameter, selector_rounds=32):
    """S(z)=[B^T(D+zI)^-1 B/G]^-1 for an exact positive Laplace z.

    z is conjugate to the algebraic repair-proposal parameter. This is not a
    continued native trajectory after terminal timeout, or a particle energy.
    """
    z = s.sympify(laplace_parameter)
    if z.is_Rational is not True or z.is_positive is not True:
        raise ValueError("a strictly positive exact rational Laplace parameter required")
    data = _data(selector_rounds)
    full = (s.Matrix(data["D"])+z*s.eye(12)).inv(method="DM")
    reduced = _exact(data["B"].T*full*data["B"]/data["G"])
    return _inverse(reduced)


def uniform_repair_control():
    """Actual uniform all30-seam first moment: a closed, degenerate slow band."""
    data = _data()
    laplacian = s.Matrix(slow_geometry.source().laplacian().tolist())
    d = laplacian/60
    coefficient = (5-s.sqrt(5))/60
    if not _zero(d*data["B"]-data["B"]*coefficient):
        raise ArithmeticError("the uniform native repair slow band changed")
    return {"D12": s.ImmutableMatrix(d), "static": s.ImmutableMatrix(coefficient*s.eye(3)),
            "derivative": s.eye(3), "coefficient": coefficient}


def higgs_mass_kernel(selector_rounds=32):
    """Actual K45(v=Omega) using the retained Higgs lines and new S0 texture."""
    static = static_family_kernel(selector_rounds)
    replacement = s.kronecker_product(s.eye(15), action.family_feature().inv()*static)
    # The retained kernel is sum C_channel tensor coefficient_channel*F.
    # Multiplication changes only F to S0 and preserves each Higgs line.
    return s.ImmutableSparseMatrix(action.mass_kernel()*replacement)


def higgs_attachment_report(selector_rounds=32):
    """Verify the actual source Higgs/current identities, not charge names alone."""
    _, primitive, currents, _, channels = action._representation()
    count = 0
    for channel in channels.values():
        components = channel["components"]
        for source, current in zip(primitive, currents, strict=True):
            h = source[3:, 3:]
            if channel["conjugate_Higgs"]:
                h = -h.T
            for component in (0, 1):
                correction = sum((h[j, component]*components[j] for j in (0, 1)), s.zeros(15))
                identity = current.T*components[component]+components[component]*current+correction
                if identity.applyfunc(s.expand) != s.zeros(15):
                    raise ArithmeticError("the retained source Higgs/current tensor identity failed")
                count += 1
    kernel = higgs_mass_kernel(selector_rounds)
    q = action.electric_charge()
    if kernel != kernel.T or q.T*kernel+kernel*q != s.zeros(45):
        raise ArithmeticError("the new static texture failed all-left charge/symmetry checks")
    # Exactly one nonzero neutral-Higgs partner per occupied internal row;
    # all entries are products of strictly positive source factors up to sign.
    entries = kernel.todok()
    nonzero_rows = {i for i, _ in entries}
    nonzero_cols = {j for _, j in entries}
    if len(entries) != 14 or len(nonzero_rows) != 14 or len(nonzero_cols) != 14:
        raise ArithmeticError("the expected rank-one-family Higgs matching changed")
    delta = s.kronecker_product(action.SPIN_EPSILON, kernel)
    if delta.T != -delta:
        raise ArithmeticError("the physical spin-pairing coefficient is not antisymmetric")
    return {
        "retained_internal_current_Higgs_identities_rechecked_exact": count,
        "family_lift_reason": "each verified internal identity tensors with the same symmetric coefficient*S0; gauge currents act as T15 tensor I3",
        "new_family_Yukawa_kernels_symmetric": True,
        "K45_symmetric_exact": True, "electric_charge_QT_K_plus_KQ_zero_exact": True,
        "K45_sparse_entries": [[i, j, str(value)] for (i, j), value in sorted(entries.items())],
        "K45_rank_exact_by_disjoint_nonzero_matching": 14,
        "K45_zero_dimension": 31,
        "Delta90_antisymmetric_exact": True, "Delta90_rank_exact": 28,
        "canonical_normalization_preserves_kernel_ranks": True,
        "massive_Dirac_fields_counting_color": 7,
        "massless_charged_Dirac_fields_counting_color": 14,
        "massless_left_Weyl_neutrino_fields": 3,
        "physical_gauge_field_or_actual_native_Yukawa_identification_proved": False,
    }


def report(selector_rounds=32):
    data = _data(selector_rounds)
    uniform = uniform_repair_control()
    coefficients = action.source_coefficients()
    gap = data["canonical_gap"]
    return {
        "schema": "oph.full_repair_response.harmonic_coupling.v1",
        "source": {"update": "existing finite mixed-cap completion proposal including its terminal selector-timeout flag",
                   "selector_rounds": selector_rounds,
                   "algebraic_form": "D=I-T12; T12 is the one-proposal unread work marginal",
                   "time_parameter": "algebraic repair proposal, not modular flow or calibrated physical time",
                   "terminal_timeout_probability": str(data["source"]["timeout"]),
                   "timeout_is_not_restarted": True,
                   "all_port_complement_modes_used": True,
                   "empirical_mass_or_coupling_inputs": []},
        "harmonic_elimination": {
            "definition": "S(z)=[G^(-1) B^T (D+zI)^(-1) B]^(-1)",
            "coordinate_constraint": "B^T U_harm=G I; 1^T U_harm=0",
            "stationarity": "D U_harm=B S0",
            "full_D_kernel": "span(e0,e1,ones_on_outer10); rank(D)=9",
            "static_family_kernel": _entries(data["static"]),
            "static_rank": 1,
            "static_kernel": "span(native family x,native family y)",
            "static_positive_coefficient": _scalar(data["static"][2, 2]),
            "harmonic_extension_12_by3": _entries(data["harmonic"]),
            "family_frequency_derivative_Z": _entries(data["derivative"]),
            "derivative_identity": "S'(0)=Z=U_harm^T U_harm/G=I+(U_harm-B)^T(U_harm-B)/G",
            "positive_derivative_at_least_identity_exact": True,
            "frequency_derivative_lifts_static_zero_modes": False,
            "nonzero_generalized_coefficient_of_S0_v_equals_lambda_Z_v": _scalar(gap),
            "rank_one_is_completion_protection_consequence": True},
        "declared_Higgs_coupling_control": {
            "assignment": "Y_u=a S0,Y_d=b S0,Y_e=J S0; common family kinetic Z on left and conjugated singlets; v=Omega",
            "assignment_forced_by_native_repair": False,
            "uses_the_existing_gauge_invariant_Higgs_tensors": True,
            "canonical_normalization": "Ycanonical=Z^(-1/2) Y Z^(-1/2), rank preserved",
            "source_predicted_family_gaps_over_clock": {
                name: ["0", "0", str(coefficients[key]*gap)]
                for name, key in (("up", "a"), ("down", "b"), ("charged_lepton", "J"))},
            "source_predicted_nonzero_gap_over_clock_float_diagnostic": {
                name: float(coefficients[key]*gap)
                for name, key in (("up", "a"), ("down", "b"), ("charged_lepton", "J"))},
            "three_nonzero_family_gaps": False,
            "canonical_singular_mass_ratios_normalized_to_nonzero_mass": [0, 0, 1],
            "all_three_charged_family_ratios_defined": False,
            "common_massive_family_axis_across_channels": True,
            "new_successful_observed_mass_spectrum_claimed": False},
        "actual_Higgs_tensor_attachment": higgs_attachment_report(selector_rounds),
        "uniform_all_seam_control": {
            "update": "uniform choice among all30 native seams and fair repair coin, one continuable native attempt",
            "D12": _entries(uniform["D12"]),
            "closed_band_relation": "D B = B*(5-sqrt5)/60",
            "static_family_kernel": _entries(uniform["static"]),
            "family_frequency_derivative_Z": _entries(uniform["derivative"]),
            "all_three_family_coefficients_equal_exact": True,
            "family_ratios_if_assigned_as_Yukawa": [1, 1, 1],
            "nonzero_coefficient": _scalar(uniform["coefficient"])},
        "interpretation": {
            "full_repair_complement_changes_the_static_family_kernel": True,
            "rank_seed_and_kinetic_coefficient_source_computed": True,
            "Dirichlet_or_Laplace_coefficients_are_automatically_Lorentzian_masses": False,
            "native_coherent_left_right_Higgs_vertex_derived": False,
            "frequency_derivative_alone_generates_missing_mass_scales": False,
            "reflection_positive_reconstruction_or_chiral_field_attachment_proved": False,
            "actual_interacting_charged_pole_response_implemented": False,
            "remaining_obligation": "derive the actual coherent Higgs/current coupling and its physical continuation; static protection gives two exact zeros, while unprotected uniform repair gives degeneracy, without adding another law"},
    }
