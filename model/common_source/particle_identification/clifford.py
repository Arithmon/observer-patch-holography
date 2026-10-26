"""Classify the separate free continuum candidate without assigning SM names.

The eight-component regulator in mass_response.dispersion has a continuum
Clifford-four action with multiplicity two.  This module constructs its spin,
flavor commutant, energy projectors and a kinematic charge conjugation.
Its finite matrix identities do not establish physical CAR, electromagnetic
normalization, an interacting field net or an observed particle identity.

Exact statements are derived algebraically in CLIFFORD.md. Floating residuals
and singular-value ranks below are checks, not outward numerical certificates.
"""
from __future__ import annotations

import itertools
import math
from numbers import Real

import numpy as np

from ..mass_response.dispersion import bloch_hamiltonian, gamma_matrices
from ..mass_response.model import parameters


def _momentum(value):
    try:
        result = tuple(itertools.islice(iter(value), 4))
    except TypeError as exc:
        raise ValueError("three finite real momentum components are required") from exc
    if len(result) != 3:
        raise ValueError("three finite real momentum components are required")
    out = []
    for coordinate in result:
        if isinstance(coordinate, (bool, np.bool_)) or not isinstance(coordinate, Real):
            raise ValueError("three finite real momentum components are required")
        try:
            number = float(coordinate)
        except (OverflowError, ValueError) as exc:
            raise ValueError("momentum must have a finite execution representation") from exc
        if not math.isfinite(number):
            raise ValueError("three finite real momentum components are required")
        out.append(number)
    return tuple(out)


def continuum_matrices():
    """Return independent copies of the declared continuum matrix dictionary."""
    gamma = gamma_matrices()
    alpha = tuple(gamma[index] for index in (1, 3, 5))
    beta = (gamma[0] + gamma[2] + gamma[4]) / math.sqrt(3)
    eta1 = (gamma[0] - gamma[2]) / math.sqrt(2)
    eta2 = (gamma[0] + gamma[2] - 2 * gamma[4]) / math.sqrt(6)
    grading = alpha[0] @ alpha[1] @ alpha[2] @ beta
    flavor = (grading @ eta1, grading @ eta2, -1j * eta1 @ eta2)
    spin = tuple(-.5j * alpha[j] @ alpha[k] for j, k in ((1, 2), (2, 0), (0, 1)))
    lorentz_gamma = (beta, *(beta @ value for value in alpha))
    chirality = -1j * alpha[0] @ alpha[1] @ alpha[2]
    return {"alpha": alpha, "beta": beta, "transverse": (eta1, eta2),
            "grading": grading, "flavor": flavor, "spin": spin,
            "lorentz_gamma": lorentz_gamma, "chirality": chirality,
            "charge_conjugation_matrix": eta1}


def spin_generators():
    return continuum_matrices()["spin"]


def commutant_basis():
    """Exact-form basis I,T1,T2,T3 of the continuum internal commutant."""
    return (np.eye(8, dtype=complex), *continuum_matrices()["flavor"])


def flavor_projectors():
    """One chosen decomposition into the two equivalent Dirac copies.

    The chosen flavor axis is not a physical generation label. All U(2)
    changes of multiplicity basis describe the same continuum operator.
    """
    t3 = continuum_matrices()["flavor"][2]
    return ((np.eye(8) + t3) / 2, (np.eye(8) - t3) / 2)


def source_scales():
    """Rest-energy and speed coefficients from the actual reference producer."""
    p = parameters()
    return {"a": p.intra, "b": p.inter, "imbalance": p.intra - p.inter,
            "clock_frequency": p.clock,
            "candidate_rest_energy": math.sqrt(3) * abs(p.intra - p.inter),
            "candidate_speed": p.clock / 2,
            "candidate_rest_energy_over_clock": math.sqrt(3) * abs(p.intra - p.inter) / p.clock,
            "units": "declared continuum coordinate and programmed modular clock",
            "physical_mass_or_light_speed_identified": False}


def continuum_hamiltonian(momentum):
    """H=c alpha.p+M beta on the separate free continuum candidate."""
    p = _momentum(momentum)
    matrices = continuum_matrices()
    scales = source_scales()
    result = scales["candidate_rest_energy"] * matrices["beta"]
    for value, alpha in zip(p, matrices["alpha"], strict=True):
        result = result + scales["candidate_speed"] * value * alpha
    return result


def energy_projectors(momentum):
    p = _momentum(momentum)
    scales = source_scales()
    energy = math.sqrt(scales["candidate_rest_energy"] ** 2
                       + scales["candidate_speed"] ** 2 * math.fsum(x * x for x in p))
    h = continuum_hamiltonian(p)
    return ((np.eye(8) + h / energy) / 2, (np.eye(8) - h / energy) / 2)


def clifford_words(generators):
    """Ordered monomial basis, including identity, for the displayed generators."""
    generators = tuple(np.asarray(value, complex) for value in generators)
    if not generators or any(value.shape != generators[0].shape for value in generators):
        raise ValueError("a nonempty equal-size matrix family required")
    if generators[0].ndim != 2 or generators[0].shape[0] != generators[0].shape[1]:
        raise ValueError("square generators required")
    out = []
    for mask in range(1 << len(generators)):
        value = np.eye(len(generators[0]), dtype=complex)
        for index, generator in enumerate(generators):
            if mask >> index & 1:
                value = value @ generator
        out.append(value)
    return tuple(out)


def commutant_nullity(generators, tolerance=1e-10):
    """Independent floating complex-linear commutant check.

    Exact dimensions used in report() also have the constructive algebraic
    proof in CLIFFORD.md; they are not inferred solely from this threshold.
    """
    generators = tuple(np.asarray(value, complex) for value in generators)
    if not generators or any(value.shape != generators[0].shape for value in generators):
        raise ValueError("a nonempty equal-size matrix family required")
    if generators[0].ndim != 2:
        raise ValueError("finite square matrices required")
    size = generators[0].shape[0]
    if generators[0].shape != (size, size) or any(not np.isfinite(value).all() for value in generators):
        raise ValueError("finite square matrices required")
    if not isinstance(tolerance, Real) or isinstance(tolerance, bool) or not 0 < tolerance < 1:
        raise ValueError("a rank threshold strictly between zero and one required")
    identity = np.eye(size)
    constraints = np.concatenate([np.kron(value.T, identity) - np.kron(identity, value)
                                  for value in generators], axis=0)
    values = np.linalg.svd(constraints, compute_uv=False)
    return size * size - int(np.count_nonzero(values > tolerance))


def lattice_probe_matrices():
    """Seven momenta already generate the full lattice Clifford-six algebra.

    This is a universal momentum-family check, not the commutant of one
    eight-by-eight Hamiltonian or every special small-volume momentum set.
    """
    out = [bloch_hamiltonian((0., 0., 0.))]
    for axis in range(3):
        for sign in (-1, 1):
            momentum = [0., 0., 0.]
            momentum[axis] = sign * math.pi / 2
            out.append(bloch_hamiltonian(momentum))
    return tuple(out)


def _norm(matrix):
    return float(np.linalg.norm(matrix, 2))


def _gram_error(matrices):
    vectors = np.array([value.ravel() for value in matrices])
    gram = vectors.conj() @ vectors.T / len(matrices[0])
    return float(np.max(abs(gram - np.eye(len(matrices)))))


def report():
    """Compact JSON-ready kinematic classification with explicit boundaries."""
    data = continuum_matrices()
    alpha, beta, flavor, spin = data["alpha"], data["beta"], data["flavor"], data["spin"]
    identity = np.eye(8)
    continuum_generators = (*alpha, beta)
    words = clifford_words(continuum_generators)
    tensor_basis = tuple(word @ factor for word in words for factor in (identity, *flavor))
    casimir = sum(value @ value for value in spin)
    p = (.31, -.27, .43)
    h = continuum_hamiltonian(p)
    positive, negative = energy_projectors(p)
    flavors = flavor_projectors()
    charge_c = data["charge_conjugation_matrix"]
    lattice = lattice_probe_matrices()
    scales = source_scales()
    gamma_lorentz = data["lorentz_gamma"]
    signs = (1., -1., -1., -1.)
    lorentz_error = max(_norm(a @ b + b @ a - (2 * signs[i] * identity if i == j else 0))
                        for i, a in enumerate(gamma_lorentz) for j, b in enumerate(gamma_lorentz))
    return {
        "schema": "oph.particle_identification.clifford.v1",
        "candidate": "separate FREE eight-component three-axis continuum extension",
        "source_scales": scales,
        "kinematic_classification": {
            "lorentz_metric_convention": "+---",
            "spin_casimir": "3/4 I8",
            "spin": "1/2 representation; physical spin identification is separate",
            "continuum_dirac_multiplicity": 2,
            "mass_degeneracy": 2,
            "positive_energy_states": int(round(np.trace(positive).real)),
            "negative_energy_states": int(round(np.trace(negative).real)),
            "chosen_flavor_projector_ranks": [int(round(np.trace(f).real)) for f in flavors],
            "spin_states_per_flavor_and_energy_sign": [
                int(round(np.trace(f @ e).real)) for f in flavors for e in (positive, negative)],
            "continuum_clifford_algebra_dimension": 16,
            "continuum_commutant_complex_dimension_exact": 4,
            "continuum_commutant_nullity_numeric": commutant_nullity(continuum_generators),
            "continuum_internal_commutant": "M2(C), unitary group U(2)",
            "physical_generation_or_flavor_assignment": False,
        },
        "particle_antiparticle_contract": {
            "energy_branches": "+E and -E, each with two spin states per Dirac copy",
            "charge_conjugation": "C=eta1 times complex conjugation, C^2=I",
            "global_phase_symmetry": "psi -> exp(i q theta) psi; q has arbitrary normalization",
            "fock_interpretation": "CAR plus a filled negative-band reference gives particles of +q and holes of -q",
            "car_field_and_vacuum_joined_to_executed_ring": False,
            "electromagnetic_charge_identified": False,
            "statistics_derived_from_spin_matrices": False,
            "positive_and_negative_bands_are_not_particles_before_vacuum_attachment": True,
            "zero_temperature_neutral_bilinear_threshold_lower_bound": 2 * scales["candidate_rest_energy"],
            "bilinear_threshold_scope": "nonzero connected particle-hole response above the filled negative-band vacuum; matrix elements can suppress the threshold, and total density at zero momentum is conserved",
            "isolated_scalar_or_vector_bound_state_proved": False,
        },
        "regulator_boundary": {
            "eight_components_retained": True,
            "lattice_all_momentum_algebra_dimension_exact": 64,
            "lattice_all_momentum_commutant_dimension_exact": 1,
            "lattice_commutant_nullity_numeric": commutant_nullity(lattice),
            "lattice_commutant_scope": "momentum-independent 8x8 internal matrices commuting with every H_lat(k); not the full finite-space spectral commutant or momentum-dependent nonlocal symmetries",
            "continuum_U2_flavor_is_exact_on_finite_lattice": False,
            "special_small_volume_momentum_sets_may_have_larger_commutants": True,
            "finite_lattice_exact_continuous_Lorentz_symmetry": False,
            "continuum_kinematic_lorentz_spinor_representation": True,
            "massless_single_particle_pole_in_this_free_gapped_continuation": False,
        },
        "standard_model_boundary": {
            "observed_species_identified": [],
            "photon_or_higgs_identified": False,
            "nontrivial_unbroken_su3_color_in_this_dirac_multiplicity": False,
            "su3_reason": "A nonzero homomorphism from simple su(3) into u(2) would be injective, but 8>4. This excludes only an internal color action commuting with this continuum Dirac operator on this two-dimensional multiplicity.",
            "remaining": ["physical current and charge normalization", "chiral gauge representation and statistics attachment",
                          "observed nondegenerate masses and family content", "interacting shared-history propagator and refinement",
                          "physical clock and observer-field identification"],
        },
        "numerical_identity_residuals": {
            "lorentz_gamma_anticommutation": lorentz_error,
            "spin_casimir": _norm(casimir - .75 * identity),
            "spin_su2_commutators": max(_norm(spin[i] @ spin[j] - spin[j] @ spin[i] - 1j * spin[k])
                                        for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1))),
            "spin_commutes_with_beta": max(_norm(s @ beta - beta @ s) for s in spin),
            "flavor_commutes_with_continuum_generators": max(_norm(t @ a - a @ t) for t in flavor for a in continuum_generators),
            "flavor_pauli_products": max(_norm(flavor[i] @ flavor[j] - 1j * flavor[k])
                                         for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1))),
            "clifford4_word_gram": _gram_error(words),
            "complete_64_operator_basis_gram": _gram_error(tensor_basis),
            "clifford6_word_gram": _gram_error(clifford_words(gamma_matrices())),
            "charge_conjugation_square": _norm(charge_c @ charge_c.conj() - identity),
            "charge_conjugation_energy_reversal": _norm(charge_c @ h.conj() @ charge_c.conj().T
                                                        + continuum_hamiltonian(tuple(-x for x in p))),
        },
        "physical_particle_identification_complete": False,
        "numerical_checks_are_outward_certificates": False,
    }
