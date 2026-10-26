"""Exact connected-occupation covariance as a native family-texture control.

The reference distribution is faithful on all binary total-two configurations.
Its covariance is an equal-reference statistical observable, not a propagating
mass, a coherent interaction Hamiltonian, or a selected Yukawa map. Compressing
it onto the existing native rank-three band yields no new mixing basis.
"""
from __future__ import annotations

from functools import lru_cache

import sympy as s

from .. import slow_geometry
from ..particle_identification import family_attachment
from ...native_modular_reference import local_preparation
from ...observer import native


def _clean(matrix):
    return s.ImmutableMatrix(matrix).applyfunc(s.simplify)


def _zero(matrix):
    return all(s.simplify(value) == 0 for value in matrix)


def _entries(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


@lru_cache(maxsize=1)
def _data():
    states, probabilities = native.states(binary=True), local_preparation.reference()
    if len(states) != 66 or len(probabilities) != 66:
        raise ArithmeticError("the complete 66-state reference is required")
    probabilities = tuple(s.Rational(value) for value in probabilities)
    if sum(probabilities) != 1 or any(p <= 0 for p in probabilities):
        raise ArithmeticError("a faithful normalized reference is required")
    if any(len(x) != 12 or sum(x) != 2 or any(bit not in (0, 1) for bit in x) for x in states):
        raise ArithmeticError("binary fixed-total-two source domain changed")
    means, second = s.zeros(12, 1), s.zeros(12)
    for state, probability in zip(states, probabilities, strict=True):
        vector = s.Matrix(state)
        means += probability*vector
        second += probability*vector*vector.T
    covariance = _clean(second-means*means.T)
    basis = s.Matrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a, b in row]
                      for row in slow_geometry.raw_frame()])
    gram = _clean(basis.T*basis)
    if gram != (10+2*s.sqrt(5))*s.eye(3) or not _zero(s.ones(1, 12)*basis):
        raise ArithmeticError("the existing metric-orthogonal zero-sum family frame changed")
    left = gram.inv()*basis.T
    feature = _clean(left*s.diag(*means)*basis)
    compressed = _clean(left*covariance*basis)
    retained = family_attachment.report()["compression"]["matrix_Qsqrt5"]
    expected_feature = s.Matrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a, b in row]
                                for row in retained])
    if feature != expected_feature:
        raise ArithmeticError("source family chart does not match the retained feature")
    expected = s.diag(s.Rational(37, 225)-17*s.sqrt(5)/2250,
                      s.Rational(3283, 22500)+433*s.sqrt(5)/112500,
                      s.Rational(19, 150))
    if compressed != expected or not _zero(feature*compressed-compressed*feature):
        raise ArithmeticError("the exact connected covariance control changed")
    return (s.ImmutableMatrix(means), s.ImmutableMatrix(second), covariance,
            s.ImmutableMatrix(basis), gram, feature, compressed)


def port_covariance():
    """C_ij=E[n_i n_j]-E[n_i]E[n_j] on the actual twelve source ports."""
    return _data()[2]


def family_covariance():
    """G^-1 B^T C B, in the existing metric-orthogonal source family frame."""
    return _data()[6]


def report():
    """JSON-ready exact negative control; no empirical targets are consulted."""
    means, second, covariance, basis, gram, feature, compressed = _data()
    # Full rank/positivity also follow from faithfulness on all two-subsets:
    # zero variance means z_i+z_j is constant for every pair, hence all z_i equal.
    full_spectrum = {s.Rational(0): 1, s.Rational(19, 150): 9,
                     s.Rational(644, 3375): 1, s.Rational(5, 18): 1}
    if covariance.eigenvals() != full_spectrum:
        raise ArithmeticError("the exact port covariance spectrum changed")
    values = sorted(compressed.diagonal(), key=lambda value: float(value))
    if not all(value.is_positive is True for value in values):
        raise ArithmeticError("the family covariance must be strictly positive")
    spread = s.simplify(values[-1]/values[0])
    if (s.Rational(5, 4)-spread).is_positive is not True:
        raise ArithmeticError("the exact covariance spread bound changed")
    return {
        "schema": "oph.hierarchy_diagnosis.native_covariance.v1",
        "source": {"binary_fixed_total": 2, "states": 66, "ports": 12,
                   "protected_core": list(local_preparation.CORE),
                   "ideal_faithful_reference_used": True,
                   "finite_preparation_error_propagated": False,
                   "empirical_mass_or_coupling_inputs": []},
        "definition": "C_ij=E_sigma[n_i n_j]-E_sigma[n_i]E_sigma[n_j]",
        "occupation_means": [str(value) for value in means],
        "second_moment": _entries(second),
        "port_covariance": _entries(covariance),
        "port_covariance_spectrum": [{"eigenvalue": str(value), "multiplicity": multiplicity}
                                     for value, multiplicity in sorted(full_spectrum.items())],
        "port_rank": 11,
        "port_kernel": "span of the all-ones vector; total occupation is exactly two",
        "family_compression": {
            "definition": "C_family=(B^T B)^(-1) B^T C B",
            "basis_Gram": _entries(gram), "basis_column_sum_is_zero_exact": True,
            "canonical_family_normalization_leaves_matrix_unchanged": True,
            "matrix": _entries(compressed), "rank": 3,
            "ascending_eigenvalues": [str(value) for value in values],
            "ascending_eigenvalues_float_diagnostic": list(map(float, values)),
            "positive_definite_exact": True,
            "largest_over_smallest": str(spread),
            "largest_over_smallest_float_diagnostic": float(spread),
            "largest_over_smallest_less_than_5_over_4_exact": True},
        "comparison_to_existing_feature": {
            "F": _entries(feature), "commutator": _entries(feature*compressed-compressed*feature),
            "commutes_with_F_exact": True,
            "same_diagonal_source_coordinate_axes": True,
            "every_separate_spectral_function_of_F_and_C_family_commutes": True},
        "scope": {
            "connected_reference_statistic_computed": True,
            "new_coherent_interaction_or_mass_selected": False,
            "raw_covariance_linear_texture_has_large_family_hierarchy": False,
            "covariance_and_F_alone_supply_nontrivial_quark_mixing": False,
            "arbitrary_interacting_native_histories_excluded": False,
            "reason": "Both supplied family matrices are diagonal; a linear covariance texture has spread below5/4. Assigning covariance to a Yukawa or interaction channel remains a separate architecture choice."},
    }
