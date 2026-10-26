"""Independent occupation, projection and basis checks of native anisotropy."""
from fractions import Fraction as F
import hashlib
import json
import math

import numpy as np
import pytest
import sympy as s

from . import family_attachment as f
from .. import slow_geometry
from ...native_modular_reference import local_preparation
from ...observer import native


@pytest.fixture(scope="module")
def packet():
    return f.report()


def _float(pair):
    return float(F(pair[0])) + math.sqrt(5)*float(F(pair[1]))


def test_first_moments_follow_actual_sixty_six_state_preparation(packet):
    states, probabilities = native.states(binary=True), local_preparation.reference()
    sector_masses = [sum((p for x,p in zip(states,probabilities,strict=True)
                          if sum(x[i] for i in local_preparation.CORE)==q), F()) for q in range(3)]
    assert sector_masses == [F(2,5), F(5,9), F(2,45)]
    core_mean = (sector_masses[1]+2*sector_masses[2])/2
    outside_mean = (2-2*core_mean)/10
    assert core_mean == F(29,90) and outside_mean == F(61,450)
    assert packet["source"]["occupations"] == [str(core_mean)]*2+[str(outside_mean)]*10


def test_three_positive_feature_eigenvalues_match_full_twelve_port_compression(packet):
    c = slow_geometry.source()
    p = c.q5_matrix_to_float(c.slow_band_projector_exact())
    means = np.array(list(map(float, map(F,packet["source"]["occupations"]))))
    values = np.linalg.eigvalsh(p @ np.diag(means) @ p)
    expected = [61/450, 41/225-7*math.sqrt(5)/750, 41/225+7*math.sqrt(5)/750]
    np.testing.assert_allclose(values[:9],0,atol=2e-16)
    np.testing.assert_allclose(values[9:],expected,rtol=0,atol=2e-16)
    np.testing.assert_allclose(list(map(_float,packet["compression"]["eigenvalues_Qsqrt5"])),
                               values[9:],rtol=0,atol=2e-16)
    assert np.ptp(values[9:]) > .06


def test_explicit_chart_binding_really_transports_edges_and_projector(packet):
    c = slow_geometry.source()
    chart = packet["family_band"]["native_to_retained_family_chart"]
    carrier = json.loads(f.MANIFEST.read_text())["carrier"]
    expected_edges = {tuple(sorted(int(p[1:]) for p in edge)) for edge in carrier["edges"]}
    assert {tuple(sorted((chart[i],chart[j]))) for i,j in c.seams()} == expected_edges
    retained = json.loads(f.FAMILY.read_text())["source_band_projector"]
    r = np.array([[float(s.sympify(x)) for x in row] for row in retained])
    p = c.q5_matrix_to_float(c.slow_band_projector_exact())
    np.testing.assert_allclose(p,r[np.ix_(chart,chart)],rtol=0,atol=3e-17)
    # Equal numeric port labels are not this chart binding.
    assert not np.array_equal(p,r)


def test_compression_eigenvalues_survive_nonorthogonal_family_coordinates(packet):
    c = slow_geometry.source()
    b = c.q5_matrix_to_float(slow_geometry.raw_frame())
    means = np.array(list(map(float,map(F,packet["source"]["occupations"]))))
    transform = np.array([[1.,2.,0.],[0.,1.,1.],[1.,0.,2.]])
    transformed = b @ transform
    metric = transformed.T @ transformed
    operator = np.linalg.solve(metric,transformed.T @ np.diag(means) @ transformed)
    np.testing.assert_allclose(metric@operator,operator.T@metric,rtol=0,atol=2e-14)
    expected = packet["compression"]["eigenvalue_floats_diagnostic"]
    np.testing.assert_allclose(np.sort(np.linalg.eigvals(operator)),expected,rtol=0,atol=3e-16)


def test_covariant_background_is_not_an_unbroken_family_intertwiner(packet):
    symmetry = packet["symmetry"]
    assert symmetry["proper_rotations_checked"] == 60
    assert symmetry["all_covariance_identities_exact"]
    assert symmetry["reference_stabilizer_size"] == 2
    assert symmetry["feature_stabilizer_size"] == 4
    assert not symmetry["full_A5_invariance"]
    assert not symmetry["unbroken_family_naturality_theorem_contradicted"]
    twirl = np.array([[ _float(x) for x in row] for row in symmetry["full_group_twirl_Qsqrt5"]])
    np.testing.assert_allclose(twirl,np.eye(3)/6,rtol=0,atol=0)
    # Complete uniform preparation has only the scalar invariant feature.
    c = slow_geometry.source()
    b = c.q5_matrix_to_float(slow_geometry.raw_frame())
    uniform = np.linalg.solve(b.T@b,b.T@(np.eye(12)/6)@b)
    np.testing.assert_allclose(uniform,np.eye(3)/6,rtol=0,atol=1e-16)


def test_conditional_yukawa_control_is_neither_empirical_fit_nor_particle_identification(packet):
    singular = list(map(_float,packet["if_Yukawa_endomorphism_is_declared_equal_to_F"]["relative_singular_values_Qsqrt5"]))
    assert singular[0] == 1 and 1 < singular[1] < singular[2] < 1.5
    interpretation = packet["interpretation"]
    assert interpretation["source_covariant_symmetry_breaker_computed"]
    assert not interpretation["source_selected_physical_Yukawa_intertwiner"]
    assert not interpretation["charged_propagating_pole_representation_identified"]
    assert not interpretation["absolute_scale_or_measured_masses_derived"]
    assert packet["source"]["empirical_mass_or_coupling_inputs"] == []
    for name,digest in packet["source_bindings"].items():
        assert hashlib.sha256((f.ROOT/name).read_bytes()).hexdigest() == digest
    json.dumps(packet,allow_nan=False)
