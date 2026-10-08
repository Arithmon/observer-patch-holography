"""Original scalar validation precedes kinetic array coercion."""
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_interacting_quantum as quantum
import whitney_quantum_packet as packet


@pytest.mark.parametrize("field", ["edge", "matter", "charge", "mass_squared", "quartic", "state"])
@pytest.mark.parametrize("bad", [True, "1", 2**54 + 1, Fraction(1, 10),
                               Decimal("0.1"), Fraction(1, 2**1075)])
def test_original_scalars_cannot_change_before_the_kinetic_calculation(field, bad):
    edges, matter, coordinates = [0]*42, [0]*13, [0]*56
    parameters = {}
    if field == "edge":
        edges[0] = bad
    elif field == "matter":
        matter[0] = bad
    elif field == "state":
        coordinates[0] = bad
    else:
        parameters[field] = bad
    with pytest.raises(ValueError):
        if field == "state":
            quantum.gaussian_half_density_log(coordinates)
        else:
            quantum.coefficients(edges, matter, **parameters)


@pytest.mark.parametrize("field", ["edge", "matter", "state"])
def test_masked_entries_do_not_turn_into_supplied_coefficients(field):
    edge = np.ma.array(np.zeros(42), mask=False)
    matter = np.ma.array(np.zeros(13), mask=False)
    state = np.ma.array(np.zeros(56), mask=False)
    {"edge": edge, "matter": matter, "state": state}[field].mask[0] = True
    with pytest.raises(ValueError):
        if field == "state":
            quantum.gaussian_half_density_log(state)
        else:
            quantum.reduced_coefficients(edge if field == "edge" else edge.data,
                                         matter if field == "matter" else matter.data)


@pytest.mark.parametrize("value", [1, 1., np.int64(1), Fraction(1), Decimal(1)])
def test_exactly_representable_original_units_keep_the_same_metric(value):
    plain = quantum.reduced_kinetic(np.zeros(42), np.ones(13), charge=.25)
    supplied = quantum.reduced_kinetic([0]*42, [value]*13, charge=Fraction(1, 4))
    np.testing.assert_array_equal(supplied.factor, plain.factor)
    assert supplied.logdet() == plain.logdet()


@pytest.mark.parametrize("value", [2**54, Fraction(1, 4), Decimal("0.25"),
                                  Fraction(1, 2**1074)])
def test_nearby_representable_configuration_values_survive_validation(value):
    edges, matter = [0]*42, [0]*13
    edges[0], matter[1] = value, value
    result_edges, result_matter = quantum._configuration(edges, matter)
    assert Fraction(float(result_edges[0])) == Fraction(value)
    assert Fraction(float(result_matter[1].real)) == Fraction(value)


def test_wider_real_precision_is_not_silently_narrowed():
    if np.finfo(np.longdouble).nmant <= np.finfo(float).nmant:
        pytest.skip("longdouble has no precision beyond binary64 on this platform")
    supplied = np.longdouble(1)+np.finfo(np.longdouble).eps
    with pytest.raises(ValueError, match="input precision"):
        quantum.reduced_kinetic([0]*42, [supplied]+[0]*12)


@pytest.mark.parametrize("phase", [1, 1j])
def test_constant_field_has_no_spurious_gradient_energy(phase):
    # Partition of unity makes this field exactly constant. The baseline
    # leaves a ~6.6e63 gradient residual and ~9.6e128 false positive energy.
    reduction = quantum.reduced_kinetic(np.zeros(42), np.full(13, phase*1e80),
        charge=0, mass_squared=0, quartic=0)
    assert reduction.potential == 0
    _, energy = quantum.coefficients(np.zeros(42), np.full(13, phase*1e80),
        charge=0, mass_squared=0, quartic=0)
    assert energy == 0


@pytest.mark.parametrize("phase", [1, 1j])
def test_nearby_nonconstant_field_retains_actual_gradient_energy(phase):
    matter = np.full(13, complex(1e15))
    matter[0] += 1
    # Integral |grad(lambda_center)|² = 30-10sqrt(5) on the original cone.
    _, energy = quantum.coefficients(np.zeros(42), phase*matter,
                                     charge=0, mass_squared=0, quartic=0)
    assert energy == pytest.approx(30-10*np.sqrt(5), rel=2e-13, abs=0)


@pytest.mark.parametrize("evaluate", [quantum.coefficients, quantum.reduced_kinetic])
def test_required_overflowing_potential_is_refused(evaluate):
    with pytest.raises(ValueError, match="precision"):
        evaluate(np.zeros(42), np.full(13, 1e80), charge=0)


@pytest.mark.parametrize("delta", [0., 1e-100])
def test_public_cotangent_refuses_unresolved_gauge_tangent_cancellation(delta):
    mesh = quantum.geometry()
    gauge = np.r_[0., 1., -1., np.zeros(10)]
    velocity = np.r_[mesh.d@gauge, np.zeros(13), gauge/4]
    velocity[55] = delta
    q = np.r_[np.zeros(42), np.ones(13), np.zeros(13)]
    # The exact horizontal tangent is only the supplied center entry delta.
    # Binary64 gauge solving leaves O(1e-16) residuals instead. A well-resolved
    # reduced metric cannot certify this already damaged tangent.
    with pytest.raises(ValueError, match="kinetic precision.*tangent"):
        packet.phase_space(q, velocity, mesh)


@pytest.mark.parametrize("delta", [0., 1e-100, 1e-200])
def test_public_cotangent_preserves_resolved_small_tangent(delta):
    q = np.r_[np.zeros(42), np.ones(13), np.zeros(13)]
    velocity = np.zeros(68)
    velocity[55] = delta
    actual = packet.phase_space(q, velocity)
    expected = np.zeros(68)
    expected[55] = delta
    np.testing.assert_array_equal(actual['full_velocity'], expected)
    if delta:
        # The map is linear in this fixed-configuration tangent, including
        # when a squared-norm diagnostic would underflow to zero.
        unit = velocity/delta
        control = packet.phase_space(q, unit)
        np.testing.assert_allclose(actual['momentum']/delta, control['momentum'],
                                   rtol=2e-12, atol=1e-15)
