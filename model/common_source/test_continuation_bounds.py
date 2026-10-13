"""Independent references and negative controls for stagewise continuation."""
from dataclasses import dataclass
from fractions import Fraction as F
import json
import math

import mpmath as mp
import numpy as np
import pytest

from . import continuation_bounds as b
from . import scalar_coupling, smooth_refinement


CANONICAL = (("evolve", F(1)), ("gauge_read", 0), ("feedback", 0),
             ("evolve", F(1, 2)), ("smooth_read", 1),
             ("evolve", F(1, 2)), ("smooth_read", 1), ("feedback", 2),
             ("evolve", F(1, 2)), ("gauge_read", 0))


def _mp(value):
    value = F(value)
    return mp.mpf(value.numerator) / value.denominator


def _reference(d, stages):
    """Unrounded analytic integrals, independent of certificate helpers."""
    amplitude = mp.sqrt(mp.mpf(143748821) / 8192000)
    pair = d // 2
    tail = mp.mpf(81) / 80 * math.comb(2 * pair, pair) / (mp.mpf(4)**pair * 81**pair)
    projected = mp.sqrt(2 * tail)
    for kind, value in stages:
        if kind == "evolve":
            growth = mp.exp(mp.mpf(13) * _mp(value) / 8)
            projected += amplitude * (growth - 1) * mp.mpf(8) / 13 * (
                mp.mpf(9) / (32 * (d + 2)) + mp.mpf(1) / (4 * mp.mpf(d + 2)**mp.mpf("1.5")))
            amplitude *= growth
        elif kind == "smooth_read":
            projected += amplitude * (mp.e - 1) / (2 * mp.mpf(d + 2)**mp.mpf("1.5"))
            amplitude *= mp.e
    outside = amplitude / (d + 2)**2
    return projected, outside, mp.sqrt(projected**2 + outside**2), amplitude**2


@pytest.mark.parametrize("d", (2, 4, 16, 128, 1024))
def test_same_history_repeated_reads_and_later_clock_dominate_unrounded_integrals(d):
    certificate = b.occupation_bound(d, CANONICAL)
    with mp.workdps(100):
        expected = _reference(d, CANONICAL)
        keys = ("final_projected_history_vector_error_upper", "final_outside_vector_norm_upper",
                "complete_stinespring_vector_error_upper", "final_fourth_number_moment_upper")
        for key, target in zip(keys, expected):
            actual = _mp(certificate[key])
            assert target <= actual
            assert actual - target < mp.mpf("1e-22")
    assert certificate["total_hamiltonian_model_time"] == "5/2"
    assert certificate["committed_binary_records"] == 4
    assert certificate["smooth_pointer_count"] == 2
    assert len(certificate["stages"]) == len(CANONICAL)
    json.dumps(certificate)


def test_prior_one_pointer_theorem_is_recovered_without_changing_its_scope():
    stages = CANONICAL[:5]
    for d in (4, 64, 256):
        extended = b.occupation_bound(d, stages)
        prior = smooth_refinement.smooth_occupation_bound(d)
        assert abs(F(extended["complete_stinespring_vector_error_upper"])
                   - F(prior["complete_stinespring_vector_error_upper"])) < F(1, 10**25)
    assert not smooth_refinement.status()["old_sharp_readout_gap_closed_by_this_result"]


def test_checkpoint_recursion_carries_post_pointer_moment_into_later_evolution():
    cert = b.occupation_bound(64, CANONICAL)
    rows = cert["stages"]
    for old, new in zip(rows, rows[1:]):
        assert old["outgoing_sqrt_fourth_moment_upper"] == new["incoming_sqrt_fourth_moment_upper"]
        assert F(new["projected_history_vector_error_upper"]) >= F(old["projected_history_vector_error_upper"])
    # Equal-duration later segments cost more after the actual smooth reads.
    assert F(rows[8]["projected_duhamel_increment_upper"]) > F(rows[5]["projected_duhamel_increment_upper"]) > F(rows[3]["projected_duhamel_increment_upper"])
    assert F(cert["final_projected_history_vector_error_upper"]) == (
        F(cert["preparation_vector_error_upper"])
        + sum(F(row["projected_duhamel_increment_upper"]) for row in rows))
    # Only the final orthogonal tail is added, not a tail per checkpoint.
    squared = (F(cert["final_projected_history_vector_error_upper"])**2
               + F(cert["final_outside_vector_norm_upper"])**2)
    assert F(cert["complete_stinespring_vector_error_upper"])**2 >= squared
    neutral = b.occupation_bound(64, (("gauge_read", 0), ("feedback", 0), ("gauge_read", 0)))
    empty = b.occupation_bound(64, ())
    for key in ("final_projected_history_vector_error_upper", "final_fourth_number_moment_upper",
                "complete_stinespring_vector_error_upper"):
        assert neutral[key] == empty[key]


@pytest.mark.parametrize("d", (2, 4, 6))
def test_pointer_fourth_moment_generator_inequality_on_the_whole_finite_carrier(d):
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    q = (a + a.T) / np.sqrt(2)
    generator = np.kron(np.diag([1., -1.]), np.kron(np.eye(d), q)) / 2
    weight = np.tile([(n0 + n1 + 2)**4 for n1 in range(d) for n0 in range(d)], 2)
    normalized = 1j * generator * (weight[None, :] - weight[:, None])
    normalized /= np.sqrt(weight[:, None] * weight[None, :])
    assert max(abs(np.linalg.eigvalsh(normalized))) < 2


def test_multiple_smooth_reads_keep_summed_branch_moment_but_not_each_conditional_moment():
    from scipy.linalg import expm
    d = 4
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    q = (a + a.T) / np.sqrt(2)
    minus, plus = expm(-.5j * q), expm(.5j * q)
    kraus = ((minus - 1j * plus) / 2, (minus + 1j * plus) / 2)
    weight = np.diag((np.arange(d) + 2)**4)
    rng = np.random.default_rng(713)
    initial = rng.normal(size=d) + 1j * rng.normal(size=d)
    initial /= np.linalg.norm(initial)
    branches = {(): initial}
    for _ in range(3):
        before = sum(np.vdot(v, weight @ v).real for v in branches.values())
        branches = {record + (bit,): k @ vector for record, vector in branches.items()
                    for bit, k in enumerate(kraus)}
        after = sum(np.vdot(v, weight @ v).real for v in branches.values())
        assert after <= math.exp(2) * before
        assert sum(np.vdot(v, v).real for v in branches.values()) == pytest.approx(1, abs=2e-14)
    assert len(branches) == 8
    # Dividing every branch by its probability destroys the ensemble contract.
    incorrectly_normalized = sum(np.vdot(v, weight @ v).real / np.vdot(v, v).real for v in branches.values())
    assert incorrectly_normalized > 2 * after


def test_fixed_protocol_cutoff_limit_and_exact_temporal_scaling():
    values = [F(b.occupation_bound(d, CANONICAL)["complete_stinespring_vector_error_upper"])
              for d in (64, 256, 1024, 4096)]
    assert all(a > c for a, c in zip(values, values[1:]))
    assert values[-1] < F(1, 10)
    for d in (2, 4):
        coarse = b.complete_error_budget(d, CANONICAL, steps=64)
        fine = b.complete_error_budget(d, CANONICAL, steps=128)
        assert len(coarse["ideal_temporal_segment_certificates"]) == 4
        for first, second in zip(coarse["ideal_temporal_segment_certificates"], fine["ideal_temporal_segment_certificates"]):
            assert F(first["unclipped_operator_error_upper"]) == 4 * F(second["unclipped_operator_error_upper"])
            source = scalar_coupling.time_error_bound(d, 64, F(first["model_time_exact"]))
            assert first["unclipped_operator_error_upper"] == source["unclipped_operator_error_upper"]
        assert not coarse["floating_execution_or_physical_service_error_included"]
        assert not coarse["nonvacuous_combined_bound"]
    analytic = b.complete_error_budget(64, CANONICAL, steps=64**2)
    assert not analytic["source_compiler_available_at_this_dimension"]


@dataclass(frozen=True)
class _Operation:
    kind: str
    time: F = F(0)
    record: int | None = None
    enabled: bool = True


def test_runtime_attribute_interface_and_disabled_feedback_are_explicit():
    operations = (_Operation("evolve", F(2)), _Operation("smooth_read"),
                  _Operation("feedback", record=0, enabled=False), _Operation("gauge_read"))
    assert b.validated_stages(operations) == (("evolve", F(2)), ("smooth_read", F(1)),
                                            ("feedback", 0), ("gauge_read", F(0)))


def test_actual_runtime_program_and_sharp_alternative_have_distinct_certificate_scope():
    from . import protocol
    program = protocol.continuation_program()
    certificate = b.program_bound(program)
    assert certificate["occupation_certificate"]["total_hamiltonian_model_time"] == "5/2"
    assert certificate["occupation_certificate"]["committed_binary_records"] == 4
    assert b.validated_stages(program.operations) == CANONICAL
    with pytest.raises(ValueError, match="sharp"):
        b.program_bound(protocol.sharp_program())


@pytest.mark.parametrize("stages", (
    (("sharp_read", 0),), (("evolve", .5),), (("evolve", True),),
    (("evolve", -1),), (("smooth_read", 0),), (("gauge_read", 1),),
    (("feedback", 0),), (("smooth_read", 1), ("feedback", 1)),
    (("smooth_read", 1), ("feedback", True)), (("reset", 0),),
    (_Operation("gauge_read", F(1)),), (_Operation("evolve", record=0),),
    (_Operation("evolve", enabled=False),),
))
def test_rejects_unsupported_semantics_instead_of_silently_certifying(stages):
    with pytest.raises(ValueError):
        b.occupation_bound(4, stages)


@pytest.mark.parametrize("d,coupling", ((True, F(1, 4)), (3, F(1, 4)), (4, .25), (4, F(1, 8))))
def test_rejects_inputs_outside_declared_occupation_family(d, coupling):
    with pytest.raises(ValueError):
        b.occupation_bound(d, CANONICAL, coupling=coupling)


def test_service_budget_never_converts_missing_or_partial_premises_into_zero_error():
    keys = b.required_service_keys((("smooth_read", 1),))
    absent = b.service_error_budget(F(1, 10), required_stage_keys=keys)
    assert absent["conditional_complete_error_upper"] is None
    assert not absent["service_bounds_supplied"]
    errors = {"preparation": F(1, 100), "operation_0:smooth_read": F(1, 50), "storage": F(1, 1000)}
    partial = b.service_error_budget(F(1, 10), {"preparation": F(1, 100)}, required_stage_keys=keys)
    assert partial["conditional_complete_error_upper"] is None
    full = b.service_error_budget(F(1, 10), errors, required_stage_keys=keys)
    assert F(full["conditional_complete_error_upper"]) == F(131, 1000)
    assert not full["service_bounds_independently_established_by_this_function"]
    assert b.service_error_budget(F(4, 5), dict.fromkeys(keys, F(1, 2)), required_stage_keys=keys)["conditional_complete_error_upper"] == "1"
    for errors in ({"preparation": -.1}, {"preparation": .001}, {"preparation": True}, {"": F(1, 10)}, {"unknown_stage": F(0)}):
        with pytest.raises(ValueError):
            b.service_error_budget(F(1, 10), errors, required_stage_keys=keys)
    empty = b.service_error_budget(F(1, 10), {}, required_stage_keys=keys)
    assert empty["conditional_complete_error_upper"] is None
    assert empty["missing_service_stage_keys"] == list(keys)
    complete_keys = b.required_service_keys(CANONICAL)
    assert len(complete_keys) == len(CANONICAL) + 2
    assert complete_keys[9] == "operation_8:evolve"
