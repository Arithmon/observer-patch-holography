"""Complete mixed-state, reference and hostile-input checks for local clocks."""
from fractions import Fraction as F
from copy import deepcopy
import math

import numpy as np
import pytest

from . import local_history as h
from ..observer import native


def test_exact_pair_controls_fix_every_other_native_and_padding_address():
    codes = h.full.support()[0]
    states = h.common.regional.source()[0].states
    for index, (a, b) in enumerate(h.PAIRS):
        for sign in (-1, 1):
            actual, target = h.rotation(index, sign), h.ideal_rotation(index, sign)
            assert h.certify_rotation(index, sign) < h.ROTATION_GATE
            np.testing.assert_allclose(actual, target, atol=2e-12)
            wrong = actual.copy(); wrong[511, 511] += .125
            assert h.certify_rotation(index, sign, wrong) > F(1, 9)
        i, j = codes.index(a), codes.index(b)
        assert h.local.reference()[j]/h.local.reference()[i] == F(25, 8)
        assert all(states[i][p] == states[j][p] for p in range(12)
                   if p not in h.common.CAPS[index])
    assert h.exact_scope()["source_control_requires_exterior_address_selection"]


def test_complete_mixed_state_protocol_and_zero_clock_countercontrol():
    sigma = np.diag(list(map(float, h.local.reference()))).astype(complex)
    bound = h.source_budget()
    for law in ("C", "E"):
        raw, ideal = h.execute(sigma, law), h.execute(sigma, law, raw=False)
        assert h.distance(raw, ideal) < float(bound["complete_history_half_diamond_upper"])
        probabilities = h.probabilities(raw)
        assert abs(sum(probabilities.values())-1) < 1e-10
        assert min(probabilities.values()) >= -1e-14
        assert any(key[0] == "read" and key[3] == 0 for key in raw)
        assert h.distance(ideal, h.execute(sigma, law, raw=False, s=F())) > 1e-4
        # Comparing only separate marginals would lose the continuing work.
        assert all(v.shape == (66, 66) for key, v in raw.items() if key[0] == "read")


def _reference_marginal(outputs):
    result = np.zeros((2, 2), complex)
    for key, value in outputs.items():
        result += value if key[0] == "abort" else np.einsum("aras->rs", value.reshape(66, 2, 66, 2))
    return result


def _primitive_instrument(law, seam):
    """Group the native primitive directly, without the interface map builder."""
    states = tuple(native.states(binary=True))
    positions = {state: i for i, state in enumerate(states)}
    groups = {}
    for coin in (0, 1):
        for i, state in enumerate(states):
            out, (_, current) = native.native(state, seam, coin)
            label = current if law == "C" else state[seam[0]]-state[seam[1]]
            key = (coin, label, current)
            groups.setdefault(key, np.zeros((66, 66), complex))[positions[out], i] = 1/math.sqrt(2)
    return groups


def _independent_history(density, law, *, feedback=True):
    """Native-space spectral clock and literal sequential Kraus composition.

    This oracle uses neither production stage/rotation helpers nor exported
    interface matrices. The common factor60 in the source clock contributes
    only a global phase, which cancels from every density block.
    """
    codes = h.full.support()[0]
    reference = len(density)//66
    phase = np.diag(np.exp(-1j*np.log(list(map(float, h.local.reference())))))

    def stage(index, sign):
        a, b = (codes.index(code) for code in h.PAIRS[index])
        rotation = np.eye(66, dtype=complex)
        rotation[np.ix_([a, b], [a, b])] = np.array([[1, sign], [-sign, 1]])/math.sqrt(2)
        return np.kron(rotation@phase@rotation.conj().T, np.eye(reference))

    first = stage(0, 1)
    ready = first@density@first.conj().T
    result = {("abort", 0): np.zeros((reference, reference), complex)}
    for record, matrix in _primitive_instrument(law, h.SEAMS[0]).items():
        k = np.kron(matrix, np.eye(reference))
        middle = k@ready@k.conj().T
        second = stage(1, -1 if feedback and record[2] < 0 else 1)
        middle = second@middle@second.conj().T
        result[("abort", 1, *record)] = np.zeros((reference, reference), complex)
        for later, matrix in _primitive_instrument(law, h.SEAMS[1]).items():
            k = np.kron(matrix, np.eye(reference))
            result[("read", *record, *later)] = k@middle@k.conj().T
    return result


@pytest.mark.parametrize("law", ("C", "E"))
def test_complete_history_matches_independent_native_kraus_oracle_and_needs_feedback(law):
    rng = np.random.default_rng(1127)
    vector = rng.normal(size=132)+1j*rng.normal(size=132)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    original = density.copy()
    expected = _independent_history(density, law)
    for raw, tolerance in ((False, 2e-14), (True, 2e-11)):
        actual = h.execute(density, law, raw=raw)
        assert actual.keys() == expected.keys()
        for key in actual:
            np.testing.assert_allclose(actual[key], expected[key], rtol=0, atol=tolerance)
        np.testing.assert_allclose(_reference_marginal(actual),
                                   np.einsum("aras->rs", density.reshape(66, 2, 66, 2)),
                                   rtol=0, atol=tolerance)
    np.testing.assert_array_equal(density, original)
    # A complete second read with its feedback disabled is a different
    # quantum-reference history, despite having the same record alphabet.
    fixed_sign = _independent_history(density, law, feedback=False)
    assert fixed_sign.keys() == expected.keys()
    difference = sum(np.linalg.norm(actual[key]-fixed_sign[key], "nuc") for key in actual)/2
    assert difference > .005


@pytest.mark.parametrize("law", ("C", "E"))
def test_genuine_export_leak_retains_offdiagonal_external_reference(monkeypatch, law):
    codes = h.full.support()[0]
    a, b = h.PAIRS[0]
    outside = next(i for i in range(512) if i not in set(codes))
    work = np.zeros(512)
    work[a] = work[b] = 1/math.sqrt(2)
    padding = np.zeros(512)
    padding[outside] = 1
    # A deliberately uncertified stage swaps a native superposition with
    # one padding address. It forces an order-one export abort and tests the
    # instrument itself, independently of the small certified source error.
    leak = (np.eye(512)-np.outer(work, work)-np.outer(padding, padding)
            +np.outer(work, padding)+np.outer(padding, work))
    monkeypatch.setattr(h, "_stage", lambda index, sign, raw, s:
                        leak if index == 0 else np.eye(512))
    vector = np.zeros(132, complex)
    vector[2*codes.index(a)] = vector[2*codes.index(b)+1] = 1/math.sqrt(2)
    density = np.outer(vector, vector.conj())
    actual = h.execute(density, law, raw=False)
    # Erasing the work leaves a coherent |+><+|/2 block on R. A scalar
    # failure probability, reset reference, or discarded abort cannot pass.
    np.testing.assert_allclose(actual[("abort", 0)], np.ones((2, 2))/4,
                               rtol=0, atol=2e-14)
    assert np.trace(actual[("abort", 0)]).real > .49
    for key, value in actual.items():
        if key[:2] == ("abort", 1):
            np.testing.assert_array_equal(value, np.zeros((2, 2)))
    np.testing.assert_allclose(_reference_marginal(actual), np.eye(2)/2,
                               rtol=0, atol=2e-14)
    assert abs(sum(np.trace(value).real for value in actual.values())-1) < 2e-14


@pytest.mark.parametrize("law", ("C", "E"))
def test_entangled_continuation_and_terminal_abort_keep_external_reference(law):
    # The first pair is initially correlated with R, not an independent probe.
    codes = h.full.support()[0]
    v = np.zeros(132, complex)
    v[2*codes.index(328)] = .6
    v[2*codes.index(73)+1] = .8j
    density = np.outer(v, v.conj())
    ideal, raw = h.execute(density, law, raw=False), h.execute(density, law)
    expected = np.diag([.36, .64])
    np.testing.assert_allclose(_reference_marginal(ideal), expected, atol=2e-14)
    np.testing.assert_allclose(_reference_marginal(raw), expected, atol=2e-11)
    assert h.distance(ideal, raw) < float(h.source_budget()["complete_history_half_diamond_upper"])
    # The second cap is fed the disturbed first output, not the initial state.
    first, _ = h._read(density, 0, 1, law, False, h.PARAMETER)
    for record, branch in first.items():
        second, _ = h._read(branch, 1, -1 if record[2] < 0 else 1, law, False, h.PARAMETER)
        for later, value in second.items():
            np.testing.assert_allclose(ideal[("read", *record, *later)], value, atol=0)


def test_finite_preparation_failures_and_success_are_one_complete_history():
    budget = h.local.budget()
    source = h.source_budget()
    upper = float(source["channel_norm_upper"]*F(budget["complete_flagged_reference_halftrace_upper"])
                  +F(5, 33)*source["complete_history_half_diamond_upper"])
    for law in ("C", "E"):
        actual = h.prepared_history(law)
        assert all(("preparation", flag) in actual for flag in
                   ("herald_failure", "warmup_timeout", "completion_timeout"))
        assert abs(sum(h.probabilities(actual).values())-1) < 2e-10
        assert h.distance(actual, h.target_history(law)) < upper
        success = sum(p for key, p in h.probabilities(actual).items() if key[0] != "preparation")
        assert .15 < success < .152


@pytest.mark.parametrize("bad", (True, -1, 2, 0., "0"))
def test_pair_address_types_are_not_coerced(bad):
    with pytest.raises(ValueError):
        h.pair_word(bad, 1)


def test_density_screen_checks_aggregate_negative_mass_and_original_arrays():
    d = np.eye(66, dtype=complex)/66
    for bad in (d.astype(object), np.full((66, 66), np.nan), d*1.0001,
                np.eye(66)*1e308):
        with pytest.raises(ValueError):
            h.execute(bad, "C", raw=False)
    negative = np.diag([1+65e-14]+[-1e-14]*65).astype(complex)
    with pytest.raises(ValueError):
        h.execute(negative, "C", raw=False)
    wrong = d.copy(); wrong[1, 0] = 1e-8j
    with pytest.raises(ValueError):
        h.execute(wrong, "E", raw=False)


@pytest.mark.parametrize("bad", (np.nan, np.inf, 1e308))
def test_distance_cannot_hide_nonfinite_or_overflowing_outputs(bad):
    with pytest.raises(ValueError):
        h.distance({"read": np.array([[bad]])}, {"read": np.array([[-bad]])})


def test_no_cached_type_alias_or_uncertified_parameter_enters_raw_history():
    density = np.eye(66, dtype=complex)/66
    for bad in (1, 0., "raw"):
        with pytest.raises(ValueError):
            h.execute(density, "C", raw=bad)
    with pytest.raises(ValueError, match="retained parameter"):
        h.execute(density, "C", s=F(2))


def test_retained_source_costs_include_every_feedback_rotation():
    budget = h.source_budget()
    pulses = sum(2*len(h.pair_word(i, 1))+len(h.word()) for i in (0, 1))
    assert budget["maximum_source_pulses"] == pulses
    assert budget["maximum_proper_code_service_events"] == 2*pulses
    assert budget["continuing_native_read_attempts"] == 2
    assert not budget["selected_configuration_initialization"]
    assert budget["projective_work_bit_readouts"] == 0


def test_receipt_rejects_negative_diagnostics_and_inconsistent_mass_inside_float_tolerance():
    receipt = h.strict_json.loads(h.full.preparation.canonical(h.run()))
    h.compare_receipt(receipt, receipt)
    negative = deepcopy(receipt)
    negative["histories"][0]["complete_quantum_halftrace_error_diagnostic"] = -1e-10
    with pytest.raises((ValueError, AssertionError)):
        h.compare_receipt(negative, receipt)
    lost = deepcopy(receipt)
    lost["histories"][0]["complete_preparation_and_history_mass_diagnostic"] += 1e-11
    with pytest.raises((ValueError, AssertionError)):
        h.compare_receipt(lost, receipt)
    probability = deepcopy(receipt)
    zero = next(r for r in probability["histories"][0]["complete_output"] if r["probability"] < 1e-11)
    zero["probability"] = -1e-10
    with pytest.raises((ValueError, AssertionError)):
        h.compare_receipt(probability, receipt)
