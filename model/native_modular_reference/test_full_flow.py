"""Complete source, omitted-address, phase and reference controls."""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
import pytest
from mpmath import MPContext
from scipy.linalg import expm

from . import full_flow as f
from ..observer.source_scalar.test_compiler import independent_source_generators


def test_full_native_support_and_exact_gray_parity_walk():
    row = f.support_certificate()
    assert len(row["previously_omitted_degree8_addresses"]) == 22
    assert row["native_dimension"] == 66 and row["identity_extension_addresses"] == 446
    assert max(row["degree8_addresses"]) > 255
    word = f.word(F(1))
    assert len(word) == 2042
    assert sum(event["kind"] == "I" for event in word) == 1
    assert sum(event["kind"] == "Z" for event in word) == 511
    assert sum(event["kind"] == "J" for event in word) == 510
    assert max(max(event["targets"]) for event in word) == 8
    # Independently track every parity chart; every cyclic walk returns home.
    masks = {0}
    for pivot in range(9):
        parity = 1 << pivot
        for index in range(1 << pivot):
            masks.add(parity)
            following = (index+1) % (1 << pivot)
            gray = index ^ (index >> 1)
            next_gray = following ^ (following >> 1)
            parity ^= gray ^ next_gray
        assert parity == 1 << pivot
    assert masks == set(range(512))


@pytest.mark.parametrize("s", (F(0), F(1), F(-13, 7), F(16)))
def test_transcendental_target_encloses_independent_high_precision_phase(s):
    ctx = MPContext(); ctx.dps = 180
    def mp(x): return ctx.mpf(x.numerator)/x.denominator
    lo, hi = f.logarithm_interval()
    true_log = ctx.log(ctx.mpf(5)/4)
    assert mp(lo) < true_log < mp(hi)
    real, imag, error = f.phase_enclosure(s)
    actual = ctx.exp(ctx.j*mp(s)*true_log)
    assert abs(actual-(mp(real)+ctx.j*mp(imag))) <= mp(error)
    assert error < F(1, 10**45)


@lru_cache(maxsize=None)
def _address_groups(targets):
    mask = sum(1 << bit for bit in targets)
    return np.array([[base + sum(((entry >> (len(targets)-1-pos)) & 1) << bit
                                for pos, bit in enumerate(targets))
                      for entry in range(1 << len(targets))]
                     for base in range(512) if not base & mask])


def _independent_source_word(tape, vectors):
    generators = independent_source_generators()
    gates = {}
    out = vectors.copy()
    for event in tape:
        key = event["kind"], event["angle"]
        if key not in gates:
            coefficients = f.compiler.integrated_controls(*key)
            full = expm(np.einsum("p,pij->ij", coefficients, generators))
            size = 4 if key[0] == "J" else 2
            assert np.linalg.norm(full[size:, :size], 2) < 1e-13
            # The same complete M6 export explicitly keeps complement failure.
            exports = f.compiler.transfer_kraus(size)
            np.testing.assert_array_equal(sum(k.conj().T@k for k in exports), np.eye(6))
            gates[key] = full[:size, :size]
        groups = _address_groups(tuple(event["targets"]))
        out[groups] = np.einsum("ij,bjr->bir", gates[key], out[groups])
    return out


def test_complete_source_execution_matches_independent_geometric_oracle_with_reference():
    # Native and padding coherences share untouched reference columns.
    rng = np.random.default_rng(6151)
    probe = rng.normal(size=(512, 3))+1j*rng.normal(size=(512, 3))
    probe /= np.linalg.norm(probe)
    independent = _independent_source_word(f.word(F(1)), probe)
    source = f.matrix(F(1))@probe
    assert np.linalg.norm(source-independent) < 1e-12
    marked = list(f.support()[1])
    ideal = probe.copy()
    ideal[marked] *= np.exp(1j*math.log(5/4))
    assert np.linalg.norm(source-ideal) < float(f.WORD_ERROR_GATE)
    assert np.linalg.norm(source.conj().T@source-probe.conj().T@probe) < 2e-12
    # Cached execution cannot be altered into a different claimed source word.
    with pytest.raises(ValueError):
        f.matrix(F(1))[0, 0] = 0


def test_omitted_native_coherences_are_measurable_in_the_executed_Ramsey_protocol():
    row = f.witness(F(16))
    assert row["initial_native_code"] == 16 and row["coherent_partner_native_code"] == 144
    for law in ('C', 'E'):
        probabilities = row["executed_complete_probabilities_by_law"][law]
        assert len(probabilities) == 4
        assert 0 < probabilities[0] < 2e-8 and 0 < probabilities[1] < 2e-8
        assert probabilities[2] < .06 and probabilities[3] > .94
        assert abs(sum(probabilities)-1) < 2e-11
    assert row['native_attempts_max_including_both_resets_and_pointer'] == 5357
    assert row['total_source_pulses'] == 2050
    assert row['proper_code_service_events'] == 4100
    assert F(row['complete_native_reset_half_diamond_upper']) == F(140, 2**32)
    assert F(row['whole_Ramsey_half_diamond_upper']) < F(34, 10**9)
    assert row["outside_native_support_probability_diagnostic"] >= 0
    assert not row["declared_projective_register_bit_readout"]
    assert row["native_initialization_and_readout_attempts_charged"]
    assert row["private_reset_tapes_discarded"]
    assert not row["physical_service_error_authenticated"]
    assert 16 not in f.corner.reference()[3] and 144 not in f.corner.reference()[3]
    # The old twenty-address operation is the identity on both witness codes.
    assert f.corner.ideal_local(1)[16, 16] == 1
    assert f.corner.ideal_local(1)[144 % 128, 144 % 128] == 1


def test_complete_history_rejects_a_relabelled_or_suppressed_native_pointer(monkeypatch):
    original = f.pointer.read_bit
    def swapped(*args, **kwargs):
        return tuple(reversed(original(*args, **kwargs)))
    monkeypatch.setattr(f.pointer, 'read_bit', swapped)
    with pytest.raises(ArithmeticError, match='complete-word comparison'):
        f.witness(F(16))


def test_native_resets_cover_collision_inputs_and_both_code_targets():
    for target in (f.compiled.circuit.decode(16), f.pointer.ready_state()):
        probabilities = f.configuration.success_probabilities(full=True, mode='reset', target=target)
        bounds = f.configuration.budget(full=True, mode='reset', target=target)
        assert len(probabilities) == 78
        assert 0 < 1-min(probabilities) <= bounds['failure_upper']
        assert bounds['native_attempts_max'] == 2678
    # The six-pulse pointer has a norm slightly above one in its portable
    # envelope; dropping this factor is not a valid raw-array telescope.
    bound = f.pointer.certificate()
    assert bound['channel_norm_upper'] > 1
    assert bound['half_diamond_upper'] > 0


def test_raw_certificate_rejects_corner_only_padding_inverse_and_amplification():
    s = F(1)
    phase = np.exp(1j*math.log(5/4))
    correct = np.eye(512, dtype=complex)
    correct[list(f.support()[1]), list(f.support()[1])] = phase
    corner_only = np.eye(512, dtype=complex)
    old = list(f.corner.reference()[3]); corner_only[old, old] = phase
    assert f.certify(s, corner_only)["operator_error_upper"] > F(1, 2)
    wrong_padding = correct.copy(); wrong_padding[511, 511] += F(1, 8)
    assert f.certify(s, wrong_padding)["operator_error_upper"] >= F(1, 8)
    assert f.certify(s, correct.conj())["operator_error_upper"] > 1
    assert f.certify(s, -correct)["operator_error_upper"] > 1
    amplified = (1+2**-20)*correct
    assert f.certify(s, amplified)["operator_error_upper"] > F(1, 2**16)


def test_intermediate_native_projection_is_not_an_admitted_shortcut():
    native = set(f.support()[0])
    probe = np.zeros(512, complex); probe[16] = 1
    # The first CNOT has H pulses that use nonnative temporary addresses.
    # A hypothetical native-shell projection loses a positive amount of norm.
    for event in f.word(F(1)):
        probe = f.compiler.execute([event], 9, probe)
        outside = sum(abs(value)**2 for code, value in enumerate(probe) if code not in native)
        if outside > F(1, 4):
            break
    else:
        raise AssertionError("the control no longer exposes temporary nonnative support")
    assert sum(abs(probe[code])**2 for code in native) < .8


def test_prepared_modular_composition_keeps_failure_cost_and_conditional_domain():
    s, attempts = F(16), 369
    row = f.composed_bound(s, attempts)
    prepared = f.preparation.budget(attempts, s)
    source = f.checked_bound(s)
    assert row["source_to_prepared_modular_half_diamond_upper"] == (
        source["half_diamond_upper"]+prepared["modular_channel_half_diamond_error_upper"])
    assert row["source_to_prepared_modular_half_diamond_upper"] < F(971, 10**6)
    assert row["complete_flagged_preparation_and_source_half_trace_upper"] == (
        prepared["warmup_half_trace_error_upper"]
        +prepared["herald_success_interval"][1]*source["half_diamond_upper"])
    assert row["expected_native_preparation_attempts_upper"] > 2442
    with pytest.raises(ValueError, match="faithful"):
        f.composed_bound(F(1), 0)


def test_source_channel_against_noncommuting_herald_reference_with_entangled_probe():
    # Adversarial state inside the certified warm-up ball; no claim that this
    # selected perturbation is a retained native-history output.
    degree, sigma = f.preparation.anchor()
    columns = f.preparation.moving_columns("C")
    i, j = next((i, j) for (i, j), outputs in sorted(columns.items())
                if i < j and any(degree[a] != degree[b] for a, b in outputs))
    perturbation = F(1, 10**8)
    assert perturbation < f.preparation.budget(369)["warmup_half_trace_error_upper"]
    rho = np.eye(66)/66
    rho[i, j] = rho[j, i] = float(perturbation)
    heralded = np.zeros((66, 66))
    for (a, b), outputs in columns.items():
        for (x, y), coefficient in outputs.items():
            heralded[x, y] += float(coefficient)*rho[a, b]
    probability = np.trace(heralded)
    state = heralded/probability
    target = np.diag([float(x) for x in sigma])
    assert np.linalg.norm(state@target-target@state) > 1e-14
    values, vectors = np.linalg.eigh(state)
    # Cancel precisely the source target's native scalar gauge and the
    # preparation theorem's normalization scalar, never reverse the dual sign.
    centered = np.log(values)+np.log(probability/float(f.preparation.Q))+np.log(60)
    state_unitary = (vectors*np.exp(-16j*centered))@vectors.T
    probe = np.zeros((66, 2), complex)
    probe[i, 0] = probe[j, 1] = 1/np.sqrt(2)
    source_output = f.matrix(F(16))[:, list(f.support()[0])]@probe
    actual_output = np.zeros((512, 2), complex)
    actual_output[list(f.support()[0])] = state_unitary@probe
    # QR evaluates the complete two-vector signed density, including padding,
    # without subtractive cancellation of almost equal squared overlaps.
    _, reduced = np.linalg.qr(np.column_stack((source_output.ravel(), actual_output.ravel())))
    difference = reduced@np.diag([1, -1])@reduced.conj().T
    half_trace = np.abs(np.linalg.eigvalsh(difference)).sum()/2
    bound = f.composed_bound(F(16))["source_to_prepared_modular_half_diamond_upper"]
    assert half_trace < float(bound)


@pytest.mark.parametrize("s", (True, False, 1., "1", None, F(1601, 100), F(-17)))
def test_parameter_rejects_ambiguous_or_undeclared_values(s):
    with pytest.raises(ValueError):
        f.word(s)


def test_hostile_arrays_rejected_before_dispatch_and_unknown_mixer_refused():
    class Hostile(np.ndarray):
        def __getattribute__(self, name):
            if name in ("shape", "dtype", "copy"):
                raise AssertionError("untrusted array dispatch")
            return super().__getattribute__(name)
    with pytest.raises(ValueError, match="base ndarray"):
        f.certify(F(1), np.eye(512).view(Hostile))
    with pytest.raises(ValueError):
        f.certify(F(1), np.full((512, 512), np.nan))
    with pytest.raises(ValueError):
        f.mixer_word(True)
    with pytest.raises(ValueError):
        f.mixer_bound(True)
