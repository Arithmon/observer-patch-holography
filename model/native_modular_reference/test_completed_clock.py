"""Same-state entropy/clock join, exact cocycle and raw source controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from mpmath import MPContext

from . import completed_clock as c
from . import regional as r


def test_native_completion_has_one_source_state_and_exact_bare_cap_restriction():
    sigma = c.reference()
    assert set(sigma) == {F(5, 336), F(19, 1200), F(37, 2400), F(1, 75), F(1, 60)}
    palm = r.source()[2]
    for sector in c.cut.cut_sectors(c.CAP):
        for row in sector.indices:
            assert sum(sigma[i] for i in row) == sum(palm[i] for i in row)
        ratios = {sigma[a[0]]/sigma[b[0]] for a in sector.indices for b in sector.indices}
        expected = ({F(1)}, {F(1), F(38, 37), F(37, 38)}, {F(1), F(5, 4), F(4, 5)})[sector.charge]
        assert ratios == expected
    # Each actual C/E cut attempt fixes every cap matrix unit, including
    # off-diagonal probes. This is stronger than preserving occupation means.
    for law in ("C", "E"):
        counts, scale = c.cut.native_cut_counts(c.CAP, law)
        for unit in r.basis_supports(c.CAP, completed=False):
            vector = np.zeros(4356, np.int64)
            for i, j in unit:
                vector[66*i+j] = 1
            np.testing.assert_array_equal(counts.T@vector, scale*vector)


def test_exact_cocycle_ratios_and_complete_walsh_reconstruction():
    clock = c.ratios("clock")
    cocycle = c.ratios("cocycle")
    codes = c.full.support()[0]
    for i, code in enumerate(codes):
        assert cocycle[code]*60*r.source()[2][i] == clock[code]
    for code in set(range(512))-set(codes):
        assert clock[code] == cocycle[code] == 1
    for kind in ("clock", "cocycle"):
        values = c.ratios(kind)
        for ratio, counts in c.walsh_counts(kind):
            for code in range(512):
                value = sum(n*(-1)**((code & mask).bit_count()) for mask, n in enumerate(counts))
                assert value == 512*int(values[code] == ratio)
        tape = c.word(kind, F(16))
        assert len(tape) == 2042 and sum(e["kind"] == "J" for e in tape) == 510


def test_rational_log_and_phase_bounds_against_independent_high_precision():
    ctx = MPContext(); ctx.dps = 150
    mp = lambda v: ctx.mpf(v.numerator)/v.denominator
    for ratio in sorted(set(c.ratios("clock")) | set(c.ratios("cocycle")) | {F(1, 2), F(2)}):
        lo, hi = c.logarithm_interval(ratio)
        assert mp(lo) <= ctx.log(mp(ratio)) <= mp(hi)
        for s in (F(-13, 7), F(0), F(16)):
            real, imag, error = c.phase_enclosure(s, ratio)
            actual = ctx.exp(-ctx.j*mp(s)*ctx.log(mp(ratio)))
            assert abs(actual-(mp(real)+ctx.j*mp(imag))) <= mp(error)


def test_raw_clock_and_cocycle_execute_on_full_register_without_hidden_projection():
    s = F(16)
    clock = c.matrix("clock", s)
    cocycle = c.matrix("cocycle", s)
    # Fresh certificates compare every original binary entry, not just probes.
    for kind, array in (("clock", clock), ("cocycle", cocycle)):
        assert c.certify(kind, s, array)["operator_error_upper"] < c.WORD_ERROR_GATE
    palm = c.full.matrix(s)
    c.full.checked_bound(s)
    assert np.linalg.norm(cocycle@palm-clock, "fro") < 1e-9
    rng = np.random.default_rng(472)
    probe = rng.normal(size=(512, 3))+1j*rng.normal(size=(512, 3))
    probe /= np.linalg.norm(probe)
    target = np.array([float(w)**(-1j*float(s)) for w in c.ratios("clock")])[:, None]*probe
    assert np.linalg.norm(clock@probe-target) < float(c.WORD_ERROR_GATE)
    # Original binary values are immutable, and padding has not been dropped.
    with pytest.raises(ValueError):
        clock[0, 0] = 0
    wrong = np.eye(512, dtype=complex)
    wrong[np.diag_indices(512)] = np.array([float(w)**(-1j*float(s)) for w in c.ratios("clock")])
    wrong[511, 511] += .125
    assert c.certify("clock", s, wrong)["operator_error_upper"] >= F(1, 8)
    assert c.certify("clock", s, clock.conj())["operator_error_upper"] > 1


def test_one_completed_word_gives_bare_cap_modular_flow_for_arbitrary_probe_reference():
    s = F(16)
    rng = np.random.default_rng(629)
    amplitude = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    amplitude /= np.linalg.norm(amplitude)
    encoded = np.zeros((512, 2), complex)
    codes = c.full.support()[0]
    encoded[list(codes)] = amplitude
    actual = c.matrix("clock", s)@encoded
    total_distance = 0.
    sigma = c.reference()
    for sector in c.cut.cut_sectors(c.CAP):
        out, target = [], []
        phases = np.array([float(sigma[row[0]])**(-1j*float(s)) for row in sector.indices])
        for j in range(len(sector.complement_states)):
            native_indices = [row[j] for row in sector.indices]
            encoded_indices = [codes[i] for i in native_indices]
            out.append(actual[encoded_indices].ravel())
            target.append((phases[:, None]*amplitude[native_indices]).ravel())
        rho_out = sum(np.outer(v, v.conj()) for v in out)
        rho_target = sum(np.outer(v, v.conj()) for v in target)
        total_distance += np.linalg.svd(rho_out-rho_target, compute_uv=False).sum()/2
    outside = sorted(set(range(512))-set(codes))
    failure_reference = actual[outside].T@actual[outside].conj()
    total_distance += np.linalg.svd(failure_reference, compute_uv=False).sum()/2
    assert total_distance < float(c.checked_bound("clock", s)["half_diamond_upper"])


def test_finite_preparation_budget_keeps_same_actual_state_flags_and_region():
    s = F(16)
    row = c.budget(s)
    palm = c.preparation.budget(369, s)
    cut = c.cut.cut_certificate(c.CAP, precision_bits=24)
    eta = F(cut["preparation"]["trace_distance_upper"])
    delta = palm["warmup_half_trace_error_upper"]
    t, u = 66*delta, 75*eta
    assert row["combined_reference_operator_sandwich_radius"] == t+u+t*u
    assert row["complement_native_attempts_after_success"] == 7440
    assert row["source_to_actual_completed_modular_half_diamond_upper"] < F(1, 1000)
    assert row["source_then_cap_restrict_vs_actual_cap_modular_half_diamond_upper"] < F(971, 10**6)
    assert row["complete_flagged_reference_half_trace_upper"] == delta+c.preparation.Q*eta
    assert row["expected_native_preparation_attempts_upper"] > 9882
    assert row["minimum_work_eigenvalue_lower"] > F(1333, 100000)
    compact = c.preparation._compact(row)
    for key, value in row.items():
        if isinstance(value, F) and key.endswith("_lower"):
            assert F(compact[key]) <= value
        if isinstance(value, F) and key.endswith("_upper"):
            assert F(compact[key]) >= value
    # Explicit adversarial noncommuting matrix sandwiches test the composition
    # algebra; actual native cut maps are separately checked on all units.
    rng = np.random.default_rng(128)
    target = np.diag(list(map(float, c.reference())))
    root = np.sqrt(target)
    h = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    h = h+h.conj().T; h /= np.linalg.norm(h, 2)
    v = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    v = v+v.conj().T; v /= np.linalg.norm(v, 2)
    tt, uu = .03, .02
    intermediate = root@(np.eye(66)+uu*h)@root
    eigenvalues, vectors = np.linalg.eigh(intermediate)
    sq = (vectors*np.sqrt(eigenvalues))@vectors.conj().T
    actual = sq@(np.eye(66)+tt*v)@sq
    radius = tt+uu+tt*uu
    assert np.linalg.eigvalsh(actual-(1-radius)*target).min() > 0
    assert np.linalg.eigvalsh((1+radius)*target-actual).min() > 0
    with pytest.raises(ValueError, match="faithful"):
        c.budget(F(1), attempts=0)


def test_actual_coherent_native_herald_and_cut_preserve_operator_sandwich():
    rng = np.random.default_rng(17)
    h = rng.normal(size=(66, 66))+1j*rng.normal(size=(66, 66))
    h = h+h.conj().T; h -= np.trace(h)*np.eye(66)/66
    h /= np.linalg.norm(h, 2)
    t = .03
    rho = (np.eye(66)+t*h)/66
    success = np.zeros((66, 66), complex)
    for (i, j), entries in c.preparation.moving_columns("C").items():
        for (a, b), value in entries.items():
            success[a, b] += float(value)*rho[i, j]
    no_counts = np.array(c.preparation.marked("C").no_current_counts)/60
    failure = no_counts*rho
    assert abs(np.trace(success)+np.trace(failure)-1) < 3e-16
    ideal = np.diag(list(map(float, r.source()[2])))
    counts, scale = c.cut.native_cut_counts(c.CAP, "C")
    actual = success/float(c.preparation.Q)
    for _ in range(3):
        actual = (counts@actual.ravel()).reshape(66, 66)/scale
        ideal = (counts@ideal.ravel()).reshape(66, 66)/scale
        assert np.linalg.eigvalsh(actual-(1-t)*ideal).min() > 0
        assert np.linalg.eigvalsh((1+t)*ideal-actual).min() > 0


def test_executed_cap_only_Ramsey_readout_uses_existing_source_services():
    row = c.witness(F(16))
    assert row["source_pulses"] == 2056
    assert .95 < row["complete_executed_bit_probabilities"][0] < .96
    assert .04 < row["complete_executed_bit_probabilities"][1] < .05
    assert abs(sum(row["complete_executed_bit_probabilities"])-1) < 2e-11
    assert row["same_outside_configuration_and_cap_charge_one"]
    assert row["outside_native_probability_diagnostic"] >= 0
    assert not row["padding_projected_or_outcomes_renormalized"]


def test_labelled_native_cut_instrument_commutes_with_clock_on_entangled_probe():
    audit = c.branch_covariance_audit()
    assert audit["ideal_completed_modular_action_commutes_with_each_Kraus_operator"]
    states = r.source()[0].states
    index = {x: i for i, x in enumerate(states)}
    seams = [e for e in r.source()[0].seams if not set(e).intersection(c.CAP)]
    rng = np.random.default_rng(89)
    probe = rng.normal(size=(66, 3))+1j*rng.normal(size=(66, 3))
    probe /= np.linalg.norm(probe)
    phases = np.array([float(w)**(-1j*16) for w in c.reference()])[:, None]
    for law in ("C", "E"):
        checked = 0
        for seam in seams:
            for coin in (0, 1):
                groups = {}
                for i, state in enumerate(states):
                    output, (_, current) = c.cut.native.native(state, seam, coin)
                    label = current if law == "C" else state[seam[0]]-state[seam[1]]
                    groups.setdefault(label, []).append((i, index[output]))
                for pairs in groups.values():
                    before, after = np.zeros_like(probe), np.zeros_like(probe)
                    for i, j in pairs:
                        before[j] += phases[i]*probe[i]
                        after[j] += probe[i]
                    after *= phases
                    np.testing.assert_allclose(before, after, atol=3e-16, rtol=0)
                    checked += 1
        assert checked == audit["laws"][law]["labelled_Kraus_branches"]


def test_clock_rejects_wrong_reference_and_parameter_contracts():
    with pytest.raises(ValueError, match="clock or cocycle"):
        c.word("selected", F(1))
    for bad in (F(0), F(3), 1.):
        with pytest.raises(ValueError, match="ratio"):
            c.logarithm_interval(bad)
    for bad in (17, 1., True):
        with pytest.raises(ValueError, match="parameter"):
            c.word("clock", bad)
