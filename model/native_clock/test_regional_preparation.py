"""Local stopping order, quantum contraction, complete herald and clock gaps."""
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import regional_preparation as r


@pytest.fixture(scope="module")
def source_channels():
    channels = {}
    for law in ("C", "E"):
        incident = remote = None
        for edge in r.mixing.seams():
            a = r.seam_counts(edge, law)
            if edge in r.incident_edges():
                incident = a if incident is None else incident+a
            else:
                remote = a if remote is None else remote+a
        channels[law] = incident/60, remote/60  # Per global attempt weights.
    return channels


def test_exact_positive_channels_and_complete_local_herald():
    certificate = r.exact_audit()
    assert certificate["both_herald_outputs_retained"]
    assert all(x["exact_unital_trace_preserving_positive_self_adjoint"]
               for x in certificate["laws"].values())


@pytest.mark.parametrize("law", ("C", "E"))
def test_stopped_channel_order_and_coherent_similarity_bound(source_channels, law):
    incident, remote = source_channels[law]
    basis = r.native.states(binary=True)
    groups = {}
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            delta = tuple(x-y for x, y in zip(a, b, strict=True))
            groups.setdefault(delta, []).append(66*i+j)
    # The population block and a genuinely coherent transported block.
    selected = [(d, ix) for d, ix in groups.items() if not any(d)]
    selected += [next((d, ix) for d, ix in groups.items() if any(d) and len(ix) > 1)]
    for delta, ix in selected:
        a = incident[ix][:, ix].toarray()
        b = remote[ix][:, ix].toarray()
        identity = np.eye(len(ix))
        m = identity-b
        s = np.linalg.solve(m.T, a.T).T
        # Remote records precede the stopping repair; reversal is detectable.
        partial, power = np.zeros_like(a), identity.copy()
        for _ in range(60):
            partial += a@power
            power = power@b
        np.testing.assert_allclose(s, partial, atol=2e-14)
        if not any(delta):
            assert np.linalg.norm(s-np.linalg.solve(m, a)) > 1e-3
            dep = np.ones_like(s)/66
            np.testing.assert_allclose(s.sum(axis=0), 1, atol=2e-14)
            np.testing.assert_allclose(s.sum(axis=1), 1, atol=2e-14)
        else:
            dep = np.zeros_like(s)
        vals, vecs = np.linalg.eigh(m)
        invroot = (vecs/np.sqrt(vals))@vecs.T
        bb = invroot@a@invroot
        np.testing.assert_allclose(bb, bb.T, atol=2e-14)
        assert np.linalg.eigvalsh(bb).min() > -2e-14
        for n in (1, 8, 32):
            bound = np.sqrt(2)*float(r.R_UPPER)**n
            assert np.linalg.norm(np.linalg.matrix_power(s, n)-dep, 2) <= bound+2e-14


def test_local_reference_is_distinct_faithful_only_after_completion():
    raw, completed, q = r.reference()
    assert q == F(5, 33)
    assert raw.count(0) == 1
    assert min(completed) == F(1, 120)
    assert sum(completed) == 1
    assert completed != r.cut.completed_palm()
    assert r.Counter(completed) == {
        F(1, 120): 28, F(11, 600): 16, F(23, 1200): 16,
        F(2, 75): 5, F(1, 30): 1,
    }
    gaps = []
    for sector in r.cut.cut_sectors():
        values = {completed[row[0]] for row in sector.indices}
        gaps.append({b/a for a in values for b in values if b > a})
        for row in sector.indices:
            assert len({completed[i] for i in row}) == 1
    assert gaps == [set(), {F(23, 22)}, {F(5, 4)}]


def test_bounds_keep_failure_and_control_conditioning_cost():
    p = r.preparation()
    assert p["local_warmup_ticks_including_zeros"] == 496
    assert not p["global_warmup_attempt_counter_required"]
    lo, hi = map(F, p["success_probability_interval"])
    flo, fhi = map(F, p["failure_probability_interval"])
    assert lo+fhi == hi+flo == 1
    assert lo < F(5, 33) < hi
    assert r.local_bound(496) <= F(1, 2**24)
    assert r.local_bound(0) == 1
    assert r.local_bound(496) >= 33*F(2)*r.R_UPPER**496
    assert F(p["conditional_output_reference_uniform_halftrace_upper"]) < F(1, 100000)
    assert F(p["relative_pre_normalization_radius"]) < F(1, 10000)
    assert F(p["flow_halfdiamond_upper_per_absolute_modular_parameter"]) < F(1, 10000)


def test_next_tick_success_is_not_an_unconditional_stationary_sample(source_channels):
    raw, _, _ = r.reference()
    # Even under a uniform input, postselecting local movement changes state.
    assert raw != (F(1, 66),)*66
    assert max(raw) == F(1, 30)
    # Remote evolution is not generally identity before the local tick.
    incident, remote = source_channels["C"]
    ix = np.arange(66)*67
    a, b = incident[ix][:, ix].toarray(), remote[ix][:, ix].toarray()
    stopped = np.linalg.solve((np.eye(66)-b).T, a.T).T
    assert np.linalg.norm(stopped-2*a) > 0.1


@pytest.mark.parametrize("value", (-1, True, 1.5, 1000001))
def test_invalid_local_counts_fail_closed(value):
    with pytest.raises(ValueError):
        r.local_bound(value)


@pytest.mark.parametrize("value", (True, 1, 11, 41, 24.0))
def test_invalid_precision_fails_closed(value):
    with pytest.raises(ValueError):
        r.preparation(value)


def test_invalid_domains_fail_closed():
    for cap in ((), (1, 0), (0, 0), (False,), list(r.CAP)):
        with pytest.raises(ValueError):
            r.reference(cap)
    for law in (None, True, "CE", "F"):
        with pytest.raises(ValueError):
            r.seam_counts(r.mixing.seams()[0], law)
    with pytest.raises(ValueError):
        r.preparation(cap=tuple(range(10)))


def test_retained_receipt_replays_exactly():
    assert json.loads(r.OUTPUT.read_text(encoding="utf-8")) == r.run()
