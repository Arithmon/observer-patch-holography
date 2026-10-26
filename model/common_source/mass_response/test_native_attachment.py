"""Complete native path checks with independent integer updates and references."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from model.observer import native
from . import experiment as e, native_attachment as a


def fake_pair(abort_mass=0.):
    """Tiny retained quantum columns exercise the classical join cheaply."""
    results = []
    for enabled in (True, False):
        branches = {record: np.zeros((2, 2), complex) for record in product((0, 1), repeat=5)}
        for weak, mass, x, y in product((0, 1), repeat=4):
            record = (weak, mass, 0 if enabled else 1, x, y)
            # Same old prefix/clock probabilities, distinct later receiver.
            branches[record] = np.eye(2, dtype=complex) * math.sqrt((1 - abort_mass) / 32)
        aborts = {("scratch_abort", 0, 1, 0): np.eye(2, dtype=complex) * math.sqrt(abort_mass / 2)}
        prefix = {(weak, mass): (1 - abort_mass) / 4 for weak, mass in product((0, 1), repeat=2)}
        results.append(e.Result(branches, aborts, prefix, {},
                                {"kick_enabled": enabled, "record_order": a.RECORD_ORDER,
                                 "preparation": "one shared known preparation"}))
    return tuple(results)


def independent_update(loads, seam, coin):
    out = list(loads)
    low, remainder = divmod(loads[seam[0]] + loads[seam[1]], 2)
    out[seam[0]], out[seam[1]] = (low, low + remainder) if coin == 0 else (low + remainder, low)
    return tuple(out)


def test_every_native_event_continues_actual_loads_and_retains_zero_outcomes():
    kick, no_kick = fake_pair()
    paths = list(a.iter_paths(kick, no_kick))
    normal = [path for path in paths if not path.aborted]
    aborts = [path for path in paths if path.aborted]
    assert len(normal) == 66 * 2 * 32 * 2
    assert len(aborts) == 66 * 2  # The explicit zero abort remains represented.
    assert sum(path.probability == 0 for path in normal) == len(normal) / 2
    for path in paths:
        assert path.after_first_loads == independent_update(path.initial_loads, native.E, path.first_record[0])
        assert path.first_record[1] == path.initial_loads[0] - path.after_first_loads[0]
        assert path.quantum_case == ("kick" if path.first_record[1] else "no_kick")
        if path.aborted:
            assert path.quantum_scale_squared == F(1, 132)
            assert path.second_seam is path.second_record is None
            assert path.final_loads == path.after_first_loads
        else:
            assert path.quantum_scale_squared == F(1, 264)
            assert path.second_seam == (native.F if path.quantum_record[2] == 0 else native.E)
            assert path.final_loads == independent_update(path.after_first_loads, path.second_seam,
                                                         path.second_record[0])
        assert sum(path.final_loads) == 2


def test_exact_movement_weights_and_feedback_changes_the_final_native_law():
    summary = a.attach(*fake_pair())
    assert summary["kick_probability_exact"] == "5/33"
    assert summary["no_kick_probability_exact"] == "28/33"
    assert summary["actual"]["normal_paths"] == 8448
    assert summary["actual"]["total_probability"] == pytest.approx(1, abs=2e-14)
    assert summary["actual"]["first_record_marginal_residual"] < 2e-14
    assert summary["quantum_mixture_residual"] < 2e-14
    assert summary["final_native_classical_total_variation"] > .01
    assert summary["public_native_receiver_transcript_total_variation"] > .01
    assert all(row["between_settings_residual"] < 2e-14
               for row in summary["pre_kick_mass_and_clock_invariance"].values())
    fixed = a.attach(*fake_pair(), feedback=False)
    assert fixed["final_native_classical_total_variation"] == 0
    assert fixed["public_native_receiver_transcript_total_variation"] == 0


def test_abort_mass_has_only_the_first_coin_and_full_quantum_reference():
    kick, no_kick = fake_pair(abort_mass=.2)
    paths = list(a.iter_paths(kick, no_kick))
    for path in (paths[0], next(path for path in paths if path.aborted)):
        output = a.quantum_output(path, kick, no_kick)
        source = {"kick": kick, "no_kick": no_kick}[path.quantum_case]
        expected = (source.aborts if path.aborted else source.branches)[path.quantum_record]
        np.testing.assert_allclose(output, math.sqrt(float(path.quantum_scale_squared)) * expected,
                                   rtol=0, atol=0)
        assert np.vdot(output, output).real == pytest.approx(path.probability, abs=2e-16)
        assert output.shape == (2, 2)
    summary = a.attach(kick, no_kick)
    assert summary["actual"]["abort_probability"] == pytest.approx(.2, abs=2e-14)
    assert summary["actual"]["normal_probability"] == pytest.approx(.8, abs=2e-14)
    assert summary["actual"]["total_probability"] == pytest.approx(1, abs=2e-14)


def test_actual_continuous_history_keeps_clock_mass_and_native_feedback():
    kick = e.execute(backend="continuous")
    no_kick = e.execute(backend="continuous", kick=False)
    summary = a.attach(kick, no_kick)
    assert summary["actual"]["normal_paths"] == 8448
    assert summary["actual"]["abort_paths"] == 0
    assert summary["final_native_classical_total_variation"] > 1e-4
    assert summary["public_native_receiver_transcript_total_variation"] > 1e-4
    assert summary["actual"]["total_probability"] == pytest.approx(1, abs=2e-13)
    assert all(row["between_settings_residual"] < 2e-13
               for row in summary["pre_kick_mass_and_clock_invariance"].values())
    assert not summary["native_initial_coherence_claimed"]
    assert not summary["physical_particle_or_clock_identification"]


def test_mismatched_incomplete_or_hidden_state_inputs_are_rejected():
    kick, no_kick = fake_pair()
    missing = dict(kick.branches)
    missing.pop(next(iter(missing)))
    with pytest.raises(ValueError):
        a.attach(replace(kick, branches=missing), no_kick)
    with pytest.raises(ValueError):
        a.attach(replace(kick, specification={**kick.specification, "preparation": "different"}), no_kick)
    with pytest.raises(ValueError):
        a.attach(replace(kick, branches={key: value * .9 for key, value in kick.branches.items()}), no_kick)
    with pytest.raises(ValueError):
        a.attach(kick, no_kick, feedback=1)
    with pytest.raises(ValueError):
        a.case_from_current(np.array([0] * 12))
    with pytest.raises(ValueError):
        a.seam_from_receiver({"hidden_loads": [0] * 12})
