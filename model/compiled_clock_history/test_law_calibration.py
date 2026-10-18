"""Independent scheduler, carrier restriction and calibration-bound checks."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import subprocess
import sys

import pytest

from ..repair_identification import calibration, core
from ..record_laws import family
from . import history, law_calibration as bridge, receipt, source


def schedules(cap):
    # Enumerate all Bernoulli incidence words, including first-wait stops.
    one = [(F(2, 5)*F(3, 5)**i, i+1, True) for i in range(cap)]
    one.append((F(3, 5)**cap, cap, False))
    rows = []
    for weight, count, complete in one:
        if not complete:
            rows.append((weight, count, False))
        else:
            rows.extend((weight*w, count+1+n, success) for w, n, success in one)
    return rows


@pytest.mark.parametrize('cap', (1, 2, 3))
def test_exposure_is_complete_independent_path_sum(cap):
    rows, result = schedules(cap), bridge.exposure(cap)
    assert sum(w for w, _, _ in rows) == 1
    assert result['native_attempts_mean'] == sum(w*n for w, n, _ in rows)
    assert result['complete_probability'] == sum(w for w, _, done in rows if done)
    assert result['native_attempts_maximum'] == max(n for _, n, _ in rows)
    assert result['native_attempts_mean'] < result['native_attempts_maximum']
    costs = receipt.costs(cap=cap)
    assert result['native_attempts_mean'] == F(costs['scheduled_native_attempts_ideal_expectation']) + F(
        costs['feedback_native_attempts_ideal_expectation'])


@pytest.mark.parametrize('cap', (1, 2, 3))
def test_expected_exposure_perturbation_bound_has_correct_first_order_coefficient(cap):
    # An independent visible error flag per call has complete TV equal to
    # P(at least one error). This is a generic sharpness control for the
    # adaptive-channel lemma, not another member of the native law family.
    rows = schedules(cap)
    delta = F(1, 10**6)
    actual = sum(w*(1-(1-delta)**n) for w, n, _ in rows)
    upper = bridge.exposure(cap)['native_attempts_mean']*delta
    second_order = sum(w*F(n*(n-1), 2)*delta**2 for w, n, _ in rows)
    assert upper-second_order <= actual <= upper
    assert actual > sum(w*(n-1)*delta for w, n, _ in rows)


def test_full_binary_carrier_is_the_calibrated_current_instrument_restriction():
    table = core.controls()['C']
    coefficients = family.multipliers(table)
    states, codes = source.STATES, source.codes()
    index = dict(zip(states, codes, strict=True))
    native_groups = dict(history._native_maps())[(0, 1)]
    actual = defaultdict(F)
    for current, pairs in native_groups:
        for (i, a), (j, b) in product(pairs, repeat=2):
            actual[current, i, j, a, b] += F(1, 2)
    expected = {}
    for x, y in product(states, repeat=2):
        for current in (-1, 0, 1):
            weight = coefficients.get((current, family.key(x), family.key(y)), family.S())
            if weight:
                assert not (weight.b or weight.c or weight.d)
                a, b = family.move(x, current), family.move(y, current)
                assert a in index and b in index
                expected[current, index[x], index[y], index[a], index[b]] = weight.a
    assert dict(actual) == expected
    assert len(states) == 66 and set(states) < set(family.STATES)


def centers():
    return tuple(p.a for p in core.predictions(core.controls()['C']))


def test_nominal_witness_uses_radius_not_diameter_and_improves_worst_branch(monkeypatch):
    monkeypatch.setattr(receipt, 'bound', lambda **kwargs: F(1, 10**12))
    result = bridge.from_tubes(centers(), (F(1, 10000),)*9)
    bounds = result['prediction_budget']
    assert result['calibration']['full_family_feasibility'] == 'witnessed_by_locked_candidate'
    assert result['nominal_probability_difference_upper'] == F(1, 10000)
    assert bounds['per_attempt_half_diamond_upper'] == F(7, 10000)
    assert bounds['native_attempts_mean'] == F(408, 125)
    assert bounds['ideal_complete_half_trace_upper'] == F(357, 156250)
    assert bounds['worst_branch_ideal_upper'] == F(7, 2000)
    assert result['physical_law_selected'] is False
    with pytest.raises(ValueError, match='outside'):
        bridge.from_tubes((F(0),)*9, (F(0),)*9)


def test_raw_non_tp_event_bound_cannot_be_halved_or_capped_at_one(monkeypatch):
    eps = receipt.ERROR
    monkeypatch.setattr(receipt, 'bound', lambda **kwargs: eps/10)
    bounds = bridge.complete_bound(F(1))
    assert bounds['raw_to_calibrated_ideal_half_trace_upper'] == 1+eps
    assert bounds['raw_event_weight_error_upper'] == 1+2*eps
    # The scalar CP map rho -> (1+2eps)rho saturates this factor two.
    bounds = bridge.complete_bound(F(0))
    assert bounds['raw_event_weight_error_upper'] == 2*eps
    assert bounds['raw_compilation_half_trace_upper'] == eps


def test_portable_envelope_still_recertifies_actual_emitted_arrays(monkeypatch):
    calls = []
    def certify(**kwargs):
        calls.append(kwargs)
        return receipt.ERROR/3
    monkeypatch.setattr(receipt, 'bound', certify)
    assert bridge.complete_bound(F(0))['raw_compilation_half_trace_upper'] == receipt.ERROR
    assert calls == [{'cap': 2, 'gamma': 1, 'reference': False}]
    monkeypatch.setattr(receipt, 'bound', lambda **kwargs: 2*receipt.ERROR)
    with pytest.raises(ArithmeticError, match='envelope'):
        bridge.complete_bound(F(0))


def test_conditional_forecasts_charge_rare_selection():
    assert bridge.event_interval(F(1, 4), F(1, 10)) == (F(3, 20), F(7, 20))
    assert bridge.conditional_interval(F(1, 4), F(1, 2), F(1, 10)) == (F(1, 4), F(7, 8))
    with pytest.raises(ValueError):
        bridge.conditional_interval(F(1, 100), F(1, 100), F(1, 100))
    with pytest.raises(ValueError):
        bridge.conditional_interval(F(3, 4), F(1, 2), F(0))


def test_complete_count_provenance_and_abort_mass_are_preserved(monkeypatch):
    monkeypatch.setattr(receipt, 'bound', lambda **kwargs: F(0))
    attempted = 10000
    rows = []
    for probe, probability in zip(core.PROBES, centers(), strict=True):
        bins = {key: 0 for key in calibration.OUTCOMES}
        target = probability*attempted
        assert target.denominator == 1
        bins[f'current_{probe.current}_plus'] = int(target)
        bins['abort'] = attempted-int(target)
        rows.append({'probe': probe.name, 'attempted': attempted, 'counts': bins})
    evidence = {'schema': calibration.SCHEMA, 'settings': rows}
    design = calibration.sampling_design((attempted,)*9, (F(1, 10),)*9)
    result = bridge.from_counts(evidence, design, (F(0),)*9)
    count = result['count_evidence']
    assert count['abort_counts'] == [row['counts']['abort'] for row in rows]
    assert count['sampling_design_prior_declaration_verified'] is False
    assert count['premises_established_by_counts'] is False
    assert count['simultaneous_sampling_failure_bound'] == '9/400'
    evidence['settings'][0]['attempted'] -= 1
    with pytest.raises(ValueError, match='sum'):
        bridge.from_counts(evidence, design, (F(0),)*9)


@pytest.mark.parametrize('bad', (True, 0.1, '0', F(-1), F(2)))
def test_invalid_bounds_rejected(bad):
    with pytest.raises(ValueError):
        bridge.complete_bound(bad)


@pytest.mark.parametrize('cap', (True, 0, 4, 2.0))
def test_invalid_caps_rejected(cap):
    with pytest.raises(ValueError):
        bridge.exposure(cap)


def test_retained_forecast_replays_without_measured_or_physical_claim():
    root = Path(__file__).resolve().parents[2]
    result = bridge.run()
    retained = json.loads((root/'model/compiled_clock_history/law_calibration_results.json').read_text())
    assert result == retained
    assert result['physical_dictionary_validated'] is False
    assert result['intervention_contrast']['nominal_probability_difference'] == '1056/15625'
    assert F(result['intervention_contrast']['calibrated_interval'][0]) > 0
    assert [r['nominal_complete_probability'] for r in result['executions']] == ['256/625']*2


def test_complete_producer_imports_are_bound():
    root = Path(__file__).resolve().parents[2]
    script = '''
import json,sys
from pathlib import Path
from model.compiled_clock_history import law_calibration
root=Path.cwd().resolve()
paths=[]
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if not name: continue
    try: path=Path(name).resolve().relative_to(root).as_posix()
    except ValueError: continue
    if path.endswith('.py') and path.startswith(('model/','code/','tools/')):
        paths.append(path)
print(json.dumps(sorted(set(paths))))
'''
    observed = json.loads(subprocess.check_output([sys.executable, '-c', script], cwd=root, text=True))
    retained = json.loads((root/'model/compiled_clock_history/law_calibration_results.json').read_text())
    assert set(observed) <= set(retained['source_sha256'])
