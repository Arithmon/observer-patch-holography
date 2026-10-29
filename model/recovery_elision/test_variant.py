"""Independent native-column checks of the admitted image and service variant."""
from fractions import Fraction as F

import numpy as np

from . import variant as v
from ..moving_defect import transport as t


def test_all_eligible_flags_are_exact_identity_recoveries():
    census = v.certificate(resources=False)['census']
    assert census['all_raw_slots'] == 43200
    assert census['positive_flags'] == 11580
    assert census['identity_recovery_flags'] == 3120
    assert census['identity_flags_with_core_motion'] == 0
    # Rebuild from original recorded transitions, independently of branch().
    checked = 0
    for center in t.LIVE:
        for law in ('C', 'E'):
            for seam in t.SEAMS:
                for coin in (0, 1):
                    for raw in range(1 << t.circuit.RECORD_BITS[law]):
                        label = center, law, seam, coin, raw
                        if not v.identity_eligible(*label):
                            continue
                        source_columns = t.supports(center)
                        images = []
                        for column in source_columns:
                            output = []
                            for pair in column:
                                state = tuple(int(i in pair) for i in range(12))
                                target, record = t.records.recorded_step(law, state, seam, coin)
                                if record[-1] == t.circuit.record_values(law, coin).get(raw):
                                    output.append(tuple(i for i, n in enumerate(target) if n))
                            images.append(output)
                        assert tuple(map(frozenset, images)) == tuple(map(frozenset, source_columns))
                        assert t.branch_probability(*label) == t.local.coin_probability(t.parameters(), seam, coin)
                        checked += 1
    assert checked == 3120


def test_reference_entangled_image_uses_no_recovery_word():
    center = t.INITIAL
    label = next((center, law, seam, coin, raw)
                 for law in ('C', 'E') for seam in t.SEAMS for coin in (0, 1)
                 for raw in range(1 << t.circuit.RECORD_BITS[law])
                 if v.identity_eligible(center, law, seam, coin, raw))
    logical = np.array([[1, 1j], [1j, 1]], complex)/2
    image = np.sqrt(float(t.branch_probability(*label)))*(t.ideal_code(center)@logical)
    actual = v.recover(*label, image)
    assert np.array_equal(actual, image)
    assert all(value == 0 for value in v.recovery_resources(*label).values())
    assert np.linalg.norm(t.recover(*label, image)-actual) < 1e-8


def test_moving_flags_keep_existing_recovery_and_costs():
    center = t.INITIAL
    label = next((center, law, seam, coin, raw)
                 for law in ('C', 'E') for seam in t.SEAMS for coin in (0, 1)
                 for raw in range(1 << t.circuit.RECORD_BITS[law])
                 if (target := t.branch(center, law, seam, coin, raw)[1]) in t.LIVE
                 and target != center)
    assert not v.identity_eligible(*label)
    assert v.recovery_resources(*label) == t.recovery_resources(*label)


def test_stationary_pulse_savings_are_exact_positive_rationals():
    result = v.certificate()
    assert F(result['stationary_saved_source_pulses_per_native_attempt']) == F(
        4737999067212019707204144087023030301441659146,
        97909756773858989776401485931436511086245)
    assert result['largest_omitted_recovery_pulses'] == 196460
