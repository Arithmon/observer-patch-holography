"""Check the complete-table event bound and the frozen experiment binding."""
from fractions import Fraction as F

import pytest

from . import occupation


def table(probabilities, budget, width=F(0)):
    return {'branches': [
        {'bits': list(bits), 'polynomial_effect_probability_interval': [value-width, value+width]}
        for bits, value in zip(((0, 0), (0, 1), (1, 0), (1, 1)), probabilities)],
        'raw_polynomial_mass': sum(probabilities, F(0)),
        'raw_probability_half_l1_error_upper': budget}


def test_event_needs_actual_mass_correction_without_renormalizing():
    # p=(.25,.25,.25,.25), q=(.25,.35,.25,.25). Half-L1=.05,
    # but the fixed event error is .1: using half-L1 alone would fail.
    proof = occupation.terminal_event(table([F(1, 4), F(7, 20), F(1, 4), F(1, 4)], F(1, 20)))
    assert proof['raw_mass_correction'] == F(1, 20)
    assert proof['raw_event_probability_interval'] == [F(3, 5), F(3, 5)]
    assert proof['unrestricted_event_probability_interval'] == [F(1, 2), F(7, 10)]
    assert not proof['branches_renormalized']


def test_event_keeps_directed_terminal_effect_enclosure_width():
    proof = occupation.terminal_event(table([F(1, 4)]*4, F(1, 100), F(1, 1000)))
    assert proof['unrestricted_event_probability_interval'] == [F(61, 125), F(64, 125)]


def test_incomplete_record_table_cannot_supply_event_bound():
    proof = table([F(1, 4)]*4, F(0))
    proof['branches'].pop()
    with pytest.raises(ValueError, match='complete four-record'):
        occupation.terminal_event(proof)


def test_native_event_weight_and_discriminator_do_not_refit_sign():
    high = {'unrestricted_event_probability_interval': [F(2, 5), F(1, 2)]}
    low = {'unrestricted_event_probability_interval': [F(1, 5), F(3, 10)]}
    positive = occupation.signed_event_separation(high, low)
    assert positive['native_positive_current_weight_exact'] == F(5, 66)
    assert positive['complete_native_history_tv_lower'] == F(1, 132)
    reversed_actions = occupation.signed_event_separation(low, high)
    assert reversed_actions['signed_conditional_probability_separation_lower'] == -F(3, 10)
    assert reversed_actions['conditional_history_tv_lower'] == 0
    assert not reversed_actions['strict_positive_discriminator_established']


def test_original_contract_is_hash_bound(tmp_path, monkeypatch):
    assert occupation.retained_source()['all_pass']
    for name in ('protocol.py', 'results.json'):
        (tmp_path/name).write_bytes((occupation.FOLDER/name).read_bytes())
    monkeypatch.setattr(occupation, 'FOLDER', tmp_path)
    with (tmp_path/'results.json').open('ab') as stream:
        stream.write(b'\n')
    with pytest.raises(ValueError, match='binding changed: results.json'):
        occupation.retained_source()
