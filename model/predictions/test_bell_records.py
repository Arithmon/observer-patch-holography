"""Independent probability, record-selection and native-boundary controls."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import subprocess
import sys

import pytest

from . import bell_records as b


def valid_row():
    return (1, 2, 5426350, 0, 5675700, 1, 0, 0, 10375, 10370,
            11000, 0, 0, 0, 0, 0)


def change(row, index, value):
    values = list(row)
    values[index] = value
    return tuple(values)


def test_every_local_strategy_and_shared_randomness_bound():
    strategies = b.local_strategies()
    assert len(set(strategies)) == 16
    assert sorted(sum(w) for _, w in strategies) == [1]*8+[3]*8
    # An arbitrary nonuniform mixture cannot improve the extremal bound.
    assert sum(F(i+1, 136)*sum(w)/4 for i, (_, w) in enumerate(strategies)) <= F(3, 4)
    # Joint access to BOTH settings invalidates the local-response premise.
    assert sum(int((-1)**(a*c)*1*((-1)**(a*c)) == 1) for a, c in b.SETTINGS) == 4


def test_native_local_example_against_uniform_sector_marginal():
    result = b.native_local_example()
    # Two distinct ports of a uniform size-two subset of twelve.
    assert result['counts'] == [[4, 40, 40, 180]]*4
    assert result['CHSH'] == '26/33'
    assert result['winning_probability'] == '79/132'
    assert result['complete_source_histories'] == 66*4*4


def test_published_eligibility_edges_and_outcome_independence():
    row = valid_row()
    assert b.eligible(row)
    assert not b.eligible(change(row, 2, 5426349))
    assert b.eligible(change(row, 2, 5478799))
    assert not b.eligible(change(row, 2, 5478800))
    assert not b.eligible(change(row, 4, 5675699))
    assert not b.eligible(change(row, 4, 5728150))
    for i in (14, 15):
        assert not b.eligible(change(row, i, 250))
        assert b.eligible(change(row, i, 251))
    for i in (12, 13):
        assert not b.eligible(change(row, i, 1))
    for a, c, t, u in product((0, 1), (0, 1), (0, 11000), (0, 11000)):
        other = list(row)
        other[6:8] = a, c
        other[10:12] = t, u
        assert b.eligible(tuple(other))


def test_readout_endpoints_keep_every_eligible_trial():
    rows = []
    for a, c in b.SETTINGS:
        for t in (10620, 10621, 14320, 14321):
            row = list(valid_row())
            row[6:8] = a, c
            row[10:12] = t, t
            rows.append(tuple(row))
    result = b.score_data(rows)
    assert result['eligible_trials'] == 16
    assert result['counts_by_setting_and_outcome'] == [[2, 0, 0, 2]]*4


def test_parser_and_source_selection_counts():
    result = b.score_data(b.parse_data(b.DATA.read_bytes()))
    assert result['deposited_candidate_rows'] == 4746
    assert result['eligible_trials'] == 245
    assert result['wins'] == 196
    assert result['setting_totals'] == [53, 79, 62, 51]
    assert result['counts_by_setting_and_outcome'] == [
        [23, 3, 4, 23], [33, 11, 5, 30], [22, 10, 6, 24], [4, 20, 21, 6]]
    assert F(result['CHSH_estimate']) == F(16036094, 6619647)
    with pytest.raises(ValueError):
        b.parse_data(b'')
    first = b.DATA.read_bytes().splitlines()[0]
    with pytest.raises(ValueError):
        b.parse_data(first+b',0')
    parts = first.split(b',')
    parts[7] = b'2'
    with pytest.raises(ValueError):
        b.parse_data(b','.join(parts))


def test_exact_tail_matches_independent_dynamic_program():
    q = F(7, 9)
    law = [F(1)]
    for n in range(1, 15):
        law = [(law[k] if k < len(law) else 0)*(1-q)
               +(law[k-1] if k else 0)*q for k in range(len(law)+1)]
        for threshold in range(n+1):
            assert b.tail(n, threshold, q) == sum(law[threshold:])
    # A process whose conditional success chance varies with prior results.
    def adversarial(prefix, mass):
        if len(prefix) == 8:
            return mass if sum(prefix) >= 7 else F(0)
        p = q if sum(prefix) % 2 else q/2
        return adversarial(prefix+(1,), mass*p)+adversarial(prefix+(0,), mass*(1-p))
    assert adversarial((), F(1)) <= b.tail(8, 7, q)
    assert b.tail(3, 0, F(0)) == 1
    assert b.tail(3, 1, F(0)) == 0
    assert b.tail(3, 3, F(1)) == 1
    for args in ((True, 0, q), (3, True, q), (3, 4, q), (3, 2, 0.75)):
        with pytest.raises(ValueError):
            b.tail(*args)


def test_exact_pvalue_enclosure_and_scope():
    result = b.run()
    test = result['conditional_null_test']
    lo, hi = map(F, test['fixed_n_tail_interval'])
    q = F(3, 4)+3*(F(27, 2500000)+F(27, 2500000)**2)
    assert lo <= b.tail(245, 196, q) <= hi
    assert F('0.03917464355') < lo <= hi < F('0.03917464357')
    assert test['optional_stopping_certificate_supplied'] is False
    assert result['scope']['whole_OPH_excluded'] is False
    assert result['scope']['quantum_C_or_E_selected'] is False
    assert result['scope']['new_native_empirical_confirmation'] is False


def test_cli_replays_and_rejects_changed_receipt(tmp_path):
    receipt = Path(__file__).with_name('bell_results.json')
    command = [sys.executable, '-m', 'model.predictions.bell_records', '--check']
    assert subprocess.run(command+[str(receipt)], cwd=b.ROOT, capture_output=True).returncode == 0
    data = json.loads(receipt.read_text())
    data['published_selection_analysis']['wins'] += 1
    bad = tmp_path/'bad.json'
    bad.write_text(json.dumps(data))
    assert subprocess.run(command+[str(bad)], cwd=b.ROOT, capture_output=True).returncode != 0
