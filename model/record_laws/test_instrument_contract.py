"""Record-basis changes must not be mistaken for current readouts."""
from fractions import Fraction as F
import json
from pathlib import Path

import pytest

from . import family as f
from . import instrument_contract as ic


def test_exact_retained_instrument_counterexample():
    result = ic.run()
    assert result == json.loads(Path(__file__).with_name("instrument_results.json").read_text())
    assert result["matrix_units_checked"] == 6084
    assert all(result["checks"].values())
    witness = result["witness"]
    assert witness["basis_input_original_moving_record_probability"] == "1/2"
    assert witness["basis_input_remixed_plus_probability"] == "1/4"
    assert witness["coherent_input_remixed_plus_probability"] == "1/2"
    assert witness["coherent_input_remixed_minus_probability"] == "0"
    assert witness["incoherent_input_remixed_plus_probability"] == "1/4"
    assert witness["incoherent_input_remixed_minus_probability"] == "1/4"


def test_mixed_record_cannot_be_decoded_as_one_definite_current():
    mixed = ic.remixed_native_c()
    branch = mixed[("record_superposition", 1)]
    x, y = f.representative((1, 0)), f.representative((0, 1))
    # The same outcome can denote a moving or a waiting input.
    assert branch[x][0] == branch[y][0] == y
    assert branch[x][1] == branch[y][1] == F(1, 2)
    with pytest.raises(AssertionError, match="mixes currents"):
        f.extract_table(mixed)


def test_selective_map_fails_covariance_that_unread_channel_retains():
    branch = ic.remixed_native_c()[("record_superposition", 1)]
    x, y = f.representative((1, 0)), f.representative((0, 1))
    output = f.apply_instrument({"plus": branch}, {(x, y): f.S(1)})
    assert output == {(y, y): f.S(F(1, 4))}
    assert f.STATES[x] != f.STATES[y]
