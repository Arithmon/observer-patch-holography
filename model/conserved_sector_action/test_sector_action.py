from fractions import Fraction as F

import pytest

from .sector_action import certify_kernel, native_kernel


def fixture():
    return [{0: F(1, 2), 1: F(1, 2)}, {0: F(1, 2), 1: F(1, 2)},
            {2: F(1)}], (0, 0, 1), (F(1, 4), F(1, 4), F(1, 2))


def test_disconnected_kernel_is_checked_sectorwise():
    result = certify_kernel(*fixture())
    assert result["connected_sector_sizes"] == {"0": 2, "1": 1}
    assert result["reference_sector_masses"] == {"0": "1/2", "1": "1/2"}


def test_a_charge_changing_transition_is_rejected():
    rows, charges, reference = fixture()
    rows[0] = {0: F(1, 2), 2: F(1, 2)}
    with pytest.raises(ValueError, match="changes the declared sector"):
        certify_kernel(rows, charges, reference)


def test_hidden_disconnection_inside_one_declared_sector_is_rejected():
    rows, _, reference = fixture()
    with pytest.raises(ValueError, match="not connected"):
        certify_kernel(rows, (0, 0, 0), reference)


def test_nonstationary_reference_is_rejected():
    rows, charges, _ = fixture()
    with pytest.raises(ValueError, match="reciprocity"):
        certify_kernel(rows, charges, (F(1, 8), F(3, 8), F(1, 2)))


def test_nonstochastic_row_is_rejected():
    rows, charges, reference = fixture()
    del rows[0][0]
    with pytest.raises(ValueError, match="stochastic"):
        certify_kernel(rows, charges, reference)


def test_physical_law_must_be_explicit():
    with pytest.raises(ValueError, match="expected C or E"):
        native_kernel("other")
