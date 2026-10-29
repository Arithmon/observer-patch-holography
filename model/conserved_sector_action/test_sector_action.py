from fractions import Fraction as F

import pytest

from .sector_action import certify_kernel, native_kernel
from . import sector_action as action


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


def test_frozen_predecessor_receipt_has_an_explicit_exact_capture():
    bindings = action.predecessor_bindings()
    assert set(bindings) == set(action.BASE_SOURCES)
    assert bindings[action.PREDECESSOR_ENERGY_SOURCE] == action.PREDECESSOR_ENERGY_SHA256
    current = action.source_bindings()
    assert action.PREDECESSOR_ENERGY_SOURCE not in current
    assert current[action.PREDECESSOR_ENERGY_CAPTURE] == action.PREDECESSOR_ENERGY_SHA256


def test_changed_predecessor_capture_is_rejected(monkeypatch, tmp_path):
    original = (action.ROOT/action.PREDECESSOR_ENERGY_CAPTURE).read_bytes()
    changed = tmp_path/action.PREDECESSOR_ENERGY_CAPTURE
    changed.parent.mkdir(parents=True)
    changed.write_bytes(original+b'\n')
    monkeypatch.setattr(action, 'ROOT', tmp_path)
    with pytest.raises(ValueError, match="predecessor capture changed"):
        action.report()


def test_live_energy_receipt_does_not_replace_the_frozen_input(monkeypatch):
    original = action.digest
    monkeypatch.setattr(action, 'digest', lambda path: '0'*64
                        if path == action.PREDECESSOR_ENERGY_SOURCE else original(path))
    assert action.predecessor_bindings()[action.PREDECESSOR_ENERGY_SOURCE] == action.PREDECESSOR_ENERGY_SHA256
    assert action.PREDECESSOR_ENERGY_SOURCE not in action.source_bindings()
