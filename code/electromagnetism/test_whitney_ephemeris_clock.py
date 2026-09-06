"""Independent geometry, invariance and false-clock controls."""
from pathlib import Path, PurePosixPath, PureWindowsPath
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_ephemeris_clock as producer
import verify_whitney_ephemeris_clock as audit


def test_actual_clock_receipt():
    result = audit.verify(audit.load())
    assert result["accepted"]
    assert .001 < result["comparison_error"] < .002
    assert not result["calibrated_physical_clock"]


def test_action_from_independent_full_tetrahedron():
    rng = np.random.default_rng(85102)
    for _ in range(12):
        q, v = rng.normal(size=(2, 5))
        np.testing.assert_allclose(producer.action_at(q, v), audit.action_at(q, v), rtol=2e-14, atol=2e-13)


def test_regrading_and_subdivision_on_same_curve():
    q = [[.1, .2, .3, -.1, .1], [.2, .3, .1, -.1, .2]]
    expected = sum(producer.durations(q, order=12))
    assert abs(sum(producer.durations(q, order=12, warped=True))-expected) < 1e-13
    middle = ((np.array(q[0])+q[1])/2).tolist()
    assert abs(sum(producer.durations([q[0], middle, q[1]], order=12))-expected) < 1e-13
    assert abs(sum(audit.increments(q))-expected) < 1e-13


def test_configuration_consumer_cannot_read_time_or_velocity():
    class ConfigurationOnly(dict):
        def __getitem__(self, key):
            assert key == "q_exact", "clock accessed an unauthorized field"
            return super().__getitem__(key)
    frames = [ConfigurationOnly(q_exact=[str(i+j) for j in range(5)]) for i in range(2)]
    assert producer.read_configurations({"frames": frames}) == [[float(i+j) for j in range(5)] for i in range(2)]


@pytest.mark.parametrize("coordinates,energy", [
    ([[0]*5, [0]*5], 20), ([[0]*5, [1]*5], -1),
    ([[0]*5, [float('nan')]*5], 20), ([[0]*4, [1]*4], 20),
    ([[0]*5, [1]*5], float('inf'))])
def test_inadmissible_paths_rejected(coordinates, energy):
    with pytest.raises(ValueError):
        producer.durations(coordinates, energy=energy)


@pytest.mark.parametrize("key", ["calibrated_physical_clock", "quantum_clock_operator",
                                "rigorous_numerical_enclosure", "input_clock"])
def test_scope_promotions_rejected(key):
    packet = audit.load()
    packet["contract"][key] = True
    with pytest.raises(ValueError, match="contract"):
        audit.verify(packet)


@pytest.mark.parametrize("mutation", ["wrong_formula", "wrong_energy", "extra_field", "pin", "source"])
def test_contract_and_custody_mutations(mutation):
    packet = audit.load()
    if mutation == "wrong_formula":
        packet["contract"]["formula"] = "d_tau=one_repair_cycle"
    elif mutation == "wrong_energy":
        packet["contract"]["energy_exact"] = "1"
    elif mutation == "extra_field":
        packet["empirical_prediction"] = True
    elif mutation == "pin":
        packet["source_pins"][next(iter(packet["source_pins"]))] = "0"*64
    else:
        packet["source"]["consumed_fields"].append("frames[].action_time_exact")
    with pytest.raises(ValueError):
        audit.verify(packet)


@pytest.mark.parametrize("text", ['{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":1e9999}'])
def test_invalid_json(tmp_path, text):
    path = tmp_path/"bad.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        audit.load(path)


def test_fabricated_perfect_clock_rejected(monkeypatch):
    # Source provenance has its own exhaustive tests; here isolate the
    # independent metric integration against a coherently rescaled clock.
    import verify_whitney_charged_instrument as instrument
    monkeypatch.setattr(audit, "instrument_verifier", lambda: SimpleNamespace(load=instrument.load, verify=lambda _: None))
    packet = audit.load()
    for row in packet["refinements"]:
        factor = 2/row["duration"]
        row["increments"] = [x*factor for x in row["increments"]]
        row["duration"] = 2.0
    packet["cumulative_clock"] = np.linspace(0, 2, 81).tolist()
    with pytest.raises(ValueError, match="segment"):
        audit.verify(packet)


def test_verifier_loads_from_an_external_working_directory(tmp_path):
    # The ledger loads verifiers by file path, without sibling sys.path setup.
    import subprocess
    code = "import importlib.util; s=importlib.util.spec_from_file_location('clock',"+repr(str(Path(audit.__file__)))+"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); assert m.instrument_verifier().OUTPUT.is_file()"
    subprocess.run([sys.executable, "-c", code], cwd=tmp_path, check=True, capture_output=True)


@pytest.mark.parametrize("path_type", [PurePosixPath, PureWindowsPath])
def test_source_attachment_is_portable_across_path_flavors(monkeypatch, path_type):
    # Actual file reads remain local; relative path serialization uses the
    # selected platform flavor, so Windows behavior is exercised on Unix too.
    class SourcePath:
        def __init__(self, actual):
            self.actual = actual

        def relative_to(self, root):
            return path_type(*self.actual.relative_to(root).parts)

        def __fspath__(self):
            return str(self.actual)

        def __getattr__(self, name):
            return getattr(self.actual, name)

    source = SourcePath(producer.SOURCE)
    monkeypatch.setattr(producer, "SOURCE", source)
    monkeypatch.setattr(audit, "SOURCE", source)
    packet = producer.build()
    assert packet["source"]["path"] == "code/electromagnetism/runtime/whitney_charged_instrument_receipt.json"
    assert audit.verify(packet)["accepted"]


def test_noncanonical_source_path_rejected():
    packet = producer.build()
    packet["source"]["path"] = packet["source"]["path"].replace("/", "\\")
    with pytest.raises(ValueError, match="source attachment"):
        audit.verify(packet)
