"""Hostile controls for the scope and arithmetic of the gauge scale contract."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


def module(name):
    spec = importlib.util.spec_from_file_location(name, HERE/(name+".py"))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


verifier = module("verify_k7_gauge_scale_contract_independent")


@pytest.fixture(scope="module")
def receipt():
    return verifier.load(verifier.OUT)


@pytest.fixture(scope="module")
def reconstructions():
    return verifier.reconstruct_native(), verifier.reconstruct_plane()


def test_live_independent_reconstruction(receipt, reconstructions, monkeypatch):
    native, plane = reconstructions
    monkeypatch.setattr(verifier,"reconstruct_native",lambda: native)
    monkeypatch.setattr(verifier,"reconstruct_plane",lambda: plane)
    assert verifier.verify(receipt)


@pytest.mark.parametrize("section,key,value", [
    ("native_signature","rg_subtraction_scale_identified",True),
    ("native_signature","physical_continuum_length_identified",True),
    ("native_signature","laboratory_energy_identified",True),
    ("native_signature","countermodel_proves_universal_physical_nonselection",True),
    ("conditional_plane_transport","requires_common_scheme_branch_and_matching",False),
    ("conditional_plane_transport","matter_kinetic_branch_physically_selected",True),
    ("conditional_plane_transport","arbitrary_threshold_or_scheme_shift_preserves_plane",True),
    ("k7_source_contract","three_couplings_common_scale_claimed",True),
    ("k7_source_contract","three_thirteenths_gut_boundary_and_dependent_rge_chain","valid"),
    ("trust_boundary","public_data_used_as_selection_input",True),
    ("trust_boundary","target_driven_bridge_choices",True),
    ("trust_boundary","k7_coupling_vector_constructed",True),
])
def test_semantic_promotions_rejected(receipt, section, key, value):
    changed = deepcopy(receipt); changed[section][key] = value
    with pytest.raises(ValueError, match="semantics"):
        verifier.verify(changed)


@pytest.mark.parametrize("section,key,value", [
    ("native_signature","counter_diffusion","4/99"),
    ("native_signature","one_attempt_law",["1/2","1/2"]),
    ("native_signature","branches",False),
    ("conditional_plane_transport","b",["11","1","-3"]),
    ("conditional_plane_transport","k",["2","2","2"]),
    ("conditional_plane_transport","k_annihilator","1"),
    ("conditional_plane_transport","mixed_component_transport_delta_D","0"),
    ("conditional_plane_transport","wrong_branch_transport_delta_D","0"),
    ("conditional_plane_transport","threshold_basis_delta_D",["0","0","0"]),
])
def test_arithmetic_and_branch_mutations_rejected(receipt, reconstructions, monkeypatch,
                                                 section, key, value):
    native, plane = reconstructions
    monkeypatch.setattr(verifier,"reconstruct_native",lambda: native)
    monkeypatch.setattr(verifier,"reconstruct_plane",lambda: plane)
    changed = deepcopy(receipt); changed[section][key] = value
    with pytest.raises(ValueError, match="reconstruction"):
        verifier.verify(changed)


def test_source_mutation_rejected(receipt, monkeypatch):
    original = Path.read_bytes
    def read(path):
        raw = original(path)
        return raw+b"\n" if str(path).endswith("model/observer/native.py") else raw
    monkeypatch.setattr(Path,"read_bytes",read)
    with pytest.raises(ValueError, match="source pin"):
        verifier.verify(receipt)


def test_duplicate_keys_rejected():
    with pytest.raises(ValueError, match="duplicate"):
        verifier.pairs([("x",1),("x",2)])
