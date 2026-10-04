#!/usr/bin/env python3
"""Smoke tests for the empirical e+e- hadron closure policy."""

from __future__ import annotations

import json
import pathlib

import yaml


ROOT = pathlib.Path(__file__).resolve().parents[2]
POLICY = ROOT.parents[0] / "docs" / "policies" / "HADRON.md"
REGISTRY = ROOT / "particles" / "hadron" / "empirical_ee_hadrons_sources.yaml"
SCHEMA = ROOT / "particles" / "hadron" / "empirical_ee_hadronic_spectral_measure.schema.json"


def test_hadron_policy_names_empirical_output_class() -> None:
    text = POLICY.read_text(encoding="utf-8")

    assert "oph_plus_empirical_hadron_closure" in text
    assert "e^+e^-\\to\\mathrm{hadrons}" in text
    assert "promotable_as_OPH_source_theorem = false" in text


def test_empirical_source_registry_keeps_source_only_guard_false() -> None:
    payload = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_empirical_ee_hadrons_source_registry"
    assert payload["fine_structure_empirical_payload"]["source_only"] is False
    assert payload["fine_structure_empirical_payload"]["empirical_hadron_closure"] is True
    assert payload["fine_structure_empirical_payload"]["promotable_as_oph_source_theorem"] is False
    source_ids = {source["id"] for source in payload["candidate_sources"]}
    assert "pdg_hadronic_xsections_r_ratio" in source_ids
    assert "hepdata_ee_hadrons" in source_ids
    assert "alphaqed_jegerlehner" in source_ids


def test_empirical_spectral_measure_schema_requires_guards_and_covariance() -> None:
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))

    assert schema["properties"]["artifact"]["const"] == "oph_empirical_ee_hadronic_spectral_measure"
    assert schema["properties"]["row_class"]["const"] == "oph_plus_empirical_hadron_closure"
    required = set(schema["required"])
    assert "guards" in required
    assert "covariance" in required
    assert "correction_policy" in required
    guard_props = schema["properties"]["guards"]["properties"]
    assert guard_props["source_only"]["const"] is False
    assert guard_props["empirical_hadron_closure"]["const"] is True
    assert guard_props["promotable_as_oph_source_theorem"]["const"] is False

