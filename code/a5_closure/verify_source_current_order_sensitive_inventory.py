#!/usr/bin/env python3
"""Independent verifier for the Stage-2 order-sensitive source inventory.

This module deliberately does not import the producer.  It checks custody,
classification semantics, the one-packet qualification rule, decisive source
facts, and the fixture-import firewall independently.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_INVENTORY = HERE / "manifests" / "source_current_order_sensitive_inventory.json"
PRODUCER = HERE / "source_current_order_sensitive_inventory.py"

SCHEMA = "oph.source_current_order_sensitive_inventory.v2"
VERDICT = "SOURCE_CURRENT_ORDER_SENSITIVE_OBJECT_NOT_PRESENT"
PREVIOUS_INVENTORY_BASE_SHA = '543298e06ce41d47f44c7707fdb2c32723cca10e'
UPSTREAM_MAIN_SHA = 'ba84976ad1195892c68fdd0d74a27ab0c899aeb8'

AUDITED_DIRECTORIES = (
    "code/source_feedback_transport",
    "code/source_scalar_instruments",
    "code/source_routing",
    "code/a5_closure",
    "code/consensus",
    "code/refinement",
    "code/angular_sprint",
    "Lean/Dynamics",
    "Lean/Time",
    "Lean/Screen",
    "Lean/Geometry",
    "Lean/InformationProjection",
    "Lean/Variational",
    "code/source_read_routing",
    "code/source_routing_refinement",
    "code/source_native_programs",
    "code/source_native_accumulator",
    "code/source_temporal_acceptance",
    "code/source_read_acceptance",
    "code/source_checkpoint_selection",
    "code/source_read_selection",
    "code/source_encoded_memory",
    "code/source_passive_memory",
    "code/source_reusable_bus",
    "code/source_native_updates",
    "code/source_scalar_finite_instrument",
    "code/source_operator_join",
    "code/sm_fermion_current",
    "code/maxwell_measurement",
    "evidence/source_net_causal_poset/routed_read_law",
    "code/source_publication_selection",
    "code/source_operation_reads",
    "code/source_selection_model",
    "code/source_record_gluing",
    "code/rg_principle",
)

EXPECTED_SCOPE_EXTENSION = (
    "code/source_read_routing",
    "code/source_routing_refinement",
    "code/source_native_programs",
    "code/source_native_accumulator",
    "code/source_temporal_acceptance",
    "code/source_read_acceptance",
    "code/source_checkpoint_selection",
    "code/source_read_selection",
    "code/source_encoded_memory",
    "code/source_passive_memory",
    "code/source_reusable_bus",
    "code/source_native_updates",
    "code/source_scalar_finite_instrument",
    "code/source_operator_join",
    "code/sm_fermion_current",
    "code/maxwell_measurement",
    "evidence/source_net_causal_poset/routed_read_law",
    "code/source_publication_selection",
    "code/source_operation_reads",
    "code/source_selection_model",
    "code/source_record_gluing",
    "code/rg_principle",
)

CONTENT_SNAPSHOT_EXACT_EXCLUSIONS = (
    "code/a5_closure/manifests/source_current_order_sensitive_inventory.json",
)
CONTENT_SNAPSHOT_EXCLUDED_DIRECTORY_NAMES = ("__pycache__",)
CONTENT_SNAPSHOT_EXCLUDED_SUFFIXES = (".pyc",)

QUALIFICATION_FIELDS = (
    "source_native",
    "port_indexed",
    "reversible",
    "both_composition_orders_recorded",
    "raw_histories_serialized",
    "target_free",
    "same_twelve_port_carrier",
    "refinement_provenance_present",
    "twelve_reversible_perturbation_families",
)

EXPECTED_CLASSIFICATIONS = {
    "registered_adjacency_response": "COMMUTATIVE_ONLY",
    "source_seam_path_tomography": "IRREVERSIBLE_ONLY",
    "source_feedback_transport": "IRREVERSIBLE_ONLY",
    "source_scalar_sequential_instrument": "IRREVERSIBLE_ONLY",
    "record_counting_repair_628": "IRREVERSIBLE_ONLY",
    "directed_seam_repair": "IRREVERSIBLE_ONLY",
    "consensus_settling_packets": "IRREVERSIBLE_ONLY",
    "carrier_and_discrete_refinement_scaffold": "STATIC_ONLY",
    "issue_566_bracket_search_space": "MISSING_ORDER_INFORMATION",
    "b14_oriented_face_bracket": "STATIC_ONLY",
    "static_response_grammar_and_poles": "STATIC_ONLY",
    "spin_deck_transport": "STATIC_ONLY",
    "declared_port_current_fixture": "DOWNSTREAM_CONTAMINATED",
    "downstream_kinetic_form_selector": "DOWNSTREAM_CONTAMINATED",
    "generic_transport_from_gluing": "UNRELATED",
    "lean_pair_mean_memory_and_reusable_bus": "IRREVERSIBLE_ONLY",
    "lean_selected_pair_mean_histories": "IRREVERSIBLE_ONLY",
    "lean_source_history_abstractions": "PARTIAL",
    "lean_conditional_a2_words": "COMMUTATIVE_ONLY",
    "lean_reversible_flow_models": "MISSING_ORDER_INFORMATION",
    "lean_declared_fin12_transport_refinement": "PARTIAL",
    "lean_irreversible_provenance_histories": "IRREVERSIBLE_ONLY",
    "lean_transport_and_current_boundaries": "STATIC_ONLY",
    "source_read_routing_full_family": "IRREVERSIBLE_ONLY",
    "native_stored_programs_and_accumulator": "IRREVERSIBLE_ONLY",
    "native_temporal_tomography_and_checkpoint_selection": "IRREVERSIBLE_ONLY",
    "finite_source_operator_join": "STATIC_ONLY",
    "fermionic_hypercharge_current_histories": "DOWNSTREAM_CONTAMINATED",
    "maxwell_measurement_adapter": "DOWNSTREAM_CONTAMINATED",
    "native_eventual_publication": "IRREVERSIBLE_ONLY",
    "finite_state_process_selection": "STATIC_ONLY",
    "registered_operation_read_histories": "IRREVERSIBLE_ONLY",
    "declared_flagged_recovery_extensions": "PARTIAL",
}


def expected_properties(*true_fields: str) -> dict[str, bool]:
    return {field: field in true_fields for field in QUALIFICATION_FIELDS}


EXPECTED_PROPERTIES = {
    "registered_adjacency_response": expected_properties(
        "source_native",
        "port_indexed",
        "reversible",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "source_seam_path_tomography": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "source_feedback_transport": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "source_scalar_sequential_instrument": expected_properties("target_free"),
    "record_counting_repair_628": expected_properties(
        "source_native",
        "port_indexed",
        "raw_histories_serialized",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "directed_seam_repair": expected_properties(
        "source_native", "target_free", "same_twelve_port_carrier"
    ),
    "consensus_settling_packets": expected_properties(
        "source_native", "raw_histories_serialized", "target_free"
    ),
    "carrier_and_discrete_refinement_scaffold": expected_properties(
        "source_native",
        "port_indexed",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "issue_566_bracket_search_space": expected_properties(
        "port_indexed", "target_free", "same_twelve_port_carrier"
    ),
    "b14_oriented_face_bracket": expected_properties(
        "source_native", "port_indexed", "target_free", "same_twelve_port_carrier"
    ),
    "static_response_grammar_and_poles": expected_properties(
        "reversible", "target_free", "same_twelve_port_carrier"
    ),
    "spin_deck_transport": expected_properties(
        "source_native",
        "reversible",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "declared_port_current_fixture": expected_properties(
        "port_indexed",
        "reversible",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "downstream_kinetic_form_selector": expected_properties(
        "same_twelve_port_carrier"
    ),
    "generic_transport_from_gluing": expected_properties("target_free"),
    "lean_pair_mean_memory_and_reusable_bus": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "lean_selected_pair_mean_histories": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "lean_source_history_abstractions": expected_properties("target_free"),
    "lean_conditional_a2_words": expected_properties(
        "port_indexed", "reversible", "target_free", "same_twelve_port_carrier"
    ),
    "lean_reversible_flow_models": expected_properties("reversible", "target_free"),
    "lean_declared_fin12_transport_refinement": expected_properties(
        "port_indexed", "reversible", "target_free", "same_twelve_port_carrier"
    ),
    "lean_irreversible_provenance_histories": expected_properties(
        "source_native",
        "raw_histories_serialized",
        "target_free",
        "refinement_provenance_present",
    ),
    "lean_transport_and_current_boundaries": expected_properties(
        "source_native", "target_free", "same_twelve_port_carrier"
    ),
    "source_read_routing_full_family": expected_properties(
        "source_native",
        "port_indexed",
        "raw_histories_serialized",
        "target_free",
        "refinement_provenance_present",
    ),
    "native_stored_programs_and_accumulator": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "native_temporal_tomography_and_checkpoint_selection": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "finite_source_operator_join": expected_properties("target_free"),
    "fermionic_hypercharge_current_histories": expected_properties(
        "reversible", "raw_histories_serialized"
    ),
    "maxwell_measurement_adapter": expected_properties(),
}
EXPECTED_PROPERTIES.update({
    'native_eventual_publication': expected_properties('source_native', 'port_indexed', 'raw_histories_serialized', 'target_free'),
    'finite_state_process_selection': expected_properties('target_free'),
    'registered_operation_read_histories': expected_properties('source_native', 'port_indexed', 'raw_histories_serialized', 'target_free'),
    'declared_flagged_recovery_extensions': expected_properties('port_indexed', 'reversible', 'target_free'),
})


EXPECTED_PRIOR_BASE_REVIEW = {
    "Lean/Geometry/SourceAccumulatorAxiomAudit.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceAccumulatorDecoder.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceAccumulatorProgram.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankBusWitness.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankCompiler.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankControl.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankExecution.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankInvariant.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankLowering.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankMachine.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankPrecision.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceBankTrace.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceCheckpointAxiomAudit.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceCheckpointEntropy.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceCheckpointMeaning.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceCheckpointNative.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceCheckpointPipeline.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceCheckpointPolicy.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceNativeAccumulator.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceNativeCore.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceNativeProgramBudget.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceNativeProgramError.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceNativeProgramsAxiomAudit.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceNativeShuttle.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceNativeStoredProgram.lean": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "Lean/Geometry/SourceReadAcceptance.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceReadAcceptanceAxiomAudit.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceReadAcceptanceSchedule.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceReadRouting.lean": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "Lean/Geometry/SourceRoutingBudget.lean": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "Lean/Geometry/SourceRoutingHierarchy.lean": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "Lean/Geometry/SourceRoutingRefinementAxiomAudit.lean": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "Lean/Geometry/SourceRoutingStorage.lean": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "Lean/Geometry/SourceTemporalAcceptanceAxiomAudit.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalAttempts.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalConnected.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalErasure.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalEssential.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalGuard.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalMenu.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalNative.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalObservation.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalSelection.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "Lean/Geometry/SourceTemporalTomography.lean": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/maxwell_measurement/README.md": ("NEW_CANDIDATE", 'maxwell_measurement_adapter'),
    "code/maxwell_measurement/contract.py": ("NEW_CANDIDATE", 'maxwell_measurement_adapter'),
    "code/sm_fermion_current/README.md": ("NEW_CANDIDATE", 'fermionic_hypercharge_current_histories'),
    "code/sm_fermion_current/coupled_README.md": ("NEW_CANDIDATE", 'fermionic_hypercharge_current_histories'),
    "code/sm_fermion_current/coupled_receipt.json": ("NEW_CANDIDATE", 'fermionic_hypercharge_current_histories'),
    "code/sm_fermion_current/current_receipt.json": ("NEW_CANDIDATE", 'fermionic_hypercharge_current_histories'),
    "code/sm_fermion_current/quantum_link_receipt.json": ("NEW_CANDIDATE", 'fermionic_hypercharge_current_histories'),
    "code/source_checkpoint_selection/CONTRACT.md": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_checkpoint_selection/README.md": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_checkpoint_selection/controls.json": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_checkpoint_selection/receipt.json": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_encoded_memory/CONTRACT.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_encoded_memory/README.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_encoded_memory/controls.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_encoded_memory/receipt.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_feedback_transport/README.md": ("EXTENDS_EXISTING_CANDIDATE", 'source_feedback_transport'),
    "code/source_native_accumulator/CONTRACT.md": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_accumulator/README.md": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_accumulator/controls.json": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_accumulator/receipt.json": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_programs/CONTRACT.md": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_programs/README.md": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_programs/controls.json": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_programs/families_receipt.json": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_programs/receipt.json": ("NEW_CANDIDATE", 'native_stored_programs_and_accumulator'),
    "code/source_native_updates/CONTRACT.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_native_updates/README.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_native_updates/controls.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_native_updates/receipt.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_operator_join/README.md": ("NEW_CANDIDATE", 'finite_source_operator_join'),
    "code/source_operator_join/operator_join_packet.json": ("NEW_CANDIDATE", 'finite_source_operator_join'),
    "code/source_passive_memory/CONTRACT.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_passive_memory/README.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_passive_memory/controls.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_passive_memory/receipt.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_read_acceptance/CONTRACT.md": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_read_acceptance/README.md": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_read_acceptance/controls.json": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_read_acceptance/receipt.json": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_read_routing/CONTRACT.md": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_read_routing/README.md": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_read_routing/controls/q3_baseline.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_read_routing/controls/q3_branch.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_read_routing/controls/q3_scratch.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_read_routing/controls/q3_source.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_read_routing/specification.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_read_selection/CONTRACT.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_selected_pair_mean_histories'),
    "code/source_read_selection/README.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_selected_pair_mean_histories'),
    "code/source_read_selection/controls.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_selected_pair_mean_histories'),
    "code/source_read_selection/receipt.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_selected_pair_mean_histories'),
    "code/source_reusable_bus/CONTRACT.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_reusable_bus/README.md": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_reusable_bus/controls.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_reusable_bus/receipt.json": ("EXTENDS_EXISTING_CANDIDATE", 'lean_pair_mean_memory_and_reusable_bus'),
    "code/source_routing_refinement/CONTRACT.md": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_routing_refinement/README.md": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_routing_refinement/receipt.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "code/source_scalar_finite_instrument/README.md": ("EXTENDS_EXISTING_CANDIDATE", 'source_scalar_sequential_instrument'),
    "code/source_scalar_finite_instrument/finite_instrument_receipt.json": ("EXTENDS_EXISTING_CANDIDATE", 'source_scalar_sequential_instrument'),
    "code/source_temporal_acceptance/CONTRACT.md": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_temporal_acceptance/README.md": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_temporal_acceptance/controls.json": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "code/source_temporal_acceptance/receipt.json": ("NEW_CANDIDATE", 'native_temporal_tomography_and_checkpoint_selection'),
    "evidence/source_net_causal_poset/routed_read_law/q13_baseline.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "evidence/source_net_causal_poset/routed_read_law/q13_source.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "evidence/source_net_causal_poset/routed_read_law/q21_baseline.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
    "evidence/source_net_causal_poset/routed_read_law/q21_source.json": ("NEW_CANDIDATE", 'source_read_routing_full_family'),
}


EXPECTED_REVIEW_543298E0 = {'Lean/Geometry/SourceOperationReads.lean': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'Lean/Geometry/SourceOperationReadsAxiomAudit.lean': ('NEW_CANDIDATE',
                                                       'registered_operation_read_histories'),
 'Lean/Geometry/SourcePublicationAxiomAudit.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourcePublicationLaw.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourcePublicationMass.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourcePublicationMenu.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourcePublicationNative.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourcePublicationPath.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourcePublicationPathInvariant.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourcePublicationTail.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourceRadiusRemoval.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourceRadiusSelection.lean': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'Lean/Geometry/SourceRepairInstrument.lean': ('NEW_CANDIDATE', 'declared_flagged_recovery_extensions'),
 'Lean/Geometry/SourceStateSelection.lean': ('NEW_CANDIDATE', 'finite_state_process_selection'),
 'code/source_operation_reads/DERIVATION.md': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/INSTRUMENT_DERIVATION.md': ('NEW_CANDIDATE',
                                                          'declared_flagged_recovery_extensions'),
 'code/source_operation_reads/README.md': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/__init__.py': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/capture.json': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/check_live.py': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/receipt.json': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/simulator_verifier.py': ('NEW_CANDIDATE',
                                                       'registered_operation_read_histories'),
 'code/source_operation_reads/source_snapshot.json': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/test_reads.py': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_operation_reads/verify.py': ('NEW_CANDIDATE', 'registered_operation_read_histories'),
 'code/source_publication_selection/CONTRACT.md': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/DERIVATION.md': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/README.md': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/STATE_SELECTION.md': ('NEW_CANDIDATE', 'finite_state_process_selection'),
 'code/source_publication_selection/__init__.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/build.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/check_exterior.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/check_locality.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/check_quantum.py': ('NEW_CANDIDATE', 'finite_state_process_selection'),
 'code/source_publication_selection/codec.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/controls.json': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/exterior.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/locality.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/quantum.py': ('NEW_CANDIDATE', 'finite_state_process_selection'),
 'code/source_publication_selection/receipt.json': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/test_quantum.py': ('NEW_CANDIDATE', 'finite_state_process_selection'),
 'code/source_publication_selection/test_selection.py': ('NEW_CANDIDATE', 'native_eventual_publication'),
 'code/source_publication_selection/verify.py': ('NEW_CANDIDATE', 'native_eventual_publication')}


EXPECTED_PRIOR_INTEGRATED_REVIEWS = (
    (
        "2b03a95caf5030272f7b426b964e816f650153f8",
        "2d9bd11bc47c56d88a2fbbca22e3cc1be171d3f9",
        {
        "Lean/Geometry/FlatDiamondError.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/FlatDiamondNormalization.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/FlatDiamondPairIntegral.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/FlatDiamondVolume.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/FlatLorentzVolume.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/GoldenSourceCausalLimit.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/GoldenSourcePairLimit.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/GoldenSourceVolumeLimit.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/OrderingFractionFourDimensional.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/SourceBusScaling.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceCausalBoundary.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/SourceConstrainedRead.lean": ("NEW_CANDIDATE", "lean_selected_pair_mean_histories"),
        "Lean/Geometry/SourceConstrainedSelection.lean": ("NEW_CANDIDATE", "lean_selected_pair_mean_histories"),
        "Lean/Geometry/SourceCountLimitAxiomAudit.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/SourceCountTransport.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/SourceEncodedMemory.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceEncodedMemoryAxiomAudit.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceNativeRecords.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceNativeUpdatesAxiomAudit.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceNetOrderLimit.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/SourceNetVolumeError.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Geometry/SourceNonlinearRecord.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourcePassiveMemoryBudget.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourcePassiveReset.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceReadSelection.lean": ("NEW_CANDIDATE", "lean_selected_pair_mean_histories"),
        "Lean/Geometry/SourceReadSelectionAxiomAudit.lean": ("NEW_CANDIDATE", "lean_selected_pair_mean_histories"),
        "Lean/Geometry/SourceReusableBus.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceReusableBusAxiomAudit.lean": ("NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus"),
        "Lean/Geometry/SourceSelectionControls.lean": ("NEW_CANDIDATE", "lean_selected_pair_mean_histories"),
        "Lean/Screen/OPHScreen.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyConeMass.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyConeMassNaturality.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyConeModes.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyConeNaturality.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyFrameMorphism.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyFrameNaturality.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyMassPremiseInstance.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyNormalModeConstruction.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        "Lean/Screen/WhitneyRotationWitness.lean": ("OUTSIDE_REVIEWED_SOURCE_SCOPE", None),
        },
    ),
    (
        "2d9bd11bc47c56d88a2fbbca22e3cc1be171d3f9",
        "afed734528edff214c34d4038a64310544df922d",
        EXPECTED_PRIOR_BASE_REVIEW,
    ),
    (
        "afed734528edff214c34d4038a64310544df922d",
        "543298e06ce41d47f44c7707fdb2c32723cca10e",
        EXPECTED_REVIEW_543298E0,
    ),
)

EXPECTED_INTEGRATED_REVIEW = {
    'Lean/Geometry/RecordGluingPrinciple.lean': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'Lean/Geometry/RecordGluingPrincipleAxiomAudit.lean': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'Lean/Geometry/SourceRecordGluing.lean': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'Lean/Geometry/SourceRecordGluingAxiomAudit.lean': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'Lean/Geometry/SourceSelectionLocality.lean': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'Lean/Geometry/SourceSelectionLocalityAxiomAudit.lean': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/CONTRACT.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/DERIVATION.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/README.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/__init__.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/certificates.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/checker.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/compiler.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/receipt.json': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/receipt.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/rg_principle/test_rg_principle.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/sm_fermion_current/test_spin_exchange_integration.py': ('EXTENDS_EXISTING_CANDIDATE', 'fermionic_hypercharge_current_histories'),
    'code/source_record_gluing/README.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/__init__.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/build.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/capture.json': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/check_process.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/process.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/receipt.json': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/test_process.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_record_gluing/verify.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/CONTRACT.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/DERIVATION.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/README.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/RECORD_GLUING.md': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/__init__.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/channels.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/finite_model.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/geometry.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/grammar.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/receipt.json': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/records.json': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/records.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/response.json': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/response.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/test_model.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/verify.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
    'code/source_selection_model/verify_response.py': ('OUTSIDE_REVIEWED_SOURCE_SCOPE', None),
}


class VerificationError(ValueError):
    def __init__(self, message: str, code: str = "VERIFICATION_ERROR") -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def check_typed(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise VerificationError(message, code)


def strict_load(path: Path) -> dict[str, Any]:
    def object_pairs(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            check(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise VerificationError(f"non-finite JSON constant: {value}")

    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=object_pairs,
        parse_constant=reject_constant,
    )
    check(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot_exclusion_reason(relative: str) -> str | None:
    path = Path(relative)
    if relative in CONTENT_SNAPSHOT_EXACT_EXCLUSIONS:
        return "self_referential_inventory_manifest"
    if any(part in CONTENT_SNAPSHOT_EXCLUDED_DIRECTORY_NAMES for part in path.parts):
        return "repository_cache_directory"
    if path.suffix in CONTENT_SNAPSHOT_EXCLUDED_SUFFIXES:
        return "repository_cache_or_build_output"
    return None


def independently_snapshot_audited_files(repo_root: Path) -> list[dict[str, Any]]:
    paths: set[str] = set()
    for relative in AUDITED_DIRECTORIES:
        directory = repo_root / relative
        check(directory.is_dir(), f"missing audited directory: {relative}")
        for path in directory.rglob("*"):
            if not path.is_file():
                continue
            repo_relative = path.relative_to(repo_root).as_posix()
            if snapshot_exclusion_reason(repo_relative) is not None:
                continue
            paths.add(repo_relative)
    return [
        {
            "path": relative,
            "byte_count": (repo_root / relative).stat().st_size,
            "sha256": file_sha256(repo_root / relative),
        }
        for relative in sorted(paths)
    ]


def verify_audited_file_snapshot(snapshot: Any, repo_root: Path) -> set[str]:
    check(isinstance(snapshot, Mapping), "audited file snapshot missing")
    check(snapshot.get("directories") == list(AUDITED_DIRECTORIES), "audit directory drift")
    expected_policy = {
        "exact_paths": list(CONTENT_SNAPSHOT_EXACT_EXCLUSIONS),
        "directory_names": list(CONTENT_SNAPSHOT_EXCLUDED_DIRECTORY_NAMES),
        "suffixes": list(CONTENT_SNAPSHOT_EXCLUDED_SUFFIXES),
        "reasons": {
            "exact_paths": "the generated inventory would otherwise contain its own byte digest",
            "directory_names": "repository cache directories cannot supply reviewed scientific source",
            "suffixes": "interpreter cache and build outputs cannot supply reviewed scientific source",
        },
    }
    check(snapshot.get("exclusion_policy") == expected_policy, "audit exclusion policy drift")
    files = snapshot.get("files")
    check(
        isinstance(files, list)
        and all(
            isinstance(row, Mapping)
            and set(row) == {"path", "byte_count", "sha256"}
            and isinstance(row.get("path"), str)
            and bool(row.get("path"))
            and isinstance(row.get("byte_count"), int)
            and row.get("byte_count") >= 0
            and isinstance(row.get("sha256"), str)
            for row in files
        ),
        "audited file content rows",
    )
    committed_paths = [row["path"] for row in files]
    check(committed_paths == sorted(set(committed_paths)), "audited file paths not canonical")
    check(snapshot.get("file_count") == len(files), "audited file count")
    check(snapshot.get("files_sha256") == canonical_sha256(files), "audited content-list hash")
    current_files = independently_snapshot_audited_files(repo_root)
    current_paths = [row["path"] for row in current_files]
    check_typed(
        committed_paths == current_paths,
        "AUDITED_PATH_DRIFT",
        "audited file snapshot has added, removed, or renamed paths",
    )
    for committed, current in zip(files, current_files, strict=True):
        check_typed(
            committed["byte_count"] == current["byte_count"],
            "AUDITED_CONTENT_DRIFT",
            f"byte-count drift: {committed['path']}",
        )
        check_typed(
            committed["sha256"] == current["sha256"],
            "AUDITED_CONTENT_DRIFT",
            f"digest drift: {committed['path']}",
        )
    return set(committed_paths)


def verify_pin_rows(rows: Any, repo_root: Path, label: str) -> set[str]:
    check(isinstance(rows, list) and rows, f"{label} pins missing")
    observed: set[str] = set()
    for row in rows:
        check(isinstance(row, Mapping), f"{label} pin row")
        relative = row.get("path")
        check(isinstance(relative, str) and relative, f"{label} pin path")
        path = Path(relative)
        check(not path.is_absolute() and ".." not in path.parts, f"unsafe pin path: {relative}")
        check(relative not in observed, f"duplicate pin path: {relative}")
        observed.add(relative)
        absolute = repo_root / path
        check(absolute.is_file(), f"missing pinned file: {relative}")
        check(row.get("bytes") == absolute.stat().st_size, f"byte pin drift: {relative}")
        check(row.get("sha256") == file_sha256(absolute), f"SHA pin drift: {relative}")
    return observed


def verify_import_firewall(producer_path: Path) -> None:
    tree = ast.parse(producer_path.read_text(encoding="utf-8"), filename=str(producer_path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {
        "__future__",
        "argparse",
        "collections",
        "hashlib",
        "json",
        "pathlib",
        "typing",
    }
    check(imports <= allowed, f"producer import firewall failed: {sorted(imports - allowed)}")
    check(
        "port_current_inner_certificate" not in imports,
        "conditional current implementation imported by producer",
    )


def independently_check_sources(repo_root: Path) -> dict[str, int]:
    projection = strict_load(
        repo_root / "code/a5_closure/manifests/source_current_capability_projection.json"
    )
    capability = strict_load(
        repo_root / "code/a5_closure/receipts/source_current_capability.receipt.json"
    )
    recurrence = projection.get("registered_recurrence", {})
    algebra = capability.get("audit", {}).get("response_word_algebra", {})
    rechart = capability.get("audit", {}).get("recharting_audit", {})
    check(recurrence.get("raw_runtime_history_serialized") is False, "raw-history scope")
    check(recurrence.get("source_count") == 12, "source count")
    check(algebra.get("exact_dimension") == 4, "response algebra dimension")
    check(algebra.get("commutative") is True, "response algebra commutativity")
    check(algebra.get("commutator_nonzero_count") == 0, "response commutator count")
    check(rechart.get("proper_recharting_count") == 60, "proper recharting count")
    check(
        rechart.get("proper_rechartings_in_response_word_algebra") == 1,
        "response/recharting intersection",
    )

    routing = strict_load(repo_root / "code/source_routing/runtime/path_tomography_receipt.json")
    routing_ops = {
        event.get("op")
        for episode in routing.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    check("pair_mean" in routing_ops, "path-tomography pair mean")
    check(
        routing.get("scope", {}).get("full_metric_neighbor_refinement") is False,
        "path-tomography refinement scope",
    )

    feedback = strict_load(repo_root / "code/source_feedback_transport/transport_receipt.json")
    feedback_ops = {
        event.get("op")
        for episode in feedback.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    check({"mean", "reset", "export"} <= feedback_ops, "feedback destructive operations")

    scalar = strict_load(
        repo_root / "code/source_scalar_instruments/sequential_instrument_receipt.json"
    )
    check(
        scalar.get("scope", {}).get("sampled_quantum_outcome_histories") is False,
        "scalar sampled-history scope",
    )
    check(
        scalar.get("scope", {}).get("source_selected_quantum_operations") is False,
        "scalar source-selection scope",
    )

    counting = strict_load(
        repo_root / "code/a5_closure/manifests/record_counting_mechanism_reference.json"
    )
    theorem = counting.get("a2_clauses", {}).get("presentation_invariance", {}).get(
        "theorem", ""
    )
    check("addition is commutative" in theorem, "record-counting order erasure")
    check(
        "append-only" in counting.get("mechanism", {}).get("event_grammar", ""),
        "record-counting append-only grammar",
    )

    directed = strict_load(
        repo_root / "code/a5_closure/manifests/directed_seam_repair_reference.json"
    )
    check(
        directed.get("verdict", {}).get("refinement_semigroup_naturality") == "open",
        "directed repair refinement boundary",
    )

    bracket1 = strict_load(
        repo_root
        / "code/a5_closure/issue_566_bracket_space_stage1/"
        "a5_alternating_bracket_space_stage1.receipt.json"
    )
    bracket2 = strict_load(
        repo_root
        / "code/a5_closure/issue_566_bracket_space_stage2/"
        "a5_jacobi_stage2.receipt.json"
    )
    check(not any(bracket1.get("later_gates", {}).values()), "stage1 later gates")
    check(
        bracket2.get("not_attained", {}).get("preferred_bracket_selected") is False,
        "stage2 bracket selection boundary",
    )

    b14 = strict_load(
        repo_root / "code/b14_jacobi/oriented_face_bracket_selector.certificate.json"
    )
    check(b14.get("jacobi_failure", {}).get("nonzero_count") == 240, "B14 Jacobi count")
    check(b14.get("jacobi_failure", {}).get("squared_norm") == 240, "B14 Jacobi norm")

    fixture = strict_load(
        repo_root / "code/a5_closure/manifests/port_current_response_reference.json"
    )
    check(fixture.get("construction_model") == "charged_double_triplet", "fixture typing")
    kinetic = strict_load(
        repo_root / "code/angular_sprint/runtime/kinetic_form_selection_receipt.json"
    )
    check(
        kinetic.get("physical_selection_boundary", {}).get("current_lift_source_selected")
        is False,
        "kinetic current source boundary",
    )

    source_history_text = (
        repo_root / "Lean/Variational/SourceToHamiltonianComposed.lean"
    ).read_text(encoding="utf-8")
    check(
        "theorem here claims that the source selects" in source_history_text,
        "Lean source-selection boundary",
    )
    response_text = (repo_root / "Lean/Screen/A5ResponseWordAlgebra.lean").read_text(
        encoding="utf-8"
    )
    check("commutative" in response_text.lower(), "Lean response-word boundary")
    a2_words = (
        repo_root / "Lean/ObserverPatchHolography/A2EndpointCommutator.lean"
    ).read_text(encoding="utf-8")
    check("FifteenCommute" in a2_words, "Lean conditional A2 word boundary")
    provenance = (
        repo_root
        / "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean"
    ).read_text(encoding="utf-8")
    check("two serializations of the diamond" in provenance, "Lean serialization boundary")

    encoded_memory = (repo_root / "Lean/Geometry/SourceEncodedMemory.lean").read_text(
        encoding="utf-8"
    )
    check("Copying halves both" in encoded_memory, "encoded-memory attenuation boundary")
    passive_reset = (repo_root / "Lean/Geometry/SourcePassiveReset.lean").read_text(
        encoding="utf-8"
    )
    check("positive_register_cannot_reset" in passive_reset, "passive-reset boundary")
    passive_budget = (
        repo_root / "Lean/Geometry/SourcePassiveMemoryBudget.lean"
    ).read_text(encoding="utf-8")
    check("closed_cycle_quiescent" in passive_budget, "closed-cycle boundary")
    native_records = (repo_root / "Lean/Geometry/SourceNativeRecords.lean").read_text(
        encoding="utf-8"
    )
    check("schedule-selection obstruction" in native_records, "schedule-order boundary")
    selected_histories = (
        repo_root / "Lean/Geometry/SourceReadSelection.lean"
    ).read_text(encoding="utf-8")
    check(
        "not impossibility of eventual" in selected_histories,
        "selected-history scope boundary",
    )

    routing_spec = strict_load(repo_root / "code/source_read_routing/specification.json")
    routing_scope = routing_spec.get("scope", {})
    check(routing_scope.get("full_q13_q21_routed_execution") is True, "routing execution scope")
    check(routing_scope.get("axiomatic_read_law_derived") is False, "routing read-law derivation boundary")
    check(
        routing_spec.get("M1", {}).get("status") == "retained supplied structural rule",
        "routing M1 status",
    )
    routing_events: dict[str, int] = {}
    for name in ("q13_baseline", "q13_source", "q21_baseline", "q21_source"):
        run = strict_load(
            repo_root / f"evidence/source_net_causal_poset/routed_read_law/{name}.json"
        )
        check(run.get("schema") == "oph.source_read_routing.run.v1", f"routing run schema: {name}")
        check(run.get("variant") == name.split("_")[1], f"routing run variant: {name}")
        routing_events[name] = int(run.get("costs", {}).get("events"))
    check(routing_events["q13_baseline"] == routing_events["q13_source"], "q13 event census")
    check(routing_events["q21_baseline"] == routing_events["q21_source"], "q21 event census")
    routing_lean = (repo_root / "Lean/Geometry/SourceReadRouting.lean").read_text(encoding="utf-8")
    check("resetReceiver" in routing_lean and "resetSource" in routing_lean, "routing hop resets")
    check("closed_sum_preserving_word_cannot_reset" in routing_lean, "routing reset obstruction")
    check("schedule_independent_readouts" in routing_lean, "routing schedule independence")
    refinement = strict_load(repo_root / "code/source_routing_refinement/receipt.json")
    check(refinement.get("scope", {}).get("M1_source_selected") is False, "refinement M1 boundary")

    programs = strict_load(repo_root / "code/source_native_programs/receipt.json")
    accumulator = strict_load(repo_root / "code/source_native_accumulator/receipt.json")
    check(programs.get("dynamic_writes") == "supported_scalar_pair_means_only", "stored-program writes")
    check(programs.get("m1_derived") is False, "stored-program M1 boundary")
    check(accumulator.get("dynamic_scalar_writes") == "pair_means_only", "accumulator writes")
    check(accumulator.get("m1_derived") is False, "accumulator M1 boundary")
    bus_witness = (repo_root / "Lean/Geometry/SourceBankBusWitness.lean").read_text(encoding="utf-8")
    check("declared graph data" in bus_witness, "bank bus W12 embedding boundary")

    temporal_controls = strict_load(repo_root / "code/source_temporal_acceptance/controls.json")
    captured_plan = temporal_controls.get("tomography", {}).get("captured_plan", {})
    check(temporal_controls.get("full_axiom_instantiation") is False, "temporal axiom boundary")
    check(captured_plan.get("payload_replay") is False, "temporal plan payload replay")
    read_acceptance = strict_load(repo_root / "code/source_read_acceptance/receipt.json")
    check(
        read_acceptance.get("source_contract") == "conditional_not_source_selected",
        "read acceptance source contract",
    )
    checkpoint = strict_load(repo_root / "code/source_checkpoint_selection/controls.json")
    check(
        checkpoint.get("scope", {}).get("canonical_axiom_selection") is False,
        "checkpoint axiom selection boundary",
    )
    erasure = (repo_root / "Lean/Geometry/SourceTemporalErasure.lean").read_text(encoding="utf-8")
    check(
        "Permanent loss before the first discriminating observation" in erasure,
        "temporal erasure boundary",
    )
    tomography = (repo_root / "Lean/Geometry/SourceTemporalTomography.lean").read_text(
        encoding="utf-8"
    )
    check("no metric menu is selected" in tomography, "temporal tomography menu boundary")

    read_selection_controls = strict_load(repo_root / "code/source_read_selection/controls.json")
    check(read_selection_controls.get("source_ports") == list(range(12)), "read selection ports")
    reusable_bus = strict_load(repo_root / "code/source_reusable_bus/receipt.json")
    check(len(reusable_bus.get("histories", [])) == 16, "reusable-bus tape census")

    finite_instrument = strict_load(
        repo_root / "code/source_scalar_finite_instrument/finite_instrument_receipt.json"
    )
    check(
        "code/source_scalar_instruments/sequential_instrument_receipt.json"
        in str(finite_instrument.get("parents")),
        "finite instrument parent pin",
    )

    operator_join = strict_load(repo_root / "code/source_operator_join/operator_join_packet.json")
    join_interpretation = operator_join.get("interpretation", {})
    check(
        join_interpretation.get("overlapping_pair_algebras_commute") is False,
        "operator join hinge",
    )
    check(join_interpretation.get("tensor_assembly") == "declared", "operator join assembly")
    check(join_interpretation.get("observed_quantum_outcomes") is False, "operator join outcomes")

    fermion = strict_load(repo_root / "code/sm_fermion_current/current_receipt.json")
    fermion_scope = fermion.get("scope", {})
    check(
        fermion_scope.get("hypercharge_edge_factors_on_prepared_q5_sites") is True,
        "fermion current label",
    )
    check(
        fermion_scope.get("source_selected_action_population_or_clock") is False,
        "fermion source-selection boundary",
    )
    maxwell_readme = (repo_root / "code/maxwell_measurement/README.md").read_text(encoding="utf-8")
    check(
        "No apparatus, physical raw data, frozen experiment or observed result is shipped."
        in maxwell_readme,
        "Maxwell adapter shipping boundary",
    )

    publication = strict_load(repo_root / "code/source_publication_selection/receipt.json")
    operation_capture = strict_load(repo_root / "code/source_operation_reads/capture.json")
    operation_receipt = strict_load(repo_root / "code/source_operation_reads/receipt.json")
    check_typed(publication.get("M1_derived") is False and publication.get("complete_prefix_consistency") is True,
            "PUBLICATION_BOUNDARY", "conditional publication boundary drift")
    check_typed(publication.get("noncommuting_segments_checked") == 3 and publication.get("quantum_channels_reconstructed") == 4,
            "STATE_SELECTION_CONTROLS", "finite state/process controls drift")
    check_typed(len(operation_capture.get("cases", [])) == 7 and operation_capture.get("scope", {}).get("complete_A1_A3") is False,
            "NATIVE_READ_SCOPE", "native read capture scope drift")
    check_typed(operation_capture.get("instruments", {}).get("extensions_selected_by_source") is False,
            "INSTRUMENT_SCOPE", "supplied instrument promoted to native source")
    check_typed(operation_receipt.get("M1_derived") is False and operation_receipt.get("complete_A1_A3_model") is False,
            "NATIVE_READ_BOUNDARY", "native read physical scope drift")

    return {
        "response_algebra_dimension": int(algebra["exact_dimension"]),
        "response_commutator_nonzero_count": int(algebra["commutator_nonzero_count"]),
        "proper_recharting_count": int(rechart["proper_recharting_count"]),
        "b14_jacobi_nonzero_count": int(b14["jacobi_failure"]["nonzero_count"]),
        "routing_q13_events": routing_events["q13_baseline"],
        "routing_q21_events": routing_events["q21_baseline"],
    }


def verify(inventory_path: Path, repo_root: Path = REPO_ROOT) -> dict[str, Any]:
    inventory = strict_load(inventory_path)
    check(inventory.get("schema") == SCHEMA, "inventory schema")
    body = {key: value for key, value in inventory.items() if key != "inventory_sha256"}
    check(inventory.get("inventory_sha256") == canonical_sha256(body), "inventory hash")
    check(inventory.get("verdict") == VERDICT, "inventory verdict")
    check(
        inventory.get("previous_inventory_base_sha") == PREVIOUS_INVENTORY_BASE_SHA,
        "previous inventory base",
    )
    check(inventory.get("audited_upstream_main_sha") == UPSTREAM_MAIN_SHA, "audited base")
    check(
        "explicitly enumerated candidates and reviewed source scope"
        in inventory.get("verdict_scope", ""),
        "bounded verdict scope",
    )
    check(inventory.get("positive_stages_authorized") is False, "positive-stage stop")
    check(
        inventory.get("physical_current_source_bridge_attained") is False,
        "false physical source gate",
    )
    check(
        inventory.get("fixture_used_as_reconstruction_oracle") is False,
        "fixture oracle firewall",
    )
    check(inventory.get("target_labels_used") is False, "target-label firewall")
    check(
        inventory.get("qualification_fields") == list(QUALIFICATION_FIELDS),
        "qualification field schema",
    )

    snapshotted_paths = verify_audited_file_snapshot(
        inventory.get("audited_file_snapshot"), repo_root
    )
    extension = inventory.get("audit_scope", {}).get("directories_added_since_initial_inventory")
    check(isinstance(extension, list), "audit scope extension missing")
    check(
        [row.get("directory") for row in extension] == list(EXPECTED_SCOPE_EXTENSION),
        "audit scope extension drift",
    )
    check(
        all(isinstance(row.get("reason"), str) and row.get("reason") for row in extension),
        "audit scope extension reasons",
    )
    check(
        all(
            any(path.startswith(directory + "/") for path in snapshotted_paths)
            for directory in EXPECTED_SCOPE_EXTENSION
        ),
        "audit scope extension not bound by the content snapshot",
    )

    integrated_review = inventory.get("integrated_tree_review")
    check(isinstance(integrated_review, Mapping), "integrated-tree review missing")
    check(
        integrated_review.get("from_upstream_main_sha") == PREVIOUS_INVENTORY_BASE_SHA,
        "integrated-tree review start",
    )
    check(
        integrated_review.get("through_upstream_main_sha") == UPSTREAM_MAIN_SHA,
        "integrated-tree review end",
    )
    surfaces = integrated_review.get("surfaces")
    check(isinstance(surfaces, list), "integrated-tree review surfaces")
    check(
        integrated_review.get("surface_count") == len(EXPECTED_INTEGRATED_REVIEW),
        "integrated-tree review count",
    )
    observed_review: dict[str, tuple[Any, Any]] = {}
    for row in surfaces:
        check(isinstance(row, Mapping), "integrated-tree review row")
        path = row.get("path")
        check(isinstance(path, str) and path, "integrated-tree review path")
        check(path not in observed_review, f"duplicate integrated-tree review: {path}")
        check(
            isinstance(row.get("reason"), str) and bool(row.get("reason")),
            f"integrated-tree review reason: {path}",
        )
        observed_review[path] = (row.get("decision"), row.get("candidate_id"))
    check(
        list(observed_review) == sorted(observed_review),
        "integrated-tree review order",
    )
    check(observed_review == EXPECTED_INTEGRATED_REVIEW, "integrated-tree semantic decisions")
    check(set(observed_review) <= snapshotted_paths, "reviewed surface missing from snapshot")
    prior_reviews = integrated_review.get("prior_reviews")
    check(isinstance(prior_reviews, list), "prior integrated-tree reviews")
    check(len(prior_reviews) == len(EXPECTED_PRIOR_INTEGRATED_REVIEWS), "prior review count")
    for prior, (start, end, expected_map) in zip(
        prior_reviews, EXPECTED_PRIOR_INTEGRATED_REVIEWS, strict=True
    ):
        check(isinstance(prior, Mapping), "prior review row")
        check(prior.get("from_upstream_main_sha") == start, "prior review start")
        check(prior.get("through_upstream_main_sha") == end, "prior review end")
        prior_surfaces = prior.get("surfaces")
        check(isinstance(prior_surfaces, list), "prior review surfaces")
        check(prior.get("surface_count") == len(expected_map), "prior review count")
        observed_prior = {
            row.get("path"): (row.get("decision"), row.get("candidate_id"))
            for row in prior_surfaces
            if isinstance(row, Mapping)
        }
        check(observed_prior == expected_map, "prior integrated-tree semantic decisions")
        check(set(observed_prior) <= snapshotted_paths, "prior reviewed surface missing from snapshot")

    source_paths = verify_pin_rows(inventory.get("source_pins"), repo_root, "source")
    implementation_paths = verify_pin_rows(
        inventory.get("implementation_pins"), repo_root, "implementation"
    )
    check(
        implementation_paths
        == {
            "code/a5_closure/source_current_order_sensitive_inventory.py",
            "code/a5_closure/verify_source_current_order_sensitive_inventory.py",
        },
        "implementation pin set",
    )
    check(
        {
            path
            for path in source_paths | implementation_paths
            if any(path == root or path.startswith(root + "/") for root in AUDITED_DIRECTORIES)
        }
        <= snapshotted_paths,
        "pinned audited file missing from snapshot",
    )
    verify_import_firewall(repo_root / "code/a5_closure/source_current_order_sensitive_inventory.py")

    candidates = inventory.get("candidates")
    check(isinstance(candidates, list), "candidate list")
    candidate_map: dict[str, Mapping[str, Any]] = {}
    for row in candidates:
        check(isinstance(row, Mapping), "candidate row")
        candidate_id = row.get("candidate_id")
        check(isinstance(candidate_id, str), "candidate id")
        check(candidate_id not in candidate_map, f"duplicate candidate: {candidate_id}")
        candidate_map[candidate_id] = row
    check(set(candidate_map) == set(EXPECTED_CLASSIFICATIONS), "candidate coverage")

    qualifying: list[str] = []
    for candidate_id, expected_classification in EXPECTED_CLASSIFICATIONS.items():
        row = candidate_map[candidate_id]
        check(row.get("classification") == expected_classification, f"classification: {candidate_id}")
        props = row.get("properties")
        check(props == EXPECTED_PROPERTIES[candidate_id], f"properties: {candidate_id}")
        failures = [field for field in QUALIFICATION_FIELDS if not props[field]]
        check(row.get("qualification_failures") == failures, f"failure list: {candidate_id}")
        check(
            isinstance(row.get("rejection_reason"), str) and row.get("rejection_reason"),
            f"rejection reason: {candidate_id}",
        )
        row_paths = row.get("source_paths")
        check(isinstance(row_paths, list) and row_paths, f"source paths: {candidate_id}")
        check(set(row_paths) <= source_paths, f"unpinned source path: {candidate_id}")
        if all(props[field] for field in QUALIFICATION_FIELDS):
            qualifying.append(candidate_id)
    check(not qualifying, f"unexpected qualifying candidates: {qualifying}")

    counts = Counter(EXPECTED_CLASSIFICATIONS.values())
    summary = inventory.get("summary", {})
    check(summary.get("candidate_count") == len(EXPECTED_CLASSIFICATIONS), "candidate count")
    check(summary.get("qualifying_candidate_count") == 0, "qualifying count")
    check(summary.get("classification_counts") == dict(sorted(counts.items())), "class counts")
    adjacency = summary.get("registered_adjacency_response_summary", {})
    check(adjacency.get("candidate_id") == "registered_adjacency_response", "adjacency-response summary")
    check(adjacency.get("response_word_dimension") == 4, "reported response dimension")
    check(adjacency.get("commutator_nonzero_count") == 0, "reported commutator count")
    check(adjacency.get("proper_rechartings_total") == 60, "reported rechartings")
    check(adjacency.get("proper_rechartings_in_response_words") == 1, "reported intersection")
    check(
        summary.get("logged_repair_candidate_summary", {}).get("candidate_id")
        == "record_counting_repair_628",
        "logged repair summary",
    )
    integrated_near = summary.get("routed_read_candidate_summary", {})
    check(
        integrated_near.get("candidate_id") == "source_read_routing_full_family",
        "routed-read summary",
    )
    check(
        isinstance(integrated_near.get("fatal_failures"), list)
        and len(integrated_near["fatal_failures"]) == 3,
        "integrated near-candidate fatal failures",
    )
    check("different primitives" in summary.get("cross_packet_noncomposition", ""), "no compositing rule")
    check("routing packet" in summary.get("cross_packet_noncomposition", ""), "routing packet non-compositing")

    missing = inventory.get("minimal_missing_fields")
    check(isinstance(missing, list), "minimal missing fields")
    check(
        [row.get("field") for row in missing]
        == [
            "perturbation_families",
            "ordered_mixed_histories",
            "raw_history_custody",
            "common_source_binding",
            "reversibility_witness",
            "source_firewall",
        ],
        "minimal missing-field set",
    )
    boundary = inventory.get("claim_boundary", {})
    required_false_promotions = (
        "bounded_repository_inventory_only",
        "verdict_limited_to_enumerated_reviewed_source_scope",
        "does_not_assert_a_no_go_for_future_source_data",
        "does_not_weaken_the_abstract_forced_lie_type_theorem",
        "does_not_reject_the_conditional_port_current_fixture",
        "does_not_make_a_physical_current_claim",
        "stages_3_through_6_not_run",
    )
    check(all(boundary.get(key) is True for key in required_false_promotions), "claim boundary")

    independent = independently_check_sources(repo_root)
    check(independent["response_algebra_dimension"] == 4, "independent dimension")
    return {
        "verdict": VERDICT,
        "candidate_count": len(candidate_map),
        "qualifying_candidate_count": len(qualifying),
        "classification_counts": dict(sorted(counts.items())),
        "response_algebra_dimension": independent["response_algebra_dimension"],
        "proper_recharting_count": independent["proper_recharting_count"],
        "b14_jacobi_nonzero_count": independent["b14_jacobi_nonzero_count"],
        "routing_q13_events": independent["routing_q13_events"],
        "routing_q21_events": independent["routing_q21_events"],
        "audited_file_count": len(snapshotted_paths),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args()
    print(json.dumps(verify(args.inventory, args.repo_root), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
