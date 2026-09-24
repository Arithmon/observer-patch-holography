#!/usr/bin/env python3
"""Stage-2 inventory for the source-current order-sensitive object.

This producer audits existing committed source surfaces.  It does not build a
new perturbation model and it never imports the conditional PORT-CURRENT-INNER
implementation.  A candidate qualifies only if one and the same source packet
contains twelve port-indexed reversible perturbation families, both mixed
composition orders, raw serialized histories, target-free provenance, and the
registered twelve-port refinement ancestry.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
INVENTORY_PATH = HERE / "manifests" / "source_current_order_sensitive_inventory.json"
INDEPENDENT_VERIFIER = HERE / "verify_source_current_order_sensitive_inventory.py"

SCHEMA = "oph.source_current_order_sensitive_inventory.v2"
VERDICT = "SOURCE_CURRENT_ORDER_SENSITIVE_OBJECT_NOT_PRESENT"
PREVIOUS_INVENTORY_BASE_SHA = 'ba84976ad1195892c68fdd0d74a27ab0c899aeb8'
UPSTREAM_MAIN_SHA = '9b527a4f8d07a21944b56bc3967a739339aff2c8'

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

AUDIT_SCOPE_EXTENSION = (
    ("code/source_read_routing",
     "runtime compiler and retained level-three tapes of the full-family routed read histories"),
    ("code/source_routing_refinement",
     "refinement companion of the routed read histories across support levels three to five"),
    ("code/source_native_programs",
     "runtime controls and receipts of the stored-program bank compiler formalized under Lean/Geometry"),
    ("code/source_native_accumulator",
     "runtime replay of native accumulation formalized under Lean/Geometry"),
    ("code/source_temporal_acceptance",
     "runtime controls of temporal acceptance and native tomography plans on the captured support"),
    ("code/source_read_acceptance",
     "runtime replay of sound partial readout on captured seams"),
    ("code/source_checkpoint_selection",
     "runtime census and policy controls of native checkpoint selection on captured seams"),
    ("code/source_read_selection",
     "runtime controls of the selected-history candidate already classified from its Lean modules"),
    ("code/source_encoded_memory",
     "runtime controls of the pair-mean memory candidate already classified from its Lean modules"),
    ("code/source_passive_memory",
     "runtime controls of the pair-mean memory candidate already classified from its Lean modules"),
    ("code/source_reusable_bus",
     "retained fixed-point tapes of the pair-mean memory candidate already classified from its Lean modules"),
    ("code/source_native_updates",
     "retained histories of the pair-mean memory candidate already classified from its Lean modules"),
    ("code/source_scalar_finite_instrument",
     "finite-duration extension of the sequential scalar instrument candidate"),
    ("code/source_operator_join",
     "operator-algebra join of committed carriers with a noncommuting hinge"),
    ("code/sm_fermion_current",
     "register-level histories that carry a physical current label and therefore test the source firewall"),
    ("code/maxwell_measurement",
     "measurement adapter that identifies registers with electromagnetic potentials and therefore tests the source firewall"),
    ("evidence/source_net_causal_poset/routed_read_law",
     "production run receipts and phase manifests of the routed read histories"),
    ("code/source_publication_selection", "Conditional publication histories and supplied state/process selection controls."),
    ("code/source_operation_reads", "Registered all-port histories and separately supplied instrument/recovery extensions."),
    ("code/source_selection_model",
     "declared scalar-seam response witness and record interventions whose explicit generator formula and algebraic jets test the source-native and raw-history requirements"),
    ("code/source_record_gluing",
     "reference execution of the proposed classical record-transport grammar with an additive diagnostic"),
    ("code/rg_principle",
     "local-clock reduction controls and finite classical compiler of the proposed record-transport law"),
)

CONTENT_SNAPSHOT_EXACT_EXCLUSIONS = (
    "code/a5_closure/manifests/source_current_order_sensitive_inventory.json",
)
CONTENT_SNAPSHOT_EXCLUDED_DIRECTORY_NAMES = ("__pycache__",)
CONTENT_SNAPSHOT_EXCLUDED_SUFFIXES = (".pyc",)

PRIOR_BASE_INTEGRATED_TREE_REVIEW = (
    ("Lean/Geometry/SourceAccumulatorAxiomAudit.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Transitive standard-axiom gate for the accumulator theorems; it adds no operation or history."),
    ("Lean/Geometry/SourceAccumulatorDecoder.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Bounded integer interval decoder for native accumulation publication; a readout policy, not a perturbation family."),
    ("Lean/Geometry/SourceAccumulatorProgram.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Finite accumulation programs with retirement and a whole-program error estimate over supplied request lists; pair means only."),
    ("Lean/Geometry/SourceBankBusWitness.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Declared finite bus discharging the compiler's route premises; its star connections are declared graph data without a W12 embedding."),
    ("Lean/Geometry/SourceBankCompiler.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "End-to-end compilation of finite bank programs from geometric route data; no perturbation family, inverse law or history custody."),
    ("Lean/Geometry/SourceBankControl.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Generated kernel check of one captured 42-instruction control trace; scalar-port embedding is not asserted."),
    ("Lean/Geometry/SourceBankExecution.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Finite noisy publication for bank programs on the declared bus; preparation and controller remain hypotheses."),
    ("Lean/Geometry/SourceBankInvariant.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Workspace, retirement and confinement invariants of the integer bank machine; no ordered perturbation data."),
    ("Lean/Geometry/SourceBankLowering.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Native lowering of bank instructions to retained-source shuttles with error allowances; retirement and cleanup are destructive."),
    ("Lean/Geometry/SourceBankMachine.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Integer semantics of a payload-independent instruction stream executed in one supplied request order per program."),
    ("Lean/Geometry/SourceBankPrecision.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Existence of finite precision for compiled words; no source operation is added."),
    ("Lean/Geometry/SourceBankTrace.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Row checker for the captured control trace; a finite comparison, not a runtime history producer."),
    ("Lean/Geometry/SourceCheckpointAxiomAudit.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Transitive standard-axiom gate for checkpoint selection."),
    ("Lean/Geometry/SourceCheckpointEntropy.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Unique information projection onto a checkpoint face for a named full-history cover; a selection theorem, not an operation family."),
    ("Lean/Geometry/SourceCheckpointMeaning.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Complete record test grammar and native continuation viability for a supplied record map."),
    ("Lean/Geometry/SourceCheckpointNative.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Composition of observability, finite feasibility and checkpoint selection; the record interface and deadline remain source obligations."),
    ("Lean/Geometry/SourceCheckpointPipeline.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Commutation of disjoint means and exact two-sample recovery on the shortest path pipeline; two records on one path, not mixed-order port families."),
    ("Lean/Geometry/SourceCheckpointPolicy.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Finite checkpoint transition law from backward masses; a policy on a supplied move grammar."),
    ("Lean/Geometry/SourceNativeAccumulator.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Reusable accumulation in which only scalar pair means write the state and the previous value is retired at each start."),
    ("Lean/Geometry/SourceNativeCore.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Reset and start words that clear the accumulator to a residual by finite pair means; clearing is irreversible."),
    ("Lean/Geometry/SourceNativeProgramBudget.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Precision margin and stationary repeated-word evaluation; arithmetic bounds only."),
    ("Lean/Geometry/SourceNativeProgramError.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Composition of native words with explicit residual bounds; comparison errors, not perturbation families."),
    ("Lean/Geometry/SourceNativeProgramsAxiomAudit.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Transitive standard-axiom gate for the stored-program theorems."),
    ("Lean/Geometry/SourceNativeShuttle.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Native transfer into a reusable record with bus cleanup by scalar means; the source is retained at half amplitude."),
    ("Lean/Geometry/SourceNativeStoredProgram.lean", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Composable stored computation from start, addition, retirement and transfer; neither the word nor the layout is source-selected."),
    ("Lean/Geometry/SourceReadAcceptance.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Sound partial publication from observational distinguishability and its native threshold realization; no source grammar is derived."),
    ("Lean/Geometry/SourceReadAcceptanceAxiomAudit.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Transitive standard-axiom gate for read acceptance."),
    ("Lean/Geometry/SourceReadAcceptanceSchedule.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Finite history counting and stopped native cost on a supplied alphabet; uniform counting laws."),
    ("Lean/Geometry/SourceReadRouting.lean", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Six-event transport word with constant resets and a destructive pair mean; readouts are proved schedule-independent and closed sum-preserving words cannot reset."),
    ("Lean/Geometry/SourceRoutingBudget.lean", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Event, read and hop counting bounds for the routing compiler; resource arithmetic without histories."),
    ("Lean/Geometry/SourceRoutingHierarchy.lean", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Walk diameters, refinement lifting and the base icosahedral diameter three; a support-distance bound, not a perturbation family."),
    ("Lean/Geometry/SourceRoutingRefinementAxiomAudit.lean", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Transitive standard-axiom gate for the routing refinement theorems."),
    ("Lean/Geometry/SourceRoutingStorage.lean", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Live logical versions in a 23-slot local store with a retirement discipline; no canonical pair-mean selection follows."),
    ("Lean/Geometry/SourceTemporalAcceptanceAxiomAudit.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Transitive standard-axiom gate for temporal acceptance and native completion."),
    ("Lean/Geometry/SourceTemporalAttempts.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Proposal-count failure bounds under a supplied uniform proposal law; counting, not histories."),
    ("Lean/Geometry/SourceTemporalConnected.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Existence of a complete calibration plan on finite connected support; exact readout existence, not schedule selection."),
    ("Lean/Geometry/SourceTemporalErasure.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Permanent loss of a prepared difference after its first mean; the separate and reversed words on one seam path are a two-order witness for one pair of unknowns."),
    ("Lean/Geometry/SourceTemporalEssential.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Essential coordinates of a finite product meaning; the minimal direct-read set."),
    ("Lean/Geometry/SourceTemporalGuard.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Maximal preservation guard with one universal completing word; proposals remain attempt history under a separate law."),
    ("Lean/Geometry/SourceTemporalMenu.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Essential direct reads versus aggregate observability for a supplied linear meaning."),
    ("Lean/Geometry/SourceTemporalNative.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Linear response forms of pair-mean words under a declared preparation; prefixes retain samples and no target enters."),
    ("Lean/Geometry/SourceTemporalObservation.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Linear acceptance grammar from a finite observation record; spans and kernels, not perturbation composition."),
    ("Lean/Geometry/SourceTemporalSelection.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Exact temporal grammar for a constrained information projection on scored full histories."),
    ("Lean/Geometry/SourceTemporalTomography.lean", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Native calibrated-path continuation from scalar seam means with no payload in the schedule; paths and incidence are explicit inputs and no metric menu is selected."),
    ("code/maxwell_measurement/README.md", "NEW_CANDIDATE", "maxwell_measurement_adapter",
     "Capture schema and preregistration for a Maxwell emulator; no capture or observed result is shipped."),
    ("code/maxwell_measurement/contract.py", "NEW_CANDIDATE", "maxwell_measurement_adapter",
     "Draft-and-freeze preregistration producer; a changed implementation invalidates a frozen contract."),
    ("code/sm_fermion_current/README.md", "NEW_CANDIDATE", "fermionic_hypercharge_current_histories",
     "Fermionic hypercharge current histories on the 64-site golden carrier with gauge and multiplet labels."),
    ("code/sm_fermion_current/coupled_README.md", "NEW_CANDIDATE", "fermionic_hypercharge_current_histories",
     "Coupled electric-link execution of the same labelled current construction."),
    ("code/sm_fermion_current/coupled_receipt.json", "NEW_CANDIDATE", "fermionic_hypercharge_current_histories",
     "Receipt of the coupled electric-link execution."),
    ("code/sm_fermion_current/current_receipt.json", "NEW_CANDIDATE", "fermionic_hypercharge_current_histories",
     "Authenticated register replay of baseline, phase intervention and gauge copy runs with hypercharge labels."),
    ("code/sm_fermion_current/quantum_link_receipt.json", "NEW_CANDIDATE", "fermionic_hypercharge_current_histories",
     "Receipt of the quantized electric link history."),
    ("code/source_checkpoint_selection/CONTRACT.md", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Contract of the checkpoint selection packet."),
    ("code/source_checkpoint_selection/README.md", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Native checkpoint selection on captured W12 seams with census commitments."),
    ("code/source_checkpoint_selection/controls.json", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Twelve controls and ten pipelines; canonical_axiom_selection = false."),
    ("code/source_checkpoint_selection/receipt.json", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Independent receipt of the checkpoint census and pipelines."),
    ("code/source_encoded_memory/CONTRACT.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Contract of the encoded memory packet."),
    ("code/source_encoded_memory/README.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Runtime controls of encoded copy, clear and reuse with attenuating amplitude."),
    ("code/source_encoded_memory/controls.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Captured copy and clear executions on the declared support."),
    ("code/source_encoded_memory/receipt.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Independent receipt; physical M1 and comparator are not derived."),
    ("code/source_feedback_transport/README.md", "EXTENDS_EXISTING_CANDIDATE", "source_feedback_transport",
     "Adds a cross-reference paragraph to the full-family compiler; the transport receipt and its classification are the same packet."),
    ("code/source_native_accumulator/CONTRACT.md", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Contract of the accumulator packet."),
    ("code/source_native_accumulator/README.md", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Reusable native accumulation with retirement at each start and a supplied request sequence."),
    ("code/source_native_accumulator/controls.json", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Exhaustive request words, three complete examples and reuse controls; source_selected = false."),
    ("code/source_native_accumulator/receipt.json", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Independent replay of 1,280 histories; dynamic scalar writes are pair means only."),
    ("code/source_native_programs/CONTRACT.md", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Contract of the stored-program packet."),
    ("code/source_native_programs/README.md", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Stored computation from pair means with two payload banks; conditional finite execution without a W12 embedding."),
    ("code/source_native_programs/controls.json", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Captured signed control with nine cases and eight routes; source_selected = false."),
    ("code/source_native_programs/families_receipt.json", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Family receipt separating complete small native executions from large-population route certificates."),
    ("code/source_native_programs/receipt.json", "NEW_CANDIDATE", "native_stored_programs_and_accumulator",
     "Independent replay receipt; dynamic writes are supported scalar pair means only and m1_derived = false."),
    ("code/source_native_updates/CONTRACT.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Contract of the native updates packet."),
    ("code/source_native_updates/README.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Native sum records with finite-chain cleanup on captured ports."),
    ("code/source_native_updates/controls.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Captured write and reread histories."),
    ("code/source_native_updates/receipt.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Independent receipt of histories, chain scaling and read bounds."),
    ("code/source_operator_join/README.md", "NEW_CANDIDATE", "finite_source_operator_join",
     "Deterministic operator-join packet of two committed pair algebras; a static noncommuting hinge with declared tensor assembly."),
    ("code/source_operator_join/operator_join_packet.json", "NEW_CANDIDATE", "finite_source_operator_join",
     "Mathematical evidence packet; overlapping_pair_algebras_commute = false and observed_quantum_outcomes = false."),
    ("code/source_passive_memory/CONTRACT.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Contract of the passive memory packet."),
    ("code/source_passive_memory/README.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Runtime controls of the passive-memory, reset and ancillary-budget obstructions."),
    ("code/source_passive_memory/controls.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Four-port equal-image witnesses and connectivity census of the level-three capture."),
    ("code/source_passive_memory/receipt.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Independent receipt of the executed means and fibre controls."),
    ("code/source_read_acceptance/CONTRACT.md", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Contract of the read acceptance packet."),
    ("code/source_read_acceptance/README.md", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Sound partial readout and a noisy two-payload read on captured seams."),
    ("code/source_read_acceptance/controls.json", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Exhaustive and example read histories with declared path laws."),
    ("code/source_read_acceptance/receipt.json", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Independent replay of stopped histories; source_contract = conditional_not_source_selected."),
    ("code/source_read_routing/CONTRACT.md", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Full-family routing contract; the exit is a conditional compiler theorem with M1 retained as a supplied rule."),
    ("code/source_read_routing/README.md", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Package front door for the q=13/q=21 routed read compiler and its regeneration commands."),
    ("code/source_read_routing/controls/q3_baseline.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Level-three baseline control receipt with a retained complete compressed tape."),
    ("code/source_read_routing/controls/q3_branch.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Level-three later-commit branch control receipt with a retained complete compressed tape."),
    ("code/source_read_routing/controls/q3_scratch.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Level-three initial-scratch control receipt with a retained complete compressed tape."),
    ("code/source_read_routing/controls/q3_source.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Level-three centre-source intervention control receipt with a retained complete compressed tape."),
    ("code/source_read_routing/specification.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Declares M1 as a retained supplied structural rule and records axiomatic_read_law_derived = false."),
    ("code/source_read_selection/CONTRACT.md", "EXTENDS_EXISTING_CANDIDATE", "lean_selected_pair_mean_histories",
     "Contract of the read selection packet."),
    ("code/source_read_selection/README.md", "EXTENDS_EXISTING_CANDIDATE", "lean_selected_pair_mean_histories",
     "Runtime controls of the constrained information-projection read obstruction on the captured level-three support."),
    ("code/source_read_selection/controls.json", "EXTENDS_EXISTING_CANDIDATE", "lean_selected_pair_mean_histories",
     "Captured and toy controls with twelve source ports and two receiver ports."),
    ("code/source_read_selection/receipt.json", "EXTENDS_EXISTING_CANDIDATE", "lean_selected_pair_mean_histories",
     "Independent receipt; status is a conditional obstruction, not a full-axiom countermodel."),
    ("code/source_reusable_bus/CONTRACT.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Contract of the reusable bus packet."),
    ("code/source_reusable_bus/README.md", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Repeated record reads between captured carriers 0 and 3 with complete fixed-point tapes."),
    ("code/source_reusable_bus/controls.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Sixteen initial payload combinations with complete tapes."),
    ("code/source_reusable_bus/receipt.json", "EXTENDS_EXISTING_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
     "Independent receipt of the histories, retained means and uniform read bounds."),
    ("code/source_routing_refinement/CONTRACT.md", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Contract of the routing refinement companion; M1 is not derived from the canonical source."),
    ("code/source_routing_refinement/README.md", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Routing refinement companion: support-distance bound across levels, 23-slot local store and raw-mass obstruction."),
    ("code/source_routing_refinement/receipt.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Independent receipt of the refinement hierarchy on the pinned level-three to level-five captures; M1_source_selected = false."),
    ("code/source_scalar_finite_instrument/README.md", "EXTENDS_EXISTING_CANDIDATE", "source_scalar_sequential_instrument",
     "Finite-duration readout certificate replaying the pinned sequential instrument; it samples no quantum outcomes."),
    ("code/source_scalar_finite_instrument/finite_instrument_receipt.json", "EXTENDS_EXISTING_CANDIDATE", "source_scalar_sequential_instrument",
     "Receipt of 21 readout slots, 210 earlier-readout pairs and 2,709 declared control operations."),
    ("code/source_temporal_acceptance/CONTRACT.md", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Contract of the temporal acceptance packet."),
    ("code/source_temporal_acceptance/README.md", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Temporal acceptance and complete native continuation on the captured 15,360-port support."),
    ("code/source_temporal_acceptance/controls.json", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Captured, toy, erasure, guarded and tomography controls; the large plan is checked structurally with payload_replay = false."),
    ("code/source_temporal_acceptance/receipt.json", "NEW_CANDIDATE", "native_temporal_tomography_and_checkpoint_selection",
     "Independent receipt of 27 captured histories and the native completion plan."),
    ("evidence/source_net_causal_poset/routed_read_law/q13_baseline.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Production run receipt binding inputs, producer source, event census and decoded-stream commitments for the q=13 baseline."),
    ("evidence/source_net_causal_poset/routed_read_law/q13_source.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Production run receipt for the q=13 centre payload +1 intervention."),
    ("evidence/source_net_causal_poset/routed_read_law/q21_baseline.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Production run receipt for the q=21 baseline."),
    ("evidence/source_net_causal_poset/routed_read_law/q21_source.json", "NEW_CANDIDATE", "source_read_routing_full_family",
     "Production run receipt for the q=21 centre payload +1 intervention."),
)


INTEGRATED_TREE_REVIEW_543298E0 = (('Lean/Geometry/SourceOperationReads.lean',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('Lean/Geometry/SourceOperationReadsAxiomAudit.lean',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('Lean/Geometry/SourcePublicationAxiomAudit.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourcePublicationLaw.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourcePublicationMass.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourcePublicationMenu.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourcePublicationNative.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourcePublicationPath.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourcePublicationPathInvariant.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourcePublicationTail.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourceRadiusRemoval.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourceRadiusSelection.lean',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('Lean/Geometry/SourceRepairInstrument.lean',
  'NEW_CANDIDATE',
  'declared_flagged_recovery_extensions',
  'Supplied swap branches, retained flags and constrained real-coordinate archive permit conditional global '
  'recovery; they are additional interfaces, not registered native operations.'),
 ('Lean/Geometry/SourceStateSelection.lean',
  'NEW_CANDIDATE',
  'finite_state_process_selection',
  'Finite compatible-state optimization and supplied channel controls; no source-selected current, temporal '
  'instrument, or native perturbation producer.'),
 ('code/source_operation_reads/DERIVATION.md',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/INSTRUMENT_DERIVATION.md',
  'NEW_CANDIDATE',
  'declared_flagged_recovery_extensions',
  'Supplied swap branches, retained flags and constrained real-coordinate archive permit conditional global '
  'recovery; they are additional interfaces, not registered native operations.'),
 ('code/source_operation_reads/README.md',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/__init__.py',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/capture.json',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/check_live.py',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/receipt.json',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/simulator_verifier.py',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/source_snapshot.json',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/test_reads.py',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_operation_reads/verify.py',
  'NEW_CANDIDATE',
  'registered_operation_read_histories',
  'Registered all-port histories with one prescribed native unitary and destructive matching means; full '
  'probe algebra does not supply reversible perturbation families.'),
 ('code/source_publication_selection/CONTRACT.md',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/DERIVATION.md',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/README.md',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/STATE_SELECTION.md',
  'NEW_CANDIDATE',
  'finite_state_process_selection',
  'Finite compatible-state optimization and supplied channel controls; no source-selected current, temporal '
  'instrument, or native perturbation producer.'),
 ('code/source_publication_selection/__init__.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/build.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/check_exterior.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/check_locality.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/check_quantum.py',
  'NEW_CANDIDATE',
  'finite_state_process_selection',
  'Finite compatible-state optimization and supplied channel controls; no source-selected current, temporal '
  'instrument, or native perturbation producer.'),
 ('code/source_publication_selection/codec.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/controls.json',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/exterior.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/locality.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/quantum.py',
  'NEW_CANDIDATE',
  'finite_state_process_selection',
  'Finite compatible-state optimization and supplied channel controls; no source-selected current, temporal '
  'instrument, or native perturbation producer.'),
 ('code/source_publication_selection/receipt.json',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/test_quantum.py',
  'NEW_CANDIDATE',
  'finite_state_process_selection',
  'Finite compatible-state optimization and supplied channel controls; no source-selected current, temporal '
  'instrument, or native perturbation producer.'),
 ('code/source_publication_selection/test_selection.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'),
 ('code/source_publication_selection/verify.py',
  'NEW_CANDIDATE',
  'native_eventual_publication',
  'Conditional publication law and scalar pair-mean histories on a declared connected grammar; no twelve '
  'reversible source perturbation families.'))

INTEGRATED_TREE_REVIEW_BA84976A = (
    ('Lean/Geometry/RecordGluingPrinciple.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Proves quadratic isotropy under finite carrier symmetry, graded path and refinement clock bounds, transcript capacity, reversible evaluation and interior-cube inequalities; it defines no port-indexed source perturbation, composition word, inverse or history.'),
    ('Lean/Geometry/RecordGluingPrincipleAxiomAudit.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Audits the axiom dependencies of the local-clock lemmas and adds no source operation or history.'),
    ('Lean/Geometry/SourceRecordGluing.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Formalizes the flight/wait grammar, cone, covering join, lowering and localized-intervention ancestry of classical retained records under the proposed RG law; classical register transport supplies no reversible twelve-port perturbation families or mixed-order response histories.'),
    ('Lean/Geometry/SourceRecordGluingAxiomAudit.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Audits the axiom dependencies of the record-gluing theorems and adds no source operation or history.'),
    ('Lean/Geometry/SourceSelectionLocality.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Proves the all-program local-transcript obstruction, the exact channel-decoder criterion and the population arithmetic of the scalar-seam model; it contains no port-indexed perturbation family or executed response history.'),
    ('Lean/Geometry/SourceSelectionLocalityAxiomAudit.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Audits the axiom dependencies of the locality theorems and adds no source operation or history.'),
    ('code/rg_principle/CONTRACT.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/DERIVATION.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/README.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/__init__.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/certificates.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/checker.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/compiler.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/receipt.json',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/receipt.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/rg_principle/test_rg_principle.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Local-clock reduction, finite Boolean compiler, copy-channel and resource controls for the proposed RG law; classical circuits and symbolic controls supply no port-indexed source perturbation histories.'),
    ('code/sm_fermion_current/test_spin_exchange_integration.py',
     'EXTENDS_EXISTING_CANDIDATE',
     'fermionic_hypercharge_current_histories',
     'Integration test binding the labelled fermionic current histories to the spin-exchange boundary; it carries the same downstream physical-current label and supplies no source-native perturbation family.'),
    ('code/source_record_gluing/README.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/__init__.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/build.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/capture.json',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/check_process.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/process.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/receipt.json',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/test_process.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_record_gluing/verify.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Reference execution of the proposed RG process grammar on classical retained records with an additive diagnostic; it carries no reversible port perturbation families or mixed-order response histories.'),
    ('code/source_selection_model/CONTRACT.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/DERIVATION.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/README.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/RECORD_GLUING.md',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/__init__.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/channels.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/finite_model.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/geometry.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/grammar.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/receipt.json',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/records.json',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/records.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/response.json',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/response.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/test_model.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/verify.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
    ('code/source_selection_model/verify_response.py',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Declared scalar-seam witness model: its twelve response generators are supplied by an explicit formula on an internal six-dimensional carrier and its ordered brackets are algebraic jets, not raw histories executed by a registered source producer on the twelve-port carrier; the non-entailment construction lies outside the source-native scope.'),
)

PRIOR_INTEGRATED_TREE_REVIEWS = (
    (
        "2b03a95caf5030272f7b426b964e816f650153f8",
        "2d9bd11bc47c56d88a2fbbca22e3cc1be171d3f9",
        (
        ("Lean/Geometry/FlatDiamondError.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Bounds geometric quadrature error; it defines no source perturbation, composition word, inverse, or history."),
        ("Lean/Geometry/FlatDiamondNormalization.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Normalizes causal-diamond pair volumes; it contains no port-indexed source operation or history custody."),
        ("Lean/Geometry/FlatDiamondPairIntegral.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Relates geometric pair sets to integrals; it supplies no perturbation family or ordered execution."),
        ("Lean/Geometry/FlatDiamondVolume.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Proves flat causal-diamond volume identities rather than source operations or response histories."),
        ("Lean/Geometry/FlatLorentzVolume.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Proves Lorentz-volume and boost identities with no source-current operation packet."),
        ("Lean/Geometry/GoldenSourceCausalLimit.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Establishes convergence of generated causal order and counts; it does not define reversible port perturbations."),
        ("Lean/Geometry/GoldenSourcePairLimit.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Establishes strict-pair count limits, not mixed-order source response histories."),
        ("Lean/Geometry/GoldenSourceVolumeLimit.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Establishes weighted source-volume limits with no perturbation composition or inverse law."),
        ("Lean/Geometry/OrderingFractionFourDimensional.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "The change connects geometric count limits to a scalar ordering fraction and adds no source operation."),
        ("Lean/Geometry/SourceBusScaling.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Defines conditional copy, clear, sweep, and cleanup words compiled to destructive scalar pair means."),
        ("Lean/Geometry/SourceCausalBoundary.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Proves null-boundary and fixed-diamond count facts; no source perturbation family is present."),
        ("Lean/Geometry/SourceConstrainedRead.lean", "NEW_CANDIDATE", "lean_selected_pair_mean_histories",
         "Connects constrained finite word laws to pair-mean transcripts, but only as an abstract irreversible read obstruction."),
        ("Lean/Geometry/SourceConstrainedSelection.lean", "NEW_CANDIDATE", "lean_selected_pair_mean_histories",
         "Supplies finite information-projection support theorems used by the selected-history packet."),
        ("Lean/Geometry/SourceCountLimitAxiomAudit.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Audits theorem dependencies for count limits and introduces no source operation or history."),
        ("Lean/Geometry/SourceCountTransport.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Transports selected finite counts to measure limits; predicates are classifiers rather than perturbation histories."),
        ("Lean/Geometry/SourceEncodedMemory.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Defines encoded copy, clear, hop, route, and read words whose amplitudes attenuate under pair means."),
        ("Lean/Geometry/SourceEncodedMemoryAxiomAudit.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Audits the encoded-memory theorem chain and adds no runtime histories or W12 binding."),
        ("Lean/Geometry/SourceNativeRecords.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Defines native arithmetic records and proves a two-order schedule-dependence example on Fin 3, not the W12 mixed-order suite."),
        ("Lean/Geometry/SourceNativeUpdatesAxiomAudit.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Audits native bus and record theorems without adding a source-selected runtime packet."),
        ("Lean/Geometry/SourceNetOrderLimit.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Proves limiting agreement of generated read order away from a null cone, not reversible perturbation histories."),
        ("Lean/Geometry/SourceNetVolumeError.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Bounds weighted Alexandrov volumes of generated orders and contains no current-tomography operation family."),
        ("Lean/Geometry/SourceNonlinearRecord.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Classifies state-only records protected by an irreversible pair-mean law."),
        ("Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Audits passive-memory, reset, and nonlinear-record theorems without supplying histories."),
        ("Lean/Geometry/SourcePassiveMemoryBudget.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Its quadratic defect proves nontrivial pair means dissipative and excludes nontrivial closed cycles."),
        ("Lean/Geometry/SourcePassiveReset.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Proves positive registers cannot reset exactly under finite nonnegative pair means."),
        ("Lean/Geometry/SourceReadSelection.lean", "NEW_CANDIDATE", "lean_selected_pair_mean_histories",
         "Defines finite pair-mean words and full prefix transcripts, but not W12 reversible operations or serialized executions."),
        ("Lean/Geometry/SourceReadSelectionAxiomAudit.lean", "NEW_CANDIDATE", "lean_selected_pair_mean_histories",
         "Audits the constrained-selection and transcript obstruction chain."),
        ("Lean/Geometry/SourceReusableBus.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Compiles copy and clear instructions to attenuating pair-mean programs with supplied schedules."),
        ("Lean/Geometry/SourceReusableBusAxiomAudit.lean", "NEW_CANDIDATE", "lean_pair_mean_memory_and_reusable_bus",
         "Audits the reusable-bus theorem chain without adding a W12 source instance."),
        ("Lean/Geometry/SourceSelectionControls.lean", "NEW_CANDIDATE", "lean_selected_pair_mean_histories",
         "Provides finite controls showing cover-dependent support selection, not a reversible current producer."),
        ("Lean/Screen/OPHScreen.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "The change only imports Whitney modules into the umbrella and adds no theorem or source operation."),
        ("Lean/Screen/WhitneyConeMass.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Defines cone mass, constraint, stiffness, and normal frames; these are static field-mode structures."),
        ("Lean/Screen/WhitneyConeMassNaturality.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Proves frame naturality for a supplied mass and contains no source perturbation histories."),
        ("Lean/Screen/WhitneyConeModes.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Constructs static cone gradients, curls, constraints, and thirty normal modes."),
        ("Lean/Screen/WhitneyConeNaturality.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Proves basis-change naturality of cone modes rather than source-operation order sensitivity."),
        ("Lean/Screen/WhitneyFrameMorphism.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Defines isometries between modal frames, not port-indexed perturbations or runtime histories."),
        ("Lean/Screen/WhitneyFrameNaturality.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Defines overlap-coordinate changes among complete modal frames; the order is basis order, not operation order."),
        ("Lean/Screen/WhitneyMassPremiseInstance.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Instantiates a positive counting mass and thirty-mode frame with no source-current producer."),
        ("Lean/Screen/WhitneyNormalModeConstruction.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Applies finite-dimensional spectral theory to supplied mass and stiffness forms."),
        ("Lean/Screen/WhitneyRotationWitness.lean", "OUTSIDE_REVIEWED_SOURCE_SCOPE", None,
         "Provides a two-mode rotation control for frame naturality, not source-current histories."),
        ),
    ),
    (
        "2d9bd11bc47c56d88a2fbbca22e3cc1be171d3f9",
        "afed734528edff214c34d4038a64310544df922d",
        PRIOR_BASE_INTEGRATED_TREE_REVIEW,
    ),
    (
        "afed734528edff214c34d4038a64310544df922d",
        "543298e06ce41d47f44c7707fdb2c32723cca10e",
        INTEGRATED_TREE_REVIEW_543298E0,
    ),
    (
        "543298e06ce41d47f44c7707fdb2c32723cca10e",
        "ba84976ad1195892c68fdd0d74a27ab0c899aeb8",
        INTEGRATED_TREE_REVIEW_BA84976A,
    ),
)

INTEGRATED_TREE_REVIEW = (
    ('Lean/Geometry/M1Interfaces.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Bounds finite crossing changes, binary volume correction and thin connecting paths and checks the sparse interface scale algebra; it contains no port-indexed source perturbation family or executed response history.'),
    ('Lean/Geometry/M1InterfacesAxiomAudit.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Audits the axiom dependencies of the interface-scale reductions and adds no source operation or history.'),
    ('Lean/Geometry/M1Necessity.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Checks the binary route subroutine, layer obstruction, finite cut algebra and positive-action normalization bounds of the sparse read-family comparison; it defines no port-indexed source perturbation, composition word, inverse or history.'),
    ('Lean/Geometry/M1NecessityAxiomAudit.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Audits the axiom dependencies of the sparse read-family reductions and adds no source operation or history.'),
    ('Lean/Geometry/M1QuantumTransport.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Proves the dot-product deficit, weighted trace budget and the other finite reductions of the native-time quantum transport theorem; it defines no port-indexed source operation or raw history.'),
    ('Lean/Geometry/M1QuantumTransportAxiomAudit.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Audits the axiom dependencies of the quantum-transport reductions and adds no source operation or history.'),
    ('Lean/Geometry/M1VacuumFidelity.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Proves finite symplectic-mode algebra identities for the full-spectrum vacuum fidelity theorem; scalar mode arithmetic supplies no twelve-port perturbation family or mixed-order history.'),
    ('Lean/Geometry/M1VacuumFidelityAxiomAudit.lean',
     'OUTSIDE_REVIEWED_SOURCE_SCOPE',
     None,
     'Audits the axiom dependencies of the vacuum-fidelity reductions and adds no source operation or history.'),
)


ALLOWED_CLASSIFICATIONS = {
    "QUALIFIES",
    "PARTIAL",
    "STATIC_ONLY",
    "COMMUTATIVE_ONLY",
    "IRREVERSIBLE_ONLY",
    "DOWNSTREAM_CONTAMINATED",
    "MISSING_ORDER_INFORMATION",
    "UNRELATED",
}

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

AUDITED_SOURCE_PATHS = (
    "code/a5_closure/manifests/echosahedral_federation_reference.json",
    "code/a5_closure/manifests/charged_response_semantic_artifact.json",
    "code/a5_closure/manifests/source_current_capability_projection.json",
    "code/a5_closure/receipts/source_current_capability.receipt.json",
    "code/source_routing/README.md",
    "code/source_routing/runtime/path_tomography_receipt.json",
    "code/source_feedback_transport/transport_receipt.json",
    "code/source_scalar_instruments/sequential_instrument_receipt.json",
    "code/a5_closure/manifests/record_counting_mechanism_reference.json",
    "code/a5_closure/manifests/directed_seam_repair_reference.json",
    "code/consensus/manifests/compiled_lattice_settling_reference.json",
    "code/refinement/runtime/discrete_refinement_theorem_receipt.json",
    "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json",
    "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_stage2.receipt.json",
    "code/b14_jacobi/oriented_face_bracket_selector.certificate.json",
    "code/a5_closure/manifests/response_grammar_completeness_reference.json",
    "code/a5_closure/receipts/response_grammar_completeness_reference.receipt.json",
    "code/a5_closure/manifests/charged_response_pole_residue_artifact.json",
    "code/a5_closure/manifests/spin_statistics_semantic_artifact.json",
    "code/a5_closure/manifests/port_current_response_reference.json",
    "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
    "code/angular_sprint/runtime/kinetic_form_selection_receipt.json",
    "code/consensus/transport_from_gluing.py",
    "Lean/InformationProjection/SourceHistoryPacket.lean",
    "Lean/Variational/SourceToHamiltonianComposed.lean",
    "Lean/QFT/SourceOperatorGeneration.lean",
    "Lean/QFT/SourceCorrelationCapstone.lean",
    "Lean/QFT/SourceHistoryGNSDynamics.lean",
    "Lean/QFT/SourceDerivedEventPrecedence.lean",
    "Lean/EventAlgebra/SourceBoundInstrumentInterface.lean",
    "Lean/ObserverPatchHolography/A2EndpointCommutator.lean",
    "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean",
    "Lean/ObserverPatchHolography/Provenance/RefinementNaturality.lean",
    "Lean/Geometry/SourceFeedbackTransport.lean",
    "Lean/Geometry/SourceSeamPathTomography.lean",
    "Lean/Geometry/CarrierDynamicsCompatibility.lean",
    "Lean/Geometry/WorldlineHopTransport.lean",
    "Lean/Screen/A2HolonomyBridge.lean",
    "Lean/Screen/A5ResponseWordAlgebra.lean",
    "Lean/Screen/CarrierEvolutionFlow.lean",
    "Lean/Screen/OrientedFaceBracketSelector.lean",
    "Lean/Screen/PrimitivePortTranslationBridge.lean",
    "Lean/Screen/RepairWordCarrierReadout.lean",
    "Lean/Geometry/SourceBusScaling.lean",
    "Lean/Geometry/SourceConstrainedRead.lean",
    "Lean/Geometry/SourceConstrainedSelection.lean",
    "Lean/Geometry/SourceEncodedMemory.lean",
    "Lean/Geometry/SourceEncodedMemoryAxiomAudit.lean",
    "Lean/Geometry/SourceNativeRecords.lean",
    "Lean/Geometry/SourceNativeUpdatesAxiomAudit.lean",
    "Lean/Geometry/SourceNonlinearRecord.lean",
    "Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean",
    "Lean/Geometry/SourcePassiveMemoryBudget.lean",
    "Lean/Geometry/SourcePassiveReset.lean",
    "Lean/Geometry/SourceReadSelection.lean",
    "Lean/Geometry/SourceReadSelectionAxiomAudit.lean",
    "Lean/Geometry/SourceReusableBus.lean",
    "Lean/Geometry/SourceReusableBusAxiomAudit.lean",
    "Lean/Geometry/SourceSelectionControls.lean",
    "Lean/Geometry/SourceAccumulatorAxiomAudit.lean",
    "Lean/Geometry/SourceAccumulatorDecoder.lean",
    "Lean/Geometry/SourceAccumulatorProgram.lean",
    "Lean/Geometry/SourceBankBusWitness.lean",
    "Lean/Geometry/SourceBankCompiler.lean",
    "Lean/Geometry/SourceBankControl.lean",
    "Lean/Geometry/SourceBankExecution.lean",
    "Lean/Geometry/SourceBankInvariant.lean",
    "Lean/Geometry/SourceBankLowering.lean",
    "Lean/Geometry/SourceBankMachine.lean",
    "Lean/Geometry/SourceBankPrecision.lean",
    "Lean/Geometry/SourceBankTrace.lean",
    "Lean/Geometry/SourceNativeAccumulator.lean",
    "Lean/Geometry/SourceNativeCore.lean",
    "Lean/Geometry/SourceNativeProgramBudget.lean",
    "Lean/Geometry/SourceNativeProgramError.lean",
    "Lean/Geometry/SourceNativeProgramsAxiomAudit.lean",
    "Lean/Geometry/SourceNativeShuttle.lean",
    "Lean/Geometry/SourceNativeStoredProgram.lean",
    "Lean/Geometry/SourceReadRouting.lean",
    "Lean/Geometry/SourceRoutingBudget.lean",
    "Lean/Geometry/SourceRoutingHierarchy.lean",
    "Lean/Geometry/SourceRoutingRefinementAxiomAudit.lean",
    "Lean/Geometry/SourceRoutingStorage.lean",
    "Lean/Geometry/SourceCheckpointAxiomAudit.lean",
    "Lean/Geometry/SourceCheckpointEntropy.lean",
    "Lean/Geometry/SourceCheckpointMeaning.lean",
    "Lean/Geometry/SourceCheckpointNative.lean",
    "Lean/Geometry/SourceCheckpointPipeline.lean",
    "Lean/Geometry/SourceCheckpointPolicy.lean",
    "Lean/Geometry/SourceReadAcceptance.lean",
    "Lean/Geometry/SourceReadAcceptanceAxiomAudit.lean",
    "Lean/Geometry/SourceReadAcceptanceSchedule.lean",
    "Lean/Geometry/SourceTemporalAcceptanceAxiomAudit.lean",
    "Lean/Geometry/SourceTemporalAttempts.lean",
    "Lean/Geometry/SourceTemporalConnected.lean",
    "Lean/Geometry/SourceTemporalErasure.lean",
    "Lean/Geometry/SourceTemporalEssential.lean",
    "Lean/Geometry/SourceTemporalGuard.lean",
    "Lean/Geometry/SourceTemporalMenu.lean",
    "Lean/Geometry/SourceTemporalNative.lean",
    "Lean/Geometry/SourceTemporalObservation.lean",
    "Lean/Geometry/SourceTemporalSelection.lean",
    "Lean/Geometry/SourceTemporalTomography.lean",
    "code/source_read_routing/CONTRACT.md",
    "code/source_read_routing/README.md",
    "code/source_read_routing/specification.json",
    "code/source_read_routing/controls/q3_baseline.json",
    "code/source_read_routing/controls/q3_source.json",
    "code/source_read_routing/controls/q3_branch.json",
    "code/source_read_routing/controls/q3_scratch.json",
    "evidence/source_net_causal_poset/routed_read_law/q13_baseline.json",
    "evidence/source_net_causal_poset/routed_read_law/q13_source.json",
    "evidence/source_net_causal_poset/routed_read_law/q21_baseline.json",
    "evidence/source_net_causal_poset/routed_read_law/q21_source.json",
    "code/source_routing_refinement/README.md",
    "code/source_routing_refinement/CONTRACT.md",
    "code/source_routing_refinement/receipt.json",
    "code/source_native_programs/README.md",
    "code/source_native_programs/CONTRACT.md",
    "code/source_native_programs/controls.json",
    "code/source_native_programs/receipt.json",
    "code/source_native_programs/families_receipt.json",
    "code/source_native_accumulator/README.md",
    "code/source_native_accumulator/CONTRACT.md",
    "code/source_native_accumulator/controls.json",
    "code/source_native_accumulator/receipt.json",
    "code/source_temporal_acceptance/README.md",
    "code/source_temporal_acceptance/CONTRACT.md",
    "code/source_temporal_acceptance/controls.json",
    "code/source_temporal_acceptance/receipt.json",
    "code/source_read_acceptance/README.md",
    "code/source_read_acceptance/CONTRACT.md",
    "code/source_read_acceptance/controls.json",
    "code/source_read_acceptance/receipt.json",
    "code/source_checkpoint_selection/README.md",
    "code/source_checkpoint_selection/CONTRACT.md",
    "code/source_checkpoint_selection/controls.json",
    "code/source_checkpoint_selection/receipt.json",
    "code/source_read_selection/README.md",
    "code/source_read_selection/CONTRACT.md",
    "code/source_read_selection/controls.json",
    "code/source_read_selection/receipt.json",
    "code/source_encoded_memory/README.md",
    "code/source_encoded_memory/CONTRACT.md",
    "code/source_encoded_memory/controls.json",
    "code/source_encoded_memory/receipt.json",
    "code/source_passive_memory/README.md",
    "code/source_passive_memory/CONTRACT.md",
    "code/source_passive_memory/controls.json",
    "code/source_passive_memory/receipt.json",
    "code/source_reusable_bus/README.md",
    "code/source_reusable_bus/CONTRACT.md",
    "code/source_reusable_bus/controls.json",
    "code/source_reusable_bus/receipt.json",
    "code/source_native_updates/README.md",
    "code/source_native_updates/CONTRACT.md",
    "code/source_native_updates/controls.json",
    "code/source_native_updates/receipt.json",
    "code/source_scalar_finite_instrument/README.md",
    "code/source_scalar_finite_instrument/finite_instrument_receipt.json",
    "code/source_operator_join/README.md",
    "code/source_operator_join/operator_join_packet.json",
    "code/sm_fermion_current/README.md",
    "code/sm_fermion_current/coupled_README.md",
    "code/sm_fermion_current/current_receipt.json",
    "code/sm_fermion_current/coupled_receipt.json",
    "code/sm_fermion_current/quantum_link_receipt.json",
    "code/maxwell_measurement/README.md",
    "code/maxwell_measurement/contract.py",
    "code/source_feedback_transport/README.md",
    "Lean/QFT/TripleCarrierOperatorJoin.lean",
    'Lean/Geometry/SourceOperationReads.lean',
    'Lean/Geometry/SourceOperationReadsAxiomAudit.lean',
    'Lean/Geometry/SourcePublicationAxiomAudit.lean',
    'Lean/Geometry/SourcePublicationLaw.lean',
    'Lean/Geometry/SourcePublicationMass.lean',
    'Lean/Geometry/SourcePublicationMenu.lean',
    'Lean/Geometry/SourcePublicationNative.lean',
    'Lean/Geometry/SourcePublicationPath.lean',
    'Lean/Geometry/SourcePublicationPathInvariant.lean',
    'Lean/Geometry/SourcePublicationTail.lean',
    'Lean/Geometry/SourceRadiusRemoval.lean',
    'Lean/Geometry/SourceRadiusSelection.lean',
    'Lean/Geometry/SourceRepairInstrument.lean',
    'Lean/Geometry/SourceStateSelection.lean',
    'code/source_operation_reads/DERIVATION.md',
    'code/source_operation_reads/INSTRUMENT_DERIVATION.md',
    'code/source_operation_reads/README.md',
    'code/source_operation_reads/__init__.py',
    'code/source_operation_reads/capture.json',
    'code/source_operation_reads/check_live.py',
    'code/source_operation_reads/receipt.json',
    'code/source_operation_reads/simulator_verifier.py',
    'code/source_operation_reads/source_snapshot.json',
    'code/source_operation_reads/test_reads.py',
    'code/source_operation_reads/verify.py',
    'code/source_publication_selection/CONTRACT.md',
    'code/source_publication_selection/DERIVATION.md',
    'code/source_publication_selection/README.md',
    'code/source_publication_selection/STATE_SELECTION.md',
    'code/source_publication_selection/__init__.py',
    'code/source_publication_selection/build.py',
    'code/source_publication_selection/check_exterior.py',
    'code/source_publication_selection/check_locality.py',
    'code/source_publication_selection/check_quantum.py',
    'code/source_publication_selection/codec.py',
    'code/source_publication_selection/controls.json',
    'code/source_publication_selection/exterior.py',
    'code/source_publication_selection/locality.py',
    'code/source_publication_selection/quantum.py',
    'code/source_publication_selection/receipt.json',
    'code/source_publication_selection/test_quantum.py',
    'code/source_publication_selection/test_selection.py',
    'code/source_publication_selection/verify.py',
)


class InventoryError(ValueError):
    """Typed fail-closed inventory error."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise InventoryError(code, message)


def load_json(relative: str) -> dict[str, Any]:
    path = REPO_ROOT / relative
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), "JSON_OBJECT", f"{relative} is not an object")
    return value


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def canonical_sha256(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def source_pin(relative: str) -> dict[str, Any]:
    path = REPO_ROOT / relative
    require(path.is_file(), "SOURCE_PATH", f"missing audited source: {relative}")
    return {
        "path": relative,
        "bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def snapshot_exclusion_reason(relative: str) -> str | None:
    path = Path(relative)
    if relative in CONTENT_SNAPSHOT_EXACT_EXCLUSIONS:
        return "self_referential_inventory_manifest"
    if any(part in CONTENT_SNAPSHOT_EXCLUDED_DIRECTORY_NAMES for part in path.parts):
        return "repository_cache_directory"
    if path.suffix in CONTENT_SNAPSHOT_EXCLUDED_SUFFIXES:
        return "repository_cache_or_build_output"
    return None


def audited_file_records(repo_root: Path = REPO_ROOT) -> list[dict[str, Any]]:
    """Return canonical content records for every included audited-tree file."""

    paths: set[str] = set()
    for relative in AUDITED_DIRECTORIES:
        directory = repo_root / relative
        require(
            directory.is_dir(),
            "AUDIT_DIRECTORY",
            f"missing audited directory: {relative}",
        )
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


def audited_file_paths(repo_root: Path = REPO_ROOT) -> list[str]:
    return [row["path"] for row in audited_file_records(repo_root)]


def audited_file_snapshot() -> dict[str, Any]:
    files = audited_file_records()
    return {
        "directories": list(AUDITED_DIRECTORIES),
        "file_count": len(files),
        "files": files,
        "files_sha256": canonical_sha256(files),
        "exclusion_policy": {
            "exact_paths": list(CONTENT_SNAPSHOT_EXACT_EXCLUSIONS),
            "directory_names": list(CONTENT_SNAPSHOT_EXCLUDED_DIRECTORY_NAMES),
            "suffixes": list(CONTENT_SNAPSHOT_EXCLUDED_SUFFIXES),
            "reasons": {
                "exact_paths": "the generated inventory would otherwise contain its own byte digest",
                "directory_names": "repository cache directories cannot supply reviewed scientific source",
                "suffixes": "interpreter cache and build outputs cannot supply reviewed scientific source",
            },
        },
        "staleness_rule": (
            "any added, removed, renamed, or byte-changed included file below "
            "an audited directory invalidates the committed inventory"
        ),
    }


def properties(**values: bool) -> dict[str, bool]:
    require(
        set(values) == set(QUALIFICATION_FIELDS),
        "PROPERTY_SCHEMA",
        "candidate property set drift",
    )
    require(
        all(isinstance(value, bool) for value in values.values()),
        "PROPERTY_TYPE",
        "candidate properties must be booleans",
    )
    return {key: values[key] for key in QUALIFICATION_FIELDS}


def candidate(
    candidate_id: str,
    classification: str,
    source_paths: Sequence[str],
    source_ancestry: str,
    history_custody: str,
    observed_properties: Mapping[str, bool],
    rejection_reason: str,
    evidence_facts: Sequence[str],
) -> dict[str, Any]:
    require(
        classification in ALLOWED_CLASSIFICATIONS,
        "CLASSIFICATION",
        f"invalid classification for {candidate_id}",
    )
    props = dict(observed_properties)
    qualifies = all(props[field] for field in QUALIFICATION_FIELDS)
    require(
        (classification == "QUALIFIES") == qualifies,
        "QUALIFICATION",
        f"classification/property mismatch for {candidate_id}",
    )
    return {
        "candidate_id": candidate_id,
        "classification": classification,
        "source_paths": list(source_paths),
        "source_ancestry": source_ancestry,
        "history_custody": history_custody,
        "properties": props,
        "qualification_failures": [
            field for field in QUALIFICATION_FIELDS if not props[field]
        ],
        "rejection_reason": rejection_reason,
        "evidence_facts": list(evidence_facts),
    }


def audit_source_semantics() -> dict[str, Any]:
    """Check the decisive source facts before any classification is emitted."""

    carrier = load_json("code/a5_closure/manifests/echosahedral_federation_reference.json")
    ports = carrier.get("carrier", {}).get("ports", [])
    require(len(ports) == 12 and len(set(ports)) == 12, "CARRIER", "twelve ports required")

    projection = load_json(
        "code/a5_closure/manifests/source_current_capability_projection.json"
    )
    capability = load_json(
        "code/a5_closure/receipts/source_current_capability.receipt.json"
    )
    recurrence = projection.get("registered_recurrence", {})
    algebra = capability.get("audit", {}).get("response_word_algebra", {})
    rechart = capability.get("audit", {}).get("recharting_audit", {})
    require(
        projection.get("schema") == "oph.source_current_capability_projection.v1"
        and recurrence.get("raw_runtime_history_serialized") is False
        and recurrence.get("source_count") == 12
        and recurrence.get("readback_count") == 12,
        "CAPABILITY_SCOPE",
        "registered recurrence scope drift",
    )
    require(
        algebra.get("exact_dimension") == 4
        and algebra.get("commutative") is True
        and algebra.get("commutator_nonzero_count") == 0
        and rechart.get("proper_recharting_count") == 60
        and rechart.get("proper_rechartings_in_response_word_algebra") == 1
        and capability.get("physical_current_source_bridge_attained") is False,
        "CAPABILITY_OBSTRUCTION",
        "bounded response-algebra obstruction drift",
    )

    path_tomography = load_json(
        "code/source_routing/runtime/path_tomography_receipt.json"
    )
    path_scope = path_tomography.get("scope", {})
    path_ops = {
        event.get("op")
        for episode in path_tomography.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    routing_readme = (REPO_ROOT / "code/source_routing/README.md").read_text(
        encoding="utf-8"
    )
    require(
        path_tomography.get("schema") == "oph.source_routing.path_tomography.v1"
        and path_scope.get("all_intermediate_events_retained") is True
        and path_scope.get("canonical_scalar_pair_means_only") is True
        and path_scope.get("full_metric_neighbor_refinement") is False
        and "pair_mean" in path_ops
        and "destructive readout construction" in routing_readme,
        "PATH_TOMOGRAPHY",
        "path tomography typing drift",
    )

    feedback = load_json("code/source_feedback_transport/transport_receipt.json")
    feedback_ops = {
        event.get("op")
        for episode in feedback.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    require(
        feedback.get("schema") == "oph.source_feedback_transport.receipt.v1"
        and "SUPPLIED_COPY_RESET_LAW" in str(feedback.get("scope"))
        and {"mean", "reset", "export"} <= feedback_ops,
        "FEEDBACK_TRANSPORT",
        "feedback transport operation typing drift",
    )

    scalar = load_json(
        "code/source_scalar_instruments/sequential_instrument_receipt.json"
    )
    scalar_scope = scalar.get("scope", {})
    require(
        scalar.get("schema") == "oph.source_scalar.centered_sequential_instrument.v1"
        and scalar_scope.get("sampled_quantum_outcome_histories") is False
        and scalar_scope.get("source_selected_quantum_operations") is False
        and scalar_scope.get("physical_clock_or_observed_outcomes") is False
        and len(scalar.get("regional_gadget", {}).get("detector_sites", [])) == 32,
        "SCALAR_INSTRUMENT",
        "sequential instrument scope drift",
    )

    counting = load_json(
        "code/a5_closure/manifests/record_counting_mechanism_reference.json"
    )
    counting_mechanism = counting.get("mechanism", {})
    presentation = counting.get("a2_clauses", {}).get("presentation_invariance", {})
    require(
        counting.get("schema") == "oph.record_counting_mechanism_certificate.v2"
        and "append-only" in str(counting_mechanism.get("event_grammar"))
        and "transfer" in str(counting_mechanism.get("repair_move"))
        and "addition is commutative" in str(presentation.get("theorem"))
        and presentation.get("permutations_checked") == 99,
        "RECORD_COUNTING",
        "record-counting irreversibility/order-erasure drift",
    )

    directed = load_json(
        "code/a5_closure/manifests/directed_seam_repair_reference.json"
    )
    directed_verdict = directed.get("verdict", {})
    require(
        directed.get("schema") == "oph.directed_seam_repair_certificate.v1"
        and directed_verdict.get("protected_record_and_checkpoint_instrument") == "open"
        and directed_verdict.get("refinement_semigroup_naturality") == "open"
        and directed.get("integer_macro_bridge", {}).get("expectation_equals_pair_average")
        is True,
        "DIRECTED_REPAIR",
        "directed repair boundary drift",
    )

    consensus = load_json(
        "code/consensus/manifests/compiled_lattice_settling_reference.json"
    )
    consensus_scope = consensus.get("scope", {})
    require(
        consensus.get("schema") == "oph.compiled_lattice_settling_certificate.v2"
        and "finite software patch dynamics" in str(consensus_scope.get("realized"))
        and "physical-hardware attachment" in str(consensus_scope.get("open")),
        "CONSENSUS_SCOPE",
        "consensus claim boundary drift",
    )

    refinement = load_json(
        "code/refinement/runtime/discrete_refinement_theorem_receipt.json"
    )
    require(
        refinement.get("schema") == "oph.discrete_refinement_theorem_packet.v1"
        and refinement.get("certified_exit")
        == "THEOREM_PACKET_AND_MESH_SCAFFOLD_ATTAINED__SOURCE_NATIVE_PHYSICAL_BIREFINEMENT_OPEN",
        "REFINEMENT_SCOPE",
        "refinement scaffold boundary drift",
    )

    bracket1 = load_json(
        "code/a5_closure/issue_566_bracket_space_stage1/"
        "a5_alternating_bracket_space_stage1.receipt.json"
    )
    bracket2 = load_json(
        "code/a5_closure/issue_566_bracket_space_stage2/"
        "a5_jacobi_stage2.receipt.json"
    )
    require(
        bracket1.get("status") == "EXACT_TARGET_FREE_SEARCH_SPACE_STAGE1_ONLY"
        and not any(bracket1.get("later_gates", {}).values())
        and bracket2.get("status")
        == "EXACT_JACOBI_SYSTEM_AND_FIXED_LINE_REDUCTION__FULL_CLASSIFICATION_OPEN"
        and bracket2.get("not_attained", {}).get("preferred_bracket_selected") is False,
        "BRACKET_SPACE",
        "issue-566 search-space boundary drift",
    )

    b14 = load_json(
        "code/b14_jacobi/oriented_face_bracket_selector.certificate.json"
    )
    require(
        b14.get("schema") == "oph.b14.oriented_face_bracket_selector.v1"
        and b14.get("jacobi_failure", {}).get("nonzero_count") == 240
        and b14.get("jacobi_failure", {}).get("squared_norm") == 240,
        "B14_CONTROL",
        "B14 Jacobi control drift",
    )

    grammar = load_json(
        "code/a5_closure/manifests/response_grammar_completeness_reference.json"
    )
    grammar_receipt = load_json(
        "code/a5_closure/receipts/response_grammar_completeness_reference.receipt.json"
    )
    require(
        grammar.get("schema") == "oph.response_grammar_completeness_certificate.v1"
        and grammar_receipt.get("schema")
        == "oph.response_grammar_completeness_receipt.v1"
        and grammar_receipt.get("verdict", {}).get("bounded_exit")
        == "independence_limited"
        and "continuum dynamics"
        in grammar_receipt.get("claim_boundary", {}).get("does_not_close", []),
        "RESPONSE_GRAMMAR",
        "response grammar boundary drift",
    )

    spin = load_json(
        "code/a5_closure/manifests/spin_statistics_semantic_artifact.json"
    )
    require(
        spin.get("schema") == "oph.spin_statistics_semantic_artifact.v1"
        and spin.get("physical_source_gate", {}).get("passed") is True
        and spin.get("physical_source_gate", {}).get(
            "continuum_spin_statistics_theorem"
        )
        is False
        and spin.get("physical_source_gate", {}).get(
            "laboratory_exchange_measurement"
        )
        is False,
        "SPIN_TRANSPORT",
        "spin/deck source gate drift",
    )

    fixture = load_json(
        "code/a5_closure/manifests/port_current_response_reference.json"
    )
    require(
        fixture.get("schema") == "oph.port_current_response_manifest.v5"
        and fixture.get("construction_model") == "charged_double_triplet"
        and "fixture" in str(
            fixture.get("response_declaration_contract", {}).get("status")
        ),
        "FIXTURE_BOUNDARY",
        "conditional port-current fixture typing drift",
    )

    kinetic = load_json(
        "code/angular_sprint/runtime/kinetic_form_selection_receipt.json"
    )
    require(
        kinetic.get("schema") == "oph.kinetic_form_selection_receipt.v1"
        and kinetic.get("physical_selection_boundary", {}).get(
            "current_lift_source_selected"
        )
        is False
        and "charged-double-triplet" in str(kinetic.get("current_fixture")),
        "KINETIC_FIXTURE",
        "kinetic-form fixture boundary drift",
    )

    routing_spec = load_json("code/source_read_routing/specification.json")
    routing_scope = routing_spec.get("scope", {})
    require(
        routing_spec.get("schema") == "oph.source_read_routing.specification.v1"
        and routing_scope.get("full_q13_q21_routed_execution") is True
        and routing_scope.get("axiomatic_read_law_derived") is False
        and routing_scope.get("universal_routing_impossibility") is False
        and routing_spec.get("M1", {}).get("status") == "retained supplied structural rule"
        and routing_spec.get("production_variants") == ["baseline", "source"],
        "ROUTING_SCOPE",
        "full-family routing specification drift",
    )
    routing_runs = {
        name: load_json(f"evidence/source_net_causal_poset/routed_read_law/{name}.json")
        for name in ("q13_baseline", "q13_source", "q21_baseline", "q21_source")
    }
    require(
        all(run.get("schema") == "oph.source_read_routing.run.v1" for run in routing_runs.values())
        and routing_runs["q13_baseline"].get("variant") == "baseline"
        and routing_runs["q13_source"].get("variant") == "source"
        and routing_runs["q21_baseline"].get("variant") == "baseline"
        and routing_runs["q21_source"].get("variant") == "source"
        and routing_runs["q13_baseline"].get("costs", {}).get("events")
        == routing_runs["q13_source"].get("costs", {}).get("events")
        and routing_runs["q21_baseline"].get("costs", {}).get("events")
        == routing_runs["q21_source"].get("costs", {}).get("events"),
        "ROUTING_RUNS",
        "full-family routing run receipt drift",
    )
    refinement = load_json("code/source_routing_refinement/receipt.json")
    require(
        refinement.get("scope", {}).get("M1_source_selected") is False
        and refinement.get("hierarchy", {}).get("base_diameter") == 3,
        "ROUTING_REFINEMENT",
        "routing refinement receipt drift",
    )
    programs = load_json("code/source_native_programs/receipt.json")
    accumulator = load_json("code/source_native_accumulator/receipt.json")
    require(
        programs.get("m1_derived") is False
        and programs.get("dynamic_writes") == "supported_scalar_pair_means_only"
        and accumulator.get("m1_derived") is False
        and accumulator.get("dynamic_scalar_writes") == "pair_means_only",
        "NATIVE_PROGRAMS",
        "native stored-program or accumulator receipt drift",
    )
    temporal_controls = load_json("code/source_temporal_acceptance/controls.json")
    temporal = load_json("code/source_temporal_acceptance/receipt.json")
    captured_plan = temporal_controls.get("tomography", {}).get("captured_plan", {})
    require(
        temporal_controls.get("full_axiom_instantiation") is False
        and temporal_controls.get("m1_derived") is False
        and captured_plan.get("payload_replay") is False
        and captured_plan.get("all_ports_covered") is True
        and temporal.get("m1_derived") is False,
        "TEMPORAL_ACCEPTANCE",
        "temporal acceptance control or receipt drift",
    )
    read_acceptance = load_json("code/source_read_acceptance/receipt.json")
    checkpoint = load_json("code/source_checkpoint_selection/controls.json")
    require(
        read_acceptance.get("source_contract") == "conditional_not_source_selected"
        and read_acceptance.get("m1_derived") is False
        and checkpoint.get("scope", {}).get("canonical_axiom_selection") is False
        and checkpoint.get("scope", {}).get("metric_radius_selection") is False,
        "ACCEPTANCE_SELECTION",
        "read acceptance or checkpoint selection scope drift",
    )
    read_selection_controls = load_json("code/source_read_selection/controls.json")
    read_selection = load_json("code/source_read_selection/receipt.json")
    require(
        read_selection_controls.get("source_ports") == list(range(12))
        and read_selection.get("status") == "conditional_obstruction_not_full_axiom_countermodel",
        "READ_SELECTION",
        "read selection control or receipt drift",
    )
    reusable_bus = load_json("code/source_reusable_bus/receipt.json")
    native_updates = load_json("code/source_native_updates/receipt.json")
    require(
        isinstance(reusable_bus.get("histories"), list)
        and len(reusable_bus["histories"]) == 16
        and isinstance(native_updates.get("histories"), list)
        and len(native_updates["histories"]) >= 1,
        "PAIR_MEAN_MEMORY_RUNTIME",
        "reusable-bus or native-update history drift",
    )
    finite_instrument = load_json(
        "code/source_scalar_finite_instrument/finite_instrument_receipt.json"
    )
    require(
        finite_instrument.get("scope", {}).get("exact_rectangular_field_pointer_pulse") is True
        and "code/source_scalar_instruments/sequential_instrument_receipt.json"
        in str(finite_instrument.get("parents")),
        "FINITE_INSTRUMENT",
        "finite instrument receipt drift",
    )
    operator_join = load_json("code/source_operator_join/operator_join_packet.json")
    join_interpretation = operator_join.get("interpretation", {})
    require(
        operator_join.get("schema") == "oph.source_operator_join.v1"
        and join_interpretation.get("overlapping_pair_algebras_commute") is False
        and join_interpretation.get("observed_quantum_outcomes") is False
        and join_interpretation.get("tensor_assembly") == "declared",
        "OPERATOR_JOIN",
        "operator join packet drift",
    )
    fermion = load_json("code/sm_fermion_current/current_receipt.json")
    fermion_scope = fermion.get("scope", {})
    require(
        fermion_scope.get("hypercharge_edge_factors_on_prepared_q5_sites") is True
        and fermion_scope.get("authenticated_register_read_write_replay") is True
        and fermion_scope.get("source_selected_action_population_or_clock") is False
        and fermion_scope.get("physical_spin_signal_order_or_laboratory_identification") is False,
        "FERMION_CURRENT",
        "fermionic current receipt drift",
    )
    maxwell_readme = (REPO_ROOT / "code/maxwell_measurement/README.md").read_text(
        encoding="utf-8"
    )
    require(
        "No apparatus, physical raw data, frozen experiment or observed result is shipped."
        in maxwell_readme,
        "MAXWELL_ADAPTER",
        "Maxwell adapter shipping boundary drift",
    )

    required_text = {
        "Lean/InformationProjection/SourceHistoryPacket.lean": (
            "sourceHistoryPacket_receipt",
        ),
        "Lean/Variational/SourceToHamiltonianComposed.lean": (
            "theorem here claims that the source selects",
        ),
        "Lean/Geometry/SourceFeedbackTransport.lean": (
            "supplied operations",
        ),
        "Lean/Geometry/SourceSeamPathTomography.lean": (
            "not a selected routing",
        ),
        "Lean/Screen/A5ResponseWordAlgebra.lean": (
            "commutative",
        ),
        "Lean/Screen/OrientedFaceBracketSelector.lean": (
            "Jacobi",
        ),
        "Lean/ObserverPatchHolography/A2EndpointCommutator.lean": (
            "fifteen strict-index commutators",
        ),
        "Lean/Screen/RepairWordCarrierReadout.lean": (
            "does not emit that complete word-action packet",
        ),
        "Lean/Screen/CarrierEvolutionFlow.lean": (
            "one-parameter group",
        ),
        "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean": (
            "two serializations of the diamond",
        ),
        "Lean/Geometry/SourceEncodedMemory.lean": (
            "copying halves both",
            "no physical gain",
        ),
        "Lean/Geometry/SourceNativeRecords.lean": (
            "schedule-selection obstruction",
            "different records",
        ),
        "Lean/Geometry/SourcePassiveReset.lean": (
            "positive_register_cannot_reset",
        ),
        "Lean/Geometry/SourcePassiveMemoryBudget.lean": (
            "closed_cycle_quiescent",
            "no_catalytic_balanced_copy",
        ),
        "Lean/Geometry/SourceReadSelection.lean": (
            "finite-word selection and remote-read obstruction",
            "not impossibility of eventual",
        ),
        "Lean/Geometry/SourceConstrainedRead.lean": (
            "prevents a guaranteed correct two-input read",
        ),
        "Lean/Geometry/SourceReadRouting.lean": (
            "Custody/control order is not erased",
            "closed_sum_preserving_word_cannot_reset",
            "schedule_independent_readouts",
        ),
        "Lean/Geometry/SourceRoutingStorage.lean": (
            "No canonical pair-mean selection or arbitrary historical-read service follows.",
        ),
        "Lean/Geometry/SourceRoutingHierarchy.lean": (
            "base_diameter",
        ),
        "Lean/Geometry/SourceNativeAccumulator.lean": (
            "Only scalar pair means write the state",
        ),
        "Lean/Geometry/SourceBankBusWitness.lean": (
            "Its star connections are declared graph data",
        ),
        "Lean/Geometry/SourceBankCompiler.lean": (
            "The backend receives geometric route data.",
        ),
        "Lean/Geometry/SourceNativeStoredProgram.lean": (
            "Neither the temporal word nor that layout is selected from the source axioms.",
        ),
        "Lean/Geometry/SourceAccumulatorProgram.lean": (
            "clearing does not refresh its error allowance",
        ),
        "Lean/Geometry/SourceTemporalErasure.lean": (
            "Permanent loss before the first discriminating observation",
            "reversedWord",
        ),
        "Lean/Geometry/SourceTemporalTomography.lean": (
            "no metric menu is selected",
            "with no record payload in the schedule",
        ),
        "Lean/Geometry/SourceTemporalConnected.lean": (
            "This is an existence theorem for exact scalar readout",
        ),
        "Lean/Geometry/SourceReadAcceptance.lean": (
            "they do not derive the source grammar or an A3",
        ),
        "Lean/Geometry/SourceCheckpointPolicy.lean": (
            "the law never samples a failed trajectory and retries it",
        ),
        "Lean/Geometry/SourceCheckpointNative.lean": (
            "Selecting the actual public record interface and physical deadline remains",
        ),
        "Lean/QFT/TripleCarrierOperatorJoin.lean": (
            "the two overlapping pair algebras do not commute",
            "No physical region, clock, quantum outcome, continuum",
        ),
    }
    for relative, needles in required_text.items():
        text = (REPO_ROOT / relative).read_text(encoding="utf-8")
        require(
            all(needle.lower() in text.lower() for needle in needles),
            "LEAN_BOUNDARY",
            f"Lean boundary drift: {relative}",
        )

    publication = load_json("code/source_publication_selection/receipt.json")
    operation_capture = load_json("code/source_operation_reads/capture.json")
    operation_receipt = load_json("code/source_operation_reads/receipt.json")
    require(publication.get("M1_derived") is False and publication.get("complete_prefix_consistency") is True,
            "PUBLICATION_BOUNDARY", "conditional publication boundary drift")
    require(publication.get("noncommuting_segments_checked") == 3 and publication.get("quantum_channels_reconstructed") == 4,
            "STATE_SELECTION_CONTROLS", "finite state/process controls drift")
    require(len(operation_capture.get("cases", [])) == 7 and operation_capture.get("scope", {}).get("complete_A1_A3") is False,
            "NATIVE_READ_SCOPE", "native read capture scope drift")
    require(operation_capture.get("instruments", {}).get("extensions_selected_by_source") is False,
            "INSTRUMENT_SCOPE", "supplied instrument promoted to native source")
    require(operation_receipt.get("M1_derived") is False and operation_receipt.get("complete_A1_A3_model") is False,
            "NATIVE_READ_BOUNDARY", "native read physical scope drift")

    return {
        "port_count": len(ports),
        "response_algebra_dimension": algebra["exact_dimension"],
        "response_commutator_nonzero_count": algebra["commutator_nonzero_count"],
        "proper_recharting_count": rechart["proper_recharting_count"],
        "path_episode_count": len(path_tomography.get("episodes", [])),
        "feedback_event_count": feedback.get("summary", {}).get("events"),
        "record_presentation_permutations": presentation["permutations_checked"],
        "b14_jacobi_nonzero_count": b14["jacobi_failure"]["nonzero_count"],
        "routing_q13_events": routing_runs["q13_baseline"]["costs"]["events"],
        "routing_q21_events": routing_runs["q21_baseline"]["costs"]["events"],
        "routing_q13_hops": routing_runs["q13_baseline"]["costs"]["hops"],
        "stored_program_histories": programs.get("histories"),
        "accumulator_histories_replayed": accumulator.get("histories_replayed"),
        "temporal_captured_histories": temporal.get("captured_histories"),
        "temporal_plan_native_means": captured_plan.get("native_means"),
        "read_acceptance_histories_replayed": read_acceptance.get("histories_replayed"),
        "reusable_bus_histories": len(reusable_bus["histories"]),
        "read_selection_captured_means": read_selection.get("captured_scalar_means"),
    }


def build_candidates(facts: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return the adjudicated list.  No candidate data are combined."""

    rows = [
        candidate(
            "registered_adjacency_response",
            "COMMUTATIVE_ONLY",
            (
                "code/a5_closure/manifests/charged_response_semantic_artifact.json",
                "code/a5_closure/manifests/source_current_capability_projection.json",
                "code/a5_closure/receipts/source_current_capability.receipt.json",
            ),
            "Pinned twelve-port carrier impulse readback and its exact adjacency recurrence.",
            "STATIC_RECONSTRUCTION",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "This candidate's registered reversible response is one functional calculus in A: exact dimension 4, zero nonzero commutators, and no raw runtime history.",
            (
                f"response-word dimension = {facts['response_algebra_dimension']}",
                f"nonzero commutators = {facts['response_commutator_nonzero_count']}",
                "registered_recurrence.raw_runtime_history_serialized = false",
            ),
        ),
        candidate(
            "source_seam_path_tomography",
            "IRREVERSIBLE_ONLY",
            (
                "code/source_routing/runtime/path_tomography_receipt.json",
                "code/source_routing/README.md",
            ),
            "Captured W12 level-three support with declared scalar pair-mean calibration episodes.",
            "RAW_DESTRUCTIVE_EVENT_LOG",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It serializes exact histories, but pair averaging is destructive, paths cross many carriers, one declared schedule is used, and full routing/refinement is false.",
            (
                f"serialized episodes = {facts['path_episode_count']}",
                "all intermediate events retained",
                "README types the construction as destructive",
            ),
        ),
        candidate(
            "source_feedback_transport",
            "IRREVERSIBLE_ONLY",
            (
                "code/source_feedback_transport/transport_receipt.json",
                "code/source_feedback_transport/README.md",
            ),
            "Classical feedback episodes over the captured W12 routing support.",
            "RAW_DESTRUCTIVE_EVENT_LOG",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Mean, reset, and export overwrite state under one fixed serial schedule; the logs are transport histories, not reversible current perturbations.",
            (f"serialized events = {facts['feedback_event_count']}", "reset and mean operations present"),
        ),
        candidate(
            "source_scalar_sequential_instrument",
            "IRREVERSIBLE_ONLY",
            (
                "code/source_scalar_instruments/sequential_instrument_receipt.json",
                "code/source_scalar_finite_instrument/README.md",
                "code/source_scalar_finite_instrument/finite_instrument_receipt.json",
            ),
            "Supplied 64-site centered quantum instrument, analytic cross-time brackets, and the finite-duration pointer readout that replays it.",
            "STATIC_INSTRUMENT_WITHOUT_SAMPLED_HISTORIES",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The full instrument includes reset/readout, is 64-site rather than W12, and explicitly has no sampled outcome histories or source-selected operations; the finite-duration extension certifies pulse bounds on the same supplied action and samples no quantum outcomes.",
            (
                "sampled_quantum_outcome_histories = false",
                "source_selected_quantum_operations = false",
                "finite instrument receipt pins the sequential instrument as parent",
            ),
        ),
        candidate(
            "record_counting_repair_628",
            "IRREVERSIBLE_ONLY",
            ("code/a5_closure/manifests/record_counting_mechanism_reference.json",),
            "Operational signed record events and strict-descent repair on the pinned twelve-port carrier.",
            "RAW_APPEND_ONLY_RECORD_LOG",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "This logged candidate's append-only/strict-descent dynamics are irreversible and presentation order is explicitly quotiented by commutative addition.",
            (
                "event grammar is append-only",
                f"presentation permutations checked = {facts['record_presentation_permutations']}",
                "three coherent refinement maps are pinned",
            ),
        ),
        candidate(
            "directed_seam_repair",
            "IRREVERSIBLE_ONLY",
            ("code/a5_closure/manifests/directed_seam_repair_reference.json",),
            "Expected pair-average seam repair and representative integer microhistories.",
            "REPRESENTATIVE_REPAIR_HISTORY",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Directed placements are seam repairs, not per-port reversible perturbations; protected-record instrumentation and refinement naturality remain open.",
            ("expectation equals pair average", "refinement_semigroup_naturality = open"),
        ),
        candidate(
            "consensus_settling_packets",
            "IRREVERSIBLE_ONLY",
            ("code/consensus/manifests/compiled_lattice_settling_reference.json",),
            "Compiler scatter/project settling trajectory and its pinned upstream implementation.",
            "REPRESENTATIVE_SETTLING_TRAJECTORY",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Scatter/project settling erases state and its local arities are not the twelve-port current carrier.",
            ("settling trajectory is serialized", "scope excludes a physical current construction"),
        ),
        candidate(
            "carrier_and_discrete_refinement_scaffold",
            "STATIC_ONLY",
            (
                "code/a5_closure/manifests/echosahedral_federation_reference.json",
                "code/refinement/runtime/discrete_refinement_theorem_receipt.json",
            ),
            "Exact twelve-port carrier, port maps, and abstract mesh/refinement theorems.",
            "STATIC_COMBINATORIAL_SCAFFOLD",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "The carrier/refinement packet supplies ancestry and commuting maps, but no perturbation, inverse, or ordered source history.",
            ("twelve carrier ports", "source-native physical birefinement remains open"),
        ),
        candidate(
            "issue_566_bracket_search_space",
            "MISSING_ORDER_INFORMATION",
            (
                "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json",
                "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_stage2.receipt.json",
            ),
            "Target-free equivariant bracket coefficient space and exact Jacobi equations.",
            "STATIC_ALGEBRAIC_SEARCH_SPACE",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It parameterizes possible brackets but records no perturbation compositions and selects no bracket or physical current.",
            ("Stage 1 is search-space only", "Stage 2 full classification remains open"),
        ),
        candidate(
            "b14_oriented_face_bracket",
            "STATIC_ONLY",
            ("code/b14_jacobi/oriented_face_bracket_selector.certificate.json",),
            "Incidence-defined equal-weight bracket on twenty oriented faces.",
            "STATIC_INCIDENCE_CONSTRUCTION",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It is not extracted from histories and its Jacobi tensor is exactly nonzero, so it is a negative control rather than a current.",
            (f"Jacobi nonzero coordinates = {facts['b14_jacobi_nonzero_count']}", "Jacobi squared norm = 240"),
        ),
        candidate(
            "static_response_grammar_and_poles",
            "STATIC_ONLY",
            (
                "code/a5_closure/manifests/response_grammar_completeness_reference.json",
                "code/a5_closure/receipts/response_grammar_completeness_reference.receipt.json",
                "code/a5_closure/manifests/charged_response_pole_residue_artifact.json",
            ),
            "Bounded involution grammar, projectors, and pole/residue response tables.",
            "STATIC_OPERATOR_TABLES",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Reversible static operators are not twelve source perturbation families and no two ordered histories are serialized.",
            ("grammar forcing is independence-limited", "operator grammar is bounded and static"),
        ),
        candidate(
            "spin_deck_transport",
            "STATIC_ONLY",
            ("code/a5_closure/manifests/spin_statistics_semantic_artifact.json",),
            "Finite deck-lift/group transport tables on the carrier and refinement maps.",
            "STATIC_GROUP_TABLE",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "Group elements are reversible algebraically, but there are no per-port perturbation histories and the physical source gate is false.",
            ("finite source-model gate passes", "continuum and laboratory rows remain false"),
        ),
        candidate(
            "declared_port_current_fixture",
            "DOWNSTREAM_CONTAMINATED",
            (
                "code/a5_closure/manifests/port_current_response_reference.json",
                "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
            ),
            "Conditional charged-double-triplet matrix-current fixture.",
            "DOWNSTREAM_MATRIX_FIXTURE",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=False,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "Its generators, bracket, block decomposition, and implementers are the forbidden comparison fixture, not reconstruction inputs.",
            ("construction_model = charged_double_triplet", "response status is declared fixture"),
        ),
        candidate(
            "downstream_kinetic_form_selector",
            "DOWNSTREAM_CONTAMINATED",
            ("code/angular_sprint/runtime/kinetic_form_selection_receipt.json",),
            "Kinetic-form analysis rebuilt from the conditional matrix-current receipt.",
            "DOWNSTREAM_FIXTURE_CONSUMER",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=False,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It imports the charged-double-triplet current and explicitly leaves current-lift source selection open.",
            ("current_lift_source_selected = false", "current_fixture is charged-double-triplet"),
        ),
        candidate(
            "generic_transport_from_gluing",
            "UNRELATED",
            ("code/consensus/transport_from_gluing.py",),
            "Generic caller-supplied edge-matrix path products and holonomy obstruction.",
            "GENERIC_ABSTRACTION_NO_SOURCE_INSTANCE",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The caller supplies the matrices; there is no source-native twelve-port instance, history custody, or refinement ancestry.",
            ("generic path-product API", "no committed current packet"),
        ),
        candidate(
            "lean_pair_mean_memory_and_reusable_bus",
            "IRREVERSIBLE_ONLY",
            (
                "Lean/Geometry/SourceEncodedMemory.lean",
                "Lean/Geometry/SourceEncodedMemoryAxiomAudit.lean",
                "Lean/Geometry/SourceNonlinearRecord.lean",
                "Lean/Geometry/SourcePassiveReset.lean",
                "Lean/Geometry/SourcePassiveMemoryBudget.lean",
                "Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean",
                "Lean/Geometry/SourceReusableBus.lean",
                "Lean/Geometry/SourceReusableBusAxiomAudit.lean",
                "Lean/Geometry/SourceBusScaling.lean",
                "Lean/Geometry/SourceNativeRecords.lean",
                "Lean/Geometry/SourceNativeUpdatesAxiomAudit.lean",
                "code/source_encoded_memory/README.md",
                "code/source_encoded_memory/CONTRACT.md",
                "code/source_encoded_memory/controls.json",
                "code/source_encoded_memory/receipt.json",
                "code/source_passive_memory/README.md",
                "code/source_passive_memory/CONTRACT.md",
                "code/source_passive_memory/controls.json",
                "code/source_passive_memory/receipt.json",
                "code/source_reusable_bus/README.md",
                "code/source_reusable_bus/CONTRACT.md",
                "code/source_reusable_bus/controls.json",
                "code/source_reusable_bus/receipt.json",
                "code/source_native_updates/README.md",
                "code/source_native_updates/CONTRACT.md",
                "code/source_native_updates/controls.json",
                "code/source_native_updates/receipt.json",
            ),
            "Formal copy, clear, route, cleanup, memory, and arithmetic programs compiled to the scalar pair-mean law, with their runtime controls on captured carriers 0, 1 and 3 of the level-three support.",
            "FORMAL_IRREVERSIBLE_PAIR_MEAN_PROGRAMS_WITH_RETAINED_CONTROL_TAPES",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The programs attenuate, average, clear, and dissipate. Their two-order Fin3 control proves schedule dependence, and the runtime tapes execute one supplied request sequence per history across three carriers; they supply neither twelve reversible port families nor the 132 serialized mixed-order histories.",
            (
                "encoded copy and routing attenuate amplitudes by powers of two",
                "positive nonnegative registers cannot reset exactly under finite pair means",
                "nontrivial closed pair-mean cycles are excluded by the quadratic defect",
                "the schedule-dependence witness covers two words on Fin3 rather than all W12 pairs",
                f"reusable-bus fixed-point tapes = {facts['reusable_bus_histories']} on captured ports of carriers 0, 1 and 3",
            ),
        ),
        candidate(
            "lean_selected_pair_mean_histories",
            "IRREVERSIBLE_ONLY",
            (
                "Lean/Geometry/SourceConstrainedRead.lean",
                "Lean/Geometry/SourceConstrainedSelection.lean",
                "Lean/Geometry/SourceReadSelection.lean",
                "Lean/Geometry/SourceReadSelectionAxiomAudit.lean",
                "Lean/Geometry/SourceSelectionControls.lean",
                "code/source_read_selection/README.md",
                "code/source_read_selection/CONTRACT.md",
                "code/source_read_selection/controls.json",
                "code/source_read_selection/receipt.json",
            ),
            "Formal finite pair-mean words, prefix transcripts, information-projection selection, and remote-read obstructions, with captured controls on twelve source ports and two receiver ports of the level-three support.",
            "FORMAL_SELECTED_IRREVERSIBLE_TRANSCRIPTS_WITH_CAPTURED_CONTROLS",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Theorems quantify over abstract finite words and transcripts of destructive pair means, and the captured controls execute two source interventions under one cut. They provide no inverse laws, refinement ancestry, or complete mixed-order family.",
            (
                "full-prefix transcripts are formal functions of supplied words",
                "cut-avoiding selected words obstruct guaranteed remote reads",
                "cover choice can change selected history support",
                f"captured scalar means = {facts['read_selection_captured_means']} with status conditional_obstruction_not_full_axiom_countermodel",
            ),
        ),
        candidate(
            "lean_source_history_abstractions",
            "PARTIAL",
            (
                "Lean/InformationProjection/SourceHistoryPacket.lean",
                "Lean/Variational/SourceToHamiltonianComposed.lean",
                "Lean/QFT/SourceOperatorGeneration.lean",
                "Lean/QFT/SourceCorrelationCapstone.lean",
                "Lean/EventAlgebra/SourceBoundInstrumentInterface.lean",
            ),
            "Finite source-history packets, observed paths, and supplied Markov/instrument interfaces.",
            "FORMAL_INTERFACE_WITHOUT_RUNTIME_INSTANCE",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "These theorems type supplied or post-processed histories, but use Fin8/Fin13/abstract carriers and instantiate neither a W12 perturbation family nor mixed-order response data.",
            ("formal source-history structures exist", "source selection/runtime custody is explicitly not claimed"),
        ),
        candidate(
            "lean_conditional_a2_words",
            "COMMUTATIVE_ONLY",
            (
                "Lean/ObserverPatchHolography/A2EndpointCommutator.lean",
                "Lean/Screen/RepairWordCarrierReadout.lean",
            ),
            "Conditional six-axis/twelve-orientation step words and cumulative port-load readback.",
            "FORMAL_CONDITIONAL_WORD_ALGEBRA",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The files explicitly lack emitted step maps/diamonds and quotient the declared words by commuting relations; no runtime source packet instantiates them.",
            ("formal steps have inverses", "fifteen commutations force an abelian quotient", "source step maps are absent"),
        ),
        candidate(
            "lean_reversible_flow_models",
            "MISSING_ORDER_INFORMATION",
            (
                "Lean/QFT/SourceHistoryGNSDynamics.lean",
                "Lean/Screen/CarrierEvolutionFlow.lean",
            ),
            "Promoted unitary/GNS and one-parameter carrier-flow models.",
            "HISTORY_RECONSTRUCTED_FROM_DECLARED_GENERATOR",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The flows are mathematically reversible but reconstructed from supplied generators on M8/Fin5 or one aggregate Fin12 mode, with no twelve P_p families or mixed-order histories.",
            ("one-parameter/unitary group laws are formalized", "no per-port source perturbation packet"),
        ),
        candidate(
            "lean_declared_fin12_transport_refinement",
            "PARTIAL",
            (
                "Lean/Geometry/CarrierDynamicsCompatibility.lean",
                "Lean/Geometry/WorldlineHopTransport.lean",
                "Lean/Screen/PrimitivePortTranslationBridge.lean",
            ),
            "Declared Fin12 bijections, seam steps, translations, histories, and geometric transport/refinement interfaces.",
            "DECLARED_TRANSPORT_WITH_EXTERNAL_HISTORIES",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Only a few global symmetries/declared steps are present; histories are arguments, and neither a run-pinned source producer nor all 132 mixed orders is supplied.",
            ("Fin12 reversible maps exist", "source/readout attachment remains open"),
        ),
        candidate(
            "lean_irreversible_provenance_histories",
            "IRREVERSIBLE_ONLY",
            (
                "Lean/QFT/SourceDerivedEventPrecedence.lean",
                "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean",
                "Lean/ObserverPatchHolography/Provenance/RefinementNaturality.lean",
            ),
            "Read-from/Boolean commit histories with two serializations and a small provenance/refinement certificate.",
            "RAW_SMALL_CARRIER_IRREVERSIBLE_HISTORY",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "The Boolean commits/read-from edges are irreversible, order-independent where both serializations occur, and live on Fin2/Fin4-to-Fin3 rather than the twelve-port carrier.",
            ("two Boolean serializations are formalized", "refinement map is non-injective and supplied"),
        ),
        candidate(
            "lean_transport_and_current_boundaries",
            "STATIC_ONLY",
            (
                "Lean/Geometry/SourceFeedbackTransport.lean",
                "Lean/Geometry/SourceSeamPathTomography.lean",
                "Lean/Screen/A2HolonomyBridge.lean",
                "Lean/Screen/A5ResponseWordAlgebra.lean",
                "Lean/Screen/OrientedFaceBracketSelector.lean",
            ),
            "Kernel-checked finite statements mirroring destructive transports and static current obstructions.",
            "FORMALIZATION_OF_EXISTING_STATIC_OR_IRREVERSIBLE_OBJECTS",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Formal proofs do not create the absent runtime packet; they preserve the same destructive/static boundaries.",
            ("response-word algebra is commutative", "oriented-face bracket fails Jacobi"),
        ),
        candidate(
            "source_read_routing_full_family",
            "IRREVERSIBLE_ONLY",
            (
                "code/source_read_routing/CONTRACT.md",
                "code/source_read_routing/README.md",
                "code/source_read_routing/specification.json",
                "code/source_read_routing/controls/q3_baseline.json",
                "code/source_read_routing/controls/q3_source.json",
                "code/source_read_routing/controls/q3_branch.json",
                "code/source_read_routing/controls/q3_scratch.json",
                "evidence/source_net_causal_poset/routed_read_law/q13_baseline.json",
                "evidence/source_net_causal_poset/routed_read_law/q13_source.json",
                "evidence/source_net_causal_poset/routed_read_law/q21_baseline.json",
                "evidence/source_net_causal_poset/routed_read_law/q21_source.json",
                "code/source_routing_refinement/README.md",
                "code/source_routing_refinement/CONTRACT.md",
                "code/source_routing_refinement/receipt.json",
                "Lean/Geometry/SourceReadRouting.lean",
                "Lean/Geometry/SourceRoutingBudget.lean",
                "Lean/Geometry/SourceRoutingHierarchy.lean",
                "Lean/Geometry/SourceRoutingRefinementAxiomAudit.lean",
                "Lean/Geometry/SourceRoutingStorage.lean",
            ),
            "Full-family q=13 and q=21 metric-read compiler on captured W12 level-four and level-five supports under the retained M1 read, memory, feedback and control law, with level-three controls and a kernel-checked six-event transport word.",
            "REGENERATED_EVENT_STREAMS_AGAINST_RETAINED_COMMITMENTS",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "Every hop executes constant resets and a destructive pair mean, and the module proves that closed sum-preserving words cannot reset, so the transport is irreversible by construction. One intervention family, the centre payload +1, is executed in one order per run and the value readouts are proved schedule-independent; no mixed-order pair of perturbations is recorded and no twelve reversible families exist. The histories span thousands of carriers rather than one twelve-port carrier.",
            (
                f"q13 retained events = {facts['routing_q13_events']}",
                f"q21 retained events = {facts['routing_q21_events']}",
                f"q13 multicast seam hops = {facts['routing_q13_hops']}",
                "six-event hop word contains resetReceiver and resetSource",
                "specification.scope.axiomatic_read_law_derived = false and M1 status = retained supplied structural rule",
                "level-three controls retain complete tapes; production streams are regenerated against phase and whole-stream SHA-256 commitments",
                "refinement receipt: base diameter 3 and M1_source_selected = false",
            ),
        ),
        candidate(
            "native_stored_programs_and_accumulator",
            "IRREVERSIBLE_ONLY",
            (
                "code/source_native_programs/README.md",
                "code/source_native_programs/CONTRACT.md",
                "code/source_native_programs/controls.json",
                "code/source_native_programs/receipt.json",
                "code/source_native_programs/families_receipt.json",
                "code/source_native_accumulator/README.md",
                "code/source_native_accumulator/CONTRACT.md",
                "code/source_native_accumulator/controls.json",
                "code/source_native_accumulator/receipt.json",
                "Lean/Geometry/SourceAccumulatorAxiomAudit.lean",
                "Lean/Geometry/SourceAccumulatorDecoder.lean",
                "Lean/Geometry/SourceAccumulatorProgram.lean",
                "Lean/Geometry/SourceBankBusWitness.lean",
                "Lean/Geometry/SourceBankCompiler.lean",
                "Lean/Geometry/SourceBankControl.lean",
                "Lean/Geometry/SourceBankExecution.lean",
                "Lean/Geometry/SourceBankInvariant.lean",
                "Lean/Geometry/SourceBankLowering.lean",
                "Lean/Geometry/SourceBankMachine.lean",
                "Lean/Geometry/SourceBankPrecision.lean",
                "Lean/Geometry/SourceBankTrace.lean",
                "Lean/Geometry/SourceNativeAccumulator.lean",
                "Lean/Geometry/SourceNativeCore.lean",
                "Lean/Geometry/SourceNativeProgramBudget.lean",
                "Lean/Geometry/SourceNativeProgramError.lean",
                "Lean/Geometry/SourceNativeProgramsAxiomAudit.lean",
                "Lean/Geometry/SourceNativeShuttle.lean",
                "Lean/Geometry/SourceNativeStoredProgram.lean",
            ),
            "Native accumulation, shuttle, stored-program and two-bank compiler words built from scalar pair means on captured seams, with a kernel-checked captured control and replayed request families.",
            "RETAINED_CONTROL_TAPES_WITH_REGENERATED_FAMILIES",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Start, addition, retirement, transfer and cleanup are pair means with retirement and residual clearing, and both receipts restrict dynamic writes to pair means only. The instruction stream is payload-independent and executes one supplied request sequence per history; no inverse law, mixed-order pair family or twelve port families is present, and the declared bus is not embedded in the captured W12 graph.",
            (
                "stored-program receipt.dynamic_writes = supported_scalar_pair_means_only",
                "accumulator receipt.dynamic_scalar_writes = pair_means_only",
                f"accumulator histories replayed = {facts['accumulator_histories_replayed']}",
                f"stored-program histories = {facts['stored_program_histories']}",
                "m1_derived = false in both receipts",
                "SourceBankBusWitness: star connections are declared graph data",
            ),
        ),
        candidate(
            "native_temporal_tomography_and_checkpoint_selection",
            "IRREVERSIBLE_ONLY",
            (
                "code/source_temporal_acceptance/README.md",
                "code/source_temporal_acceptance/CONTRACT.md",
                "code/source_temporal_acceptance/controls.json",
                "code/source_temporal_acceptance/receipt.json",
                "code/source_read_acceptance/README.md",
                "code/source_read_acceptance/CONTRACT.md",
                "code/source_read_acceptance/controls.json",
                "code/source_read_acceptance/receipt.json",
                "code/source_checkpoint_selection/README.md",
                "code/source_checkpoint_selection/CONTRACT.md",
                "code/source_checkpoint_selection/controls.json",
                "code/source_checkpoint_selection/receipt.json",
                "Lean/Geometry/SourceCheckpointAxiomAudit.lean",
                "Lean/Geometry/SourceCheckpointEntropy.lean",
                "Lean/Geometry/SourceCheckpointMeaning.lean",
                "Lean/Geometry/SourceCheckpointNative.lean",
                "Lean/Geometry/SourceCheckpointPipeline.lean",
                "Lean/Geometry/SourceCheckpointPolicy.lean",
                "Lean/Geometry/SourceReadAcceptance.lean",
                "Lean/Geometry/SourceReadAcceptanceAxiomAudit.lean",
                "Lean/Geometry/SourceReadAcceptanceSchedule.lean",
                "Lean/Geometry/SourceTemporalAcceptanceAxiomAudit.lean",
                "Lean/Geometry/SourceTemporalAttempts.lean",
                "Lean/Geometry/SourceTemporalConnected.lean",
                "Lean/Geometry/SourceTemporalErasure.lean",
                "Lean/Geometry/SourceTemporalEssential.lean",
                "Lean/Geometry/SourceTemporalGuard.lean",
                "Lean/Geometry/SourceTemporalMenu.lean",
                "Lean/Geometry/SourceTemporalNative.lean",
                "Lean/Geometry/SourceTemporalObservation.lean",
                "Lean/Geometry/SourceTemporalSelection.lean",
                "Lean/Geometry/SourceTemporalTomography.lean",
            ),
            "Linear observation grammar, sound partial publication, connected native completion plans and checkpoint selection on the captured 15,360-port level-three support, built from scalar seam means.",
            "RETAINED_CONTROLS_WITH_STRUCTURAL_PLAN_CHECK",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The words are destructive pair means, and the erasure theorem proves permanent loss of a prepared difference after its first mean. The separate and reversed four-mean words on one seam path are a two-order witness for one pair of unknowns, not 132 ordered compositions of twelve reversible families; the large calibration plan is checked structurally without payload replay.",
            (
                "temporal controls: full_axiom_instantiation = false",
                f"captured plan native means = {facts['temporal_plan_native_means']} with payload_replay = false",
                f"captured histories = {facts['temporal_captured_histories']}",
                f"read-acceptance histories replayed = {facts['read_acceptance_histories_replayed']} with source_contract = conditional_not_source_selected",
                "checkpoint scope.canonical_axiom_selection = false",
                "SourceTemporalErasure: permanent loss before the first discriminating observation",
            ),
        ),
        candidate(
            "finite_source_operator_join",
            "STATIC_ONLY",
            (
                "code/source_operator_join/README.md",
                "code/source_operator_join/operator_join_packet.json",
                "Lean/QFT/TripleCarrierOperatorJoin.lean",
            ),
            "Operator embeddings of two committed pair algebras into the triple-carrier matrix algebra, intertwining 31 retained source-row transitions.",
            "DETERMINISTIC_MATHEMATICAL_PACKET",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The packet is a static algebra join: the noncommuting hinge is an operator-algebra fact about supplied embeddings, tensor assembly and path-read transpositions are declared postprocessing, and no perturbation family, ordered history or inverse law is executed.",
            (
                "interpretation.overlapping_pair_algebras_commute = false",
                "interpretation.tensor_assembly = declared",
                "interpretation.observed_quantum_outcomes = false",
            ),
        ),
        candidate(
            "fermionic_hypercharge_current_histories",
            "DOWNSTREAM_CONTAMINATED",
            (
                "code/sm_fermion_current/README.md",
                "code/sm_fermion_current/coupled_README.md",
                "code/sm_fermion_current/current_receipt.json",
                "code/sm_fermion_current/coupled_receipt.json",
                "code/sm_fermion_current/quantum_link_receipt.json",
            ),
            "Declared CAR hopping factors on the 64-site golden carrier with registered multiplets, retaining register read and write histories.",
            "AUTHENTICATED_REGISTER_REPLAY",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=False,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The histories carry gauge, multiplet and hypercharge labels under a supplied action, so the source firewall excludes them; they live on the 64-site carrier rather than the twelve-port carrier and record one baseline, one phase intervention and one gauge copy rather than mixed-order port perturbations.",
            (
                "scope.hypercharge_edge_factors_on_prepared_q5_sites = true",
                "scope.authenticated_register_read_write_replay = true",
                "scope.source_selected_action_population_or_clock = false",
            ),
        ),
        candidate(
            "maxwell_measurement_adapter",
            "DOWNSTREAM_CONTAMINATED",
            (
                "code/maxwell_measurement/README.md",
                "code/maxwell_measurement/contract.py",
            ),
            "Capture schema and preregistration contract for a fixed Maxwell emulator with destructive pair-average probes and restoring feedback.",
            "NO_CAPTURE_SHIPPED",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=False,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "No capture, apparatus data or observed result is shipped, and the adapter identifies encoded registers with model electromagnetic potentials, so any capture would enter with a downstream target label.",
            (
                "README: no apparatus, physical raw data, frozen experiment or observed result is shipped",
                "three ordered controls and two gauge executions are schema fields, not retained histories",
            ),
        ),
        candidate(
            'native_eventual_publication',
            'IRREVERSIBLE_ONLY',
            ('Lean/Geometry/SourcePublicationAxiomAudit.lean',
             'Lean/Geometry/SourcePublicationLaw.lean',
             'Lean/Geometry/SourcePublicationMass.lean',
             'Lean/Geometry/SourcePublicationMenu.lean',
             'Lean/Geometry/SourcePublicationNative.lean',
             'Lean/Geometry/SourcePublicationPath.lean',
             'Lean/Geometry/SourcePublicationPathInvariant.lean',
             'Lean/Geometry/SourcePublicationTail.lean',
             'Lean/Geometry/SourceRadiusRemoval.lean',
             'Lean/Geometry/SourceRadiusSelection.lean',
             'code/source_publication_selection/CONTRACT.md',
             'code/source_publication_selection/DERIVATION.md',
             'code/source_publication_selection/README.md',
             'code/source_publication_selection/__init__.py',
             'code/source_publication_selection/build.py',
             'code/source_publication_selection/check_exterior.py',
             'code/source_publication_selection/check_locality.py',
             'code/source_publication_selection/codec.py',
             'code/source_publication_selection/controls.json',
             'code/source_publication_selection/exterior.py',
             'code/source_publication_selection/locality.py',
             'code/source_publication_selection/receipt.json',
             'code/source_publication_selection/test_selection.py',
             'code/source_publication_selection/verify.py'),
            ('Eventual-publication and radius selection on a declared connected scalar pair-mean grammar with '
             'captured-support controls; this is not the registered all-port matching driver.'),
            'PINNED_FINITE_CONTROLS_AND_REVIEWED_ANALYTIC_SCOPE',
            properties(source_native=True, port_indexed=True, reversible=False, both_composition_orders_recorded=False, raw_histories_serialized=True, target_free=True, same_twelve_port_carrier=False, refinement_provenance_present=False, twelve_reversible_perturbation_families=False),
            'The conditional connected pair-mean grammar can be order-sensitive but its means are destructive. Publication conditioning and a supplied demand/read law do not instantiate twelve native reversible perturbation families or same-current overlap holonomy.',
            ('complete_prefix_consistency = true; M1_derived = false', 'connected scalar mean grammar and demand interface are declared', 'uniform native block instantiation is analytic, not a finite-census theorem'),
        ),
        candidate(
            'finite_state_process_selection',
            'STATIC_ONLY',
            ('Lean/Geometry/SourceStateSelection.lean',
             'code/source_publication_selection/STATE_SELECTION.md',
             'code/source_publication_selection/check_quantum.py',
             'code/source_publication_selection/quantum.py',
             'code/source_publication_selection/test_quantum.py',
             'code/source_publication_selection/receipt.json'),
            ('Finite-algebra A3 support theorems, exact noncommuting state segments and supplied identity, '
             'depolarization and record-preserving channels.'),
            'PINNED_FINITE_CONTROLS_AND_REVIEWED_ANALYTIC_SCOPE',
            properties(source_native=False, port_indexed=False, reversible=False, both_composition_orders_recorded=False, raw_histories_serialized=False, target_free=True, same_twelve_port_carrier=False, refinement_provenance_present=False, twelve_reversible_perturbation_families=False),
            'The finite A3 theorem identifies maximal feasible local support and positive failure floors. The process controls use a second declared Choi-state optimization problem. Neither constructs the source feasible grammar or records mixed-order native perturbations.',
            ('noncommuting_segments_checked = 3; quantum_channels_reconstructed = 4', 'selected state does not select a record-preserving process', 'no global extending state is assumed by the finite support theorem'),
        ),
        candidate(
            'registered_operation_read_histories',
            'IRREVERSIBLE_ONLY',
            ('Lean/Geometry/SourceOperationReads.lean',
             'Lean/Geometry/SourceOperationReadsAxiomAudit.lean',
             'code/source_operation_reads/DERIVATION.md',
             'code/source_operation_reads/README.md',
             'code/source_operation_reads/__init__.py',
             'code/source_operation_reads/capture.json',
             'code/source_operation_reads/check_live.py',
             'code/source_operation_reads/receipt.json',
             'code/source_operation_reads/simulator_verifier.py',
             'code/source_operation_reads/source_snapshot.json',
             'code/source_operation_reads/test_reads.py',
             'code/source_operation_reads/verify.py'),
            ('Exact-pinned simulator histories of the registered all-port driver, its native unitary, retained '
             'local ledger reads, and disjoint port-pair repair.'),
            'PINNED_FINITE_CONTROLS_AND_REVIEWED_ANALYTIC_SCOPE',
            properties(source_native=True, port_indexed=True, reversible=False, both_composition_orders_recorded=False, raw_histories_serialized=True, target_free=True, same_twelve_port_carrier=False, refinement_provenance_present=False, twelve_reversible_perturbation_families=False),
            'The matching means are commuting idempotent projections and are irreversible. The icosahedral Laplacian unitary is one prescribed native evolution, not twelve source-derived port flows. Coordinate probes and a time-translated probe generate M12(C), but those noninvertible probes and their algebraic closure do not supply reversible perturbation families. History extension is not spatial refinement.',
            ('seven native captured cases; exact source closure and local numeric records', 'even N >= 4 full-sweep normalized rank is 6N-1, with kernel 5N+1', 'accessible probe algebra is M12(C); this differs from the four-dimensional icosahedral adjacency-response candidate'),
        ),
        candidate(
            'declared_flagged_recovery_extensions',
            'PARTIAL',
            ('Lean/Geometry/SourceRepairInstrument.lean',
             'code/source_operation_reads/INSTRUMENT_DERIVATION.md',
             'code/source_operation_reads/capture.json',
             'code/source_operation_reads/receipt.json'),
            ('Explicitly supplied swap-twirl/dephased channels, fixed recorded unitary branches and a '
             'constrained linear archive extension of mean repair.'),
            'PINNED_FINITE_CONTROLS_AND_REVIEWED_ANALYTIC_SCOPE',
            properties(source_native=False, port_indexed=True, reversible=True, both_composition_orders_recorded=False, raw_histories_serialized=False, target_free=True, same_twelve_port_carrier=False, refinement_provenance_present=False, twelve_reversible_perturbation_families=False),
            'Here reversible=true refers only to each fixed recorded swap word, which is unitary and has a two-sided inverse on the quantum state. The complete flagged channel has recovery only on its image with retained state and 6N bits; it is not an automorphism of the enlarged output algebra. The deterministic archive needs 5N+1 extra real coordinates and global recovery on its constrained normalized image. Unconditional twirl and dephased repair remain irreversible. None of these extra interfaces is selected by the registered source or supplies twelve mixed-order native families.',
            ('extensions_selected_by_source = false', 'conditional flag recovery retains quantum state and 6N bits', 'linear archive dimension 5N+1 is a real-coordinate count, not a finite-bit or local decoder'),
        ),
    ]
    require(
        len({row["candidate_id"] for row in rows}) == len(rows),
        "CANDIDATE_ID",
        "duplicate candidate id",
    )
    return rows


def integrated_tree_review(
    table: Sequence[tuple[str, str, str | None, str]] = INTEGRATED_TREE_REVIEW,
) -> list[dict[str, Any]]:
    rows = [
        {
            "path": path,
            "decision": decision,
            "candidate_id": candidate_id,
            "reason": reason,
        }
        for path, decision, candidate_id, reason in table
    ]
    require(
        [row["path"] for row in rows]
        == sorted({row["path"] for row in rows}),
        "INTEGRATED_REVIEW_ORDER",
        "integrated-tree review paths must be unique and canonical",
    )
    require(
        all(
            row["decision"]
            in {
                "NEW_CANDIDATE",
                "EXTENDS_EXISTING_CANDIDATE",
                "OUTSIDE_REVIEWED_SOURCE_SCOPE",
            }
            for row in rows
        ),
        "INTEGRATED_REVIEW_DECISION",
        "invalid integrated-tree semantic decision",
    )
    return rows


def build_inventory() -> dict[str, Any]:
    facts = audit_source_semantics()
    candidates = build_candidates(facts)
    qualifying = [
        row
        for row in candidates
        if all(row["properties"][field] for field in QUALIFICATION_FIELDS)
    ]
    require(not qualifying, "STOP_CONDITION", "a qualifying candidate requires Stage-2 review")
    counts = Counter(row["classification"] for row in candidates)
    inventory: dict[str, Any] = {
        "schema": SCHEMA,
        "program": "SOURCE-CURRENT-TOMOGRAPHY-0",
        "stage": 2,
        "previous_inventory_base_sha": PREVIOUS_INVENTORY_BASE_SHA,
        "audited_upstream_main_sha": UPSTREAM_MAIN_SHA,
        "verdict": VERDICT,
        "verdict_scope": (
            "The verdict applies only to the explicitly enumerated candidates "
            "and reviewed source scope at the audited upstream base."
        ),
        "positive_stages_authorized": False,
        "physical_current_source_bridge_attained": False,
        "fixture_used_as_reconstruction_oracle": False,
        "target_labels_used": False,
        "audit_scope": {
            "rule": "one candidate packet must satisfy every qualification field; facts from distinct packets are not composited",
            "directories": list(AUDITED_DIRECTORIES),
            "directories_added_since_initial_inventory": [
                {"directory": directory, "reason": reason}
                for directory, reason in AUDIT_SCOPE_EXTENSION
            ],
            "repository_wide_search_terms": [
                "port-indexed perturbation",
                "reversible perturbation",
                "ordered response/history/composition",
                "mixed-order",
                "raw runtime history serialized",
                "sequential instrument",
                "commutator",
                "overlap word",
                "recharting",
                "refinement provenance",
            ],
        },
        "audited_file_snapshot": audited_file_snapshot(),
        "integrated_tree_review": {
            "from_upstream_main_sha": PREVIOUS_INVENTORY_BASE_SHA,
            "through_upstream_main_sha": UPSTREAM_MAIN_SHA,
            "surface_count": len(INTEGRATED_TREE_REVIEW),
            "surfaces": integrated_tree_review(),
            "prior_reviews": [
                {
                    "from_upstream_main_sha": start,
                    "through_upstream_main_sha": end,
                    "surface_count": len(table),
                    "surfaces": integrated_tree_review(table),
                }
                for start, end, table in PRIOR_INTEGRATED_TREE_REVIEWS
            ],
        },
        "qualification_fields": list(QUALIFICATION_FIELDS),
        "source_pins": [source_pin(path) for path in AUDITED_SOURCE_PATHS],
        "implementation_pins": [
            source_pin("code/a5_closure/source_current_order_sensitive_inventory.py"),
            source_pin("code/a5_closure/verify_source_current_order_sensitive_inventory.py"),
        ],
        "candidates": candidates,
        "summary": {
            "candidate_count": len(candidates),
            "qualifying_candidate_count": len(qualifying),
            "classification_counts": dict(sorted(counts.items())),
            "registered_adjacency_response_summary": {
                "candidate_id": "registered_adjacency_response",
                "response_word_dimension": facts["response_algebra_dimension"],
                "commutator_nonzero_count": facts[
                    "response_commutator_nonzero_count"
                ],
                "proper_rechartings_total": facts["proper_recharting_count"],
                "proper_rechartings_in_response_words": 1,
            },
            "logged_repair_candidate_summary": {
                "candidate_id": "record_counting_repair_628",
                "fatal_failures": [
                    "irreversible strict-descent repair",
                    "presentation order erased by commutative addition",
                    "no twelve reversible perturbation families",
                ],
            },
            "routed_read_candidate_summary": {
                "candidate_id": "source_read_routing_full_family",
                "fatal_failures": [
                    "constant resets and destructive pair means in every hop",
                    "one intervention family in one order per run with schedule-independent readouts",
                    "no twelve reversible perturbation families",
                ],
            },
            "cross_packet_noncomposition": (
                "The reversible four-dimensional response packet, the logged "
                "record-repair packet and the regenerated full-family routing packet "
                "have different primitives and semantics. "
                "Combining their complementary fields would be a new producer, not "
                "an inventory finding."
            ),
        },
        "minimal_missing_fields": [
            {
                "field": "perturbation_families",
                "requirement": "exactly twelve reversible source-native families P_0,...,P_11, one per primitive port",
            },
            {
                "field": "ordered_mixed_histories",
                "requirement": "both p-then-q and q-then-p for all 66 unordered distinct-port pairs (132 ordered compositions)",
            },
            {
                "field": "raw_history_custody",
                "requirement": "canonical serialized raw histories with provenance and tamper-evident hashes, not reconstruction from one static generator",
            },
            {
                "field": "common_source_binding",
                "requirement": "all histories share the same twelve-port carrier, orientation, normalization, source family, and refinement stage/tower",
            },
            {
                "field": "reversibility_witness",
                "requirement": "an inverse/two-sided law for every perturbation family, excluding reset, averaging, measurement, and strict-descent repair",
            },
            {
                "field": "source_firewall",
                "requirement": "no current fixture, gauge/particle label, measurement target, or downstream selector enters the producer",
            },
        ],
        "claim_boundary": {
            "bounded_repository_inventory_only": True,
            "verdict_limited_to_enumerated_reviewed_source_scope": True,
            "does_not_assert_a_no_go_for_future_source_data": True,
            "does_not_weaken_the_abstract_forced_lie_type_theorem": True,
            "does_not_reject_the_conditional_port_current_fixture": True,
            "does_not_make_a_physical_current_claim": True,
            "stages_3_through_6_not_run": True,
            "missing_packet": "No enumerated candidate supplies the complete reversible ordered-history contract; existing static and irreversible packets do not supply it by composition.",
        },
        "verifier_command": (
            "python3 code/a5_closure/verify_source_current_order_sensitive_inventory.py "
            "--inventory code/a5_closure/manifests/source_current_order_sensitive_inventory.json"
        ),
    }
    inventory["inventory_sha256"] = canonical_sha256(inventory)
    return inventory


def verify_inventory(inventory: Mapping[str, Any]) -> None:
    require(inventory.get("schema") == SCHEMA, "INVENTORY_SCHEMA", "wrong schema")
    body = {key: value for key, value in inventory.items() if key != "inventory_sha256"}
    require(
        inventory.get("inventory_sha256") == canonical_sha256(body),
        "INVENTORY_HASH",
        "inventory self-hash failed",
    )
    require(
        dict(inventory) == build_inventory(),
        "INVENTORY_REPLAY",
        "committed inventory does not replay exactly",
    )


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("all", "show", "verify"))
    parser.add_argument("--inventory", type=Path, default=INVENTORY_PATH)
    args = parser.parse_args()

    if args.command == "all":
        inventory = build_inventory()
        write_json(args.inventory, inventory)
        print(VERDICT)
        return 0
    if args.command == "show":
        print(json.dumps(build_inventory(), indent=2, sort_keys=True))
        return 0
    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    require(isinstance(inventory, dict), "JSON_OBJECT", "inventory is not an object")
    verify_inventory(inventory)
    print("SOURCE_CURRENT_ORDER_SENSITIVE_INVENTORY_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
