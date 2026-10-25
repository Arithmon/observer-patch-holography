# Particle Provenance Ledger

Generated: `2026-10-25T12:33:45Z`

This ledger records target-use and convention-sensitivity status for the public quantitative particle rows.

## Scientific Boundary

- Promotion allowed: `False`
- Public rows classified: `True`
- Numeric sensitivity intervals supplied: `False`
- Boundary: The ledger classifies target ancestry, blindness, and declared convention sensitivity. It does not supply the source spectral measure or interval-composition certificates required for a numerical sensitivity sweep.

## Rows

| Particle | Value | Class | Blind status | Target use | Promotable | Convention sensitivity |
| --- | ---: | --- | --- | --- | --- | --- |
| `higgs` | `125.1995304097179 GeV` | `conditional_declared_surface_candidate` | `conditionally_blind_on_declared_surface` | candidate_upstream_d10_repair_not_source_promoted | `False` | depends_on_declared_D10_D11_running_matching_threshold_surface |

## Withheld Non-Prediction Rows

These rows have source artifacts but no public prediction value in the output tables.

| Particle | Claim label | Blind status | Target use | Reason |
| --- | --- | --- | --- | --- |
| `electron` | `exact_target_anchored_current_family_witness` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `muon` | `exact_target_anchored_current_family_witness` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `tau` | `exact_target_anchored_current_family_witness` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `up_quark` | `selected_class_target_anchored_mixed_convention_mass_texture_audit` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `charm_quark` | `selected_class_target_anchored_mixed_convention_mass_texture_audit` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `top_quark` | `selected_class_target_anchored_mixed_convention_mass_texture_audit` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `down_quark` | `selected_class_target_anchored_mixed_convention_mass_texture_audit` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `strange_quark` | `selected_class_target_anchored_mixed_convention_mass_texture_audit` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `bottom_quark` | `selected_class_target_anchored_mixed_convention_mass_texture_audit` | `withheld_not_blind` | target_values_or_target_derived_datum_used | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `electron_neutrino` | `rejected_target_informed_weighted_cycle_candidate` | `withheld_not_blind_rejected_candidate` | target_ranked_selector_development_and_correlated_profile_rejection | target_informed_candidate_rejected_by_correlated_profile |
| `muon_neutrino` | `rejected_target_informed_weighted_cycle_candidate` | `withheld_not_blind_rejected_candidate` | target_ranked_selector_development_and_correlated_profile_rejection | target_informed_candidate_rejected_by_correlated_profile |
| `tau_neutrino` | `rejected_target_informed_weighted_cycle_candidate` | `withheld_not_blind_rejected_candidate` | target_ranked_selector_development_and_correlated_profile_rejection | target_informed_candidate_rejected_by_correlated_profile |

## Separated Classical Carrier Modes

These zero hard quadratic parameters are branch-conditional mode statements, not public quantum-particle mass predictions.

| Carrier | Hard parameter squared | Classical gate | Quantum gate | Particle promotion |
| --- | ---: | --- | --- | --- |
| `photon` | `0` | `conditional_pass_on_declared_action_phase_branch` | `not_passed` | `False` |
| `gluon` | `0` | `conditional_pass_on_declared_action_phase_branch` | `not_passed` | `False` |
| `graviton` | `0` | `conditional_pass_on_declared_action_phase_branch` | `not_passed` | `False` |

## Prospective Comparison Workflows

- `new_quantity_pre_reference_provenance`: `protocol_defined__not_exercised`. For any quantitative row outside construction inputs, timestamp and hash the source artifacts, record allowed conventions, then fetch or reveal the external reference. Required evidence: `source_artifact_hashes`, `forbidden_target_inputs`, `convention_set`, `pre_reference_runtime_output`, `post_reference_comparison_only_delta`.
- `convention_sensitivity_sweep`: `taxonomy_declared__certified_intervals_required`. Vary only declared scheme, matching, and threshold choices inside certified intervals; report induced intervals for every public quantitative row. Required evidence: `scheme_record`, `threshold_map`, `matching_interval_composition_certificate`, `rowwise_sensitivity_intervals`.

## Convention Sensitivity

- Classification: `declared_taxonomy__numeric_sweep_not_performed`
- RG contract classification: `open_source_rg_frontier_partial`
- Endpoint contract classification: `source_spectral_reduction_closed_measure_payload_absent`
- Endpoint package classification: `None`
- Missing evidence: populated source spectral measure payload and interval-composition certificates

## Restricted finite matter and field contracts

The response-derived internal multiplets, invariant interactions and shared field operators below are finite witnesses under declared contracts. They supply no observed mass value and promote none of the quantitative rows.

| Witness | Premise connections | Finite supply | Receipt |
|---|---|---|---|
| matter_decoder | PR-59 | Exact conditional fifteen-state charge/Casimir decoder and all source-current brackets. | `model/observer_particles/results/representation_identification.json` |
| one_higgs_vertices | PR-50, PR-67 | One-doublet invariant lines and all 648 family/current/component identities; coefficients are unselected. | `model/observer_particles/results/unified_yukawa_action.json` |
| spin_matter_controls | PR-47 | Independent finite Spin and matter controls with complete readout; physical chirality is unselected. | `model/observer_particles/results/spin_matter_join.json` |
| shared_spatial_operator | PR-53, PR-71 | Complete conditional CAR field shares the ideal scalar spatial generator and incident coordinate; physical field and refinement are separate. | `model/observer_particles/results/shared_spatial_matter.json` |
| full_cochain_control | PR-53, PR-71 | Complete cone cochains preserve the scalar block and remove truncated edge cycles; this cone is distinct from the radial cylinder. | `model/observer_particles/results/cochain_completion.json` |
| metric_and_null_dictionary | PR-25, PR-53 | Common pointwise tetrad, field symbols and rank-nine null dictionary; independent Einstein balances are not supplied. | `model/observer_particles/results/geometric_matter_join.json` |
| native_entropy_and_area | PR-26, PR-27 | Finite constrained entropy and geometric area/volume identities; physical vacuum fixed-volume stationarity is not supplied. | `model/observer_particles/results/geometric_stationarity.json` |
| entangled_pair_scalar | PR-17, PR-54 | Wissner-Gross's finite spectral scalar, independently reconstructed with complete outcome controls; physical Thomson and same-cell identifications are separate. | `model/observer_particles/results/entangled_alpha.json` |
