# Conditional Candidate Outputs

Generated: `2026-10-26T06:30:05Z`

This surface renders the numeric conditional candidates that the strict
source-only promotion policy withholds from the public OPH-value columns.
Each row carries its recorded claim label and its open gates. No row is
promotable, and this surface must not feed back into promotion decisions.
The public claim boundary stays in `RESULTS_STATUS.md`; target-anchored
exact fits stay in `EXACT_FITS_ONLY.md`.

The electroweak rows are the conditional D10 value law. No measured mass
enters, but every chart value is evaluated at the working pixel
`P_C = 1.630968209403959`, which is located from the measured
fine-structure constant; the forward closure root
`P_fwd = 1.630972095858897` is a different number. The chain runs the
gauge couplings with the declared one-loop coefficients `(33/5, 1, -3)`
and `beta_EW = 4`. The QT1-QT5 quotient-path source-entailment
certificate is not emitted, so the value follows from those premises
only conditionally. The diagnostic base-running chart and the
measured-reference inverse adapter stay in the technical audit README
and off this surface.

## Conditional GeV Chart Candidates

| Observable | Value | Chart | Open gates | Source |
| --- | ---: | --- | --- | --- |
| `m_H` | `125.1995304097179 GeV` | `declared_d10_d11_running_matching_threshold_surface` | `d10_target_free_repair_closure, qt1_qt5_quotient_path_source_entailment, ds1_ds5_d11_split_character_and_rigidity, physical_complex_pole_conversion` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |
| `m_t` | `172.3523553288312 GeV` | `declared_d10_d11_running_matching_threshold_surface` | `d10_target_free_repair_closure, qt1_qt5_quotient_path_source_entailment, ds1_ds5_d11_split_character_and_rigidity, physical_complex_pole_conversion, top_threshold_control` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |
| `M_W` | `80.37700001539531 GeV` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `M_Z` | `91.18797807794321 GeV` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |

## Conditional Scale-Free Chart Ratios

These ratios cancel the chart scale. The value-law triple shares one
emitted quintet, so the three coordinates stand or fall together. The
Higgs/top ratio shares the declared-surface split pair.

| Observable | Value | Chart | Open gates | Source |
| --- | ---: | --- | --- | --- |
| `M_W / M_Z` | `0.8814429457652062` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `sin2_theta_W_eff` | `0.22305833336075578` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `alpha_em_eff_inv` | `132.7524088191668` | `archived_d10_value_law_candidate` | `qt1_qt5_quotient_path_source_entailment, physical_complex_pole_attachment` | `code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json` |
| `m_H / m_t` | `0.726416126839982` | `declared_d10_d11_running_matching_threshold_surface` | `d10_target_free_repair_closure, ds1_ds5_d11_split_character_and_rigidity, physical_complex_pole_conversion` | `code/particles/runs/calibration/d11_live_exact_split_pair_theorem.json` |

## Conditional Dimensionless Candidates

| Observable | Value | Chart | Open gates | Source |
| --- | ---: | --- | --- | --- |
| `v_over_E_star` | `2.0199803239725553e-17` | `closed_form_hierarchy_law_on_public_endpoint_packet` | `e_star_physical_normalization, clock_scale_certificate, frozen_running_threshold_convention` | `code/particles/hierarchy/computations/hierarchy_numeric_witness.json` |
