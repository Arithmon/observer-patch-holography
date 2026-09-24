import Mathlib.Util.AssertNoSorry
import Geometry.M1Necessity

/-!
# M1 necessity proof audit

All public declarations of Geometry.M1Necessity are checked transitively.
The stencil and continuum arguments remain analytic, as stated in that module.
-/

open Lean Elab Command in
elab "audit_m1_necessity " n:ident : command => do
  let name ← liftCoreM <| Lean.Elab.realizeGlobalConstNoOverloadWithInfo n
  let axioms ← Lean.collectAxioms name
  let permitted := #[``propext, ``Classical.choice, ``Quot.sound]
  let unexpected := axioms.filter (fun ax => !permitted.contains ax)
  unless unexpected.isEmpty do
    throwError "M1 necessity axiom audit rejected {unexpected.toList}"
  logInfo m!"'{name}' depends on axioms: {axioms.toList}"

/-- error: M1 necessity axiom audit rejected [sorryAx] -/
#guard_msgs in
audit_m1_necessity sorryAx

/-- error: M1 necessity axiom audit rejected [Lean.ofReduceBool, Lean.trustCompiler] -/
#guard_msgs in
audit_m1_necessity Lean.ofReduceBool

audit_m1_necessity OPH.M1Necessity.cube_face_count
audit_m1_necessity OPH.M1Necessity.binary_route
audit_m1_necessity OPH.M1Necessity.axis_route
audit_m1_necessity OPH.M1Necessity.route_budget
audit_m1_necessity OPH.M1Necessity.strict_layer_increase
audit_m1_necessity OPH.M1Necessity.adjacent_layers_require_edge
audit_m1_necessity OPH.M1Necessity.sum_lengths
audit_m1_necessity OPH.M1Necessity.sum_squared_lengths
audit_m1_necessity OPH.M1Necessity.axial_coordinate_cut
audit_m1_necessity OPH.M1Necessity.axial_diagonal_cut
audit_m1_necessity OPH.M1Necessity.sparse_dense_cut_separation
audit_m1_necessity OPH.M1Necessity.scalar_calibration_cancels
audit_m1_necessity OPH.M1Necessity.isotropic_moment_normalization
audit_m1_necessity OPH.M1Necessity.normalized_fourth_remainder
audit_m1_necessity OPH.M1Necessity.full_radius_tick_instability
audit_m1_necessity OPH.M1Necessity.balanced_crossing_scale
audit_m1_necessity OPH.M1Necessity.removed_fraction_enclosure
