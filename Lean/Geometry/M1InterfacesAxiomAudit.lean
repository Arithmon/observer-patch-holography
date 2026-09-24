import Mathlib.Util.AssertNoSorry
import Geometry.M1Interfaces

/-! Transitive standard-axiom audit of every public interface reduction. -/

open Lean Elab Command in
elab "audit_m1_interfaces " n:ident : command => do
  let name ← liftCoreM <| Lean.Elab.realizeGlobalConstNoOverloadWithInfo n
  let axioms ← Lean.collectAxioms name
  let permitted := #[``propext, ``Classical.choice, ``Quot.sound]
  let unexpected := axioms.filter (fun ax => !permitted.contains ax)
  unless unexpected.isEmpty do
    throwError "M1 interface axiom audit rejected {unexpected.toList}"
  logInfo m!"'{name}' depends on axioms: {axioms.toList}"

/-- error: M1 interface axiom audit rejected [sorryAx] -/
#guard_msgs in
audit_m1_interfaces sorryAx

/-- error: M1 interface axiom audit rejected [Lean.ofReduceBool, Lean.trustCompiler] -/
#guard_msgs in
audit_m1_interfaces Lean.ofReduceBool

audit_m1_interfaces OPH.M1Interfaces.crossing_cover
audit_m1_interfaces OPH.M1Interfaces.cut_card_le
audit_m1_interfaces OPH.M1Interfaces.absolute_cut_change
audit_m1_interfaces OPH.M1Interfaces.normalized_volume_correction
audit_m1_interfaces OPH.M1Interfaces.compactness_scale_identity
audit_m1_interfaces OPH.M1Interfaces.alias_action_scale_identity
audit_m1_interfaces OPH.M1Interfaces.sum_of_saturations
audit_m1_interfaces OPH.M1Interfaces.explicit_repair_factors
