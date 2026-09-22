import Mathlib.Util.AssertNoSorry
import Geometry.RecordGluingPrinciple

/-! Transitive whitelist audit, including negative controls for proof shortcuts. -/
open Lean Elab Command in
elab "audit_rg_principle " n:ident : command => do
  let name ← liftCoreM <| Lean.Elab.realizeGlobalConstNoOverloadWithInfo n
  let axioms ← Lean.collectAxioms name
  let permitted := #[``propext, ``Classical.choice, ``Quot.sound]
  let unexpected := axioms.filter (fun ax => !permitted.contains ax)
  unless unexpected.isEmpty do
    throwError "RG principle axiom audit rejected {unexpected.toList}"
  logInfo m!"'{name}' depends on axioms: {axioms.toList}"

/-- error: RG principle axiom audit rejected [sorryAx] -/
#guard_msgs in
 audit_rg_principle sorryAx

/-- error: RG principle axiom audit rejected [Lean.ofReduceBool, Lean.trustCompiler] -/
#guard_msgs in
 audit_rg_principle Lean.ofReduceBool

audit_rg_principle OPH.RecordGluingPrinciple.finite_symmetry_quadratic
audit_rg_principle OPH.RecordGluingPrinciple.graded_path_bound
audit_rg_principle OPH.RecordGluingPrinciple.refinement_clock_bound
audit_rg_principle OPH.RecordGluingPrinciple.subdivided_time_formula
audit_rg_principle OPH.RecordGluingPrinciple.finite_positive_overhead
audit_rg_principle OPH.RecordGluingPrinciple.exists_subdivision_before
audit_rg_principle OPH.RecordGluingPrinciple.lossless_injective
audit_rg_principle OPH.RecordGluingPrinciple.transcript_capacity
audit_rg_principle OPH.RecordGluingPrinciple.independent_records_capacity
audit_rg_principle OPH.RecordGluingPrinciple.reversible_evaluation_involutive
audit_rg_principle OPH.RecordGluingPrinciple.reversible_evaluation_blank
audit_rg_principle OPH.RecordGluingPrinciple.dense_inner_cube
