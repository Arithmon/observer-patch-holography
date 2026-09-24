import Mathlib.Util.AssertNoSorry
import Geometry.M1VacuumFidelity

/-! Transitive standard-axiom audit of the finite vacuum-fidelity reductions. -/

open Lean Elab Command in
elab "audit_m1_vacuum " n:ident : command => do
  let name ← liftCoreM <| Lean.Elab.realizeGlobalConstNoOverloadWithInfo n
  let axioms ← Lean.collectAxioms name
  let allowed := #[``propext, ``Classical.choice, ``Quot.sound]
  let unexpected := axioms.filter (fun ax => !allowed.contains ax)
  unless unexpected.isEmpty do
    throwError "M1 vacuum axiom audit rejected {unexpected.toList}"
  logInfo m!"'{name}' depends on axioms: {axioms.toList}"

/-- error: M1 vacuum axiom audit rejected [sorryAx] -/
#guard_msgs in
audit_m1_vacuum sorryAx

/-- error: M1 vacuum axiom audit rejected [Lean.ofReduceBool, Lean.trustCompiler] -/
#guard_msgs in
audit_m1_vacuum Lean.ofReduceBool

audit_m1_vacuum OPH.M1VacuumFidelity.one_step_excitation
audit_m1_vacuum OPH.M1VacuumFidelity.second_order_symplectic
audit_m1_vacuum OPH.M1VacuumFidelity.second_order_excitation
audit_m1_vacuum OPH.M1VacuumFidelity.invariant_form
audit_m1_vacuum OPH.M1VacuumFidelity.composition_duration
audit_m1_vacuum OPH.M1VacuumFidelity.cubic_cancellation
audit_m1_vacuum OPH.M1VacuumFidelity.bulk_fraction
audit_m1_vacuum OPH.M1VacuumFidelity.repair_energy_scale
audit_m1_vacuum OPH.M1VacuumFidelity.repair_number_scale
audit_m1_vacuum OPH.M1VacuumFidelity.second_order_energy_scale
