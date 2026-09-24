import Mathlib.Util.AssertNoSorry
import Geometry.M1QuantumTransport

/-! Transitive standard-axiom audit of the finite quantum-transport reductions. -/

open Lean Elab Command in
elab "audit_m1_transport " n:ident : command => do
  let name ← liftCoreM <| Lean.Elab.realizeGlobalConstNoOverloadWithInfo n
  let axioms ← Lean.collectAxioms name
  let allowed := #[``propext, ``Classical.choice, ``Quot.sound]
  let unexpected := axioms.filter (fun ax => !allowed.contains ax)
  unless unexpected.isEmpty do
    throwError "M1 transport axiom audit rejected {unexpected.toList}"
  logInfo m!"'{name}' depends on axioms: {axioms.toList}"

/-- error: M1 transport axiom audit rejected [sorryAx] -/
#guard_msgs in
audit_m1_transport sorryAx

/-- error: M1 transport axiom audit rejected [Lean.ofReduceBool, Lean.trustCompiler] -/
#guard_msgs in
audit_m1_transport Lean.ofReduceBool

audit_m1_transport OPH.M1QuantumTransport.dot_deficit
audit_m1_transport OPH.M1QuantumTransport.weighted_trace_budget
audit_m1_transport OPH.M1QuantumTransport.isotropic_speed
audit_m1_transport OPH.M1QuantumTransport.shared_speed_impossible
audit_m1_transport OPH.M1QuantumTransport.charged_three_flights
audit_m1_transport OPH.M1QuantumTransport.overhead_lowers_speed
audit_m1_transport OPH.M1QuantumTransport.three_axis_remainder
audit_m1_transport OPH.M1QuantumTransport.charged_error
audit_m1_transport OPH.M1QuantumTransport.pauli_current
audit_m1_transport OPH.M1QuantumTransport.zero_interaction_coefficients
audit_m1_transport OPH.M1QuantumTransport.saturated_weighted_bound
audit_m1_transport OPH.M1QuantumTransport.finite_stencil_gap
audit_m1_transport OPH.M1QuantumTransport.direction_resource_bound
