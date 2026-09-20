import Mathlib.Util.AssertNoSorry
import Geometry.SourceReusableBus

/-! Every theorem is audited transitively, including proof dependencies. -/
open Lean Elab Command in
elab "audit_reusable_bus_axioms " n:ident : command => do
  let name ← liftCoreM <| Lean.Elab.realizeGlobalConstNoOverloadWithInfo n
  let axioms ← Lean.collectAxioms name
  let permitted := #[``propext, ``Classical.choice, ``Quot.sound]
  let unexpected := axioms.filter (fun ax => !permitted.contains ax)
  unless unexpected.isEmpty do
    throwError "Reusable bus axiom audit rejected {unexpected.toList}"
  logInfo m!"'{name}' depends on axioms: {axioms.toList}"

/-- error: Reusable bus axiom audit rejected [sorryAx] -/
#guard_msgs in
 audit_reusable_bus_axioms sorryAx

/-- error: Reusable bus axiom audit rejected [Lean.ofReduceBool, Lean.trustCompiler] -/
#guard_msgs in
 audit_reusable_bus_axioms Lean.ofReduceBool

audit_reusable_bus_axioms OPH.SourceReusableBus.instruction_native
audit_reusable_bus_axioms OPH.SourceReusableBus.program_native
audit_reusable_bus_axioms OPH.SourceReusableBus.compile_append
audit_reusable_bus_axioms OPH.SourceReusableBus.scrub_formula
audit_reusable_bus_axioms OPH.SourceReusableBus.forward_receiver
audit_reusable_bus_axioms OPH.SourceReusableBus.scrub_native
audit_reusable_bus_axioms OPH.SourceReusableBus.scrub_length
audit_reusable_bus_axioms OPH.SourceReusableBus.forward_length
audit_reusable_bus_axioms OPH.SourceReusableBus.scrub_contraction
audit_reusable_bus_axioms OPH.SourceReusableBus.scrub_archives
audit_reusable_bus_axioms OPH.SourceReusableBus.cleaning_length
audit_reusable_bus_axioms OPH.SourceReusableBus.readCycle_length
audit_reusable_bus_axioms OPH.SourceReusableBus.execute_append
audit_reusable_bus_axioms OPH.SourceReusableBus.cleaning_native
audit_reusable_bus_axioms OPH.SourceReusableBus.clean_bound
audit_reusable_bus_axioms OPH.SourceReusableBus.clean_archives
audit_reusable_bus_axioms OPH.SourceReusableBus.export_error
audit_reusable_bus_axioms OPH.SourceReusableBus.retained_error
audit_reusable_bus_axioms OPH.SourceReusableBus.receiver_error
audit_reusable_bus_axioms OPH.SourceReusableBus.contrast_error
audit_reusable_bus_axioms OPH.SourceReusableBus.program_readout_error
audit_reusable_bus_axioms OPH.SourceReusableBus.separated_readouts
audit_reusable_bus_axioms OPH.SourceReusableBus.roundedHalf_integer_error
audit_reusable_bus_axioms OPH.SourceReusableBus.roundedHalf_error
audit_reusable_bus_axioms OPH.SourceReusableBus.roundedHalf_between
