# Encoded native memory operations toward M1

## Objective

Derive a concrete, reusable encoded copy/reset protocol from the existing
scalar pair mean, without protected export, assignments after initialization,
payload-dependent ancillas or a supplied raw-zero reset. Charge attenuation
and give a finite noise margin. This investigates the restricted-code opening;
selection of the physical source law is outside its scope.

## Deliverables

1. Kernel-checked native words for copying an unknown balanced input to blank
   rails and clearing a used pair to its encoded blank state. Track every
   register, including the attenuated source.
2. Composition theorems for repeated reads of one retained logical sign and
   for moving that sign along a declared two-rail walk, including revisits.
3. A bound connecting errors in native operations and terminal readout to the
   decreasing sign margin; a sharp terminal ambiguity control. Do not equate
   exact sign preservation with unlimited finite-precision memory.
4. Independently replayed rational executions, actual writer/value custody,
   operation and storage counts, mutation tests, transitive standard-axiom
   audit and CI coverage.

## Exit

A finite classical encoded protocol, with the actual native operations and
their limits proved and checked. A two-rail cell is prepared as `[b+a,b-a]`;
blank is `[b,b]`. Means halve amplitude when copying to a blank cell. A native
mean clears the old balanced pair. This is logical reset, not raw-zero reset,
and logical copying, not copying an analog amplitude unchanged.

Locality is conditional on the declared two-rail edges. A captured W12 control
must check its actual seams; it must not be called a full support embedding or
a distant route. Initial encoding, accepted schedule, isolation, comparator,
addresses, versions, precision and physical time remain supplied. No source
population, physical energy, M1 archive service, amplitude refresh, quantum
instrument or common-world realization is derived.
