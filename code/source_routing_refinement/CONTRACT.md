# M1 routing refinement and live storage

## Objective

Derive useful implementation requirements of the supplied routed family instead
of extrapolating the q=13 and q=21 receipts. Separate what follows from its
program and support construction from what requires a separate source law.

## Deliverables

1. A constructive routing bound assembled analytically from four
   kernel-checked components, with captured W12 levels compared against that
   construction. The all-level graph identifications and composed conclusion
   are not Lean theorems.
2. A proof and independently replayed realization of safe local storage reuse
   for the layered read program. Account separately for scalar slots, scalar
   precision, controller state and retained audit history. The kernel layout
   bound is `23C`; `14C+9n` is the finite occupied-host census and a Lean
   arithmetic corollary.
3. An explicit raw-operation-count comparison with the logical population.
   Do not assign the logical count-volume theorem to added routing events.
4. Kernel axiom checks, adversarial controls, a reproducible finite receipt,
   and an exact statement of the remaining M1 assumptions.

## Exit

Proved bounds and a checked finite realization, or a named obstruction to a
specific proposed interface. A conditional implementation is not source
selection. No law, geometry, physical time, constant-bit memory bound or
physical event measure is inferred from a software storage allocation.

The live-storage interface serves the specified preceding-layer reads. An
arbitrary request for an older retired version is outside that interface and
must be rejected, not silently redirected to a newer value. Every original
event and consumed-writer identifier remains in the audit comparison.
