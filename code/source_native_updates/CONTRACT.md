# Native writes and finite-chain scaling

## Objective

Remove two concrete restrictions of the existing native-memory construction:
records must all be prepared in advance, and cleanup is proved only for a
four-cell bus. These sub-results are not a derivation of a native record
service.

## Deliverables

1. Kernel proofs that a native write produces a sum record while preserving
   the two logical inputs at declared smaller scales; native halving and
   scale identities; all-word raw-range and unbounded-growth obstructions.
2. General finite-chain cleanup, including its actual compiled word,
   untouched archives, zero-edge boundary, weighted contraction, and a
   polynomial sufficient scalar-mean count for binary cleanup accuracy.
3. A captured-support execution of seven means creating `C = A+B`, followed
   by reads `C,A,C,B` across carriers, each with 80 bus cleanup sweeps. All
   16 inputs in `{-3/2,-1/2,1/2,3/2}^2` must decode; old values and the new
   version must survive subsequent operations in the declared logical sense.
4. Independent exact coefficient, fixed-point, consumed-writer, ancestry,
   resource and noise-bound checks. Retain one complete representative trace
   (`A=1/2,B=-3/2`) and compact commitments and results for the other fifteen;
   both producer and verifier independently reconstruct every event of every
   history. A missing or substituted history must fail.
5. A supported two-schedule control showing that the candidate local mean
   law alone does not determine the read history. This is not an A1--A3
   countermodel. State precisely what is not established for source selection.

## Exit

The native sum and rereads pass for all 16 input pairs with positive certified
receiver margins, including independently bounded preparation error `2^-18`,
per-mean disturbance `2^-28`, grid rounding at `Q=2^20`, and local readout
error `2^-16`. Retained histories execute rounding, not the other disturbances.
The source support and every touched seam are checked. Altered programs,
records, operations, writers, optimistic bounds, histories, pins and invalid
JSON must be rejected even if their digests are recomputed. Every theorem in this packet
must pass a transitive audit allowing only standard Lean axioms. No `sorry`,
compiler-trust proof, invented premise or weakened old test is permitted.

## Boundary

This is a finite native write/commit/read experiment with three records, not
an arbitrary mutable record server. Commit is a declared publication boundary
after the write word; its controller and version metadata are supplied.
The sum is computed by the physical means, not by the receiver or preparation.
Raw amplitudes change; decoding uses the declared record scale and read count.
The two inputs, codebooks, scalar mean law, isolation, controller, embedding
and error budgets are supplied. The physical origin of bounded errors is open.

The general chain theorem assumes two rail edges per hop and one reset rung
at the first cell. It does not produce such an embedding for every route in
the captured graph. More distant routes and general positive affine record
programs, full q=13/q=21 native histories, intervention-preserving compilation,
and the derivation of the metric-neighbour read law are not supplied by this
package.
