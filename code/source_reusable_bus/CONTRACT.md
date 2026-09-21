# Captured reusable record transport

## Objective

Construct a finite-precision classical record protocol between distinct
carriers of the captured W12 support. Derive transport and approximate bus
cleanup from pair means, retain two prepared records, and reread one after
transporting a different record. No protected export or assigned reset is
an operation in the ideal protocol.

## Deliverables

1. A kernel-checked paired-operation compiler, the four-cell bus cleanup
   formula and contraction, repeated cleanup, record-retention error and
   receiver error bounds, with a transitive standard-axiom audit.
2. An injective embedding of every participating rail into the captured
   support. Source records are on carrier 0; receiver rails are on carrier 3.
   Every executed mean must be an actual captured seam.
3. Complete fixed-point executions for all 16 combinations of two prepared
   records from `{-3/2,-1/2,1/2,3/2}`, with request sequence `[0,1,0]`.
   Baseline is 2, grid denominator is `2^20`, and each read is followed by
   80 cleanup sweeps. Retain every mean, input, output and consumed writer.
4. Independent exact verification of the ideal linear maps and finite
   rounding errors, decoding margins, intervention comparisons, storage and
   operation counts, mutation controls, CI coverage and an explicit audit.

## Exit

All requested payloads decode from the receiver's two local rails and the
declared codebook/read-count interface, with a strictly positive certified
margin. The complete protocol reuses the same 12 scalar registers, executes
1704 means per history, and retains the cleanup operations and their ancestry.
The ideal cleanup residual is at most `(3/2)*(7/8)^80`; finite implementation
error is charged separately. An insufficient-precision control must fail the
margin, and changed edges, writers, operations and optimistic bounds must be
rejected independently of the producer.

## Boundary

The selected support, preparation, finite codebook, isolation, schedule,
addresses, version labels, comparator and error limits are supplied. Rounded
means are an explicitly bounded implementation model of the exact mean law;
no physical rounding mechanism or measured noise bound is derived. The
analysis allows initial max-coordinate error `2^-18`, additional disturbance
`2^-28` per mean on every register, and terminal readout error `2^-16` per
rail. Retained histories execute rounding, not those additional disturbances.

Stored analog amplitudes attenuate and acquire bounded residual errors.
Logical payload reuse is not immutable scalar storage or unit-gain export.
Approximate cleanup preserves raw cross-record influence and causal ancestry.
No full record-metric compiler, physical causal-order identification, clock,
population production, quantum implementation or common-world join is claimed.
The M1 acceptance contracts are not closed by this packet.
