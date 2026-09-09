# Closure loop at carrier scale

This package mirrors the simulator receipt
`data/exact/closure_loop_receipt.json` (producer `oph_exact/closure_loop.py`
in `oph-physics-sim`, with an independent verifier and a test file) and the
fifteen stored event logs of the single-carrier loops.

## What the loop does

A generator runs an exact carrier under a repair law from seeded integer loads
and emits an event log only: the initial readings, each event's changed ports
with their readings before and after, the final readings, and a probe log of
readbacks after repeated repair from one-hot states. A fixed recovery
algorithm reads that log and nothing else. From it the recovery returns the
port count, the seam set (the port pairs that change together), the
inverse-port pairing at graph distance three, the incidence automorphism group
and its orientation-preserving subgroup, the repair rule class, the schedule
law, the normalized response Gram with its rank, the Lie-type dimension count,
the descent statistics and the terminal law. The recovered specification is
instantiated, run with a fresh seed, and recovered again. Equality of the
invariant vectors of the inhabited structure, the constructed realization and
the second iteration is the fixed-point receipt at that scale.

## Result

For the canonical source (icosahedral carrier, seam-mean law, uniform schedule)
the recovered specification is twelve ports, thirty seams, degree five
throughout, the antipode at distance three, automorphism order 120, rotation
order 60, rule class `seam_mean` (every event conservative and symmetric),
seam chi-square 181/6 against the uniform threshold at significance 0.001,
Gram eigenvalues `(4, 4, 4, 0, ...)` at the largest probe step with rank
three, the Lie split `1 + 3 + 8`, zero strict-descent violations of the
squared norm, and a terminal state equal to the initial mean across eight
schedules. The three invariant vectors agree:
`CLOSURE_FIXED_POINT_AT_CARRIER_SCALE`. The isolated federation of twenty
carriers recovers as twenty disjoint copies with the same component
invariants.

Negative controls: the tetrahedron and the octahedron recover their own
specifications (four and six ports, groups of order 12 and 24, Gram rank
three for both, Lie splits `1 + 3` and none), so the Gram rank alone selects
nothing and the pairing with the rotation order does; the overwrite law is
recovered as nonconservative with a schedule-dependent terminal state and no
fixed point; the integer nearest-agreement law reaches the fixed point on the
quotient invariants with the tie placed uniformly.

The descent potential is the flagship's `V = sum N_i^2`, which every
conservative unit transfer across a seam with oriented mismatch `d >= 2`
lowers by exactly `2(d - 1)`; the seam-sum form rises on 75 events of the
canonical log and is recorded as a diagnostic.

## Independent check

```bash
python3 evidence/closure_loop/verify_closure_loop_archive.py
```

The checker imports no simulator code. It verifies every manifest digest,
strict and canonical JSON, the receipt's self digest and schema, the log pins,
the invariant-vector equality flags of every loop, the negative-control flags,
and recomputes from the canonical source log alone the seam set, the degree
sequence, the antipode, the rule class of every event, the seam chi-square,
the conserved total, the strict descent of the centered squared norm, the
seam-form increases, and the probe readback as the exact power of the
recovered expectation operator. Recorded output:

```text
VERIFIED_CLOSURE_LOOP_ARCHIVE, 16 archived files verified, 15 logs recovered and compared
```

## Boundary

One carrier and the isolated federation of twenty carriers. The orientation
class, the uniform schedule, the initial loads and the Lie classification
table are declared. Universe-level closure, existence or uniqueness of the
cosmic fixed point, the physical identification of any invariant, and the
glued federation are outside this package; the glued federation is work in
progress.
