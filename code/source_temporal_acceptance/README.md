# Temporal acceptance and complete native continuation

This package derives an exact continuation criterion and a native protocol
that completes every recoverable read. On any finite connected scalar
support, one receiver can reconstruct the entire initial state using only
supported pair means and its own retained samples. Connected topology
constructs the protocol; a successful route or decoder is not a premise.
The [contract](CONTRACT.md) states the objective, deliverables and exit;
the [derivation](DERIVATION.md) connects the result to the A3 support theorem
and audits the remaining M1 source-selection obligations.

For N ports, the general Lean construction uses N receiver samples and at
most N(N-1)/2 native means. The captured 15,360-port support has a completely
checked topology plan using 425,979 means, with maximum path length 51.
This large plan is checked structurally; payload execution is separately
replayed on eight five-port histories, including every basis preparation.
No claim of large-plan payload replay or attainable precision is made.

At an existing checkpoint, combine the retained samples L with the current
complete state response M. The connected completion theorem proves:

```
a native completion can recover f
    iff ker [L; M] is contained in ker f.
```

This distinguishes a record that has not arrived yet from one that has
been irreversibly lost. It supplies sufficiency as well as an obstruction.
No future target payload enters the test or completion schedule. Exact
sampling, retained classical records and decoder arithmetic are explicit
parts of this scalar interface, not a full A1--A3 instantiation.

The same criterion gives a maximal preservation guard: admit a proposed
mean exactly when the proposed joint kernel preserves the source records.
The Lean proof derives one bounded completing word that this guard accepts
from every protected checkpoint. An occurrence amid arbitrary other
proposals therefore completes all records. With an explicitly supplied
uniform full proposal law, a block of B proposals from M event identities
gives failure at most `(1-M^(-B))^k` after k blocks. Rejections remain in the
denominator and proposal cost; no successful-history conditioning is used.
This derives progress for that specified guard/proposal interface, without
claiming that full canonical A3 selects its proposal law or physical guard.

For response matrix L and requested linear meaning f:

```
sound publication on every feasible sample
    iff ker L is contained in ker f
    iff f belongs to the row span of L
    iff a decoder c exists with cL=f.
```

Every negative certificate supplies an invisible intervention z with Lz=0
and f(z)=1. Observation forms are constructed from preparation, scalar
pair means and local sample ports without the requested answer. The
general Lean theorem covers all real linear meanings. Coefficient
identities in the finite executable evidence are checked independently.

Under convex finite full-history KL selection with faithful reference and
positive weights, certain publication requires the same criterion on
every history admitted by any feasible law. A favorable optimizer cannot
erase an admitted ambiguity. For direct record observations, the theorem
forces inclusion of every nonzero coefficient of a specified affine
meaning. It does not select that meaning or forbid additional reads.
For arbitrary nonlinear meanings on a finite independent product domain,
the unique minimal direct-read set consists of coordinates whose isolated
variation can change the meaning. That stronger statement is also kernel
checked. Correlated preparation domains require their own fiber test.

The captured W12 controls distinguish three native words:

| Word | Receiver result |
| --- | --- |
| Merge sources, then transport | Their sum is identifiable; neither separate source is identifiable |
| Transport the second source before merging | Both sources are identifiable from two retained samples |
| Execute the path edges in reverse temporal order | Neither source nor their sum is identifiable |

The seven-mean positive word has receiver deviations y3=x1/8 and
y7=x0/16+5*x1/32. Thus x1=8*y3 and x0=16*y7-20*y3. These identities are
kernel checked for arbitrary real payloads on 15,360 scalar coordinates;
actual edge membership is verified against the pinned captured support.
Decoder arithmetic, sample retention and physical error limits are supplied.

An additional kernel theorem covers arbitrary infinite continuations. If
a first mean erases the whole state distinction before it is observed,
every subsequent receiver transcript agrees. A sound policy cannot repair
this loss by waiting. On the stated preparation, the first mean of ports
0 and 1 erases their difference. Under the separately supplied uniform
first-seam law this has probability 1/46050. Encodings or observations that
preserve the distinction change the hypotheses; protected-bank execution
is not excluded.

The compact receipt authenticates 27 captured histories, 135 scalar means
and 72 prefix/target certificates. A separate three-edge chain exhausts
9,841 attempt words through horizon eight. At horizon eight, the first
source is identifiable in 673/6561 words and their sum in 2929/6561.
Expected stopped work is 17281/2187 and 15454/2187 means respectively.
The complete census includes failed attempts, initial sampling and every
executed prefix. It is an exact conditional schedule calculation, not a
physical frequency measurement.
Each horizon also separates currently identified, recoverable by native
completion, and permanently lost meanings against the same full denominator.
The guarded horizon-eight census completes both records in 857/6561 words
and permanently loses them in none. It charges 51,750 stopped proposals,
45,445 executed means, 6,305 rejected proposals and 52,006 samples across all
6,561 attempts. The six-proposal completing block gives the conservative
analytic expectation bound of 4,374 proposals under the product proposal
law; this bound is not a measured physical time or a tight mean.

M1's metric radius and particular public meaning are not selected by this
result. Native recoverability on a connected finite support extends to every
port. The existing causal-limit theorem also admits a range of mesoscopic
radii. The derivation documents why neither fact selects the specified
metric menu; it does not claim a countermodel to the full axioms.

## Reproduction

From the repository root, with `code` on `PYTHONPATH`:

```
python -m source_temporal_acceptance.build
python -m source_temporal_acceptance.verify --write-receipt
python -m source_temporal_acceptance.verify
python -m pytest -q code/source_temporal_acceptance
```

From `Lean/`:

```
lake build Geometry.SourceTemporalAcceptanceAxiomAudit
```

The producer performs row elimination. The verifier checks primal/dual
identities against independently replayed basis preparations; the exhaustive
two-source controls use an independent determinant test. Source custody
covers the local transitive proof imports, source specification and capture.
Adversarial tests exercise the actual verifier CLI and resealed semantic
forgeries. No expanded attempt-word tapes are committed.
