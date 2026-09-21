# Reusable native accumulation

A mutable accumulator computes `1 + sum(requested input values)` using scalar
pair means after one preparation. A four-mean rectangle clears the accumulator
and workspace, a two-mean export reuses the unit seed, and each addition aligns
the two dyadic scales before combining them. The input values remain decodable
at their updated scales. The accumulator's previous value is retired at each
start.

[Stored computation](../source_native_programs/README.md) uses two payload
banks so that computed results can be inputs to subsequent layers, with native
transfer, retirement and a whole-history precision budget.

The [contract](CONTRACT.md) states the objective, deliverables, exit and retained
premises. The [derivation](DERIVATION.md) connects the instruction sequence to
the Lean theorems; the [receipt](receipt.json) records independent replay.

The captured instance covers 1,280 complete histories and 78,255 scalar means,
including all request words of lengths zero through seven over the unit and
payload inputs, five signed integer payloads, and 15-episode reuse. Full tapes
are regenerated; the checked-in evidence contains compact commitments and
three complete examples. The final receiver is on a neighbouring carrier.

Run from the repository root:

```
python code/source_native_accumulator/build.py
python code/source_native_accumulator/verify.py --write-receipt
python code/source_native_accumulator/verify.py
python -m pytest -q code/source_native_accumulator
cd Lean
lake build Geometry.SourceAccumulatorAxiomAudit
```

The verifier independently constructs each word, executes its means on exact
rationals and on the rounded grid, checks the captured seams and current
writers, and checks the public interval decoder at every declared checkpoint.
The CLI requires the complete fixed fixture inventory. It is not a validator
for arbitrary submitted programs. The Lean result quantifies over arbitrary
finite valid request lists and input values.

This realizes the local arithmetic interface used by the conditional M1
compiler. The preparation, supported placement, request sequence, metadata
and publication instrument are supplied. The source-selection criterion and
cover obstruction of [the selection package](../source_read_selection/README.md)
apply independently. A mean-only implementation of a request does not select
the request from A1--A3. No M1 classification is promoted.
