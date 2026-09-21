# Reading only when the local evidence distinguishes the answer

A receiver can publish a value with certainty when every source state
compatible with its observations gives that value. This package derives that
partial readout rule and implements a noisy two-payload read using native
means on captured seams. Waiting and aborted histories retain their full cost.
The general path result relates delivery probability to a declared random
schedule and finite precision. Deriving the preparation, permitted schedules
and publication instrument from the source axioms is the main remaining step.

[Objective and exit](CONTRACT.md), [derivation and boundaries](DERIVATION.md),
[independent receipt](receipt.json).

Run from the repository root:

```
python code/source_read_acceptance/build.py
python code/source_read_acceptance/verify.py --write-receipt
python code/source_read_acceptance/verify.py
python -m pytest -q code/source_read_acceptance
cd Lean
lake build Geometry.SourceReadAcceptanceAxiomAudit
```

The verifier reconstructs every stopped history. The retained examples and
aggregate commitments are not substitutes for replay. The receipt establishes
a conditional finite read service; it does not classify M1 as derived.
