# Sampled manifold reproduction interface

The four files named in the original sampled-run receipt are preserved with
matching SHA256 hashes: the source net, sampled producer, independent verifier
and producer tests. A small arithmetic dependency and namespace initializers
make them executable without a simulator checkout. `manifest.json` identifies
which hashes come directly from the historical run and which package files are
new. These are records of observer-like patches with local ports, read events,
causal ordering and bounded observation regions under a supplied read law.

From the RER root:

```sh
python3 evidence/source_net_causal_poset/verify_sampled_manifold.py
python3 evidence/source_net_causal_poset/verify_sampled_manifold.py --replay-references --write-receipt
```

The first command verifies pins and large-receipt arithmetic, checks the
sampled estimators against exact small orders and exercises corruption tests.
The second also regenerates the seeded continuum reference samples and writes
a separate verification receipt. No large source-net run is performed.

The original q55/q89 event samples are absent. Their standard errors are
retained reported quantities, not independently reconstructed here. Joint
C2/C3/C4 sample covariance, weighted-CDF uncertainty and the discarded second
C3 estimator cannot be recovered from the available aggregates. The weighted
CDF comparison uses its declared 25-point grid. These boundaries survive a
successful arithmetic and small-order verification.
