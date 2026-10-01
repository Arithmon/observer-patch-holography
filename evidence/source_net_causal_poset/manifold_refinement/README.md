# Finite manifold refinement interface

This package reconciles retained observer-like patches with local port state,
read events, records and finite observation regions under the supplied
source-net read law. It performs arithmetic on existing receipts; it runs no
new source population or large sampled experiment.

From the research repository root:

```sh
python3 evidence/source_net_causal_poset/manifold_refinement/build.py
python3 -O evidence/source_net_causal_poset/manifold_refinement/verify.py
python3 -m pytest -q evidence/source_net_causal_poset/manifold_refinement/test_verify.py
```

`build.py` generates `receipt.json` and `refinement.csv`. The CSV is the
per-statistic comparison table; the JSON additionally retains every signed
CDF-grid difference, three-level descriptive comparison and all 120 outcomes
of the separate small clock control. There are 96 regions across the three
sampled levels and all three spatial dimensions. Centre, moving and every
off-centre direction remain separate; absent one-dimensional moving regions
are not fabricated. `verify.py` imports no generator code. It checks frozen
input hashes, complete region/statistic coverage, independently evaluated
flat coefficients, arithmetic, table agreement and the stated uncertainty
boundary. Corruption tests include omitted controls, discarded clock failures
and invented precision.

Chain coefficients and the interval mean have analytic flat references.
CDFs and second moments use the retained seeded Monte Carlo references; a
25-point maximum is not a continuous KS statistic. Marginal chain errors are
reported historical estimates, not reconstructed sampling errors. Missing
joint covariance and self-normalized mean/CDF uncertainty stay missing.
Finite count-volume deviations are geometric discretization readings, not
sampling error estimates. Empty intervals remain lattice-sensitive controls.

The verdict describes the retained point estimates. It rejects a statement
of uniform monotone improvement across these readings; it neither establishes
nor statistically rejects continuum convergence. Ratios compare absolute
deviations only when the preceding deviation is nonzero. No independent-sample
or joint-significance interpretation is attached to them. The receipt also
specifies the data needed for a future uncertainty-qualified bounded control,
without supplying those samples.

Input receipts and their original sampled draws are different objects. The
draws and the q144 chunk-parallel driver are unavailable. Graph digest ties
and recorded small-level identities do not recover them. This interface
does not select the population/read law, identify a physical clock or remove
the preparation and volume-normalization inputs.
