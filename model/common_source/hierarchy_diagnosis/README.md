# Why the current chiral candidate misses the mass hierarchy

This diagnostic audits the retained
[Weyl–Higgs candidate](../chiral_candidate/README.md). It does not modify that
action or its quantum records. The
[research diagnosis](../../../docs/research/OPH_PARTICLE_MASS_RESPONSE.md#hierarchy-diagnosis)
gives the derivations and the next physical obligations.

The calculation separates exact source algebra, retained observer readout,
native statistical/dynamical controls, and empirical comparison. Its main
exclusion applies to **every complex affine family law** `Y=uI+vF`, with
canonical kinetic normalization and a common Higgs factor. It does not
exclude arbitrary OPH dynamics. Arbitrary functions of a single feature can
change its eigenvalues but retain its family axes; that freedom is not a
mass prediction.

The adjacent-sector calculation applies the existing native endpoint rule
to binary loads zero through three. It restores the two CAR contractions
about the total-two reference but proves that the recorded C/E dynamics
still has no cross-port charged response. An exact family of faithful
references agrees on all total-two data while changing adjacent-sector
modular thresholds. The chosen reference, coherent continuation and physical
field map therefore have distinct selection obligations. These are finite
operator and source-domain results; no Standard Model mass is inferred.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/common_source/hierarchy_diagnosis
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.common_source.hierarchy_diagnosis.audit --output /tmp/fresh-hierarchy-diagnosis.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.common_source.hierarchy_diagnosis.audit --verify model/common_source/hierarchy_diagnosis/results.json
```

The producer refuses to overwrite a report. Verification recomputes this
diagnosis and checks source custody and the prior retained field/clock
weights. It **does not replay the prior source evolution**; that replay
remains available through the original candidate's commands. No empirical
mass is used to change a source coefficient, and no observed hierarchy is
reported as derived.
Exact claims, rational strings and source hashes must match exactly. Floating
diagnostics are compared with a relative tolerance of `2e-9` and absolute
tolerance of `5e-11` across the pinned local and CI numerical dependencies;
these are replay tolerances, not physical error bounds.
