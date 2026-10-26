# Why the current chiral candidate misses the mass hierarchy

This is a new diagnostic layer over the retained
[Weyl–Higgs candidate](../chiral_candidate/README.md). It does not modify that
action, its quantum records or the original research note. The
[research diagnosis](../../../docs/research/OPH_HIERARCHY_DIAGNOSIS_2026-10-26.md)
gives the derivations and the next physical obligations.

The calculation separates exact source algebra, retained observer readout,
native statistical/dynamical controls, and empirical comparison. Its main
exclusion applies to **every complex affine family law** `Y=uI+vF`, with
canonical kinetic normalization and a common Higgs factor. It does not
exclude arbitrary OPH dynamics. Arbitrary functions of a single feature can
change its eigenvalues but retain its family axes; that freedom is not a
mass prediction.

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
