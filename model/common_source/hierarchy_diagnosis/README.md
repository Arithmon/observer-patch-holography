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

The collective calculation constructs `X_j=|0,j><5,j|` on the original
total-two carrier. Its density `Q_j=n_0 n_j` obeys exact native continuity,
and the common occupation transports a core-relative charged response.
A complete first-incident-clock experiment has signed Ramsey response
`1/2048` at a local detector, retaining every timeout and detector output.
Two overlapping local observer instruments have the same complete
postmeasurement record law. The operators have matrix-unit statistics and
a diffusive transfer response; the tested direct sharp-pole dictionary fails.
The exact characteristic polynomial of that ten-mode accepted transfer has
one double root and eight simple roots for both C and E. Hence an unbroken
continuous color `SU(3)` action confined to these linear modes and commuting
with their transfer is necessarily trivial: each invariant eigenspace has
dimension at most two. This is a scoped direct-color-attachment obstruction,
not a restriction on composite, enlarged or differently encoded responses.
Physical fermion/boson and Lorentzian mass identifications remain open.
Preparation and coherent readout are declared apparatus, without another
pulse-level execution claim.
Independent ideal-reference copies give an ordered-Weyl thermal bosonic
fluctuation limit with occupation `8/17`, under an explicit independence
and control contract. Exact finite-copy moments and direct tensor controls
certify that statement without claiming finite operator-norm CCR or an
executed growing-copy experiment. Copy normalization preserves the two-point
poles. The collective field has only two positive modular frequencies;
exact bounds keep every tested transfer-log hierarchy below 71. These
conditional rate assignments cannot reproduce the charged-lepton hierarchy,
and they do not establish that the rates are physical masses.

The [same-carrier exterior test](native_exterior_attachment.py) uses the
native binary two-load space as the exterior square of twelve port
coefficients. Its constant and rank-five bands contain a 15-dimensional
subspace with a `10+5` split under a specified exterior vector-space map.
Exterior orientation signs are an additional representation choice; the
native unordered-pair repair does not inherit them. The protected core
projects to a genuine two-plane in that rank-five band. A natural additive
charge action reproduces
the decoded `10` charges but gives `5` rather than the required conjugate
`bar5`; shifting the constant-line charge cannot fix both parts. The actual
one-round mixed-cap repair also leaks out of this subspace, with exact
projector commutator rank 30; retaining selector timeout gives the same
commutator. Its actual C/E quantum branches give nonzero outside probability
from the maximally mixed rank-15 input. More broadly, the complete C/E work
channel has fixed-algebra
projection ranks `0,1,10,11,20,21,45,46,55,56,65,66`; no 15-dimensional
work-state code remains invariant without active correction or an enlarged
carrier. The exterior vector-space map alone supplies no fermion statistics.

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
