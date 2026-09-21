# Coupled fermion and gauge-field transport

This interface extends the pinned semiclassical current history with electric
link drift, Wilson plaquette kicks and a further matter-current readout. The
64 prepared golden sites carry bounded software patches with local orbital
state, link ports, retained readback, electric feedback and public records.
The verifier checks actual consumed writers, versions and values.

From the repository root:

```sh
python3 code/sm_fermion_current/build_coupled.py --write
python3 code/sm_fermion_current/verify_coupled.py
python3 -m pytest -q code/sm_fermion_current/test_coupled.py
```

The builder writes `coupled_receipt.json`. It preserves the parent execution
inside each run's `prefix`; it does not rewrite `current_receipt.json`.
The continuation performs all 144 electric drifts and 108 plaquette kicks,
then applies the matter factor on edge 14 and records its result. That final
factor is a declared partial second matter sweep. It consumes the evolved
link phase. The baseline, phase intervention and gauge copy follow the same
factor word. `drift_disabled` is a separate counterfactual control whose link
values are retained unchanged in place of electric drift.

The action dictionary distinguishes the fixed kinetic-shape parameter
`kinetic_lambda` from `electric_duration`. Their values are both `1/2` in
this fixture. The rational electric phase is an exact drift of the declared
convex nonquadratic kinetic law over that duration. It is not the drift of
the quadratic Maxwell kinetic law. Matter factors use the variational CAR
midpoint rule; Wilson kicks use the same declared gauge action. The composed
word preserves mean Gauss and charge continuity, without asserting exact
summed-energy conservation or a continuous-trajectory error estimate.

The action has three declared families and 27 arbitrary complex Yukawa
coefficients. Family 0 is occupied; the other two are explicit Fock vacua.
Internal covariance is proportional to the identity, anomalous pair
expectations vanish, and Higgs coordinates and momenta are zero. The omitted
nonabelian and Higgs/Yukawa forces vanish in expectation under this supplied
restriction. Their operator fluctuations are not set to zero. The gauge
action uses the direct-product cover; the charge-one plaquette is not
asserted to descend to its diagonal central quotient.

Checkpoint `instantaneous_current` differs from the return factor's
`return_midpoint_current`. Exact rational local energy terms are retained;
their sum intervals use outward rounding on the declared `2^40` dyadic grid.
Electric kinetic values are represented by their rational arguments and the
explicit action function. No binary64 value defines a scientific result.

The receipt is capped at 20 MB and each retained run at 4.5 MB. It provides a
finite semiclassical action and authenticated software transport. Source
selection, an operational physical clock, quantum operator Gauss, a chiral
continuum theory and laboratory current identification are outside its scope.
