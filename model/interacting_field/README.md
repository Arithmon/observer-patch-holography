# A source-compiled interacting Weyl lattice

This CP-06 / [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142)
advance executes a local, non-Gaussian fermion interaction over 27, 125 and
343 spatial sites, retains complete continued detector records, and bounds
their comparison with the infinite lattice. It builds on the existing
native CAR compiler while treating draft [#1151](https://github.com/FloatingPragma/observer-patch-holography/pull/1151)
as available for development. The recovered clock does not supply this
field action or its physical calibration.

The [contract](campaign.json) states the supplied cubic geometry, Weyl
kinetic term, density interaction, preparation and coherent services.
The [derivation](DERIVATION.md) proves the exact-sector, interaction and
history bounds and maps the remaining premises. No integrated gate is
accepted.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.interacting_field --check model/interacting_field/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/interacting_field
```

The action is

```
H = dGamma(h) + U sum_x n_(x,0)n_(x,1),
h(k) = -c sigma.sin(k),    c=1/8, U=1.
```

The exact two-fermion sector contains **all** 54, 250 and 686 CAR modes of
the three cubes. Its dimensions are 1,431, 31,125 and 234,955. This is an
invariant sector, not a Gaussian approximation or a spatial mode cutoff.
The coordinate time is supplied Hamiltonian model time. The interaction
and kinetic evolution share one action throughout each protocol.

A compiled local preparation, two evolution intervals, two complete Ramsey
reads and record-controlled feedback produce four retained outcome branches
with their disturbed quantum outputs. There is no reset or postselection.
All local hopping, contact, preparation and detector matrices are obtained
by contracting native words on their full local Fock spaces. Exact CAR
routing carries the FSWAP sign and its full requested resource cost; the
large run does not replay every repeated routing pulse numerically. Tiny
forbidden floating entries are compressed in the global sector replay; their
full local residuals are reported. Exact invariance refers to the ideal
completed words. The full quantum vectors are regenerated during replay;
the JSON retains probabilities and their complete-output comparisons.

For 27, 125 and 343 sites, respectively, the interaction-deletion history
TV distances are 0.0326546, 0.0326541 and 0.0326541. Source/sparse-reference
history TV distances are about 2e-9 at each size. The ideal full-output
bounds to the infinite lattice are 0.437504, 0.0478550 and 0.00349258.

The exact rational bounds concern ideal completed compilation. Floating
local/source residuals, their repeated-call estimates and physical hardware
errors are separate. Numerical comparisons do not become rigorous physical
error bars. The [receipt](results.json) retains both kinds of evidence.
Replay compares exact counts, bounds and source hashes exactly and floating
observables within 2e-10. The accumulated evolution residuals sum local
4-by-4 floating matrix errors over millions of calls; backend roundoff can
therefore exceed that absolute tolerance. For those two diagnostics alone,
replay checks a cap of 128 binary64 eps per counted hopping, contact or
routing call and the defining `expm1` telescoping identity. This numerical
replay cap is not an outward certificate, and changes no observable or
ideal-model acceptance bound.

Two positive interaction results strengthen the observation contract:

- An initially Gaussian Slater state acquires `2-Tr(gamma^2)=sin^2(1/4)`
  and connected four-point **-0.01530218**. Sixteen complete source detector
  settings reconstruct gamma, and a separate complete paired read measures
  the four-point. Seventeen prepared copies and every detector complement
  are accounted. The U=0 control remains Gaussian.
- With kinetic and contact terms both running, a separate short complete
  read has interaction-minus-free probability at least
  `141294691/21579694080 > 0.00654` after ideal source splitting error.
  This bound is uniform in volume and remains valid on the infinite lattice.

The source-specific clock boundary is also sharper. For every positive
product-weight reference, including #1151's record-inferred table,
`-log(pi)=-sum log(theta_i)n_i+constant` on its occupation CAR dictionary.
Its modular unitary therefore preserves Slater states. It cannot generate
the observed non-Gaussian contact response through that dictionary. Both
target and actual inferred-weight clock controls are executed. A new joint
state/algebra map or additional field controls remain necessary.

This is a nontrivial interacting **three-spatial-dimensional lattice field
with continuous model time**. The kinetic symbol has eight Weyl nodes,
four of each chirality, all retained. The one-particle low-momentum symbol
has a quantitative Weyl comparison, but contact scattering can leave its
momentum patch. There is no interacting chiral continuum claim or taste
removal. The vacuum, physical particle identity, gauge/Higgs sector,
calibrated stress/energy and gravity joins remain open. The result supplies
executable interacting-field machinery and precise obstructions for the
completion program; it does not establish the Standard Model.
