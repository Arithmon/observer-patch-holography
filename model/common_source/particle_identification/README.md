# Identifying the measured excitations

This audit identifies the charged operator response of the
[continuing repair/clock experiment](../mass_response/README.md), classifies
its separate free continuum candidate, and connects the results to the
existing OPH Standard Model representation theorems. It also constructs a
native family symmetry breaker. These are particle-content and attachment
tests; no observed Standard Model species or mass is claimed as identified.

## What the finite observer actually detects

[response.py](response.py) computes the complete finite Lehmann response of
the same four-mode interacting law. It retains all matrix residues, including
dark channels, and verifies CAR zeroth/first moments, density covariance and
charge sum rules. Its lowest charged addition and removal costs are both
`1.028595953717032`. This is the gap measured in the prior clock experiment.
The local charged response is nonzero there. A globally charge-neutral
transfer between the two cells has pair frequency `2.057191907434065` and
nonzero spectral weight `0.061838316443786` at the chosen sites.

Thus the earlier clocked gap really does couple to an encoded fermionic
field response. The neutral density threshold is instead
`1.305496913774670`. It cannot be substituted for that charged gap.
The exact spectral decomposition is finite-dimensional; displayed numerical
energies and residues are floating diagnostics, not outward certificates.

The [retained replay](results.json) goes further: it reexecutes the source
history and appends an actual even CAR transfer read to **every** surviving
branch, keeping all earlier records, clocks and purification references.
There is still one work preparation and no reset. It retains 64 normal
six-bit records and 12 terminal aborts. The companion
[quantum outputs](results.npz) contain their complete unnormalized
work/reference vectors. An independent continuous instrument supplies the
comparison; the ideal splitting budget contracts through the appended read.
The final transfer mean after the earlier instruments is a continuing field
observation, rather than a second independent pole fit.

The charge is global integer occupation. If it is promoted to a local gauge
charge, a transfer requires its corresponding link or dressing. The present
apparatus supplies neither electromagnetic identification nor spacetime
locality for that inter-cell transfer. The finite ring is spinless. Its
nondegenerate elementary charged levels cannot acquire a color-triplet
representation merely through a basis change. These restrictions apply to
this finite candidate, not to all OPH constructions.

Read the [response derivation](RESPONSE.md) for spectral and source-word
conventions.

The [inverse-response diagnostic](kinetic_response.py) computes the inverse
finite CAR matrix and its frequency derivative from those same complete
Lehmann residues. Independent tests compare the free inverse with the actual
one-particle action, differentiate the interacting response numerically, and
change operator coordinates nonunitarily. Generalized free poles remain
fixed only when the kinetic and constant terms transform together; ordinary
eigenvalues of the transformed constant matrix do not. The interacting
frequency dependence is retained. This is a check using a supplied full
response, not a reconstruction of that matrix from observer records, a
spatial-momentum derivative or a physical positive kinetic metric.

The transfer half-gap is now admitted only at the source's symmetric
chemical-potential convention. A nonzero common chemical potential leaves
the neutral pair frequency unchanged while splitting addition/removal costs;
the interface then reports those costs and withholds the single-gap readout.

## The free continuum candidate contains two Dirac copies

[clifford.py](clifford.py) constructs the actual continuum spin generators
of the separate eight-component three-axis extension. Their Casimir is
`3/4 I8`. The four Clifford generators act as two equivalent four-component
Dirac representations. Each energy sign has two spin states per copy.

The common rest-energy coefficient is `sqrt(3)*log(125/64)`, with ratio
`0.7204303414` to the programmed clock frequency. These are coordinate/clock
units, not measured particle masses. A particle/antiparticle interpretation
requires the CAR field and filled negative-band reference. Matrix signs
alone are not a vacuum or a spin-statistics derivation.

The continuum internal commutant is `M2(C)`. It cannot carry a nontrivial
commuting `SU(3)` color action on that two-dimensional multiplicity. The
full finite regulator has only a scalar universal momentum-independent
internal commutant. A free gapped continuation has no massless photon pole;
neutral bilinear thresholds do not establish a Higgs or vector bound state.
The [Clifford derivation](CLIFFORD.md) proves these statements constructively.

## Existing Standard Model results and an available family feature

The prior OPH representation decoder already recovers
`Q,u^c,d^c,L,e^c`, with their color, weak and hypercharge representations,
under its declared lift. Those internal types remain valid. Their physical
chiral Spin/CAR action must be attached to a propagating response to identify
observed particles. Simply tensoring identity spatial evolution with their
labels does not supply that attachment.

The new [family feature](family_attachment.py) combines the existing native
reference and rank-three slow band. In its actual source chart,

```text
F = (B^T B)^(-1) B^T diag(E_sigma[n_i]) B
eigenvalues(F) = 61/450, 41/225 - 7 sqrt(5)/750, 41/225 + 7 sqrt(5)/750.
```

All three eigenvalues are distinct and positive. The producer verifies all
60 exact covariance identities and explicitly transports the native and
retained-family port charts. The full group twirl is `I3/6`. Therefore the
prepared core supplies a real covariant symmetry breaker: the earlier
unbroken-family Schur conclusion does not apply to this state. That theorem
remains valid under its original hypotheses.

This closes the availability of a nonscalar native family operator. It does
not make that operator a physical Yukawa interaction. The direct declaration
`Y=F` with common kinetic normalization gives ratios within a factor 1.5,
and supplies no observed hierarchy or distinct quark/lepton textures.
The specific remaining map is how the coherent action couples this feature
to the already verified Higgs interaction lines and kinetic forms, and what
charged propagators result. The exact ideal reference, protected-core choice
and physical family attachment retain their stated contracts.

[SM_MATCH.md](SM_MATCH.md) maps the existing gauge, matter, Spin, family and
field theorems to these new responses, including public ledgers and primary
benchmark references. The small ring, its free extension and the separate
Hubbard law are distinct candidates throughout the audit.

## Replay

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/common_source/particle_identification
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -u -m model.common_source.particle_identification.run --output /tmp/particle-audit.json
.venv/bin/python -m model.common_source.particle_identification.run --verify model/common_source/particle_identification/results.json
```

Replay creates a fresh JSON receipt and a companion NPZ. It refuses to
overwrite either. Verification checks source custody and retained quantum
weights, without rerunning the dynamics. Public representation producers
are read and checked; their frozen artifacts are never regenerated here.
All source, preparation, clock, numerical and physical boundaries of the
underlying [mass-response experiment](../mass_response/README.md) remain in
force. This package is unreleased and changes no r2053 artifact.
