# Chiral regulator routes and the current construction's domain

The current executable candidate supplies an all-left Weyl symbol and a
finite spectral menu. It is not a local nonperturbative lattice realization
of the gauged Standard Model. The naive periodic sine derivative is retained
as a negative control: it has eight nodes with balanced handedness. A host
cutoff or a projection that drops these states is not a decoupling theorem.

The [source action](ACTION.md) has the exact internal module
`M15=Lambda²(C5)+Lambda⁴(C5)=10+bar5`, with integer hypercharge `q=6Y`.
A potential regulator must transport this representation and its Higgs
vertices, rather than substitute a vectorlike field with the same dimension.
It must also preserve the observer clock and actual charged response used by
the executable construction. These are compatibility requirements for this
project, independent of any particular regulator update.

## Ginsparg–Wilson and overlap

[Lüscher (1998)](https://arxiv.org/abs/hep-lat/9802011) establishes an exact
lattice chiral symmetry from the Ginsparg–Wilson relation. Its chiral
transformation differs from the naive anticommuting one; the usual doubling
argument therefore does not directly apply. This is a constructive route
for testing the kinetic projector, not a selection principle for OPH
couplings or a proof of a complete nonabelian fermion measure.

For a suitable overlap operator `D`, the finite identities to check are

```text
gamma5 D + D gamma5 = a D gamma5 D
hat_gamma5 = gamma5 (I-a D),  hat_gamma5^2=I
hat_P_left=(I-hat_gamma5)/2.
```

Gauge-dependent chiral projectors also require consistent measure phases.
Checking these algebraic identities at one gauge configuration does not
establish locality, global consistency, or a Lorentzian Hamiltonian.

## A recent measure construction matches the internal module

[Craig, 6 October 2026, section 10](https://arxiv.org/html/2610.09093v1#S10)
proposes a nonperturbative overlap measure on admissible gauge fields. Its
Standard Model restriction uses exactly `10+bar5` and
`diag(z^(-2) g3,z^3 g2)`, the same representation pattern as the OPH lift.
The auxiliary `SU(5)` construction does not introduce extra Standard Model
gauge fields. The paper leaves reflection positivity, an interacting
continuum definition and some phase choices open. This is a relevant recent
preprint to audit; its theorem and supplementary formalization have not
been independently replayed here. It does not select masses or supply the
continuing OPH Hamiltonian.

The representation match supplies a concrete comparison map. An OPH reuse
would still need actual gauge links in the stated admissible domain, its
global sector/measure data and a justified route from the Euclidean measure
to the observer's continuing quantum dynamics. Existing finite current
matrices alone do not meet those conditions.

## Hamiltonian alternatives

[Thorngren, Preskill and Fidkowski (2026)](https://arxiv.org/html/2601.04304v2)
give symmetry disentanglers and an exactly solvable 1+1-dimensional chiral
example. Their 3+1-dimensional hypercharge discussion is a proposed route,
not a completed local Standard Model Hamiltonian. It is relevant because
finite local unitary circuits could be mapped to source operations, provided
the required interface and low-energy modes are actually established.

[Sen (October 2026)](https://arxiv.org/abs/2610.11035) proposes a Hamiltonian
framework with a disk boundary and supplies an Abelian bulk-field extension.
It is another candidate architecture for comparison, not an implemented
nonabelian OPH model. A bulk or mirror sector must remain in the accounting
until its decoupling and residual effects are controlled.

The current spectral candidate can test charge assignments, Higgs pairing,
source execution, spin/helicity and mass readout while these regulator
questions remain explicit. None of the regulator routes would make its
present trial Yukawa spectrum agree with measured masses automatically.
