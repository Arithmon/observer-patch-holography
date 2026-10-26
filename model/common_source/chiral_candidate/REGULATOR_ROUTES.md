# Chiral regulator routes and the current construction's domain

The protected Weyl--Higgs candidate supplies an all-left symbol and a finite
spectral menu. The separate [overlap calculation](../overlap_regulator.py)
now builds a finite Euclidean, gauge-covariant chiral kinetic operator on the
same retained matter representation. It is not a quantum gauged Standard Model
or a continuing Lorentzian field. The naive periodic sine derivative is retained
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

The [retained calculation](../overlap_regulator_results.json) assembles all
four Wilson directions on a periodic $2^4$ lattice, with the actual
$\Lambda^2\mathbb C^5\oplus\Lambda^4\mathbb C^5$ link matrices. The
$960$-dimensional one-family Dirac kernel has one free Dirac zero at the
origin; the other 15 corners have $aD=2$. At a nonconstant nonabelian
$SU(3)\times SU(2)\times U(1)$ background, the full operator and modified
left projector obey the Ginsparg--Wilson relation, the Weyl-map identity and
site-dependent gauge covariance within numerical admission. This implements
a kinetic projector for three independent copies of the declared generation;
it does not select their family couplings.

The declared small-link domain has represented plaquette norm below $1/40$.
The [Hernandez--Jansen--Lüscher overlap gap bound](https://arxiv.org/html/hep-lat/9808010v2)
then gives $1/4\le H^2\le64$ on that domain.
Its degree-$m$ Legendre approximation has range at most $2m+1$ lattice links
and ideal operator tail at most $16(15/17)^{m+1}$, independent of volume;
floating evaluation error is separate. The free Euclidean symbol satisfies
$\|D_a(p)-i\gamma\cdot p\|\le a|p|^2$ when $a|p|\le1/4$.
Two complete native source queries execute the free four-spinor overlap
kernel with both success and failure amplitudes retained. They are kernel
queries, not execution of the three-family interacting Fock field. The
spacing, gauge background, source service accuracy and whole-field state are
still supplied. Global chiral measure phases, reflection positivity,
Lorentzian continuation, interacting gauge/Higgs dynamics and physical
mass/clock identification remain separate tests.

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
