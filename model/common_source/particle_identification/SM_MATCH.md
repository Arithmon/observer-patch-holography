# Existing Standard Model identification and the new measured excitations

The public OPH construction already recovers the internal quantum numbers
of one Standard Model generation under its stated lift. It also supplies
exact Higgs interaction tensors, conditional Spin and spatial-field joins,
and a conditional three-family band. These results can constrain the new
response experiment. They do not yet identify its charged gap as an electron,
quark, Higgs, or vector-boson mass.

This audit uses the public [gauge paper](../../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex),
[particle paper](../../../paper/deriving_the_particle_zoo_from_observer_consistency.tex),
[observer-particle contracts](../../observer_particles/DERIVATION.md),
[selection ledger, rows 13–14](../../../docs/registers/SELECTION_LEDGER.md),
[observation ledger, OL-G1–G9, OL-K1–K7 and OL-M1–M3](../../../docs/registers/OBSERVATION_LEDGER_V3.md),
[premise register](../../../docs/registers/PREMISE_REGISTER_V3.md), and
[particle theorem-gap register](../../../code/particles/THEOREM_GAP_REGISTER.md).
Their distinctions between exact internal representation, declared field
action, conditional numerical coordinate and physical observation are retained.

## The strongest established internal match

Complete twelve-port compact response with endogenous inner transport has
local Lie type `u(1) + su(2) + su(3)` under the gauge paper's A1/A2 response
clauses. This is stronger than a dimension coincidence. The inverse-port
recurrence by itself does not construct that nonabelian current: the paper
also records exact abelian and oversized-current controls. The finite
six-level geometric response is an explicit conditional current realization.

The [representation decoder](../../observer_particles/representation_identification.py)
and its [exact receipt](../../observer_particles/results/representation_identification.json)
lift the response to a five-component carrier with primitive central
character `(-2,-2,-2,3,3)`. Its nonvacuum even exterior component
`M15 = Lambda²(C⁵) + Lambda⁴(C⁵)` has the following joint charge/Casimir
eigenspaces. The decoder recovers these blocks without using a field-name
lookup table; the names below identify the resulting representation types.

| Internal type, all left Weyl convention | `(SU(3), SU(2), Y)` | Dimension | `(C2_color, C2_weak)` | `Q=T3+Y` weights |
| --- | --- | ---: | --- | --- |
| `Q` | `(3,2,1/6)` | 6 | `(4/3,3/4)` | `2/3,-1/3`, each three colors |
| `u^c` | `(bar3,1,-2/3)` | 3 | `(4/3,0)` | `-2/3` |
| `d^c` | `(bar3,1,1/3)` | 3 | `(4/3,0)` | `1/3` |
| `L` | `(1,2,-1/2)` | 2 | `(0,3/4)` | `0,-1` |
| `e^c` | `(1,1,1)` | 1 | `(0,0)` | `1` |

Highest weights distinguish color triplets from antitriplets even though
their quadratic Casimirs coincide. `u^c,d^c,e^c` are the left-handed
charge-conjugated singlets, rather than the right-handed fields with the
same charge signs. Fifteen counts internal Weyl components; it is not a
list of fifteen independent massive Dirac particles. The minimal table
contains a left neutrino and no right-neutrino singlet.

The exterior census cancels the five perturbative gauge/mixed anomaly
forms and contains four weak doublets per generation. It has the diagonal
`Z6` kernel, with exact descent criterion `2t+3d+6Y=0 mod6`. The maximal
faithful image of this representation is therefore the familiar quotient.
The physical current, complete character category and physical global-form
selection remain the separate OL-G2/PR-35/PR-46 attachments. Exterior
parity and this internal lift do not establish physical spacetime chirality
or fermion statistics (OL-G3/G5, PR-47/PR-59).

For comparison with the physical chiral benchmark, see the
[PDG electroweak review, introduction and equations 10.2–10.6](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf).
Its hypercharge convention uses `Q=T3+Y_PDG/2`; the table here uses
`Y=Y_PDG/2`. The comparison imports no measured mass or coupling into a
source calculation.

## Bosons, interactions, families and masses already available

| Existing result | Established content | Remaining attachment relevant to the new experiment |
| --- | --- | --- |
| [Gauge response decomposition](../../observer_particles/results/gauge_particle_identification.json) | Exact current directions with adjoint dimensions `8,3,1`, the internal types of gluon, weak and hypercharge currents | A common propagating vector-field response, physical current normalization, transverse helicities and the relevant phase/poles; internal adjoints are not those poles |
| Same complete 36-coordinate operator domain | Weak Casimirs `0,2,6`; no `3/4` doublet eigenspace | This rejects a doublet inside that particular `End(C6)` domain, not the separate exterior matter construction or all OPH sectors |
| [One-Higgs action](../../observer_particles/results/unified_yukawa_action.json) | One declared `H=(1,2,1/2)` admits `Q u^c H`, `Q d^c H*`, `L e^c H*`; all twelve internal current identities and 648 three-family/component checks hold | Physical scalar/vacuum selection, kinetic normalization and Yukawa coefficients; three declared families have 27 complex charged-fermion coefficients before further constraints |
| [Spin/matter join](../../observer_particles/results/spin_matter_join.json) | A source-executed 30-dimensional `C²_Spin tensor M15` representation, commuting factor actions and uniform Spin central sign `−I` | The source rotation data allow both Weyl boost extensions; the physical chirality, connection, CAR attachment and charged propagator remain separate |
| [Shared spatial CAR field](../../observer_particles/results/shared_spatial_matter.json) and [metric join](../../observer_particles/results/geometric_matter_join.json) | Actual shared spatial operators and conditional scalar/Dirac/Maxwell principal cone; extra modes retained | A common charged pole carrying the decoded internal representation; Dirac/cochain identities alone do not remove unwanted species or identify physical chirality |
| [Rank-three band theorem](../../../code/a5_closure/family_band_attachment_certificate.py) | Under single-complete-band and cost-order premises, `5−sqrt(5)<6<5+sqrt(5)` selects the faithful rank-three response band; tensoring `M15` gives rank 45 | Physical family interpretation and a shared matter action remain open (OL-G4, PR-36/PR-47/PR-54) |
| [Family naturality](../../observer_particles/results/family_naturality.json) | On that actual band, the full unbroken `A5` commutant has complex dimension one; identical invariant left/right family attachments give `Y=yI` and `1:1:1` singular values | A covariant noninvariant background can evade these hypotheses; a concrete native one is computed below |
| [Optional neutrino action](../../observer_particles/results/neutrino_action_extension.json) | Same-module dimension-five `LH LH` channel, with six complex three-family coefficients and a scale | Neither coefficient matrix, scale, physical mechanism nor neutrino masses are selected |

There is no nonzero **physical pole mass** supplied as a source-only result
by these packets. This is also the particle paper's explicit abstract
boundary. The [numeric status](../../../code/particles/RESULTS_STATUS.md),
[closure matrix](../../../code/particles/DERIVATION_CHAIN_CLOSURE_MATRIX.md)
and [postdiction ledger](../../../docs/POSTDICTION_LEDGER.md) retain separate
conditional coordinates and empirical comparisons: running-chart `W/Z`
coordinates about `(80.330,91.119)` GeV, a conditional Higgs/top criticality
coordinate about `(125.77,172.63)` GeV, target-anchored charged-lepton
comparisons, and withheld source-only quark masses. These are not masses
measured by the new interferometer. Their chart, selection, scale, scheme,
transport and physical pole premises cannot be removed by matching a number.
The [conditional electroweak mass-form receipt](../../observer_particles/results/electroweak_mass_form.json)
likewise leaves its vacuum scale and physical pole identification open.

## A native symmetry breaker is already present

[family_attachment.py](family_attachment.py) adds a concrete use of two
existing source results: the [prepared protected-core reference](../../native_modular_reference/local_preparation.py)
and the [native slow response band](../slow_geometry.py). It does not change
either construction or the unbroken-family theorem.

In the ideal binary two-occupation reference, the core `I=(0,1)` has sector
probabilities `Pr(q=0,1,2)=(2/5,5/9,2/45)`. Consequently

```text
E[n_i] = 29/90   for i in I,
         61/450  otherwise.
M_sigma = diag(E[n_i]) = (61/450) I12 + (14/75) P_I.
```

This reference distinguishes the core from the other ports. Full `A5`
invariance, a hypothesis of the family degeneracy theorem, is absent.
The protected core and preparation protocol are existing operational choices;
the resulting anisotropy is computed, rather than a spontaneous selection
of that core from a fully symmetric vacuum.
Let `P3` be the actual native slow-band projector, `B` the existing native
coordinate frame spanning its range, and `G=B^T B`. The source feature

```text
F_sigma = G^(-1) B^T M_sigma B
```

is `G`-self-adjoint and has three distinct positive eigenvalues, exactly

```text
61/450,
41/225 - 7 sqrt(5)/750,
41/225 + 7 sqrt(5)/750.
```

Their approximate values are `0.135556,0.161352,0.203092`. The producer
checks all sixty exact covariance identities
`F_(g sigma)=R_g F_sigma R_g^(-1)`. The reference stabilizer has size two;
the compressed feature stabilizer has size four because compression loses
some information. Its full group twirl is exactly `I3/6`. Thus it is a
covariant source feature, not an invariant endomorphism contradicting Schur's
lemma. The exact source arithmetic uses `Q(sqrt(5))`; float displays are
diagnostics. The exact ideal reference is consumed; finite preparation
errors are not propagated by this feature calculation.

Native and retained family artifacts use different port charts. The
producer constructs and checks an oriented incidence isomorphism, then
checks all 144 transported projector entries. It does not silently identify
ports carrying the same integer label. The recovered chart and source hashes
are exported by `family_attachment.report()`.

This **supplies an actual available family symmetry breaker**. The next map
is concrete: determine how the coherent action couples this endomorphism to
the already verified `Q u^c H`, `Q d^c H*`, and `L e^c H*` interaction lines,
and to their kinetic forms. Assigning the endomorphism as `Y=F_sigma` is an
explicit possible architecture, rather than a consequence of compression.
With common canonical kinetic normalization its relative singular values
would be

```text
1 : (82/61 - 21 sqrt(5)/305) : (82/61 + 21 sqrt(5)/305),
```

all within a factor 1.5. That declaration supplies neither the observed
hierarchies nor three different up/down/lepton textures or mixing; no
empirical fit is performed. Nonlinear or interaction-generated spectral maps
could produce other scales, but they require an actual selected action and
response calculation. The missing object is therefore the mass-generating
coupling/readout of an available native family feature, not the existence of
any possible native anisotropy.

## Attaching the representations to a propagating spinor

The public Spin/matter join and the metric Dirac construction permit the
explicit finite tensor architecture `S tensor M15 tensor F3`: a genuine
declared spatial Spin module, decoded internal matter module and native
family band. Independent factor actions preserve the known charge projectors
and invariant tensors. This uses the existing positive results directly.

Tensoring the new free eight-component field with `M15` is weaker than
identifying a Standard Model chiral field. A spatial operator
`D_spatial tensor I15` gives every internal component the same spatial
dispersion. Appending `I3` adds three degenerate copies. A four-component
Dirac factor with the same internal gauge action on both chiralities is a
vectorlike construction unless a further chiral assignment is made. Neither
tensor identity derives that assignment or the physical currents. A bare
mass inserted across unrelated left Weyl charge sectors is not a substitute
for the verified Higgs interaction lines. Coupling `F_sigma` through those
lines is a useful specified next architecture; its pole representations,
statistics, kinetic normalization and observer readout must then be computed
from that action. Operator identities on a tensor product do not themselves
prove that an observed charged pole transforms in that product.

## What the new excitations actually match

| Construction | Response established | Match to the physical Standard Model |
| --- | --- | --- |
| [Finite interacting SSH ring](../mass_response/DERIVATION.md) | Spinless CAR occupation sectors; identical-copy neutral transfer frequency `2 Delta4`, and a measured ratio `Delta4/Omega≈0.63910260` with a later local density response | A source-clocked finite charged excitation. Its occupation charge is not yet the decoded hypercharge/electric charge; no color, weak doublet, Weyl chirality or bulk pole is identified |
| [Separate free three-axis continuation](../mass_response/dispersion.py) | Eight retained components, positive-band gap `sqrt(3) abs(J)` and an analytic refinement to a Lorentz-form free dispersion | A propagation/spinor-regulator control. Internal gauge quantum numbers, physical chiral assignment and family splitting are not inferred from its Clifford identities |
| [Alternative Hubbard law](../mass_response/HUBBARD_GAP.md) | Uniform spinful chain, exact nonperturbative charge threshold and half-threshold holon/antiholon energies; gapless spin sector | An interaction-generated 1+1-dimensional charge-sector example. A holon has charge but no spin; a spinon is neutral and gapless. Neither is a Standard Model quark/lepton, and the electron response is not an isolated relativistic electron pole |
| Existing conditional `M15` and Spin join | The exact internal types in the first table, with independent finite Spin action | These already match the required representation types. Their attachment to the newly measured charged excitation is not yet executed or selected |

The Hubbard charge theorem remains valuable without supplying observed
particles: its positive threshold vanishes when the interaction is disabled
while hopping and clock remain fixed. It demonstrates a mechanism for
creating a small scale without a bare mass input. It is a different coherent
law from the executed finite SSH ring and is not an interacting 3+1 theorem.

The most specific integration target is to combine the decoded `M15`, the
actual native covariant family feature and the exact one-Higgs tensors in a
single propagating Spin/CAR history. Charges and masses would then be read
from the **same charged correlation functions and clock**, including weak
current selection, kinetic residues and any extra modes. The current results
supply substantial pieces of that target. They do not yet show that the
physical particles or their masses have been identified, and they do not
invalidate the larger OPH construction because the small ring has fewer
quantum numbers.

Run the exact feature and its independent checks from the repository root:

```sh
python -m model.common_source.particle_identification.family_attachment
python -m pytest -q model/common_source/particle_identification/test_family_attachment.py
```
