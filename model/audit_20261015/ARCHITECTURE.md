# Architecture and identification audit — 15 October 2026

The finite constructions have substantive, explicitly scoped execution and
error results. The principal architectural risk is spending additional work
on universal implementation while leaving the physical selection problem
unchanged. The current documentation usually distinguishes those obligations
correctly. This audit found concrete ways to test and simplify that distinction;
it does not invalidate the existing results within their stated families.

The review covered [source selection](../SOURCE_SELECTION.md), the
[construction contract](../CONSTRUCTION.md), the native scalar compiler,
common-source histories, scalar observations, the three reconstruction packets
and the relativistic scalar regulator. Frozen evidence was not changed.

## 1. Exact calibration ambiguity across occupation cutoffs

**Finding: four-level calibration cannot identify an unrestricted oscillator
action without the supplied action-family premise.** This is stronger than
finite-shot uncertainty or an imperfect short-time derivative estimate.

The [coupled calibration](../coupled_reconstruction/source.py) uses four levels
at each of eight sites. Consider adding the same nonnegative onsite operator

\[
 B(N)=\frac{N(N-1)(N-2)(N-3)}{24}
\]

to the Hamiltonian at every site. Its first eight diagonal entries are exactly
`[0,0,0,0,1,5,15,35]`. Its compression to the entire calibration code is zero.
Every calibration preparation, evolution interval, detector and reference
therefore has exactly the same law, even with perfect probabilities and
arbitrarily many calibration times. This is an exact statement about the
implemented compressed calibration, not an assertion that physical
unrestricted evolution cannot leave the four-level subspace.

The two Hamiltonians differ once the read site retains eight levels. Applying
the existing two-read instrument, feedback and model intervals `1/128` to the
`8 × 2^7` reference calculation gives:

| Pointer sign | TV between the original and extended two-read laws |
| --- | ---: |
| −1 | `0.0000635812943` |
| +1 | `0.0000635407035` |

These are independently calculated sparse-Hamiltonian diagnostics, not new
native source executions or outward numerical certificates. The calibration
identity itself is exact. The alternative is outside the declared canonical
quartic family; the family assumption remains legitimate, but the records do
not infer it. The added term is local, positive and homogeneous, so these
generic structural properties do not exclude it.

**Consequence:** predicting the larger cutoff is a test of the chosen
extension, rather than a unique reconstruction of that extension. More
precision on the current probes cannot resolve the ambiguity. The inexpensive
next diagnostic is a preparation/phase read involving level four, followed by
an explicit wider family containing diagonal and sextic controls. A bounded
energy/observable target can also make unidentifiable higher-level terms
irrelevant, provided their effect is independently bounded.

The exact identity and numerical contrast are retained in
[architecture_counterexamples.py](architecture_counterexamples.py) and
[architecture_results.json](architecture_results.json).

## 2. The cube graph does not determine a three-dimensional metric

**Finding: the recovered rank-three Gram uses an extra geometric convention,
not only the recovered coupling graph or unit edge lengths.**

Take the three planar unit vectors

\[
 v_1=(1,0),\quad v_2=(3/5,4/5),\quad v_3=(-5/13,12/13).
\]

The eight points `b1 v1 + b2 v2 + b3 v3`, with each `bj` zero or one, are
distinct. Connect points whose bit labels differ in one coordinate. This is
the same cube graph, with all twelve edges having exactly unit Euclidean
length, embedded in a plane. The graph Hamiltonian, preparations and detector
records in the reconstruction packets are unchanged.

[cube_geometry](../observer_reconstruction/recovery.py) instead interprets
graph distance as **squared** Euclidean distance. This is correct for an
orthogonal unit cube and gives rank three. It fails for the planar embedding;
the retained rational check reports a maximum discrepancy of `194/65` between
the two distance assignments. The existing rank-three result is correct under
its declared cube interpretation. It is not a measurement that excludes the
planar realization, an internal-channel graph, or an extra ambient dimension.

**Consequence:** a new spatial-dimension claim needs observations constraining
the metric/translation meaning beyond the graph. Suitable targets include
independently calibrated distance-sensitive interventions and a growing
family with a common observable definition, comparing volume growth and
low-energy propagation. A single finite adjacency matrix or its spectral
multiplicities cannot supply that inference by itself. The same issue must
be resolved before identifying native rank-three response with field space.

## 3. Simplify the fixed-band spatial calculation

The [seven-mode scalar action](../scalar_observations/lattice.py) already
supports its exact spectral gradient using `grid=None`. All five retained
grids have the same seven oscillator modes. The quartic integral is exact
from grid size five onward; only the common first-harmonic gradient
coefficient `kappa_n = (sin(pi/n)/(pi/n))^2` changes. The audit verifies that
all onsite and pair quartic coefficients are identical across those grids.

For future fixed-band predictions, use the exact-gradient endpoint directly,
retaining the existing grid runs as discretization controls. Increasing the
grid beyond 65 adds neither independent modes nor an ultraviolet test. The
new work with physical value is controlled **mode enlargement**, including
the induced parameter changes and omitted-mode response for the same detector.
Occupation, spatial volume, spacing and momentum bandwidth are four different
limits and should remain separate receipt fields.

## 4. Attach universal implementation results once

The native compiler checks exact source generators and builds proper-code
one- and two-register controls. Given coherent code transfer, blanks,
arbitrary admitted control angles and addressed gates, this has the structure
covered by the standard universal-gate result: arbitrary one-qubit operations
and a controlled NOT suffice for finite qubit unitaries. This is a statement
about the ideal admitted service family; it does not claim every dimension
has been compiled by today's bounded routines.
[Barenco et al., *Elementary gates for quantum computation*](https://arxiv.org/abs/quant-ph/9503016)
provides the reusable theorem.

For a finite instrument with Kraus operators `K[r,a]`, the retained map has
isometry `V psi = sum_(r,a) |r,a> K[r,a] psi` when
`sum K[r,a]^* K[r,a] = I`. Retaining the record, tracing the declared
environment and extending the isometry to a unitary supplies the standard
implementation route. The relevant existence theorem is
[Stinespring's dilation theorem](https://community.ams.org/journals/proc/1955-006-02/S0002-9939-1955-0069403-4/S0002-9939-1955-0069403-4.pdf).
Its map here requires the actual complete instrument, including rejected
outputs and arbitrary references; matching one input's probabilities is
insufficient. The theorem supplies an ideal dilation, not its physical services.

**Consequence:** once a law and interface are supplied, native compilation is
strong implementation evidence but normally not a discriminator of that law.
An alternative supplied finite Hamiltonian or channel can be implemented by
the same admitted control family. The experimental postdictions therefore
test the calibrated apparatus laws that were compiled. To test an OPH-specific
restriction, exhibit a law or intervention excluded by a source-derived
constraint while allowed by the comparison family, and compare that difference
with data. If no difference exists within the admitted interface, the
implementations belong to the same operational class for that test.

An economical workflow is to maintain a tested semantic channel evaluator,
the verified source compiler and their explicit error contract. Explore
candidate laws with the semantic evaluator; execute and hash source words
for locked predictions, changed primitives and compiler regressions. Repeating
millions of already verified gates for every exploratory variant adds cost
without independently establishing a new physical law. Cached contractions
already move the implementation in this direction.

## 5. Treat the continuing observer as a quantum instrument with memory

[Chiribella, D'Ariano and Perinotti's quantum-network framework](https://arxiv.org/abs/0904.4483)
already gives a finite causal network a positive comb representation and
realization by a sequence of memory channels. Here the slots are declared
preparations/interventions, the memory is the surviving work register plus
committed records, and the final output retains the quantum state and all
classical branches. Those objects are explicit in the common-source protocol.

This supplies a useful reusable attachment: compare implementations as
multitime processes under the **same allowed testers**, not only their final
scalar means. Complete operational equality is stronger than equality on the
few retained probes; it must be proved or tomographically controlled. No
additional observer-specific operation is required once the process/interface
is fixed. A controller reproducing that complete process is an operational
implementation of it. Memory bounds or restrictions on available operations
can still discriminate particular controller classes.

For a finite set of classical record experiments, a classical generator can
also reproduce the joint probabilities by drawing successive conditional
distributions, retaining the entire history if needed. This does not establish
efficient classical simulation, preserve arbitrary reference outputs or
reproduce arbitrary future interventions. It demonstrates why record fit alone
does not establish a uniquely quantum or observer-specific mechanism.

A useful next empirical target is a causal-break test: after different past
interventions, prepare the same present system state and compare future laws.
Persistent differences detect operational memory in the chosen system access.
[Pollock et al., *Operational Markov condition for quantum processes*](https://arxiv.org/abs/1801.09811)
provides the theorem and clarifies why matching a one-time channel or fitting
a Lindblad curve does not establish the full Markov property. For OPH, the
initial system/environment split and allowed causal breaks must be specified.
This offers a more direct test of memory/feedback than additional passive
trajectories, although it requires suitable experimental data.

## 6. Apply identification theorems with the actual control menu

[Burgarth and Yuasa, *Quantum System Identification*](https://arxiv.org/abs/1104.0583)
provides an input/output equivalence framework. Its stronger identification
statements require their controllability/observability and model hypotheses;
they do not follow from the two times and few preparation/detector families
currently used. The concrete interface map is the opaque channels, trusted
number preparations, known detector phases, calibrated clock and allowed
action family. Basis/measurement gauge must be fixed by trusted calibration
or retained as an equivalence rather than counted as new observable physics.

The most economical next attachment is a rank/nullspace calculation of the
response map for an enlarged finite action family, followed by experiments
chosen to detect its null directions. The exact invisible term in section 1
already proves one null direction when the action extension is unrestricted.
For accessible directions, optimize times and readouts against finite-shot
uncertainty. The existing derivative estimator's probability error is divided
by its tiny model interval; more symbolic precision cannot remove that
sampling penalty. Longer-time likelihood fitting with independent held-out
interventions is worth comparing against the short-time estimator, while
retaining its honest family and uncertainty assumptions.

## 7. Do not repeat reconstruction where the missing hypothesis is selection

The current GNS and fixed-Gibbs modular attachments are appropriate finite
mathematics. They are already documented as representations of supplied
objects. Choosing `sigma = exp(-beta H)/Z` makes its modular flow agree with
Hamiltonian evolution at `t=-beta s`; the identity cannot identify an unknown
`H` independently when the reference was constructed from that same `H`.
The actionable missing attachment is a source-selected reference/access pair
and an independently readable clock, not a further check of this identity.

For a genuinely relativistic reconstruction, the
[corrected Osterwalder–Schrader theorem](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-42/issue-3/Axioms-for-Euclidean-Greens-functions-II-with-an-Appendix-by/cmp/1103899050.pdf)
is the relevant existing result once a common Euclidean correlation hierarchy
has the required growth/regularity, symmetry, Euclidean covariance and
reflection positivity (with clustering for the corresponding vacuum property).
A finite Hamiltonian or GNS representation alone does not provide that
hierarchy. The current finite-band and finite-cell results do not establish
these hypotheses, especially the interacting spatial limit and covariance.
An attachment must carry the actual observable maps and limit, not merely
identify another positive matrix.

The scalar route also has a known restriction already cited by the regulator:
[Aizenman–Duminil-Copin](https://arxiv.org/abs/1912.07973) establishes Gaussian
scaling limits for the stated critical four-dimensional Ising/phi4 families.
This is not a blanket no-go for every OPH construction or a finite-cutoff
effective theory. It does mean that simply increasing the ordinary positive
scalar regulator's cutoff cannot be presented as an established route to a
nontrivial four-dimensional ultraviolet theory. Parameter matching and an
explicit effective validity domain are the useful immediate targets.

## Priority after this audit

1. Keep the empirical apparatus comparison and native source-selection tests
   distinct. Quantify what new constraint, if any, the latter imposes on the
   former before expanding the physical claim.
2. Retain the exact calibration-extension and planar-embedding controls. Use
   them to design the next measurement, rather than add precision to probes
   they defeat identically.
3. Use the exact seven-mode gradient endpoint for new fixed-band work. Spend
   refinement effort on omitted modes, matched parameters and a shared
   observation domain.
4. Reuse universal dilation/compiler and quantum-comb results as explicit
   interface contracts. Concentrate new proofs on the unproved attachment:
   a constrained source law, accessible clock or physically meaningful metric.

Reproduce the new audit controls from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.audit_20261015.architecture_counterexamples
```

No theorem above supplies empirical validity merely because its mathematical
implementation succeeds. Equally, the counterexamples do not undo the
completed finite predictions under the premises they were designed to test.
