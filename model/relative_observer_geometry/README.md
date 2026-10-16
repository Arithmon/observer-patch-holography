# Relative observer fields, records and clocks

This packet closes a specific observer-access gap: on a **known diagonal native
preparation**, absolute occupation fields can be read from retained signed
currents. They need not remain simulator-only variables. Two regional observers
can reconstruct compatible parts of the same field history, continue from
their own checkpoints, and compare local event clocks on shared records.

The result is finite and conditional on declared access. It does not identify
the fields as physical matter or the local counters as physical time.

## What executes

The carrier is the existing binary total-two sector: 66 configurations of two
occupied ports among twelve. The native endpoint rule acts on an oriented
admitted seam `(a,b)`, commits its fair coin and signed current `j`, and changes
the loads by

\[
n_i(k)-n_i(k-1)=j_k(\mathbf1_{i=b_k}-\mathbf1_{i=a_k}).
\]

The rule, carrier, preparation service, record access and scheduler remain
premises. The packet first checks this equality on **all 3,960** native
configuration/seam/coin branches. Summing these one-step identities gives an
exact reconstruction for every finite prefix, without importing a field action:

\[
n_i(N)=n_i(0)+\sum_{k\le N}j_k(\mathbf1_{i=b_k}-\mathbf1_{i=a_k}).
\]

Each `LocalCounter` receives only its initial local loads and every incident
seam record. Boundary seams are essential, including seams whose other endpoint
is outside the region. Nonincident records are absent; incident zero-current
records remain. Common event identifiers authenticate correspondence and give a
source-order cut. They carry no duration, physical simultaneity or global clock.
Authenticity and completeness of the incident record service are premises:
omitting a boundary current can be locally undetectable.

The retained finite example supplies the seam tape
`(0,1), (0,1), (0,5), (0,5), (0,1)`. It enumerates all **2,112** histories from
66 equally weighted, **disclosed** initial configurations and 32 fair-coin words.
There is no successful-outcome selection. These are exact laws on a declared
tape, not the stationary uniform-seam experiment. Regional observers access
ports `0..6` and `5..11`; their overlap agrees on every history and their union
covers all twelve ports. At the same retained prefix their values therefore
glue uniquely. Local checkpoint restoration and continuation agree for every
history. No hidden load vector occurs in a local checkpoint.

Separately, the primitive CQ runtime executes the first three operations from
the disclosed preparation with ports 0 and 11 occupied. Each operation has six
coin/current outcomes. All **216** outcome words and full terminal density
blocks remain, including **208 zero-probability words**. Exact enumeration
checks the eight positive words, and runtime continuation checks the complete
serialized checkpoint state. The record-only `n_0` expectation equals the separately
audited terminal occupation expectation. Zero-probability impossible words are
kept with an arbitrary zero public-field value and are not declared physically
consistent trajectories. Floating runtime comparisons are diagnostics at
`2e-11`, not directed roundoff certificates; the diagonal enumeration and
continuity identities are exact rational/integer computations.

## Local fields and common state

For a known basis preparation, the current prefix determines each local field
value and its products. Averaging those record functions gives the same
occupation correlators as the diagonal source law. A randomly chosen initial
configuration is usable in this way only when its local preparation label is
disclosed. An unknown initial law instead permits conditional estimates; it
does not generally reveal the realized occupation. This packet does not
implement a Bayesian filter or extend the pathwise assertion to arbitrary
coherent preparations.

The existing effective-field construction can predict **current increments**
without this initial disclosure. Absolute occupation readback and increment
readback thus have different access requirements. No extra direct occupation
measurement is silently granted here.

At the classical record level, a local view is a deterministic map `v(h)` from
the retained common history. Its state is the pushforward of the common law.
Pullback gives the explicit finite GNS isometry

\[
(Jf)(h)=f(v(h)),\qquad
\sum_h p(h)|f(v(h))|^2=\sum_v p(v)|f(v)|^2.
\]

The complete finite enumeration records the local fibers and checks their
exact masses. Distinct indicator fibers are disjoint, so their off-diagonal
Gram entries vanish. This is a commutative record-algebra inclusion and a
compatible expectation, not a constructed quantum local net. In particular,
state marginalization does not supply the homomorphic algebra retractions in
[`FiniteCausalObserverNet`](../../Lean/QFT/FiniteCausalObserverNet.lean).
[`primitive_observer.algebra`](../primitive_observer/algebra.py) supplies the
actual CQ public-record expectation used for the runtime readback check.

Pathwise value gluing requires actual pointwise coverage and a shared retained
cut. It does **not** reconstruct a joint probability law from regional
marginals. The negative control uses the laws with equally likely occupied
pairs `(0,10)/(1,11)` and `(0,11)/(1,10)`: both regional marginal laws on `0..5`
and `6..11` agree, while `E[n_0 n_10]` is respectively `1/2` and `0`. This is
the kind of failure ruled out by an explicit descent receipt, not by isotony
alone; see [`ObserverNetDescent`](../../Lean/QFT/ObserverNetDescent.lean).

## A successful and a failed relative clock comparison

A local incident-event counter increments on each nonzero movement incident
to that observer. On the actual coin word `0,1,0,1,0`, shared seam `(0,1)`
movements occur at record identifiers 1, 2 and 5. The local `(port 0, port 1)`
tick readings are

| Shared event identifier | Port 0 ticks | Port 1 ticks |
|---|---:|---:|
| 1 | 1 | 1 |
| 2 | 2 | 2 |
| 5 | 5 | 3 |

The first two shared events fix the affine comparison `v=u`. The predesignated
third event has residual **−2** because port 0 also recorded the excursion
through port 5. Its probability is exactly **1/32 given the disclosed initial
preparation**. Local event counters need not be affine on every shared-event
domain; this does not preclude other valid comparisons or clock constructions.

A positive control reads the shared-seam movement count `u` and a second
declared display `v=2u+5`. The same three shared events give `(1,7)`, `(2,9)` and
`(3,11)` with zero third-event residual. The offset and scale are a supplied
display convention, not independently recovered physical clock calibration.

These computations use the exact two-anchor/third-event algebra isolated in
[`Time/ClockComparison`](../../Lean/Time/ClockComparison.lean). That Lean
interface additionally assumes timelike histories; no such geometric
assumption is instantiated here. The relation of these native local tick
records to observer-relative modular flow remains to be established. A
preferred universal physical clock is not required.

## Geometry connection and remaining boundary

The existing record-only native response calibration is executed again. Its
recovered unweighted graph gives exactly the **same thirty seam interfaces**
used by the field counters; its exact slow response Gram has rank three.
This reuses [`observer_reconstruction/native_geometry.py`](../observer_reconstruction/native_geometry.py),
including its supplied uniform-seam calibration and twelve preparations with
both loads at one port. Those calibrations live in the larger total-two sector,
not the binary history sector. They connect two uses of the same endpoint rule
and labeled interfaces; they are not measurements on these same five-step
histories. No target coordinate array enters the reconstruction. The exact
unit Gram matrix `X=4P` supplies a finite Euclidean spatial interpretation up
to orthogonal isometry. Its 66 unordered port pairs fall into exact shells:

| Graph shell | Pair count | Squared chord length in the recovered unit Gram |
|---|---:|---:|
| Native seam | 30 | `2 - 2 sqrt(5)/5` |
| Two seams | 30 | `2 + 2 sqrt(5)/5` |
| Antipodal | 6 | `4` |

The executable checks every Gram pair against the graph-distance class and
confirms that all thirty native exchange seams are exactly the shortest
metric shell. The rank-three projector is symmetric and idempotent, so `X` is
a positive semidefinite rank-three Gram with unit diagonal and zero row sums.
Port relabeling preserves the exact shell census; a symmetric change to one
edge entry is rejected. Together with
[`native_effective_fields`](../native_effective_fields/README.md), this means
the source-derived hard-core exchange terms are local on this reconstructed
finite metric. The logarithm of the transfer can still contain longer-range
configuration couplings.

This is an exact spatial interpretation of one finite twelve-port carrier,
not a measurement of the dimension of physical space. It supplies no
physical length unit, global spacetime, Lorentz cone, signal-speed law or
continuum limit. The time readouts below remain record parameters, and the
metric reconstruction uses a separately disclosed calibration preparation.

### Transport through source refinement

The response-recovered carrier Gram can be transported into the actual level-3
and level-5 source screens without changing its twelve root-port labels. The
screens have exactly the same induced native-seam graph through nine hops from
the corresponding root carrier: 222 ports and 561 seams. The exact unitless
chord metric on the root twelve ports is therefore identical in both screens.
The root graph-hop ball counts at radii zero through ten are
`12, 15, 30, 45, 54, 84, 114, 132, 186, 222, 243`; the first induced adjacency
difference appears at radius ten. The separate
[`native_action` refinement packet](../native_action/INCIDENT_REFINEMENT.md)
bounds how much this finite boundary difference can affect the complete local
repair output at a declared incident-event horizon.

An exact breadth-first readout covers every retained federation level, 0
through 5. The largest root-port distances are 13, 31, 60, 121, 242 and 484
native seam hops. At level 5, the port-ball counts at radii 16, 32, 64, 128 and
256 are 576, 2,250, 8,805, 37,224 and 123,926; the corresponding exact
doubling ratios are `125/32`, `587/150`, `12408/2935` and `5633/1692`. This is
roughly quadratic growth across part of this finite range, followed by
saturation on the finite sphere. It is a property of native seam-hop balls,
not evidence for a calibrated three-dimensional physical volume law. The
rank-three metric applies to each twelve-port carrier; it has not been glued
into a bulk metric.

This transports an operational neighborhood and one local metric, not a
Euclidean coordinate chart over the 222-port ball. The existing federation
ancestry map, supplied by the source rather than recovered from observer
records, copies each local port label to every child carrier. It preserves
the carrier-cell adjacency, but it does not map every fine port seam to a
coarse port seam: between levels 0 and 1, 74 fine seams map to non-seams and
two collapse; between levels 1 and 2, 280 map to non-seams and four collapse.
The exact per-edge Laplacian intertwining defect is one at both steps. Thus
this map supplies neither a path-independent carrier-frame transport nor a
scale law for a spatial metric.

The federation builder also has an inter-carrier routing convention: it uses
the source sphere's cell-centre directions and a declared local tangent-frame
rule to choose which ports are joined. The builder marks that routing
`source_derived: false`; it is not inferred from native observer records. At
levels 3 and 5, respectively, only `113/160` and `5951/7680` of inter-carrier
endpoint pairs are antipodal-consistent with the common local port labels.
This convention constructs the seam list, but it does not supply a
record-derived frame-transition law or a scale factor. No path-independent
physical coordinates follow from it.

The obstruction is also visible without extrapolating the ancestry map. As a
non-derivation control, take an actual inter-carrier seam at the level-3 root
and additionally require its endpoint distance to equal the local
shortest-shell length. Two independent orientations of its neighboring
carrier can satisfy that same link-length constraint, preserve both carriers'
exact internal Gram matrices, and still change other cross-carrier Gram
entries. This equality is used only in the control; it is not a model input.
The current repair records identify the seam endpoints but do not record the
missing cross-carrier inner products. So the substrate supplies exact
native-hop neighborhoods and a local unitless carrier metric; it does not yet
supply a global physical space, inter-carrier ruler, or continuum geometry.

The source-order/frame packet in
[`Geometry/SourceOrderFrameCompatibilityPacket`](../../Lean/Geometry/SourceOrderFrameCompatibilityPacket.lean)
needs its spatial placement, cone attachment and frame premises. This
packet adds no Lorentz cone, spacetime metric, signal-speed law or continuum
limit. Similarly, the finite representation in
[`QFT/SourceHistoryGNSDynamics`](../../Lean/QFT/SourceHistoryGNSDynamics.lean)
has its own supplied operator extension and flow parameter; record ordering
does not identify that parameter with these local ticks.

The concrete advance is operational: absolute native fields have an
observer-visible readback construction under an explicit preparation/access
contract, compatible local restrictions, complete continuing records, and a
testable clock comparison domain. This is model-level evidence, not a new
empirical comparison. The next step is to connect these record-accessible
fields and clock domains to the model's effective propagation and modular
descriptions while preserving their joint history and stated premises.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.relative_observer_geometry --check model/relative_observer_geometry/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/relative_observer_geometry
```

`results.json` binds its producer, runtime, native source, response recovery and
the cited Lean interfaces by SHA-256. `manifest.json` binds every packet file.
