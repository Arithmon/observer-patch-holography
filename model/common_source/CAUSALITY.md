# Exact native causal cones

The existing native seam primitive has a finite propagation bound on its
declared 12-port, 30-edge carrier. A perturbation can cross at most one edge
per causal layer of local seam events. The actual source rule supplies
examples that attain this bound. This is an exact statement about native
loads and local placement coins, before any identification with an observer
clock or a continuum signal.

Run or replay the retained computation from the repository root:

```sh
python -m model.common_source.causality --check model/common_source/causality_results.json
python -m pytest -q model/common_source/test_causality.py
```

## Primitive and dependency proof

For an oriented active seam `(a,b)`, let `s=x_a+x_b` and placement coin
`c∈{0,1}`. The actual [native rule](../observer/native.py) replaces the two
endpoint loads by

\[
 x'_a=\lfloor(s+c)/2\rfloor,\qquad x'_b=s-x'_a,
\]

and leaves every other port unchanged. The rule conserves total load and
preserves nonnegative integer loads. Fix the entire ordered seam schedule
and coin word before comparing inputs.

Starting from a selected output port, scan the schedule backwards. Whenever
a seam touches the current dependency set, include both its endpoints and
its local coin. Induction on the number of events proves that the output is
a function only of these initial loads and coins. Consequently, two initial
states agreeing on this set give identical selected outputs when the coins
are shared. The construction gives a sufficient dependency set: inclusion
does not assert an independently variable input in every constrained sector.

Equivalently, dependencies follow time-ordered paths that wait at vertices
and cross activated edges. If an initial difference supported on `D` can
influence port `v`, some such path starts in `D` and ends at `v`. Therefore

\[
 \operatorname{dist}(D,v)
 \leq \text{edge crossings on that causal path}
 \leq \text{number of elapsed seam events}.
\]

The first inequality is the static graph metric bound; the required edge
ordering makes it strictly more informative than connectivity alone. A
sequential sweep through all 30 edges is 30 events, not one parallel round.
If simultaneous disjoint events are grouped into a causal layer, at most
one edge can be crossed per layer. There is no simultaneous-matching
execution claim in this receipt.

## Executed signals and exhaustive checks

The carrier and labels are read from the existing pinned source through
[slow_geometry.py](slow_geometry.py). [causality.py](causality.py) exhausts
all 78 total-load-two states for each of 27 declared schedule/coin pairs,
checking every output port against its dependency set. All coin words are
enumerated for the five short schedules; four specified words are checked
for each forward and reverse full sweep. The general dependency proof does
not depend on this finite selection of coin words. The receipt retains the
actual schedules, groups and comparison counts.

Two concrete source executions establish delayed influence:

* States with occupied ports `{1,2}` and `{1,3}` have difference support
  `{2,3}` at minimum distance two from receiver `0`. The ordered seams
  `(2,10),(10,0)`, both with coin zero, produce receiver traces `(0,0,1)`
  and `(0,0,0)`. Reversing the seams gives `(0,0,0)` in both cases.
* From the shared state `{1,3}`, run `(3,2),(2,10),(10,0)`. Changing only
  the first coin from zero to one produces traces `(0,0,0,1)` and
  `(0,0,0,0)`. Immediately after that first local perturbation, its changed
  endpoints are still distance two from the receiver. Two more properly
  ordered seams are required.

The mass constraint matters. Two different total-two states cannot differ
at only one port. The graph has diameter three, with one antipodal port per
receiver, but this does not supply a pair of admissible states differing
only at that antipode. The first example uses two admissible binary states;
the second perturbs one local coin and follows its two-endpoint influence.
Neither replaces that constraint with an impossible single-port input.

## Units and remaining operational identification

Taking one graph edge as the distance unit and one causal layer as the
duration unit gives native signal speed at most one; these witnesses
saturate the bound along their paths. Observers may choose units with a
reference signal speed equal to one. There is no requirement to derive the
names or numerical sizes of metres and seconds, or to pick unique units.

The additional physical question concerns relations between measurements:
whether an actual observer clock, propagation over inferred distances,
different directions and different matter sectors agree on the required
signal laws. A host loop counter or a history-repair iteration is not
automatically such a clock. This receipt does not establish a relation
between these seam events and the modular clock readouts, isotropic
continuum propagation, Lorentz symmetry or the speed measured in nature.
The carrier's diameter three is not a claim of three spatial dimensions.

The locality statement is for the native primitive with a fixed external
schedule and shared coins. An adaptive controller that reads remote data
and changes a seam selection can itself communicate information; its
communication must be modeled before extending this proof. Shared record
storage, the quantum compiler's transfer/load/export services and the
entire observer protocol have not been localized by this calculation.
Thus native causal cones are one concrete consistency check, not a proof
that every operation in the full apparatus has this propagation bound.
