# Native exchange histories have a nontrivial finite class

The binary two-load carrier has a positive history-topology result: **its
configuration complex has first integral homology `Z/2`** under the explicit
local deformations below. A three-move history exchanges the two native
worldlines inside one triangular face. Complete C and E records reconstruct
that exchange, and it cannot be removed by those deformations. No physical
particle-label register or assigned seam phase enters this construction.

This is a class of complete paths. It is not a conserved charge of the
instantaneous source, a fermionic exchange phase or a protected particle.
The [particle-readiness test](#what-is-missing-for-a-particle) makes those
distinctions executable. The [successor contract](joins_campaign.json) fixes
the domain, and [`defect_history.py`](defect_history.py) retains the finite
certificates and complete output laws.

## Compatibility and the native configuration complex

A vertex is an unordered pair of distinct occupied ports, giving 66 vertices.
An edge moves one occupied port along a native seam to an empty port, giving
300 edges. This is specified from carrier incidence before applying an
update. The complete C/E instruments have exactly this movement graph:
every edge has positive branch weight, and every nonmovement is a stutter.
All 240 C and 480 E raw maps, including the zero maps, enter the check.

The admitted elementary deformations are independently specified as:

- insertion or deletion of a stutter or an edge followed by its reverse;
- movement of one particle around a native triangular face while the other
  particle is at a port outside that face: `20*9=180` triangular cells;
- interchange of two moves along vertex-disjoint native seams: 315 square
  cells.

These moves never permit the two particles to coincide. The declaration does
not identify every endpoint-preserving rewrite with a physical equivalence,
nor does it identify the configuration complex with the twelve-port spatial
complex. Its higher-dimensional cells and full fundamental group are not
computed here.

## The invariant comes from lifting actual paths

At each configuration use the sorted occupied pair as a chart. Lift a native
movement by following the particle that its recorded current moves. Comparing
the resulting ordered pair with the target's sorted chart gives an edge
label in `{+1,-1}`. A zero current gives `+1`. This procedure derives the
labels from native trajectories; it does not attach a sign to a Kraus
amplitude.

The labels reverse by group inverse. Their product is `+1` around every
declared triangle and square. Relabeling the two particles separately at
each configuration multiplies an edge by its endpoint chart factors; these
cancel around a closed loop. Thus the resulting closed-loop parity is
chart independent and invariant under the declared local deformations.

For the native face `(0,11,5)`, consider

```text
{0,11} -> {5,11} -> {0,5} -> {0,11}.
```

The particle initially at 0 ends at 11 and the particle initially at 11
ends at 0. The loop's holonomy is `-1`. All configurations remain on the
same triangular face; no occupied port or seam is inserted into the source.
The opposite labeling of the initial pair gives the same closed parity.

The mathematical transport and chart arguments instantiate the generic
theorems `transport_append`, `transport_reverse`,
`transportRatio_eq_closedLoopHolonomy` and `holonomy_rechart_invariant` in
[`HolonomyInterference.lean`](../../Lean/Screen/HolonomyInterference.lean).
Here its vertex type is the **66 configurations**, not the twelve ports.
The theorem's group and reversal predicates are supplied by the actual lift.
For its total edge-label function, extend the sign by `+1` on nonedges;
physical paths remain restricted to native edges and stutters. This extension
satisfies the universal reversal predicate without adding physical moves.
Its character theorem gives the algebraic sign character; an implemented
interference experiment remains a separate requirement. No fixed Lean seam
numbering or the nineteen-cycle spatial seam classification is imported.

## Exact integral certificate

The configuration graph is connected, with a retained 65-edge spanning tree.
Its integral cycle lattice therefore has rank `300-66+1=235`. Projection
onto the non-tree edges is an integral coordinate system for this lattice.
All 495 local-cell boundaries have zero vertex boundary.

The receipt verifies three exact facts:

1. The local-cell boundary matrix has rank 234 over `GF(2)`.
2. A specified 235-cell minor in the integral cycle coordinates has
   determinant 2.
3. A retained integral combination of 62 actual cells has boundary exactly
   twice the triangular exchange loop, checked on all 300 edges.

The determinant-two sublattice makes the full quotient finite with order
dividing two. The rank defect modulo two excludes order one. Hence the
integral quotient is exactly `Z/2`; the odd loop detects its nonzero element.
Independent tests also compute the full rectangular Smith form, with 234
unit entries and one entry two. This is not a numerical rank estimate.

Parity is complete for this **abelian homology quotient**. No completeness
claim is made for path homotopy. Declaring the exchange triangle itself an
additional filled cell raises the boundary rank modulo two to 235 and
kills this invariant. That negative control changes the deformation
contract. Native binary repair does not create double occupations; extension
to a collision domain would need its own path and deformation maps.

## Complete native loop experiments

The fixed local schedule uses seams `(0,5)`, `(0,11)`, `(5,11)`. From the
known pair `{0,11}`, all eight fair-coin words are retained:

| Three-attempt output | Exact probability |
| --- | --- |
| Closed, even exchange parity | `1/4` |
| Closed, odd exchange parity | `1/8` |
| Open path | `5/8` |

The complete three-attempt quantum kernels have 38 nonzero C blocks and
56 nonzero E blocks on arbitrary binary inputs. The raw menus include all
512 C and 4096 E record slots. Their partial-injection Gram sums preserve
trace, including untouched references. The diagonal loop probabilities
agree, while the full C/E records and surviving coherences differ.

That fixed schedule does not explore every seam. A separate complete run
uses each preferred seam with weight `15/16` plus uniform exploration with
weight `1/16`. Every native seam and coin has positive probability. The same
specified three-move exchange word has probability `(451/960)^3`; the full
unconditioned run additionally retains other exchange paths, open paths and
every excursion outside the face. The receipt stores their exact rational
masses. The calculation starts from a declared known pair; it does not
execute or waive the source's preparation service.

Only the retained records and known initial preparation are needed to
reconstruct the ordered path. Its 132 possible ordered pairs need at most
eight bits in an offline tracker. This is a classical calculation, not a
new physical label register, a protected qubit or an implemented online
controller. Full records retain their own storage cost.

## A finite signed compatibility form

The derived labels also instantiate
`kernel_trivial_of_negative_cycle` and `gram_posDef_of_negative_cycle` in
[`SignedGap.lean`](../../Lean/OphGap/SignedGap.lean). Its actual incidence
has 300 rows and 66 columns, with row `-u_x+sign(x,y)u_y`. The verified
spanning tree supplies connectedness, and the exchange loop supplies the
negative cycle. The signed Gram `L_sign=D_sign^T D_sign` is positive definite.
An independent exact Sylvester certificate gives `L_sign>2I`; the degree
bound gives `L_sign<=20I`.

This is a positive finite compatibility cost at fixed holonomy. It is not
the native quantum Hamiltonian or an action minimized by C/E repair.
The unsigned native population form has a constant zero mode. Choosing the
signed form as a physical action would require a separate dynamical map.

## What is missing for a particle

The same source endpoint admits both even and odd complete histories. A
future native exchange loop can change the accumulated parity and return
the source to the same configuration. The invariant therefore classifies
histories under specified deformations; it does not protect a live
source-only charge. This preserves the complete-branch commutant result.

The canonical native Kraus amplitude on the odd three-move word is
`+sqrt(1/8)`. A mathematical sign does not turn it into a fermionic phase.
Multiplying any **committed** history block by a scalar sign leaves its
complete source/reference output unchanged:

```text
(-K_h) rho (-K_h)* = K_h rho K_h*.
```

The record blocks are classical alternatives. A proposed phase-sensitive
transport must specify an actual coherent operation and readout. An
input-dependent phase can change a quantum branch, but that is a changed
operation, not a consequence of naming the path class. No fermionic
statistics, spin or species assignment follows from `Z/2` alone.

There is a further exact native history observable. Under uniform repair,
let `S^n[i,j]` be the probability of reaching `j` with even lifted parity
minus the probability with odd lifted parity, starting in chart `i`. The
actual complete current records give

```text
S = I - L_sign/60,       ||S||_2 <= 29/30.
```

At three attempts the unconditioned odd return probability at the local
starting pair is `1/54000`. The receipt retains exact even, odd and open
masses through 64 repairs. This character transfer has signed entries and
is **not** a native quantum source channel. Its contraction is in native
repair count. It neither erases a retained classical history nor makes its
gap a rest energy or a particle decay width.

The [adjacent-pair test](LOCALIZED.md) supplies the remaining localization
control: under uniform repair an adjacent pair exits adjacency after 15
attempts on average and fails the tested reference-like exterior condition.
The closed loop does not fix that ongoing-state failure.

A particle construction therefore needs a live localized sector with the
reference exterior, a native representation/readout of any claimed exchange
phase, and protection under the complete continuing dynamics. Rest action,
moving response and stress must then come from the same process. The first
concrete test is to compare a proposed phase-carrying transport against the
unchanged complete C/E quantum outputs; decorating a classical history
with its sign cannot pass that test by itself.

```sh
python -m model.clock_defect_energy.defect_history --check
python -m pytest -q -W error model/clock_defect_energy/test_defect_history.py
```
