# A local action for native quantum observations

The native C/E instruments fix an interacting local generator on the binary
operator algebra. No field coupling or Gaussian action needs to be fitted to
obtain that generator. This packet also bounds how much an observer's **whole
local quantum and record output** changes when distant ports are omitted.
The bound includes arbitrary untouched references and interventions controlled
by the observer's own local records.

This adds a quantum and complete-history extension to the existing
[density-only propagation bound](../federated_observer/SCALING.md). It reuses
the [signed-difference theorem](../observer/QUANTUM_MIXING.md) and
[actual federation instruments](../federated_quantum_sector/README.md);
the general block formula is not a new result. The new attachment supplies
a local action representation, an observational truncation theorem and a
conditional local-limit construction for the full quantum process.

## What the declared inputs supply

There are binary factors at the ports of a finite graph, and the same
native endpoint repair acts on each seam. C and E remain two specified
candidate instruments; this calculation does not select between them.
The existing architecture supplies the twelve-port carriers, their seam
wiring and the allowed quantum interface. Uniform, independent seam
selection is a declared scheduling protocol. Coherent preparation and local
phase interventions in the finite controls are declared services.

For a graph with `M` seams, let `Phi_e` be one native repair with its records
discarded, and let `Phi=(1/M) sum_e Phi_e`. These maps are constructed from
the source endpoint function and the actual C/E record equality. Their
extension by identity on other ports is exactly the existing finite tensor
split. A federation needs no larger elementary carrier.

## The source fixes a local Liouville action

On the Hilbert space of operators with the Hilbert--Schmidt inner product,
define

```text
h_e = I-Phi_e,                 H_native = sum_e h_e,
Phi = I-H_native/M.
```

The endpoint audit constructs orthogonal projectors such that
`h_e=P_1+(1/2)P_half` for C and `h_e=P_1` for E. Each term is positive and
acts only on the two endpoint operator factors. These facts hold on the
full binary tensor algebra, including coherences between number sectors.
The identity is a common zero vector, so the local operator is frustration
free. Its kernel need not be one dimensional: vacuum/full coherence remains
fixed on the unrestricted algebra. The earlier fixed-number uniqueness
theorem has a narrower domain.

The signed-difference decomposition makes the action explicit. For a matrix
unit `|x><y|`, write `z=x-y`, and `Z={v:z_v=0}`. The common occupied sites
form a subset `S` of `Z`. On each fixed common-number block,

```text
H_native[z] = (1/2) L_exclusion(G[Z]) + V_L(z) I,
V_C(z) = (1/2) sum_edges |z_a-z_b|,
V_E(z) = sum_edges 1[z_a != z_b].
```

Thus its quadratic form on block amplitudes `f(S)` is

```text
<f,H_native[z]f>
 = (1/4) sum_(S,ab in G[Z]) |f(S)-f(S^ab)|^2
   + V_L(z) sum_S |f(S)|^2.
```

An exchange of equal occupations contributes zero. These coefficients come
from the native fair placement coin and record rule. The hard occupation
constraint and seam exchange couple distinct ports: for example the local
generator takes the population `|10><10|` to half the difference between
that population and `|01><01|`. It cannot be a sum of independent on-site
updates. Coherence differences additionally carry the explicitly derived
C/E killing terms. A Gaussian field with an independently chosen potential
is not inserted.

`H_native` is a positive **Liouville** operator acting on operators. It is
not a Hamiltonian acting on the source's state vectors. In particular the
native CP dynamics mixes some pure states and is not unitary. For the
original discrete scheduler, the exact transfer remains `I-H_native/M`;
replacing it by `exp(-H_native/M)` would change predictions. The logarithm
of a transfer, when defined, also differs from this local sum.

## Exact path weights and the optional continuous parameter

An optional bookkeeping construction gives each seam an independent
rate-one Poisson attempt process. Its parameter `s` counts expected
attempts **per seam**. The total count is `N~Poisson(Ms)`, and conditional
on `N=n` the ordered seam word is exactly the original uniform scheduler.
Consequently the unread channel is

```text
T(s) = sum_(n>=0) exp(-Ms) (Ms)^n/n! Phi^n
     = exp[-s H_native].
```

This is also a local Lindblad generator: writing each native Kraus map as
`K_(e,a)`, completeness gives
`-H_native(X)=sum_(e,a)[K_(e,a)XK_(e,a)^* - {K_(e,a)^*K_(e,a),X}/2]`.
No Hamiltonian term or fitted rate is required for this particular
scheduler extension. The Poisson count distribution itself is declared;
it is not derived as an observer's physical clock.

More explicitly, one ordered seam-and-coin path of length `n`, after
integrating its ordered attempt locations within `[0,s]`, has coefficient
`exp(-Ms)(s/2)^n/n!`. A matrix unit contributes only if its left and right
native paths give the same hidden Kraus label at every step. Visible
currents are then retained, while hidden alternatives are summed as CP
maps. This is a local path-weight action with hard record-compatibility
constraints. Complex input coherences and final quantum readouts still
contract with this kernel; the construction is not a classical
hidden-variable replacement for arbitrary quantum experiments.

For protocols depending on event order and records, every joint record law
at a deterministic native count remains available by conditioning this
auxiliary count construction on `N=n`. A protocol using externally assigned
Poisson locations needs those locations retained as additional inputs. Clock readings
formed from the native records can therefore be retained directly.
A physical duration, modular-flow identification, Gaussian measure and
Lorentzian action do not follow from this reparameterization.

## Complete local observation theorem

Let `A` be the fixed set of observed ports, `Delta` the maximum graph degree,
and `B_R(A)` its closed graph-radius-`R` ball. Compare the full process with
a simulation retaining this ball and its initial reduced state, including
any untouched reference. Delete seams crossing the ball boundary. In a
fixed-count discrete simulation, keep deleted attempts as identity holds
under the original scheduler, preserving its seam probabilities.

The observer may keep every native current/coin record of seams wholly
inside `A`, execute local instruments on `A`, and choose later such
instruments from those local records. It may also retain the final quantum
state on `A` jointly with all those records and an arbitrary reference.
Local memories storing these records are included. The additional intervention
protocol is finite, or locally finite and well posed on every bounded horizon;
an accumulation of infinitely many control operations is excluded. The controller may not
consult external records, act outside `A`, or observe the global attempt
tape as an extra signal. The native scheduler is state independent.

For the **unconditional complete output**, in trace distance, the error is
at most

```text
discrete n attempts:
  min{1, |A| binom(n,R+1) (Delta/M)^(R+1)},

optional per-edge Poisson parameter s:
  min{1, |A| (Delta*s)^(R+1)/(R+1)!}.
```

The binomial coefficient is zero when `n<R+1`. The same bounds control all
measurements of the retained output, and hence joint field and local-clock
readouts within the declared horizon. Selected rare outcomes can amplify
conditional errors; the bound does not silently become a postselected
bound. An unbounded stopping time requires an additional stopping-tail
estimate. Finite caps must retain their timeout outcomes.

**Proof.** Couple the full and truncated source schedules inside the ball.
For a fixed seam word, trace the dependency cone of `A` backwards. An edge
can enlarge that cone only when it meets the existing cone. Record memories
and their allowed feedback remain in `A`. If the cone never crosses the
retained boundary, the whole local instrument is supported within the ball.
The original initial state and its restriction then give identical outputs,
including arbitrary references. This is the actual endpoint instance of
[DependencyCone](../../Lean/ObserverPatchHolography/Locality/DependencyCone.lean)
and the nonselective finite tensor identity in
[SlotLocalitySurface](../../Lean/EventAlgebra/SlotLocalitySurface.lean).
The causal-access requirement follows the existing
[AdaptiveScheduler](../../Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean)
contract; no new Lean proof is claimed.

Crossing the boundary requires a temporally ordered walk of at least
`R+1` edges starting in `A`. There are at most `|A| Delta^(R+1)` possible
walk prefixes. In the uniform discrete schedule, a particular ordered
edge list occurs at specified `R+1` attempt positions with probability
`M^(-(R+1))`; union over their `binom(n,R+1)` choices. For independent
rate-one seam Poisson processes, the expected number of ordered occurrences
is `s^(R+1)/(R+1)!`, including repeated edges by the Poisson factorial
moment identity. A union bound gives the two displayed probabilities.
On every good schedule the complete outputs coincide; on a bad schedule
the trace distance of normalized outputs is at most one. Averaging proves
the result without any dimension factor or reference-size restriction.

The idea is related to general dissipative
[Lieb--Robinson locality](https://arxiv.org/abs/1003.3675).
Here the explicit native Poisson unraveling permits the elementary proof
above, including local records. Its bound does not assume that a
single-particle swap coupling also controls coherent histories.

## Finite executions and conditional federation limit

The runtime uses the actual source path `(0,1,82,72)` across two elementary
carriers, with its three seams uniformly selected for three attempts. This
is a declared three-seam experiment, not the full-federation scheduler.
All binary states on the four active ports are included; all other ports
are empty and are not touched by this declared schedule. The initial probe
is entangled with a reference qubit, and a remote port starts occupied.
Current records on the first seam control a local phase intervention.

For each candidate, all 27 schedules are executed both in full and with
the outermost seam replaced by identity. The runtime retains 125 complete
local record branches across these schedules, full terminal matrices and
a serialized checkpoint continuation. Exactly `7/27` of schedules allow
a boundary-crossing cone. All others agree to roundoff. An independent
rational test checks complete local output equality for **all 256 input
matrix units on every good schedule**, including the feedback operation.
The bad-schedule controls have a nonzero discrepancy, so truncation is not
assumed exact at arbitrary horizons.

A second executed control reads an external port and flips a local port
from that record. Two initial states agreeing on the radius-two ball then
produce perfectly distinguishable local outputs with zero native attempts.
That controller violates the local-access hypothesis. Local native gates
alone do not constrain information supplied through a global controller.

The receipt also evaluates the uniform bound for the existing 240, 960 and
3,840-port federations, whose maximum degree is six. For two observed ports,
`R=8` and `s=1/4`, the complete quantum-history error allowance is
`243/1146880 < 0.000212` at every size. The discrete attempt allowances and
actual retained ball sizes are reported separately. This is a simulation
error guarantee under the stated scheduler and local-access assumptions,
not an empirical error bar or a measured physical propagation speed.

For a **fixed bounded-degree infinite graph** with consistent local initial
states, use independent per-edge Poisson processes with the auxiliary parameter
`s`. A uniform choice among all edges of an infinite graph is unavailable.
The Poisson bounds make finite-ball quantum observation laws Cauchy for every
bounded `s` horizon. The two-ball comparison uses a common
inner ball and the sum of their truncation allowances. Positivity,
normalization and compatibility pass to the limit. This constructs a
continuing local quantum process and all specified local instrument
histories from federations of the same elementary factors. The limiting maps compose to give CP dynamics on the quasilocal algebra.
An individual observer subalgebra need not evolve autonomously: its reduced
description can retain the memory exhibited by the joint clock process.

This conditional thermodynamic construction does not require a new larger
elementary carrier. It does require the graph neighborhoods, source
instruments, initial restrictions and observer access maps to agree. The
pinned screen refinements have not been shown to obey those hypotheses;
their earlier coarse-graining obstruction remains. The result supplies no
Lorentz symmetry, continuum dimension, or physical identification of `s`.
Those require a further observational correspondence and scaling argument.

```sh
python -m model.native_action --check model/native_action/results.json
python -m pytest -q -W error model/native_action
```
