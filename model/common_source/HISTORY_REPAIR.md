# One complete history from local constraint repair

The existing common observer history can be specified as twelve vector-valued
event nodes with local compatibility constraints. For the fixed preparation,
source operators and event graph, there is exactly one compatible complete
history. Updating any incompatible child from its current parent recovers
that history under every node-fair schedule, from arbitrary initial work
vectors. Five covering rounds suffice; a topological pass needs only one.

This gives the finite construction a concrete timeless-history description
and a local algorithm for finding it. It does not identify the auxiliary
constraint solver with a quantum hardware process or derive the preparation,
action, event graph or meaning of its records. The supplied causal instrument
is being solved as a global system of local equations.

## The history is defined independently of its repair rule

The work space at each node is the same `C^(4d^2)` used by
[experiment.py](experiment.py), with `d=2` or `4`. Every vector is
unnormalized. The boundary `x_root` is the source-prepared scalar vacuum and
gauge blank followed by the source pulse selected by the fixed native current
`j in {-1,0,1}`. The remaining nodes and equations are:

| Depth | Nodes | Compatibility equation |
| --- | --- | --- |
| 0 | Fixed boundary | `x_root = prepared(j)` |
| 1 | First evolution | `x_1 = U_1 x_root` |
| 2 | First records `a=0,1` | `x_(a) = P_a x_1` |
| 3 | Feedback branches `a=0,1` | `y_(a) = F_a x_(a)` |
| 4 | Second evolution branches `a=0,1` | `z_(a) = U_2 y_(a)` |
| 5 | Terminal records `(a,b)` | `w_(a,b) = B_b z_(a)` |

`U_1` and `U_2` are the same compiled interacting source product operators at
model-time intervals 1 and 1/2. `P_a` are the source pointer maps for the first
gauge occupation record. `F_a` are the same source conditional gauge pulses;
with feedback disabled both become the declared zero-angle controls. `B_b`
are the actual source scalar-sign pointer maps. No target state, fitted Born
probability or precomputed completed history enters these equations.

This solved object is the complete branching instrument law, retaining all
four record pairs and their unnormalized quantum outputs. Uniqueness of that
law does not select one realized record pair or eliminate Born uncertainty.

In ideal arithmetic all these maps are contractions: evolution and feedback
are unitary, while pointer outcome maps are Kraus contractions. The numerical
implementation uses the existing guarded source contractions without altering
their matrices to force a desired repair behavior.

Define `C` to be the set of assignments satisfying **all** these equations and
the fixed boundary. This compatibility set is defined before choosing an
update order. It is nonempty: recursively composing the given maps produces
an assignment. It has one element: agreement at depth zero forces agreement
at depth one, and induction forces agreement at all five depths. This is
uniqueness conditional on the apparatus inputs; it is not selection of one
physical universe from OPH's axioms.

The scalar and gauge work are coupled within `U_1,U_2`. The first retained
record really selects the feedback operation before the second evolution;
the resulting scalar record remains available to the native continuation.
Repair does not remove this operational feedback by averaging its branches.
The event semantics and record interpretation are declared, not learned by
the solver.

## Local repair, preservation and schedule independence

For a nonboundary node `v`, write its equation as
`x_v = A_v x_parent(v)`. One repair replaces only

\[
 x_v\leftarrow A_vx_{\operatorname{parent}(v)}.
\]

Every other node, particularly the boundary and already written ancestor
record labels, is left untouched. This rule reads one adjacent parent and
writes one child. It is local on the **supplied event graph**. The work
vectors of incomplete candidate histories are solver variables, not mutable
past records inside an already realized physical history.

Each element of `C` is fixed by every local update. Conversely, an assignment
fixed by every local update satisfies every equation and therefore belongs
to `C`. Thus the common fixed set is exactly the independently declared
compatibility set.

Call a schedule node-fair if every nonboundary node is selected infinitely
often. Once a node's ancestors have become correct, its next selection makes
it correct permanently. The boundary is correct from the start, and fairness
ensures that this happens at a finite update index at every depth. Since there
are finitely many nodes, every fair schedule reaches the unique completed
history after finitely many updates. It remains there thereafter.

A covering round is any consecutive finite block that selects every
nonboundary node at least once; repeated selections are allowed. After round
`k`, every node of depth at most `k` is correct, regardless of the order within
that round. The maximum depth is five, so five covering rounds suffice from
any initial assignment. A forward topological pass propagates the boundary
through all depths in eleven updates. A reverse order may require all five
rounds, or 55 updates. An unfair schedule has no such guarantee and the code
does not mark one successful merely because its requested update budget ends.

This argument is a direct finite acyclic dependency result. It does not need
an appeal to Newman confluence or an assumption that all arbitrary reduction
sequences terminate. In particular, a schedule can repeat a compatible node
forever while leaving another constraint wrong.

## A strictly decreasing scalar residual

Put `r_v=x_v-A_v x_parent(v)`. Recursively assign

\[
 w_v=1+\sum_{u:\operatorname{parent}(u)=v}w_u,
 \qquad
 \Phi(x)=\sum_{v\ne\mathrm{root}}w_v\|r_v\|.
\]

These are subtree sizes. The weights at successive levels are 12 at the fixed
root, then 11, 5, 4, 3 and 1; the root's residual is not included in `Phi`.
All weights are positive, so `Phi=0` exactly on `C`.

Updating `v` sets its own residual to zero. For each immediate child `u`,

\[
 r'_u=r_u+A_u r_v,
 \qquad \|r'_u\|-\|r_u\|\le\|r_v\|,
\]

by contractivity and the triangle inequality. No other residual changes.
Consequently

\[
 \Phi(x')-\Phi(x)
 \le\left(-w_v+\sum_{u:\operatorname{parent}(u)=v}w_u\right)\|r_v\|
 =-\|r_v\|.
\]

Every active repair strictly decreases this scalar potential. This is a
real-valued norm potential, not an integer rank: strict decrease alone is
not the finite-termination argument. The finite-DAG/fairness induction above
supplies finite recovery. In particular this construction does not supply
an integer `BareConsensusTower` witness, a new general consensus theorem, or
a certificate that the native M78 integer-repair hardware implements this
global history solver.

In floating arithmetic an executed edge may have norm slightly above one
because of rounding. The exact argument concerns ideal source maps. The
receipt separately records the maximum executed edge norm, final constraint
residual and maximum observed defect in the potential-drop inequality, with
an explicit numerical tolerance. A numerical check is not a replacement for
the exact proof.

## Timeless completion and observations inside the history

The compatible assignment is one complete history: it contains preparation,
different quantum work states at different events, both committed record
labels, and their feedback dependencies. Its being a fixed point of history
repair does not mean that all those event states equal one instantaneous
stationary state. The encoded unitary evolutions remain nontrivial.

The repair update count describes a host solver seeking that assignment.
The represented Hamiltonian uses the declared model-time intervals 1 and 1/2.
The causal ordering of recorded events is another part of the supplied
history. None of these is silently identified with a universally preferred
physical clock or with an observer's calibrated elapsed time. Different
repair schedules can use different update counts and intermediate guesses
while representing the same completed quantum observations.

Every deterministic extraction from the unique completed assignment is also
independent of the repair schedule. This includes its branch probabilities,
unnormalized terminal density operators, and any defined common observer
state/GNS representation or modular parameters obtained from those operators
using the same declared conventions. When such extracted parameters feed the
same declared clock experiment, its resulting law likewise does not depend
on which fair solver schedule found the history.

This inheritance is conditional on the extraction being well-defined on the
completed state. It does not add a new modular-domain theorem or turn an
ensemble-state-derived control angle into an online observable available to
one individual observer. Any state extraction, clock attachment and meaning
dictionary remain specified operations of the construction.

## Executable checks and retained scope

[history_repair.py](history_repair.py) builds its edges directly from the
existing source preparation, evolution, pointer and feedback APIs. It does
not call `experiment.history` to initialize or solve the constraints. The
audit calls that independent execution path only **after** the repair results
exist, to compare every terminal unnormalized vector and density operator.
The tests also make the terminal-history API raise an exception while building
and solving a history, checking that no such oracle is required.

```python
from model.common_source import history_repair as h

system = h.build_constraints(4, steps=512, current=0)
schedule = h.cover_schedule(system, "reverse", rounds=5)
solution = h.repair(system, schedule)
branches = solution["terminal_vectors"]  # (a,b) -> unnormalized work vector
assert solution["closed_within_tolerance"]
```

From the repository root:

```bash
python -m model.common_source.history_repair
python -m pytest -q model/common_source/test_history_repair.py
```

The compact audit retains schedule-dependent update counts, residuals,
potential checks, all four terminal probabilities and full-state comparison
errors. It also extracts the same reduced scalar state directly from each
repaired ensemble, retains its spectrum and the `d=2` modular frequency, and
checks equality across schedules. These are numerical state diagnostics,
with no added physical clock calibration. Tests cover arbitrary complex
inconsistent assignments, reverse and
shuffled fair schedules, all native currents, both feedback settings,
boundary preservation, unfair schedules and deterministic state extraction.

The added result is a concrete local repair formulation of this same finite
apparatus history. The occupation theorem and source temporal certificates
remain their separately stated results; their errors do not become zero
because an apparatus constraint system has been solved. The construction
does not derive the graph or dynamics, implement the global repair on the
native quantum hardware, certify an unrestricted physical-history registry,
recover all continuum sectors, or establish empirical agreement with nature.
