# Established mathematics for observer agreement

Research method, 11 October 2026. This document governs theorem reuse in the
[model search](MODEL_SEARCH.md), [candidate register](MODEL_CANDIDATES.md), and
their existing [premise obligations](MODEL_PREMISE_MAP.md). The target remains
one physical process with compatible carrier, response, record and clock
readouts. Established results about distributed recovery and sheaf diffusion
are tools for that construction. They do not select its physical law.

For each proposed result, identify the nearest theorem below, map its objects
and hypotheses to the actual candidate, and isolate the remaining obligation.
Reuse the conclusion once those hypotheses hold. A renamed consensus result
or a longer simulation of it is not an independent theoretical contribution.

## 1. Recovery is a property of a system and its execution assumptions

[Dijkstra (1974), *Self-stabilizing systems in spite of distributed control*](https://www.cs.utexas.edu/~EWD/transcriptions/EWD04xx/EWD426.html)
provides local finite-state constructions that enter a legitimate set from
arbitrary initial states and remain legitimate. His central daemon selects
one enabled move at a time; the guarantee covers every such selection.
A distributed daemon is explicitly outside that paper's scope. Legitimate
states have enabled moves, so legitimate operation continues after recovery.
This is not a theorem that arbitrary local rules stabilize, or that every
legitimate state is a fixed point.

**OPH application contract.** Before invoking recovery language, supply:

- The complete state and an independently specified legitimate set `L`.
  Defining `L` as whatever the chosen solver eventually reaches is insufficient.
- The initial domain, permitted perturbations and conserved sectors. Identify
  a fault span closed under the modeled fault and protocol moves. If faults
  can alter supposedly protected records, a proof that assumes those records
  are intact covers a narrower fault model.
- The enabled moves and execution class: central or concurrent scheduling,
  atomicity, stale reads, delivery, crashes, random choices and any fairness.
- **Closure:** each allowed fault-free move from `L` stays in `L`.
- **Recovery:** every maximal fault-free execution from the admitted fault span
  reaches `L`, with deterministic, almost-sure or expected guarantees stated
  separately. State the unit of a recovery bound.

This is the programme's contract for recovery after a finite disturbance.
Continued forcing requires a tracking or noise-tube result. Recovery does not
promise safety throughout the recovery interval or restoration of information
that the allowed disturbance erased. Fair scheduling alone gives no wall-clock
bound; a time interpretation needs additional timing premises.

**Immediate use:** apply this contract to the native candidate's invariant
binary shell, the larger integer-load domain, and protected record registers
separately. Keep ensemble convergence distinct from convergence of individual
recorded trajectories. Reuse the existing finite completion proof where its
domain and scheduler match; do not assert Dijkstra's ring theorem for this
different transition system.

## 2. Termination, confluence and stabilization answer different questions

[Newman (1942)](https://doi.org/10.2307/1968867) supplies the established
termination-plus-local-confluence route to confluence. The
[finite consensus paper](../paper/reality_as_consensus_protocol.tex),
`thm:quotient-repair-normal-form`, uses that route on its accepted quotient
relation, with separate completeness and protected-boundary premises.

Continue to verify those concrete premises rather than reproving the abstract
lemma. A unique terminal state from one initial quotient does not imply that
different histories share one terminal state. That stronger statement needs
the boundary-fibre uniqueness condition. Conversely, a legitimate operating
set may contain cycles and need no terminal normal form. Preserve the
distinction when specifying an observer's continued record-producing activity.

The existing noisy fair-block branch also has its own contraction, stopping-time
and occupation assumptions. It is not a corollary of the word “asynchronous.”

## 3. Use cellular sheaves for a declared linear overlap model

For a finite graph, let each vertex carry a finite-dimensional inner-product
space `F(v)`, each edge a comparison space `F(e)`, and each incidence a fixed
linear restriction `R_ve`. With an orientation on each edge `e=(u,v)`, set

```text
(delta x)_e = R_ve x_v - R_ue x_u,
L_F = delta* delta,
E(x) = ||delta x||^2 / 2.
```

[Hansen and Ghrist, *Opinion Dynamics on Discourse Sheaves*](https://arxiv.org/html/2005.12798),
§§2 and 4, gives `ker L_F = H^0(G;F) = ker delta` and convergence of
`dx/dt = -alpha L_F x`, `alpha>0`, to the orthogonal projection of the initial
state onto that space. Agreement compares restricted expressions, so it need
not make distinct vertex states identical. The theorem assumes the fixed
linear sheaf and this evolution equation; it does not establish arbitrary
nonlinear, asynchronous, adaptive or quantum-channel convergence.

**Immediate use:** factor the existing weighted scalar disagreement into its
actual restriction maps and adjoints. Record the mass/inner-product convention,
nullspace, protected coordinates and boundary data before comparing spectra.
If a native quantum candidate admits a suitable operator-space description,
prove that identification and its positivity, invariant-sector and dynamical
conditions first. Calling density matrices “stalks” does not discharge them.
Keep conservative Hamiltonian evolution distinct from dissipative heat flow.

Changing the restrictions can reduce the reported discrepancy by changing what
is compared. Therefore an adaptive proposal must preserve a declared observation
map or expose and test its changed meaning. A zero residual under a modified
readout cannot validate the original physical interpretation.

## 4. Boundary uniqueness does not imply exact agreement

[Hansen and Ghrist, *Toward a Spectral Theory of Cellular Sheaves* (2019)](https://arxiv.org/html/1808.01513),
§4.1, establishes harmonic extension results. For 0-cochains with values pinned
on vertices `B`, an extension minimizing `||delta x||^2` always exists; it is
unique exactly when `H^0(G,B;F)=0`. With interior vertices `I`, it solves

```text
L_II x_I = -L_IB x_B.
```

This interior harmonicity is weaker than `delta x=0`. Boundary data must lie
in the image of `H^0(G;F) -> C^0(B;F)` to admit a global section. When the
boundary data are a section on a boundary subcomplex, the connecting map into
relative first cohomology expresses the extension obstruction. For the
programme's graph computations, explicitly testing the linear image is enough.

For example, an ordinary three-vertex path pinned to 0 and 1 at its endpoints
has the unique harmonic interior value `1/2` and positive disagreement.
Do not label it complete agreement. Compute both boundary uniqueness and the
remaining residual. On homogeneous unconstrained linear sheaves, the zero
section always exists; nonzero `H^1` is not a blanket assertion that no global
section exists.

**Immediate use:** before adding a boundary selector, test whether the current
protected data extend and whether a nonzero global section vanishes on the
boundary. This separates an inconsistent target from a free interior mode or
a failing solver. It also tests the hypothesis needed by OPH's existing
boundary-conditioned uniqueness result.

## 5. Keep relational, probabilistic and quantum extension separate

[Abramsky and Carù (2019)](https://arxiv.org/html/1911.03521), Eq. (4),
Theorem 4.1 and Propositions 5.1/6.2, formalize local agreement and global
extension. For relations, a nonempty natural join supplies a satisfying
assignment, while an exact global relation must recover every local relation
under projection. Probabilistic extension additionally requires the specified
weights. Full local support can hide a Bell obstruction. The natural-join
criterion for relations does not justify multiplying probability marginals
and calling the product their extension.

**Immediate use:** name the valuation type of each candidate's agreement
target before choosing a solver. Check probability compatibility separately
from support compatibility. A positive joint law for all counterfactual
measurement outcomes is not the same object as a global quantum density
operator with context-dependent measurements. A classical extension failure
rejects that classical interpretation; it does not reject every quantum
realization. A2 itself
[does not imply global extension](AXIOM_REFERENCE.md#a2-observer-agreement).

## 6. Apply distributed impossibility results only to matching models

[Fischer, Lynch and Paterson (1985), Theorem 1](https://groups.csail.mit.edu/tds/papers/Lynch/jacm85.pdf)
rules out guaranteed termination of deterministic nontrivial binary consensus
in its fully asynchronous message-passing model with a possible crash, even
with eventual reliable delivery to nonfaulty processes. It is a relevant
obstruction if an OPH implementation makes those commitments. It is not a
no-go theorem for fault-free finite rewriting, asymptotic diffusion, or a
protocol with genuinely supplied additional synchrony assumptions.

Before invoking this result, map processes, decisions, validity, communication,
failures and admissible schedules. Before claiming to evade it, identify the
changed hypothesis and its implementation cost. The
[BFT appendix](../paper/appendix_B_bft_qecc_extensions.tex) separates
protocol safety from the partial synchrony needed by its liveness branch.
No network theorem supplies a physical clock merely by naming an event a tick.

## 7. Work products and scientific contribution

The following bounded outputs use the existing issue ownership. They add no
new issue queue or mechanism and retain the existing candidate evidence.
Reuse a completed certificate when the mapping is unchanged; calculate only
what the actual candidate's missing premise or changed interface requires.

| Work | Required output | Existing owner |
| --- | --- | --- |
| Specify native recovery | One closure/recovery contract with fault span, scheduler and protected records; distinguish deterministic and probabilistic claims. | [#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076), [#1026](https://github.com/FloatingPragma/observer-patch-holography/issues/1026) |
| Reuse the linear theory | Actual stalks/restrictions, inner products, `delta`, nullspace, boundary extension residual and uniqueness test for the scalar comparison. | [#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076) |
| Preserve one implementation | Accepted operation grammar and same-source preparation, evolution and record history; unchanged reference checks plus the relevant theorem premises. | [#1077](https://github.com/FloatingPragma/observer-patch-holography/issues/1077) |
| Discriminate physical laws | Two viable laws with the same preparation, effects, boundary and clock; a predicted transient or response difference with declared calibration freedom. | [#1025](https://github.com/FloatingPragma/observer-patch-holography/issues/1025), [#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076) |
| Review the inference | Nearest primary theorem, premise mapping, unchanged conclusion, candidate-specific remaining work, and scope of any obstruction. | [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033), [#1078](https://github.com/FloatingPragma/observer-patch-holography/issues/1078) |

Record each reusable result as **applicable with checked premises**,
**conditional on named missing premises**, or **not applicable to this model**.
These are research-workflow labels, not changes to the scientific claim registry.
Keep links to the supplied premise evidence; a citation alone is not that evidence.

Recovery, linear diffusion, harmonic extension, confluence and the distinction
between compatibility and extension are established mathematics. An OPH-specific
contribution must identify the common physical realization or a necessary
additional result absent from the cited theorem, and ultimately supply a
discriminating observable relation. Shared equilibria do not choose transition
rates, instruments, reference states or physical clocks. The existing geometry,
gauge, continuum and empirical obligations remain in the premise map.

The sources above have different conventions and scopes. In particular, the
arXiv HTML renderings can renumber theorem counters relative to the PDFs;
the linked sections and displayed equations identify the intended results.
