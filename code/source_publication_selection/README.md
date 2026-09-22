# Eventual publication and source selection

The [derivation](DERIVATION.md) removes an arbitrary publication deadline
and derives a consistent native history law. It also tests the proposed
selection mechanism on a common demand algebra and on radius/population
alternatives. The [contract](CONTRACT.md) states the declared interface and
scientific boundary.

The [finite-state derivation](STATE_SELECTION.md) addresses A3's actual
noncommutative local-state objective. It proves maximal feasible support,
the exact criterion for forcing a positive public failure effect to zero,
quantitative failure floors, and an explicit bound against a false
approximate minimum. These are analytic matrix results with kernel-checked
real reductions. Exact noncommuting matrices and independently reconstructed
channels test the distinction between a selected state and a selected
record-preserving process. Neither a global state nor a shared eigenbasis
is assumed in the general theorem.

From the repository root, with `PYTHONPATH=code`:

```text
python -m source_publication_selection.verify
python -m pytest -q code/source_publication_selection
```

To regenerate the compact control packet, run
`python -m source_publication_selection.build`; the independent verifier
returns the receipt through `verify.verify(codec.load(...))`. Both artifacts
are canonical JSON. Source custody covers the contract, derivation, Python
implementation and transitive local Lean dependencies. The fixed experiment
specification lives in `codec.py`; it contains no response or geometry solver.

The producer traverses rational response checkpoints. The verifier enumerates
words independently and runs two integer-scaled scalar preparations. Geometry
uses integer radical signs in the producer and certified rational enclosures
in the verifier. No numerical tolerance decides a metric boundary.
The path-family producer uses positive exterior-power dynamics; its verifier
replays scalar columns and recomputes every maximal minor independently.
Locality controls compare closed formulas with independent DAG closure and
grid breadth-first search. They show why deleting all spatial restrictions
or keeping only a fixed microscopic stencil fails the causal target.

The formal gate is `Geometry.SourcePublicationAxiomAudit`. The native ordered
two-row minor invariant for arbitrary words and the generic
harmonic law, finite-marginal convergence, support, tail contraction, demand
marginal, radius/counting formulas, complete-read order, interval-fraction
bounds, and the balance of the actual volume-error estimate are kernel checked. Infinite-path
measure/entropy arguments, the uniform native block instantiation, and the
probability classification and general k-column Cauchy--Binet argument are
analytic; finite census checks do not
replace those proofs. Neither the paper nor the receipt claims full M1.
