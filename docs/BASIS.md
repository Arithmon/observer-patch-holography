# The OPH basis

Observer Patch Holography has one hypothesis, one rule that follows from it,
one condition that defines what is shared, and one executable specification.
The conditional theorems of the [theorem library](../paper/README.md) keep
their displayed premises; the [premise map](../model/PREMISES.md) records
which of those premises the executable specification supplies.

## The closure hypothesis

The Universe is a timeless, self-referential structure. Its observers recover
the specification of their own world from completed records and build
hardware whose execution realizes that same history, including its builders
and their hardware. The structure reads itself into existence.

Write `U` for the complete inhabited history, `D(U)` for an account of its
architecture recoverable through records inside it, and `R` for an admissible
realization of that account constructed inside `U`. The closure condition is
`R(D(U)) ≅ U`, where `≅` is agreement in the observable and interventional
structure. As a fixed point of a universe-level closure operator it reads
`T(U) = U`. No external host, external time or global clock is supplied;
construction and execution happen in emergent, mutually calibrated local
times. Fixedness of the complete history permits changing states within it.

The equation is schematic until a candidate specifies the history domain and
the action of `T`. Existence, uniqueness and the identification of our
Universe with such a fixed point are open questions of the programme. The
[working terms](INTERPRETATION_AND_TERMS.md) distinguish complete history,
instantaneous fixed state, stationary distribution and continuing regime, and
[Essay B](../essays/B-why-does-anything-exist.tex) develops the hypothesis
philosophically.

## The no-outside rule

A self-contained structure has no external vantage point from which a premise
could be imposed. The rule has three applications.

- **Justification.** A premise is justified when observers recover it and the
  specification that implements it runs and reproduces their measurements.
  The carrier, the observer algebra, the repair law and the closure
  coordinates P and N owe no further a-priori argument, and construction
  choices need no uniqueness theorem.
- **State selection.** No information may be assumed beyond the recovered
  constraints. The selected state is the least informative one relative to
  the recovered reference under every observer-visible constraint.
- **Transformations.** Every change of description that observers can make is
  implemented by the structure's own operations: each proper recharting of a
  carrier is realized from inside its complete reversible response.

## The agreement condition

Observers assign the same operational meaning to accepted data on shared
boundaries, naturally under every restriction, recharting and refinement. A
record is public when it survives every authorized comparison, and "our
Universe" is the world that public record describes. This is a definition and
an acceptance criterion. The [interpretation contract](INTERPRETATION_AND_TERMS.md#observer-measurements-and-operational-objectivity)
states what inter-observer agreement must preserve.

## The executable specification

The specification consists of the twelve-port carrier, the native repair
law, the observer algebra of states, preparations, measurements, retained
records and feedback, the declared services, and the closure inputs P and N.
One carrier supports the local observer, clock and geometry results; the
spherical federation, whose cross-carrier routing is a declared convention,
is the construction used for larger observer domains and refinement. The
specification is implemented and run in the
[model directory](../model/README.md). Its status is a candidate, judged by
the [acceptance criteria](../model/CONSTRUCTION.md#completion-target-and-acceptance-criteria)
of the completion target.

## Entry point and order

The exposition enters at the executable specification, because it is
concrete and because observers stand there, holding records and inferring.
The order of justification puts the hypothesis first, because it is what
exempts the specification from outside justification. The structure itself
has no first element in either sense, and there is no starting point in a
timeless loop.

## Theorem-library premises

The library papers state their theorems with three labelled premises. They
name parts of the basis above.

| Library label | Registry name | Part of the basis |
| --- | --- | --- |
| A1 | `observer_screen_architecture` | The carrier, federation and complete-response layer of the executable specification |
| A2 | `observer_meaning_agreement` | The agreement condition, with the transformation application of the rule as its endogeneity clause |
| A3 | `conditional_maximum_randomness` | The state-selection application of the rule |

The [premise reference](AXIOM_REFERENCE.md) gives the formal statements of
A1, A2 and A3, and the machine-readable source is
[the premise registry](../claims/axiom_registry.yaml). A theorem that uses them
states an implication; the [premise map](../model/PREMISES.md) says whether
the executable specification supplies its antecedents.

## The paper set

- The [flagship synthesis](../flagship/from_observer_consensus_to_standard_physics.tex)
  gives the account, develops selected results and states the closure
  hypothesis in full.
- The [executable observer paper](../paper/executable_observer_construction.tex)
  states the specification, its acceptance criteria, its certified components
  and the premise map.
- The theorem library holds the long proofs: consensus, observers and records, gauge structure, matter, screen microphysics, and spacetime with Einstein dynamics.
