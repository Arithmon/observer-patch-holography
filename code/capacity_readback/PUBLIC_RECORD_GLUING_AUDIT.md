# Preserve the exact public record set

This repair under issue #1033 corrects the finite construction of public
record sections before checkpoint capacity is evaluated. It implements the
existing atom-gluing equation in Pro5. It supplies no physical source,
capacity selector, temperature or clock.

## Reproduced failures

On main `7f99c8cc`, the optimized section enumerator can disagree with direct
Cartesian enumeration. With one observer, atoms `0,1`, and a self-interface
whose two readouts are respectively `(x,y)` and `(y,x)`, neither atom is
compatible. The optimized enumerator nevertheless returns both. Appending
the analogous contradictory self-interface to an otherwise valid two-record
packet still yields `PASS` and saturated capacity two.

The search tested an interface only when its other endpoint was already
assigned. A self-interface's endpoint was never assigned at that test, so
the constraint was omitted entirely. Rejecting every self-interface would
also be wrong: an interface with readouts `(x,y)` and `(x,x)` admits exactly
atom `0`, and equal readouts admit every atom.

A separate defect at the same construction stage merged distinct sections:

```text
{"a": "x|b=y", "b": "z"}
{"a": "x",     "b": "y|b=z"}
```

Both formerly became `a=x|b=y|b=z`. On the unconstrained product of these
two two-atom alphabets there are four global records but only three old
identifiers. A dictionary keyed by those identifiers loses a real record.
The repaired pipeline preserves all four through an identity-checkpoint
capacity calculation.

Further regressions cover readouts outside declared interface alphabets,
non-string/coerced labels, string values used as entire atom alphabets,
recursion failure on a 1,200-observer singleton diagram, and exceptions
instead of packet failure classifications for invalid reachability IDs.
All seventeen cases in `test_public_record_gluing_regressions.py` fail on
the baseline before the repair.

## The gluing law and the solver

For local atom sets `X_O` and readouts `r_Le`, `r_Re`, the public sections are

```text
{(x_O) in product_O X_O : r_Le(x_L) = r_Re(x_R) for every interface e}.
```

For distinct endpoints this is a fiber-product condition. For a
self-interface it is the equalizer `r_Le(x) = r_Re(x)` on **one** atom.
Introducing two independently selected copies of that atom would implement
a different condition. Noninjective readouts are permitted: two distinct
local atoms with the same readout stay distinct local atoms.

`public_record_csp.py` validates the complete diagram before solving it:

- Observer labels, local atoms and interface atoms are nonempty strings;
  each alphabet is a finite nonempty sequence without duplicates.
- Every interface has known endpoints and total readouts. When an
  `interface_atoms` alphabet is declared, both readout images must lie in it.
  Older declarations without that field use equality of the supplied image
  labels; they do not assert an additional interface alphabet.
- Parallel and self-interfaces are all conjunctive constraints. A malformed
  later interface cannot be hidden by an inconsistent earlier one.

Self-interfaces first restrict their own domains. Directed binary interfaces
then remove atoms having no matching readout in the current neighboring
domain. A queue propagates every domain reduction to affected interfaces.
Remaining multivalued domains are partitioned by an explicit stack.

The correctness argument is finite and exact:

1. A global section cannot contain an atom rejected by a self-interface.
2. If an atom has no compatible atom across one binary interface, no global
   section extending the current domains contains it. Repeated propagation
   therefore preserves every global extension.
3. Choosing one atom of a remaining multivalued domain partitions all global
   extensions into disjoint cases. The search explores all surviving cases.
4. At a leaf every domain is a singleton. Stable support across each binary
   interface then implies its readout equality, and the self-interface
   restrictions still hold. The emitted section is globally compatible.

Together these prove soundness, completeness and absence of duplicates.
Arc consistency alone is insufficient: an odd cycle of binary inequality
constraints has two solutions on each edge but no global section. The
retained control ensures the solver branches and rejects that obstruction.

No Python recursion is used. In a connected 1,200-observer equality chain,
one binary choice propagates through the entire chain; the complete two-record
set needs three allocated search nodes. The retained twelve-observer,
twenty-four-atom source diagram needs twenty-five nodes. These are explicit
budget controls, not timing-based benchmarks or polynomial-time claims for
arbitrary finite constraint systems.

The default limits are 100,000 allocated search nodes and 100,000 output
sections. Callers of the solver may set either explicitly. Exhaustion raises
an exception instead of returning a partial list or certifying nonexistence.
The packet evaluator returns its existing failure status with the refusal
reason, never a capacity result. Output list order is not a mathematical
contract; record identity and membership are.

## Record identity

The common encoder escapes `%`, `|` and `=` inside each observer/atom label,
then joins observer-sorted pairs using the historical separators. Escaping
`%` first distinguishes a literal `%7C` from an encoded `|`. The decoder
splits only on unescaped separators, reverses exactly those three escapes,
and checks canonical re-encoding. Duplicate keys, alternate escape
spellings, empty labels and out-of-order observer pairs are rejected.

This decoder is a left inverse of the encoder, proving injectivity. Every
existing identifier whose component labels avoid the three reserved
characters keeps its exact bytes. Labels containing those characters must
be re-encoded from their explicit section mappings; ambiguous historical
strings cannot determine which record was intended and are not guessed.
The reversible packet producer now uses the shared encoder as well.

Reachability is checked against the encoded compatible section set. Unknown
or malformed identifiers produce `NO_PUBLIC_RECORD_REACHABILITY` through the
packet evaluator, including when a newly enforced interface removes every
previously declared record.

## Independent completeness evidence

`verify_public_record_sections.py` imports neither the solver nor the packet
producer. It independently validates the diagram and checks the **complete**
section set, including omissions and duplicates, by two separate procedures:

- If the Cartesian candidate count is within its budget, it enumerates all
  assignments and checks the original readout equations directly.
- For a larger diagram with injective binary readouts, it chooses one root
  per connected component. Each root atom determines at most one assignment
  throughout that component by inverse readouts. Every possible root tuple
  is tried, and all interfaces, including self-interfaces, are checked.

Every global section has one of the enumerated root tuples and is uniquely
determined by it, so the second procedure is also a completeness check.
It independently replays all twenty-four records of the retained source
diagram with twenty-four candidates, rather than enumerating `24^12` tuples.
For larger noninjective diagrams or too many root tuples it raises
`ReplayBudgetExceeded`; it does not claim that their record sets are empty.

The independent test controls include:

- All 4,096 binary-readout triangles and all 512 pairs of parallel or
  self-interfaces, compared with direct Cartesian enumeration.
- A further 5,184 injective diagrams with unequal local alphabet sizes and
  a disconnected component. All two-edge chain maps are exhausted, with
  separate cycle, parallel-edge and self-interface controls. The verifier's
  four-root-candidate budget forces reconstruction instead of its 24-candidate
  Cartesian path; direct enumeration independently checks completeness,
  omitted records and incompatible extras.
- Relabeling/orientation invariance and 343 identifier round-trip cases,
  including literal escape strings, separators and Unicode labels.
- Nine forged-list controls, including constant `PASS`, missing records,
  incompatible extra records, duplicates, malformed maps and wrong codomains.
- The full retained source diagram replayed under `python -O` with the
  producer absent, accepting its complete set and rejecting an omitted record.
- Seven isolated implementation mutations: omitted self constraints,
  omitted propagation, ignored codomains, colliding identifiers, truncated
  output, membership-only verification and omitted replay self constraints.
  Each mutation must cause an ordinary assertion or missing-exception failure;
  import and collection failures do not count as detection.

## Downstream and provenance impact

The section solver feeds the exact capacity evaluator, reversible reference
packet and source-derived packet's reachability, marginal, refinement and
decoder checks. Contradictory self-interfaces now remove impossible records;
reserved label characters no longer merge distinct records. Well-formed
ordinary labels and compatible diagrams keep the same record set.

The three retained source artifacts were regenerated separately using main's
implementation and the repaired implementation. Both regenerate exactly the
tracked bytes, with these SHA-256 values:

```text
packet      3604147f715327abdc96ad56160fa5ff5384af767d0a57c49beae5b92f3ed96a
manifest    34e22b6c554c657ca0573516a8c543c2bb995e581cd072c5592f3a4b36c94f69
certificate bf48d67735360ca3d54d4c2b58bb14f0884bb917cc0f4d8d66e7229a75a3ff93
```

The existing Pro5 atom-gluing equation already imposes every interface
equality. No paper, public claim, registry payload or frozen evidence needs
changing. In particular, the frozen readback specification and mandatory
runner remain unchanged. This PR is independent of #1054's channel-arithmetic
repair and addresses an earlier stage in the same finite record pipeline.
Their shared capacity workflow is byte-identical so the branches can combine
without competing test configurations.

Reproduce the entire affected suite, including independent and mutation
controls, on either platform with:

```sh
python -m pytest -q code/capacity_readback -W error
```
