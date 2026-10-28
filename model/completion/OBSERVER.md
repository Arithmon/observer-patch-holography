# Continuing C/E observer with accounted record export

This is a bounded CP-01 component under the [campaign contract](campaign.json).
It extends the existing two-event observer to an indefinitely iterable ideal
instrument with a finite classical FIFO. Exact retained executions cover six
model attempts for C, E, and the declared delayed C/E mixture. A model attempt
is the discrete operation parameter; no physical duration is assigned.

The result does **not** accept CP-01: preparation, random selection, commitment,
record storage/export and feedback are supplied services in this experiment.
It does not source-compile those services or realize remote record replication.
The original q=3 common-code preparation/history remains unexecuted and is not
used or discharged here.

## Actual instrument and independent compatibility contract

The source is the existing total-two space `H_78`. An arbitrary untouched
reference `R` is allowed. The classical memory is a FIFO of at most the
campaign's declared capacity; each record includes its hold/repair tag,
exploration/feedback selection tag, seam, C/E law, coin, native record label
and derived current. The record current is a decoder of the selected C/E
instrument, not an additional measurement of the loads.

The independently specified source obligations are nonnegative integer loads,
conserved total two, unchanged spectator ports and, on an actual repair,
floor/ceiling agreement at its endpoints. An identity hold preserves loads
without claiming to repair a disagreement. Record compatibility requires
exact agreement on the event-labelled restriction of a common committed tape
to shared seams. The implementation verifies this shared-tape property; it
does not use it as evidence for separately constructed communicating observers.

The controller proposes `(0,5)` when the last retained current is positive,
and `(0,1)` otherwise. It never reads the source configuration. Each repair
chooses the proposed seam with weight `1-epsilon` or an independently uniform
native seam with total weight `epsilon`. A declared intervention can change
the proposal while retaining this exploration floor. The C/E mixture tag is
committed, so it is a classical mixture of instruments, not a coherent sum.
Both law endpoints and all semantic zero branches remain represented.

For a branch `b`, the source Kraus operator is

```text
K_b = sqrt(w_b) M_b,
w_b = repair_probability * selection_probability * law_probability / 2.
```

`M_b` is the existing native partial injection selected by its coin and C/E
record. Its untouched-reference action is `K_b tensor I_R`. The hold branch
has operator `sqrt(1-repair_probability) I`. Committed memory is classical:
the corresponding source/memory algebra is a direct sum over memory words.
No preservation of superpositions of different controller memories is claimed.

For each allowed classical memory input, every `M_b` is a partial injection
and the sum of its weighted diagonal effects is exactly `I_78`. Consequently
the complete instrument is trace preserving on every complex source matrix
and on arbitrary reference extensions. This reuses the existing native
instrument certificate; the new work is its explicit composition with the
declared scheduler and bounded-memory/export law. The tests reconstruct the
native branches independently from floor/ceiling endpoint arithmetic.

Induction supplies the complete finite-history channel: append each branch
record, update its classical memory, and compose its partial injections.
`history()` computes any recorded word's unnormalized output and probability,
including its exported records. Zero branches remain zero; there is no
postselection renormalization. This kernel specifies the whole exponential
history tree without storing it explicitly.

## Finite memory and continuing operation

With `export_oldest`, the oldest FIFO record is transferred to a declared
environmental output whenever a new record replaces it. A hold is itself
recorded and has current zero. Capacity-one operation therefore retains one
record and exports `n-1` records after `n` attempts. Concatenating exports and
the remaining FIFO exactly reconstructs the committed history.

This has bounded working memory and unbounded total history storage as
`n` increases. Export is not free erasure, an implicit reset, or a proof of
closed finite-resource continuation. The report explicitly counts exports.
With `halt`, a full memory instead emits the timeout outcome, preserves its
records and source, and enters an absorbing halt control. It provides a
finite-storage control against the exported-history continuation.

The efficient exact runner executes the marginal on source, reference and
the last complete record for a one-slot export policy. Earlier records are
traced out in that **particular marginal**. The full recorded channel remains
available through the transition kernel and prefix replay. Grouping old
memories by their proposed seam is exact because the old slot is exported
and all future control depends only on that proposal. The two-step tests
compare this grouping against explicit recorded paths.

## Retained checks and limits

- Exact six-attempt executions preserve total probability and the
  unconditional reference marginal for C, E and the delayed mixture.
  Native current conditioning changes the reference branch state, which is
  retained rather than discarded as an irrelevant outcome.
- A coherent source/reference matrix unit survives a C branch and is killed
  by the corresponding E contrast commitment. The two laws remain distinct.
- Exact replay of a six-hold word has probability `1/4096`, five exported
  records and one retained record. Continuing operation exceeds the original
  two-event tape without hiding its storage cost.
- On the uniform binary preparation, the two-attempt adaptive scheduler
  changes 30 source probabilities, with maximum deviation `1/960` from
  `1/66`. This independently reproduces the
  [adaptive-reference countercontrol](adaptive_reference.py): the old uniform
  source reference is not automatically stationary for the joint controller
  process. The separate exact source/bit stationary construction is lifted
  into this actual complete-record observer below.
- Shared-tape regional restrictions agree on overlap event labels. Native
  acquisition, transmission, remote storage and independent observer clocks
  remain unconstructed here.

The record-selected source operators, full instruments and arbitrary-reference
extension reuse the existing native observer result. This extension adds a
specified continuing feedback/memory composition, its complete history
semantics, accounted overwrite exports, timeout control and exact executions.
It does not derive fields, topological matter, physical energy or clocks.

## Exact stationary reference on the actual retained record

`observer_reference.py` takes the exact stationary source/current-sign-bit
law from `stationary_reference.py` and applies **this observer's complete
record instrument** once. The output is a diagonal density on the source and
the full one-slot FIFO, including seam, selector, law, coin and record labels.
The construction is performed separately for the C and E endpoints with
repair probability one and exploration probability one half.

Projecting the new record to `1[current>0]` gives exactly the original solved
source/bit law. Applying the actual full-record observer again leaves the
entire new source/record density unchanged, with exact rational residual zero.
Thus the bit quotient is a genuine intermediate construction; its two bit
labels are not silently identified with the much larger complete record.

The full lifted reference is retained by a canonical hash, exact per-record
probabilities and support dimensions, and the deterministic rational lift
from the retained source/bit solution. Positive support is checked in every
record block. Faithfulness holds on the direct sum of these explicitly
compressed support algebras, not on full `M_78` or on unrestricted full-record
blocks. Access to every operator in these supports is not supplied by the
mathematical compression.

This stationary object is a finite source/last-record marginal. Its exported
history continues to grow. Preparing this reference, implementing its modular
action, comparing local clocks and fixing physical energy remain CP-02
obligations. Holds in this observer overwrite the last current with zero;
the delayed member therefore needs a separate stationary construction and
does not inherit the repair-only endpoint lift.

Run from the public repository root:

```sh
python3 -m model.completion.observer --check model/completion/observer_results.json
python3 -m pytest -q -W error model/completion/test_observer.py
python3 -m model.completion.observer_reference --check
python3 -m pytest -q -W error model/completion/test_observer_reference.py
```

The receipt binds the frozen campaign and the executed source files. Refresh
it only after a deliberate source change; it is a new component result, not a
rewrite of an earlier retained packet.
