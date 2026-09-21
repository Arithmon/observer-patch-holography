# Full-family routing contract

## Objective

Realize the complete q=13 and q=21 metric read menus on the captured W12
support, conditional on the supplied record, local feedback and routing law
in `specification.json`. Establish value/intervention preservation, induced
semantic order and explicit costs while distinguishing physical control
ancestry from the logical projection.

## Deliverables

| Obligation | Evidence |
| --- | --- |
| Local transport and finite history theorem | `SourceFeedbackTransport.lean` and `SourceReadRouting.lean`: exact six-operation hop, 6-event/7-read/7-write census, conditional full-history equality and both induced-order inclusions. |
| Complete production histories | Four retained q13/q21 run receipts and every phase checksum; `regenerate.py` executes and independently verifies every event against these commitments. |
| Independent semantics and provenance | `verify.py` reconstructs every metric decision and completed value; `verify.cpp` checks each event and consumed writer. `check_control.py` exhausts all logical pairs of four retained q3 controls. |
| Explicit costs | Preparation, transport, accumulation, commits, register lifetimes, routing-controller work and host allocation observations are separately counted. |
| Support comparison | `support_comparison.py` checks all L3/L4/L5 captures; full q13 baseline replay compares the diagnostic's read menu and every logical value. |

## Exit boundary

The result is a conditional compiler theorem with M1 retained as a supplied
structural rule. The specification records the separate source-derivation
obligation in the premise register. The result supplies neither an A1--A3
derivation of that law nor a universal impossibility theorem, physical clock,
quantum instrument or four-volume interpretation of routing primitives.

## Construction and proof

`capture_support.py` joins precisely the triangular cells sharing a vertex
and uses the public simulator's geometric port-assignment function. Those
public primitives reproduce every face and glued port pair of the L3
fixture before L4/L5 are captured. Each larger mesh has twelve degree-five
vertices, all other vertices degree six, and `6*C-30` glued seams. The sixty
unused port slots at pentagonal cells are retained. The unpublished
`oph_exact/support_wiring.py` is not a dependency. Geometric floating-point
assignment is captured as declared finite combinatorics; exact replay does
not treat those coordinates as a physical geometry or an interval proof.

Golden-coordinate Morton order assigns the q^3 source sites injectively to
the first q^3 host cells. This placement is supplied. Every host initially
has the twelve-port signed antipodal split of its prepared golden record;
six protected address scalars are also initialized. Repair can change live
loads while the prepared address remains protected, as in the bounded protected-address
construction. Unoccupied hosts are relays. Source-payload values are separate
immutable registers, initially `site+1`.

For each source version, an ascending-neighbour BFS tree is pruned to the
union of paths to its complete metric-recipient set. Each used edge performs
export, receiver reset, pair mean, factor-two capture, source reset, receiver
reset. All three feedback operations read only their owner carrier's state.
A carrier may access its own different ports through its protected local
memory; that local access is part of M1 and is charged as an export/capture,
not silently treated as an extra intra-carrier mean. Intermediate protected
versions branch to all needed children. A receiver accumulates each required
version once. Only after every read is complete are new immutable layer
versions committed. No route or read interface depends on payload values.

The exact real hop and sharp additive error bound are inherited from
`Lean/Geometry/SourceFeedbackTransport.lean`. The
`Lean/Geometry/SourceReadRouting.lean` additionally evaluates the concrete
six-operation word in integer half-units, proving source preservation,
capture, cleanup and its 6-event/7-read/7-write cost. This word discharges the
transport hypothesis in `six_event_compiledHistory_eq`, which proves
full-history value equality for every initial integer state and integer
local-commit intervention. It also proves equality of induced logical reachability
from the local-edge and lifted-read certificates. The native verifier checks
the order hypotheses against every consumed writer. The C++ verifier itself
is reviewed executable code, not a Lean-verified interpreter. This covers the complete
order without enumerating billions of logical pairs. A separate Python
traversal checks all 6,561 logical pairs in each of four q=3 controls.

The full q=13/q=21 `source` runs change the centre's initial payload by +1;
every completed value is checked against an independent exact recurrence.
The q=3 controls additionally change a later commit and the initial scratch
loads. Immutable registers are never overwritten in any replay. The
4,928-event `source_feedback_transport` receipt additionally executes old
version rereads after newer versions and intersecting scratch reuse; its
31-test suite remains part of validation. This is an inherited transport
control, not a q=21 old-version experiment.

The closed-register sum invariant gives a limited obstruction: no finite
composition of total-load-preserving operations can implement a reset that
changes that total. The Lean module checks this finite-word statement.
It explains why this compiler's reset is additional feedback. It does not
exclude tomography, other decoders, changed ancillary state, other laws, or
every operational realization of the metric order.

## Resource and custody model

For C hosts, n sites, K rounds, R logical reads and H multicast hops:

- Preparation uses `13*C+8*n` events: twelve scratch ports and one zero per
  host, one accumulator per site, six protected address scalars per site,
  and one initial payload version per site.
- Transport uses `6*H` events, `7*H` scalar reads and `7*H` scalar writes.
- Layer start/commit uses `2*K*n` events; accumulation uses R events and
  `2*R` scalar reads. Each start reads its own zero, and each commit its own
  accumulator. All writes and read-from edges are counted.
- Protected storage is `C+6*n+(K+1)*n+H` scalar registers; mutable storage is
  `12*C+n`. No relay archive is discarded from that accounting.
- A requested depth-d read traverses `6*d` transport events plus its local
  accumulation; shared multicast prefixes are charged once in H. Total
  serial duration in the declared unit-event schedule equals the total
  event count. It is not the distinguished layer count K or a physical clock.

Each scalar in the retained executions has an exact signed 64-bit half-unit
representation. Independent replay checks mean exactness, arithmetic bounds
and every result; all supported runs fit. The arbitrary-intervention theorem
uses unbounded integers: it does not promise that arbitrary interventions or
inputs fit the native binary. Each explicit
event row occupies 64 bytes, including opcode, owner, input/output register
identifiers, both consumed writer identifiers and the result. Two-output
mean writes share the recorded result. Absent reads have a reserved sentinel;
event IDs are row offsets. Register versions are their checked successive
writer IDs, and immutable cells have only one writer.

Input bytes include the complete read menu, port wiring, host mapping and
golden addresses. `implementation_resources` records actual native vector
allocation requests (including simultaneous old/new register allocations),
auxiliary input/BFS buffers, the 16-byte wide arithmetic temporary, the
64-byte event row, and measured native peak RSS including runtime/allocator
overhead. These are host-software measurements, not a physical patch-memory
certificate. The producer's register structure is 24 bytes; the verifier's
extra origin certificate makes its structure 32 bytes. Codec work is in
65,536-event blocks, with a 4 MiB raw block and at most three differential
dependency levels for the archived controls. Stored manifests and compressed
part byte sizes account for the archive metadata too; compression is not a
reduction in executed events or scalar accesses. The raw-block size is not
the codec's total memory use; the RSS observation covers the native producer,
not the whole Python preparation/codec pipeline.

Routing computation is also charged separately. The implementation
recomputes each source BFS in each round: `K*n*C` vertex dequeues,
`2*K*n*(6*C-30)` directed-edge examinations, `K*n*(C-1)` discoveries and
tree-mark examinations, and `2*K*n*C` parent/mark initializations. Pruning
traverses H parent links and makes R recipient-path queries. These control
counts are recorded in `routing_control_work`; at q=21 there are
11,377,138,500 directed-edge examinations. No optimal compiler-cost claim is
made. Metric-menu preparation is a separate declared setup calculation;
the independent input verifier checks all `n*n` metric decisions. These
control operations and host storage/compression work are not assigned a
physical clock by the unit-event repair schedule.

The complete production event streams are deterministically regenerated and
independently replayed. The retained run receipts bind the input bytes,
producer source, whole decoded stream, event census and resource census.
The retained segment manifests bind each decoded phase. Regeneration checks
all these commitments without modifying them. A producer process emits every
64-byte event through a bounded pipe into the independent native verifier;
there is no event sampling or summary-only acceptance.

The q=3 controls retain their losslessly compressed event files. For optional
production materialization, later layers store field-by-field uint64
differences from the first layer; intervention phases store differences from
the corresponding baseline phase. Addition modulo 2^64 recovers each field.
The retained segment manifests describe those compressed parts, including
their byte sizes and checksums. Production compressed parts are not stored
in Git, and regeneration authenticates decoded events rather than claiming
to retrieve the original compressed bytes. No external custody, signature
or laboratory attestation is inferred.

Semantic read-from order, serial control order, resource/write hazards and
archive custody remain distinct. Constant resets do not consume the previous
scratch payload, but they remain physical resource operations in the tape.
Routing events are not assigned four-volume. Physical clock units,
header/address transmission, source selection, and quantum copying or
instruments are outside this finite classical implementation.

The two production levels are separate finite regulators. No cross-regulator
state/provenance refinement or uniform physical clock/capacity limit is
claimed. Glued-carrier distance counts inter-carrier means with the supplied
local archive-to-port access; it is not the old port-graph distance when
intra-carrier motion is restricted to means alone.

## Support comparison

After authenticating every q13 baseline primitive and its consumed writers,
`verify.py` additionally compares the entire metric read menu and all **10,985
logical preparation/commit values** with the support-wiring diagnostic's separately executed q13 tape.
That tape has **1,176,764 reads**. Thus the checked logical projection has the
same order intervals as the support-wiring diagnostic's record-metric comparison. This does not assert
equality with the support-wiring diagnostic's canonical-only W12 trajectory or with the full order of
this compiler's **19,113,548 primitive events**.

The diagnostic also limits the interpretation of the contrast:

* Its one-level W12 intervals contain only 9 and 75 events, with MM estimates
  1.18605 and 2.24830. They do not establish a stable 2+1 dimension. The two
  L3-to-L5 intervals give 3.38441 and 3.52077 under the declared selections.
* The q13 lag-three record-metric interval is interior and gives 2.82788.
  Lag four gives 4.14933 but is **clipped**. Routing that logical order does
  not turn its lag-four estimate into an interior result.
* The q13 one- and two-dimensional scientific controls belong to the the support-wiring diagnostic
  diagnostic; their lag-four rows are also clipped. These controls are not routed here. The routed intervention controls remain the declared q=3
  baseline/source/branch/scratch executions.

The shared graph therefore does not make the operation laws identical. the support-wiring diagnostic
executes canonical means on every intra-carrier and glued seam under its fixed
sweep schedule. This compiler adds protected memory, export, reset, capture,
accumulation and control under the retained M1 premise. It constructs the
record-metric read law **conditional on those additions**; it does not derive
that law from the canonical-only dynamics or identify a physical dimensional
transition. The comparison does not supply a premise of the routing theorem.

| Quantity, per baseline or source-intervention run | q=13 | q=21 |
| --- | ---: | ---: |
| Source sites | 2,197 | 9,261 |
| W12 host carriers | 5,120 (L4) | 20,480 (L5) |
| Declared rounds | 4 | 5 |
| Logical reads, including self | 1,176,764 | 14,559,225 |
| Executed multicast seam hops | 2,972,512 | 46,288,315 |
| Retained events, including preparation/accumulation/commit | 19,113,548 | 292,722,053 |
| Scalar register reads | 23,178,688 | 353,229,265 |
| Scalar register writes | 22,086,060 | 339,010,368 |
| Protected registers, including every retained relay | 3,001,799 | 46,419,927 |
| Mutable registers, including accumulators | 63,637 | 255,021 |
| Maximum shortest glued-carrier distance of a requested read | 44 | 91 |
| Sum of requested read distances, including all rounds | 12,599,152 | 294,030,780 |
