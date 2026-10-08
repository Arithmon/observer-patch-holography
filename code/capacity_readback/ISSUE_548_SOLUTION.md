# Issue #548: Source-Derived Public-Checkpoint Capacity Packet

## Result

The first fixed-cutoff physical packet is the oriented edge-center register

```text
X_reach = P_12 x {write, check},
D = |X_reach| = 24.
```

The word *physical* is used in the repository's typed fixed-cutoff sense: the record atoms and checkpoint laws are generated from the declared finite screen carrier rather than arbitrary `record_i` labels. The artifact does not claim a laboratory realization, a capacity-indexed cosmic family, a unique finite-size slack zero, horizon saturation, or an electroweak load identification.

The exact certificate emits

```text
|Omega_tilde(r,D)| = 1,
|X_pub(q)| = |X_reach(q)| = 24,
|K(q)| = 40,
|<support(K)>| = 40,
E(G_q) = empty,
M_0(q) = alpha(G_q) = 24 = D.
```

## Acceptance mapping

| Acceptance item | Executable receipt |
|---|---|
| Frozen carrier; unclosed trials; output-blind membership | `build_terminal_fiber_manifest`: exact world plus all 30 single-edge deletions, 24 single-slot deletions, and 12 inverse fixed-point faults. `is_terminal_world` reads only structural carrier fields. |
| Complete terminal-fiber manifest | 67 fully materialized declared trials, SHA-256 constructor/candidate/manifest receipts, exactly one terminal ID, and `terminal_fiber_complete=true`. |
| Observer/interface atoms and total readouts | 12 observers, 24 local atoms each, 30 interfaces, 24 interface atoms each, total endpoint readout maps. |
| Endogenous histories and preregistered publicness | One semantic propagation/check/commit history per public section; universal twelve-port policy frozen before evaluation. |
| Complete joint kernels and support compositions | Each of the 40 named `D5 x C2_antipodal x C2_orientation` actions is checked on all 24 actual public sections. Every entry of the supplied 40-by-40 table is replayed as left-after-right on those sections: 38,400 pointwise composition identities. |
| Local-marginal checks | Complete continuation, observer and source domains are required. Marginals are derived from every normalized global row and compared exactly against original local probabilities: 11,520 row checks, without float coercion. |
| Injective reversible generators and exact model count | Every continuation is deterministic, injective, and surjective; the CSP backend counts exactly 24 public sections. |
| Compound graph, MIS, exact decoders | Empty 24-vertex graph, complete 24-vertex independent set, and inverse decoder for every continuation. |
| Approximate branch and TV robustness | Exact worst-input success 1; for rowwise TV distance `delta`, the same decoder gives success at least `1-delta`, hence `M_delta=24` by the carrier upper bound. |
| Carrier bound and rank-one saturation | Sparse projections `P_x=|i><i|`, 276 pairwise orthogonality checks, rank sum 24, identity resolution. |
| Empty/incomplete/ambiguous/singleton fibers | Explicit `classify_terminal_fiber` controls for all four cases. |
| Required controls | Isomorphic relabeling; cyclic permutation; same-marginal/different-joint coupling; tiny full-support noise; circular-definition rejection; target taint; identity family; erasure family; equal-finite-suffix nonpromotion. |
| Extension and refinement injections | Separate 24-to-48 capacity extension and fixed-24 refinement embeddings, both checked by `no_new_confusability`, plus negative controls that deliberately add a new edge and are rejected. |

## Why the CSP change is required

The previous global-section evaluator enumerated the Cartesian product of local atom sets. The source carrier has 24 local atoms at each of 12 observers, so that strategy starts from `24^12` candidates. `public_record_csp.py` performs the same exact enumeration with early interface propagation and a most-constrained-variable order. The existing finite tests verify extensional equality with Cartesian enumeration on a generic noninjective record diagram.

## Reproduction

```bash
cd code/capacity_readback
python3 source_derived_public_checkpoint_packet.py --output-dir runtime
python3 -m pytest -q
```

The generated JSON files are:

```text
runtime/source_derived_terminal_fiber_manifest.json
runtime/source_derived_public_checkpoint_packet.json
runtime/source_derived_public_checkpoint_certificate.json
```

## Source-binding audit

The original checker regenerated the canonical multiplication table and
compared it with the supplied table, without composing the supplied kernels.
Swapping the identity and rotation kernels, together with their local
marginals, retained `PASS` despite 4,620 failed pointwise composition
identities. Inverting every kernel retained identity, permutation closure
and capacity 24, yet failed 19,200 of the 38,400 identities.

Composition is also insufficient by itself. Replacing every named rotation
by its inverse through a group automorphism preserves the multiplication
table, but sends `upper_0/write` to `upper_4/write` under the named `r1`
operation, instead of the source's `upper_1/write`. Certification therefore
checks both the actual composition law and the named source action.
The record-to-slot map comes from the supplied section atoms and interfaces;
an alias declaration cannot substitute a different source meaning.
The actual observer domain and each authorized read set must also realize
the universal twelve-port source policy. A one-observer policy cannot retain
that source certificate merely because its continuation permutations agree.

The retained tests first reproduce false acceptance on main `a20d4736`,
then check the repaired public certificate with recomputed packet hashes.
The independent scalar oracle checks all 960 named input/output mappings;
separate positive controls preserve valid presentation changes. The existing
near-unit normalization policy for global relative row weights is unchanged.
Local marginals are compared to those normalized probabilities exactly.
Zero entries may be omitted or placed anywhere in a row. The downstream
noise control now decodes the positive support, rather than the first map
key: a leading zero previously produced a success probability of zero while
the control still reported `PASS`. Its returned status now checks its noise
and decoding identities, and nonreversible source kernels are refused.

The follow-up audit challenged the other evidence carried by the same source
certificate. At PR head `67f9a9f`, a valid hash still admitted an empty actual
terminal fiber, reversed or target-reading histories, a changed interface
topology, missing sparse projections, and an identity-only generating set.
The generated false certificates also survived serialization into the direct-N
consumer's fixed-cutoff result. The universe-level result stayed unavailable;
this was a false finite-source claim, not a new cosmological result.

The checker now reconstructs the 67 declared source faults from the ports,
edges and oriented slots, validates each supplied candidate and recomputes
membership. Neither the trial enumerator nor candidate materializer supplies
the verifier's expected census or candidate. Histories are replayed from birth
and seed through eleven valid propagations, followed by public audit and commit;
every propagation must reach a new child from an already reached parent along
a source edge. The supplied witness, semantic metadata and target-use flags
must agree. Alternative valid spanning trees and trial enumeration orders pass.

The actual thirty interfaces must implement the declared slot readouts, with
complete matching atom and section inventories. Every sparse carrier projector
is checked against its declared basis and support. Coherent basis changes,
interface atom renamings and reversed edge presentation remain valid. The
declared generators must generate all forty supplied operations. Fixed-source
identity and recipe fields are checked; free-form input prose cannot change the
certificate's scientific boundary. Refused evidence remains refused after JSON
serialization into the direct-N consumer.

This is a repair of verification for the fixed finite source. The canonical
packet, terminal manifest and certificate retain their existing outputs;
downstream capacity and physical classifications do not change. It supplies
no new source-selection law or physical realization.

Validation uses the complete capacity directory on Windows/Python 3.13 and
Linux/Python 3.12 with warnings treated as errors. Isolated faulty variants
exercise removed source binding, table replay, domain validation, exact marginal
comparison, diagram and carrier checks, administrative flags, restored first-key
decoding, unchecked publicness, constant success and blanket refusal. Two mutants
initially survived: omitted table-domain checking and omitted local output-domain
checking. Retained public-path controls now reject both, including a vacuous
empty table and an unknown zero-weight read atom. Independent downstream lift
replay and byte-exact canonical receipt regeneration remain required.
The existing capacity workflow executes this entire directory on both
operating systems; these bounded checks are not a formal verification claim.

## Remaining theorem boundary

This issue closes the first source-derived reversible packet at frozen `D=24`. It does not supply the source family `D -> packet(D)` or the exact slack law selecting one regulator-stable cosmic dimension. Consequently it does not by itself establish the universe-level equation `N=log M_0(U_N)` at a selected `N`.
