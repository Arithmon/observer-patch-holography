# Exchange history, sharp sectors and active memory

This CP-03/#1139 diagnostic supplies a live **relational** sector and an
optimal actively corrected qubit inside one such sector. It also gives two
exact boundaries: every sharp sector differs from the declared reference on
the full exterior, and every closed odd exchange under uninterrupted native
C/E repair has a rank-one source input fiber. A localized physical particle,
fermionic exchange phase, rest energy, inertia and gravitational stress remain
unestablished. CP-03 and CP-05 remain open.

The carrier is the existing binary total-two shell, with 66 source states.
The classical core is an unordered pair of ports, stored in seven bits.
For every committed nonzero native current, exchange the seam endpoints in
both source and core. Then `q(x,c)=sum(x[p] for p in c)` is exactly preserved.
Core update, complete coin/current or contrast records, and any further core
reset are explicit parts of the state. This is neither a bare-source reducing
projection nor a new unrecorded phase operation. The core need not remain an
edge or geometrically small cap after transport.

## A sharp sector has a measurable exterior defect

The transported ideal reference assigns per-state weights
`(2/225, 1/36, 2/45)` according to q. It is the declared faithful extension,
not an assertion that a finite physical preparation has exactly that state.
For a fixed core the sector dimensions are `(45,20,1)` and their reference
probabilities are `(2/5,5/9,2/45)`.

Every density matrix supported in q has exterior particle number `2-q`.
Measuring the projection onto that number gives probability one, whereas the
reference gives `p_q`. Therefore the distance of its restriction to the
**full ten-port exterior algebra** from the reference is at least `1-p_q`:

| q | Source/core sector dimension | Minimum full-exterior distance |
| --- | ---: | ---: |
| 0 | 2970 | 3/5 |
| 1 | 1320 | 4/9 |
| 2 | 66 | 43/45 |

The reference conditioned on q attains each bound. This proves the exact
minimum, including arbitrary coherent source states in a sharp sector. It
does not bound every smaller cap, exclude screened/reference-dependent
particle constructions, or exclude superpositions/mixtures of these sectors.
For a reference known only within trace distance epsilon, subtract epsilon
from the bound. Changing the reference is a changed comparison contract.

For the uniform-core ensemble conditioned on q, the source marginal is
`I_66/66` for all three q. Thus no source-only channel with core-independent
controls can distinguish these ensembles, even with a complete native tape.
The joint core and terminal source computational read can distinguish q.
The implemented read retains all 512 source addresses, including padding and
arbitrary reference vectors. It is a supplied destructive read service,
already used by the source observer; it is not an additional native seam
measurement or a quantum nondemolition read.

## A closed odd native exchange resolves its input

The existing [configuration-complex calculation](../clock_defect_energy/DEFECT_HISTORY.md)
establishes `H_1(Z)=Z/2` for its declared collision-free local cells. Its odd
class is a label derived from actual ordered particle paths. The new result
concerns the **complete quantum instrument** carrying those paths.

Take any two distinct binary total-two configurations compatible with the
same complete C/E history. A shared signed current translates both loads
identically, so their difference is constant. Wherever the difference is
nonzero, binary occupancy forces each load to stay fixed. These are their
distinct spectators. A common nonzero move can move only a particle occupied
in both configurations, into a site empty in both. They share at most one
particle. Thus if either trajectory closes, its distinguished spectator never
moved and the other particle returns to itself: the exchange is even.

Consequently, if a positive history closes with odd exchange for one input,
it has **only that input**. Its Kraus map is
`K_h=sqrt(w_h)|y><x|`, of rank one. The local induction premises are checked
against every raw C/E map; the proof applies to any finite record-adaptive
native schedule, with holds and stopping allowed and no intervening quantum
controls. Classical random choices must be retained or conditioned upon.
The implementation also checks all triangular exchange orientations.

For the three-seam word `(0,5),(0,11),(5,11)`, each law has one odd closed
flag, with initial pair `{0,11}` and probability `1/8` on that input. All
38 C or 56 E nonzero blocks are retained, with 512 or 4096 total raw slots.
The positive rank-one flag is entanglement breaking, even with an untouched
reference. Multiplying its committed Kraus operator by a scalar minus sign
changes no output channel.

Flag-aware exact recovery of an unknown k-dimensional code requires
`V* K_h* K_h V = p_h I_k`. For k greater than one, a rank-one compression can
meet this only with `p_h=0`. Therefore an exactly recoverable code for the
whole uninterrupted native history has zero probability for every such odd
flag. This applies to a final decoder; it does **not** rule out intermediate
active recovery, coherent extra order registers, different record commitment,
or refinement. Those change the tested grammar and need their own account.

## The optimal q=1 code is an active memory control

The existing [flagged-code proof](../common_source/ENCODED_MEMORY.md)
establishes that any exact all-seam C/E correctable code compresses every
port occupation to a scalar. The actual moving effects are
`n_a(1-n_b)` and `n_b(1-n_a)`; their difference is `n_a-n_b`.
Connectivity and total occupation two then give scalar compression of all
`n_a`. On sharp q=1, every exterior `n_j` has rank two and their sum is the
identity. At least one scalar compression is nonzero, so **k is at most two**.

The existing source-executed code attains this bound for core `{0,1}`:

```
|0_L> = (|{0,2}> + |{1,3}>)/sqrt(2)
|1_L> = (|{0,3}> + |{1,2}>)/sqrt(2).
```

Its complete scalar branch effects are verified on all thirty seams for both
laws. This packet executes the original source encoder, the routed `(0,5)`
repair, every record-selected full-register recovery, and two separate
terminal read settings. It does not declare new all-seam source compilation.
Before recovery, native source plus moved core has q=1 exactly. Recovery
returns the source to the fixed code, while the moved core can be `{1,5}`.
Keeping that old core gives aggregate q=1 probability **7/8**, for either law
and any logical input. Explicitly resetting the core to `{0,1}` restores q=1.
Its value changes with probability **1/4**.

The reset appends one status/old-core/new-core record (15 bits) per non-abort
attempt; it does not overwrite the native record. It is a counted classical
control, not a consequence of native conservation. Recovery needs no fresh
quantum bits. The initial unknown logical qubit and eight blank code bits,
native coin/record preparations, all source pulses and proper-code load/export
events, core storage/reset records, and terminal read services are declared.
Append-only export capacity and physical service accuracy are not supplied
by this finite control.

The inverse encoder and one source Y rotation distinguish logical plus/minus
with probabilities one/zero, retaining every complementary address. The
boundary computational read and logical read are separate terminal
experiments. A boundary basis read destroys the logical coherence; no claim
of simultaneous preserved quantum memory follows from coarse-graining its
classical output. All source execution comparisons are bounded numerical
diagnostics against exact ideal map identities, not hardware error bounds.

This is the smallest same-object interface presently supplied to CP-04:
an identified code, sharp relational charge, native branch probabilities,
recorded recovery/core resets, and explicit recovery work counts. Pulse
counts are not rest energy. Under uniform seam selection, the complete native
moving-current probability on this code is **1/6 for every logical state**,
for both C and E. The uniform-source background is `5/33`, giving excess
**1/66 per native attempted repair**. The fixed-core ideal reference instead
has activity `341/2250`, giving excess `17/1125`. These are alternative
declared baselines; neither statistic includes the coherent recovery work.
The all-seam activity is an exact native effect calculation, while the new
source execution above uses the single specified seam.

Relating this object's excess maintained action
to independently measured dispersion, clock response and stress is still
required; neither uniform nor transported-reference averaging performs that
identification. Its unavoidable full-exterior mismatch also prevents calling
the current q=1 code a localized reference-preserving excitation.

## Reused finite results and new applications

| Existing result | Checked interface and application here |
| --- | --- |
| [Complete native history transport](../native_modular_reference/history_transport.py) | Binary66, full retained seam/coin/value tape, partial injective Kraus maps, no intervening coherent control. New: an odd closed lift forces fiber rank one. |
| [Adaptive record recovery](../observer/adaptive_record_recovery.py) | Equal-record differences stay fixed; binary occupancy freezes distinct spectators. This local invariant gives the finite-history argument without assuming stationary or uniform scheduling. |
| [Native exchange topology](../clock_defect_energy/defect_history.py) | Ordered lifts label the admitted paths. Existing homology is retained; it does not insert a quantum exchange phase. |
| [Encoded memory](../common_source/encoded_memory.py) and [optimal record code](../common_source/record_code.py) | All-seam flagged KL, source encoder/recovery, untouched reference and complete padding. New: sharp q=1 optimal dimension two, explicit core-reset failure/control, exterior read. |
| [Record information](../common_source/record_information.py) | The corrected tape alone carries no unknown logical information; quantum logical output and destructive read are separate objects. |
| [Lean balanced-rail memory](../../Lean/Geometry/SourceEncodedMemory.lean) | Checked but not applied: its scalar pair-average/sign encoding is a different state/operation interface from this quantum C/E code. It proves no quantum correction or particle claim here. |

No abstract theorem is extended or silently specialized across those
different interfaces. The local exact proof, finite certificate and executed
code control are the new work.

## Reproduce

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/exchange_sector/test_sector.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.exchange_sector.sector --check
```

The [campaign](campaign.json) binds the predecessor interfaces. The
[receipt](sector_results.json) binds the static repository Python import
closure, this derivation and cited source texts. It records exact claims and
the execution tolerance, avoiding platform-specific floating residual bytes.
