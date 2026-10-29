# Proper-code compiler attachment probe

Verdict: **SOURCE_INTERFACE_UNDERDETERMINED**. Work baseline:
`a44330802c726037f88356f9f62044184ed429a5`; scientific source bytes:
`FloatingPragma/observer-patch-holography@53e365b248a4611bb911e75d7afae0f52bfd8675`.
The accepted Phase-0 STOP remains in force. This finite probe tests an existing
export service; it does not implement a spatial/internal operational bridge.

## Objects and exact attachment

The spatial object is the real slow subspace of twelve-component load
readouts, with projector `P_s`. Native total-two configurations have 66 binary
or 78 full basis states. `slow_geometry.load_readout` audits integer loads;
it does not prepare a quantum state on the rank-three quotient. Native
configuration reset is a measure-and-prepare map with retained timeouts;
readout alone retains a path-dependent work state and is not a Lueders
instrument on that work. Coherent synthesis is a separate admitted service. Full native
tomography identifies states on the 66-dimensional sector; it does not
identify that sector with the slow real vector space.

The existing `model/predictions/prepared_clock_source.py` also supplies a
declared `B(C12)` comparison workspace inside sixteen protected addresses,
with a source-Laplacian spectral basis, finite sampler and heralded prepared
reference. This is an actual additional quantum domain, not an absent
capability. Reusing its incidence spectrum does not itself identify native
load/current observations with amplitudes or density operators on that
workspace. This export probe supplies no such preparation/observable map.

The internal object is `H_upper = span_C(e0,e1,e2)` in the declared six-level
response processor. Full internal matrix effects and calibrated preparations
belong to that abstract source architecture. The compiler uses the same
twelve source generators and exact control rows for `X,Y,Z,I,H,J`; coherent
proper-code load/export and routing remain declared services. None of these
inputs identifies a real load vector with a density operator.

Let `phi=(1+sqrt(5))/2`. Transposing the existing native raw frame gives the
source-bound linear map `T : R^12 -> R^3`:

```text
T = [-1, 1,-1, 1, 0, 0, 0, 0, phi, phi,-phi,-phi
     phi,phi,-phi,-phi,-1,1,-1,1, 0, 0, 0, 0
      0, 0, 0, 0, phi,phi,-phi,-phi,-1,1,-1,1].
```

Exact arithmetic gives `T T* = c I3`, `T* T = c P_s`, `c=10+2sqrt(5)`.
Thus `T*/c` is its inverse on the slow subspace. All 60 source actions satisfy
`T P(g)=U_upper(g) T`. The independent verifier reconstructs the frame via
response vertices in `Q(phi)` and recomputes character inner products:
`dim Hom(slow,upper)=1`, `dim Hom(slow,lower)=0`. This fixes the mathematical
attachment and its source normalization; it supplies no quantum encoding.

The marker is transported from the frozen generator **before** any character
is evaluated. It selects path index 4, with slow/upper trace
`(1+sqrt(5))/2` and lower trace `(1-sqrt(5))/2`. Historical path index 1 has
the opposite pair. An inverse element with the same character is not accepted
as a replacement marker. No sign criterion chooses a block or a twist.

## The actual export and its bounded obstruction

`compiler.transfer_kraus(d, coherent)` maps `M6` to `M_(d+1)`, for the source's
admitted `d=2,4`. In the coherent case its successful Kraus map is
`K_s=sum_(j<d)|j><j|`; each complemented column has a Kraus map
`K_j=|failure><j|`, `j>=d`. Dephased success instead sums separate
`|j><j|` CP maps. The two public sectors are success and failure. The
individual complementary Kraus indices are not added as public read access.
The receiver overwrite and sender reset remain the fixed outputs omitted by
the source's rectangular-map convention. This is not their joint retained
apparatus instrument or a timeout law for physical services.

| Source service | Upper-sector state transfer | Success effect under all 60 upper actions |
| --- | --- | --- |
| Coherent `d=2` | Truncates the third upper coordinate | Commutes with only 4 actions |
| Dephased `d=2` | Truncates and dephases | Commutes with only 4 actions |
| Coherent `d=4` | Exact upper-sector identity in the existing first three code coordinates, with zero failure there | Commutes with all 60 actions |
| Dephased `d=4` | Same success effect, but upper coherences are destroyed | Commutes with all 60 actions |

For the fixed source marker and the first source upper basis input `|0>`,
direct `d=2` export has success/failure probabilities `(1,0)` before the
action and `(3/4,1/4)` after it. This rules out a label-preserving covariant
direct `d=2` export on the whole upper triplet. This is an internal-domain
unit control, not a preparation of a spatial vector or a physical experiment.
The source also admits `d=4` and multi-pulse register operations; therefore
this witness is not a universal compiler impossibility theorem.

The coherent `d=4` positive control preserves all nine upper matrix units,
including their surviving quantum states. It therefore preserves arbitrary
untouched-reference extensions by complex linearity. Each of the four
service cases is checked on all 36 input matrix units, with every failure
branch and trace preservation retained. The independent route constructs
complete Choi blocks and checks their output partial trace. A rational Bell
input on the existing internal sector and a two-level untouched reference
also tests correlations: dephasing keeps the success effect and reference
marginal but loses the off-diagonal joint state. This reference is a unit-test
input, not an ancillary triplet introduced into the source.

The `d=4` statement is restricted to `H_upper`. On full `M6`, this export
retains only one lower-block coordinate; its full success effect commutes
with only 4 actions, and the receipt includes a separate lower-input witness.
No full-processor covariance is inferred from perfect upper transfer.

## What remains missing

The positive export control identifies what the *declared internal code
service* can preserve. It does not derive a source identification of the
spatial load quotient with an admitted quantum domain, its preparation and
observable dictionary, or an executable operation encoding that domain into
`H_upper`. The generic real-orthogonal
compiler can synthesize requested logical matrices, but this probe freezes
no marked-action pulse family or rigorous complete compiled-word error budget.
An induced algebraic action on a code subspace is not reported as an executed
operation. In particular no arbitrary isometry, triplet register, clock,
Galois selector or target-fitted transport fills the gap.

A falsifier of this underdetermination is a source-derived native preparation/observation
contract identifying the spatial observation domain with an admitted quantum
workspace and providing a coherent encoding, followed by source-bound marked words and complete
service/error comparison. Losing the `d=2` probability witness under the
fixed source data, or a nonidentity upper `d=4` matrix-unit image, falsifies
the corresponding finite result. Physical golden selection is false throughout.

## Historical drift and scientific inventory

The whole-law checker still rejects `DERIVATION.md`. Its pinned bytes were
independently located at commit `a84438d99a355dfd31e39aeb667e00d8ea60e445`.
The diff is one paragraph: the older reference to RFC A1-R/A2-R is replaced
by a statement about complete primitive grammars and completed reconciliation
not being silently added to the canonical axioms. The other three pins
(`geometry.py`, `response.json`, `verify_response.py`) match. The receipt binds
the old Git object and the exact diff hash. This scope wording deserves its
own review; no historical hash is updated and no whole-law replay is claimed
to have passed its guard.

The new files audit a fixed declared export and existing representations.
They provide no twelve reversible port perturbation families, mixed-order
raw histories or refinement ancestry. They are not a new Stage-2 source-current
candidate. Regenerating its canonical snapshot must retain the 33-candidate,
zero-qualifying verdict; snapshot refresh is not scientific promotion.
Source-selection, construction and observation ledger statuses remain unchanged.

## Reproduction and evidence

The producer executes the source's export API and uses its exact integral
entries. The independent verifier imports neither that producer nor the
compiler: it checks the export function's restricted AST and reconstructs
Choi blocks from its index rules. Both bind Git commit/blob, SHA-256 and byte
counts. All arithmetic claims use rationals and `Q(sqrt(5))`/`Q(phi)`, with
no floating tolerance. Tests reject changed marking, Galois-swapped frame,
wrong normalization, missing failure, dimension mismatch, discarded reference
coherence, diagonal-only coverage, counterfeit pins and physical/operational
promotions. No heavy pool, dependency installation or physical data is used.

```sh
python3 -B code/a5_closure/proper_code_attachment_probe.py --check
python3 -B code/a5_closure/verify_proper_code_attachment_independent.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B -m pytest -q -W error -o cache_dir=/tmp/oph-probe-cache code/a5_closure/tests/test_proper_code_attachment_probe.py
```

Frozen output: `receipts/proper_code_attachment_probe.json`. To generate a
fresh comparison without overwriting it, use the producer's `--write` with a
new path, then pass that path to the independent verifier.

The next bounded task is to audit the existing protected `C12` preparation
as a candidate domain and derive or refute, within its admitted grammar,
the missing native-current preparation/observable identification with its
slow spectral sector. This precedes any new spatial-to-internal pulse family.

Local validation used the existing CPython 3.14.4 environment; no dependencies
were installed. The probe's exact arithmetic and independent Choi reconstruction
pass, with 30 probe/hostile tests, 14 existing marking tests plus 13 subtests,
13 targeted compiler/export tests (56 unrelated compiler cases deselected),
and 25 existing scale-contract tests. The standard mandatory shard at index 8
of 10 passes 426 tests and 17 subtests and separately collects 20,778 tests;
collection does not execute that complete set. The canonical inventory has
914 files and retains 33 candidates, 0 qualifying. All five register checks
pass without ledger edits. The whole-law guard's known prose-pin failure is
tested as an expected negative. Full mandatory/model and cross-platform pinned
CI replay are not claimed by these local results.
