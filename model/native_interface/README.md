# Native preparation and phase readout with one shared mixer

The [quantum-history construction](../native_quantum_history/README.md)
supplied a coherent density matrix and a phase detector independently. This
packet replaces both with **actual native records and one declared balanced
real mixer**. It executes a complete preparation, evolution and readout
history through the primitive observer runtime. The mixer remains an extra
coherent service: the committed native grammar cannot construct it from
diagonal preparations and classical feedback alone.

Both signs of a moving native current can herald a useful preparation. The
probability of a useful trial is `5/33`, with `5/66` for each sign. Failed
preparations, probability `28/33`, remain explicit outcomes. There is no new
experimental dataset or claim that these model probabilities have been
measured in nature.

## The end-to-end operation

The carrier is the existing 66-dimensional binary sector: two occupied sites
among twelve. Begin with the declared stationary diagonal state `I/66`. Choose
the first source seam `(u,v)` and execute one native C or E repair with its
fair coin `c`. Read its signed current `j=x_u-x'_u`.

1. If `j=0`, retain a failed-preparation record and stop using native repairs
   or mixers in that trial. The controller's remaining slots act as identity.
2. If `j=-1`, the repair has prepared the uniform mixture of the ten
   configurations `|uw>`, with `w` outside the seam. If `j=+1`, it has instead
   prepared the corresponding mixture of `|vw>`. The sign remains recorded.
3. Apply the same declared mixer `B` to each of these ten orientation pairs:
   `B|uw>=(|uw>+|vw>)/sqrt(2)` and
   `B|vw>=(|uw>-|vw>)/sqrt(2)`. It acts as identity on the other 46 states.
   Thus the two herald signs prepare the previous `rho_+` and `rho_-` exactly.
4. Execute `n` unread native attempts with the existing uniform seam scheduler,
   then apply `B` again.
5. Execute a final fixed-seam native repair. Declare its raw pointer outcome
   plus exactly when `2*j+2*c-1>0`. This uses only the recorded native coin
   and current. For the minus preparation, reverse the final bit, so the
   interpreted readout always means **matches/opposes the prepared phase**.

The runtime receives the diagonal initial state, actual native Kraus maps and
`B`. It receives no independently supplied target coherent state or target
phase POVM. The old apparatus is used only as a comparison oracle. Successful
trials cost `n+2` native attempts and two applications of the supplied mixer;
failed trials cost one native attempt and no mixer. The parameter `n` counts
native attempts, not seconds or a continuous-time Hamiltonian flow.

For each candidate and each `n=1,2,4`, the executable retains all six records:

| Herald | Interpreted final bit | Exact probability |
| --- | --- | --- |
| Failed | 0 | `28/33` |
| Failed | 1 (impossible) | `0` |
| Plus preparation | Match | `(5/66) p_L(n)` |
| Plus preparation | Opposite | `(5/66) [1-p_L(n)]` |
| Minus preparation | Match | `(5/66) p_L(n)` |
| Minus preparation | Opposite | `(5/66) [1-p_L(n)]` |

The old exact source transfer gives

```text
p_L(n) = [1+r_L^n]/2,       r_C=5/6,       r_E=7/10.
```

| Unread native attempts | C: conditional match | E: conditional match |
| --- | --- | --- |
| 1 | `11/12` | `17/20` |
| 2 | `61/72` | `149/200` |
| 4 | `1921/2592` | `12401/20000` |

At four attempts, the conditional difference is `12259/101250`, about
`0.12108`. Including all failed preparations, the total variation distance
between the two complete six-record distributions is **`12259/668250`**, about
`0.01834`. Using both herald signs doubles the usable yield relative to keeping
just one. Four is the maximizing integer for this particular conditional
phase contrast, as established by the existing quantum-memory calculation;
no optimum over other preparations, detectors or quantum channel tests is
claimed.

All terminal densities are compared with the independently reconstructed
source transfer blocks, including every failure and impossible record. A
checkpoint after preparation is serialized, independently restored and
continued. A separate entangled reference-qubit input tests every complete
132-dimensional terminal density and the untouched reference marginal.
Entanglement survives its successful herald, so this check is not reduced to
a classical input test. Its herald probabilities differ from those of `I/66`;
the stationary `5/66` law is not asserted for arbitrary inputs.

## Why native records implement the phase effects

Let `P_u` and `P_v` project onto the two ten-dimensional orientation subspaces,
let `Z=P_u-P_v`, and let
`S=sum_w (|uw><vw|+|vw><uw|)`. Exact native branch counts show that the raw
coin/current rule has effects

```text
D_+ = (I+Z)/2,       D_- = (I-Z)/2.
```

The effects are deterministic on either orientation and give a fair bit on
the remaining 46 states. Conjugation by the balanced mixer gives
`B* D_± B=(I±S)/2`: exactly the earlier phase effects. This equality holds on
every input, including with an untouched reference. Omitting the native coin
and assigning a fair bit to `j=0` instead gives the weaker effect
`I/2±S/4`. Retaining the coin correlation is therefore essential.

Equal effects do not identify the full postmeasurement operation. On `rho_+`,
the native detector reports plus with certainty but changes its state. Its
terminal density has trace distance **`1/2`** from the former positive-square-
root detector output. The producer explicitly checks that distinction. The
effects are projectors only after restriction to the active 20-dimensional
orientation corner; globally their inactive eigenvalue is `1/2`.

The existing finite theorem
[`instrument_component_fixing_certain_states_is_lueders`](../../Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean)
requires a projective effect and invariance of every state certain for that
outcome. Applied to the active corner, its invariance premise fails here.
Thus that theorem does not identify this native detector with the old
instrument. It is not applied to the nonprojective whole-carrier effects.

## The exact reachability boundary

For every source seam, coin and retained candidate record, the native Kraus
matrix is a weighted **partial injection**: each nonzero column and row has
one entry. The producer checks all 3,960 source branch addresses per
candidate, grouped into 120 C or 180 E matrices across the 30 seams. This
proves, with `Delta` denoting work-basis dephasing,

```text
Delta(K rho K*) = K Delta(rho) K*,       K* K is diagonal.
```

Products and tensor products of partial injections have the same property.
Coarse-graining outcomes sums their CP maps without adding amplitudes between
records. Hence every finite path built from these native maps, classical
record conditioning, classical feedback and diagonal ancillary preparations
preserves diagonal states. Conditional normalization cannot create a missing
off-diagonal entry. Every accessible classical record effect in this grammar
is diagonal, so those record probabilities depend only on the input diagonal.
The same argument retains an arbitrary untouched reference and its conditional
marginals. It does **not** equate the complete quantum terminal states of
coherent and dephased inputs: native maps can preserve coherence that was
already present. An adaptive two-step reference test checks that distinction.

The [source repair instrument theorems](../../Lean/Geometry/SourceRepairInstrument.lean)
`dephase_kraus`, `dephase_cptp` and `dephase_diagonal` supply the finite
dephasing operation. The displayed intertwining identity is the direct
partial-injection attachment for these actual native branches; no new Lean
theorem is claimed. This is a reachability obstruction for the declared
operational grammar, not for every possible OPH interface.

The mixer has two nonzero entries in each active column and creates coherence
from a heralded diagonal state. It therefore lies outside that native-only
closure. Native support derives its pairing, but neither membership in a
generated matrix algebra nor the existence of coherent finite-theorem
constructions supplies this operation as an executable native service. The
positive reduction is to **one shared coherent service**; it is not a proof
that this service is physically minimal or source-selected.

Even adding `B` leaves every matrix real. All finite adaptive record effects
remain real symmetric and therefore cannot distinguish a state from its
complex conjugate. In particular the two mixtures of
`(|uw>±i|vw>)/sqrt(2)` are orthogonal (trace distance one), while both give
the same half/half native pointer law. The executable checks this control.
The existing
[`rhoY_bornWeight_equal_on_real_symmetric`](../../Lean/QFT/SourcePhaseLiftBridge.lean)
theorem gives this identity in dimension two; compression to each orientation
pair applies it here. The 66-dimensional real-effect closure is the direct
finite-matrix argument above, not a claim that this larger apparatus already
has a new Lean instantiation or belongs to the theorem's original generator
set. Full complex tomography still needs a resource breaking this real-effect
restriction.

## What is still supplied and how to reproduce

The binary carrier sector, `I/66` input, chosen seam, uniform evolution
scheduler, candidate C or E, balanced mixer and classical protocol controller
remain premises. The construction does not derive their physical availability,
select C over E, emit source pulse words, define a physical clock, establish
an interacting continuum field theory or match a new laboratory observation.
It removes independent coherent preparation and phase-effect assumptions
conditional on a single clearly identified mixer. The previously derived
first-current phase probabilities can use these conditional preparations and
the same effects, but this packet executes fixed native attempt counts only.

From the repository root, with `model/requirements.txt` installed:

```sh
python -m pytest model/native_interface/test_construction.py -q
python -m model.native_interface --check model/native_interface/results.json
python -m model.native_interface --write /tmp/native-interface-new.json
```

`--write` refuses to overwrite existing files. `--check` rejects malformed
JSON and compares source hashes, exact fractions, complete outcome dictionaries
and numerical diagnostics using the shared portable receipt comparison.
`manifest.json` binds the packet. Exact rational source counts establish the
stated fractions; floating density calculations are residual diagnostics with
declared tolerance `2e-10`, not outward roundoff certificates.
