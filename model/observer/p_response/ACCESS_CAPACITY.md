# Detector access and the conditional capacity budget

Interface contract, 12 October 2026. This note separates the apparatus
that transcribes a native current, the actuator that consumes its detector
bit, and the richer audit record.

The native first event produces `r1=(c1,j1)`. The declared probe uses `j1` to
route the phase `g(P) j1 chi`, then writes a detector value `y in {-1,+1}`.
The second-seam selector accepts **only y**. Therefore the actuator's record
algebra is `C({-1,+1})`, even though the probe router uses a native current
and the full simulation can audit it. A binary actuator does not imply that
the complete apparatus has only one bit of state.

For the second native event, a precise preservation contract is

```
(loads, committed y, coin) -> (native(loads, seam(y), coin), same y).
```

Equivalently, condition the native working channel on the classical detector
block and apply identity to its committed label. The projectors onto y are
fixed by that continuation. It is unnecessary, and generally false, to infer
preservation of coherence between different record labels from this classical
property. The first native label must not be silently supplied to the seam
selector as another input.

The implemented `decoder_seam(detector_bit)` and
`controlled_native_step(native, loads, detector_bit, coin)` satisfy this access
split. A separate endpoint floor/ceiling oracle checked all 312 combinations
of 78 source states, two detector values and two coins; every returned native
state, current letter and preserved detector bit agreed. The separate P-coupling calculation and its statistical experiment are
covered by the [calculation receipt](results.json).

## Two different record capacities

At the declared nonsaturated phase, both detector outcomes have strictly
positive probability for every native current. Independent floor/ceiling
enumeration of the 66-state binary source gives:

| Retained stage | Distinct reachable records |
| --- | ---: |
| Empty prefix | 1 |
| First native letter `r1` | 4 |
| First letter and detector `(r1,y)` | 8 |
| Terminal full audit `(r1,y,r2)` | 24 |

The full audit has 37 prefixes in total and thus admits a six-bit fixed-width
classical history encoding. Its 24 terminal records are not all visible to
the bit-only actuator. With exact reread and halt continuation, the existing
correctable-record criterion gives `M0=2` for the detector-bit interface and
`M0=24` for the full terminal audit interface. These are different declared
algebras. The phase, probe, source state and program need their own accounts.
The six-bit count is not a quantum-space compression theorem or a proof of
minimal hardware.

At this exposure, changing P through its specified coupling map changes the
detector probabilities and record/response correlations. It does not change
the two-symbol detector alphabet or the 24-word support, because all detector
branches remain positive. Response sensitivity and zero-error storage
capacity are therefore different uses of the same candidate parameters.

## A reusable conditional P/N area inequality

The existing OPH dictionary uses `P = a_cell / ell_star^2` and records
capacity in nats. If an **identified cut enclosing the retained record**
obeys the area-capacity law, and its area is K equal geometric chart cells,
then

```
N_cut = A_cut / (4 ell_star^2) = K P / 4,
ln M0 <= N_cut,
therefore K >= 4 ln(M0) / P.
```

Using the pinned Thomson-branch P enclosure, direct high-precision numerical
evaluation gives `4 ln 2 / P` approximately `1.6999608573` and
`4 ln 24 / P` approximately `7.7942567832`. Both enclosure endpoints stay
well inside their respective integer intervals. This yields the conditional
necessary area budgets **K >= 2** for a retained detector bit and **K >= 8**
for the full 24-record terminal audit. This arithmetic is not a new interval
certificate or a hardware sufficiency result.

The one-cell sensing aperture in this design is a different object from the
cut enclosing its memory and controller. It must not be assigned the entire
observer's record budget by relabeling it. Likewise, a chart cell carrying
`P/4` nats of area density is not an autonomous Hilbert factor of fractional
dimension. The [existing discussion of cell/edge capacity](../../../paper/observers_are_all_you_need.tex)
already makes that distinction.

The [public capacity theorem](../../../paper/tex_fragments/OBSERVERS_N_CLOSURE_PRO5.tex)
provides the correctable-record bound. Its cosmic relation
`Lambda ell_star^2 = 3 pi / N_star` additionally identifies the same globally
selected capacity with a de Sitter horizon. Accepting that global relation
does not by itself identify this apparatus's local record cut. The useful
simplification here is to reuse the capacity inequality for the actual
declared record interface, while keeping the needed cut attachment explicit.
