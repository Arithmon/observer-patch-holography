# A native reference state supplies a modular frequency

This packet connects the existing native event statistics to a programmed,
observer-relative modular readout. The faithful Palm reference is
`sigma = diag(degree/600)`, derived in
[`event_circulation`](../event_circulation/README.md). Its two eigenvalues have
ratio `5/4`. The reference therefore supplies a nonzero modular frequency
`log(5/4)` without choosing a Gibbs state from an externally supplied
Hamiltonian. This is a finite reference-state attachment, not a recovery of
interacting spacetime or a new comparison with experimental data.

The reference is **not the probe preparation**. The executed probe begins at
`I66/66`, uses the existing fixed-seam signed-current herald, and prepares a
coherent mixture with the compiled native balanced mixer. Each herald sign
has probability `5/66`; failed preparation has probability `28/33` and remains
in the raw record. The ten orientation pairs have degree census
`(8,10):2`, `(10,8):2`, `(8,8):2`, `(10,10):4`. Thus four pairs see a
nontrivial relative modular phase and six do not.

## Exact relation and executed computation

With Heisenberg convention `sigma^(is) A sigma^(-is)`, the dual state action
uses `sigma^(-is)`. On the full native algebra,

```
-log(sigma) = log(60) I + log(5/4) P_degree8.
```

The twenty-dimensional active corner commutes with the reference. There the
implemented diagonal target is `D(theta)=exp(i theta P_degree8)`, up to the
common scalar phase, with `theta=s log(5/4)`. Normalizing the restricted
reference changes only that scalar. Its full-reference weight is `23/75`.
The chosen extension assigns zero phase to the remaining native addresses;
it omits 22 other degree-eight states. **It does not implement the modular
flow on the full 66-dimensional native algebra.**

Each observer program executes:

1. Native signed-current herald, retaining failure and both signs.
2. On success, the 34-pulse compiled balanced mixer `B`.
3. A fresh compiled diagonal source word `D` on seven work bits, lifted by
   two spectator bits to all 512 encoded and padding coordinates.
4. The same 34-pulse mixer `B`, then the native pointer with sign-aware
   decoding. Failure receives its explicit read record and no phase cost.

The diagonal compiler acts on eight declared degree-eight addresses among
the twenty active addresses. All source matrices are raw binary pulse
products, with no polar projection. Padding remains in the complete runtime;
native herald padding is failure and pointer padding has an unbiased bit.
There are six joint herald/read records, including the zero-probability
failed-preparation/read-one record. The C and E candidate laws are both run.
The reference preparation itself reuses the separately executed Palm-herald
construction; this packet does not rerun that experiment as a probe.

All settings consume the actual `runtime.Program` events, controls and
Kraus arrays with a factor executor: the 66 initial basis columns retain
their original binary weight `1/66`, and every hidden Kraus alternative
adds a column. Final densities are weighted sums of outer products on all
512 coordinates. This avoids subtractive cancellation at exactly zero
ideal probabilities without discarding small amplitudes, projecting states
or normalizing outcomes. At phase `pi/2`, the generic dense observer runtime
independently reproduces the full output and serialized checkpoint
continuation. Its relative positivity guard can reject cancellation noise
at nominally zero outcomes; it is not relaxed for this packet.

The exact ideal conditional sign-match law is

```
P(match | successful herald, s) = 4/5 + cos(s log(5/4))/5.
```

| Programmed phase | Modular parameter | Conditional match | Joint match, each sign |
|---|---|---|---|
| 0 | 0 | 1 | 5/66 |
| pi/2 | pi/(2 log(5/4)) | 4/5 | 2/33 |
| pi | pi/log(5/4) | 3/5 | 1/22 |
| 2 pi | 2 pi/log(5/4) | 1 | 5/66 |

These settings use exact rational-complex target entries `1`, `i` and `-1`.
No binary logarithm or exponential is treated as an exact source target.
The displayed decimal frequency is diagnostic only. Uniform reference
`I66/66` gives identity modular action and conditional match one at every
setting. Conversely, the native unread dissipative attempt has a different,
law-dependent contraction. Programmed modular control is not asserted to be
that natural waiting evolution.

## Error boundary and operational meaning

The certificate converts every original binary component of the local
128-by-128 phase map to an exact rational, compares every coordinate against
the exact diagonal target, and encloses its Frobenius error `epsilon`.
Consequently its operator error is at most `epsilon`, its half-diamond error
at most `delta=epsilon(2+epsilon)/2`, and its channel norm at most
`g=(1+epsilon)^2`. The spectator lift and arbitrary untouched quantum
references preserve these operator bounds. The existing eight-dimensional
mixer certificate is reused. For the whole raw `B D B` sequence the bound is

```
delta_B (1 + g_D g_B) + delta_D g_B.
```

Both mixer uses and numerical amplification are charged. Exact ideal native
instruments and adaptive failed branches contract this comparison. This
bound covers the compiled gates against their stated complete extensions;
it does **not** certify binary native Kraus production, subsequent floating
state propagation, or hardware service errors. Full terminal density and
probability comparisons and checkpoint replay are separately labeled
numerical diagnostics. Physical service discrepancy remains unbounded.
For portable replay the fresh exact phase comparison must pass the fixed
operator-error envelope `1e-10`; the mixer must pass `1e-11`. The exported
whole-sequence half-diamond allowance is outward rounded from those fixed
envelopes, rather than depending on platform-specific last bits.

The law is an ensemble clock candidate on the monotone window
`0 <= s log(5/4) <= pi`; one bit is not an exact clock reading, and sign and
period aliases remain. A noncommuting prepared probe is necessary: the
reference itself is stationary under its modular flow. Coherent services,
the probe design and selected phase settings remain supplied. Calling `s`
an observer's candidate time needs no universal clock or conversion to SI
seconds. What remains is to establish consistency with continuing clock
records, other accessible correlations and interventions, and effective
field dynamics. This packet supplies a concrete reference-derived frequency
and testable finite readout law, rather than that larger identification.

The finite modular-inner attachment is the existing
[`PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean); the current
construction uses explicit finite spectral calculus. The prior native
interface supplies the exact herald and pointer laws, while
[`compiled_native_interface`](../compiled_native_interface/README.md) supplies
the complete source mixer and its numerical comparison theorem. No additional
axiom or newly claimed abstract theorem is needed.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference --check model/native_modular_reference/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_modular_reference
```

Controls cover the faithful reference and 4+6 pair census, omitted native
coordinates, all padding columns, inverse/common phase, scalar amplification
that defeats an unweighted telescope, retained preparation failures, uniform
reference, native waiting, checkpoint continuation and hostile ndarray
subclasses. The retained receipt is software evidence, not measured data.
