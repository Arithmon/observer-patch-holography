# Scalar preparation, evolution and detection through source controls

This runnable construction compiles a two-site quantum scalar apparatus into
the existing twelve-port source operations. It starts with blank registers,
prepares an entangled field, applies the first native current as a momentum
impulse, evolves the field, reads a quantum pointer and passes its committed
bit to native repair. Both conditional field states are retained.

The new connection is an **emitted and executed finite source program** for
the field apparatus. Its scalar action, site interpretation, coherent code
services and model clock are supplied construction choices. This is a
completed finite software construction, with no claim that the action was
derived from integer repair or established experimentally.

## Run and inspect

From the repository root after `python3 -m pip install -r model/requirements.txt`:

```sh
python3 -m model.observer.source_scalar --check model/observer/source_scalar/results.json
python3 -m pytest -q -W error model/observer/source_scalar
# Optional complete control routines and their repeat recipe; refuses overwrite.
python3 -m model.observer.source_scalar --emit /tmp/source-scalar-program.json
```

`--write PATH` retains a new result without overwriting one. Each primitive in
the emitted program specifies a control kind, angle and register addresses;
`compiler.integrated_controls(kind, angle)` gives its twelve source-port
coefficients. The retained source summary contains their exact algebraic
definitions and source hashes. Neither target exponentials nor a requested
vacuum vector replace the executed source gates.

## What ran

The chosen Hamiltonian is
`H=(P0²+P1²+Q0²+Q1²)/2 + (9/64)(Q0-Q1)²`, with normal frequencies 1 and 5/4.
Each site uses four, eight or sixteen number states. The prepared target retains
total occupation through `d-2` from the analytic infinite vacuum and is
compiled from the all-zero state. The evolution uses 64 or 128 symmetric
split steps to reach the declared model coordinate `t=1`. These are neither
physical seconds nor native repair counts.

For a positive unit native current, the retained results are:

| Onsite levels | Steps | Source-executed `Pr(y=+1)` | Exact rational finite splitting bound, rounded upward |
| --- | --- | --- | --- |
| 4 | 64 | 0.516680958 | 0.000103957 |
| 4 | 128 | 0.516679563 | 0.000025990 |
| 8 | 64 | 0.516424174 | 0.000693013 |
| 8 | 128 | 0.516422803 | 0.000173254 |
| 16 | 64 | 0.516423330 | 0.005076875 |
| 16 | 128 | 0.516421959 | 0.001269219 |

The independent continuous **finite** Hamiltonian gives approximately
`0.516422345746` at `d=8` and `0.516421501748353` at `d=16`;
the independent infinite two-oscillator reference
gives `0.516421501748`. The latter comparison improves as the occupation
code grows; that numerical diagnostic is separate from the rigorous
whole-history cutoff certificate below.
The source program and independent finite split calculation agree within
`2e-10` in the finest detector probability. Complete branch-density tests
also compare the source instrument with independent direct matrices.

The ideal prepared target has a separate rigorous trace-distance bound from
the infinite coupled vacuum: below `0.007576` at `d=4`, `0.000079899` at
`d=8` and `0.000000010324` at `d=16`. That bound cannot be used alone to
compare two different subsequent evolution or detector channels.

Only the detector bit reaches the controller. The existing exact 528-path
enumeration gives `E[j2]=5r/132`, where `r=Pr(y=+1|j1=1)-1/2`. Three
200,000-trial native runs retain the actual detector response, terminal
branches, complete record counts, reversal, disconnected and independent
donor controls. A new trial freshly prepares its field; measured branches
are not reused as undisturbed vacua.

## Source operations and their cost

Every gate is reconstructed as a full six-level exponential of the twelve
retained source generators. Exact arithmetic checks the base control
identities. Phase-preserving CNOTs, Gray-order basis rotations and diagonal
Walsh expansions supply every basis change and pointer operation. The
larger state preparation uses a binary amplitude tree with commuting
Walsh-Y factors: 5,379 source pulses for the eight field qubits. The
smaller preparations retain their original source words.
Zero-angle requests and scalar phases remain in the schedule.

The `d=16`, 128-step protocol uses nine work qubits and **18,612,125 native
pulses**, with **37,224,250 proper-code load/export service events**. The
executor contracts each generated subroutine once and reuses its matrix
in the declared repeated schedule. These are logical operation counts,
not host iterations, elapsed physical duration or optimized hardware costs.
Initialization, storage and final readout costs are additional.

Before repetition, a guarded polar decomposition restores numerical
unitarity to each matrix computed from the source pulses. It uses no target
matrix and preserves the word's phase. The receipt records each raw defect,
correction and final defect. A raw unitarity defect above `1e-10` or a
correction above `5.1e-11` fails. Matrices whose nonunitarity exceeds the
declared rounding tolerance are rejected. This prevents backend-dependent
rounding bias from accumulating across the long program. Outcome probabilities and conditional
states still pass the original independent checks. These diagnostics are
not an interval arithmetic certificate or a physical error correction service.

The full export channel includes all complement failure outcomes and is
tested on the entire six-level operator algebra and with a reference system.
The executed histories stay in invariant proper codes, so those failure
outcomes have zero probability in exact arithmetic. Coherent transport,
reset, storage, routing and readout remain specified apparatus services;
their physical implementation and noise budget are not established here.

## Error boundaries and next connection

[DERIVATION.md](DERIVATION.md) proves the finite target, complete ideal
occupation-history bound, pointer law and uniform splitting bound.
[bounds.py](bounds.py) evaluates
the nested-commutator bound with exact rational matrix entries and directed
integer square roots. It controls all inputs, including entangled references,
for the same finite Hamiltonian. Under parity, the bound transfers to
classical native-history total variation by `5/33` and mean current by
`5/132`. Floating compilation error and physical service error are separate.

[cutoff.py](cutoff.py) uses exact rational Gaussian number-moment bounds and
projected Duhamel comparison throughout the impulse, continuous evolution
and both pointer arms. Its trace-distance bound includes the terminal
quantum branches without dividing by their probabilities. The bound is
state-specific for this prepared history, not an all-input diamond norm.

| Onsite levels | Complete ideal occupation bound, rounded upward | Source word executed |
| --- | --- | --- |
| 4 | 0.799447 | Yes |
| 8 | 0.077171 | Yes |
| 16 | 0.000196 | Yes |

Adding the finite splitting bound gives a common ideal instrument error
budget relative to the infinite two-oscillator reference. The occupation
bound tends to zero with the cutoff. For the executed `d=16`, 128-step design,
the combined ideal error is below `0.001464436`. A separate exact rational
enclosure of the infinite Gaussian probability then gives

```
0.514957066204 <= ideal Pr(y=+1 | j1=1) <= 0.517885937293.
```

The endpoints displayed here are rounded outward; the receipt retains
exact rational endpoints. This interval excludes the fair-detector response
`1/2`, so the declared ideal coupled apparatus predicts a distinguishable
positive-current response. The enclosure uses bounded sine and exponential
series, the full occupation bound and the finite splitting bound; it does
not use the floating execution result as a proof. Numerical execution and
physical service errors require additional terms. This is a prediction of
the specified apparatus, with no claim that it identifies a measurement in
nature.

A field comparison additionally needs a spatial refinement family for the
same implementation. The earlier [scalar continuum adapter](../scalar_response/README.md)
uses a different spatial construction; its bounds do not automatically
apply to these two sites. An empirical prediction additionally needs a
physical identification and calibration with independent comparison data.

This construction adds an implemented field apparatus to the source and
observer interfaces. It does not complete continuum physics or OPH.
The separate [C/E prediction test](../../predictions/README.md) distinguishes
two finite native instruments with a fixed complete-count protocol; it does
not use this scalar cutoff bound or assert a laboratory result.
