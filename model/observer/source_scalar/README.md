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
Each site uses four or eight number states. The prepared target retains
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

The independent continuous **finite** Hamiltonian gives approximately
`0.516422345746` at `d=8`; the independent infinite two-oscillator reference
gives `0.516421501748`. The latter comparison improves as the occupation
code grows, but is a numerical diagnostic, not a certified cutoff bound.
The source program and independent finite split calculation agree within
`3e-11` in the finest detector probability. Complete branch-density tests
also compare the source instrument with independent direct matrices.

The ideal prepared target has a separate rigorous trace-distance bound from
the infinite coupled vacuum: below `0.007576` at `d=4` and `0.000079899` at
`d=8`. That bound cannot be used alone to compare two different subsequent
evolution or detector channels.

Only the detector bit reaches the controller. The existing exact 528-path
enumeration gives `E[j2]=5r/132`, where `r=Pr(y=+1|j1=1)-1/2`. Two
200,000-trial native runs retain the actual detector response, terminal
branches, complete record counts, reversal, disconnected and independent
donor controls. A new trial freshly prepares its field; measured branches
are not reused as undisturbed vacua.

## Source operations and their cost

Every gate is reconstructed as a full six-level exponential of the twelve
retained source generators. Exact arithmetic checks the base control
identities. Phase-preserving CNOTs, Gray-order basis rotations and diagonal
Walsh expansions then supply every basis change and pointer operation.
Zero-angle requests and scalar phases remain in the schedule.

The full six-level pulse exponentials use 50 decimal digits before conversion
to binary64, reducing drift accumulated through repeated pulse contractions.
The assembled generators and subsequent matrix products remain binary64;
this numerical safeguard supplies no certified roundoff or physical error bound.

The `d=8`, 128-step protocol uses seven work qubits and **1,774,010 native
pulses**, with **3,548,020 proper-code load/export service events**. The
executor contracts each generated subroutine once and reuses its matrix
in the declared repeated schedule. These are logical operation counts,
not host iterations, elapsed physical duration or optimized hardware costs.
Initialization, storage and final readout costs are additional.

The full export channel includes all complement failure outcomes and is
tested on the entire six-level operator algebra and with a reference system.
The executed histories stay in invariant proper codes, so those failure
outcomes have zero probability in exact arithmetic. Coherent transport,
reset, storage, routing and readout remain specified apparatus services;
their physical implementation and noise budget are not established here.

## Error boundaries and next connection

[DERIVATION.md](DERIVATION.md) proves the finite target, preparation tail,
pointer law and uniform splitting bound. [bounds.py](bounds.py) evaluates
the nested-commutator bound with exact rational matrix entries and directed
integer square roots. It controls all inputs, including entangled references,
for the same finite Hamiltonian. Under parity, the bound transfers to
classical native-history total variation by `5/33` and mean current by
`5/132`. Floating compilation error and physical service error are separate.

A full field comparison still needs certified occupation-tail errors during
the impulse, evolution and readout, followed by a spatial refinement family
for the same implementation. The earlier [scalar continuum adapter](../scalar_response/README.md)
uses a different spatial construction; its bounds do not automatically
apply to these two sites. An empirical prediction additionally needs a
physical identification and calibration with independent comparison data.

This construction adds an implemented field apparatus to the source and
observer interfaces. It does not complete continuum physics or OPH.
