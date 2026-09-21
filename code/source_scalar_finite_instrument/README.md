# Finite-duration scalar quantum readout

This package certifies rectangular field-pointer interactions on the supplied
64-site scalar action. Its observer-like detector has 32 local pointer sites,
field-coupling ports, retained bit and parity slots, and explicit backaction.
The receipt contains all 21 readout slots, their 210 earlier-readout pairs,
and 2,709 declared control operations. It does not sample quantum outcomes.

From the repository root:

```sh
python3 code/source_scalar_finite_instrument/finite_instrument.py --write
python3 code/source_scalar_finite_instrument/verify_finite_instrument.py
python3 -m pytest -q code/source_scalar_finite_instrument/test_finite_instrument.py
```

The producer writes `finite_instrument_receipt.json`; `--check` checks byte
parity. The independent verifier replays the pinned sequential instrument
and its quantum/clock parents. Existing parent files are unchanged.

Every pulse is centered within its inherited reference-clock interval.
The free field evolves during preparation, coupling and readout. Local field
couplings act simultaneously; their common duration is counted once in the
time budget. Pointer preparation and read/decode operations have separate
positive durations. Completed record availability follows the pulse center
by half its duration plus the retained read/decode latency.

The parameter region includes pulse durations from `1/1000` to `1/100`,
pointer-operation duration `1/10000`, and a supplied half-diamond error at
most `1e-7` at each declared fault location. Initial baseline and intervention
states each have an allowed trace-distance preparation error of `1e-5`.
The paired-response bound includes errors in both noisy preparations and both
noisy records. These are unconditional marginals of one joint instrument;
the 21 slots are not independent repetitions.

Pulse-endpoint faults are specified channels, not errors inferred from
Hamiltonian coefficient precision. The ideal energy-work bound does not
bound noisy energy or total apparatus power. The effective instantaneous
smear spreads beyond the hardware coupling region. Initial state preparation,
Born readout and quantum controls are supplied; the 16 historical coherent
preparation factors have no timed implementation in this certificate.
Reference-clock intervals do not certify noisy branch clocks, physical units,
native quantum routing or a regional continuum realization.
