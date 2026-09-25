# Classical common-history feedback packet

The packet extends the authenticated 64-site golden scalar execution by
making its detector an actual sequence of local averaging probes, retained
records, restoration feedback and routed accumulation. Subsequent action
steps consume the restored field values and their actual writers. This is a
bounded observer-like software patch with local state, ports, readback,
records and feedback, on declared classical operations and model time.

The new comparison includes feedback error propagated into later field
states, both preparation errors, detector error, the same finite action's
discretization error, clock readback uncertainty and causal hold error across
the entire interval `[15*tau,17*tau]`. The generic bound permits arbitrary
bounded record errors; deterministic alternating errors provide retained
examples, not an estimated noise law. Six complete histories retain baseline,
intervention, doubled input, two noisy preparations and disabled feedback.
Routing, pointer resets, averaging, restoration arithmetic and accumulation
are exact operations in this contract. The error budget is not a bound for
unmodeled channel corruption or processor faults at those operations.

Run from the research repository root:

```sh
python3 code/common_history_packet/packet.py --check
python3 code/common_history_packet/verify.py --write evidence/common_history_packet_20260925/verification.json
python3 -m pytest -q code/common_history_packet/test_common_history_packet.py
```

The verifier imports no new producer code. It rebuilds the action from the
existing independent edge-energy assembler, replays the immutable full
parent, interprets every consumed writer/value, and checks all restored
configurations. Resealed mutations test substituted action, population,
clock, readout, foreign records, forged writers and disabled feedback.

The clock advances once per discrete action step. Every auxiliary operation
is retained and counted; none is assigned a laboratory duration. Routing,
classical arithmetic, source addresses, preparation, action, record access
and bounded errors remain supplied. There is no quantum outcome, physical
clock, native action selection, spatial continuum or experimental claim.
The completed histories fit 2,456 retained scalar slots each; their largest
rational coefficient uses 185 bits, within the declared 4,096-bit capacity.
This counts exact stored values, not physical memory density or execution time.
