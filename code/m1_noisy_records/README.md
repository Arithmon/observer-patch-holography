# Noisy control and live public histories

This extends the merged fixed-rate source result by removing reliable runtime
classical control. The [proof](../../extra/NOISY_SOURCE_RECORDS.md) constructs
protected causal decisions, disjoint export islands and actual noisy classical
archive maintenance. The [contract](CONTRACT.md) fixes the scope and exit.

The finite evidence includes full adaptive-instrument Choi comparisons, exact
Toffoli phases, and all 8,100 single faults of a 515-location archive circuit,
with an existing input error. All outcomes and live-wire idles are included.
It also constructs the [[2047,1,63]] punctured Reed--Muller computation code,
checks its complete generator algebra, and supplies an analytic distance proof
and a conservative spread margin for the imported universal construction.
The arbitrary-size archive lemma and universal source protection are analytic;
the measurement-free universal gadgets are credited Aharonov--Ben-Or results,
not certified by the small archive or the parent's different recovery census.

From the repository root:

```sh
PYTHONPATH=code python -m m1_noisy_records.build
PYTHONPATH=code python -m m1_noisy_records.verify
python -m pytest -q code/m1_noisy_records
```

PowerShell: set `$env:PYTHONPATH='code'` first. CI replays on Linux and Windows.
`independent.py` reconstructs embedded gates and compares against independently
formed branch projectors. `archive_check.py` uses Boolean arrays and validates
every layer's ownership and idle coverage, independently of the producer's
integer-bit execution. Neither imports a producer. The compact receipt retains
instructions and counts, not a large expanded fault tape. Verification under
`python -O`, producer-disabled replay, malformed inputs and altered physical
circuits are regression tests. Illustrative noise constants are never called
a measured or derived full-machine threshold.

The archive increases the safe active location bound to O(q^5 log^16 q).
Retained passive protection diagnostics bring total physical storage-time
volume to O(q^6 log^25 q); they do not feed back into the computation.
Five constant quantum protection levels and logarithmic
classical words yield joint operational error O(q^-20 log^16 q), with paid
cofinal controls. Bare terminal bits, unprotected external controllers and
unknown raw-input encoders do not inherit that bound. The fixed source and
external circuit clock remain specified capabilities, not A1--A3 deductions.
