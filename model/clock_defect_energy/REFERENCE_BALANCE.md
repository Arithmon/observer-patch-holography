# Fixed clock-reference balance under bare C/E repair

The all-core clock service transports a faithful reference along each
recorded branch. An energy/action interpretation also needs to say which
measure the unconditioned native repair law preserves. This packet tests that
question on the same binary H66 source, before adding a new repair or memory
primitive.

The [exact receipt](reference_balance_results.json),
[source](reference_balance.py) and [independent Kraus check](test_reference_balance.py)
retain every C and E seam/coin output, including zero-current events. The
population kernel is symmetric: each of the 300 native movement edges has
directed probability `1/60` per attempt. The connected graph therefore has a
unique normalized reversible measure, uniform on 66 configurations. This is
an actual-model instance of
[`normalized_symmetrizer_unique`](../../Lean/Geometry/SourceActionMeasure.lean).

Holding the transported clock reference **fixed** during unconditioned
repair gives the following exact results:

| Core | One-attempt total variation drift | Edges violating detailed balance | Positive mode form / variance |
| --- | --- | --- | --- |
| `{0,5}` | `17/750` | `80/300` | `710/6983` |
| `{0,8}` | `17/600` | `100/300` | `8515/83796` |

The mode is the same native current observable used in
[MODE_ENERGY.md](MODE_ENERGY.md), with port profile
`(1,-1,-1,1,0,0,0,0,0,0,0,0)`. Its fixed-reference quadratic change stays
positive, but the reference drifts. The two ratios differ from the uniform
stationary value `1/10`. Consequently, a proof using the bare C/E reversible
action measure cannot silently use the fixed transported clock reference as
its equilibrium measure. A common physical action would need an explicitly
changed law, a co-moving reference/history construction with its own action,
or a separately justified dictionary. The calculation does not select one.

The predecessor's branchwise covariance remains intact: it updates the core
metadata, transports the reference and retains each record. This test holds
one reference fixed while averaging all outcomes, so it is a separate
compatibility check. No clock calibration, particle rest energy or physical
stress follows from either measure alone.

## A stationary co-moving source and recorded core

There is a positive joint result under **uniform seam selection**. Include the
already counted classical two-port core label `c` and update it by the same
recorded seam transposition as the source pair `x` whenever current is nonzero.
On the `66 × 66 = 4,356` source/core states, set

```text
pi(x,c) = sigma_c(x) / 66,
q(x,c) = number of occupied source ports in c.
```

Every movement and its reverse have the same native probability `1/60` and
preserve `q`, hence preserve `pi`. Zero-current events are self-loops. The
receipt checks all `300 × 66 = 19,800` movement/core pairs, both uniform
marginals, and every partial-injection pair of the complete C/E record maps.
Each raw record block preserves the three `q` projectors on this enlarged
source-plus-controller domain. Their stationary weights are exactly `2/5`,
`5/9` and `2/45`. For a source-only current observable, the joint source
marginal recovers the uniform form `D(v)=2/33` for the chosen mode.

This gives an exact finite way for the transported clock reference and native
current action **measure** to coexist on one recorded source/core process.
The core is counted classical metadata, so `q` is a relational controller
charge, not a new autonomous source charge. Its conservation does not keep
the occupied pair spatially localized, select a fermionic exchange phase or
derive physical rest energy. Stationarity under adaptive feedback or a
laboratory clock has not been shown.

```sh
python -m model.clock_defect_energy.reference_balance --check
python -m pytest -q -W error model/clock_defect_energy/test_reference_balance.py
```
