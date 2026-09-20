# Native sum records and finite-chain cleanup

This packet adds a native write to the captured paired-rail memory of #917
and proves cleanup bounds for arbitrary finite chain length. It does **not**
finish M1. The exact objective and exit conditions are in [CONTRACT.md](CONTRACT.md).

## A new record is computed, then reused

The ideal raw values of a balanced pair are `(2+a, 2-a)`. Pair means on the
two corresponding rails implement `(a,b) -> ((a+b)/2,(a+b)/2)`; averaging
the two rails of one cell clears its amplitude. These are the operations
proved in `SourceEncodedMemory.lean`, with the composition and arbitrary
per-step error theorem in `SourceReusableBus.lean`.

Begin with amplitudes `(A,B,0,0)` in two archives and two blank cells:

1. Export A to the first blank by two means: both become A/2.
2. Export B to the second blank by two means: both become B/2.
3. Average the two new pairs: both become (A+B)/4, using two means.
4. Clear the first temporary pair by one mean.

The result is `(A/2,B/2,0,(A+B)/4)`. Publication after those seven means
declares the last cell to be record C, with decoded value A+B and scale
exponent 2. A and B have scale exponent 1. No preparation write, host-side
addition, protected export, assigned reset or factor-two physical amplifier
creates C. Publication and the scale/version metadata are supplied control.
The decoder only classifies its own two receiver samples against the
specified codebook and exponent; it does not receive A or B.

The injective embedding uses these actual captured ports:

| Role | Positive, negative ports |
| --- | --- |
| A | 0, 9 |
| B | 7, 4 |
| C | 8, 11 |
| Bus 0 | 1, 5 |
| Bus 1 | 14, 40 |
| Bus 2 | 23, 38 |
| Bus 3 / receiver | 45, 39 |

The verifier checks every scalar mean against
`code/source_routing/support_w12_l3.json`. C occupies the second blank used
in the write; Bus 0 is the first blank and is cleared before transport.
The same bus then serves `C,A,C,B`, cleaning for 80 sweeps after each read.
The sources are on carrier 0 and the receiver on neighbouring carrier 3;
this experiment is not a long-distance routing demonstration.

For all 16 pairs of A,B from `{-3/2,-1/2,1/2,3/2}`, all four reads decode
correctly, including the second C read after A has used the same bus.
Changes to B leave decoded A unchanged, changes to A leave decoded B
unchanged, and different input pairs with the same sum give the same decoded
C. These comparisons are exhaustive for this finite input set. The ideal
linear maps are independently computed for arbitrary real input amplitudes.

The receiver exponents for the four reads are `6,5,7,5`. Half spacings are
`1/128,1/64,1/256,1/64`; every certified margin is positive. The error budget
includes the ideal cleanup residual `(3/2)*(7/8)^80`, initial coordinate
error `2^-18`, `K*(2^-28+1/(2*2^20))` for K preceding means, and readout
error `2^-16`. Those three extra error allowances are hypotheses, not
measured hardware guarantees. The retained histories execute nearest-even
fixed-point means at Q=`2^20`.

Raw causal ancestry still includes the other records and cleanup operations.
This is a decoded functional correspondence, not equality of physical causal
orders or of all induced writer relations. No intermediate operation is
deleted from the replay or cost.

## Polynomial cleanup on a supplied chain

For n edges and n+1 bus cells, clear cell 0 and sweep paired copies from
0 to n. Only the first cell needs a reset rung; every hop needs two rail
edges. The scalar cost of one sweep is `1+2n`. Cells outside the bus are
unchanged. The zero-edge bus is cleared exactly in one mean.

For n >= 1 put

```
M = (n+1)^2
w_i = M - (n-i)^2,    0 <= i <= n
beta = 1 - 1/M
```

The actual sweep maps the nonnegative weight vector to at most `beta*w`.
The proof uses `w_(j-1)+w_(j+1)=2*w_j-2`, the first-cell estimate, and the
last-cell duplication. Positivity and the triangle inequality give the
same weighted bound for arbitrary signed inputs. Since `min(w)=2n+1` and
`max(w)=M`, initial amplitudes bounded by A have residual at most

```
A * M/(2n+1) * beta^R
```

after R sweeps. Lean also proves `beta^r <= M/(M+r)`, hence each block of
M sweeps reduces the weighted norm by at least one half. With `R=M*k`:

```
residual <= A * M/(2n+1) * 2^(-k)
scalar means = (1+2n) * (n+1)^2 * k
```

This gives polynomial work in chain length and requested binary cleanup
accuracy. It replaces an exponentially loose sufficient-sweep estimate
obtained by generalizing the one-sweep unweighted bound. For the existing
three-edge bus, the older specialized `(7/8)^80` bound remains sharper and
is still used in the captured experiment.

This theorem does not supply the rail embedding in a general captured graph.
Forward transport still attenuates its signal by `2^-n`, and repeated export
further attenuates the source. The general forward-read error theorem and
`program_readout_error` retain those gains and absolute error terms. Increasing
cleanup also increases the accumulated noise allowance. Unbounded reuse at
fixed grid and fixed nonzero noise is not established.

## Two precise boundaries

`raw_growth_obstruction` proves that every finite word of scalar means stays
inside its initial uniform scalar range, including words selected adaptively.
An unscaled sequence growing by at least one per layer eventually leaves
that range. In particular, nonnegative `1+sum` updates with a self-read cannot
run indefinitely as raw values on a fixed finite closed scalar store.
Scaled encodings (including this packet), finite horizons, fresh inputs and
growing representations are outside that obstruction. It is not an
impossibility theorem for M1 or logical memory.

The captured path `0--1--14` with loads `(1,3,5)` admits two supported,
strictly quadratic-energy-decreasing words using the same two edges:
`(0,1);(1,14)` reads 2 at port 0, whereas `(1,14);(0,1)` reads 5/2 there,
with different consumed-input ancestry. This is a counterexample to selection
by the candidate local mean law alone. It does not verify the complete A1--A3
clauses or their accepted-record domain and is **not** an A1--A3 countermodel.

The canonical [axiom reference](../../docs/AXIOM_REFERENCE.md) requires an
A1-generated state/observable space, complete A2-visible constraints, an
exact reference and a map from the A3 optimizer to the claimed output.
Those objects and that selection theorem have not been supplied for the
metric-neighbour read controller. `CoreAxioms.lean` explicitly contains typed
shadows; inhabiting them would not settle this missing derivation.

## Evidence and reproduction

Run from the repository root:

```
python code/source_native_updates/verify.py
python -m pytest -q code/source_native_updates
cd Lean
lake build Geometry.SourceNativeUpdatesAxiomAudit
```

`build.py` expands the write and transport controller independently of the
verifier. `controls.json` retains one full 2,279-event representative trace
and the initial, committed, read, final and trace-commitment data for all
16 input pairs. Each event records its endpoints, consumed values, consumed
writer IDs and output. Scalar arrays occupy one line. Every verification
reconstructs and checks all 36,464 events; the fifteen compact histories are
not accepted on their hashes alone. The verifier never imports the producer.
`receipt.json` contains the independently derived read bounds, ancestry,
resources, chain controls and schedule witness.

The modeled per-history resources include 14 active scalar registers, three
logical records after commit, 2,279 means (972 crossing carriers), 4,558
native scalar reads and writes each, 14 preparation writes, and eight local
receiver samples. The receipt separately counts register width, mean and
receiver arithmetic, expanded schedule, port map, counters, record IDs,
codebook entries and publication flag. These are explicit representation
bounds, not Python memory measurements. They exclude simulator/evidence
storage, immutable controller code, physical isolation costs and the inactive
ports of the 15,360-port support; those exclusions prevent treating 14
registers as the capacity of the whole physical realization.

Regenerate with `build.py`, then `verify.py --write-receipt`. Source pins and
strict canonical bytes bind the result. Production and verification use
different control expansion, rounding and decoding implementations. Mutations
are tested after resealing their commitments, including a different valid
native program that computes the same sum.

## M1 work still owed

| Obligation | This packet | Remaining work |
| --- | --- | --- |
| Native record service | One new sum version and preservation/rereads of old logical values on captured seams | Arbitrary finite record programs, bias generation, many live versions, allocation and commit/abort policy |
| Scaling | Any supplied finite chain has a polynomial cleanup bound and explicit forward attenuation | General captured rail embeddings, realistic precision/noise/capacity for complete histories |
| Full-family compilation | Sixteen finite write/read episodes | Complete q=13/q=21 native executions and intervention-preserving operational refinement with all operations counted |
| Source selection | Exact ambiguity under the candidate mean law; restricted raw-growth obstruction | Derive the controller/read law from the full axioms, or supply a full qualifying non-derivability result |

The new finite construction removes the need to prepare **every** record in
advance. It does not remove preparation altogether or prove a complete native
service. No supplied-remainder closure of #777 or transfer of the M1
derivation obligation is taken.
