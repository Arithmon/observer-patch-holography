# Native observer records reuse the existing capacity construction

The finite observer's record budget can
be expressed in exactly the mathematical language already used by OPH's
cosmological-capacity work. This packet makes that connection and executes a
smaller classical controller representation. It does not identify this local
budget with the capacity of a cosmological horizon.

## The existing equations apply directly to the records

The current [public-record closure](../../paper/tex_fragments/OBSERVERS_N_CLOSURE_PRO5.tex)
defines `M0` as the largest correctable set of reachable public records:

```
M0 = alpha(compound checkpoint confusability graph),
N_record = ln M0  (nats),
M0 <= D,  where D is the declared carrier dimension.
```

The [existing evaluator](../../code/capacity_readback/correctable_public_record_capacity.py)
already implements the required record atoms, reachability, joint checkpoint
channels, graph and carrier bound. `capacity.py` imports it directly,
with source pins covering its local dependencies. No new capacity definition
or graph solver is introduced.

The [native observer](README.md) produces
11 distinct terminal two-letter records from its declared uniform binary
66-state preparation. All 264 complete histories are retained. The full
78-state total-two input family produces 17 records from 312 histories.
The atom algebra has all 36 allocated two-letter words; only words with an
actual native history enter the reachable subset.

One observer reads its full classical tape. The declared terminal
continuations are halt ticks and exact classical rereads. They leave that tape
unchanged, so their joint kernel is identity, its compound graph has no
edges, and the existing theorem gives:

| Input family | Reachable first records | Reachable terminal records / `M0` | Terminal capacity |
| --- | ---: | ---: | --- |
| Declared binary benchmark | 4 | 11 | `ln 11` nats, about 3.459 bits |
| All 78 total-two configurations | 6 | 17 | `ln 17` nats, about 4.087 bits |

The continuation family is complete for this explicitly halting program.
It is not a census of the complete OPH terminal fiber. A single observer's
joint kernel equals its one local marginal; no multi-observer coupling is
inferred from separate marginals. Actual Shannon entropy under the benchmark
record probabilities is another quantity, about 2.023 nats, and need not equal
the maximal code capacity `ln 11`.

The distinct dimension counts are useful:

| Quantity | Value | Meaning |
| --- | ---: | --- |
| Record slots | 2 | Maximum number of writes in this program |
| Terminal tape dimension | 36 | Six letters in each allocated slot |
| Terminal source plus tape sector | 2808 | `78 * 36` |
| All-phase carrier dimension | 3354 | `78 * (1 + 6 + 36)` |
| Reachable benchmark correctable record count | 11 | The code size actually returned by the existing evaluator |

The carrier certificate uses the real full 3354-dimensional allocation.
For terminal word index `w`, its projector occupies basis coordinates
`546 + 78*w` through `546 + 78*w + 77`; each projector has rank 78.
Different words are orthogonal. All 36 terminal projectors resolve the
terminal-sector identity, not the identity on the earlier pointer sectors.
The 11 reachable projections have total rank 858. The result is unsaturated,
with slack `ln(3354/11)`: unused allocation is not extra readable information.
The adapter never replaces this with an artificial rank-one carrier of size
11 in order to manufacture closure.

## A concrete simplification: four controller bits on the benchmark

The benchmark reaches exactly `1 + 4 + 11 = 16` classical record prefixes:
one empty prefix, four one-letter prefixes and eleven two-letter prefixes.
`results/capacity.json` contains their explicit encode/decode and append tables. Four
bits therefore suffice for losslessly retaining this finite controller
history. This is the minimum fixed-width encoding of all 16 prefixes, not a
minimum-state theorem for control policies allowed to forget old records.
Its code already determines the phase, next seam and supplied forecast.

`compact_view` reads only that code. `compact_tick` attaches it to the
unchanged source state, executes the same native event and returns the new
code. Every one of the 264 benchmark histories is checked against the original
execution, including phase, action, forecast, source loads, record and halt.
An unrecognized code or a controller/source pair outside the reached support
is rejected. The wider 78-input family has 24 prefixes and needs five bits.

This compares with eight bits for a naive separate two-bit phase and two
three-bit slots, or six bits for an optimal index over all 43 allocated
prefixes. It compresses the reached **classical controller representation**;
it supplies no amplifier, physical memory device or quantum compression map.
The 3354-dimensional quantum construction and its arbitrary-input isometry
are unchanged. Program tables, source state and implementation workspace are
not included in the four controller bits.

Keeping only predictions would lose records. Reusing the same capacity
evaluator on declared diagnostic coarse readouts gives:

| Retained terminal information | Correctable code size |
| --- | ---: |
| Full two-letter tape | 11 |
| First letter alone | 4 |
| Second letter alone | 4 |
| Forecast alone | 3 |
| Chosen second seam alone | 2 |
| Forecast and chosen seam | 4 |
| Complete erasure | 1 |

These channels are controls that discard distinctions, not newly derived
native operations. If erasure is admitted as a continuation alongside halt,
the compound zero-error capacity is one. Positive full-support noise can do
the same, however small its probability; a noisy device needs the existing
approximate-capacity definition with a specified error and horizon.

## The precise connection to Lambda

The existing global proposal adds a much stronger closure condition:

```
N_star = ln D_star = ln M0(U_N_star),
F_set,r,0(D_star) = {D_star}.
```

It asks for the same direct capacity across the complete terminal fiber and a
law selecting the physical carrier. If that same carrier is then identified
with a de Sitter horizon's entropy, the existing geometric dictionary is

```
N_star = A_dS / (4 ell_star^2),
A_dS = 12*pi / Lambda,
therefore Lambda * ell_star^2 = 3*pi / N_star.
```

`N_star` is measured in nats, `Lambda` in inverse area, and `ell_star` is the
length fixed by the corresponding scale identification. Expressing it as
the Planck length additionally uses that identification. Bits would introduce
the factor `ln 2`. Two recorder slots do not mean `N_star=2`, and raw Hilbert
dimension 3354 does not mean `N_star=3354`.

The useful link is constructive: native records can now supply an actual
finite input to the same capacity evaluator. The cosmological step still
requires a capacity-indexed source family, a complete-fiber closure and
physical selector, and the same-carrier horizon/scale attachment. This
packet makes none of those identifications and computes no Lambda. Existing
[fiber-product sewing](../../paper/tex_fragments/OBSERVERS_N_CLOSURE_PRO5.tex)
also provides the next finite-network counting rule:
`|X_(A union B)| = sum_z |r_A^-1(z)| |r_B^-1(z)|`.
It avoids multiplying local record counts when shared records impose
compatibility.

## Reproduce

From the repository root:

```sh
python3 -m model.observer capacity --check model/observer/results/capacity.json
python3 -m unittest model.observer.test_observer_capacity -v
```

The independent tests reconstruct floor/ceiling events, all reachable words,
capacity changes and all compact-controller paths. They reject overlapping
carrier projectors and unsupported compact inputs. `--write` creates a new
receipt exclusively and refuses to overwrite an existing result.
