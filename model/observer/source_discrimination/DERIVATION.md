# Source implementation of the complete C/E discriminator

The [existing prediction protocol](../../predictions/DERIVATION.md) specifies
two recorded repair channels on all 78 total-load-two configurations, four
preparations and one binary source measurement. This construction supplies
finite source words for those operations in an explicit code. It preserves
the old candidate definitions and fixed statistical protocol. The code's
identification with physical load states, coherent proper-code services,
reset, classical records and pointer measurement remain supplied apparatus.

## Full carrier and independently specified repair

Write a configuration as `(a,b,s)`, where the first two entries are the
chosen seam and `s` lists the other ten ports. The seam has six possible
states, encoded in the low three source bits as

| Endpoint pair | `00` | `01` | `10` | `02` | `11` | `20` |
| --- | --- | --- | --- | --- | --- | --- |
| Code | 0 | 1 | 2 | 3 | 4 | 5 |

The remaining six bits label the spectators. Their total is `2-a-b`, so
labels may be reused between spectator totals without identifying distinct
native states. There are 55 spectator states at total two, ten at total
one and one at total zero. The number of native code states is therefore
`55 + 2*10 + 3*1 = 78`. The nine-qubit register has 434 unused addresses.
The fixed enumeration puts `X={2,3}` at address zero and `Z={1,2}` at
address one; `W={0,2}` is at address two. `encode` is injective and `decode`
recovers every native configuration. Denote this isometry by `J`.

The independent target is endpoint averaging. For endpoint total `m=a+b`,
let `l=floor(m/2)` and `h=ceil(m/2)`. Coin zero selects `(l,h)` and coin one
selects `(h,l)`, with unchanged spectators. Call the result `F_c(x)`.
The current is the transferred load `j=a-F_c(x)_0`; the contrast is the
original difference `d=a-b`. This characterization does not evaluate the
prediction oracle's current formula.

Prepare and measure a fresh fair coin first, retaining its outcome `c`.
For C the current record codes are `0->0`, `+1->1`, `-1->2`. With blank
record bits above the three endpoint bits, the local basis index is
`endpoint+8*current_code`. The following disjoint transpositions define
the entire five-qubit local permutation; all other basis states are fixed.

| Coin | C local transpositions |
| --- | --- |
| 0 | `(2,9) (3,20) (5,12)` |
| 1 | `(1,18) (3,20) (5,12)` |

For every admitted blank-record input this outputs precisely `|F_c(x)>|j>`.
E applies the same C word with one additional blank record bit, then one
CNOT from endpoint bit `c` to local bit five. After C, even seam total has
endpoint code zero or four; odd seam total has code one for coin zero or
code two for coin one. The copied bit is therefore exactly seam-total
parity `p`. E's raw record code is `current_code+4*p`. Its contrast decoder
is `d=2*j+(2*c-1)*p`, with the admitted values

| Coin | Raw record code to contrast |
| --- | --- |
| 0 | `0->0, 1->2, 2->-2, 4->-1, 5->1` |
| 1 | `0->0, 1->2, 2->-2, 4->1, 6->-1` |

C uses five local qubits and E uses six. Spectator bits are untouched.
Both words define unitaries on their entire local registers, although the
candidate instruments require blank record initialization. All admitted
record outcomes are retained; relabeling `(c,j,p)` as `(c,d)` is bijective.
The remaining record codes are padding, with zero probability on the ideal
native input code. Failure handling must retain any actual invalid outcome.

## Exact recorded-channel intertwining

Each input column of the joint coin/record/source isometry has two
coefficients `1/sqrt(2)`. For C, a fixed current translates the endpoints
injectively. For E, fixed contrast and unchanged spectators determine the
original endpoints, so its record/source map is also injective. Coin is
part of the record. Thus distinct input columns have disjoint joint support
and every column has squared norm one: `V*V=I` on all 78 inputs.

For any matrix unit `|x><z|`, recording the labels gives

```text
E_C(|x><z|) = (1/2) sum_{c: j(x,c)=j(z,c)}
    |c,j(x,c)><c,j(x,c)| tensor |F_c(x)><F_c(z)|,

E_E(|x><z|) = (1/2) sum_{c: d(x)=d(z)}
    |c,d(x)><c,d(x)| tensor |F_c(x)><F_c(z)|.
```

`certificate.channel_certificate` independently enumerates the 78 integer
configurations, derives floor/ceiling outputs, and checks the permutation
image and all `78^2=6084` matrix units for each candidate. It also checks
the complete Gram identity, with 156 nonzero Kraus entries. Consequently

```text
E_code(J rho J*) = (I_record tensor J) E_native(rho) (I_record tensor J*).
```

This identity holds for arbitrary operators, not only basis populations or
the four experimental preparations. Expanding a source/reference operator
as `sum_{x,z} |x><z| tensor R_xz` proves the same identity with any untouched
reference system. It includes all unnormalized terminal record blocks.
It does not derive a physical encoding of a pre-existing load carrier.

## Contrast is the minimal binary refinement of the current record

Let `Pi_p` project onto seam-total parity `p=(a+b) mod 2`. Endpoint averaging
preserves `a+b`, so this parity can be measured before or after repair. Its
relation to current and contrast is an identity of the declared rule:

```text
d = 2*j + (2*c-1)*p,
p = abs(d) mod 2,
j = (d-(2*c-1)*p)/2.
```

For fixed `c`, the admitted pair `(j,p)` and the contrast `d` therefore
determine each other. The fair-coin factor `1/sqrt(2)` is unchanged.
If `L_(c,j)` is the current-record Kraus operator, the parity-refined
operator is `Pi_p L_(c,j)`, whose support is precisely the contrast group
with the displayed `d`. Hence the entire E instrument is C followed by
recording seam parity and this classical relabeling. This is an equality of
complete unnormalized quantum branches, not only outcome probabilities.

On a matrix unit `|x><z|`, C retains a coin branch exactly when the two
currents agree; the added parity record retains it exactly when the two
parities also agree. By the identity above those two conditions are
equivalent to equal contrasts. `parity_refinement_certificate` checks all
6084 units, including the 2650 connecting different spectator-total shells.
For example, the coherence between `X={2,3}` (seam total zero) and `{0,1}`
(seam total two) survives: both have even parity and zero current. No
measurement of the spectator shell is introduced. Linearity extends the
identity to any reference entangled across these shells.

For a fixed C record, at most two distinct E records occur, indexed by
`p=0,1`. Two are necessary: `X={2,3}` and `Z={1,2}` share C's `(c,j)=(0,0)`
record but give E contrasts zero and minus one. Thus an additional
one-valued record could not distinguish the required branches. The maximum
refinement fiber has size two, proving that one additional binary classical
outcome is necessary and sufficient when the existing C records are
retained. This statement concerns the specified classical refinement. It
does not claim minimum physical memory, environment dimension or gate count;
the joint labels can be encoded in other alphabets.

Forgetting the added parity label gives, inside each old current block,

```text
rho_(c,j) -> Pi_0 rho_(c,j) Pi_0 + Pi_1 rho_(c,j) Pi_1.
```

It leaves parity dephasing and does not recover C. In particular the
`X,Z` off-diagonal survives C's coin-zero/current-zero branch and is lost
under E. The fixed coherent phase/readout test consequently measures the
effect of this additional parity recording, while its diagonal control
has no such coherence to lose. This factorization introduces no new
candidate hypothesis: it follows from the existing C and E definitions.
Whether an identified physical apparatus records that extra parity remains
an operational choice or a further physical selection premise.

## Source compilation, phase and complete readout

A controlled bit flip is compiled from its matching-controls projector
`P`. Since `P` commutes with the target Pauli `X`, the exact identity is

```text
controlled_X = exp(-i*pi*P/2) exp(i*pi*P*X/2).
```

The nonmatching subspace sees identity; on the matching subspace the two
factors are `-i` and `iX`, giving exactly `X`. Walsh expansion of `P` makes
both factors products of commuting tensor-Pauli rotations. The scalar
Walsh term is retained for exact unitary equality. Omitting only that term
would change a common global phase, which cancels from the standalone
instrument. Relative column phases can change the instrument. Full-unitary
tests enforce the declared phase convention; terminal-density tests accept
common global phases and reject relative-phase corruption. The existing
source compiler realizes every Pauli rotation
using its full six-level generators and declared proper-code services.

To transpose two arbitrary local basis states, choose a path of single-bit
flips between them. Apply its edge transpositions in order and then all but
the last edge in reverse. The endpoints exchange with phase `+1`; every
intermediate and unrelated basis state returns to itself. This compiles the
local permutations above without a primitive many-qubit swap. Zero-angle
and scalar-phase requests remain in the resource schedule.
E then appends the existing phase-exact native CNOT, costing three more
source pulses: 5370 for either C repair word and 5373 for either E word.
The parity record's initialization and readout retain their service costs.

Let `U=exp(-i*pi*Y_0/4)`, so `U|0...0>=(|0>+|1>)/sqrt(2)`. The preparation
then applies the bit-zero phase `diag(1,exp(i*theta))`, compiled as an
identity phase followed by a Z rotation. For `theta=0,pi,pi/2`, the prepared
native state is exactly the corresponding `(X+exp(i*theta)Z)/sqrt(2)`.
The diagonal control uses an additional independent fair randomizer and a
recorded mixture of identity and Z. Both randomizer outcomes and all
attempts remain present. Summing their unnormalized blocks gives
`(|X><X|+|Z><Z|)/2`; neither outcome is selected or renormalized away.

For readout, initialize a pointer to zero, apply `U*`, flip the pointer
only at the all-zero nine-bit source address, and apply `U` again. Measuring
only this pointer gives source Kraus operators

```text
P_plus = U |0...0><0...0| U* = |X+Z><X+Z|/2,
P_minus = I-P_plus.
```

The pointer value one denotes plus. The source basis change is restored
before its measurement. In particular the minus outcome is one rank-77
projector on the native carrier, not a measurement of its separate basis
states. Coherences throughout that space survive. The exact checker verifies
both projector columns on all 512 encoded addresses and trace preservation
on all 6084 native matrix units. Applying the same pointer channel separately
inside every old record block preserves all old record probabilities.
Classically grouping E's contrast records into the old common current
interface sums density operators and never restores the lost coherence.

These identities recover the existing predictions
`p_C(theta)=3/8+cos(theta)/4`, `p_E(theta)=3/8`, and the diagonal-control
probability `3/8` for both. Source execution is an independent numerical
implementation check of this encoded construction. Its floating residuals
are not used as bounds in the exact channel argument.

## A constructive limit on selection by the source alphabet

The two realizations give a counterexample to selecting a unique quantum
repair from the available control alphabet and the classical endpoint law.
They use the same twelve generators, nine-bit source code, placement coin,
proper-code services, fixed preparation and final pointer measurement.
If a common total carrier is required, pad C's two-bit repair record with
one untouched blank bit; both then fit E's three-bit record register.
Both channels conserve total load, leave all spectator occupations fixed
and induce exactly the same classical transition for every native basis
state. Yet the common plus observable at phase zero differs by `1/4`.

Therefore availability of these source controls, code services and
classical transitions does not determine which quantum record instrument
is executed. The missing choice is operationally measurable: retaining
contrast distinguishes inputs whose current record coincides and destroys
their interference. This is an explicit pair of implemented alternatives,
not a claim that both satisfy every premise of the theorem library, nor an obstruction
to selection by an additional law governing admissible records. The
frozen count test probes precisely that remaining choice once its physical
preparation, readout and error premises are independently established.

## Compositional apparatus contract

Use trace distance `D(rho,sigma)=||rho-sigma||_1/2` and channel distance
`d(Phi,Psi)=||Phi-Psi||_diamond/2`, on complete channels with arbitrary
reference systems. If an actual pulse and its ideal unitary satisfy
`||U_actual-U_ideal|| <= u`, their channel distance is at most `u`.
This follows by applying the vector bound to a purification and then
contracting trace distance; it includes global and relative phase errors
without relying on a fitted execution residual.

Let `N_p,N_l,N_e` be upper bounds on native pulse, proper-code load and
proper-code export invocations on every complete attempted-trial branch.
They include the emitted fair coin, the diagonal randomizer where present,
zero-angle requests and the full binary pointer word. Let `u,l,e` bound
the pulse operator error and the half-diamond error of each complete load
and export channel. Leakage and service failures belong to those channels;
they cannot be deleted or divided out.

The remaining whole-trial error allowances are:

- `p`: preparation and supplied encoding error beyond already charged words;
- `z`: resets and blank initialization;
- `c`: fair-coin generation/recording errors beyond charged pulses/services;
- `a`: the diagonal randomizer's remaining generation/recording error;
- `r`: other record acquisition, storage and classical control errors;
- `s`: source storage or transport not charged above;
- `m`: final pointer measurement and output/abort assignment;
- `g`: an independently bounded difference between exact mathematical
  coefficients and the compiled finite coefficients;
- `n`: a certified numerical execution error, if a numerical result is
  itself being bounded.

Each term is an aggregate upper bound for its stated operations. It must
include arbitrary-reference behavior for channels that act on quantum
inputs. Preparation of a fixed state may instead use its trace-distance
bound, since later common channels contract it. The partition must account
for every operation; any known overlap may be charged twice conservatively.
No term is inferred from the comparison counts.

Telescoping the actual and ideal channel compositions and using contraction
gives the sufficient condition

```text
epsilon = N_p*u + N_l*l + N_e*e + p+z+c+a+r+s+m+g+n <= 1/64.
```

For classical branches, use the maximum complete-branch resource counts and
uniform channel bounds at every allowed record value. Controlled subsequent
operations and record-preserving grouping are common channels. Thus the same
bound controls the complete final classical/quantum instrument, including
all attempted-trial failures. It requires no normalized branch comparison
and never divides by a possibly small outcome probability.

Every joint-bin event, the aggregate plus event and the abort event are
effects with operator norm at most one. Their probability differences are
therefore each at most `min(1,epsilon)`. In particular the aggregate plus
bound is the same epsilon, without multiplying it by the number of bins.
If the displayed sufficient condition holds, it supplies the existing
protocol's unchanged `eta=1/64` premise simultaneously for every setting.
Its 8192 attempts per setting, four settings, simultaneous intervals and
decision rule remain unchanged. This argument does not establish fresh-trial
independence or calibration; those remain the fixed protocol's premises.

`apparatus_error_contract` requires all rational error terms explicitly and
checks this arithmetic implication. It supplies no physical values for
them and declares neither physical calibration nor certified numerical
execution. An ideal symbolic construction can set these errors to zero as
a mathematical case; that assignment does not certify actual services.
The full source-coded apparatus is a completed encoded construction, with
physical identification and independently established error bounds still
required for empirical use.
