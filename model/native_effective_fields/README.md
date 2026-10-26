# Native exchange fields and measured four-point correlations

The existing native observer supplies a concrete finite hard-core exchange
field. Its coefficients follow the native seam rule, and its recorded currents
have nonzero connected four-point correlations. This packet identifies that
field, bounds its transfer Hamiltonian by exact rational polynomials, and
executes three complete four-read experiments. No fitted oscillator, Gaussian
field assumption or independent replica limit enters these predictions.

This is a field on the supplied twelve-port carrier, in its binary total-two
sector. It does not establish a spatial continuum, a number of physical spatial
dimensions, relativistic dynamics, or agreement with a laboratory dataset.

## Which field follows from the native rule

The 66 basis configurations have occupations `x_i in {0,1}` and
`sum_i x_i=2`. Define `n_i |x> = x_i |x>`. For each source seam `(a,b)`,
the following operator preserves this sector:

```text
G_ab = n_a + n_b - 2 n_a n_b
       - (raise_a lower_b + raise_b lower_a).
G = sum over the 30 source seams of G_ab.
K = I - G/60.
```

Here `raise` and `lower` are the hard-core occupation ladder operators on a
port. `G_ab` is exactly identity minus the endpoint-swap permutation. It is
positive with eigenvalues zero and two. Disjoint seam terms commute; terms
sharing a port need not commute. The implementation compares every entry of
`G` against the independently enumerated native transition counts. Thus `K`
is precisely the existing one-attempt transfer, with no choice of coefficients
left in this identification.

The density interaction matters at fixed hopping. Deleting `-2 n_a n_b`
changes `G` by twice the adjacent-occupation indicator. The constant vector
ceases to be a ground vector, and the resulting formal step `I-G_deleted/60`
takes the uniform probability law to total mass **65/66**. It therefore is
not a stochastic replacement for the native rule. This is an operator-level
exclusion interaction; a nonzero fourth cumulant alone would not establish
quantum interaction.

The scheduler and sector remain premises: uniform independent selection of
the thirty seams, the native fair placement coin, and the binary total-two
carrier. These results do not transfer unchanged to the production observer's
adaptive schedule.

## A controlled effective Hamiltonian, including its longer-range terms

The existing [finite transfer construction](../primitive_observer/DERIVATION_TRANSFER.md)
establishes `2I/3 <= K <= I`, full-history bond reflection positivity, and
`H=-log K`. The present addition supplies exact computable approximants:

```text
D = G/60,                 0 <= D <= I/3,
H_m = sum_(k=1..m) D^k/k,
0 <= H-H_m <= r_m I,      r_m = 1/[2(m+1)3^m].
```

The bound follows by the positive scalar logarithm series and spectral
calculus: replace each `1/k` in the tail by `1/(m+1)` and sum its geometric
majorant. All entries of `H_m` are computed as exact fractions. For real
unitary parameter `t`, the operator error is at most `abs(t) r_m`; for
nonnegative Euclidean parameter `s`, the semigroup error is at most `s r_m`.
These are analytic errors of the exact polynomials, not certificates for
later floating matrix products or physical services.

At order sixteen, `r_16=1/1463588514`. The receipt gives an exact positive
interval for `H[0,2]` although `G[0,2]=0`. Thus the logarithm contains
longer-range configuration couplings. The original seam generator is local
on the declared port graph; its exact discrete-attempt logarithm is not
restricted to single seam exchanges. Integer Euclidean gaps satisfy
`K^n=exp(-n H)`. The existing negative-entry witness for `sqrt(K)` prevents
treating every fractional Euclidean step as a stochastic native step.

The same transfer supplies a faithful mathematical reference `sigma_K=K/56`
on this 66-dimensional **configuration transfer Hilbert space**. Consequently
`-log(sigma_K)=H+log(56)I`, so its modular automorphisms use this Hamiltonian.
With the convention `A -> sigma_K^(is) A sigma_K^(-is)`, the dual state
action is `rho -> exp(isH) rho exp(-isH)`. This is not the native stationary
density `I66/66`, an executed native preparation of `sigma_K`, or a proof
that this reference is uniquely selected. Unitary continuation remains
distinct from native completely positive repair. Native attempt count,
Euclidean parameter, modular parameter and a calibrated physical clock are
not identified.

## Finite spatial support recovered from observer records

The twelve-port response calibration reconstructs an exact rank-three unit
Gram matrix from retained native currents; its edge set is the same thirty
seams used above. In that recovered metric, each source seam has squared
chord length `2-2 sqrt(5)/5`, and every local exchange term `G_ab` is exactly
identity minus the port swap on one such shortest-shell edge. The executable
checks all thirty terms and their exact sum `G`. Thus the hard-core field's
finite generator is local on the record-reconstructed twelve-point metric.

The logarithm `H=-log K` need not remain local in that metric: its certified
longer-range witness is already reported above. The rank-three metric is
unitless and finite; it does not supply calibrated lengths, an extended
spatial manifold, or a continuum limit. The calibration uses a separately
declared preparation in the same total-two sector, rather than the stationary
current-read preparation used by the four-point experiments. See the
[record reconstruction](../relative_observer_geometry/README.md) for the
complete shell census and relabeling control.

## Actual current-only measurement predictions

At each of four successive native attempts, the observer retains the signed
change of occupation at the declared read port. It obtains that change from
the selected seam and its committed current; it never measures the hidden
configuration. For seam `(a,b)` with signed current `j`, the record at port
`i` is `delta n_i=-(1[i=a]-1[i=b])j`, taking values `-1,0,1`.

The preparation is the declared stationary mixture `I66/66`, without disclosure
of the actual starting configuration. The [producer](construction.py)
constructs the existing rule-C Kraus maps, including their measurement
disturbance. The shared observer runtime executes four complete instruments.
Every one of the **81 record words**, including exact zero outcomes, and its
full unnormalized 66-dimensional terminal density are checked against an
independent integer transition calculation. A serialized two-read checkpoint
is replayed and continued to the same complete output. No successful branch
is selected and no final state is renormalized.

For four zero-mean recorded currents `J_1,...,J_4`, the connected correlation
is `kappa4=E[J1 J2 J3 J4]-C12 C34-C13 C24-C14 C23`. The comparison Gaussian
has exactly the measured two-point covariances, so it predicts zero for this
connected quantity.

| Four read ports, in source numbering | Exact connected four-point correlation |
| --- | --- |
| `(0,0,0,0)`: one port | `-200297/101616768000` |
| `(0,1,0,1)`: adjacent ports | `-3701/4064670720` |
| `(0,3,0,3)`: opposite ports | `-4537/20323353600` |

These are predictions of the specified native model, including spatial
dependence on its finite carrier. They rule out a Gaussian description of
this complete record family at the same two-point law. They do not imply
an interacting relativistic QFT or establish a quantum advantage over the
exact classical transition description of these particular records.

Summing signed currents does not reveal an unknown initial occupation.
For example, the record `(0,0,0,0)` at port zero remains compatible with
both initial occupations. Additional initial-state disclosure or a justified
inference procedure is necessary to turn current integration into a local
occupation readout.

## The observable transfer has an exact minimal three-state realization

The existing [native clock reduction](../native_clock/README.md) groups the
two occupied ports by separation: adjacent, distance two, or opposite. Strong
lumping preserves every finite absolute-current word. In this representation,

```text
P = [[14/15, 1/15,   0  ],       pi = (5/11, 5/11, 1/11).
     [ 1/15, 9/10, 1/30],
     [  0,   1/6,  5/6 ]]
```

This packet proves that three is also the **minimal linear realization
dimension** of this stationary record law: the `3x3` Hankel matrix formed
from probabilities of all-one words has determinant `1/9702990000`. A smaller
linear representation cannot reproduce all those words. This result concerns
the chosen observable family, not the full native quantum algebra or every
terminal state.

Detailed balance makes `S=sqrt(diag(pi)) P sqrt(diag(pi))^(-1)` symmetric.
Its nonconstant transfer eigenvalues and current spectral weights are

```text
r_minus, r_plus = (25 - sqrt(3))/30, (25 + sqrt(3))/30,
w_minus, w_plus = 1/7260 - sqrt(3)/16335, 1/7260 + sqrt(3)/16335.
```

Both weights are positive. For two moving-current bits with `n` unread
attempts between their retained events, the exact centered covariance is
`w_minus r_minus^n + w_plus r_plus^n`. The producer checks the spectral
recurrence symbolically, not merely a finite fit. This supplies an economical
two-pole prediction and fixes the associated Euclidean decay parameters.
The four finite record covariances at `n=0,1,2,3` recover the pole sum `5/3`
and product `311/450` without reading the hidden configurations or assuming
the eigenvalues beforehand. Their Hankel determinant is `1/9702990000`.
The resulting dimensionless ratio of logarithmic decay rates is
`log(30/(25-sqrt(3)))/log(30/(25+sqrt(3)))`, approximately `2.20334203`.
It is a **current-record relaxation ratio per native attempt**, not a
particle-mass ratio. The stationary preparation and uniform scheduler are
declared protocols, and the production E/F observer halts after two events;
neither its fixed tape nor a physical clock supplies this long-lag law.

A concrete preparation control shows why the wider readout matters. Occupied
ports `(4,7)` and `(4,8)` have the same one-bit move probability `1/6`, but
their probabilities of two successive moving bits are respectively `1/36`
and `23/900` under this protocol. The exact 66-state native transition counts
independently reproduce both values. These records resolve this finite hidden
preparation ambiguity; they do not select a particle family operator.

The current transfer is restricted to `N=sum_i n_i=2 I_66`. Therefore every
operator inside that 66-state space is neutral under its occupation `U(1)`,
and a single-port creation or annihilation operator compressed to the space
vanishes. Its two record poles cannot be charged-particle two-point poles.
The full native quantum channel does retain cross-number coherences. Direct
enumeration of all thirty seams shows that `|one token at i><vacuum|` is an
eigenoperator of the unread channel at every port, with eigenvalue `11/12`
per native attempt under coherent rule C, or `5/6` under rule E. These twelve
single-port modes are degenerate within each law, so this block supplies no
three-family mass splitting. Formally,
`Tr(b_i Phi^n(b_i^dagger |vac><vac|))` has the corresponding decay. C and E
give the same classical endpoint/current records but different coherent
decays. Measuring this expression would require a prepared vacuum/one-token
coherence, a charge reference and an analyzer outside the present current
observer. It is a dissipative Liouville coherence rate, not a state-vector
Hamiltonian gap or a selected charged-lepton mass; the hard-core ladder is
not a derived CAR electron operator. The [native action](../native_action/README.md)
keeps those mathematical objects distinct.

A number-conserving *relational* check makes the missing readout precise.
Prepare the declared total-one source/reference state
`(|one_i,0_ref>+|vacuum,1_ref>)/sqrt(2)`, retain the reference, make `n`
uniform native seam attempts on the source, then apply the declared positive
effect `(I+X_i,ref)/2`, where `X_i,ref` swaps those two states and is zero
elsewhere. Direct propagation of every current word through the original C/E
Kraus maps gives `p_+(n)=[1+r_L^n]/2`: at one attempt it is `23/24` for C and
`11/12` for E. This is an exact charged-coherence signal accessible *if* the
reference preparation and exchange analyzer are supplied. Neither is in the
current two-event observer or selected by its repair records.

The declared reference test has an encoded realization in the existing
six-level processor source words.
A five-pulse word prepares the total-one source/reference superposition;
the existing C/E repair words act on the vacuum/one-token and binary seam
columns; a 598-pulse word couples the pair-exchange effect to a pointer.
Every coin, raw record (including reserved zero-probability codes) and
pointer outcome is retained. The compiled branch densities match the
independent exact local Kraus law within the stated numerical tolerance.
For one seam incident on the source port, the exact plus probabilities are
`3/4` for C and `1/2` for E. The native degree-five, thirty-seam carrier
then gives `23/24` and `11/12` under the separately declared uniform
scheduler. The source-word execution checks the incident unit; the
25 remote seam units retain the source/reference coherence exactly. Their
compiled vacuum columns are identity with the source token as an untouched
spectator; all remote units and both coin words are checked. The
thirty-seam mixture follows these exact local units and the graph count,
while the compiled incident probabilities carry a separate floating error.
This establishes ideal encoded *realizability* using existing control
primitives. It does not make the preparation, analyzer, coherent law or
scheduler a consequence of repair, and it does not bound physical
proper-code service errors or identify a particle pole.

The records alone cannot fix an oscillation frequency. Multiplying each
native Kraus map by source charge parity after repair leaves every current
word unchanged, including adaptive words, yet changes the charge-one rate
from `r_L` to `-r_L`. The exact one-attempt plus probabilities then become
`1/24` for C and `1/12` for E. This is an **alternative quantum completion of
the classical endpoint/current rule**, not the chosen C or E instrument. The
general source-only charge rotation `U_theta=exp(i theta N_source)` gives the
same current words: each native outcome preserves source number on both sides
and commutes with `U_theta`, so a word of length `n` changes only its terminal
state by `U_theta^n`. Against the untouched reference the declared analyzer
would read `[1+r_L^n cos(n theta)]/2`. This is relative phase, not a global
gauge rotation of source and reference together. Selecting a
full quantum lift, a source/reference phase service, a physical clock and a
unitary or relativistic spectral continuation are separate requirements. A
decay pole in this declared attempt count cannot be identified with a
particle rest mass by analytic continuation alone.

The [charged spatial-pole gate](charged_reference.py) sharpens this limit using
the **same** C/E native instruments. For every port `i`, seam and retained
current, the map of `|one_i><vacuum|` is either zero or a multiple of that
same matrix unit. Nonzero output requires zero current. A seam remote from
`i` gives multiplier `1`; an incident seam gives `1/2` under C and `0` under
E. All 2,160 port/seam/current/law columns are checked. The existing
portwise signed-difference invariant of the full matrix-unit rule extends
this nontransport to any charge insertion of fixed signed difference; the
present gate computes the vacuum/one-token two-point block explicitly.
Any finite sequence of native current instruments, even with seam choice
adapted to prior current records, leaves that coherence at its original
port. Under the separately declared uniform scheduler,

```text
G_ij(n) = Tr(b_i Phi_L^n(b_j^dagger |vac><vac|))
        = delta_ij r_L^n,  r_C=11/12, r_E=5/6.
sum_(n>=0) z^n G_ij(n) = delta_ij/(1-z r_L).
```

This has a positive scalar temporal Hankel form
`sum_i |sum_n a_(i,n) r_L^n|^2`; it is not a verification of the full
Osterwalder--Schrader axioms or a physical particle reconstruction. The
charged band is flat in
**every** graph-Laplacian mode: it has no offsite response or positive-speed
dispersion. By contrast, the one-token **population** block of the same
instrument transfers across each directed native seam with probability
`1/60` per uniform attempt. Population transport cannot be substituted for a
charged two-point function. A scalar conversion of attempts to observer time
cannot create its absent spatial dispersion. This excludes interpreting this
particular C/E vacuum-to-one-token block as a propagating relativistic
particle on the recovered graph; it does not exclude a different native
charged sector or a source-selected coherent continuation. The existing
[unitary coherent walk](../observer/COHERENT_CONTINUATION.md) propagates a
neutral density amplitude after additional uncommitted-record services. It
does not supply a charged particle pole for this block.

The marked matrices inserted at the reads are transition insertions, **not
POVM effects**: their sum is `P`, and their symmetric moving insertion has a
negative determinant. A coherent probe of these transfer modes needs its own
valid preparation, instruments and service contract.

## Spatial neutral coherence in the same native clock

The [neutral transport construction](neutral_transport.py) supplies such a
complete measurement for a spatially resolved **neutral** coherence. Choose
two distinct source ports `a,b`. In the same binary total-two sector define
`E_j=|a,j><b,j|`, for the ten common occupation positions `j` outside `{a,b}`.
The signed-difference endpoints `a,b` stay fixed. The common occupation moves
on the actual induced graph `Z=G\{a,b}`, with adjacency `A_Z`, degrees `d_j`
and Laplacian `L_Z`. The existing signed-difference theorem gives

```text
Phi_L(E_j) = sum_k T_L[k,j] E_k,
T_L = I - L_Z/60 - (V_L/30) I,
V_C = 5;   V_E = 10 - 1[(a,b) is a source seam].
```

Actual C/E endpoint enumeration checks every column for all 66 choices of
`{a,b}`. Thus this coherence block has offsite response and a nonflat band
`1-V_L/30-lambda/60`, with `lambda` from the punctured source Laplacian.
For adjacent endpoints its characteristic polynomial is
`lambda (lambda-5)^3 (lambda²-10lambda+22) (lambda²-9lambda+16)
(lambda²-8lambda+10)`. These are source-derived dissipative response modes,
not charged-particle energies. The coordinate that moves is the shared
occupation; the earlier nontransport result for charged insertions remains.

The same response can be read when the **first global moving current** occurs,
using the existing [stopped native instrument](../native_relational_readout/README.md).
This protocol requires all-seam current access. It is not one port's incident
clock. The zero-current coherence multiplier is
`W_j=1-(2V_L+d_j)/60`; the moving part is `A_Z/60`. Summing every unobserved
wait therefore gives the exact first-movement response

```text
R_L[k,j] = A_Z[k,j]/(2V_L+d_j).
```

Prepare `(|a,j> +/- |b,j>)/sqrt(2)` and, after the tick, read all twenty
orthogonal projectors onto `(|a,k> +/- |b,k>)/sqrt(2)`, together with their
46-dimensional complement. The difference of plus/minus probabilities at
`k` is exactly `+/- R_L[k,j]`. Every probability is nonnegative and the full
21-outcome law sums to one. A contrast is not a probability or a herald:
its attenuation is retained in the phase and complement outcomes. For
`(a,b,j)=(0,1,2)`, each of five neighboring common positions has contrast
`1/15` for C and `1/23` for E at the first movement. All other positions have
zero contrast. No fitted time conversion enters this comparison.
For this preparation the complete measurement's C/E total variation is
`4/69`, exactly the output states' trace distance: their difference splits
into orthogonal two-dimensional phase blocks diagonalized by this analyzer.
The two candidates still have the same population and current waiting laws.

Successive global movement ticks compose as `R_L^m`. Independent inserted
identity holds cancel from the stopped law; genuine zero-current repairs
cannot be discarded. `R_L` is similar to the real symmetric matrix
`D^-1/2 A_Z D^-1/2`, with `D_jj=2V_L+d_j`. It has negative eigenvalues, witnessed
by a negative two-port principal minor, so this stopped transfer is not a
positive Euclidean transfer admitting a real self-adjoint `-log R_L`.
The change of clock cannot be used to invent a particle Hamiltonian.

The shared runtime executes both preparation signs and both laws with a
three-attempt cap, retaining all **168 flag/measurement outcomes**, including
impossible words and the unresolved timeout of probability `125/216`.
Every terminal 66-dimensional density agrees with an independent exact CP
calculation, and a serialized checkpoint gives the same continuation.
The infinite first-tick law is analytically summed, not finitely executed.
This ideal-response experiment declares coherent preparation and a joint
analyzer. Their encoded realization below uses the existing source services.
Global record access, the uniform scheduler and logical halt remain explicit;
no laboratory identification, Lorentzian limit or particle mass follows.


## Source-coded preparation and a one-bit-pointer spatial analyzer

The [neutral instrument](neutral_instrument.py) realizes the preparation and
complete spatial POVM with existing source controls and pointer services.
Two actual native repairs first prepare the required occupation: record
current `-1` on seam `(2,3)`, then either nonzero current on `(0,1)`.
Conservation of total charge two fixes the successful work state to `|0,2>`
or `|1,2>`. This purity holds for arbitrary input, including entanglement
with an untouched reference. At the uniform reference, each sign has
probability `1/264`; the failed preparation has probability `131/132`.
Every failed trial is retained. The existing 34-pulse balanced mixer turns
these two successful states into the plus and minus phase preparations.
No independent pure coherent input is supplied.

After the first global moving current, with a **two-attempt cap**, the same
34-pulse mixer converts each phase pair to its encoded occupation address.
Nine source CNOTs copy the nine work bits, one at a time, to the existing
binary pointer. Each CNOT uses three source pulses. A reset supplies zero
before each copy: one initial blank and eight subsequent resets; all nine
readout bits are retained. This executes **512 raw outcomes**, including
nonbinary and padding addresses. Twenty addresses label the ten spatial
plus/minus outcomes; every other address belongs to the complement. Full
code-basis measurement is thus decomposed into explicit source words and
the existing one-bit reset/read service.

On the native code the POVM is exactly the twenty phase projectors plus
their complement. Its terminal operation additionally dephases addresses
inside the complement: it is **not** the earlier Lüders instrument. A
complement superposition gives trace distance `1/2` between those updates.
The full finite experiment retains preparation failure, successful movement
at either attempt, and unresolved timeout. Its exact predictions are

| Complete finite output quantity | Value |
| --- | --- |
| Successful preparation followed by timeout, both signs | `25/4752` |
| C/E total variation, all preparation and timeout outcomes | `29/39600` |
| Half-L1 difference on completed-movement outcomes alone | `1/23760` |

The third row is an unnormalized subdistribution comparison; it is not a
success-conditioned probability. The ideal uncapped `4/69` distance belongs
to the separately specified pure-preparation first-movement law. The first
moving attempt has the same law under C and E; their different zero-current
coherence retention produces the second-attempt distinction.

Every successful trial uses at most **95 source pulses**, **190 proper-code
load/export events**, **nine pointer resets** and **nine pointer reads**.
Two native herald repairs precede at most two movement attempts. Discarded
uniform native warm-up removes exact stationary preparation: **369 attempts**
give a complete-output half-trace bound of at most `9.19317e-7`, including the
original binary source-word comparison and its norm amplification. The
maximum native count is therefore 373. The bound applies before conditioning
on a rare success and uniformly over untouched references. A finite native
population execution checks the warm-up attachment independently.

The source certificate encloses the original eight-dimensional mixer and
four-dimensional copy arrays, lifts them over every spectator, and composes
the two mixer calls and nine pointer copies. No array is normalized, polar
corrected or compressed to the native sector. For this comparison the native
instruments have a stated CPTP extension: dephase native code/complement,
retain complement as an absorbing failure, and use the actual native maps
on the code. This extension is a mathematical comparison contract, not a
newly derived physical code-validity detector. Runtime terminal states and
checkpoint continuation use the ideal target instrument; separate raw source
execution retains all 512 pointer probabilities. Floating propagation
residuals are diagnostics, not rigorous roundoff or physical-fidelity bounds.

The coherent proper-code service, pointer blank/reset/read operations,
carrier isolation, uniform schedule, global current access and halt control
remain explicit. Native committed repairs and classical feedback alone
preserve diagonal states and cannot create this phase preparation. The
source-word realization removes a separately supplied coherent state and
arbitrary spatial analyzer within the existing apparatus contract; it does
not derive those underlying services or identify a particle or laboratory
clock.

## Existing theorem attachments and remaining scope

The finite bond-positivity argument uses the algebra in
[`SourceLaw/ReflectionPositivity.lean`](../../Lean/SourceLaw/ReflectionPositivity.lean);
the native-history factorization is already stated in the transfer derivation.
The spectral-calculus operations align with
[`SourceLaw/TransferMatrixLift.lean`](../../Lean/SourceLaw/TransferMatrixLift.lean).
Its Gaussian covariance hypotheses are not imported into this non-Gaussian
process. The inner-unitary/modular correspondence uses the finite algebraic
structure of [`Dynamics/PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean).
None of these general theorems chooses the source scheduler, prepares the
reference state, or identifies unitary evolution with the native instrument.
The new native matrix identities and bounds are exact executable/analytic
attachments; they are not presented as new Lean proofs.

The logical resource charge is four native attempts and four retained signed
records for each complete history. The program executes native CP semantics;
it does not emit a physical coherent pulse sequence. Floating propagation
residuals are diagnostics with an explicit acceptance tolerance, not outward
roundoff or apparatus error certificates. The exact integer record law and
rational logarithm remainder are separate from those diagnostics.

The remaining physics bridge is a controlled family of larger carriers,
operational access and continuation for the relevant fields, and a scaling
limit with an identified observer clock and local spacetime structure. The
present finite model gives a concrete interacting generator and richer
readouts against which such a bridge can be checked.

## Reproduction and API

```sh
python -m model.native_effective_fields --check model/native_effective_fields/results.json
python -m model.native_effective_fields.neutral_transport --check model/native_effective_fields/neutral_transport_results.json
python -m model.native_effective_fields.neutral_instrument --check model/native_effective_fields/neutral_instrument_results.json
python -m pytest -q -W error model/native_effective_fields
```

`generator()` and `transfer_matrix()` return immutable exact matrices.
`log_polynomial(m)` and `log_remainder(m)` accept integer orders zero through
twenty-four. `exact_current_law(ports)` returns all integer terminal counts
and their common denominator; `current_moments(ports)` derives four-point
quantities. `record_spectral_reduction()` provides the exact reduced matrices,
spectral formulas and minimality witness. The manifest binds this packet;
the receipt separately binds all producer and theorem dependencies.
