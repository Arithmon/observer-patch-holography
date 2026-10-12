# A native record coupled to a scalar field instrument

This construction connects the [native observer](../DERIVATION.md) to the
existing [scalar continuum family](../../../code/source_scalar_packet/README.md).
The scalar action, vacuum, address geometry and clock are declared inputs.
The result is a composed observation and feedback experiment with a controlled
continuum comparison. It does not derive a scalar Hamiltonian from native
load repair.

## The field and its binary instrument

Use the parent's actual golden nodes, product masses, Dirichlet boundary,
dimensionless mass one and compact preparation and detector profiles. After
mass whitening, write

```text
A = I + B tensor I tensor I + I tensor B tensor I + I tensor I tensor B,
Omega = sqrt(A),       H = dGamma(Omega),       hbar = 1.
```

Here `B=W^(-1/2) K W^(-1/2)` is the parent's one-dimensional stiffness.
The profiles `f` and `g` below include square roots of the product masses.
They are the compact left preparation and right detector, respectively.
No modes are removed from this finite scalar field.

Starting from the original vacuum, the first native current `j` controls a
momentum displacement of amplitude `j gamma f`, with the parent's `gamma=4`.
At model time `t`, its coherent amplitude and detector vector are

```text
alpha_j(t) = i j gamma Omega^(-1/2) exp(-i Omega t) f / sqrt(2),
v = Omega^(-1/2) g / sqrt(2),
X = a^dagger v + v^dagger a.
```

For a committed detector record `y` in `{-1,+1}`, take

```text
K_y = [exp(-i X/2) - i y exp(i X/2)] / 2,
K_y^dagger K_y = [I + y sin(X)] / 2,
sum_y K_y^dagger K_y = I.
```

These are the centered-pointer Kraus operators used by the existing
[scalar instrument](../../../code/source_scalar_instruments/README.md).
Their boundedness makes the instrument well-defined on the oscillator
Hilbert space without an occupation cutoff. The field has only one readout
in this experiment. There is no earlier field measurement whose backaction
could be omitted.

The branch is an actual state update. With `D(beta)` the coherent
displacement, `exp(±i X/2)=D(±i v/2)`, and

```text
D(beta)|alpha> = exp[(beta.alpha* - beta*.alpha)/2] |alpha+beta>,
<alpha|beta> = exp[-(||alpha||^2+||beta||^2)/2 + alpha*.beta].
```

Consequently `K_y|alpha_j(t)>` is an unnormalized sum of two coherent
components, including their relative phase. Its squared norm is

```text
p_h(y|j) = 1/2 + y j r_h(t),
r_h(t) = (1/2) exp(-||v||^2/2) sin(mu_h(t)),
mu_h(t) = gamma <g, Omega^(-1) sin(Omega t) f>.
```

The identity uses `j` in `{-1,0,+1}` and oddness of sine. The zero-current
branch has probability one half for either record. Opposite preparation
currents give opposite response with the same covariance. A tensor spectral
evaluation can compute the coherent inner products and their overlaps
without allocating a Fock-space state vector. Floating evaluation of these
quantities is distinct from the inherited directed interval certificate.

## Local split evolution of the same apparatus

The separate local extension replaces exact continuous propagation by
canonical kick–drift–kick gates for the same whitened action `A`:

```text
p_half = p - (tau/2) A q,
q_next = q + tau p_half,
p_next = p_half - (tau/2) A q_next.
```

The position recurrence is `q_(n+1)=(2I-tau^2 A)q_n-q_(n-1)` with the
initial half kick included. Each application of `A` uses only the site and
its six axial neighbors, with the declared Dirichlet boundary. All modes
remain present. The chosen time step obeys the inherited stability bound.

Let `E_n` and `F_n` be the position blocks of the nth split propagator, so
`q_n=E_n q_0+F_n p_0`. For the detector `X=g^T q_n`, define
`a_n=E_n^T g` and `b_n=F_n^T g`. Propagating the mean and these two columns
suffices for the single terminal readout. In the original action vacuum,

```text
mu_n = g^T mean(q_n),
variance_n = [a_n^T Omega^(-1) a_n + b_n^T Omega b_n]/2,
r_split(n) = exp(-variance_n/2) sin(mu_n)/2.
```

The variance changes because the split propagator generally squeezes the
original vacuum. Gaussian characteristic functions give this probability
without allocating a full covariance matrix or Fock vector. The continuous
vacuum is an input, not a split-invariant replacement state. The final Weyl
instrument has the same effects; its general output is a superposition of
displaced squeezed states. This extension computes the terminal readout law
from its mean and variance and retains each terminal branch's coefficients
and Gaussian Gram overlaps without allocating a Fock vector.

The local recurrence is executed at q=13 and q=34 at the nearest certified
step to `39/40`; the exact attained model time is retained. The q=13
step-halving comparison fixes that attained time. An independent spectral
split propagator checks the local result. Retained traces are checkpoint
summaries and native outcome records, not a complete authenticated field
event history or native quantum-control compilation.

## One channel, retained records and native feedback

The first native event uses seam `E=(0,1)` and commits `r=(c1,j)`. The field
preparation reads that record. After the field instrument commits `y`, the
second native event uses `F=(0,5)` for `y=+1` and `E` for `y=-1`.
The actuator receives only `y`; it does not read a field amplitude or the
native source loads.

Let `L_r^E` and `L_s^a(y)` be the native Kraus operators, and let `V_j` be
the controlled field preparation displacement. On source and field, the
complete branch operator is

```text
M_(r,y,s) = (L_s^a(y) L_r^E) tensor (K_y exp(-i H t) V_j).
```

The output places that branch in the orthogonal classical record block
`|r,y,s><r,y,s|`. Summing `M^dagger M` first over `s`, then `y`, then `r`
gives identity: each native instrument is trace preserving, `V_j` and the
free field propagation are unitary, and the binary effects sum to identity.
Thus this is one CPTP instrument on arbitrary joint source/field inputs.
The stated probabilities additionally require the declared uniform binary
native preparation and independent field vacuum.

Later operations act within each committed-prefix block and append records.
Summing over every continuation recovers the previous prefix probabilities.
This preserves committed classical records; it does not require coherence
between different record values. The finite tape halts after this protocol.
There are two native repair events and one intervening field observation,
not two successive measurements of the field.

## Exact native consequences

For the uniform binary total-two shell and fair first placement coin,

```text
P(j=+1)=P(j=-1)=5/66,       P(j=0)=28/33.
```

All `66 * 2 * 2 * 2 = 528` source/coin/detector/coin histories are retained.
A history's probability is `p_h(y|j)/264`. Direct application of the native
floor/ceiling rule gives

```text
E[j2]             = 5 r_h / 132,
E[N0_final]       = 1/6 - 5 r_h / 132,
E[y j]            = 10 r_h / 33,
P(y=+1)=P(y=-1)   = 1/2.
```

Disconnecting the preparation current, or replacing the actuator input by
an independent fair bit, gives `E[j2]=0` and `E[N0_final]=1/6` with the
same two action frequencies. Reversing the preparation kick changes
`r_h` to `-r_h` and reverses both contrasts. These controls distinguish
record-dependent feedback from an action-frequency change. A shuffled bit
must be independent of the native state and record in the comparison law;
a finite empirical shuffle has its own sampling uncertainty.

## Transfer of the continuum error

Let the existing parent bound on its declared time window be
`|r_h(t)-r_infinity(t)| <= epsilon_h`. It applies to the positive current
preparation. Parity gives the same bound for negative current; zero current
has exactly zero induced sine response in both theories.

Conditioned on a native initial state and first coin, the two detector laws
have total variation distance `|j| |r_h-r_infinity|`. The remaining native
update is the same stochastic map. Retaining the entire history therefore
gives

```text
TV(full history_h, full history_infinity)
    = (5/33) |r_h-r_infinity| <= (5/33) epsilon_h,
|Delta E[j2]| = |Delta E[N0_final]| <= (5/132) epsilon_h,
|Delta E[y j]| <= (10/33) epsilon_h.
```

The same argument applies to local split propagation: opposite currents
reverse the mean while preserving covariance, and zero current has zero
sine response. The existing time-refinement certificate bounds both the
mean displacement error and the vacuum squeezing. Thus one may replace
`epsilon_h` by `min(1, epsilon_spatial+epsilon_temporal)` for a comparison
at the attained update time. At a distinct target time `s`, add
`2|s-n tau|` before capping at one, using the parent's continuum probability
Lipschitz bound. Its nearest admissible update is at distance at most `tau`.
These transferred bounds concern the ideal mathematical recurrence;
floating execution needs a separate numerical-error allowance.

The equality concerns the complete labelled histories; forgetting records
can only decrease total variation. The two sharper expectation bounds use
the exact affine formulas, rather than a generic bounded-observable bound.

At `q=233`, the retained parent certifies, uniformly for `19/20 <= t <= 1`,

```text
epsilon_h <= 12513633763 / 250000000000,
-111322580247/500000000000 <= r_h <= -34438518047/250000000000.
```

The conditional prediction therefore has negative `E[j2]` and final port-0
occupation above `1/6`; multiplying the displayed response interval by
`5/132` gives its certified current range. These are predictions of the
specified software experiment before sampling. They are not calibrated
laboratory predictions or comparisons with unseen physical observations.

The parent's convergence of finite Weyl products on the same refinement
family also applies to this finite instrument word. The explicit bound here
uses its single-readout probability certificate; no full Fock-state norm
limit, interacting continuum or physical clock is inferred. This construction
does not instantiate a field-to-native quantum compiler or provide a uniform
hardware resource bound. Preparation, field controls, the clock, record
storage and physical implementation errors remain separate supplied services.
