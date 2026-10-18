# A spatial fluctuation field from one native history

Long occupation averages of the native repair process have a Gaussian field
limit with a calculable spatial action. The same derivation works in a port's
incident clock, with an explicit correction for its random sampling. On a
**supplied periodic square graph**, both clocks give the same massless
Gaussian field on every fixed finite Fourier band after spatial refinement.
These are two successive limits: observation count first, spatial mesh second.
The spherical source screens have not met the required Green convergence
criterion. Neither result identifies native evolution with a relativistic wave.

This extends the existing [passage geometry](README.md),
[native replica fluctuations](../native_fluctuations/README.md), and
[one-history circulation](../native_circulation/README.md). The new observable
is the time sum of occupations on **one continuing history**. It does not use
independent replicas, impose a Gaussian microscopic state, or replace repair
by a chosen Hamiltonian. The
[incident-field correspondence](../observer_time_fields/README.md) supplies the
stopping construction; the calculation below extracts its occupation
fluctuation action and its spatial convergence conditions.

```sh
python -m model.physical_space_time.occupation_field --check model/physical_space_time/occupation_field_results.json
python -m pytest -q -W error model/physical_space_time
```

## Source, observer access and compatibility

On the original twelve-port carrier, prepare the uniform diagonal law on the
66 binary configurations of total two. Use the existing native endpoint rule,
uniform independent selection of its thirty seams, and fair placement coins.
Both C and E have this diagonal law. Their coherent outputs are not equated.
Stationarity, the preparation, scheduler and readout access are premises.

The vector observer is given the initial occupations and every seam/current
record. It computes each subsequent occupation from the native continuity
equation, then samples the vector either after each attempt or at attempts
incident to port zero. The second clock ignores remote attempts in its count;
the vector decoder receives their currents. A port-only observer with
its initial occupation and incident-current tape can read the scalar at that
port, but cannot reconstruct the complete vector. The execution checks a
counterexample with the same initial state and identical port-zero incident
tape, but different remote occupations at the first tick. Full-field access
is a supplied observation interface, not a hidden-state oracle relabeled as
a local measurement. No intervention or feedback occurs during these averages.

The independent compatibility requirements are conservation of the disclosed
total, continuity-based readback, normalized complete stopped laws, and
agreement between the native configuration process and field predictions.
All 3,960 native branches are checked. An independent 66-state transition
calculation checks the field closure, stationary distribution and lag laws.
All counter moments include every terminal configuration; no event is
postselected. The large torus example reuses the proved local binary identity
and checks 16,384 actual single-particle native branches on its 64-port graph.

## Finite-graph theorem and the source-derived action

The algebra applies to any finite connected simple graph with `V` vertices,
`M` seams and a uniform fixed total `0<q<V`. Its binary native step is the
average of identity and endpoint exchange. Let

```text
eta = occupations - (q/V) 1,
P0 = I - 11^T/V,
chi = q(V-q)/[V(V-1)],
L = graph Laplacian,                  E[eta eta^T] = chi P0.
```

Choose a nonempty fixed set of `m` clock seams. An attempt in that set ticks,
including attempts with zero current. Let `L_i,L_r` be the incident and remote
Laplacians. On density vectors the remote and incident subtransfers are

```text
A = (1-m/M) I - L_r/(2M),
B = (m/M) I - L_i/(2M),
D = I-A,               K = B D^(-1).
```

Every remote word is summed before the next incident attempt. Its tail mass
after `h` attempts is `(1-m/M)^h`, independently of the configuration.
`K` is the exact completed density map, not a mean waiting-time substitution.
The full configuration stopped kernel is bistochastic and keeps the uniform
fixed-total preparation stationary. It is irreducible: each incident exchange
is a one-tick transition; a remote exchange followed by an incident no-change
branch is also a possible one-tick transition. Connected exchanges generate
the fixed-number sector. Positive self-loops give aperiodicity.

On the constant-free space, write `L^+` for the inverse of `L`, zero on the
constant mode. Since `I-K=(L/(2M))D^(-1)`,

```text
H = (I-K)^(-1) = (2m I+L_r)L^+,
Sigma = chi [H+H^T-P0]
      = chi [4m L^+ + L_r L^+ + L^+ L_r - P0].                 (1)
```

The inverses in these formulas are restricted inverses, not inverses on the
conserved constant direction. For the vector sum
`Y_n(u)=n^(-1/2) sum_(j<floor(nu)) eta_j`, the finite-chain functional CLT gives

```text
Y_n -> W,        Cov(W(u),W(v)) = min(u,v) Sigma.              (2)
```

Here `u` is rescaled count of the specified observer events. To check the
imported theorem's premises, use `g(eta)=H eta`. Its Poisson equation is
`g-E[g_next|eta]=eta`; the sum differs from the martingale with increments
`H(eta_next-K eta)` by a bounded endpoint term. The finite state space bounds
these increments. The ergodic theorem makes their averaged conditional
covariance converge. Its stationary value is exactly
`H chi(P0-KK^T) H^T=Sigma`, independently checked in the executable. This
satisfies the bounded-jump case of
[Whitt, Theorem 2.1](https://www.columbia.edu/~ww2040/WWproofs2007pub.pdf).
The probabilistic limit theorem is reused; finite tests do not prove a CLT.

At `u=1`, the field density is proportional to `exp(-y^T Sigma^+ y/2)` on
the constant-free subspace. The receipt contains its exact covariance,
precision and positive rational LDL pivots. On the twelve-port source,
`chi=5/33`. The two clock formulas become

```text
Sigma_attempt  = (5/33)(120 L^+ - P0),
Sigma_incident = (5/33)(20 L^+ + L_r L^+ + L^+ L_r - P0).
```

This is an occupation-history fluctuation action, not the action of the
instantaneous occupation state. At a finite averaging count the state remains
discrete and non-Gaussian. The exact first four moments, computed independently
from complete 66-state histories, retain nonzero fourth cumulants. No finite
discrete law is asserted close to a continuous Gaussian in total variation.

The mean wait is six attempts per port-zero incident tick. Nevertheless,
`Sigma_incident != Sigma_attempt/6`. The difference after division by clock
exposure, `Sigma_incident/5-Sigma_attempt/30`, has entries
`(0,0)=-5/216`, `(1,1)=23/1080`, `(0,1)=61/11880`. Its anisotropy rules out
repairing this discrepancy by a scalar amplitude or clock conversion.
The action precision at this finite source is not a scalar multiple of `L`.

## Finite counts, caps and observation errors

Stationarity and density closure give an exact finite-count covariance:

```text
Sigma_n = chi [P0 + sum_(h=1..n-1) (1-h/n)(K^h+(K^T)^h)P0].   (3)
```

If `||K^h P0|| <= c r^h`, with `r<1`, summing the geometric series gives

```text
||Sigma_n-Sigma|| <= 2 chi c r/[n(1-r)^2].                    (4)
```

The exact source bounds `2P0<L<10P0` are checked by rational LDL certificates.
For the attempt clock use `c=1,r=29/30`. For the incident clock, `D` lies
between `I/6` and `I/3`. Its similarity to the symmetric positive matrix
`D^(-1/2) B D^(-1/2)` gives `c=2,r=9/10`. These are conservative operator
bounds. Equation (4) is a covariance error, **not a distributional CLT rate**.

A finite implementation that caps each requested interval at `h` attempts
must keep timeout. For `n` requested ticks its total timeout probability is
at most `n(1-m/M)^h`. Coupling to the completed experiment bounds the
complete-output total variation by that probability; expectations of bounded
characteristic functions differ by at most twice it. Fixed-cap success laws
must not be renormalized and treated as completed laws. Increasing caps with
`n(1-m/M)^h -> 0` preserves (2). No auxiliary Poisson clock is required.

## A quantitative spatial-limit criterion

For a graph family indexed by `a`, supply a positive spatial scale `s_a`,
fixed-total preparation `chi_a>0`, and a declared observation matrix `W_a`
whose columns lie in the constant-free space. These columns specify spatial
smearings, including their normalization; neither a plot nor a graph rank
chooses them. Normalize the history average by `1/sqrt(4 chi_a m_a s_a)`.
Its limiting covariance on those observations is exactly

```text
C_a = W_a^T (s_a L_a)^+ W_a + E_a,
E_a = W_a^T (L_r L_a^+ + L_a^+ L_r - P0) W_a/(4 m_a s_a).    (5)
```

If the first term differs from a target covariance `Gamma` by at most
`epsilon_G`, and `||E_a||<=epsilon_clock`, then
`||C_a-Gamma||<=epsilon_G+epsilon_clock`. For a Fourier test vector `theta`,
the corresponding Gaussian characteristic functions differ by at most
`||theta||^2 (epsilon_G+epsilon_clock)/2`. Equation (4), multiplied by
`||W_a||^2/(4 chi_a m_a s_a)`, also bounds the finite-count **covariance**
error. These are observable matrix quantities that can be measured or
certified; they do not assume that agreement of means establishes a field law.

Thus Green convergence and vanishing clock correction give convergence of
the limiting Gaussian finite-dimensional laws. The order is explicitly
`lim_(graph) lim_(observation count)`. No uniform simultaneous limit or
uniform finite-count distributional error has been proved. A fixed total
two on growing graphs has `chi_a -> 0`; its amplitude normalization cannot
be carried over from the twelve-port value unchanged.

## Constructive control: unchanged repair on a supplied two-dimensional torus

Choose the `N by N` periodic square graph, a chart of period `2pi` in each
direction, `s_N=N^2/(4pi^2)`, and canonical total `floor(N^2/3)`. This is a
declared graph/preparation family, **not the existing spherical federation**.
The unchanged native binary rule gives every transition. No oscillator,
Hamiltonian or different update law is inserted. The field observer reads
normalized spatial Fourier sums from its disclosed occupation/current data.

For any fixed finite set of distinct nonzero integer wavevectors, retain
normalized real sine and cosine modes with no duplicate signs or lattice
aliasing. Their Laplacian energies are

```text
kappa_N(k) = (N/pi)^2 [sin^2(pi k_1/N)+sin^2(pi k_2/N)],
0 <= |k|^2-kappa_N(k) <= pi^2 (k_1^4+k_2^4)/(3N^2).          (6)
```

The canonical susceptibility is
`q_N(N^2-q_N)/[N^2(N^2-1)]`, and is divided out in (5).
For the attempt clock (`m=M=2N^2`), the covariance on this band is diagonal:
`1/kappa_N(k)-pi^2/(2N^4)` for each real mode.

For the clock counting the four seams at the origin, let
`G_N=diag(1/kappa_N(k))` and `J_N=W_N^T L_i W_N`. Its covariance is

```text
C_N,incident = G_N + I/(16s_N) - (J_N G_N+G_N J_N)/16.       (7)
```

This retains finite-mesh cross-mode correlations. For the sum `S` of
`|k|^2` over real modes, their normalized endpoint differences give
`||J_N|| <= 16pi^2 S/N^4`. Equations (6) and (7) bound the clock correction by
`pi^2/(4N^2) + 2pi^2 S/(kappa_min N^4)`. Hence **both actual clocks converge
to covariance `diag(1/|k|^2)` on every fixed band**, with an `O(N^-2)` bound.

The execution retains wavevectors `(1,0),(0,1),(1,1),(2,0)` and their eight
real modes. It uses `pi^2<10` to store rational error bounds:

| N | Incident covariance error bound | Computed error diagnostic |
| --- | ---: | ---: |
| 8 | 0.187270 | 0.096215 |
| 16 | 0.028450 | 0.022855 |
| 32 | 0.006046 | 0.005647 |
| 64 | 0.001446 | 0.001408 |

The inequalities certify the exact mathematical matrices. Floating band
evaluations check them independently; their residuals are diagnostics.
The theorem applies to arbitrary fixed finite bands. Their limiting Gaussian
laws are compatible under deleting modes, and define the zero-mode-removed
massless Gaussian generalized field on the declared torus. This identifies
its finite-dimensional laws, not uniform convergence in a spatial path-space
topology. The target is a random distribution: for every `epsilon>0`, its
expected squared `H^(-epsilon)` norm is proportional to the convergent sum
`sum_(k!=0) (1+|k|^2)^(-epsilon)/|k|^2`. This target construction does not
prove tightness of the native field sequence. The formal action is the Dirichlet energy
`(1/2) integral |grad phi|^2` in the declared normalized spatial measure.
For the relation between inverse-Laplacian covariance, generalized fields and
two-dimensional conformal Dirichlet energy, see
[Sheffield, Sections 2.1–2.4](https://arxiv.org/pdf/math/0312099).
This is a Euclidean free-field connection under the supplied graph/chart.
No interacting CFT, stress tensor or Lorentzian continuation is certified here.

The same native mean density has the dissipative multiplier
`1-lambda/(2M)` per attempt. Its long-wavelength energy ratio at momenta
`2k` and `k` tends to four. A positive deterministic rescaling of this clock
does not turn decay into a propagating wave. The Gaussian spatial action
therefore does not identify the original repair dynamics with relativistic
field evolution. Other observables or selected coherent continuations need
their own dynamical correspondence.

## Requirements for the actual source screens

The [coherent screen result](../coherent_screen_observer/README.md) compares
complete local quantum outputs on two matched finite neighborhoods. The
[compiled clock history](../compiled_clock_history/README.md) implements the
finite local preparation/control/readout contract. Both are compatible with
the conditional interpretation here; neither supplies the spatial maps and
Green convergence in (5). A fixed-radius local history bound does not control
the inverse-Laplacian contribution from arbitrarily slow global modes.

For those actual screens the next certificate must supply compatible spatial
smearings and scales, bound their Green matrices, evaluate the incident-clock
correction, and carry preparation/readout errors through the same comparison.
This packet gives that test and a constructive source-law control, without
asserting that the screens pass it. Physical dimensionality, Lorentzian
propagation, interacting continuum fields, and identification with observed
physics remain separate obligations. No federation axiom is needed for the
finite single-carrier result.
