# A finite environment in the continuing observer

[open_environment.py](open_environment.py) attaches independently refreshed
number-phase channels to the canonical [observer program](COMPOSITION.md).
It reuses the source-executed finite environment of
[anharmonic_transfer](../anharmonic_transfer/DERIVATION_ENVIRONMENT.md).
This attachment supplies a dimensionless variance at each declared checkpoint.
It does not transfer that experiment's rate in inverse seconds, infer a noise
mechanism, or calibrate the common model's clock.

The finite environment is an unobserved complete channel. Its alternatives
are summed as density operators, without conditioning or a coherent sum of
their amplitudes. Existing visible records, disturbed work states and an
untouched reference continue through the later program. A repeated host cache
lookup requires a fresh independent draw each time; reusing a draw throughout
the history would implement a different, correlated environment.

## Schedule and working algebra

`PhaseNoise(after, site, variance, order)` specifies a cursor, scalar site
zero or one, nonnegative exact integer/Fraction variance, and quadrature order
from 1 through 64. `Environment(events)` requires an immutable tuple with
strictly increasing, distinct `(after, site)` addresses. Cursor zero means
after the single preparation and native-current gauge kick. Cursor `k>0`
means after operation `k` of the existing program. An event does not write
another observer bit or refer to future records. Each application acts on
every current visible-record block, including zero-mass blocks.

`validated(program, environment)` validates and reconstructs both the supplied
program interface and the environment, returning detached event objects.
`Environment.specification()` supplies the canonical JSON-safe schedule. Its
exact variance strings distinguish values that round to the same source float.

The work index is

\[
 a=g+4(n_0+d n_1),\qquad 0\le g<4,
\]

and an external reference is the least-significant tensor factor in the
density representation. If `r` is its dimension, full index `ar+s` has scalar
occupation `n_j`. A local Schur channel therefore acts by

\[
 \rho_{ar+s,br+t}\longmapsto S_{n_j(a),n_j(b)}\rho_{ar+s,br+t}.
\]

The adapter implements these indices directly. The source environment's
original `reference tensor oscillator` layout cannot be applied directly to
the common register, whose gauge bits lie below its scalar sites.
Independent tests instead form tensor products of every executed local Kraus
operator on both sites and on arbitrary reference inputs.

## Executed channel and complete recorded histories

Let `N_j` be local oscillator number and let `(z_m,w_m)` be the positive
standard-normal Gaussian quadrature rule of order `M`. For supplied variance
`v`, the ideal finite channel is

\[
 \mathcal Q_{M,v}(\rho)=\sum_{m=1}^M w_m
 e^{-i\sqrt v z_m N_j}\rho e^{i\sqrt v z_m N_j}.
\]

Every branch uses the existing number-phase source word: one scalar identity
phase and one Z phase for each occupation bit. The local matrix is embedded
on site wires `2+j log2(d),...,1+(j+1)log2(d)`. Global source phases remain in
the source guards. Prediction uses the Schur matrix assembled from those
executed branch phases. The analytic Gaussian multiplier appears only in the
separate comparison backend.

For an exact normalized quadrature rule,
`sum_m (sqrt(w_m) U_m)^dagger(sqrt(w_m) U_m)=I`. This proves complete positivity
and trace preservation on arbitrary work/reference operators. In a visible
record block `h`, a later instrument has operators `K_(h,b)`. A noise insertion
and that instrument compose as

\[
 \rho_h\longmapsto\bigoplus_b\sum_m w_m
 K_{h,b}U_m\rho_hU_m^\dagger K_{h,b}^\dagger.
\]

Summing `b,m` preserves the old block's trace. Repetition proves completeness
of the whole visible-record history, preservation of earlier record
probabilities after summing their continuations, and preservation of the
untouched reference marginal. It does not preserve the old conditional work
state. Treating a discarded draw as a visible bit would define a different
output algebra; averaging `sqrt(w_m) U_m` would generally fail completeness.

The [encoded-memory identity](ENCODED_MEMORY.md) already holds with an
arbitrary complete work instrument. Replacing its work instrument by this
noisy continuation therefore retains the logical identity tensor the work
channel. After summing every native/work record and tracing work, the logical
and untouched-reference marginal is preserved, including initially correlated
logical/work inputs. This statement neither protects the work state nor
corrects physical errors on the native storage register.

## Reference-stable finite error budgets

Write `R=d-1`. The source API receives `v_f=float(v)` as its duration and
rate one; these arguments denote a dimensionless phase variance here. The
ideal Gaussian comparison at exact declared variance is

\[
 (\mathcal D_v\rho)_{nm}=e^{-v(n-m)^2/2}\rho_{nm}.
\]

The reused finite-rule theorem gives

\[
 \epsilon_{\rm quad}(d,M,v_f)
 =\min\left(1,\frac{d M!}{2(2M)!}[v_f R^2]^M\right)
\]

as a uniform half-diamond bound between the ideal finite rule and
`D_(v_f)`. The factor is local occupation dimension `d`, not the dimension
of the work register or its reference. Exact Fraction arithmetic uses the
stored binary value of `v_f`.

The change of declared variance has a separate bound. Center `N_j` by its
spectral midpoint, so `B=N_j-RI/2` has norm `R/2`. Scalar shifts disappear from
commutators, and `||ad_B||_diamond<=2||B||=R`. The Gaussian variance generator
is `-ad_B^2/2`; integrating its CPTP semigroup gives

\[
 \tfrac12\|\mathcal D_v-\mathcal D_{v_f}\|_\diamond
 \le\epsilon_{\rm conversion}
 :=\min(1,|v-v_f|R^2/4).
\]

Consequently the ideal represented finite channel differs from the
exact-variance Gaussian channel by at most
`min(1,epsilon_quad+epsilon_conversion)`. This charge remains nonzero in exact
arithmetic when a positive supplied variance underflows to source float zero.

There is also a direct finite-channel bound against identity that does not
pay the quadrature error. For a phase `theta`, pair its two signs. The first
derivative of `(exp(theta L)+exp(-theta L))/2` vanishes at zero, where
`L=-i ad_B`. Its second derivative has diamond norm at most `R^2`, because
unitary conjugations have norm one. Taylor's integral remainder therefore
gives half-diamond distance at most `theta^2 R^2/4`. The exact quadrature rule
is symmetric and, for `M>=2`, has second moment one. Thus

\[
 \boxed{\tfrac12\|\mathcal Q_{M,v_f}-I\|_\diamond
       \le\min(1,v_f R^2/4).}
\]

For `M=1`, the only node is zero: the finite channel is exactly identity,
even if the requested Gaussian variance is nonzero. Its Gaussian-comparison
error is not zero. The field
`ideal_exact_variance_finite_to_identity_half_diamond_upper` describes the
finite rule at `v`; `ideal_finite_to_identity_half_diamond_upper` describes
the represented finite rule at `v_f`.

`schedule_budget` sums either local half-diamond radius over every scheduled
invocation and caps the result at one. Arbitrary intervening complete
instruments, native continuations and reference systems cannot amplify a
local radius; the telescope has no multiplier for the number of visible
branches. Different insertions may use different sites, variances and orders.
These are discrete checkpoint channels. No comparison to simultaneous noisy
Hamiltonian/Lindblad evolution is asserted. In particular the specialized
anharmonic pulse Strang bound is not a bound for the common Hamiltonian.

The previous [prepared smooth-history bounds](SMOOTH_HISTORY.md) target its
specific noise-free trajectory. They do not become noisy-history cutoff
certificates by renaming the target. One valid triangle comparison is the
new finite-to-identity schedule radius plus the previous prepared ideal
finite-to-infinite bound, with that **same noise-free infinite target** and
the previous preparation hypotheses explicitly retained. Source arithmetic
and physical service errors are separate obligations in that comparison.

## Occupation refinement through noise

The uniform continuation proof carries the summed fourth moment
`w=(N_0+N_1+2)^4` and the rectangular projection `Pi_d`. Each local phase
unitary commutes with both. Therefore

\[
 \sum_m w_m U_m^\dagger w U_m=w,\qquad
 U_{m,d}\Pi_d=\Pi_d U_m.
\]

For the common Stinespring dilation
`V psi=sum_m sqrt(w_m)|m> tensor U_m psi`, projection commutes with `V`, and
`V` is isometric. A phase insertion thus leaves both the moment amplitude
and the projected comparison error unchanged. Its row in the
[continuation recurrence](COMPOSITION.md#repeated-smooth-reads-and-projected-error-induction)
is `(e,A)->(e,A)`. Subsequent Hamiltonian and smooth-read rows retain their
existing estimates. The incoming preparation, coupling and unconditioned
history hypotheses remain necessary; arbitrary supplied inputs do not acquire
a vacuum occupation certificate.

This proves the same occupation bound for a finite mixture with fixed
schedule and quadrature order, comparing its finite and infinite occupation
versions using the same phases. A Gaussian environment admits the same
argument with its continuous phase dilation. It does not follow that fixed
`M` approximates the Gaussian uniformly as `d` grows: the displayed quadrature
bound grows with the available occupation gap.

An analytic joint limit can choose increasing order. Since
`M!/(2M)! <= (M+1)^(-M)`, choosing
`M(d)>=max(1,2 v_max)d^2` for the finitely many exact ideal variances gives
each quadrature radius at most `d*2^(-M(d))/2`, which tends to zero. Combining
this with the unchanged occupation recurrence and the existing temporal
refinement proves a joint **ideal** limit. The implemented order cap 64 and
source cutoffs 2 and 4 do not execute that unbounded family. Holding a
non-dyadic variance at one binary64 representation also does not eliminate
its separately reported conversion error or certify numerical convergence.

`open_bounds.gaussian_refinement` also admits the sufficient order
`M(d)>=max(d,ceil(2v(d-1)^2))`. It reports the geometric radius rounded
outward on the grid `1/[10^12(d+2)^6]`. The rounding allowance tends to zero;
the helper need not construct or serialize an enormous power-of-two
denominator to retain a valid convergence certificate.

## Resources, cache integrity and evidence

Every invocation costs `log2(d)+1` native pulses on a realized environment
branch, independent of the drawn node, including zero variance. The expected
and maximum branch counts are therefore equal. Contracting all `M` branches
costs `M(log2(d)+1)` pulses in the enumeration ledger. The adapter retains the
corresponding source load/export/service counts and one independent classical
draw per invocation. The schedule sum is per realized observer path; the
runtime's expanded visible-record tree has a separate application count.
Neither host caching nor discarding the draw makes its supply free.

Before each source application, the adapter copies the cached arrays and
metadata. It checks weights against the declared canonical quadrature,
requested parameters, every phase branch and its emitted-word hash, resource
counts, and agreement of the Schur matrix with the complete copied Kraus sum.
It uses those executed arrays after validation, not the target phases used
in the guard. Caller mutation of returned metadata cannot mutate the source
cache. The tolerances are rejection guards, not directed roundoff bounds.

The executable adapter has an explicit certificate-serialization budget.
Every numerator and denominator of the supplied exact variance, represented
variance, conversion charge, inherited quadrature certificates and combined
event/schedule radii must have
a conservatively estimated decimal length of at most 4,000 digits, or the
interpreter's active limit if that is smaller. The check includes the
uncapped quadrature remainder retained by the reused source producer, even
when its half-diamond radius clips to one. It runs before source execution
and never changes the process-wide integer conversion policy. Unsupported
order/variance combinations raise a specific `ValueError`; they do not
receive a zero error or a truncated exact certificate. For example, the
smallest positive binary64 variance is supported at order two but exceeds
this reporting budget at order sixteen. Its positive half-sized neighbor,
which rounds to zero, remains supported with a nonzero exact conversion
charge. These are computational-domain limits, not singularities of the
mathematical phase channel.

The inherited source channel rounds weights, stabilizes executed unitaries
and makes a small diagonal representation correction. Those changes and
floating contractions are not included in the exact ideal bounds above.
Unknown source-service, physical storage, detector, variance and independent
draw accuracy remain unknown. A supplied noise model is not empirical
evidence of that model or of OPH's physical completion.

[test_open_environment.py](test_open_environment.py) independently sums all
executed branches on both scalar sites and complex reference inputs, checks
every matrix unit at the smaller work dimension, preserves reference and
occupation marginals, separates fresh from reused or coherently added draws,
checks Gaussian and identity radii including order one and variance underflow,
and rejects malformed schedules and poisoned source caches.
