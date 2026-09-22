# What the full finite A3 objective can force

This result concerns the canonical finite-algebra A3 objective, including
noncommuting local states. It does not replace that objective by entropy on
a supplied schedule simplex. The source of an M1 failure observable and the
complete A1/A2 feasible family must be identified before applying it to M1.

## Theorem: local support and exact public requirements

Let K be a nonempty convex set of compatible normalized local states on a
finite net of finite-dimensional C*-algebras. A global extending state is
not assumed. Let G be the finite A3 observer cover, w_P > 0 its exact
weights, and tau_P faithful local reference states. Use the declared
faithful traces, natural logarithms, and the support-aware Umegaki
divergence. Suppose rho* attains the minimum of

    F(rho) = sum_(P in G) w_P D(rho_P || tau_P).

Then, for every sigma in K and every P in G,

    supp(sigma_P) <= supp(rho*_P).                         (1)

Consequently every positive effect E in a scored algebra satisfies

    rho*_P(E) = 0  iff  sigma_P(E) = 0 for every sigma in K. (2)

The same equivalence holds for a finite nonnegative sum of such
expectations, and for a record obtained through an actual positive
Heisenberg pullback into a scored algebra. Cover injectivity alone does not
supply a positive pullback. Uniqueness and cover injectivity are not needed
for (1); canonical A3 supplies them for its other purposes.

In particular, if E_bad is the positive failure effect of a required read,
the selected state makes the read certain exactly when the entire feasible
family already excludes its failure. The conclusion is about that
identified effect, not every nonlinear predicate of a state, an unscored
reconstructed history, or a property of the underlying model space.

### Proof without simultaneous diagonalization

Fix rho = rho*, sigma in K and 0 < t < 1. The mixture
rho_t = (1-t)rho + t sigma belongs to K, with the same mixing coefficient
at every patch. Local reference faithfulness makes all objectives finite.
Expansion of the trace formula gives the exact identity

    (1-t)D(rho_P||tau_P) + t D(sigma_P||tau_P)
      - D((rho_t)_P||tau_P)
    = (1-t)D(rho_P||(rho_t)_P) + t D(sigma_P||(rho_t)_P). (3)

All logarithms on the right are restricted to the support of rho_t; it
contains both endpoint supports. In particular no logarithm of a zero
eigenvalue is silently totalized. This identity uses linearity of the
trace, not commutation of rho, sigma or tau.

If (1) fails at P, let Q be the kernel projection of rho_P and put
b = Tr(Q sigma_P) > 0. Measure the two effects Q and 1-Q. Relative entropy
data processing bounds the right side of (3) below by the binary mixture
gap for distributions (0,1), (b,1-b), and (tb,1-tb). In that gap the
rho term is nonnegative, and the sigma term obeys

    D((b,1-b)||(tb,1-tb))
      >= -b log t + (1-b)log(1-b)
      >= -b log t - H(b),

where H(b) = -b log b -(1-b)log(1-b), with 0 log 0 = 0. The endpoint b=1
gives exactly -log t and is included. At every other patch the right side
of (3) is nonnegative. Multiplying by positive weights and summing gives

    F(rho_t) - F(rho)
      <= t [F(sigma)-F(rho)+w_P H(b)] + t w_P b log t.    (4)

The last coefficient is strictly positive and log t tends to minus
infinity, contradicting minimality. This proves (1). For positive E,
Tr(rho E)=0 means E annihilates the support of rho: equivalently the
positive matrix rho^(1/2) E rho^(1/2) has zero trace and is zero. Support
inclusion gives the forward implication of (2); rho* in K gives the
reverse. Positive pullback and nonnegative sums preserve this argument.

The matrix facts used here are finite-dimensional relative entropy
nonnegativity, the trace identity (3), and data processing under a binary
measurement. See John Watrous, *The Theory of Quantum Information*,
Section 5.2, particularly Theorem 5.35, in the
[author's text](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).
The weighted compatible-family argument above requires no global density
matrix and no assertion that separately optimized local states agree.

## Quantitative obstruction: a feasible witness gives a failure floor

The support result first makes all D(sigma_P||rho*_P) finite. The one-sided
derivative at t=0 can then be taken on each rho*_P support. The positive
eigenvalues there are bounded away from zero at this fixed finite cutoff.
Minimality along rho_t and the derivative of Tr(X log X) give

    sum_P w_P D(sigma_P||rho*_P) <= F(sigma)-F(rho*).     (5)

This is the information-projection Pythagorean inequality for this convex
family. It is not an equality assumption or a stationarity oracle.
Normalization cancels the derivative's identity terms. The difference
between the right and left sides is precisely
sum_P w_P Tr((sigma_P-rho*_P)(log rho*_P-log tau_P)),
the nonnegative directional derivative.

For any effect 0 <= E <= 1 in patch P, write b=sigma_P(E)>0 and
p=rho*_P(E). Measuring E and 1-E, then dropping the nonnegative divergence
terms from other patches, yields w_P d(b||p) <= Delta, where
Delta=F(sigma)-F(rho*). The binary bound d(b||p) >= -H(b)-b log p gives

    p >= exp(-(Delta/w_P + H(b))/b)
      >= exp(-(F(sigma)/w_P + H(b))/b) > 0.              (6)

For b=1 this becomes p >= exp(-F(sigma)/w_P). Thus a certified feasible
failure witness gives a positive floor without solving the optimizer. The
bound can be very small; it is not a regulator-uniform lower bound unless
the witness costs and cover weights have corresponding uniform bounds.
If p=1 the conclusion is immediate; p=0 is excluded by (1).

There is also a certificate against falsely declaring a zero-support state
approximately optimal. If rho_P(E)=0 and sigma_P(E)=b>0, use the effect
measurement directly in (3). Let a=w_P b and C=F(sigma)+w_P H(b)>=0.
Equation (4), with F(rho)>=0, and t=exp(-1-C/a) prove

    F(rho_t) <= F(rho) - a exp(-1-C/a).                 (7)

Such rho cannot be epsilon-optimal when epsilon is smaller than this
explicit improvement. A floating-point zero and an approximate solver
status therefore do not certify exact absence of failures.

## Singular references are a different constraint

If some tau_P is singular, replace K by its finite-objective face

    K_fin = {rho in K: supp(rho_P) <= supp(tau_P), P in G}.

Assume K_fin is nonempty. Compress each scored algebra to supp(tau_P).
The same proof applies, and (1)--(7) quantify over K_fin. They do not
quantify over the discarded infinite-divergence states in K. For example,
on the full qubit state space, tau=|0><0| selects |0><0| and excludes the
otherwise feasible |1><1|. The exclusion came from the reference support.
Calling that reference a source derivation of a read would be circular
unless its zero support had itself been derived from A1/A2.

## Exact noncommuting controls and a temporal distinction

For each rational unit vector v=(a,b), with 0<a,b<1, take the convex segment
between R=|0><0| and S=|v><v|, and tau=I/2. R and S do not commute. The
orthogonal unitary U=[[a,b],[b,-a]] exchanges them. Unitary invariance and
strict convexity of relative entropy give the unique midpoint optimizer
M=(R+S)/2. Its determinant is b^2/4>0, and the failure effect |1><1| has
probability b^2/2, although it is zero at R. Eigenvalues of M are (1+a)/2
and (1-a)/2. The packet uses (a,b)=(3/5,4/5),(5/13,12/13),(8/17,15/17).
No numerical optimizer or tolerance decides these claims.

A three-level face, diag(1-t,t,0), has midpoint diag(1/2,1/2,0) for the
faithful reference I/3. The last-coordinate failure is excluded throughout
the feasible family. This distinguishes maximal feasible support from an
incorrect claim that every optimizer is faithful on its whole algebra.

Finally, identity and complete-depolarization channels both preserve I/2
and are covariant under every unitary, but their classical two-time records
differ. Starting with a uniform recorded bit, identity has mismatch
probability zero; depolarization has mismatch probability 1/2. Their
normalized Choi states are respectively |Omega><Omega| and I/4. On the full
common convex set of normalized Choi states with input marginal I/2,
entropy relative to I/4 selects depolarization. Its mismatch effect is the
sum of the |01> and |10> projectors. The bit-flip channel is a feasible
deterministic-failure witness, with divergence log 4. Equation (6) bounds
selected mismatch below by 1/4; its actual value is 1/2.

Imposing the zero-mismatch face instead selects diag(1/2,0,0,1/2), the
dephasing channel, uniquely. Indeed positivity restricts the Choi support
to span{|00>,|11>}; the uniform state on that support obeys the marginal
constraint and has the largest entropy log 2. This enforces the identity
law on classical records while permitting phase destruction. It does not
select the identity quantum channel. The required record constraint is
the additional input that changes this optimizer.

This is an exact state-versus-process distinction, not a full A1--A3
countermodel: no claim is made that both channels belong to the actual A1
update grammar. It shows why a stationary selected state, even together
with full unitary covariance, cannot be substituted for a derived temporal
instrument or a preserved version record. The Choi objective is explicitly
a second, declared process-state problem, not the original ontic problem.

## Consequence for the unchanged M1 exit

For an actual M1 candidate, identify the source-generated finite record
algebra, the public meaning of a failed required read, its positive effect
in a scored algebra (or its positive restriction pullback), and the complete
A1/A2 feasible family. With faithful references, (2) is a necessary and
sufficient test for exact success under A3; (6) certifies failure if an
admitted state violates it. Apply it separately to required reads, forbidden
reads, missing population records, and overwritten version records only
when each is represented by such an effect in that same architecture.

This settles the state-selection step for those observables. It does not
show that the canonical architecture forces or admits any particular one
of these failure effects. There is no full-axiom independence theorem here.
The missing source grammar cannot be replaced by a successful-history
conditioning, a zero in the reference, or an unproved identification of an
ontic state with a channel. No issue is closed by this theorem, and literal
M1 is not claimed. The obligations remain in this PR's original exit rather
than being transferred to a newly proposed issue.

## Verification boundary

The matrix proof, its local-family application, and the entropy derivative
are analytic. `SourceStateSelection.lean` checks the real inequality
reductions, weighted aggregation, support-exclusion logic, explicit
improvement and probability-floor arithmetic. Its entropy-gap hypotheses
are discharged by the analytic matrix argument above, not by a Lean
formalization of Umegaki data processing. The exact packet checks matrix
identities, positivity/rank witnesses, unitary symmetry, and actual channel
actions independently of its producer. It is regression evidence for the
examples, not a finite enumeration proof of the general theorem.
