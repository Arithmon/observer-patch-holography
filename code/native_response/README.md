# One N4 process: port geometry and a non-native timing response

This calculation resolves the finite bare-observable interface of the existing
[N4 candidate](../../docs/MODEL_CANDIDATES.md#n4-a-reciprocal-coherent-exchange).
It adds no interaction, state factor or axiom. N4 retains the native port
Dirichlet geometry at every real coupling, but its linear occupation space is
not dynamically closed at the two certified couplings. This supplies an
internal discriminator against a native heat law with a refitted clock.

## Fixed inputs and experiment

The state space is the binary total-two shell (dimension 66), with one
six-dimensional private factor per port. The response matrices are exactly
the formula in [SOURCE_RECORD_GLUING.tex](../../paper/tex_fragments/SOURCE_RECORD_GLUING.tex),
in the chart of [verify_response.py](../source_selection_model/verify_response.py).
Here the first neighboring pair is **(0,2)**; (0,1) belongs to a different
distance class. The seam set is derived from this chart, not imported from
another labeling. The existing `response.json` independently binds all twelve
matrices and sixty proper carrier actions.

The candidate inputs are its binary quantum representation, response coupling
θ, equal quadrature weights, exponential exchange, the full normalized-trace
state τ, and independent uniform choices among thirty seams. The measurement
input is access to the unchanged commuting occupations N_p. These are not
consequences of A1–A3 or physical calibration. Neither mass, inertia, a bulk
boundary condition nor a new field is introduced by this calculation.

The two rational couplings reuse the preceding feedback test; neither is
fitted to this result. The discriminator is strict positivity of the
full centered-subspace memory defect after error enclosure, against exact
zero for the native and one-term controls. The compute budget uses local
72-dimensional matrices and three representative edges with protected
18-dimensional blocks, rather than constructing the full network matrix.

In independent trials, measure all initial binary loads, retain the outcome,
apply n unmodified N4 events, then measure all final loads. Use n=1 and n=2;
there is **no intermediate measurement, private reset or postselection**.
The first unconditioned measurement preserves τ. The Lüders joint law gives
`E[f(X_0)f(X_n)] = τ[f Φⁿ(f)]`, where Φ discards the event records. Thus the
following moments are operational moments of one declared instrument, not
three separately chosen geometry, response and timing fixtures. Branch label
j retains its transformed-occupation meaning; it is not relabeled a bare
current. Event count, or a supplied Poisson attempt rate, is model time.

## Geometry theorem, including all couplings

Use normalized Hilbert–Schmidt inner products throughout. On the odd
two-endpoint load block, `Z=diag(1,−1)`, private dimension is 36, and

```text
g=D_b−D_a, h=D_a+D_b−J/6, J=ΣD_p,
G=g⊗I+I⊗g, H=h⊗I+I⊗h, B=H−iG,
T=[[0,B],[-B†,0]], W=exp(θT),
R=Tr_load(W† Z W)/2, Q=Φ_e(Z)=W(I_2⊗R)W†,
β=Tr_36(R²)/36, α=1−β.
```

C and E coincide on this input: their odd-block action is the conditional
expectation `X ↦ I_2⊗Tr_load(X)/2`. The even block is fixed. Locality and
occupation conservation force each seam contribution to the bare Dirichlet
matrix to be an incidence outer product. The odd-shell probability is
`20/66=10/33`; `N_a=(N_a+N_b+Z)/2`. Proper source covariance makes α equal
on all thirty seams. Therefore

```text
K_pq=τ[N_p(I−Φ)N_q]=α L_ico,pq/396,
Cov_τ(N)=(5/33)(I−11ᵀ/12),
P(I−Φ)P on centered linear loads = α L_ico/60.
```

Here P is the orthogonal projection onto linear occupations. The last line
is a compression, not yet an invariant-subspace claim. Direct expansion of
the endpoint joint law gives `E[ΔX ΔXᵀ]=2K`; hence `60K/Tr K=L_ico`.
Four times its smallest-positive-eigenvalue projector is the source Gram
`v_p·v_q/||v_p||²`, of rank three. Selecting this band is the existing finite
port readout, not a derivation of physical spacetime dimension.

Each endpoint's protected 3+3 splitting leaves four private blocks of
dimension nine. On each, R is a traceless Hermitian contraction: contraction
follows from normalized partial trace of the involution `W†ZW`, and its trace
vanishes because W preserves that block. Nine eigenvalues in [−1,1] summing
to zero have sum of squares at most eight. Maximize the convex sum on this
polytope: at a vertex eight entries are endpoints and the last must be zero.
Consequently **`1/9 ≤ α ≤ 1` for every real θ**. No singular coupling loses
the finite port geometry.

There is also a useful quantitative bound. The block exponential gives
`R=[f(BB†)−f(B†B)]/2`, `f(x)=cos(2θ√x)`. For x≥0, `|f'(x)|≤2θ²`.
In eigenbases of Hermitian A,C,
`||f(A)−f(C)||_F²=Σ_ij |f(a_i)−f(c_j)|² |<u_i,v_j>|²`, which bounds this
by `4θ⁴||A−C||_F²`. The source identity
`Tr_6(i[g,h])²=49+19√5` and `Tr_6(i[g,h])=0` thus give

```text
α ≥ max{1/9, 1−(4/3)(49+19√5)θ⁴}.
```

Its quartic coefficient is also the exact leading coefficient of `1−α`;
the asymptotic expansion is not used as a finite-angle certificate.

## Exact memory reduction and the decisive counterexample

For centered v with `L_ico v=λv`, define `f_v=Σv_pN_p` and
`c_n=τ[f_v Φⁿ(f_v)]/τ[f_v²]`. Self-adjointness of Φ and orthogonality of P
give, without an approximation,

```text
c_1=1−r, r=αλ/60,
c_2−c_1²=||(I−P)(I−Φ)f_v||²_τ/||f_v||²_τ = b_v.
```

Under Poissonization at rate one, the initial decay is r and the initial
curvature is `r²+b_v`. The ratio `b_v/r²` is unchanged by a common clock
rescaling. A closed native heat mode has b_v=0. Positive b_v excludes that
law on the same observables; it does not rule out a larger closed system
with its private degrees of freedom retained.

The full-network norm reduces to local overlaps. Define the edge residual
`S_e=Q_e−βZ_e`, extended by zero outside odd occupation, and let ∂ have
oriented edge rows. Each diagonal load block of S_e has zero private trace.
Its normalized Gram matrix B on the full reference is:

* `B_ee=(10/33)Tr_72(S_e²)/72=(10/33)β(1−β)`.
* Disjoint edges have zero entry.
* For edges meeting at p, reduce each diagonal load block by normalized
  partial trace over the other endpoint. Write the resulting six-dimensional
  matrices as `M_e^1,M_e^0` according to p's occupation. Then
  `B_ef=[9 Tr_6(M_e^1 M_f^1)+Tr_6(M_e^0 M_f^0)]/(66·6)`.

The weights count full-shell configurations: when p is occupied there are
nine places for the other particle; when p is empty the other two endpoints
must both be occupied. Distinct-edge off-diagonal load blocks have zero
product trace. Same-edge coherences remain in the complete local norm. Thus
the leakage quadratic form is exactly `∂ᵀB∂/3600`, and dividing by the bare
covariance 5/33 gives b_v for unit coefficient vectors. This is a locality
reduction of the original network, not an independent mean-field model.

At θ=1/8 and 1/16, the interval calculation certifies strict diagonal
dominance of B with positive diagonal. Since the carrier is connected,
**every nonconstant bare carrier vector has strictly positive leakage**.
For its rank-three band, approximate `(α,b_v,b_v/r²)` are respectively
`(0.975986704560, 4.178223344e−5, 0.0206705310)` and
`(0.998236412623, 3.298779751e−6, 0.0015600348)`. Positivity is certified;
these printed decimals are numerical summaries. Both one-term exchange
controls have normal B, hence R=0, α=1 and zero leakage exactly.

## A corrected single-port readout

The same result admits a local scalar readout without reconstructing a carrier
eigenmode. Keep the tracial preparation and IID uniform attempt schedule, fix
a port p, and put `f=N_p−1/6`, `c_n=τ[f Φⁿ(f)]/τ[f²]`, and `A=I−Φ`.
The carrier identities `L_pp=5`, `(L²)_pp=30`, the covariance above and
`τ[f²]=5/36` give

```text
r = 1−c_1 = α/11,
||P A f||²_τ / ||f||²_τ = α²/110,
b_local = ||(I−P) A f||²_τ / ||f||²_τ
        = c_2 − c_1² − (1−c_1)²/10.
```

This is a corollary of the existing compression and leakage proof. The
uncorrected `c_2−c_1²` contains `α²/1210` from mixing carrier bands;
native repair already gives `1/1210`, despite closed mean dynamics. The
corrected statistic is zero for native repair and the one-term controls.
The certified edge bound `B≥mI`, together with `||∂e_p||²=5`, gives
`b_local≥m/100`. Its certified lower bounds are `186661676/10¹⁴` at θ=1/16
and `3778430107/10¹⁴` at θ=1/8. Numerical evaluations are respectively
`4.14290705794e−6` and `6.05048265887e−5`.

Only N_p needs measurement at the two endpoints of each independently prepared
trial. If E_a is its outcome projector and P_x the full load-pattern
projectors, then `Σ_{x:x_p=a} P_x τ P_x = E_a τ E_a = τ E_a`. Thus local
endpoint measurements give the same marginal joint law as the all-load
protocol. No other port's record, intermediate measurement, private reset or
postselection is needed. Preparation and measurement access remain inputs.

The stopping rule is still one or two **global IID attempts**. Stopping at
observer-local pulses, delayed messages or state-dependent events defines
another instrument and requires its own null calculation. With a supplied
constant Poisson rate ν>0, `c(t)=τ[f exp(−νAt)f]/τ[f²]` instead gives

```text
c′(0)=−νr,       c″(0)=ν²((11/10)r²+b_local),
c″(0)/c′(0)² − 11/10 = b_local/r².
```

The ratio is invariant under a common constant change of time units, not
arbitrary retiming. It does not calibrate physical duration. C and E still
give the same first two moments, so this local statistic cannot select their
completion or uniquely identify N4 among larger-state laws.

A measured test must propagate preparation, readout, scheduling and sampling
errors through the complete corrected statistic. Its quadratic terms make a
finite-sample plug-in estimate biased; a positive point estimate alone is not
evidence of leakage. Trial allocation, shared uncertainty, confidence rule and
adequate power remain prospective obligations. The regressions evaluate
unprojected incident-edge moments independently of the leakage formula and
retain the exact native false-positive control; they are finite-model checks.

## Reproduction and scope

Run `python -m pytest -q -W error code/native_response` from the repository
root with the pinned requirements. CI runs this on Linux and Windows in the
existing finite quantum-information job. The tests bind the source chart,
compare independent Kraus and source formulas, recover the existing Gram,
check full-shell overlap counts, retain the simpler controls, and replay
the finite-coupling interval certificate.

The floating evaluators are diagnostics on `|θ|≤1`; they do not certify
zeros, tiny differences or arbitrary-angle forward accuracy. The separate
certificate starts from the algebraic source entries and exact rational
couplings. All-angle nondegeneracy follows from the proof above.
The finite common fixed algebra, absolute quantum gap, physical clock,
source-law selection and continuum limit are not inferred. In particular
C and E share these first two occupation moments, so this test alone does
not distinguish those two completions. No M01–M12 parent or #1076 is closed
by the result.

### Finite certificate and error budget

`PYTHONPATH=code python -m native_response.certificate` prints the two
recomputed rational enclosures (set `PYTHONPATH` using the local shell).
There is no user-supplied success packet to trust. The certificate separately
reconstructs the algebraic source using directed interval operations,
including √5, in an isolated 35-digit context. It evaluates degree-24
Taylor polynomials on each protected 18-dimensional load/private block.
The interval row sums prove `||θT||∞<2`; skew-Hermiticity also gives
`||θT||₂≤||θT||∞`. The entire exponential tail is bounded by
`e²2²⁵/25! < 8·2²⁵/25! < 2⁻⁵⁰`. Each real and imaginary matrix entry is
enlarged by that error before the nonlinear trace computations. Arithmetic,
source and truncation uncertainty therefore all reach the final enclosure.
The existing broad rational endpoint export checks outward containment.

The 120 unordered pairs of incident carrier edges form two orbits of sixty
under the proper carrier action. Representatives are `{(0,2),(0,4)}` and
`{(0,2),(0,8)}`. Symmetry of the Gram and endpoint orientation signs leave two
absolute off-diagonal entries, with four of each in every row; all other
off-diagonal entries vanish. The exact permutation enumeration is tested
alongside source covariance and the independently computed thirty-edge Gram.
If m denotes the diagonal minus those eight absolute entries, the certificate
gives the strict bounds

| Exact coupling | Certified lower bound on m |
| --- | --- |
| 1/16 | 186661676 / 10¹² |
| 1/8 | 3778430107 / 10¹² |

Thus `B≥mI` and the normalized leakage matrix is at least
`(11m/6000)L_ico`. This proves positivity for the whole centered subspace,
including arbitrary superpositions, without diagonalizing an enormous
network operator. The enclosure does not assert those rational lower bounds
are optimal, nor that nonclosure holds at every nonzero coupling.
