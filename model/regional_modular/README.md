# Regional expectations and a common modular reference

The [signed CAR regional maps](../field_regions/README.md) give complete
state restrictions and ambient trace-preserving conditional expectations.
The [finite field Gibbs calculation](../field_stress/README.md) compares the
actual reduced modular generator with bare regional energy. This package
reuses finite conditional-expectation theory to state the missing common
reference condition exactly. Its main field witness is rational: boundary
hopping forbids a Gibbs-preserving UCP retraction onto the stated regional
algebra, even if ordinary trace preservation is relaxed.

This does not invalidate the regional state trace, the existing tracial
expectation, or a clock defined from the actual reduced density. It does
not exclude approximate, non-idempotent or operationally restricted maps.
The action, coordinates, reference choice and Hamiltonian model time are
supplied. No physical particle, mass, geometry or continuum is identified.

## The finite algebra and the exact criterion

Let

```
H = direct_sum_n (A_n tensor B_n),
M = direct_sum_n (B(A_n) tensor I_(B_n)) subset B(H),
```

with every displayed factor nonzero. In the actual fixed-total-number CAR
sector, `A_n = wedge^n C^U` and `B_n = wedge^(N-n) C^(U-complement)`.
The signed full-Fock split supplies this direct sum; the global fixed-N
space is not replaced by an unrestricted tensor product. Absent sectors
are omitted. The central projection `P_n` records regional particle number.

An expectation here is a **unital completely positive retraction**
`E:B(H)->M`, with the inclusion into `B(H)` understood: it fixes every
element of `M`. Its range condition implies idempotence. A faithful density
`rho>0` is preserved when `Tr(rho E(X))=Tr(rho X)` for every matrix `X`.
This is preservation of a state functional by a Heisenberg map; it is
different from ordinary ambient trace preservation.

**Finite criterion.** Such a `rho`-preserving expectation exists if and
only if

```
rho = direct_sum_n (r_n tensor tau_n),
r_n > 0, tau_n > 0, Tr(tau_n)=1, sum_n Tr(r_n)=1.
```

The sector weights `Tr(r_n)` can be arbitrary positive numbers. This is a
conditional product in each sector, not a global product of unconstrained
Fock densities. When it exists, the expectation is unique:

```
E(X)_n = [(id tensor phi_tau_n)(P_n X P_n)] tensor I_(B_n),
phi_tau_n(Y) = Tr(tau_n Y).
```

**Proof.** A UCP map obeys Kadison--Schwarz. If a projection or matrix `a`
and both `a* a` and `a a*` are fixed, equality puts `a` in the map's
multiplicative domain: `E(aX)=aE(X)` and `E(Xa)=E(X)a`. For completeness,
in a finite Stinespring form `E(X)=V* pi(X) V`, `V*V=I`, the Schwarz
difference is
`V* pi(X)* (I-VV*) pi(X) V`. Equality makes the middle defect zero;
applying this to `a` and `a*` gives both bimodule identities.

Every `a` in `M` meets these equality conditions because `M` is a
*-algebra fixed pointwise. For distinct central sectors,
`E(P_n X P_m)=P_n E(X) P_m=0`. State preservation therefore forces the
off-diagonal blocks of `rho` to vanish. Within a block, the bimodule
identities force `E(I tensor Y)` to commute with every `B(A_n)` element;
since it belongs to `M`, it equals `phi_n(Y) I`. Positivity and unitality
make `phi_n` a state, represented by one density `tau_n`.
The same identities give the displayed formula for all matrix units.
Preservation of `rho` then identifies its block as `r_n tensor tau_n`,
where `r_n=Tr_B(P_n rho P_n)`. Faithfulness makes both factors positive
definite and fixes `tau_n` uniquely.

Conversely, the formula fixes `M` and preserves the displayed state.
For `tau_n=sum_i lambda_i |v_i><v_i|`, its Heisenberg Kraus operators are
`sqrt(lambda_i) I_A tensor |v_i><j|` in block `n`, extended by zero
elsewhere. Their `K* K` sum is the full identity. They establish complete
positivity, including arbitrary untouched reference factors. This proves
the criterion without a continuum modular theorem or a fitting procedure.

## Consequences for clocks and regional composition

For a finite self-adjoint Hamiltonian `H` and
`rho_beta=exp(-beta H)/Z` with **beta nonzero**, the criterion is equivalent
to the blockwise additive form

```
[H,P_n]=0,
H restricted to P_n = H_A,n tensor I + I tensor H_B,n.
```

Sector scalars can be absorbed into either summand. This follows by taking
the logarithm of the positive product blocks, and conversely exponentiating
commuting tensor summands. The resulting regional algebra is invariant
under the full Gibbs modular flow, and its expectation intertwines that
flow. At beta zero the state is tracial and preservation imposes no
condition on a separately chosen Hamiltonian; it cannot identify that
Hamiltonian's time with a zero modular generator.

The center gives a cheap necessary test. If `N_U=sum_n n P_n`, existence
requires `[H,N_U]=0`. This condition is **not sufficient**: a cross-factor
interaction inside one number block can leave the center fixed while its
Gibbs density is correlated. The control `rho=diag(3,1,1,3)/8` on a
two-by-two factor has that property. No factor expectation preserves it.

One compatible faithful reference also simplifies tower bookkeeping.
If `M_A subset M_B` both admit expectations preserving the same `rho`,
their composition `E_A E_B` is a UCP retraction onto `M_A` preserving
`rho`. Uniqueness gives `E_A E_B=E_A`, while `E_B E_A=E_A` follows from
the inclusion. Independent regional references need not give a tower;
the tests retain a three-factor countercontrol with mismatched exterior
states. Selecting each region's density separately is insufficient.

The existing [CPRegionalNet](../../Lean/QFT/CPRestrictionNet.lean) fixes
**ordinary trace preservation**. Its actual CAR witness uses
`tau_n=I/dim(B_n)` and is the Hilbert--Schmidt orthogonal projection.
Such expectations always compose for nested regions. Requiring this same
map also to preserve `rho` restricts the criterion further to
`rho_n=r_n tensor I/dim(B_n)`. A nontracial weighted expectation is therefore
not automatically an inhabitant of that unchanged Lean interface.
The code supplies no new Lean inhabitant or kernel proof.

There is a useful distinction between preserving a Gibbs reference and
intertwining the Hamiltonian generator. The tracial expectation commutes
with `ad_H(X)=[H,X]` exactly when `H` has the blockwise additive form above;
its exterior energy need not be scalar. To see necessity, commutation
makes `ad_H` preserve `M`. Applying it to each `P_n` forces the off-diagonal
blocks of `H` to vanish. Within a block, let
`H_A=(id tensor normalized_trace_B)(H)`. If `[H,X tensor I]` lies in the
factor, its partial trace is `[H_A,X]`, so
`H-H_A tensor I` commutes with the entire local matrix factor and equals
`I tensor H_B`. Conversely, the partial trace kills the exterior
commutator and intertwines the local one. Thus a tracial expectation can
intertwine an additive Hamiltonian while failing to preserve its Gibbs
state; a control checks both properties on the same matrices.

## Exact boundary witness on the supplied field

For `N=2`, the field Hamiltonian is `H=dGamma(h)+U P_contact`. The contact
term commutes with every regional number. Put
`D_U=diag(1_(p in U))` on the one-particle space of dimension `m`. Then

```
[N_U,H] = dGamma([D_U,h]),
||[N_U,H]||_HS^2 = (m-2) ||[D_U,h]||_HS^2
                 = (m-2) sum_(p,q crossing U) |h_pq|^2.
```

The last sum is over **ordered** matrix entries: both `p,q` and `q,p`
are counted. To prove the norm identity directly, the one-body commutator
has zero diagonal. Each coefficient `A_pq`, `p!=q`, sends the N=2 basis
state `{q,r}` to `{p,r}` for each of the `m-2` distinct spectator modes.
The CAR sign has modulus one, and a pair of input/output configurations
determines `p,q,r` uniquely, so no different coefficient contributes to
the same matrix entry. Summing the squared entries proves the formula.

On the canonical open radius-one Weyl cube, take both spin modes at the
origin and its positive-x neighbor. Here `m=54`, the N=2 dimension is
1,431, and there are 36 ordered nonzero crossing entries, each of modulus
`1/16`. Thus

```
||[D_U,h]||_HS^2 = 36/256 = 9/64,
||[N_U,H]||_HS^2 = 52*36/256 = 117/16 > 0.
```

The [certificate](retraction.py) checks integer real and imaginary parts
of `16h` and `16H`, independently counts both commutators, and compares
their exact rational norms. It is identical for `U=0` and `U=1`.
Consequently no expectation onto this regional number-preserving algebra
can preserve the full-field Gibbs state at nonzero beta. Computing its
1,431-dimensional density or taking a numerical logarithm is unnecessary
for that conclusion. Every UCP retraction onto this algebra fixes `N_U`
and kills the off-diagonal number blocks of `[H,N_U]`. Its generator
intertwining defect on `N_U` therefore has the same squared norm `117/16`.
Cutting all boundary hopping for this full-site region is a positive
control; the sector-product Gibbs density then admits the weighted
expectation.

This is a finite exact boundary for one algebra/reference/action join.
The actual reduced density and its modular Hamiltonian remain legitimate.
Boundary corrections, changed references, non-idempotent channels,
different accessible domains and controlled approximation remain possible.
Nothing here supplies Lorentz stress, the physical Ward-limit hypotheses,
an Einstein-register inhabitant, or particle mass calibration.

## A merely operator-system range does not evade the exact faithful case

Suppose a UCP idempotent `E` preserves a faithful finite density, without
initially requiring its fixed range to be an algebra. For a fixed matrix
`X`, Schwarz gives `E(X*X)-X*X>=0`. Its expectation in the faithful density
is zero by invariance, so this positive difference is zero. The same holds
for `XX*`. The multiplicative-domain argument makes products of fixed
matrices fixed; their adjoints and identity are fixed too. Thus the range
is automatically a C*-algebra. Merely changing the target's name to an
operator system does not remove the exact retraction condition.

Faithfulness is essential. The UCP projection on `M3` given by
`E(X)=diag(X_[0:2,0:2], Tr(X_[0:2,0:2])/2)` preserves
`diag(1/2,1/2,0)`, but its range is not closed under ordinary products:
`diag(1,-1,0)` is fixed and its square is not. The test gives its complete
Kraus realization. Non-idempotent observations and approximate preservation
also fall outside the exact argument.

## Reproduction and scope

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/regional_modular
python -m model.regional_modular.retraction
```

The tests compare the block map with independent Kraus operators and the
actual signed CAR expectation, use complex reference controls, retain
sector-weighted positive states and crossing/within-block negative controls,
and check the zero-beta and nonfaithful exceptions. Dense eigenvalues and
matrix exponentials are diagnostics; the exact field witness and the stated
finite proofs carry the mathematical conclusions. No source apparatus or
physical density preparation is asserted.
