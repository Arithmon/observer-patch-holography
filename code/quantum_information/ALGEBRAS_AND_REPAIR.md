# When finite entropy completion is a repair channel

For a finite observable algebra and a faithful reference state, relative-entropy
completion of all its marginal states is affine **exactly when** the algebra
is invariant under the reference's modular flow. In that case the completion
is the dual of the unique reference-preserving conditional expectation.
Outside that case, no single linear channel implements all those completions.
This classifies the whole finite algebraic problem, including boundary states;
it is not a fit to selected successful inputs.

The result connects the finite A3 objective, the collar alignment criterion
and the conditional expectations in the
[canonical repair proposal](../../docs/CANONICAL_REPAIR_LAW_RFC.md).
It also identifies which premises are doing work. An observable span need
not be an algebra, a weighted projection need not be positive, and a selected
reference need not be compatible with a proposed retained algebra.
These are three different tests.

This note audits the live finite implementations on top of the
[entropy and sector-label audit](README.md), at predecessor commit
09170c3f00299bb61ffe62e1166255cb332759da. It uses ordinary matrix trace and
natural logarithms. The finite mathematical results are proved below;
the Python checks are floating-point diagnostics supplemented by exact
rational counterexamples. None of these checks identifies a simulator's
transition law merely by constructing the desired channel.

## Correct subsystem and separation calculations

Let \(H=H_R\otimes H_S\), with dimensions \(D=d_Rd_S\), where the tensor
slots in R can be noncontiguous. On the full Hilbert-Schmidt space \(HS(H)\),
use left multiplication \(L(a)v=av\) and the vector \(\Omega=\sqrt{\rho}\).
Then

\[
\|L(a_R\otimes I_S)\Omega\|^2
 =\operatorname{Tr}(\rho_R a_R^*a_R),\qquad
\|L(a_R\otimes I_S)\|_{\mathrm{HS}}^2
 =D d_S\operatorname{Tr}(a_R^*a_R).
\]

Consequently the separation modulus with the *ambient Hilbert-Schmidt
operator norm* used by the legacy helper is

\[
\inf_{\|L(a)\|_{\mathrm{HS}}=1}\|L(a)\Omega\|
 =\sqrt{\frac{\lambda_{\min}(\rho_R)}{D d_S}}.
\]

The infimum follows by choosing a rank-one \(a_R^*a_R\) in a minimum
eigenspace. In particular, the represented regional algebra separates
\(\Omega\) exactly when its actual reduced state is faithful. For the full
algebra, both cyclicity on full \(HS(H)\) and separation hold exactly when
\(\rho\) is faithful: right multiplication by \(\sqrt{\rho}\) is then
invertible. If its rank is r, its cyclic span instead has dimension Dr.

For singular \(\rho\), that proper cyclic subspace is the actual GNS
representation of the state. The GNS vector is cyclic there by construction.
Calling all of \(HS(H)\) the state's GNS space and then calling its GNS vector
noncyclic confuses two representations. The faithful case agrees with both
descriptions.

For a general independent operator basis \(b_i\), set
\(G_{ij}=\operatorname{Tr}(b_i^*b_j)\) and
\(H_{ij}=\langle b_i\Omega,b_j\Omega\rangle\).
The squared modulus is the least generalized eigenvalue of \(Hc=\lambda Gc\).
It is not generally the least eigenvalue of H divided by one basis norm.
Orthonormalizing the operator span gives the same calculation without assuming
orthogonality or equal norms.

| Reproduced old behavior | Correction |
| --- | --- |
| For \(\rho=\mathrm{diag}(1/2,0,1/2,0)\), dimensions (2,2), region [1], the helper returned separating. The nonzero operator \(I\otimes|1\rangle\langle1|\) annihilates \(\Omega\). | Use the actual region's partial trace. The A region separates; the B region does not. The reversed example is also tested. |
| Bases [I,Z] and [I,I+Z] for the same diagonal qubit algebra gave moduli 0.7071067811865474 and 0.43701602444882093 at \(\Omega=(1,1)/\sqrt2\). | Both now give \(1/\sqrt2\), also under unequal rescaling. An independent generalized eigensolver checks the result. |
| The explicit positive spectrum (1,1.5e-20) gave cyclic but not separating. | Use one resolved spectral-support criterion for both; do not discard small positive eigenvalues at an arbitrary 1e-10 threshold. |
| Full standardness constructed \(D^2\) dense left-action matrices of size \(D^2\) by \(D^2\). At D=64 their complex entries alone need one tebibyte. | Use the D by D state and its spectrum. A D=64 test disables both dense representation constructors and checks full and noncontiguous regional standardness. |

The state validator still rejects invalid densities. Dense spectra whose
support cannot be resolved at floating-point precision raise an explicit
error; a numerical answer does not certify exact rank. Explicit diagonal
positive entries, including the smallest binary64 subnormal, remain positive.

## The complete finite completion theorem

Let \(B\subset M_D(\mathbb C)\) be a unital complex *-subalgebra and
\(\rho>0\) a density matrix. A state b on B means **all** its observable
expectations, not just a chosen nonclosed list of observables. Define

\[
C_\rho(b)=\arg\min_{\tau\ge0,\ \operatorname{Tr}\tau=1,\ \tau|_B=b}
                   D(\tau\|\rho).
\]

The feasible set here is the entire fiber of states with B data b.
Every such fiber is nonempty and compact. Since \(\rho\) is faithful, the
objective is continuous and strictly convex, so its minimizer exists and
is unique, including on boundary fibers.

**Theorem.** The following are equivalent:

1. \(C_\rho\) is affine on the entire state space of B.
2. B is invariant under \(x\mapsto\rho^{it}x\rho^{-it}\) for every real t.
3. There is a \(\rho\)-preserving unital completely positive projection E
   with range B, fixing every element of B.
4. In Wedderburn coordinates
   \[
   H=\bigoplus_\alpha(H_{L,\alpha}\otimes H_{R,\alpha}),\qquad
   B=\bigoplus_\alpha(M_{n_\alpha}\otimes I_{m_\alpha}),
   \]
   the reference has the form
   \[
   \rho=\bigoplus_\alpha
            p_\alpha\eta_\alpha\otimes\xi_\alpha,
   \quad p_\alpha>0,\quad \eta_\alpha>0,\quad\xi_\alpha>0,
   \]
   with normalized sector weights and normalized factor states.

When they hold, E is unique and

\[
E(x)=\bigoplus_\alpha
  \operatorname{Tr}_{R,\alpha}[(I\otimes\xi_\alpha)x_{\alpha\alpha}]
            \otimes I_{R,\alpha},
\qquad
E_*(\sigma)=\bigoplus_\alpha b_\alpha\otimes\xi_\alpha,
\]

where \(b_\alpha=\operatorname{Tr}_{R,\alpha}\sigma_{\alpha\alpha}\) are
subnormalized left marginals. E acts on observables; its Hilbert-Schmidt
adjoint \(E_*\) is the trace-preserving state channel. E itself need not
preserve ordinary trace. The completion is \(C_\rho(\sigma|_B)=E_*(\sigma)\).

### Proof of necessity, including affinity

Any positive affine assignment of full states to **every** state of B,
consistent with those B data, must have the product form
\(b\mapsto\bigoplus_\alpha b_\alpha\otimes\xi_\alpha\) with fixed right
states. To see this, first take data supported in one central sector.
Positivity forces the output to have zero support in the other sectors,
including cross-sector coherences. A pure left marginal then forces the
output to be \(|u\rangle\langle u|\otimes\xi_u\).

For orthogonal u,v, compare the two decompositions of
\((|u\rangle\langle u|+|v\rangle\langle v|)/2\), using u,v and
\((u+v)/\sqrt2,(u-v)/\sqrt2\). Affinity and the diagonal blocks give
\(\xi_u=\xi_v=(\xi_++\xi_-)/2\); the off-diagonal block gives
\(\xi_+=\xi_-\). Within any orthonormal basis all right states therefore
agree. The assigned state of \(I/n_\alpha\) is independent of its
decomposition, so this common right state is independent of the basis.
Every pure state lies in such a basis; spectral decomposition gives the
formula for every mixed left state. For \(n_\alpha=1\) there is only one
normalized left state. Finally, affinity over central weights gives the
formula on all of B's state space.

Apply this to \(C_\rho\). At the reference's own data it returns \(\rho\),
since \(D(\tau\|\rho)\ge0\), with equality only at \(\tau=\rho\).
Thus affinity forces exactly the block product reference in item 4.
This argument assumes positivity and full-domain affinity, not complete
positivity. It is the finite assignment argument underlying
[Pechukas's assignment-map result](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.73.1060);
the extension over the central sectors is explicit here.
The arbitrary-system assignment result is also proved by
[Rodriguez-Rosario, Modi and Aspuru-Guzik](https://arxiv.org/abs/0910.5568).

For modular invariance, a continuous automorphism group of B fixes each
minimal central projection: it cannot continuously permute a finite set.
Therefore \(\log\rho\) has no cross-sector blocks. In a sector, its
commutator defines a derivation of \(M_n\otimes I_m\). Such a derivation is
inner on \(M_n\): equivalently, examining matrix-unit commutators shows that
off-diagonal m by m blocks of \(\log\rho\) are scalar, and differences
between its diagonal blocks are scalar. Hence
\((\log\rho)_\alpha=K_{L,\alpha}\otimes I+I\otimes K_{R,\alpha}\).
Exponentiation and normalization give item 4. The converse follows by
direct conjugation. This is the finite form of
[Takesaki's conditional-expectation theorem](https://www.sciencedirect.com/science/article/pii/0022123672900043).

For completeness, item 3 also forces item 4 without assuming the theorem.
Since E fixes B pointwise and is unital completely positive, equality in
the Schwarz inequality places B in its multiplicative domain:
\(E(bxc)=bE(x)c\). Its dual consequently depends only on the restriction
of the input state to B, and is a positive affine consistent assignment
on the entire state space. The preceding assignment argument applies.
Preservation of \(\rho\) then gives item 4.

### Sufficiency, uniqueness and exact entropy accounting

For the product reference, the displayed E is a composition of central
pinching, positive-state partial evaluation and inclusion. It is unital and
completely positive, fixes B, and preserves \(\rho\). For any expectation
with these properties, bimodularity gives

\[
\operatorname{Tr}(\rho b^*E(x))
 =\operatorname{Tr}(\rho E(b^*x))
 =\operatorname{Tr}(\rho b^*x),\qquad b\in B.
\]

It is therefore the orthogonal projection onto B for the faithful GNS
inner product \(\langle x,y\rangle_\rho=\operatorname{Tr}(\rho x^*y)\),
which proves uniqueness.

Write \(\Gamma=E_*(\sigma)\). For faithful marginals,
\(\log\Gamma-\log\rho\) belongs to B: the right-factor logarithms cancel.
Since \(\sigma\) and \(\Gamma\) have the same B expectations,

\[
D(\sigma\|\rho)
 =D(\sigma\|\Gamma)+D(\Gamma\|\rho).
\tag{1}
\]

The same equation holds on boundary fibers. Let Q be the direct sum of
the supports of \(b_\alpha\), tensored with each full right factor.
Positivity and the zero marginal weights force \(\sigma=Q\sigma Q\),
including when \(\sigma\) has central coherences. Since every \(\xi_\alpha\)
is faithful, Q is exactly the support of \(\Gamma\). On Q, the right
logarithms still cancel, and the compressed left logarithms define an
observable in B after extending it by zero outside Q. Equality of its
expectations proves (1). No logarithm of zero or assumed invertibility of
the marginal is needed. The reference term is \(Q(\log\rho)Q\); it must
not be replaced by the logarithm of \(Q\rho Q\). The two differ when Q
does not commute with \(\rho\), a case included in the regression suite.

Every \(\tau\) in the same B fiber has the same \(\Gamma\). Equation (1)
then proves that \(\Gamma\) is its unique minimum, since
\(D(\tau\|\Gamma)\ge0\) with equality only at \(\tau=\Gamma\). This proves
item 4 implies item 1 and the stated completion formula. It also shows
exactly how much relative entropy one completed repair removes. Conditional
expectations and their entropy identities are treated in
[Bardet, Capel and Rouze, sections 2 and 3](https://link.springer.com/article/10.1007/s00023-021-01088-3).

The earlier audit's labelled direct sum turns the fixed-weight finite A3
objective into one relative entropy. The present theorem applies to such
an objective when its feasible set is precisely a complete B-data fiber.
An arbitrary OPH feasible family K may impose further constraints or lack
a joint matrix-state realization. This theorem neither removes those
constraints nor constructs that realization. In particular, the finite
classification does not adopt the proposed A1-R or A2-R clauses.

## Three exact controls that distinguish the premises

**Nonlinear completion.** Take B to be the diagonal qubit algebra and
\(\rho=\begin{pmatrix}1/2&1/4\\1/4&1/2\end{pmatrix}\).
The diagonal data (1,0) and (0,1) each have only one positive completion.
Their mean completion is \(I/2\). At their mean data, however, the unique
relative-entropy minimum is \(\rho\), because it already has those data.
The affine candidate loses by
\(D(I/2\|\rho)=\tfrac12\log(4/3)>0\). Thus a reference alone does not make
A3 completion a linear repair channel. There is no optimizer in this test.

**A weighted projection need not be positive.** Let
\(\rho=(I+\tfrac12 X\otimes X)/4\) and \(B=M_2\otimes I_2\).
The \(\rho\)-GNS orthogonal projection is unital, idempotent, fixes B,
has range B, preserves \(\rho\), and is GNS self-adjoint. Yet it sends the
positive input \(I+Z\otimes X\), whose eigenvalues are 0 and 2, to
\(I+\tfrac12 ZX\otimes I\), which is not even Hermitian. Here X and Z are
Pauli matrices. For the projected left factor y, the orthogonality equations
are \(y\rho_L=\operatorname{Tr}_R(x\rho)\); here \(\rho_L=I/2\), giving
the stated image. The tests check those equations with exact rational
matrices against a full Pauli basis. This prevents treating weighted
least squares as an automatic quantum channel.

**An operator system need not support a quantum repair.** Retain only the
span of I,X,Z with tracial reference \(I/2\). Entropy maximization sets the
unconstrained Y coefficient to zero. The resulting map
\(T(a)=(a+a^{\mathsf T})/2\) is affine, positive, unital and idempotent,
but its unnormalized Choi eigenvalues are
\((-1/2,1/2,1/2,3/2)\). It is not completely positive. Its range is closed
under adjoint but not multiplication, so it is not B as required above.
The algebra validator rejects that range before any repair construction.

## Noncommuting primitive repairs still have a common limit

Let \(E_m\) be finitely many reference-preserving expectations on one
common finite workspace, with strictly positive supplied rates \(\gamma_m\).
Define

\[
\mathcal L=\sum_m\gamma_m(E_m-I),\qquad B_\cap=\bigcap_m B_m.
\]

Each \(E_m\) is a GNS orthogonal projection. Therefore

\[
\langle x,-\mathcal Lx\rangle_\rho
 =\sum_m\gamma_m\|(I-E_m)x\|_\rho^2,\qquad
\ker\mathcal L=B_\cap.
\tag{2}
\]

Positivity of every rate makes the kernel equality exact. No commutation
between different \(E_m\) is needed. Put
\(\gamma=\sum_m\gamma_m\) and
\(P=\sum_m(\gamma_m/\gamma)E_m\). Then
\(e^{t\mathcal L}=e^{-\gamma t}\sum_{k\ge0}(\gamma t)^kP^k/k!\) is unital,
completely positive and reference-preserving. The finite spectral theorem
and (2) show convergence to the GNS projection \(E_\cap\) onto \(B_\cap\).
That intersection is modular invariant, so the projection is its unique
reference-preserving conditional expectation.

When the orthogonal complement is nonzero, let \(\lambda>0\) be the least
eigenvalue of \(-\mathcal L\) on it. Then
\(\|(e^{t\mathcal L}-E_\cap)x\|_\rho
\le e^{-\lambda t}\|(I-E_\cap)x\|_\rho\).
When the intersection is the full algebra there is no positive gap to
report. The dual semigroup is a state channel, so relative entropy to
\(\rho\) is nonincreasing by data processing. This does not identify the
GNS gap with an entropy decay constant.

An exact control uses qubit dephasing in Z and
\(N=(3X+4Z)/5\), each at rate 1/2 with \(\rho=I/2\).
These expectations do not commute. Their common algebra is the scalars
and the generator spectrum is \(0,-1/10,-9/10,-1\).
SymPy verifies the rational characteristic spectrum independently of the
numerical gap routine. Explicit matrix exponentials check complete
positivity, entropy decrease and the corresponding norm bound.
This supplies a concrete noncommuting example, not just a commuting
special case.

## Implementation, independent checks and downstream effects

[algebras.py](algebras.py) validates a complete independent basis, its
identity, adjoints and every pairwise product. It embeds operators in the
declared tensor slots, computes intersections, and tests modular invariance
by all commutators with \(\log\rho\). This finite generator criterion is
equivalent to invariance under all modular times.

[expectations.py](expectations.py) separates the weighted projection
candidate from a validated expectation. Its recognizer accepts the **full
linear superoperator**, not a callable tested on a few random states.
It independently recomputes the full Choi matrix, unitality, reference
preservation, range, fixation of B, idempotence and GNS self-adjointness.
It never decides acceptance by comparing a supplied map with the producer's
projection formula. The observable map and its dual use explicit row-major
vectorization; complex unitary examples test the convention.

The numerical tolerance is finite, positive and bounded; it cannot be set
to infinity or a Boolean to disable checks. Invalid states, singular
references, dependent bases, empty families, omitted/zero/negative rates,
nonfinite matrices, diagnostic overflow and unresolved Gram matrices fail
explicitly. Validation remains active under optimized Python.
Vector norms are evaluated in complex128 arithmetic before checking
normalization. At audit commit df1d99e5, the int64 vector \((2^{32},1)\)
had a wrapped squared norm of one and was incorrectly accepted; float32
could also round the squared norm of \((1,10^{-4})\) to one. Regression
controls reject both, along with unsigned-integer and complex64 variants.
These guards are not interval arithmetic: deviations below the declared
tolerance can pass. In particular, intersection dimension and spectral gaps
are numerical diagnostics, not exact certificates. Almost coincident
algebras can have arbitrarily small gaps. The routine uses a machine-scale
rank threshold, rejects unresolved positive gaps, and retains a 1e-8-angle
negative control; sub-resolution distinctions still need exact input or
certified enclosures. A numerical construction is not a source-law receipt.

| Existing component | Effect of this audit |
| --- | --- |
| Null-net standardness, #524 | Fixes actual-region selection, basis normalization, small-eigenvalue inconsistency and dense allocation. Published faithful and single-sector controls remain valid. The legacy full-standardness helper can now raise on unresolved support instead of silently deciding rank. |
| Collar modular alignment, #543 | Its complete-basis criterion now has a cross-check against existence of the actual finite expectation. Aligned product blocks pass and generic incompatible blocks fail. The independent collar calculation is retained. |
| Finite MaxEnt/A3 | Identifies exactly which full-marginal problems give affine repair and proves the Pythagorean accounting; arbitrary feasible families retain their own constraints. |
| Canonical repair RFC | Supplies the finite existence/recognition calculation and noncommuting convergence proof for specified algebras and reference. Protected records, conserved functionals, checkpoint instruments, grammar completeness, refinement and clock selection are additional RFC obligations. |
| Classical consensus projection | A separately specified row-stochastic three-state transition table reproduces weighted fiber resampling. Replacing stationary weights by uniform weights fails. The existing demand for independent transition extraction is retained. |
| Source operator-join evidence | Its exact source-specific verifier remains independent. A generic constructed channel does not replace its source identification or Lean proof. |

The one-particle zero-extension helper in the null-net receipts is also a
different operation from tensor-identity embedding and remains separate.
The change does not recalculate frozen receipts, amend the paper or promote
a physical model. It supports the existing-theory audit alongside the
bounded physical-candidate issues #1025/#1026.

The regression suite uses explicit tensor indices, independent generalized
Gram eigenvalues, weighted partial traces, independent Choi blocks, SciPy
matrix logarithms, exact rational controls and adversarial supplied maps.
A central-sector test includes a pure state with cross-sector coherence
and a singular completed marginal, checking (1) against scalar logarithms.
Separate controls treat a singular marginal whose support does not commute
with the reference, and noncommuting repairs preserving a complex nontracial
reference after a global unitary change of coordinates.
Run the live consumers and the new tests with the pinned requirements:

~~~text
PYTHONPATH=code OPENBLAS_NUM_THREADS=1 python -m pytest -q \
  code/quantum_information code/geometry code/collar_alignment code/maxent
~~~

The dedicated finite-quantum-information workflow runs these tests on both
Linux and Windows.
