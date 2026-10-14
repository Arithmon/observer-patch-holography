# The declared repair law on the Euclidean history register: inertia, the golden refinement family and the reference state

This note accompanies a Lean library about a Langevin law declared on the Euclidean history register
of Observer Patch Holography (OPH), whose drift is the mean of OPH's capacity-weighted seam
agreement. It contributes three connections to rows of OPH's premise map, each a conditional
candidate under named premises, with every input labelled derived, candidate hypothesis,
calibration, target-exposed choice or supplied. For inertia, the two-time correlations of the
stationary chain fix the reconstructed acceleration matrix, whose slab-thickness limit has on OPH's
three-site control the distance-two entry OPH records; read through OPH's calibration fit
$`\alpha I+\beta A`$, this response distinguishes the same-equilibrium counterfamily exactly when
the whitened stiffness has at least three distinct eigenvalues, as at golden levels $`q\ge3`$ but
not on OPH's two-site action. For the continuum, at every golden level the reconstructed
coefficients and vacuum covariance tend to those of OPH's canonical quantization, with
vacuum-covariance error at most a quarter of the slab thickness per mode, so with OPH's continuum
proposition the limit is joint for its smearings. For the reference state, the declared dynamics'
iterated limit has, at unit $`\hbar`$, the position covariance of the vacuum OPH supplies or
targets. Nothing here selects a physical vacuum or identifies physical time.

## 1. Scope and status

In a paragraph labelled *[Lean]* every claim is a theorem of the Lean library `Lean/SourceLaw/`
(modules `SourceLaw.*`), named in parentheses, and its definitions and attributions only fix
notation and sources (Lean 4.29.1, Mathlib `5e932f97`, built from scratch with no errors, no
warnings and no `sorry`; an axiom-closure check over all constants of the library finds only
`propext`, `Classical.choice` and `Quot.sound`). A paragraph labelled *[Lean, with cited
hypotheses]* applies a named theorem whose hypotheses the preceding *[Cite]* paragraph supplies. The
library has 23 modules; sixteen of them carry the operator-level results of the declared law (its
stationary covariance, the transfer structure and reflection positivity of the stationary chain, and
what determines the inertia), and a paragraph or table entry labelled *[Lean, operator level]*
restates results of those sixteen modules under the same rule as *[Lean]*. A paragraph labelled
*[Std]* cites a standard result whose hypotheses we check against the Lean statements, and *[Cite]*
restates a result of OPH's documents unchanged. Statements about OPH refer to this repository at
revision `33a2eb0d`. Each connection follows OPH's result format: a conditional candidate under
named premises, with every input labelled derived (computed here), candidate hypothesis, calibration
or target-exposed choice (a choice the calculation leaves to the target), inputs that OPH itself
declares labelled supplied, and the physical choices the calculation does not determine stated as
such.

## 2. The declared law and the operator-level results

On a finite carrier without self-seams, with capacities $`M>0`$, seam conductances $`c\ge0`$ and
grounding $`g>0`$, the Euclidean history register has cells $`(i,t)`$, $`t\in\{0,\dots,T\}`$,
spatial seams $`\Delta\tau\,c_e`$ in every slab, grounds $`\Delta\tau\,g_i`$ and time seams
$`J_i/\Delta\tau`$; let $`P`$ be its stiffness and $`\mathbb{M}=\mathrm{diag}(M)`$ on cells. The
declared law is

```math
\mathrm{d}f=-\mathbb{M}^{-1}Pf\,\mathrm{d}s+\sqrt{2\hbar\,\mathbb{M}^{-1}}\,\mathrm{d}\mathcal{W}_s \tag{1}
```

with sampling time $`s`$, which is not physical time; here $`V`$ is the stiffness of the carrier
(its seams $`c`$ and grounds $`g`$). *[Lean, operator level]* Its drift matrix $`\mathbb{M}^{-1}P`$
is the mean repair generator, the sum over register seams of OPH's capacity-weighted agreement with
declared grounds (`repairGenerator_apply`, `register_drift`); its Lyapunov equation has the unique
solution $`\hbar P^{-1}`$ (`lyapunov_stationary`, `lyapunov_unique`, `register_stationary_unique`);
for $`J=M`$ the stationary chain with the register's interior precision has transfer eigenvalues
$`r(\lambda)`$ and capacity-whitened, unit-$`\hbar`$ equal-time covariance $`C_0=c_0(A)`$
with $`A=M^{-1/2}VM^{-1/2}`$
(`transferEig_add_inv`, `equalTimeVar_closed`, `chain_covariance_row`, `contact_identity_mode`,
`transferMat_add_inv`, `equalTimeCov_mul_one_sub`), its covariance is bond reflection positive
(`chain_bond_rp`, `block_bond_rp`), and the coefficients of the Hamiltonian reconstructed from it
tend, mode by mode, to those of OPH's canonical completion as $`\Delta\tau\to0`$
(`latticeFreq_tendsto`, `kineticCoeff_tendsto`, `potentialCoeff_tendsto`; the reconstruction itself
is the standard lattice Osterwalder–Schrader construction with Mehler's formula).

Throughout the transfer and two-time reader statements, correlations are normalized by the
supplied $`\hbar>0`$ (equivalently, those statements use $`\hbar=1`$). In the original carrier
coordinates, the raw equal-time covariance is $`\hbar M^{-1/2}C_0M^{-1/2}`$. If
$`\widetilde C_j=\hbar C_j`$ are raw stationary-chain correlations in the reader's coordinates,
its calibrated formula is
$`\hbar\Delta\tau(\widetilde C_0\widetilde C_1^{-1}\widetilde C_0-\widetilde C_1)^{-1}=J`$.
Without the factor $`\hbar`$, applying that reader to raw correlations returns $`J/\hbar`$.

With uniform capacities, seam rates $`1/N`$ for $`N`$ seams and no grounds, the mean repair
generator is $`L_s/(2N)`$ for the seam Laplacian $`L_s`$ (the corresponding case of
`repairGenerator_apply`), the generator of OPH's uniform seam average $`I-L_s/(2N)`$
([`model/common_source/NATIVE_FIELDS.md`](../../model/common_source/NATIVE_FIELDS.md)); OPH's Lean
library machine-checks the equal-capacity pair average and that uniform seam average
([`Lean/ObserverPatchHolography/ScalarSeamRepair.lean`](../../Lean/ObserverPatchHolography/ScalarSeamRepair.lean))
and, partly through `native_decide`, the diffusion kernel of the twelve-port repair mean
([`Lean/Screen/PortGramRepairCovariance.lean`](../../Lean/Screen/PortGramRepairCovariance.lean)).

Write a quadratic Hamiltonian as $`\tfrac12(p^\top\mathcal{T}p+q^\top Vq)`$ with inverse inertia
$`\mathcal{T}`$ and potential $`V`$; for real symmetric positive-definite $`G,\mathcal{T},V`$
and $`\hbar=1`$, a centred Gaussian $`\exp(-q^\top Gq/2)`$ is
its ground state exactly when $`G\mathcal{T}G=V`$ (a standard computation outside Lean), and its
position covariance is $`C=(2G)^{-1}`$. The Lean statements about $`G`$ concern only the relation
$`G\mathcal{T}G=V`$.

*[Lean, operator level]* For $`\beta\neq0`$ and invertible $`V`$, the configuration block of the
inverse of the Gibbs precision $`\mathrm{blockdiag}(\beta V,\beta\mathcal{T})`$ is
$`(\beta V)^{-1}`$ for every invertible $`\mathcal{T}`$, so configuration data do not determine the
inertia (`config_cov_indep`, which holds by construction); every $`\mathcal{T}`$ satisfies the
relation for the potential $`G\mathcal{T}G`$, so $`G`$ alone does not determine it
(`groundRel_compatible`, which holds by construction); and for invertible $`G`$ the relation
together with the potential does, $`\mathcal{T}=G^{-1}VG^{-1}=4CVC`$ (`inverse_inertia_unique`,
`inverse_inertia_formula`).

*[Cite]* OPH records that configuration equilibrium alone does not determine inertia, and that its
native counterpair shares the full reference with distinct finite kinetics
([`model/CONSTRUCTION.md`](../../model/CONSTRUCTION.md)).

So the vacuum together with the potential fixes the inertia within the class of pure Gaussian ground
states of a quadratic Hamiltonian with the given potential, and outside that class a fixed reference
does not determine the kinetics. The scalar benchmark these results attach to is OPH's S1, the
minimal canonical inertial completion of weighted scalar agreement
([`model/CANDIDATES.md`](../../model/CANDIDATES.md)).

## 3. Inertia: exclusion, equivalence, dynamical discrimination and a shared layer parameter

OPH's closure condition for inertia asks to derive the kinetic identification from admissible
operations or to constrain the dynamical response independently, and to exclude the same-equilibrium
alternative inertias ([`model/PREMISES.md`](../../model/PREMISES.md), RE-06).

**Proposition 1** (Conditional candidate for RE-06). Let the spatial stiffness $`K`$ be realized by
agreement seams of nonnegative conductance and let the slab index carry at least one time seam per
node. Assume (p1) positive capacities $`M`$; (p2) for every small enough $`\Delta\tau>0`$ the
history precision $`\Delta\tau\,K\otimes1+\Delta\tau^{-1}J\otimes L`$, with $`L`$ the Laplacian of
the time path, is realized by agreement seams of nonnegative conductance; (p3) every time seam
repairs at the rate $`2/\Delta\tau`$. Then the time coupling is $`J=\mathrm{diag}\,M`$. Under (p2)
alone $`J`$ is node-local, and when $`V=K`$ couples two distinct nodes the same-equilibrium family
$`J_\varepsilon=M(M+\varepsilon V)^{-1}M`$, with $`M+\varepsilon V`$ invertible, meets (p2) only for
$`\varepsilon=0`$. Independently of (p2) and (p3), for a positive definite coupling $`J`$ with
$`J^{-1/2}VJ^{-1/2}`$ positive definite, the two-time correlations of the stationary chain return
$`J`$ and the reconstructed acceleration matrix tends to $`J^{-1}V`$ as $`\Delta\tau\to0^+`$; on
every member of the family with $`V`$ and $`M+\varepsilon V`$ positive definite, negative
$`\varepsilon`$ included, this limit is $`M^{-1}V+\varepsilon M^{-1}VM^{-1}V`$, which is
$`A+\varepsilon A^2`$ in whitened coordinates; on OPH's three-site control that domain is exactly
$`\varepsilon>-1/4`$, and there the distance-two entry of the limit is $`\varepsilon`$. Reading the
common clock, mass and speed calibration as OPH's control script does, as a fit of the same-capacity
responses $`\alpha I+\beta A`$ with real $`\alpha`$ and $`\beta`$, the family response
$`A+\varepsilon A^2`$ is such a response exactly when $`\varepsilon=0`$ or $`A^2`$ lies in the span
of $`I`$ and $`A`$; this leaves only $`\varepsilon=0`$ on OPH's three-site control and at every
level $`q\ge3`$ of the golden family of Section 4, while on OPH's two-site action every family
response has this form, with $`\beta=1+41\varepsilon/16`$, positive exactly when
$`\varepsilon>-16/41`$.

*[Lean]* The first three statements are `eventually_repair_clock_iff` (for an arbitrary coupling
matrix $`J`$), `eventually_repair_isDiag` and `counterfamily_eventually_repair_iff`; at every
$`\Delta\tau>0`$ the two-time reader returns every positive definite coupling $`J`$ with
$`J^{-1/2}VJ^{-1/2}`$ positive definite (`history_reader_posDef`), on the family with $`V`$ positive
definite and $`\varepsilon\ge0`$ it returns $`J_\varepsilon`$ (`counterfamily_history_reader`), and
when $`V`$ couples two distinct nodes `counterfamily_reader_isDiag_iff` separates $`\varepsilon=0`$
from $`\varepsilon>0`$. For a positive definite coupling $`J`$ with $`A_J=J^{-1/2}VJ^{-1/2}`$
positive definite and $`\Delta\tau>0`$, the reconstructed coefficients, functions of $`A_J`$,
satisfy $`\mathcal{K}_{\Delta\tau}\mathcal{V}_{\Delta\tau}=\Omega_{\Delta\tau}^2`$
(`reconstructed_acceleration_eq`), and $`J^{-1/2}\,\Omega_{\Delta\tau}(A_J)^2J^{1/2}\to J^{-1}V`$
(`reconstructed_acceleration_tendsto`). On the family,
$`J_\varepsilon^{-1}V=M^{-1}V+\varepsilon M^{-1}VM^{-1}V`$ (`counterfamily_acceleration`), which is
$`A+\varepsilon A^2`$ in whitened coordinates (`counterfamily_acceleration_whitened`), and for
positive capacities with $`V`$ and $`M+\varepsilon V`$ positive definite, negative $`\varepsilon`$
included, the inertia $`J_\varepsilon`$ is positive definite, the reader returns it and the
acceleration matrix tends to that value (`counterfamilyInertia_posDef_of_shift`,
`counterfamily_dynamical_response_shift`; the case $`\varepsilon\ge0`$ is
`counterfamily_dynamical_response`). On the path $`A=I+L_{\mathrm{path3}}`$ with unit capacities,
the matrix of OPH's kinetic control
([`code/programme_obstruction_audit/kinetic_control.py`](../../code/programme_obstruction_audit/kinetic_control.py)),
the distance-two entry of $`\gamma(A+\varepsilon A^2)`$ is $`\gamma\varepsilon`$ for every
$`\varepsilon`$ (`path3_stiffness_response`), and the kinetic matrix $`I+\varepsilon A`$ is positive
definite exactly when $`\varepsilon>-1/4`$ (`path3_shift_posDef_iff`); of the cases
$`\varepsilon=0,\pm1/16`$ and the singular case $`-1/4`$ that the control script lists
([`code/programme_obstruction_audit/kinetic_control.py`](../../code/programme_obstruction_audit/kinetic_control.py)),
this domain contains the first three and excludes the last, so for each of the first three, with
$`A`$ positive definite (`path3A_posDef`), the acceleration matrix of the family tends to
$`A+\varepsilon A^2`$ (`counterfamily_dynamical_response_shift` with $`M=1`$), whose distance-two
entry is $`\varepsilon`$. For $`\varepsilon\ge0`$ the reader and that entry are also
`path3_dynamical_response`. OPH records $`Y_{13}=\gamma\varepsilon`$ for
$`Y=\gamma(A+\varepsilon A^2)`$ with actual acceleration $`-Yq`$
([`model/PREMISES.md`](../../model/PREMISES.md), RE-06); our statement is its case $`\gamma=1`$. For
every square matrix $`A`$ and every $`\varepsilon`$, some $`\alpha I+\beta A`$ equals
$`A+\varepsilon A^2`$ exactly when $`\varepsilon=0`$ or $`A^2=aI+bA`$ for some $`a,b`$
(`calibrated_response_iff`). On the path no $`\alpha I+\beta A`$ equals $`A+\varepsilon A^2`$ unless
$`\varepsilon=0`$, since the distance-two entry of $`\alpha I+\beta A`$ is $`0`$
(`path3_calibrated_separation`). At every golden level $`q\ge3`$ (Section 4) and every mass, the
whitened packet $`A_q`$ has entry $`0`$ and $`A_q^2`$ has entry $`2B_{01}^2`$ between the diagonal
neighbours $`(0,0,0)`$ and $`(1,1,0)`$, where the entry $`B_{01}`$ of the line operator is nonzero,
so again only $`\varepsilon=0`$ remains (`golden_distance_two_entry`, `lineOp_zero_one_ne`,
`golden_calibrated_separation`). On the two-site carrier of Section 5, with OPH's spring
$`\kappa=9/32`$ ([`model/PREMISES.md`](../../model/PREMISES.md), RE-05, and
[`model/observer/source_scalar/field.py`](../../model/observer/source_scalar/field.py)),
$`\gamma(V+\varepsilon V^2)=-\tfrac{25\gamma\varepsilon}{16}I+\gamma\big(1+\tfrac{41\varepsilon}{16}\big)V`$
for all $`\varepsilon`$ and $`\gamma`$ (`twoSite_calibrated_equiv`).

*[Cite]* OPH states that three-mode curvature and an exact distance-two stiffness-response
coefficient distinguish the counterfamily after common clock, mass and speed calibration
([`model/PREMISES.md`](../../model/PREMISES.md), RE-06).

*[Std]* For symmetric $`A`$, $`A^2`$ lies in the span of $`I`$ and $`A`$ exactly when $`A`$ has at
most two distinct eigenvalues (its minimal polynomial then has degree at most two), so in this
reading the calibrated response distinguishes the family from same-capacity dynamics exactly when
the whitened stiffness $`M^{-1/2}VM^{-1/2}`$ has at least three distinct eigenvalues, and this
criterion gives the carriers on which OPH's statement holds; the two-site stiffness has the two
eigenvalues $`1`$ and $`25/16`$. When the whitened stiffness has at most two distinct eigenvalues,
as on the two-site action, the calibrated acceleration response leaves $`\varepsilon`$ open, while
the two-time reader separates $`\varepsilon=0`$ from $`\varepsilon>0`$
(`counterfamily_reader_isDiag_iff` above) and $`\varepsilon\neq0`$ is excluded by (p2) whenever
$`V`$ couples two distinct nodes. If the calibration factors must be positive, the two-site members
with $`-16/25<\varepsilon\le-16/41`$, where $`\beta\le0`$, are distinguished as well.

The acceleration is that of the reconstructed Hamiltonian in model time (the Heisenberg equations of
$`\tfrac12(p^\top\mathcal{K}p+q^\top\mathcal{V}q)`$ give $`\ddot q=-\mathcal{K}\mathcal{V}q`$, a
standard computation outside Lean). The two-time correlations are those of the stationary chain; for
the finite open register Section 5 proves only that its equal-time blocks at interior slabs tend to
those of the chain as the window grows, while the first slab keeps the limit $`C_0(1+R)`$.

| Input | Status |
|---|---|
| Capacities $`M`$ | calibration |
| Correlations divided by the supplied positive action scale $`\hbar`$ before applying the reader | supplied normalization; not inferred from the reader |
| Small-step admissibility (p2), a matrix-level condition of this work | candidate hypothesis |
| Time-seam rate $`2/\Delta\tau`$ (p3) | shared layer parameter; calibration |
| Time conductances, the common rate scale | target-exposed choice |
| Common clock, mass and speed calibration, read as OPH's control script fits it (responses $`\alpha I+\beta A`$, real $`\alpha,\beta`$) | calibration (our reading of OPH's calibration sentence) |
| $`J=M`$, node locality, exclusion of $`\varepsilon\neq0`$, the limit $`J^{-1}V`$ of the reconstructed acceleration matrix, the calibration criterion | derived (Lean; equations of motion and the eigenvalue form [Std]) |

## 4. Continuum: the golden refinement family

OPH's scalar packet uses, at its levels $`q=F_n`$, the sorted points $`\{b\varphi\}`$, $`0\le b<q`$,
followed by the Dirichlet marker $`1`$, lumped weights $`W`$, Dirichlet stiffness $`K`$ and the
whitened operator $`A_q=I+B\otimes I\otimes I+I\otimes B\otimes I+I\otimes I\otimes B`$ with
$`B=W^{-1/2}KW^{-1/2}`$; its canonical quantization is the Fock space of $`\omega_q=A_q^{1/2}`$, and
its continuum proposition gives, along the golden refinements, convergence of the stated smeared
responses and Weyl-effect expectations to those of the massive Dirichlet free scalar on the cube,
with continuum operator $`\omega=(I-\Delta_D)^{1/2}`$
([`paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex`](../../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex)).
Its scalar-response adapter accepts the levels $`q=13,34,89,233`$
([`model/PREMISES.md`](../../model/PREMISES.md), RE-11, and
[`code/source_scalar_packet/source_common_scalar.py`](../../code/source_scalar_packet/source_common_scalar.py)).
The theorems below hold at every level $`q\ge2`$, which contains these.

**Proposition 2** (Conditional candidate for the S1 continuum connection). At every level $`q`$ of
the golden family the register theorems apply, and as $`\Delta\tau\to0^+`$ the register's
reconstructed frequency $`\Omega_{\Delta\tau}`$, its coefficients, its vacuum covariance and the
kernel $`\Omega_{\Delta\tau}^{-1}\sin(\theta\,\Omega_{\Delta\tau})`$ at every model-time lag
$`\theta`$ tend to $`\omega_q`$, to $`(I,A_q)`$, to $`\tfrac12\omega_q^{-1}`$ and to
$`\omega_q^{-1}\sin(\theta\,\omega_q)`$. In capacity-whitened coordinates the vacuum-covariance
error is at most $`\Delta\tau/4`$ times the squared length of the test vector, uniformly in $`q`$
and in the mode. Composed with OPH's continuum proposition, unchanged, this gives the iterated limit
$`q\to\infty`$ after $`\Delta\tau\to0`$, and for the vacuum quadratic forms of the smearings that
proposition admits, sampled as in Section 4, the limit is joint along every sequence
$`\Delta\tau_q\to0`$.

*[Lean]* For every $`q\ge2`$ the points $`\{b\varphi\}`$, $`b<q`$, are distinct (irrationality of
$`\varphi`$, the same statement as OPH's `orbit_injective`
([`Lean/Geometry/GoldenSourceAssignment.lean`](../../Lean/Geometry/GoldenSourceAssignment.lean)),
restated here because this library imports only Mathlib) and the grid is strictly increasing from
$`0`$ to the marker (`goldenFract_injective`, `goldenGrid_strictMono`, `goldenGrid_zero`,
`goldenGrid_last`); for every mass $`m^2>0`$ the packet operator and its whitened form are positive
definite at every level, the premise of the register theorems (`golden_packet_posDef`,
`golden_whitened_posDef`); the whitened packet is
$`B\otimes1\otimes1+1\otimes B\otimes1+1\otimes1\otimes B+m^2\,1`$ (`whitenedPacket_kroneckerSum`);
the $`q=13`$ carrier of the operator-level modules is this grid (`golden13_eq_goldenGrid`), and
$`q=34`$ is an instance (`golden34_whitened_posDef`). With slope $`3/2`$ two points coincide
(`rational_slope_collision`). At every level and for every mass $`m^2>0`$, as $`\Delta\tau\to0^+`$,
the limits of the proposition hold, and so does the Gaussian characteristic function
$`\exp(-\tfrac12 g^\top C_0(\Delta\tau)g)\to\exp(-\tfrac12 g^\top\tfrac12\omega_q^{-1}g)`$ for every
$`g`$ (`golden_continuum_time_limit`, from `latticeFreq_cfc_tendsto`, `kineticCoeff_cfc_tendsto`,
`potentialCoeff_cfc_tendsto`, `equalTimeCov_tendsto`, `vacuum_charFun_tendsto`,
`response_kernel_tendsto`). In one mode,
$`0\le\tfrac1{2\sqrt\lambda}-c_0(\Delta\tau,\lambda)\le\Delta\tau/4`$ for every $`\lambda>0`$ and
$`\Delta\tau>0`$ (`equalTimeVar_err_uniform`); the best uniform constant is
$`\varphi^{-5/2}/4\approx0.075`$, attained at $`\lambda=4\varphi/\Delta\tau^2`$
(`equalTimeVar_err_sharp`, `equalTimeVar_err_sharp_eq`). Hence
$`|v^\top(C_0(\Delta\tau)-\tfrac12A^{-1/2})v|\le\tfrac{\Delta\tau}{4}v^\top v`$ for every positive
definite $`A`$ and $`\Delta\tau>0`$ (`equalTimeCov_err_L2`), alongside
$`\le\tfrac{\Delta\tau^2}{16}v^\top A^{1/2}v`$ (`equalTimeCov_quadratic_err`). For a family of
positive definite $`A_k`$ on index sets that change with $`k`$ and vectors $`g_k`$ with
$`g_k^\top g_k`$ bounded and $`g_k^\top\tfrac12A_k^{-1/2}g_k\to c`$, every sequence
$`\Delta\tau_k\to0^+`$ gives $`g_k^\top C_0(\Delta\tau_k;A_k)g_k\to c`$ (`joint_limit_L2`; with
$`g_k^\top A_k^{1/2}g_k`$ bounded instead, `joint_limit`), and along the golden family with
$`m^2>0`$ this holds for every sequence of levels $`q_k\ge2`$ (`golden_joint_limit_L2`); without
such a bound the conclusion fails on one mode (`joint_limit_bound_needed`).

*[Cite]* OPH's continuum proposition admits smearings $`f`$ that are Riemann integrable with Riemann
integrable products, which includes finitely many slab-boundary jumps and so its compact slab
functions; for them the coherent-state expectations of finitely many bounded smeared Weyl effects
converge along the golden refinements to those of the massive Dirichlet scalar on the cube, and its
proof controls the squared norms of the Weyl parameters and compares sampled squared norms with
$`L^2`$ norms by product Riemann quadrature
([`paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex`](../../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex),
its continuum proposition and proof). For the iterated limit, the outer limit is that proposition
with its own argument (a fixed sine approximation first, then its error to zero); we use it
unchanged.

*[Lean, with cited hypotheses]* Let $`g_q`$ be the samples of such an $`f`$ weighted by the square
roots of the tensor masses. The cited proof gives
$`g_q^\top\omega_q^{-1}g_q\to\|\omega^{-1/2}f\|^2`$ and the quadrature gives
$`g_q^\top g_q\to\|f\|^2`$; these are the two hypotheses of `golden_joint_limit_L2`, so the
register's vacuum characteristic functions $`\exp(-\tfrac12 g_q^\top C_0(\Delta\tau_q;A_q)g_q)`$
tend to $`\exp(-\tfrac14\|\omega^{-1/2}f\|^2)`$, the continuum vacuum value, along every sequence
$`\Delta\tau_q\to0^+`$. The inner limit of the iterated limit is `golden_continuum_time_limit`.

| Input | Status |
|---|---|
| Golden coordinates, weights, stiffness, $`A_q`$ | supplied (OPH's declared inputs); positivity and Kronecker form checked in Lean |
| Mass $`m^2=1`$ | supplied (OPH's dimensionless mass one) |
| Register over each level with $`J=M`$ | candidate hypothesis (Section 3) |
| Smearings in the class of OPH's continuum proposition, sampled as in this section | target-exposed choice |
| Outer limit $`q\to\infty`$ | supplied (OPH's continuum proposition, cited unchanged) |

## 5. The reference state and the vacuum identity

OPH records that its scalar-response preparation uses a vacuum and controls supplied rather than
source selected, that a separate two-site action prepares its coupled vacuum to a stated trace
distance, and that preparation is not autonomous vacuum selection
([`model/PREMISES.md`](../../model/PREMISES.md), RE-09).

**Proposition 3** (Conditional candidate for RE-09). For the declared law (1) with $`\hbar>0`$ on
the register of a carrier with $`c\ge0`$, $`g>0`$, $`J=M>0`$ and $`T\ge1`$: in sampling time the
mean of every initial law with finite second moments tends to $`0`$ and its covariance to
$`\hbar P_T^{-1}`$; as the window grows, the equal-time block of $`\hbar P_T^{-1}`$ at interior
slabs tends to $`\hbar M^{-1/2}C_0(\Delta\tau)M^{-1/2}`$, while at the first slab it tends to
$`\hbar M^{-1/2}C_0(\Delta\tau)(1+R)M^{-1/2}`$, with $`R=r(A)`$ the transfer matrix of the chain; as
$`\Delta\tau\to0^+`$ the interior value tends to $`\hbar\,C_{\mathrm{vac}}`$, where
$`C_{\mathrm{vac}}=\tfrac12M^{-1/2}A^{-1/2}M^{-1/2}`$ is the position covariance at unit $`\hbar`$
of the ground state of $`\tfrac12\Pi^\top M^{-1}\Pi+\tfrac12f^\top Vf`$. At every level of the
golden family this limit is, in capacity-whitened coordinates, $`\hbar`$ times the position
covariance $`\tfrac12\omega_q^{-1}`$ of the vacuum of $`d\Gamma(\omega_q)`$, the vacuum that OPH's
scalar adapter supplies at its levels $`q=13,34,89,233`$, and on OPH's two-site action with
$`\kappa=9/32`$ it is $`\hbar\begin{pmatrix}9/20&1/20\\1/20&9/20\end{pmatrix}`$, $`\hbar`$ times the
position covariance at unit $`\hbar`$ of the coupled vacuum its preparation targets.

*[Lean]* For positive capacities and positive definite $`P`$, with
$`\tilde A=\mathbb{M}^{-1/2}P\mathbb{M}^{-1/2}`$ and
$`E(s)=\mathbb{M}^{-1/2}e^{-s\tilde A}\mathbb{M}^{1/2}`$: $`E(0)=1`$, $`E(s+s')=E(s)E(s')`$,
entrywise $`\tfrac{\mathrm{d}}{\mathrm{d}s}E(s)=-\mathbb{M}^{-1}P\,E(s)`$, $`E(s)\to0`$ and
$`E(s)m_0\to0`$ (`flow_zero`, `flow_add`, `flow_hasDerivAt`, `flow_tendsto_zero`, `mean_tendsto`);
for every initial matrix $`\Sigma_0`$,
$`\Sigma(s)=E(s)(\Sigma_0-\hbar P^{-1})E(s)^\top+\hbar P^{-1}`$ has $`\Sigma(0)=\Sigma_0`$, solves
$`\dot\Sigma=-\mathbb{M}^{-1}P\,\Sigma-\Sigma\,(\mathbb{M}^{-1}P)^\top+2\hbar\,\mathbb{M}^{-1}`$
entrywise and tends to $`\hbar P^{-1}`$ (`covFlow_zero`, `covFlow_hasDerivAt`, `covFlow_tendsto`).
One mode of the open register precision is inverted explicitly by the stationary kernel
$`c_0r^{|u-v|}`$ between slabs $`u`$ and $`v`$ plus the two boundary solutions $`r^u`$ and
$`r^{T-u}`$ (`greenMode_inv`); its diagonal exceeds $`c_0`$ by at most
$`c_0r\,(r^{2t}+r^{2(T-t)})/(1-r)`$ (`greenMode_diag_bound`); the register inverse is the
corresponding block matrix (`register_inv`, `declared_register_inv`); in whitened coordinates its
slab-$`t_k`$ block tends to $`C_0`$ whenever $`t_k\to\infty`$ and $`T_k-t_k\to\infty`$
(`interior_block_tendsto`) and its first-slab block tends to $`C_0(1+R)`$, which differs from
$`C_0`$ on a nonempty carrier (`end_block_tendsto`, `end_block_ne`); for the declared register these
limits are $`M^{-1/2}C_0M^{-1/2}`$ and $`M^{-1/2}C_0(1+R)M^{-1/2}`$
(`declared_interior_block_tendsto`, `declared_end_block_tendsto`), and the first-slab block of
$`\hbar P_T^{-1}`$ tends to $`\hbar M^{-1/2}C_0(1+R)M^{-1/2}`$, which for $`\hbar\neq0`$ differs
from the interior value (`first_slab_reference_limit`). The three limits and their values are
`iterated_reference_limit`, whose positivity hypothesis holds on every carrier with $`c\ge0`$,
$`g>0`$ and $`M>0`$ (`carrier_whitened_posDef`), with the vacuum covariance and its ground relation
$`GM^{-1}G=V`$ for $`G=(2C_{\mathrm{vac}})^{-1}`$ from `origCov_tendsto`,
`vacuumCov_mul_vacuumPrecision` and `vacuumPrecision_groundRel`; this $`G`$ is the only positive
definite solution of the relation (`vacuumPrecision_unique`). At every golden level the final value
in whitened coordinates is $`\tfrac{\hbar}{2}\omega_q^{-1}`$ (`golden_reference_limit`); the
two-site carrier (capacities $`1`$, one seam of conductance $`9/32`$, grounds $`1`$) has stiffness
$`\begin{pmatrix}41/32&-9/32\\-9/32&41/32\end{pmatrix}`$ and the stated limit (`twoSite_stiffness`,
`twoSite_vacuum_cov`, `twoSite_reference_limit`).

*[Std]* The moment equations of a linear SDE have unique solutions and a Gaussian initial law stays
Gaussian, so the law of the declared dynamics from a Gaussian initial law converges to the centred
Gaussian with covariance $`\hbar P_T^{-1}`$; for a general initial law the characteristic function
of the Ornstein–Uhlenbeck transition gives weak convergence [R1]. For positive definite $`G`$ the
Gaussian $`\psi(q)\propto\exp(-q^\top Gq/2\hbar)`$ solves the eigenvalue equation of
$`\tfrac12\Pi^\top M^{-1}\Pi+\tfrac12q^\top Vq`$ exactly when $`GM^{-1}G=V`$; it is then the ground
state (it has no node), its position covariance is $`\hbar(2G)^{-1}`$, and among centred pure
Gaussian wavefunctions with no position-dependent phase it is fixed by that covariance, up to
an overall phase. The ground state of the reconstructed Hamiltonian at slab
thickness $`\Delta\tau`$ has, per whitened mode, fidelity $`2\sqrt{\gamma}/(1+\gamma)`$ with the
canonical one, $`\gamma=\sqrt{1+\Delta\tau^2\lambda/4}`$ (a one-mode Gaussian overlap).

*[Cite]* For OPH's two-site action, directed rational bounds certify the initial trace distance of
the prepared state from that action's coupled vacuum
([`model/PREMISES.md`](../../model/PREMISES.md), RE-09).

*[Std]* By the identity above, that bound is a bound on the distance to the $`\Delta\tau\to0^+`$
limit of the vacua reconstructed from the declared dynamics, and at finite $`\Delta\tau`$ the
triangle inequality adds the trace distance
$`\big(1-\prod_\lambda 2\sqrt{\gamma_\lambda}/(1+\gamma_\lambda)\big)^{1/2}`$ between the
reconstructed and the canonical vacuum.

| Input | Status |
|---|---|
| Declared law (1) (a changed operation, not OPH's microscopic repair) | candidate hypothesis |
| $`\hbar`$ | supplied (OPH's absolute action scale) |
| Slab data (time conductances, $`\Delta\tau`$) | target-exposed choice |
| Two-site action ($`\kappa=9/32`$, capacities and grounds $`1`$), golden levels | supplied (OPH's declared inputs); limits checked in Lean |
| Convergence of moments, interior slab limit, vacuum covariance | derived (Lean; law-level convergence and ground-state identification [Std]) |

## 6. Non-claims and prior work

Nothing here selects a physical vacuum or OPH's common-source Gibbs reference. The sampling time and
the slab label are not identified with physical time; the premises (p2) and (p3), $`M`$, $`m`$,
$`\hbar`$ and $`\Delta\tau`$ are not derived; no source program family with records and service is
supplied, and neither interactions nor gravity are treated. The declared law differs from OPH's
microscopic repair in its noise and stationary law; for uniform capacities, seam rates $`1/N`$ and
no grounds the two agree at the level of the mean generator. A limit in a different order is not
asserted. OPH notes that vanishing kinetic perturbations can share the position-Weyl limit
([`model/PREMISES.md`](../../model/PREMISES.md), RE-11); Proposition 2 selects no inertia. The
remaining parts are standard: continuity of the functional calculus (in general form Mathlib's
`tendsto_cfc_fun`), uniqueness of positive square roots (in general form Mathlib's
`CFC.sqrt_unique`), convergence of Ornstein–Uhlenbeck moments (an abstract Lean treatment of
Ornstein–Uhlenbeck semigroups and the Mehler formula is in [R2]), Gaussian overlaps and iterated
limits. As far as our literature check found, what is specific to this work is the machine-checked
form of the *[Lean]* statements above for this declared law and its instances on OPH's carriers; we
do not claim new mathematics.

## References

[R1] H. Risken, *The Fokker–Planck Equation*, 2nd ed., Springer Series in Synergetics 18, Springer, Berlin, 1989.  
[R2] M. R. Douglas, Lean repositories `markov-semigroups` and `reflection-positivity`, <https://github.com/mrdouglasny> (accessed 11 October 2026).

## Appendix: Sentence-to-anchor table

| Statement | Anchor |
|---|---|
| Drift and mean repair generator; Lyapunov solution $`\hbar P^{-1}`$; transfer eigenvalues and equal-time covariance; bond reflection positivity; per-mode coefficient limits; what determines the inertia | [Lean, operator level] `repairGenerator_apply`, `register_drift`, `lyapunov_stationary`, `lyapunov_unique`, `register_stationary_unique`, `transferEig_add_inv`, `equalTimeVar_closed`, `chain_covariance_row`, `contact_identity_mode`, `transferMat_add_inv`, `equalTimeCov_mul_one_sub`, `chain_bond_rp`, `block_bond_rp`, `latticeFreq_tendsto`, `kineticCoeff_tendsto`, `potentialCoeff_tendsto`, `config_cov_indep`, `groundRel_compatible` (both by construction), `inverse_inertia_unique`, `inverse_inertia_formula` |
| $`J=M`$ under (p1)–(p3); node locality; $`\varepsilon=0`$ only | [Lean, operator level] `eventually_repair_clock_iff`, `eventually_repair_isDiag`, `counterfamily_eventually_repair_iff` |
| Reader returns $`J`$; acceleration $`\to J^{-1}V`$; family $`A+\varepsilon A^2`$ for $`V\succ0`$ and $`M+\varepsilon V\succ0`$; on the path $`\varepsilon>-1/4`$ exactly, distance-two entry $`\varepsilon`$ | [Lean] `history_reader_posDef`, `reconstructed_acceleration_tendsto`, `counterfamily_dynamical_response_shift`, `counterfamily_acceleration_whitened`, `path3_shift_posDef_iff`, `path3_stiffness_response` |
| Fitted responses $`\alpha I+\beta A`$ (real $`\alpha,\beta`$) reproduce the family exactly when $`\varepsilon=0`$ or $`A^2\in\mathrm{span}(I,A)`$; separation on the path and at golden $`q\ge3`$; two-site form with $`\beta=1+41\varepsilon/16`$ | [Lean] `calibrated_response_iff`, `path3_calibrated_separation`, `golden_calibrated_separation`, `twoSite_calibrated_equiv`; [Std] eigenvalue count |
| Golden family: grid, positivity, Kronecker sum, $`q=13`$, $`q=34`$ | [Lean] `goldenGrid_strictMono`, `golden_whitened_posDef`, `whitenedPacket_kroneckerSum`, `golden13_eq_goldenGrid`, `golden34_whitened_posDef` |
| The limits at each level; uniform error $`\Delta\tau/4`$ and its best constant; joint limit for bounded squared length; bound needed | [Lean] `golden_continuum_time_limit`, `equalTimeVar_err_uniform`, `equalTimeVar_err_sharp`, `equalTimeVar_err_sharp_eq`, `equalTimeCov_err_L2`, `joint_limit_L2`, `golden_joint_limit_L2`, `joint_limit_bound_needed` |
| OPH's smearings meet the joint-limit hypotheses; outer limit | [Cite] ([`paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex`](../../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex)) with [Lean, with cited hypotheses] `golden_joint_limit_L2` |
| Moments converge; Lyapunov flow | [Lean] `mean_tendsto`, `covFlow_hasDerivAt`, `covFlow_tendsto` |
| Interior slabs lose the end effect; the first slab keeps it | [Lean] `declared_interior_block_tendsto`, `first_slab_reference_limit` |
| Iterated limit on every carrier; unique positive solution of the ground relation; golden and two-site values | [Lean] `iterated_reference_limit`, `carrier_whitened_posDef`, `vacuumPrecision_unique`, `golden_reference_limit`, `twoSite_reference_limit` |
| State determined by covariance; finite-$`\Delta\tau`$ distance; transported bound | [Std] |
