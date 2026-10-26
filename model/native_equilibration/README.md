# Preparing the stationary native sector, with an error bound

The [quantum mixing theorem](../observer/QUANTUM_MIXING.md)
applies to any finite connected weighted graph, not just one elementary
observer. This packet attaches that result to the actual screen federations.
It replaces **exact stationary preparation** with **discarded native repair
attempts plus a certified finite error**, within each supplied binary
particle-number sector. It does not select that number or an initial mixture
of different sectors.

The architecture supplies the graph; the implemented C or E rule supplies
the Kraus maps. The schedule is independent uniform seam selection and a
fresh fair coin. Forgotten preparation records are not used for conditioning.
These are the hypotheses of this construction, and must not be silently
replaced by the production observer's adaptive schedule.

## Mixing theorem and its concrete attachment

For `V` ports, `M` seams and fixed total `k`, write `d=binomial(V,k)`.
The source's matrix-unit difference blocks establish self-adjoint positive
channels. Their diagonal block is symmetric exclusion. The theorem of
[Caputo, Liggett and Richthammer, J. Amer. Math. Soc. 23 (2010) 831-851, §4.1.1, equation (4.1)](https://arxiv.org/abs/0906.1238)
identifies its gap with the underlying graph's random-walk gap. The existing
OPH quantum proof bounds every coherence block by the same gap. Consequently,

```text
tau = I/d,       gamma = lambda_1(graph Laplacian)/(2M),
distance((Phi^n tensor id_R)(rho), tau tensor rho_R)
    <= min(1, d/2 (1-gamma)^n).
```

This is uniform over input densities and arbitrary untouched references.
It concerns the unread channel, not the accumulated record tape. It also
bounds the discrepancy of **any common subsequent complete instrument**,
including adaptive interventions. The total variation of the resulting
complete record laws cannot exceed the state trace distance. Conditioning
on a rare successful herald can amplify that error; the same unconditional
allowance cannot simply be assigned to the conditional law.

The graph certificate verifies the retained federation graph and gives an exact,
conservative gap certificate without relying on a numerical eigenvalue.
If `D` bounds its diameter, the range/path and variance bounds give
`lambda_1 >= 4/(V D)`. A breadth-first traversal supplies
`D <= min(V-1,2 eccentricity(root))`. Thus `gamma >= 2/(M V D)`.
This bound is deliberately conservative; it is not evidence of rapid
equilibration or a gap uniform in federation size.

The integer preparation budget avoids numerical logarithms and huge rational
powers. For `b=ceil(1/gamma_lower)`, Bernoulli's inequality gives
`(1-gamma_lower)^b <= 1/2`. Enough such blocks make `d/2 * 2^-blocks`
smaller than the requested tolerance. The receipt gives budgets for total
two and half filling at all three retained screen levels. Their large counts
are sufficient bounds, not measured minimum preparation times.

## What executes

On the elementary carrier, the exact gap is `(5-sqrt(5))/60`; the rational
bound `11/240` follows from `sqrt(5)<9/4`. Start with a specified complex
superposition of two configurations, enumerate actual source Kraus records,
and apply their signed-difference blocks repeatedly for the certified attempt
count. Both C and E are executed. The code then measures the complete four-word
law of two successive absolute-current bits and compares it with the tracial
preparation's law. No Hamiltonian or fitted state is used in that propagation.
The conservative budget is **374 discarded attempts**, with a uniform trace
distance allowance **`33/131072 < 1/1024`**. The particular executed initial
state has numerical distance about `1.04e-8`; that diagnostic does not replace
the uniform bound. The large-federation budgets are computed certificates,
not trillions of executed native attempts.

The execution uses the exact source-enumerated native channel with floating
matrix arithmetic; its numerical residuals are diagnostics, separate from
the analytic error certificate. Independent small-graph tests enumerate all
native matrix units, check the exact spectral identification numerically,
and propagate an entangled untouched reference. No fresh Lean theorem is
claimed. The executable attachment supplies preparation/observation error bounds
using the stated general mixing theorem.

## What this resolves, and what it does not

**Within this connected fixed-total unread regime**, the law has a unique
attracting stationary state.
An experiment may prepare it approximately with a declared error budget.
Exact stationarity at a finite attempt count is not asserted. Disconnected
graphs, conserved totals, preparation records and adaptive schedules change
the conclusion and remain explicit.

This result does not select a nontrivial modular clock. The modular flow of
`I/d` on its full matrix algebra is the identity. A Palm reference conditioned
on an event, or a reference on a reconstructed transfer algebra, is a different
operational construction. Native events and changing records continue at
stationarity; a stationary density does not make the history static.
Nor does convergence select a physical vacuum,
the number of spatial dimensions, or a relativistic field action. All counts
here are **discarded native attempts**, not elapsed physical or modular time.

The full unconstrained tensor algebra also contains coherences between the
all-empty and all-full configurations that the native difference records
cannot distinguish. The fixed-sector hypothesis is therefore essential;
there is no claim of a unique state on that larger algebra.

## Continuing records and constrained entropy

[history_entropy.py](history_entropy.py) extends the two-moment reference
envelope to actual recorded population histories, including adaptive seams
crossing the core and bounded record-dependent stopping. Two native E
contrast records recover the protected core from a connected three-port
census; the C-current implementation uses a bounded waiting service with
state-independent success and every timeout retained. Neither readout
receives hidden occupation bits. The instrument counts compatible endpoints
and accounts explicitly for the entropy of its retained fair coins.

At fixed initial occupation moments, its cost-corrected record entropy has
an exact relative-entropy deficit and nonpositive Fisher curvature. The
finite reference-preparation sandwich transfers to the same adaptive
record readout. Cap predictions are compatible marginals of one conditional
endpoint distribution, not separately prepared states. Unresolved C
timeouts have a quantified counting excess; forgetting informative records
and freezing an input-dependent stopping cost are tested failures.

The current derivation, service costs and physical boundary are in
[the common-source geometry discussion](../common_source/SLOW_GEOMETRY.md#entropy-stationarity-from-continuing-native-records).
The population formula and the coherent formula below use the same finite
native instrument and are retained inside `variation_results.json`. Proper
volume, geometric area and physical stress remain separate interfaces.

### Coherent history stationarity and recovery

`history_entropy.quantum_response` closes the finite **coherent recorded
entropy-stationarity** interface on this same adaptive instrument. It verifies
every primitive endpoint, complete tape, fair-coin probability and record-only
control decision. For a branch `K_b=sqrt(w_b) U_b P_b`, with
`w_b=2^(-T_b)/beta`, it computes all 4,356 exact rational entries
`G_ij=sum_b w_b 1{i,j in domain_b}`. The reference-Petz decoder pulls the
complete quantum output back to `G` entrywise multiplied by the input. The
reference factors cancel for every faithful diagonal source reference; the
identity extends to arbitrary untouched references. This applies the
established [Petz sufficiency theorem](https://arxiv.org/abs/quant-ph/0209053),
and identifies its equality algebra for the actual native histories.

For the declared two-attempt adaptive prefix and census, both C and E have
recovery algebra `M_28 + 4 M_8 + 6 C`, of dimension 1,046. These are exactly
the initial-occupation fibres on the four used ports `0,1,5,8`; every
diagonal population and arbitrary entanglement inside each block is recovered.
The complete C instrument destroys 558 matrix units; conditioning on its
state-independent successful flag destroys 2,814. E destroys 3,278. The C
timeouts therefore carry genuine coherent output, even where a population
counting completion is ambiguous. The exact counts include every source
matrix unit, rather than sampling prepared states or fitting numerical ranks.

For the existing source-selected faithful `sigma(m,c)`, put
`F(rho)=S(Gamma rho)-E_rho[T] log(2)-log(beta)`. Direct branch calculation gives

```text
F(rho) = S(sigma)+h1 delta<m>+h2 delta<c>-D(Gamma rho || Gamma sigma),
dF_sigma[X] = 0                    if Tr X=Tr(q X)=Tr(c X)=0,
d2F_sigma[X,X] = -sum_ij G_ij |X_ij|^2 L(sigma_i,sigma_j),
L(a,b)=(log(a)-log(b))/(a-b),       L(a,a)=1/a.
```

Thus all admitted coherent tangents are stationary. Output-flat directions
have `G_ij=0`; information-lossless directions instead have `G_ij=1`.
Invisible coherences can share the maximum output entropy with `sigma`;
this maximum does not determine a unique source state.
The strict fixed-moment curvature dimensions are 3,795 for complete C,
1,539 for successful C, and 1,075 for E. The same reference is stationary
after each of the 1,024 core-containing cap completions; their visible and
recoverable algebra intersections are enumerated exactly. The existing finite
preparation sandwich also gives its same gradient bound for the full quantum
output, by positivity and operator monotonicity of logarithm. These are
variations in the declared dimensionless state parameter. The geometric
metric/volume variation, area coefficient, stress normalization and continuum
tail required by the Einstein branch remain separate premises.

### Joint vacuum and geometry test

The [finite slab calculation](vacuum_geometry.py) uses the actual twelve-port
angular operator, a three-site Dirichlet radial action and its Gaussian ground
state. It varies the action-derived spatial metric at fixed proper slab
volume, then computes the central radial cell's reduced entropy and both
boundary-face areas from that same vacuum. For its uniform variation,
the [receipt](vacuum_geometry_results.json) gives entropy derivative
-1.1004136186 and two-face-area derivative -4/3; these have the
same sign for every nontrivial uniform shape variation. No positive
bare-area coefficient makes their sum stationary on this slab.

The finite regional modular kernel is not port-local: two Gaussian momentum
perturbations with identical canonical stress changes at every port give
different entropy derivatives. A varying radial coordinate length cannot
remove the sign obstruction. Fixed native site count also fails to keep one
isotropic proper ultraviolet spacing through the nontrivial variations.
This tests the supplied Gaussian slab, not a renormalized small-ball
generalized entropy. The native compiler has not selected or prepared its
vacuum, and the experiment supplies neither physical stress nor curvature.

### Proper-ball modular and geometry comparison

[ball_modular.py](ball_modular.py) restricts the full three-dimensional
Dirichlet Gaussian vacuum covariance to a proper ball with equal lattice
spacing in three Cartesian directions. A frozen native chart coefficient
defines the constant comparison metric and proper-time coordinate; it is not
identified with the unknown native continuum tensor. The regional modular
quadratic form is reconstructed from the reduced covariance, with a precision
refusal for unresolved near-pure modes.

For momentum displacements, the exact entropy derivative of a Gaussian-noise
mixture equals the coherent relative-entropy charge. Nested regions of the
same vacuum satisfy the corresponding monotonicity check. The independently
fixed [conformal-ball stress weight](https://arxiv.org/abs/1102.0440) gives
`4*pi^2*R^4/15` for a constant momentum profile and `4*pi^2*R^4/105` for a
linear profile. At `R=1`, box side `16` and spacing `2/7`, all 179 regional
modes are retained; the calculated charges are `2.687795261` and
`0.397474951`, respectively 2.124% and 5.716% above those continuum values.
The [receipt](ball_modular_results.json) retains six lattice cuts, four box
sizes, an independent sine-spectrum computation and entropy secants. The
grid errors alternate and the finest cut refuses insufficient precision;
these controls do not establish a continuum convergence rate. Equal spacing
still leaves finite cubic rotational anisotropy.

An independent space-form metric calculation gives `G_tau,tau=3*k` and the
fixed-proper-volume coefficient `dA/dG_tau,tau=-4*pi*R^4/15`. This supplies
the geometric coefficient without defining curvature from entropy or stress.
Combining that identity with the ball first law still requires a common
native small-ball limit, a renormalized area density and selection of joint
state--geometry stationarity. Neither Newton's constant nor the Einstein
equation is obtained by fitting these finite comparisons.

```sh
python -m model.native_equilibration --check model/native_equilibration/results.json
python -m model.native_equilibration.vacuum_geometry --check model/native_equilibration/vacuum_geometry_results.json
python -m model.native_equilibration.ball_modular --check model/native_equilibration/ball_modular_results.json
python -m pytest -q -W error model/native_equilibration
```

### Fermion vacuum and regional algebra

[fermion_ball.py](../common_source/fermion_ball.py) uses the same calibrated
M15 fermion action as the spectrum experiment. A region-local canonical
particle-hole map gives 21 Dirac species (six quark flavors in three colors and three charged leptons) and three massless Weyl
species. One antiperiodic Fourier cutoff and one ball define the restricted
global vacuum projector, its modular kernel and entropy. Independent entropy
secants and nested relative entropy test the state variation. At the retained
fine cut, massless and running-top comparisons with the conformal stress
weight differ by -0.1953% and +22.83%. That weight is not assumed to be the
massive modular Hamiltonian; the grid errors are not a convergence theorem.

A separate 32-CAR source preparation executes 13,616 native pulses and checks
its covariance against the all-left BdG vacuum. Larger regional grids are
covariance calculations. The region contains the full free-matter CAR algebra;
interacting gauge constraints, dressing and gauge-field edge contributions
are not supplied. Calibration, spectral cutoff, vacuum selection and joint
geometric stationarity remain explicit premises.

### Vacuum-compatible active repair and its work cost

[vacuum_repair.py](vacuum_repair.py) specifies compatibility independently of
the update: a uniform native geometry bit and conditional low excess energy
above the field Hamiltonian's vacuum. One actual C/E seam event resamples
that bit. Its committed record selects a source preparation and ordinary
SWAP reset with fresh blanks, retaining the displaced input in the environment.
The geometry family has density $1/r$ and angular scale $r^3$, hence fixed
proper slab volume branch by branch on one fixed Fock algebra.

Eight source preparations cover both geometry branches and all four native
angular bands, at oscillator cutoff eight. Exact rational Rayleigh enclosures
and the full oscillator gap bound each branch's 36-oscillator product vacuum
error by $1/256$ in half trace norm, including the binary norm defect. The
large product is not materialized. A separate cutoff-two, one-band joint
reset executes all 16 geometry/field input columns; full-field reset costs
remain tensor-product interface requirements, with fresh blanks and entropy
export charged.

The ideal finite-ground-state family's regional entropy and bare two-face
area have derivatives of the same sign. The existence of this family of joint fixed states
does not imply bare generalized-entropy stationarity. Source fidelity at
discrete geometries does not certify an entropy derivative of the source
word family. A separate exact work test is stronger: for each common normal
mode, $a_\pm=\mu_\pm\omega_\pm$ gives stationary cycle work

```math
W=\frac{\omega_++\omega_-}{16}
  \frac{(a_+-a_-)^2}{a_+a_-}\ge0.
```

The constant native angular mode makes the total strictly positive for
distinct geometries. At branches $r=0.9,1.1$, the ideal total is about
0.163501 in the shared model-coordinate Hamiltonian units. Exact Rayleigh
intervals separately certify work greater than $1/1000$ for the retained
source preparations. A static geometry potential has zero mean change in
this stationary cycle and cannot pay that cost. Reciprocal backreaction or
an explicit work reservoir is needed. This supplies a finite preparation
and work ledger; it derives neither autonomous cooling nor a renormalized
small-ball Einstein balance.

### Reciprocal thermal matter and geometry

[reciprocal_vacuum.py](reciprocal_vacuum.py) gives a separate finite thermal
construction on the same cutoff-two scalar Hamiltonians at $r=9/10,11/10$.
Compatibility is the joint Gibbs state of those declared Hamiltonians at
the supplied inverse temperature $1/\sqrt{14}$. Auxiliary native proposals
select energy-level pairs. A bank of fixed-gap thermal qubits supplies
energy-conserving exchanges, with exact pairwise detailed balance and
all rejected proposals and bath outputs retained. The full quantum state
is stationary even though the conditional scalar Hamiltonians do not commute.

Matter free energy fixes geometry probabilities $p_g=Z_g/(Z_0+Z_1)$,
approximately $0.57509,0.42491$. A uniform-geometry conditional Gibbs state
has stationarity defect about $0.05419$ in population L1 distance. The same
law acts on the first native correlation/exchange fermion chain in its
fixed-occupation-one sector; charge-sector preparation is explicit.
Energy readout, native proposals and selected complete source collisions
are executed. The full reservoir-bank controller and its thermal refresh
history are not executed. Endpoint energy balance is distinct from
physical pulse-controller work. Temperature, Hamiltonian and thermal-input
services remain premises. Free-energy stationarity under changes of the
geometry probabilities does not supply a metric variation, a vacuum
selection theorem or the Einstein entropy--area condition.

```sh
python -m model.native_equilibration.reciprocal_vacuum --check model/native_equilibration/reciprocal_vacuum_results.json
```

### Zero-temperature free-field wedge

[wedge_modular.py](wedge_modular.py) uses the uncalibrated staggered action
derived in [the common source field](../common_source/native_staggered_field.py).
The negative-energy sea is the declared vacuum. Transverse Fourier reduction
preserves the spatial cut and gives four scalar staggered chains per reduced
transverse momentum, retaining two Dirac tastes. Independently placed sites
and bonds define the lapse coupling and boost energy before the modular law.
The exact half-chain theorem gives $K_{\rm half}=\gamma B$, with
$\gamma=4\kappa K(\kappa')/c$, $\kappa=(1+\mu^2)^{-1/2}$ and $\mu=Ma/c$.
Directed rational AGM bounds prove
$0\le1-\gamma/(2\pi/c)\le\mu^2/4$. This fixed-band continuum normalization
uses no temperature fit.
Here $H$, $B$ and $K_{\rm half}$ are one-particle kernels. A finite regional
density has $-\log\rho_A=d\Gamma(k_A)+\log Z_A\,I$, with its actual kernel
$k_A=\log[(I-C_A)/C_A]$; the finite kernel is not asserted to equal the
infinite-wedge weight. The infinite quadratic operator is defined relative
to the filled sea, equivalently by particle/hole excitation normal ordering.
Scalar normalization cancels in trace-zero first-law variations.

A Legendre sign-polynomial supplies a local covariance error bound
$\Lambda r^{k+1}/m$, $r=(\Lambda-m)/(\Lambda+m)$, when degree-$2k+1$ walks
cannot reach the external boundary. The stated gap, norm and walk hypotheses
give an explicit volume schedule. This is not a global state-distance bound.
The separate eight-CAR source preparation retains its raw 256 amplitudes;
an exact binary Rayleigh quotient and algebraic Fermi gap certify its global
vacuum fidelity. Full-Fock partial traces and CAR derivatives check the
entropy first law. Finite-interval disagreement with the infinite-wedge
weight and unresolved-purity refusals remain explicit.

The staggered entanglement-level series, with its own tail bound, gives
bare planar entropy after transverse factorization. The independently
geometric area is $N_yN_za^2$; its coefficient is regulator- and
species-dependent. Larger grids are ideal vacuum controls. The source
preparation, coherent/parity services and physical accuracy remain specified
requirements. The result establishes the free half-space modular/stress
attachment, not a massive small-ball identity, physical Newton constant or
joint curved-metric entropy stationarity. See the
[mass research doc](../../docs/research/OPH_PARTICLE_MASS_RESPONSE.md#515-zero-temperature-modular-response-of-the-same-free-field)
for the common field and the distinct interacting-pole test.

```sh
python -m model.native_equilibration.wedge_modular --check model/native_equilibration/wedge_modular_results.json
```
