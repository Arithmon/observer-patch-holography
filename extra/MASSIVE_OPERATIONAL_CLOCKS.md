# Massive flights, resolved reads and observable clocks

## Question and result

The timed tetrahedral walk has isotropic field speed `v=c/3`, while its
complete record cone has speed `c`. This distinction by itself does not
decide whether clocks and measurements at fixed spatial resolution share
the field cone. The construction below decides that question for an explicit
free process: a local mass coupling gives Dirac dynamics; smooth local CAR
operations converge to the same cone uniformly over states; a two-mass
interference clock exhibits the corresponding inertial time dilation.
Exact cell-resolving records retain speed `c`. Consequently the two cones
have different operational domains, with a quantitative resolution bound
between them. No change of the native time coordinate is made.

The parent results are the charged-flight construction in
`code/m1_quantum_transport`, the read-interface distinctions in
`code/m1_interfaces`, and the full-residual comparison in
`paper/tex_fragments/EFFECTIVE_QUANTUM_COMPARISON.tex`. The supplied graph,
quantum operations, native speed, masses and apparatus are explicit inputs.
This is not a source-selection theorem for A1--A3.

Throughout, hbar=1, c>0, m>0, and a>0. Matrix norms are operator norms;
function norms are L2 norms. Fourier transforms are unitary. The spatial
lattice is `a Z^3`; finite packet experiments use its periodic quotient.

## 1. One charged flight with mass

Let the four sign vectors be (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1), and

```
             [0  1  1  1]
C = 1/sqrt(3)[1  0  i -i].
             [1 -i  0  i]
             [1  i -i  0]
```

Then C*=C and C^2=I. Set `F_a(k)=diag(exp(-i a k.s))`, `W=F_a C`,
`tau=sqrt(3)a/c`, `R=I_2 tensor C`, and `beta=sigma_x tensor I_4`.
The eight-channel massive walk is

```
U_m(k) = exp(-i tau m beta) diag(W(k), W(k)*).
```

Here * is the adjoint, not entrywise conjugation. Execute, in order:

1. the onsite coin `diag(C,I)`;
2. parallel flights `diag(F,F*)`, with displacements `a s` and `-a s`;
3. the onsite coin `diag(I,C)`;
4. the onsite mass gate `cos(tau m) I - i sin(tau m) beta`.

Each moving register travels distance sqrt(3)a in time tau. Both flavours
fly concurrently. Ideal instantaneous onsite gates have the same convention
as the parent flight grammar; a physical nonzero gate overhead must be
charged and changes the speed. This is not a bounded-control-strength
refinement claim. Second quantization lifts each onsite unitary and register
permutation to number-conserving fermionic operations. The packet experiment
executes their one-particle sector, not a many-body hardware device.

### Entire spectrum

Write `mu=m tau`, assume `0<mu<pi/2`, and define

```
sin(epsilon)^2 = (sin(a k_x)^2+sin(a k_y)^2+sin(a k_z)^2)/3,
0 <= epsilon <= pi/2,
omega = acos(cos(mu) cos(epsilon)).
```

The four eigenvalues of W are `exp(+-i epsilon)` and
`-exp(+-i epsilon)`. In an eigenvector of W with eigenvalue exp(i theta),
the flavour block of U is `exp(-i mu sigma_x) exp(i theta sigma_z)`.
Its determinant is one and half its trace is `cos(mu) cos(theta)`.
It follows that the complete eight eigenvalues are

```
exp(+-i omega), exp(+-i(pi-omega)), each twice.
```

Multiplicity is retained at crossings. There are eight low-mass valleys
`a k in {0,pi}^3`; none is discarded. Each has four low-band Dirac modes
and four high-band modes. The positive excitation magnitudes are
`omega/tau` (four times) and `(pi-omega)/tau` (four times), with lower
bounds m and pi/(2 tau), respectively. For this construction define the
positive accounting observable `E_a=|i log U|/tau`, with the principal
logarithm, and its number-conserving Fock lift `dGamma(E_a)`. Its unique
vacuum is the empty Fock state, which is fixed by the actual walk. All clock
preparations use this same reference. This explicitly declared positive
observable commutes with the walk but is not its signed logarithmic
generator; neither it nor a filled negative-energy sea is being identified
with a derived physical Hamiltonian. These are declared Floquet excitation
magnitudes, not a thermodynamic energy law from A1--A3. In particular squaring
U does not turn the high band into low-energy matter.

For distance delta to `{0,pi}^3` in angular coordinates, the parent bound
`epsilon >= 2 delta/(pi sqrt(3))` and `omega>=epsilon` give a global
coercive bound. On an even cubic torus of fixed side L, dominated convergence
therefore gives the full free excitation partition function

```
log Z_a -> 32 sum_{n in Z^3} log(1+exp(-b sqrt(m^2+v^2 |2 pi n/L|^2))).
```

Indeed each of eight valleys contributes four low excitations; for each
fixed offset omega/tau tends to the displayed massive energy. The global
coercive bound dominates their tails by a summable exponential uniformly
in a. The high-band contribution is at most `4 (L/a)^3 exp(-b pi/(2 tau))`
and tends to zero. This proves a limit of the full spectrum, not a truncation
of finite spectra. The temperature parameter b is positive and fixed.

### Sharp velocity of every band and ballistic propagation

There is a stronger distinction than low versus high energy. Every smooth
band of this walk has group speed at most `v cos(mu)`, and this supremum
is sharp. Set `x_i=sin(p_i)^2` and `s=(x_1+x_2+x_3)/3`. Direct
differentiation gives

```
|grad_p epsilon|^2 = [3s-sum_i x_i^2]/[9s(1-s)] <= 1/3,
|d omega/d epsilon| = cos(mu) sin(epsilon)/sin(omega) <= cos(mu).
```

The first inequality is the variance inequality `sum x_i^2>=3s^2`.
Multiplying by a/tau gives the claimed speed. All four eigenphase choices
have the same speed magnitude. Near each middle point
`p_i in {pi/2,3pi/2}`, omega tends to pi/2 and its directional cone slope
tends to `cos(mu)/sqrt(3)`, attaining the supremum. The middle-band
crossings are retained; they have diverging excitation energy under
refinement, and are not extra low-mass species.

At a fixed lattice spacing, this is also a statement about transport of
arbitrary fixed normalized one-particle preparations, rather than just a
group-velocity interpretation. The distribution of `X/(n tau)` under U^n
converges weakly to the spectral velocity distribution, supported in the
ball of radius `v cos(mu)`. To prove this, write the characteristic function
in momentum space: it contains
`U(k)^(-n) U(k-xi/(n tau))^n` and the translated initial wavefunction.
Away from the finitely many band-crossing points, spectral projectors and
phases are smooth. Terms with different eigenvalues have projector product
tending to zero; terms in the same eigenspace tend to the phase
`exp(i xi.grad(theta)/tau)` with the appropriate sign convention.
The integrand is bounded by an L1 function after approximating the initial
L2 wavefunction by a smooth one; unitarity bounds the approximation error.
Dominated convergence and the characteristic-function continuity theorem
give the result. The crossing set has measure zero. In particular the
probability outside any ball of radius `(v cos(mu)+d)n tau`, d>0,
tends to zero. The statement fixes the preparation and lattice first; it
is not uniform over n-dependent encoders or a joint refinement/time limit.

A perfect one-flight record at c therefore does not imply a free ballistic
channel at c. The resolution theorem below is a separate fixed-duration
refinement statement and does not interchange these limits.

## 2. Full-channel dynamics at native time

Put `B(k)=diag(k.s)` and

```
X(k) = sqrt(3)/2 (B(k)+C B(k) C),
H_m(k) = m beta + v sigma_z tensor X(k),  v=c/3.
```

The matrices X_i commute with C and satisfy `{X_i,X_j}=2 delta_ij I`.
Thus the alpha matrices anticommute with beta and
`H_m(k)^2=(m^2+v^2 |k|^2) I_8`. The full eight-channel field, including
both signs and both R eigenspaces, has Dirac current norm at most v times
its density and propagation speed v. Integration of the continuity equation
over an expanding ball proves finite propagation for compact initial data;
the mass term cancels from the current identity. A positive-energy projection
alone would not have this property and is not made here.

Here is an explicit operator error, including leakage between R eigenspaces.
For `|k|<=K`, put `b=a sqrt(3)K`, `mu=m tau`. At `t=n tau`,

```
||U_m(k)^n - R^n exp(-i t H_m(k))||
 <= D_n(a,K,m)
 := floor(n/2)(b^2+4 mu b) + (n mod 2)(b+mu b).       (1)
```

For a proof set `A=diag(B,-CBC)` and `A'=RAR=diag(CBC,-B)`.
Then `U=exp(-i mu beta) exp(-i a A) R` and

```
U^2 = exp(-i mu beta) exp(-i a A)
      exp(-i mu beta) exp(-i a A').
```

The sum of these four Hermitian exponents is `2 tau H_m`. The unitary
product/exponential difference is bounded by half the sum of pairwise
commutator norms, proved by Duhamel integration and induction. Four
mass/flight pairs contribute at most mu b each, the flight pair at most
b^2, and the mass/mass pair zero. Telescope complete pairs. For an odd
tick, compare A with `(A+A')/2`, at cost at most b, and apply the two-factor
bound to the remaining mass and averaged flight, at cost mu b. R commutes
with H_m. This proves (1), without a compressed-band residual or a hidden
renormalization of time. At K=0 the error is exactly zero, including odd n.

At valley `pi z/a`, U gains the scalar `eta=(-1)^(z_x+z_y+z_z)`; use
`(eta R)^n` in (1). Every valley obeys the same full-channel bound.
For `n tau<=T`, a useful upper bound is

```
D <= a T (sqrt(3)c K^2/2 + 2 sqrt(3)m K)
     + sqrt(3)a K + 3 a^2 m K/c.                     (2)
```

## 3. Actual local reads, uniformly over the state

A resolution assumption belongs on the operations, not on an unobserved
choice of successful states. Let f be an eight-component H2 function and
`S_a f(x)=a^(3/2) f(a x)` its actual lattice samples. Define

```
M2(f) = || |k|^2 fhat(k) ||,
A_a(f) = (a/pi)^2 sqrt(4 pi^2 + pi^4/45) M2(f).
```

Let J_a be the Fourier isometry after restriction to the Brillouin cube
`[-pi/a,pi/a)^3`. Poisson summation followed by Cauchy--Schwarz proves

```
||S_a f-J_a f|| <= A_a(f).                           (3)
```

To see the constant explicitly, a nonzero alias l has
`|k+2 pi l/a| >= (pi/a)||l||_infinity`. The weighted square sum is
`sum_{r>=1}(24 r^2+2)/r^4 = 4 pi^2+pi^4/45`. Integration of the other
Cauchy--Schwarz factor is bounded by M2(f)^2 because alias cubes partition
the exterior. Density extends the identity from Schwartz functions to H2;
point sampling is defined since H2 embeds in continuous functions in 3D.

For `K<pi/(2a)`, equations (1) and (3) give

```
||U^n S_a f - R^n S_a exp(-it H_m) f||
 <= Q_a(f) := 2 A_a(f) + D_n ||f|| + 2 M2(f)/K^2.   (4)
```

The last term can be replaced by twice the actual Fourier tail norm beyond
K. The evolution preserves M2. The estimate is on full vectors: the other
bands and valleys remain in the tail estimate, rather than being removed
from the actual state or evolution. Modulation by any of the eight valley
carriers gives the corresponding version with the envelope's M2.
For fixed f,m,T, choosing `K=a^(-1/4)` after fixing units yields
`Q_a=O(sqrt(a))`. This is a resolution limit with a stated error, not an
assertion of exact microscopic relativistic causality.

Taking adjoints in the low-momentum operator estimate gives (4) for inverse
evolution as well, with U^(-n), R^(-n) and exp(+it H_m). This is the version
used for a future Heisenberg read.

Write a(h) for the CAR annihilation operator, whose norm is ||h||. For
nonzero sampled f,g, use normalized kernels f_a,g_a and number effects
`N_f=a(f_a)* a(f_a)`, `N_g=a(g_a)* a(g_a)`. Assume f and g are compactly
supported and their supports are separated by more than vt. The continuum
solution of f has zero overlap with g, as does its sampled version; R is
onsite. Equation (4) implies

```
|<g_a,U^(-n) f_a>| <= Q_a(f)/||S_a f||,
||[N_g,N_{U^(-n) f}]|| <= 2 Q_a(f)/||S_a f||.        (5)
```

The commutator formula follows directly from CAR and Cauchy--Schwarz.
Encoding a bit by applying `exp(i theta N_g)` or the identity changes the
later N_f probability, in *any* Fock state and with any spectator ancilla,
by at most `min(1,2 |theta| Q_a(f)/||S_a f||)`. Duhamel's formula for
conjugation proves the bound. Finite products of such encoders obey the
sum of their bounds. This gives actual compactly supported preparations
and reads, not a nonlocal hard momentum-cutoff detector. It covers this
specified CAR operation class; arbitrary local encoders are not asserted
to have fixed resolution merely because one detector is smooth.

### What the bound does not erase

Before the precoin, put a single particle `C|s>` at one site in the first
flavour. After that coin it is |s>, after one flight it is exactly at a s,
and both remaining onsite gates leave its total occupation there equal to
one. An empty input gives zero there. This is a perfect speed-c record in
the massive walk too. A delta-site kernel does not form a bounded-H2
envelope under refinement, so (5) supplies no vanishing error for it.
Likewise a bound on final matter energy does not bound transient
cell-resolving encodings or resettable controller resources. No such
replacement of the operation hypothesis is used.

## 4. A readable normalizable clock

Take two copies with masses `m1<m2` and coherent internal labels 1,2.
These labels are apparatus degrees of freedom, not a new source axiom.
For a specified constant velocity u with |u|<v, let

```
gamma = 1/sqrt(1-|u|^2/v^2),
k_j = gamma m_j u/v^2,  E_j = gamma m_j.
```

These momenta are the solutions of the group-velocity equation
`grad E_m(k)=v^2 k/E_m(k)=u` for the dispersion derived in section 2.
They are not assigned by a proper-time clock rule.

There is a common unit spinor xi with `H_mj(k_j)xi=E_j xi` and `R xi=xi`,
since the normalized matrices H/E are identical for both masses and commute
with R. The common positive eigenspace in the R=+1 sector is nonempty.
Use the normalized
three-dimensional envelope

```
phi_sigma(x) = (2 pi sigma^2)^(-3/4) exp(-|x|^2/(4 sigma^2)),
Psi_j(0,x) = phi_sigma(x) exp(i k_j.x) xi.
```

No positive-energy projection of this preparation is needed. Its small
negative-energy part is included in the following bound. The actual
continuum solution differs in norm from the rigid packet

```
Phi_j(t,x) = phi_sigma(x-u t) exp(i k_j.x-i E_j t) xi
```

by at most

```
e_j(T) = sqrt(3) v/(m_j sigma)
         + sqrt(15) T v^2/(8 m_j sigma^2).           (6)
```

Proof: the positive projector `P(k)=(I+H(k)/E(k))/2` is Lipschitz with
constant at most v/m_j. The negative component contributes at most
`2 v sqrt(<|q|^2>)/m_j`. The scalar energy Hessian has norm at most
`v^2/m_j`; Taylor's remainder contributes
`T v^2 sqrt(<|q|^4>)/(2m_j)`. For this Gaussian the two moments are
`3/(4 sigma^2)` and `15/(16 sigma^4)`. This proves (6). It applies at all
times in [0,T], not just at selected phases. For the superposition
`(Psi_1, Psi_2)/sqrt(2)`, use `e=max(e_1,e_2)`.

Use a positive detector effect centered at u t:

```
F_theta(t,x) = exp(-|x-u t|^2/(2R^2)) |+_theta><+_theta| tensor I_8,
|+_theta> = (|1>+exp(i theta)|2>)/sqrt(2).
```

For the rigid packets its actual probability is

```
p_theta(t) = [p_R + V_R cos(Delta m t/gamma + theta)]/2,
p_R = (1+sigma^2/R^2)^(-3/2),
V_R = p_R exp(-s^2 |Delta k|^2/2),
s^2 = sigma^2 R^2/(sigma^2+R^2).                    (7)
```

Here V_R is the reference fringe coefficient; the conventional normalized
visibility of this reference is V_R/p_R. Neither number alone bounds the
contrast of the finite clock after its error terms are included.

This sign convention follows `<+_theta|Psi>`; changing the detector phase
convention changes the sign of theta only. Gaussian integration gives (7)
because `Delta E - Delta k.u = Delta m/gamma`. A single unobservable
global phase is not being called a clock. The relative phase is read by a
bounded positive internal interference effect. Truncating the detector to
`|x-u t|<=L_d` gives finite extent and
changes every probability by at most `exp(-L_d^2/(2R^2))`. The actual
continuum probability differs from (7) by at most e plus this truncation
error. The rest reference has frequency Delta m; the moving reference has
frequency Delta m/gamma, using the same native time t and speed v. Actual
packet reads obey the stated error relative to those laws.

### Finite lattice error and a constructive parameter choice

Initialize the lattice with the normalized samples of the two Gaussian
packets and apply the two actual walks for n ticks. Sample the same detector
at the lattice sites. The odd-tick R carrier commutes with this detector and
acts identically on the mass labels, so it cannot create a clock phase.
For a Gaussian centered at k_j,

```
M2_j^2 = |k_j|^4 + 5 |k_j|^2/(2 sigma^2) + 15/(16 sigma^4).
```

Choose `K>max|k_j|`, less than pi/(2a), and let r=K-max|k_j|. A valid
common Fourier tail norm is

```
h = sqrt(erfc(sqrt(2) sigma r)
         + 2 sqrt(2) sigma r exp(-2 sigma^2 r^2)/sqrt(pi)).
```

Let `A=max A_a(Psi_j)`, `D=max D_n(a,K,m_j)`. Put
`B(z)=sum_{l in Z^3, l!=0} exp(-z |l|^2)`. This sum can be bounded without
an infinite computation by

```
B(z) <= (1+2 exp(-z)/(1-exp(-3z)))^3 - 1,  z>0,
```

because `l^2 >= 1+3(l-1)` for positive integer l. Use a cancellation-safe
form of the cubic when evaluating a tiny bound. Define
`b_sigma=B(2 pi^2 sigma^2/a^2)` and `b_s=B(pi^2 s^2/(2a^2))`.
If `a |Delta k|<=pi` and `b_sigma<=1/2`, the lattice probability differs
from (7) by at most

```
4(e + D + 4A + 2h) + 4(b_s+b_sigma)
    + exp(-L_d^2/(2R^2)).                           (8)
```

For detail, (4) compares the sampled exact walk to the sampled continuum
solution at cost `D+2A+2h`. Sampling the continuum/rigid difference adds
at most `e+2A`. Poisson summation bounds both sampled Gaussian norms
squared within `1+-b_sigma`. Passing from vector error to an effect
expectation after normalization costs at most
`(1+sqrt(3))/sqrt(1/2)<4`. The sampled rigid Gaussian detector integral
has alias error at most b_s: for a nonzero reciprocal vector,
`|2 pi l/a-Delta k|>=pi |l|/a`. Dividing by the initial norm and comparing
to the integral adds at most `2(b_s+b_sigma)`, covered by (8).
Detector truncation costs its sup-norm tail. This proves (8) for the
specified local measurement; no assumed detector agreement is inserted.

The preparation can also have finite extent. Multiply each sampled packet
by a fixed C2 cutoff equal to one on the cube `[-P,P]^3`, tapering to zero
within the larger cube `[-P-sigma,P+sigma]^3`, and normalize. For example
use the product of the one-dimensional tapers
`1-10s^3+15s^4-6s^5` in the transition interval 0<=s<=1. This retains
R=+1 and uniformly finite Fourier moments. The lost Gaussian probability
is at most `3 erfc(a floor(P/a)/(sqrt(2)sigma))`: bound each discrete
one-dimensional tail by its integral starting at floor(P/a), and use the
union bound; the unshifted Gaussian lattice normalization is at least its
integral by Poisson summation. The normalized preparation changes by at
most twice the square root of this tail in norm. Unitarity and the effect
bound therefore add at most

```
e_prep = 2 sqrt(3 erfc((P-a)/(sqrt(2)sigma)))         (P>a)            (8a)
```

to (8), at every time, with no postselection on successful evolution.

On Fock space the bounded click effect is `I-Gamma(I-F_theta)`, where
Gamma denotes fermionic second quantization; on the one-particle sector
it restricts to F_theta. After the spatial truncation it acts only on the
finite detector region and has spectrum in [0,1]. It can be read by onsite
mass-label rotations, local occupation detection and the stated acceptance
weight. The probability is that of those local clicks at their detection
time. Collecting distributed click records can add latency; no instantaneous
global report or zero-duration communication across the detector is claimed.

For example choose c=3, v=1, m1=100, m2=100.1, u=(0.6,0,0), sigma=10,
R=100, L_d=600, P=80, a=1e-9, K=max|k_j|+1 and
`T=2 pi gamma/Delta m+tau` (through the first native tick after one moving
cycle). Equations (6)--(8a) give an outward-enclosed probability bound below
0.035153, while the reference V_R exceeds 0.7457. All
parameters are finite. The very small lattice spacing is an analytic
existence witness, not a claimed affordable lattice run. The accompanying
finite runs use coarser spacings, retain all Fourier modes and report their
actual probabilities and errors rather than claiming (8) is sharp there.
The independent checker uses outward interval arithmetic for these strict
inequalities, with `erfc(x)<=exp(-x^2)/(x sqrt(pi))`; it does not treat a
rounded point evaluation as a certified enclosure.

There is an actual observable contrast guarantee as well. With theta=0,
compare time zero and the native tick nearest `pi gamma/Delta m`. Its
reference phase misses pi by at most `d=Delta m tau/(2 gamma)`. If B is
the full probability bound, the finite clock's probability swing is at
least `V_R(1-d^2/4)-2B`, using `cos(d)>=1-d^2/2`. Outward interval
arithmetic gives a swing **greater than 0.6753** for the example, including
both finite cutoffs. This is a claim about actual click probabilities.
The finite packet signal is compared with a sinusoid with explicit error;
it is not claimed to be an exactly monochromatic clock at finite resources.

More generally any fixed number of readable cycles and error tolerance
can be achieved with finite parameters: first choose R/sigma and
L_d/R large; choose sigma Delta m small for visibility; increase m1 to
make (6) small over the chosen cycles; choose r to reduce h; finally
decrease a to reduce D,A and the aliases. This is an explicit construction,
not a premise that good clocks exist. It concerns inertial free clocks.

### Finite accounting energy and the localized record cost

The clock construction does not require diverging energy in its declared
positive observable. For `mu<=pi/3`, any normalized lattice packet whose
spinor lies in R=+1 has the full-spectrum estimate

```
<E_a> <= m + (pi c/6)<|k|> + 4 sqrt(3) pi c a <|k|^2>.               (9)
```

Indeed `||Re U-cos(mu) R|| <= a sqrt(3)|k|`. If P_hi is the negative
spectral subspace of Re U, then
`||P_hi psi|| <= a sqrt(3)|k| ||psi||/cos(mu)` for `R psi=psi`:
on that subspace `|Re U-cos(mu)|>=cos(mu)`. The high-band weight is
therefore at most `12 a^2 |k|^2`. Also `omega<=mu+epsilon` and
`epsilon<=pi a|k|/(2 sqrt(3))`, since `asin(x)<=pi x/2` for 0<=x<=1.
Bound the high energy by pi/tau and integrate to obtain (9). At isolated
middle crossings either equal-energy subspace convention gives the same
energy; the pointwise estimate can be interpreted almost everywhere.
The sampled fixed Gaussian has uniformly bounded Fourier moments as a
tends to zero, by its exponentially decaying aliases. Thus (9) is uniform
under refinement. The actual walk conserves E_a exactly. This argument
keeps the high-band leakage and charges its full cost.

In contrast any normalized single-site particle, whatever its internal
spinor, has

```
<E_a> >= asin(1/sqrt(6))/(2 tau).                                  (10)
```

Its Fourier density is uniform. On half the Brillouin cube,
`sin(p_x)^2>=1/2`, hence every excitation magnitude is at least
`epsilon/tau>=asin(1/sqrt(6))/tau`. Integration proves (10). This
divergence applies in particular to the perfect speed-c record. It does
not prove that a bound on final matter energy restricts arbitrary control
protocols. It identifies the cost of the explicit counterexample in the
same observable used for the packets. The conversion of E_a to physical
energy remains a separate identification.

## 5. Consequence for the M1 interface

The same process supplies an isotropic massive dispersion, a quantum
interference read, and a native-duration clock law at v=c/3. Equations
(4)--(8) state which operations share that cone and quantify their errors.
An exact requirement equating every microscopic record with this field
cone excludes this process, even though its fixed-resolution matter and
clocks converge together. A requirement restricted to the stated resolved
read algebra admits it with controlled error. No claim is made that the
restricted algebra exhausts all physically realizable OPH observations.

The result answers a clock and interface question left by the speed-budget
obstruction. It does not supply the A1--A3 admission of this quantum gate
family, the physical selection of its mass parameters or an interacting
matter theory. Those are existing source/physics boundaries, not additional
premises used to infer the clock or resolved-read formulas.

## Reproduction and intellectual attribution

The companion `code/m1_operational_clocks` contains independent numerical
replay, full-mode periodic packets, explicit wrong-time/high-band/phase
controls and a strict receipt verifier. The inequalities and all-level
limits above are analytic proofs. Finite numerical catalogs do not prove
their universal quantifiers, and these operator theorems are not claimed
as Lean formalizations.

Weyl-to-Dirac automaton constructions and relativistic dispersion are
standard; see D'Ariano and Perinotti, [Quantum cellular automata and free
quantum field theory](https://arxiv.org/abs/1608.02004). Quantum-clock
interference and time dilation are standard subjects; see Chiba and
Kinoshita, [Quantum Clocks, Gravitational Time Dilation, and Quantum
Interference](https://doi.org/10.1103/PhysRevD.106.124035). The use of a
restricted operational norm, instead of unconstrained channel norm, has
precedents including Winter, [Energy-constrained diamond norm with
applications to the uniform continuity of continuous variable channel
capacities](https://arxiv.org/abs/1712.10267). Here the restriction is on
spatially resolved CAR operations, not an energy-constrained state set.
The OPH-specific contribution is the charged tetrahedral extension, its
full-channel errors and spectrum, and a common constructive read/clock
interface with an explicit surviving faster record.
