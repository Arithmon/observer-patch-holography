# What the anomalous density must supply to predict lensing

## Result and relation to the OPH source

Issue #751 asks whether one collar-source model fixes dynamics and lensing
without separately fitting their gravitational responses. The answer for the
released density interface is **no**. The obstruction survives the Einstein
equation, covariant stress conservation, positive density, the dominant energy
condition, regular central geometry, finite total mass and identical exterior
geometry. It is not the earlier freedom to choose a non-Einstein spatial
response: here the same Einstein equation holds in every alternative.

We also obtain a positive result: an exact joint observable reconstructs the
missing source parameter in the flat-rotation annulus when the ray's winding
is known, with sharp bounds under specified stress restrictions. For a general
spherical system, the full rotation and unwrapped scattering profiles
determine the remaining metric function by an Abel inverse. Thus the missing information is a radial stress/density
closure, not two independently adjustable lensing and rotation laws.

These are completed mathematical results on the Einstein branch, not a new
microscopic collar model. The [dark-sector paper](../cosmology/oph_dark_matter_paper.tex)
and `DarkSector.lean` supply a scalar modular-charge interface and a
conditional linear enclosed-mass law. They do not select the relativistic
stress law or its normalization. We do not claim the alternatives below are
microscopic OPH realizations. Their purpose is to test exactly the proposed
macroscopic inference. Parent #751 remains open for physical attachment,
magnitude, localization and abundance; its bounded dynamics/lensing
underdetermination exit is answered here.

The distinction between density, pressure, rotation and lensing is established
GR, including [Faber and Visser (2006)](https://arxiv.org/abs/astro-ph/0512213)
and the exact geodesic treatment of
[Scharf and Bräunlich (2009)](https://arxiv.org/abs/0912.5058). Our contribution is the
explicit source-interface test: sharp finite-annulus bounds, a regular
finite-mass degeneracy with a rigorous nonzero ray separation, and independently
checked certificates tied to the OPH claim boundary. No novelty is claimed
for Einstein's spherical equations or Abel inversion themselves.

## 1. Conventions and exact source reconstruction

Use areal radius and geometrical units G=c=1:

    ds² = -A(r) dt² + B(r) dr² + r² dΩ²,
    A=exp(2ν), B=(1-2m/r)^(-1),
    T in a static orthonormal frame = diag(ρ,p_r,p_t,p_t).

The cosmological constant is zero in this isolated-system model.
There is no radial energy flux. Let u=rν'. Timelike circular geodesics have
local speed squared u, independently of B, for 0<u<1. Radial stability requires
(r²u/(1-u))'>0; in particular every constant-u orbit used below is stable.
The Einstein equations and conservation give

    4πr²ρ = m',
    4πr²p_r = (r-2m)ν' - m/r =: P,
    4πr²p_t = (rP' + (m'+P)u)/2.                         (1)

The last equation is the anisotropic TOV identity
p_r'=-(ρ+p_r)ν'+2(p_t-p_r)/r. Direct substitution in the angular
Einstein equation gives the same expression; the checker uses that angular
equation rather than assuming the conservation identity.

The explicit countermodels set the baryonic stress to zero and use geodesic
test probes without backreaction. Their length and strength are declared
model inputs, not estimates for a named galaxy. This isolates the scalar
dark-source inference; a comparison with a real disk must restore its
baryonic geometry and jointly model both channels. Each countermodel has
one metric and one conserved tensor for all probes.

Given A, **one function m remains**. Conversely, given A and physical ρ on
the full radial domain plus the central mass constant, m'=4πr²ρ fixes m,
then (1) fixes both pressures and B. One cannot hold both this complete
physical density and A fixed and still choose arbitrary lensing. A scalar
collar charge not yet identified with that density does not supply this
missing function. Nor does a fitted Newtonian dynamical mass automatically
equal the Einstein mass m.

In SI units m=G M/c², u=v_c²/c², ρ_geom=G ρ_mass/c² and
p_geom=G p_SI/c⁴.
Here ρ_mass denotes total mass-equivalent energy density, including kinetic
energy, rather than just rest-mass density. The coefficient κ gives
M(r)=κ c²r/G. If identified with the paper's anomalous deep profile,
κ=√(G M_b a_0)/c². That identification uses a source-supplied a_0;
the measured speed does not separately derive it. No H0 or Λ normalization
is introduced. The flat annulus is the ideal anomalous-dominated profile,
not a spherical replacement for a measured baryonic disk or a prediction
through the Newtonian transition.

## 2. Full flat-annulus classification

Fix 0<u<1, 0<r0<R, r0≤r≤R, and

    A(r)=A(R)(r/R)^(2u),  m(r)=κr,  0<κ<1/2.

Every such metric is a solution with

    4πr²ρ   = κ,
    4πr²p_r = u-(1+2u)κ,
    4πr²p_t = u²(1-2κ)/2.                               (2)

This class exhausts constant circular speed and exactly linear Einstein
enclosed mass on the annulus (the central integration constant is specified
to be zero in this local mass law). The pressures have the same r^-2
dependence and obey conservation exactly. The dominant energy condition
ρ≥|p_r|, ρ≥|p_t| is equivalent to

    u/[2(1+u)] ≤ κ < 1/2.                               (3)

Proof: the radial inequality gives the lower endpoint; its other side is
automatic. The tangential inequality gives κ≥u²/[2(1+u²)], a weaker
bound for 0<u<1. If both pressures must also be nonnegative, the sharp set is

    u/[2(1+u)] ≤ κ ≤ u/(1+2u).                          (4)

Both endpoints and all interior points occur. In particular, even positive
pressures and dominant energy do not select one density normalization from
the rotation curve. With the additional radial coldness bound
|p_r|≤ηρ, 0≤η<1, intersect (3) with

    u/(1+2u+η) ≤ κ ≤ u/(1+2u-η).                        (5)

This is a bound on radial pressure; no hidden tangential coldness is claimed.
At η=0, p_r=0 fixes κ=u/(1+2u) and p_t/ρ=u/2. Static nonzero dust
(both pressures zero) cannot support this profile. The commonly inserted
Newtonian identification κ=u instead has p_r/ρ=-2u and
p_t/ρ=u(1-2u)/2; it is a small-stress approximation at u≪1, not an
exact zero-pressure Einstein source. Under (4), the weak-field relation is
κ/u∈[1/2,1]+O(u), not a unique mass coefficient.

### What a cold particle source would add

There is a physically motivated restriction stronger than dominant energy.
For a nonnegative distribution of massive particles in the static frame,
T^(ab)=∫ f p^a p^b d³p/p^0. Hence the principal pressures are
nonnegative and

    p_r+2p_t = ∫ f p^0 |v|² d³p ≤ δρ                  (2a)

if every occupied speed obeys |v|²≤δ<1 (a bound on the corresponding
energy-weighted second moment suffices). Opposite directions remove net
flux without affecting this inequality. This is a standard kinetic-matter
restriction; it does not follow from the scalar modular charge or from
homogeneous comoving a^-3 dilution. An anomalous medium need not be a
massive-particle gas, so imposing (2a) is a substantive source test.

For (2), (2a) and nonnegative pressures are equivalent at the tensor level to

    u(1+u)/(1+2u+2u²+δ) ≤ κ ≤ u/(1+2u),   δ≥u.        (2b)

The interval is empty if δ<u. Its fractional mass width is exactly

    (κ_max-κ_min)/κ_max = (δ-u)/(1+2u+2u²+δ).           (2c)

Proof: the pressure trace scaled by 4πr² is
u+u²-(1+2u+2u²)κ. Its upper inequality gives the lower endpoint
in (2b); p_r≥0 gives the upper one. Subtracting them gives (2c).
These inequalities imply the dominant energy bounds. At δ=u they select
p_r=0 and p_t/ρ=u/2. For δ=C u with fixed C≥1, the relative mass
ambiguity is at most (C-1)u. The exact bending interval follows by (6),
using these same endpoints with no independent light-deflection parameter.

The certificate supplies rational outward bounds on the ratio
α(κ_min)/α(κ_max), obtained by squaring rational square-root brackets,
not by trusting a rounded decimal. The endpoint geometry cancels from this
ratio. At u=10^-6 and δ=2u the fractional bending spread is about
5×10^-7, while the nonnegative-pressure-only family allows about 25%
relative to its upper bending endpoint. Thus a demonstrated cold stress
law would almost eliminate this ambiguity in a weak field. Large slip is
not an automatic prediction of the anomalous-density proposal.

Equation (2b) is an exhaustive **tensor-level** restriction, not a proof that
every point is realized by a global stationary collisionless distribution.
The zero-radial endpoint has the standard
[Einstein-cluster interpretation](https://arxiv.org/abs/0705.1756) as an
isotropically oriented circular-orbit ensemble: each occupied particle has speed squared
u, radial pressure zero and two equal tangential stresses ρu/2. General
distribution-function existence and matter perturbation stability are not
needed for, and do not follow from, this necessary-bound theorem. Physical
OPH promotion requires deriving the stress restriction or directly deriving
the density; the theorem quantifies what that missing evidence would buy.

## 3. An exact finite-endpoint lensing observable and its inverse

A ray has turning radius r0 and both endpoints on the shell r=R. Its
impact parameter is b=r0/√A(r0). Define s=(r0/R)^(1-u) and
θ=acos(s). The acute angle with the radial line measured by either endpoint's
static orthonormal frame is Ψ=asin(s). Relative to the outward radial direction,
the incoming source angle is π-Ψ and the outgoing receiver angle is Ψ.
The unwrapped angular sweep along the ray is exactly

    Δφ = 2θ / [(1-u)√(1-2κ)].

Indeed dφ/dr = √B/[r√((r/r0)^(2-2u)-1)]; substitution
z=(r0/r)^(1-u) gives the integral. Define the finite-distance bending
angle by the endpoint-angle convention

    α = Δφ + 2Ψ - π
      = 2θ {1/[(1-u)√(1-2κ)] - 1}.                     (6)

This is the endpoint-angle convention of
[Ishihara et al.](https://arxiv.org/abs/1604.08308), with the angular sweep
retained along the path. This operational comparison specifies both endpoint
angles and the unwrapped central azimuthal separation; it is zero in the flat reference. No infinite
isothermal halo, asymptotic flatness of the annulus, exterior cutoff or
cosmological thin-lens distance is smuggled into (6). For astrophysical
images one must supply the exterior and observer/source geometry separately.

At fixed u,r0/R the angle increases strictly with κ. Substitution of the
endpoints of (4) or (5) gives exact sharp lensing intervals. Conversely,

    κ = (1 - [2θ/((1-u)(α+2θ))]^2)/2.                  (7)

Thus a joint exact measurement with known winding in this model fixes κ,
p_r/ρ and p_t/ρ; no second fit amplitude is available. Positions alone only
give angular separation modulo 2π. The winding qualification is essential:
at u=1/5 and r0/R=2^(-5/4), set

    κ_n = (1-25/[36(2n+1)²])/2,    n=0,1,2,... .         (7a)

These distinct dominant-energy sources all have Ψ=π/6 and the same
endpoint azimuthal difference π modulo 2π, while their unwrapped sweeps
are (2n+1)π and α=(2n+1/3)π. The receipt verifies n=0,1,2 exactly;
the displayed identity proves the complete family. Unknown winding
therefore precludes a unique inverse from those endpoint data. Conversely,
all retained nonnegative-pressure rays have u≤1/5, so
Δφ<π√(1+2u)/(1-u)≤π√35/4<2π. They have zero full winding.
The stronger global scattering inverse below likewise requires unwrapped
branch-resolved scattering data, not just an observed angle modulo 2π.

For r0/R→0 followed by u→0,
the nonnegative-pressure interval gives α/u→[3π/2,2π]. Its upper/lower
ratio tends to 4/3. This is a limit of the explicitly defined local
observable, not a claim about an isolated infinite flat-rotation galaxy.
The certificate includes finite R/r0=2,10,100 and u=10^-6,10^-2,1/5,
all endpoint/interior controls, and inversion checks. Reported decimal
values are reproducible numerical evaluations; the interval formulas and
endpoint inequalities, not floating quadrature, prove the bounds.

## 4. The degeneracy survives a regular centre and identical exterior

The annulus alone does not establish a globally admissible isolated system.
Here is a separate full-domain witness, not an invented cutoff of (2).
Fix 0<ε≤1/100 in units of an arbitrary length scale and set

    m0(r)=ε r³/(1+r³),
    u0(r)=m0/(r-2m0),
    ν(r)=-∫_r^∞ u0(s) ds/s.

This fixes A with A(∞)=1. Since r²≤1+r³, m0/r≤ε and u0≤1/98.
The centre is regular: m0=O(r³), ν'=O(r), density tends to
3ε/(4π). The metric is at least C² in regular central coordinates;
no analyticity at the centre is asserted. The total mass is ε, the lapse
integral converges at both ends, and there is no horizon. The baseline has
p_r=0, p_t=ρu0/2 and satisfies dominant energy everywhere.

Let

    h(r)=(r-1)^4(2-r)^4  for 1≤r≤2, and 0 otherwise,
    m_t(r)=m0(r)+t ε h(r),       |t|≤1/10,                (8)

with the **same A** for every t. The bump and its first three derivatives
match zero at both ends. All metrics have the same central neighbourhood,
the same lapse and entire circular-speed profile, the same exterior at r≥2,
and the same total mass. No thin shell is introduced. The source is obtained
from (1), so it is covariantly conserved and obeys every Einstein equation.
There is one static diagonal stress tensor in each model, not independently
tuned potentials for two observables.

Here is an all-radius dominant-energy proof, not a grid test. Outside [1,2]
it is the baseline. Inside, write R0=m0', δR=tεh',
δP=-tεh(1+2u0)/r and T=4πr²p_t. The elementary bounds are

    ε/27 ≤ R0 ≤ 3ε,  |h|≤1/256, |h'|≤1/16,
    |δR|≤ε/160, |δP|≤ε/2000, |δP'|≤9ε/1280,
    |T|/ε ≤ 3/196 + 9/1280 + (1/160+1/2000)/196
             < 1/27-1/160 ≤ (R0+δR)/ε.                (9)

For the derivative estimate use |u0'|≤8ε, obtained from
u0'=(rm0'-m0)/(r-2m0)^2, and
|(1+2u0)/r|≤21/20, |[(1+2u0)/r]'|≤6/5. Equation (9) also
dominates |δP|, so density is positive and both principal pressures satisfy
the dominant energy condition. The bound is deliberately conservative and
is checked as exact rational arithmetic. The pressures can include tension;
the separate annulus family already supplies a nonnegative-pressure witness.
Matter perturbation stability, a microscopic equation of state and a
distribution function are not inferred from energy conditions. These
examples meet the stated macroscopic interface, not unstated matter laws.

For a scattering ray with r0=1/2, the turning point is common and unique
because (r/√A)'=(1-u0)/√A>0. The deflection difference is supported only
on [1,2]:

    α_+ - α_- = 2∫_1^2 b√A/r² ·
      [√B_+ - √B_-]/√(1-Ab²/r²) dr > 0.                (10)

Here + and - mean t=±1/10. Since d√B/dm≥1/r,
b√A≥r0, and h≥81/65536 on [5/4,7/4],

    α_+ - α_- ≥ 81 ε / 3512320 radians > 0.             (11)

Every factor and the integration interval is retained in this lower bound.
There is no weak-field truncation, fitted curve, missing ray tail or noisy
numerical cancellation in the separation proof. At ε=1/100 the bound is
81/351232000 radians. Both the ADM mass and exterior metric are identical;
adding an exterior-mass measurement cannot remove this interior ambiguity.
The smoothness of (8), algebraic field equations and complete bounds are
independently replayed in the certificate.

A stronger two-sided enclosure uses the full positive bump integral:

    I = ∫_1^2 h(r)/r³ dr = 248 log(2)-1719/10,
    ε I/5 ≤ α_+ - α_- ≤ (539/400) ε I/5.                (11a)

For the upper bound, b√A/r0≤exp(log(4)/98)<49/48;
the turning-point denominator contributes at most 6/5 because
(49/48)²/4<1-25/36. With the compactness bound in (9), the
spatial derivative contributes at most 11/10, since
(1-(2+1/1280)/100)^3(11/10)^2>1. Their product is 539/400.
For the lower bound each corresponding factor is at least one. The
independent checker encloses log(2) using the positive atanh series with
argument 1/3 and its exact geometric tail; cancellation in I is therefore
controlled. The producer's 70-digit numerical logarithm is not trusted as
an interval. The resulting rational enclosures are supplied at ε=1/100
and ε=10^-6. This is a certified observable difference at finite strength,
including all ray tails (which cancel exactly), not merely a leading-order
coefficient or a floating-point sign test.

A separate independent numerical check integrates the actual lapse
ν(r)-ν(1/2) and both radial metrics in (10), with the difference of square
roots rationalized to avoid cancellation. It gives approximately
1.09172226528319×10^-6 and 1.06941244667967×10^-10 radians at the two
strengths, strictly inside their respective exact enclosures. Precision
refinement is tested. These numerical values are controls; the inequalities
(9)--(11a), not quadrature accuracy, certify the result.

## 5. What a joint measurement determines in the general system

For a regular asymptotically flat spherical geometry with 0≤u<1, define
the optical radius x=r/√A and

    Q(x)=√B(r)/(1-u(r)).

Assume differentiability and decay sufficient for the following Abel
integrals and their derivative to converge. This holds for (8). The exact
unwrapped scattering law and inverse are

    α(b)=2b∫_b^∞ [Q(x)-1]/[x√(x²-b²)] dx,
    Q(x)=1-x²/π ∫_x^∞ [d(α(b)/b)/db]/√(b²-x²) db.     (12)

To prove the inverse, set f=(Q-1)/x² and recognize α/b as
2∫_b^∞ f(x)x dx/√(x²-b²). The elementary Abel inversion follows by
interchanging the two convergent integrals; the inner beta integral is π.
The subtraction Q-1 removes the flat π sweep. Rotation supplies u and,
with normalization at infinity, A. Thus (12) yields B, then m and the
entire tensor (1). This is uniqueness from the **complete** scattering
profile, not a claim of stable reconstruction from noisy finite images.
The checker exercises forward and inverse integrals for four independent
power-law optical controls and a linear combination, without importing the
producer's gamma-function formula. These are transform controls, not
additional physical halo models.

The forward map smooths while the inverse differentiates α/b. Finite noisy
data require a likelihood and regularization or a source-selected profile;
the exact inverse supplies no statistical error bar. Missing large-impact
tails and the astronomical mass-sheet/geometry degeneracies cannot be set
to zero merely to use (12).

## 6. Consequence for the proposed physical comparison

The macroscopic Einstein interface alone permits the different outcomes
(6) and (10). A scalar anomaly bound, comoving a^-3 dilution, universal
coupling, or agreement with a fitted rotation curve does not select their
radial pressure relation. This is the specific missing selector. If OPH
supplies the full physical density independently, (1) and the central mass
condition remove it; if it supplies a pressure closure, the coupled system
must first be solved with boundary data. In either case the same amplitude
must carry into lensing.

This result does not exclude anomalous modular matter or establish a new
dark component. It prevents counting a pressure-blind rotation fit as a
joint dark-sector test. The constructive output (3)--(7) states precisely
which lensing ranges would follow from three distinct source classes:
dominant energy alone, nonnegative pressures, or a declared radial coldness
bound. Changing class after seeing a mismatch is a new candidate. The
additional positive-particle moment bound (2b) makes this criterion useful
for a cold-matter reference: its speed budget must come from the same
source/state and cannot be selected separately in the lensing channel. The
pressureless-radial branch has the same leading joint signal as ordinary
cold gravitating matter; that agreement by itself is not OPH-specific
evidence.

No new natural measurements are ingested. The existing SPARC sample is
already exposed and consists of disk rotation data; it supplies neither the
spherical null-scattering function nor this missing stress law. An eligible
comparison must fix a system/sample, baryon geometry and stellar population,
distances and inclinations, lens/source redshifts, external convergence and
shear, selection and joint covariance, then transport one source parameter
set unchanged across both channels. A baryons-plus-dark-matter reference
needs matched freedom and its own stress/profile assumptions. Until the
source closure is selected, a numerical joint likelihood would test an
analyst's extra halo assumption. The bounded exit is therefore an explicit
underdetermination witness and a quantitative joint-observable interface,
not a claimed ready prediction or a new prospective freeze.

## 7. Reproduction and adversarial checks

The [package](../code/dark_source_lensing/README.md) produces a compact
receipt. Exact symbolic checks compare independently reconstructed angular
Einstein components, conservation, orbit stability, the bump endpoints and
global rational bounds. Numerical controls compare closed-form ray sweeps
with a transformed geodesic integral and Abel inversion with direct
quadrature. They are high-precision controls, not interval-certified natural
observations. The positive separation (11) and the stress inequalities are
exact certificates.

The verifier rejects omitted cases, changed stresses, signs, units of the
reported normalized quantities, equality at a forbidden endpoint, broken
bump matching, wrong lower bounds, altered inversion, invalid numeric types,
duplicate fields, extra fields and false physical-promotion claims. Semantic
mutations bypass hash custody; the CLI is also exercised with assertions
disabled and producer imports blocked. Historical receipts and physical
premise classifications are not rewritten to turn this result into source
selection.

## 8. Maintainer-style audit findings

The inverse was audited against the entire stated parameter domain, not only
the weak-field numerical controls. This exposed an omitted qualification:
angles modulo a full turn do not determine the unwrapped ray sweep. Equation
(7a) now gives exact, dominant-energy counterexamples, and the theorem,
registry and receipt specify known winding. The retained weak-field catalogue
has no full winding by a separate bound. The degenerate endpoint r0=R is
excluded explicitly, and source/receiver angle orientations are fixed.

The global result was also checked directly at the metric level. The verifier
integrates the common lapse from the actual mass function and rationalizes
the difference between the two radial metrics. Both source strengths agree
with the exact enclosures, and independent precision refinement tests the
numerical control. A deliberately incompatible enclosure fails that control
even when the analytic-certificate checker is bypassed. Exact global energy
and bending inequalities remain the proof; the new quadrature supplies a
distinct implementation check.

The closure audit retains the source-identification boundary: dominant energy
is not a microscopic OPH admission theorem, a kinetic moment inequality is
not a distribution-function existence proof, and complete unwrapped data is
stronger than measured finite images. The issue's explicit bounded
underdetermination exit is satisfied; its full physical-source and likelihood
requirements remain open. The new claim's contextual M1 ancestry carries its
own assumptions without transferring dense-radius premises.
