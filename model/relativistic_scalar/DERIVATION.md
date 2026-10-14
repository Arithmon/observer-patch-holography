# Canonical scalar field, observer records and comparison domains

This packet changes the field action and occupation domain while reusing the
existing twelve-port compiler. Its nearest predecessors are the
[spatial Bose construction](../spatial_interaction/DERIVATION.md), which has an
exact two-particle sector, and the
[infinite two-oscillator history](../common_source/INFINITE_HISTORY.md), which
uses factorial-weighted polynomials and spectral remainders. The new result is
an eight-site, three-dimensional **canonical quartic** field with source
execution and exact complete-record comparison to unrestricted oscillators.
The abstract quantum-field interfaces do not supply its action or geometry.
The [SourceLaw covariance identification](../../docs/research/SOURCE_LAW_HISTORY_REGISTER.md)
is a separate free-field reference under its specified normalization and
Euclidean-history hypotheses. Its Gaussian vacuum statement is not an
interacting-vacuum theorem or an identification of its time parameter with
this packet's Hamiltonian coordinate.

## 1. Fix the target and operational meaning first

For a finite cubic graph Lambda, with positive spacing a, mass m and speed c,
and lambda >= 0, define on L2(R^|Lambda|)

    H = (1/2) sum_x p_x^2 + (m^2/2) sum_x q_x^2
        + (c^2/(2a^2)) sum_<xy> (q_x-q_y)^2
        + (lambda/(24a^3)) sum_x q_x^4,
    [q_x,p_y] = i delta_xy,  q_x = a^(3/2) phi_x.

The nonnegative closed energy form specifies its Friedrichs self-adjoint
Hamiltonian. Its confining potential bounds m^2 |q|^2/2 from below, giving
compact resolvent in every finite spatial box. Finite oscillator polynomials
times the Gaussian lie in every required Hamiltonian-power domain: polynomial
differentiation preserves the Schwartz class. The spectral Taylor estimates
below need these finite moments, not an assertion that their infinite Taylor
series converges for quartic evolution.

This is the standard spatial regulator for scalar field theory; see
[Jordan, Lee and Preskill](https://arxiv.org/abs/1112.4833). With continuous
Hamiltonian **model time**, its untruncated Heisenberg equations give

    d_t^2 phi_x = c^2 Delta_a phi_x - m^2 phi_x - (lambda/6) phi_x^3.

The nonlinear Klein--Gordon form identifies the intended relativistic
effective action. A finite spatial regulator preserves cubic symmetry and
does not make exact Lorentz or spacelike-commutator claims. Finite occupation
matrices additionally have a boundary correction to the canonical commutator.

The executed cell is {0,1}^3, with all twelve nearest-neighbor bonds, and
m=c=a=lambda=hbar=1. The local oscillator basis has frequency Omega=2.
Each vertex has degree three, so its quadratic onsite term is 2n+1. This
basis frequency is neither a measured gap nor an independent physical clock.
The occupation boxes are (8,2,2,2,2,2,2,2) and (8,4,4,4,4,4,4,4).
They are approximations, not invariant sectors. The Hamiltonian and geometry
are supplied; native repair does not select their coefficients here.

Before computing any result, define this observation protocol:

1. Prepare the product Omega=2 oscillator blank. It is not the interacting
   ground state, or the coupled free-field ground state.
2. Evolve for model time t=1/128.
3. Record the two-outcome local instrument (I-P_s, P_s), where
   e_s=(|0>+s i sqrt(24)|4>)/5, P_s=|e_s><e_s| at the origin, and s is -1 or +1.
4. If the committed bit is one, multiply the origin's fourth occupation by i.
5. Continue the same disturbed field for another t=1/128, then record the
   same two-outcome instrument. Retain all four joint outcomes.

The first native repair current selects s=+1 exactly when j1>0. The last
field bit selects the next E/F native repair. The specified 66-state prior
and both fair repair coins give 1,056 joint histories at each cutoff. This
attachment is classical on that prior; it does not replace the separate
coherent full-carrier construction. Blank preparation, routing, coherent
load/export and pointer readout remain supplied services.

Compatibility is specified independently of the executor: embed finite
occupation states into the same oscillator Hilbert space, and compare these
same instruments, retained records and feedback operations. Both P_s and
the feedback preserve every chosen occupation box. No outcome is normalized
or discarded. Model-time steps approximate the stipulated evolution;
source pulse counts and host runtime do not calibrate physical duration.

## 2. Correct compression and actual local source words

Let Pi_d be the rectangular occupation projector. The finite Hamiltonian is
H_d=Pi_d H Pi_d. Form each full polynomial before compressing it. In
particular, Pi_d q^4 Pi_d differs from (Pi_d q Pi_d)^4. The latter would
silently remove excursions through levels above the cutoff even when both
endpoints are retained. Four extra local levels suffice to calculate every
retained quartic matrix element exactly.

Write H_d=A+B, where A is the sum of onsite terms and
B=-sum_<xy> q_x,d q_y,d. Onsite terms on distinct sites commute. Every bond
term also commutes with every other bond term, including at a shared vertex,
because they use the same q matrix there. Each symmetric step therefore is

    S(dt) = exp(-i dt A/2) exp(-i dt B) exp(-i dt A/2).

The onsite eigensystem is compiled through the existing real-orthogonal and
diagonal-phase routines. Bond words compile the two local q eigensystems
separately, then the product-eigenvalue diagonal phases. Pointer words
rotate the declared e_s to a basis state, perform a controlled pointer flip,
and restore the local field basis. Every retained matrix is obtained by
contracting emitted source pulses and checked against an independent target
on all local code inputs, with absolute phase retained. The global numerical
reference exponentiates the full sparse H_d, without this split schedule.

For a conservative exact bound, subtract the onsite scalar 1+lambda/128
from every onsite term. Its total phase is identical in the exact and
product evolutions. Let

    S = sum_x ||h_x - (1+lambda/128)I||
        + sum_<xy> ||q_x,d|| ||q_y,d||.

The implementation bounds each Hermitian onsite norm by its absolute row
sums in the normalized Fock basis, with rational outward square-root
enclosures. It uses ||q_d|| <= sqrt(d-1). The exact and symmetric-product
curves agree through second order at zero; their third derivatives have
norm at most S^3. Taylor's integral remainder and a unitary telescoping sum
give

    || exp(-itH_d) - S(t/n)^n || <= S^3 t^3/(3n^2).

Errors of the two continued intervals add through the complete instrument
isometries. Thus the complete record TV error is at most twice this bound.
This argument concerns the exact compiled unitaries. Compiler arithmetic,
angle synthesis and physical service errors are separate; floating residuals
and the 2e-8 software guard are diagnostics, not interval certificates.

## 3. Exact comparison to eight unrestricted oscillators

Use Bargmann monomials z^n=(a_dagger)^n|0>, whose squared norms are prod(n_x!).
Here q=(a+a_dagger)/2 and the centered Hamiltonian is

    H' = H - (8+lambda/16) I
       = 2N - (1/4) sum_<xy> B_x B_y
         + (lambda/384) sum_x (B_x^4-3),  B=a+a_dagger.

Its action has rational coefficients. Creation increments occupation and
annihilation differentiates the monomial. Hamiltonian powers expand their
entire finite support; no occupation truncation is applied to the infinite
target. For the finite comparison, project only after each complete H'.
The same scalar is removed on both sides, so full histories differ by no
relative scalar phase.

For a finite-support input psi and Taylor degree r,

    || exp(-itH') psi - sum_(j=0)^r (-itH')^j psi/j! ||
        <= |t|^(r+1) ||(H')^(r+1) psi||/(r+1)!.

This follows from the scalar spectral remainder for real energies. Every
polynomial coefficient and squared norm is computed exactly; rational upper
square roots are verified by exact squaring. Degree two suffices for the
retained experiment. The local instrument has rational monomial coefficients
because e_s=(1+s i z_0^4)/5 and ||z_0^4||^2=24. Feedback also acts exactly.

Let e1 bound the first segment and e20,e21 bound the second segments on
their respective unnormalized polynomial branches. The complete instrument
is an isometry, hence the combined retained-vector error is bounded by

    epsilon = e1 + sqrt(e20^2+e21^2).

For a final branch with polynomial squared norm p and branch-vector error e,
the exact probability lies in [p-2sqrt(p)e-e^2, p+2sqrt(p)e+e^2], intersected
with [0,1]. Polynomial branch weights are not renormalized. Their half-L1
distance from the true record law is at most epsilon+epsilon^2/2.

Compare the finite and unrestricted polynomial histories in the same
factorial-weighted space. Their vector difference plus their two Taylor
budgets bounds the exact retained-vector difference, and therefore both
the record TV distance and the complete retained-output trace distance.
Adding the source product bound then controls the actual ideal source law.
The common record-controlled native attachment cannot increase this error.
The certificate covers **eight untruncated oscillators**, not an infinite
spatial lattice. Its finite Taylor remainder is not a claim of a vanishing
error at fixed degree as the occupation box grows.

## 4. A bounded, strictly positive quartic witness

A separate one-read protocol starts from the same product blank, evolves
for t=1/4096 and reads P_- at the origin. At t=0 its probability is 1/25.
The quartic vacuum-to-four matrix element gives

    d_t [p_lambda=1(t)-p_lambda=0(t)] at zero = 1/200.

The quadratic action alone has no vacuum-to-four matrix element at first
order. In monomial coordinates H'|vac> has origin-four coefficient 1/384;
the weighted projector matrix element supplies the stated derivative.
This is a bounded local instrument, rather than an unbounded moment used
as a proxy for an observer record.

For either centered Hamiltonian K and normalized initial psi, its probability
remainder is at most t^2 (||K psi||^2+||K^2 psi||). Indeed the absolute second
derivative is at most twice that expression, by Cauchy--Schwarz and ||P||=1.
Exact first and second moments give a positive contrast lower bound for
the unrestricted action and for both executed occupation boxes. Subtracting
both finite source product errors still leaves a positive executed margin.
The result establishes a response to this declared quartic interaction;
it is not a measured relativistic scattering amplitude or a fit to nature.

## 5. Gap of the untruncated finite-volume target

The configuration-space potential satisfies

    Hess V = m^2 I + (c^2/a^2)L_graph
             + diag(lambda q_x^2/(2a^3)) >= m^2 I.

Apply [Andrews--Clutterbuck, Theorem 1.3 and Section 5](https://arxiv.org/abs/1006.1686)
to 2H=-Delta+2V, using the one-dimensional comparison potential m^2 s^2.
Its convexity modulus is dominated by that of 2V. On expanding
configuration-space balls, their comparison bounds the gap by the
corresponding harmonic-oscillator gap. Dirichlet exhaustion gives

    E1(H)-E0(H) >= m.

The bound is uniform in finite site count, positive spacing and nonnegative
coupling. Passing it to an infinite-spatial-volume vacuum representation
requires constructing that representation and its dynamics. It also does
not transfer automatically to a finite occupation matrix. The mass m is
supplied, and this protocol neither prepares the interacting ground state
nor measures its gap. For lambda>0 an exactly Gaussian ground wavefunction
cannot solve the Schrödinger equation: its kinetic contribution has degree
at most two, leaving an uncancelled positive quartic polynomial.

## 6. Free relativistic comparison and the remaining interacting obligation

The infinite or periodic free spatial lattice has frequency

    omega_a(k)^2 = m^2 + (4c^2/a^2) sum_j sin^2(a k_j/2).

This differs from the earlier Bose action's quadratic energy above m. For a
common Fourier band |k|<=K with aK<=1,

    0 <= omega(k)-omega_a(k) <= c^2 a^2 K^4/(24m),
    omega(k) = sqrt(m^2+c^2|k|^2).

Use 0<=x^2-4sin^2(x/2)<=x^4/12 and divide the squared-frequency difference
by omega+omega_a>=2m. The bound holds in every momentum direction. On the
common band-limited free Fock representation with at most N particles,
the evolution norm difference is at most N|t| times this frequency bound.
The equal-time covariance multiplier difference is at most the frequency
bound divided by 2m^2. This is an analytic free-field comparison; it includes
neither the eight-site boundary nor an interacting mass/coupling correction.

The finite scalar regulator is a useful relativistic EFT candidate. A
controlled **interacting** validity domain still needs matched renormalized
parameters and bounds on its observable laws, including spatial boundaries,
momentum/energy restrictions and Lorentz-breaking corrections. The free
bound cannot discharge that obligation. Exact spacelike locality, a
continuum scattering theory, calibrated physical clocks, action selection
and comparison with independent measurements remain separate work.

Nor is positive scalar phi4 a promise of nontrivial ultraviolet completion:
[Aizenman--Duminil-Copin](https://arxiv.org/abs/1912.07973) prove Gaussian
scaling limits for the stated nearest-neighbor ferromagnetic four-dimensional
Ising/phi4 families. This leaves finite-cutoff effective use available; it
prevents treating ordinary cutoff removal as an established route to an
interacting four-dimensional continuum here.
